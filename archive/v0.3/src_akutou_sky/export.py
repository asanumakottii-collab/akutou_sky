"""Export any v2 run without hardcoded elevation, angle or wavelength counts."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import xarray as xr
from .colorimetry import band_mean, display_values, xyz_from_radiance

def color_values(r):
    """Legacy NPZ compatibility for independent comparisons."""
    xyz = xyz_from_radiance(r['radiance'], r['wavelength_nm'])
    xy, rgb, gamut = display_values(xyz)
    return xyz, xy, rgb, gamut

def csv_write(path, header, rows):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for row in rows:
            writer.writerow(['' if isinstance(v, (float, np.floating)) and not np.isfinite(v) else v for v in row])

def export(ds, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    d, az, e, w = [ds[k].values for k in ['depression', 'azimuth', 'elevation', 'wavelength']]
    xyz, xy = ds.total_XYZ.values, ds.total_xy.values
    rgb = ds.total_preview_sRGB.values
    solar = ds.solar_XYZ.values
    gamut = ds.total_gamut_clipped.values
    coords = [(i,j,k) for i in range(len(d)) for j in range(len(az)) for k in range(len(e))]
    csv_write(out/'colors.csv',
        ['solar_depression_deg','relative_azimuth_deg','elevation_deg','X','Y_cd_m2','Z','CIE_x','CIE_y',
         'preview_sRGB_R8','preview_sRGB_G8','preview_sRGB_B8','gamut_clipped',
         'solar_Y_cd_m2','background_Y_cd_m2','near_horizon','night_background_important'],
        ([d[i],az[j],e[k],*xyz[i,j,k],*xy[i,j,k],*np.rint(255*rgb[i,j,k]).astype(int),
          int(gamut[i,j,k]),solar[i,j,k,1],xyz[i,j,k,1]-solar[i,j,k,1],int(e[k]<5),int(d[i]>=12)] for i,j,k in coords))
    spec = ds.solar_stokes_radiance.values
    total = ds.total_radiance.values
    bounds = ds.wavelength_bounds.values
    csv_write(out/'spectra.csv',
        ['solar_depression_deg','relative_azimuth_deg','elevation_deg','wavelength_nm','bin_left_nm','bin_right_nm',
         'solar_I_W_m2_sr_nm','total_I_W_m2_sr_nm'] + [f'solar_{s}_W_m2_sr_nm' for s in ds.stokes.values[1:]],
        ([d[i],az[j],e[k],ww,*bounds[l],spec[i,j,k,l,0],total[i,j,k,l],*spec[i,j,k,l,1:]]
         for i,j,k in coords for l,ww in enumerate(w)))
    bands = []
    for lower, upper in [(0,5),(5,15),(15,35),(35,60)]:
        if e.min() <= lower < upper <= e.max():
            mean = band_mean(e, xyz, lower, upper)
            mxy, _, _ = display_values(mean)
            for i, dd in enumerate(d):
                for j, aa in enumerate(az):
                    bands.append([dd,aa,lower,upper,*mean[i,j],*mxy[i,j]])
    csv_write(out/'bands.csv', ['solar_depression_deg','relative_azimuth_deg','h_min_deg','h_max_deg',
                              'X','Y_cd_m2','Z','CIE_x','CIE_y'], bands)
    summary = {
        'model_version': ds.attrs['model_version'], 'shape': dict(ds.sizes),
        'Y_min_cd_m2': float(xyz[...,1].min()), 'Y_max_cd_m2': float(xyz[...,1].max()),
        'gamut_clipped_count': int(gamut.sum()), 'background': ds.attrs['background'],
        'convergence': ds.attrs['convergence'], 'elapsed_seconds': ds.attrs['elapsed_seconds'],
        'limitations': ['地平線付近は視線・散乱光の屈折が未対応。',
                        '太陽円盤の有限サイズ、雲、水平方向の大気変化は未モデル化。',
                        'O2・H2Oの狭い吸収線、NO2吸収は未実装。',
                        '深い薄明の総光量には実測の夜空背景が必要。',
                        '粒径・屈折率は仮定。実測校正・全天での数値収束保証なし。'],
    }
    if ds.sizes['stokes'] == 3:
        summary['limitations'].append('偏光は試験機能。粗い格子で左右対称性検査が不合格。Q/UはSASKTRAN2既定のStandard基底（全球z軸と視線の面）。')
    validation_path = out/'validation.json'
    if validation_path.exists() and (out/'sky.nc').exists():
        validation = json.loads(validation_path.read_text())
        with (out/'sky.nc').open('rb') as stream:
            current_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
        if validation.get('input_sha256') == current_hash:
            summary['validation'] = validation
    (out/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    write_viewer(ds, out/'index.html', summary)
    return out

def write_viewer(ds, path, summary):
    payload = {key: ds[key].values.tolist() for key in ['depression','azimuth','elevation','total_XYZ','total_preview_sRGB']}
    payload['config'] = json.loads(ds.attrs['config_json'])
    payload['summary'] = summary
    # No network dependencies; safe even if a metadata string contains </script>.
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    template = '''<!doctype html><html lang="ja"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>あくとう · 大気散乱モデル v2</title>
<style>
:root{color-scheme:dark;font-family:-apple-system,BlinkMacSystemFont,"Noto Sans JP",sans-serif;background:#0b1420;color:#e8edf4}
body{max-width:1100px;margin:auto;padding:32px 22px}h1{font-size:28px;margin-bottom:8px}p{line-height:1.8;color:#bdc9d8}
.controls,.cards{display:flex;gap:18px;flex-wrap:wrap;background:#172334;padding:20px;border-radius:12px;margin:20px 0}
label{display:flex;gap:8px;align-items:center;flex-wrap:wrap}select,input{font:inherit;accent-color:#7dc7d1}select{padding:6px;background:#233248;color:inherit;border:1px solid #60748c;border-radius:5px}
.card{flex:1;min-width:280px}canvas{width:100%;background:#0d1928;border-radius:8px}h2{font-size:18px}strong{color:#90d7e0}a{color:#90d7e0}pre{white-space:pre-wrap;font-size:12px}small{color:#acbbcd}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:7px;text-align:right;border-bottom:1px solid #34465e}
</style>
<h1>あくとう · 朝焼け・夕焼け・薄明</h1>
<p>球面大気の多重散乱から求めた分光放射輝度と色。太陽の伏角を変えて、空の高度方向の分布を確認できます。</p>
<div class="controls"><label>太陽の伏角 <select id="dep"></select></label><label>太陽からの方位角 <select id="az"></select></label><label>表示 <select id="mode"><option value="color">色のみ（各点を正規化）</option><option value="exposure">共通の露出で表示</option></select></label><label>表示基準 <input id="exposure" type="range" min="-5" max="5" step="0.1" value="2"><span id="exptext"></span> cd/m²</label></div>
<p id="state" aria-live="polite"></p>
<p id="numerics"></p><p id="validationStatus"></p>
<div class="cards"><div class="card"><h2>高度ごとの色</h2><canvas id="sky" width="470" height="500"></canvas><p id="colorNote"></p></div>
<div class="card"><h2>輝度 Y [cd/m²] · 対数目盛</h2><canvas id="plot" width="470" height="500"></canvas><p>色のみの表示では光量の差を取り除いています。深い薄明でも肉眼に鮮やかな色が見えるという意味ではありません。</p></div></div>
<p id="limits"></p><details><summary>代表点の数値</summary><table><thead><tr><th>高度</th><th>Y [cd/m²]</th><th>x</th><th>y</th></tr></thead><tbody id="table"></tbody></table></details>
<details><summary>計算条件・モデルの限界</summary><pre id="config"></pre></details>
<p><a href="colors.csv">色と輝度 CSV</a> · <a href="spectra.csv">分光 CSV</a> · <a href="bands.csv">4帯域平均 CSV</a> · <a href="sky.nc">全計算結果 NetCDF</a> · <a href="manifest.json">再現条件</a></p>
<small>太陽伏角は幾何学的な太陽中心の高度に負号を付けた値。朝焼けは伏角が減る向き、夕焼けは増える向きです。同じ大気条件での対応です。画面のRGBをLEDのPWMには直接使用できません。</small>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent), el=id=>document.getElementById(id);
for (const [id,values] of [['dep',D.depression],['az',D.azimuth]]) values.forEach((v,i)=>{let o=document.createElement('option');o.value=i;o.textContent=v+'°';el(id).appendChild(o)});
const encode=x=>x<=.0031308?12.92*x:1.055*Math.pow(x,1/2.4)-.055;
function draw(){let i=+el('dep').value,j=+el('az').value,ref=10**(+el('exposure').value),mode=el('mode').value;
el('exptext').textContent=ref.toPrecision(3);let order=D.elevation.map((_,k)=>k).sort((a,b)=>D.elevation[a]-D.elevation[b]);
let xyz=D.total_XYZ[i][j],rgb=D.total_preview_sRGB[i][j],sky=el('sky').getContext('2d'),plot=el('plot').getContext('2d');
for(let c of [sky,plot]){c.clearRect(0,0,470,500);c.font='14px sans-serif';c.fillStyle='#bdc9d8'}
let emin=D.elevation[order[0]],emax=D.elevation[order.at(-1)],span=Math.max(emax-emin,1),y=e=>455-(e-emin)/span*425;
order.forEach((k,q)=>{let low=q===0?emin:(D.elevation[k]+D.elevation[order[q-1]])/2,high=q===order.length-1?emax:(D.elevation[k]+D.elevation[order[q+1]])/2,c=rgb[k];
if(mode==='exposure'){let [X,Y,Z]=xyz[k];c=[3.2406*X-1.5372*Y-.4986*Z,-.9689*X+1.8758*Y+.0415*Z,.0557*X-.2040*Y+1.0570*Z].map(v=>encode(Math.max(0,Math.min(1,v/ref))))}
sky.fillStyle='rgb('+c.map(v=>Math.round(v*255)).join(',')+')';let top=Math.floor(y(high)),bottom=Math.ceil(y(low));sky.fillRect(75,top,350,Math.max(1,bottom-top)) });
for(let e of [0,5,15,30,45,60,75,90])if(e>=emin&&e<=emax){sky.fillStyle='#d3dfed';sky.fillText(e+'°',20,y(e)+4)}
let vals=xyz.map(v=>Math.log10(Math.max(v[1],1e-20))),lo=Math.floor(Math.min(...vals)),hi=Math.max(lo+1,Math.ceil(Math.max(...vals))),x=v=>65+(v-lo)/(hi-lo)*370;
for(let t=lo;t<=hi;t+=Math.max(1,Math.ceil((hi-lo)/6))){plot.strokeStyle='#263951';plot.beginPath();plot.moveTo(x(t),30);plot.lineTo(x(t),455);plot.stroke();plot.fillStyle='#bdc9d8';plot.fillText('10^'+t,x(t)-15,481)}
for(let e of [0,15,30,45,60,75,90])if(e>=emin&&e<=emax)plot.fillText(e+'°',8,y(e)+4);
plot.strokeStyle='#91d7df';plot.lineWidth=2;plot.beginPath();order.forEach((k,n)=>n?plot.lineTo(x(vals[k]),y(D.elevation[k])):plot.moveTo(x(vals[k]),y(D.elevation[k])));plot.stroke();
el('state').textContent='伏角 '+D.depression[i]+'° / 方位差 '+D.azimuth[j]+'° / 輝度 '+Math.min(...xyz.map(v=>v[1])).toExponential(2)+' ～ '+Math.max(...xyz.map(v=>v[1])).toExponential(2)+' cd/m²';
let comparisons=D.summary.validation?Object.values(D.summary.validation.checks).flatMap(c=>(c.per_depression||[]).filter(r=>r.depression_deg===D.depression[i])):[];
el('numerics').textContent=comparisons.length?'数値感度：設定を細かくした比較での最大輝度差 '+(100*Math.max(...comparisons.map(c=>c.all.max_relative_Y))).toFixed(2)+'%。検証した高度での比較差であり、誤差保証ではありません。':'この伏角の数値収束は未確認です。';
el('colorNote').textContent=mode==='color'?'各高度を別々に明るくした色比較です。':'共通の露出による表示です。表示範囲を超えた光はクリップされます。';
el('table').replaceChildren();order.filter((k,q)=>q%Math.max(1,Math.floor(order.length/12))===0||q===order.length-1).forEach(k=>{let tr=document.createElement('tr'),v=xyz[k],s=v.reduce((a,b)=>a+b,0);[D.elevation[k]+'°',v[1].toExponential(4),s?(v[0]/s).toFixed(5):'未定義',s?(v[1]/s).toFixed(5):'未定義'].forEach(t=>{let td=document.createElement('td');td.textContent=t;tr.appendChild(td)});el('table').appendChild(tr)});
}
el('config').textContent=JSON.stringify(D.config,null,2)+'\\n\\n'+D.summary.limitations.join('\\n');
if(D.summary.validation)el('config').textContent+='\\n\\n数値感度の検証：\\n'+JSON.stringify(D.summary.validation,null,2);
el('validationStatus').textContent=D.summary.validation?(D.summary.validation.passes_all_sampled_checks?'実施した数値感度の比較は基準内です。実測精度の保証ではありません。':'数値感度の比較に基準外の項目があります。定量的な予測値としての精度は未確立です。'):'数値感度の検証結果はこのデータに紐付けられていません。';
el('limits').textContent='高度5°未満は参考値。'+(D.summary.background==='none'?'夜空背景は含まれていません。':'指定された背景スペクトルを加算しています。')+' 数値収束の判定は validation.json を確認してください。';
for(let id of ['dep','az','mode','exposure'])el(id).addEventListener('input',draw);draw();
</script></html>'''
    Path(path).write_text(template.replace('__DATA__', data), encoding='utf-8')

