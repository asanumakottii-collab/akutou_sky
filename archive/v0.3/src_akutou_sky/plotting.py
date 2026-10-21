"""Reproducible standalone scientific figure from the labeled v2 NetCDF."""
import json
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import colors, font_manager

def plot_profiles(source,out,*,azimuth=0.):
    """Save a scientific PNG from a Dataset or NetCDF path; return the output path.

    Uses a private Agg canvas without changing the application's plotting backend
    or global rcParams. Install akutou-sky[plot] to use this optional module.
    """
    from .io import load_dataset
    from .colorimetry import add_colorimetry
    dataset = source if isinstance(source,xr.Dataset) else load_dataset(source)
    if 'total_XYZ' not in dataset:
        dataset = add_colorimetry(dataset)
    ds=dataset.sel(azimuth=azimuth).sortby('depression').sortby('elevation')
    fonts={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in ['Hiragino Sans','Noto Sans CJK JP','Yu Gothic','Meiryo','Arial Unicode MS'] if f in fonts),'DejaVu Sans')
    with matplotlib.rc_context({'font.family':font,'font.size':11,'axes.spines.top':False,'axes.spines.right':False}):
        d,e=ds.depression.values,ds.elevation.values
        # No invented intermediate depressions. Unequal elevation spacing uses actual cell edges.
        if len(e)<2 or len(d)<2:raise ValueError('Figure requires at least two elevations and depressions')
        yedges=np.r_[e[0],(e[:-1]+e[1:])/2,e[-1]]
        xedges=np.r_[d[0]-(d[1]-d[0])/2,(d[:-1]+d[1:])/2,d[-1]+(d[-1]-d[-2])/2]
        rgb=ds.total_preview_sRGB.values.transpose(1,0,2)
        Y=ds.total_XYZ.sel(tristimulus='Y').values.T
        fig=Figure(figsize=(13,6.8))
        FigureCanvasAgg(fig)
        axes=fig.subplots(1,2,gridspec_kw={'width_ratios':[1,1.14]})
        fig.subplots_adjust(left=.065,right=.94,bottom=.24,top=.80,wspace=.24)
        axes[0].pcolormesh(xedges,yedges,rgb,shading='flat',rasterized=True)
        positive=Y[Y>0]
        norm=colors.LogNorm(vmin=10**np.floor(np.log10(positive.min())),vmax=10**np.ceil(np.log10(positive.max())))
        im=axes[1].pcolormesh(xedges,yedges,np.ma.masked_less_equal(Y,0),shading='flat',cmap='cividis',norm=norm,rasterized=True)
        fig.colorbar(im,ax=axes[1],label='輝度 Y [cd/m²]',fraction=.045,pad=.035)
        axes[0].set_title('色の分布：各点を明るく正規化',pad=12)
        axes[1].set_title('光量の分布：共通の対数目盛',pad=12)
        for ax in axes:
            ax.set_ylim(e[0],e[-1]);ax.set_xlim(xedges[0],xedges[-1]);ax.set_xticks(d,[f'{x:g}' for x in d])
            ax.set_yticks([x for x in [0,5,15,30,45,60,75,90] if e[0]<=x<=e[-1]])
            ax.set_xlabel('太陽中心の伏角 [°]');ax.set_ylabel('地平線からの高度 [°]')
            ax.axhspan(0,5,facecolor='none',edgecolor='#404040',hatch='////',lw=.3)
            if d.min()<=12<=d.max():ax.axvline(12,color='white',ls='--',lw=1)
        config=json.loads(ds.attrs['config_json'])
        direction='太陽方向' if azimuth==0 else f'方位差{azimuth:g}°'
        physics=f"球面多重散乱・{config['aerosol_model'].upper()}エアロゾル・太陽光の屈折{'あり' if config['solar_refraction'] else 'なし'}"
        fig.suptitle(f'あくとう v2  |  晴天・{direction}の色と明るさ',x=.065,ha='left',fontsize=20,y=.96,color='#243448')
        fig.text(.065,.885,physics+f"  /  {config['wavelength_step_nm']:g} nm区間の分光計算",fontsize=11,color='#526170')
        fig.text(.065,.15,'数値収束未達の試算です。定量精度は未確立。左図は各点を正規化し、光量や暗順応の見え方を表しません。',fontsize=10)
        background='夜空背景を含まない残光です。' if ds.attrs['background']=='none' else '指定された背景スペクトルを含みます。'
        fig.text(.065,.10,'高度5°未満（斜線）は参考領域。伏角12°以降（破線）は'+background+' 雲・太陽円盤は未モデル化。',fontsize=10)
        source_label='in-memory Dataset' if isinstance(source,xr.Dataset) else Path(source).name
        fig.text(.065,.05,f'出典：{source_label}。列は選択した計算点の代表区間。区間内の連続変化を計算した図ではありません。',fontsize=9,color='#59636c')
        out=Path(out)
        out.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(out,dpi=160)
        return out

