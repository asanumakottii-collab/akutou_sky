# あくとう：大気散乱による朝焼け・夕焼け・薄明 v2

`./sky.sh` または `akutou-sky simulate` で計算し、`results/<tag>/index.html` で色と輝度を確認します。作成済みの主データは `results/improved/`（太陽方向）、3方位の軽量試算は `results/overview/` にあります。既定は晴天・海抜2 mの仮想的な大気です。実測校正済みの予報やLEDのPWM値ではありません。

**検証結果：物理過程と分光積分は拡張しましたが、主データの数値収束は未達です。** 波長・反復の比較差は小さい一方、散乱角・水平・鉛直格子の比較には基準外の点があり、伏角6°でも大きな差が残ります。偏光の粗い格子では左右対称性の検査も不合格です。現状は研究用の試算として扱ってください。詳細は `改善内容と検証.md` に記録しています。

## 今回の変更

- 球面2Dの地球の影と多重散乱を維持し、3成分の偏光計算（I, Q, U）も選択可能に。通常設定はメモリを抑えたスカラー計算。
- 霞に **Mie散乱** を導入。対数正規の粒径分布と複素屈折率から、波長ごとの消散・吸収・散乱行列を計算。従来のHG近似も選択可能。
- **太陽入射光の屈折**を追加。Ciddorの乾燥空気の屈折率を使用。視線と多重散乱光の屈折は現行エンジンの2D経路が未対応。
- 波長は380–780 nmの区間を5 nm幅で覆い、区間平均太陽スペクトルと区間積分したCIE等色関数を使用。端の波長の重複や欠落を避ける。
- 下層大気は125 m、10–30 kmは500 m刻み。視線の高度0–5°は0.25°刻み、5–90°は1°刻み。物性・光源場・視線の刻みは別々の設定。
- 通常設定の光源場は500 m / 1 km / 4 km刻み、水平41点、散乱角38点。水平格子は太陽天頂角85–110°付近に集中させ、変化の急な薄明を詳しく扱います。メモリ節約のため旧版より粗い方向もあり、各刻みの影響は数値検証で評価します。
- オゾン全量、霞の量・高さ・粒径、地表反射率、観測高度、太陽距離、成層圏エアロゾルを設定可能。
- 方位角を追加。0°は太陽側、90°は横、180°は反太陽側。朝夕は同じ大気を仮定する場合、伏角を逆順に利用。
- NetCDFに座標・単位・設定を保存し、CSVとオフラインHTMLを同時生成。旧来の固定行数への依存を解消。
- 各実行のソース・物理データ・環境のハッシュを保存。数値収束は別スクリプトで検証し、変化をそのまま報告。

## 実行

Python 3.13で検証しています。Mac用の `setup_environment.sh` は、同期フォルダでのライブラリ読み込み待ちを避けるため、`~/Library/Caches/akutou-science/python313` に専用環境を作ります。システムPythonは変更しません。

```sh
cd /Users/asanumawatarushi/Documents/東大天文部/あくとう/output/akutou_sky
sh setup_environment.sh  # 初回のみ
./sky.sh --tag my_sky
```

新しい環境では次の手順です。

```sh
sh setup_environment.sh
./sky.sh --tag my_sky
```

最初は `--preset preview` で条件を探り、通常の試算は `standard`、数値確認には `akutou-sky validate` を使ってください。偏光は試験機能で、`--preset polarized` または `--num-stokes 3` で有効になります。偏光と高密度の空間・散乱角格子の組合せは数十GB以上のメモリを使うことがあるため、16 GB環境ではpolarizedプリセットから始めてください。プリセット名だけで精度は保証されません。Q/Uの基底はSASKTRAN2既定のStandard（全球z軸と視線の面）です。3方位を手早く見渡す場合は `--preset overview --azimuths 0 90 180` を使ってください。特に地平線付近は数値感度と未実装の物理過程の両方が残ります。計算は順番に実行してください。

```sh
# 軽い試算
./sky.sh --tag quick --preset preview --d 0 6 12
# オゾン・霞・高度を変更（単位はオプション名に含む）
./sky.sh --tag high_observer --aod550 0.03 --ozone-du 330 --observer-altitude-m 2000
# 従来のHG近似の条件を指定
./sky.sh --tag hg --aerosol-model hg --aerosol-g 0.65 --aerosol-ssa 0.95 --angstrom-exponent 1.3
# 同じ物理条件でスペクトル・鉛直・水平・散乱角・反復を独立に比較
"$HOME/Library/Caches/akutou-science/python313/bin/python" -m akutou_sky validate --input results/my_sky/sky.nc --cache-dir model_cache
# 保存済みの計算から出力だけを再生成
"$HOME/Library/Caches/akutou-science/python313/bin/python" -m akutou_sky export results/my_sky/sky.nc --output results/my_sky --overwrite
# 自動テスト
"$HOME/Library/Caches/akutou-science/python313/bin/python" -m unittest discover -s tests -v
```

実行環境のPythonは `AKUTOU_PYTHON` でも指定できます。全項目は `akutou-sky simulate --help` または `src/akutou_sky/config.py` を参照。`--config scenario.json` にModelConfigのキーを記述すると繰り返し実行できます。優先順位は「標準値 → プリセット → JSON → コマンドライン」です。未知のキー、負の光学的厚さ、不正な波長刻み等は計算前に拒否します。同名結果の上書きは `--overwrite` が必要です。光源場の高度数×水平数×散乱角数×Stokes成分数に上限を設け、極端なメモリ負荷を事前に検出します。`--max-source-samples` は十分なメモリがある機械でのみ増やしてください。この個数による判定はRAMの正確な予測ではありません。

観測視線数および大気物性の高度数との積にも上限を設けています。`--max-view-source-product` と `--max-atmosphere-source-product` は資源上限の変更用で、物理的な精度設定ではありません。

## 出力

| ファイル | 内容 |
|---|---|
| `preview.png` | 主データの色と輝度を一覧する図（`akutou-sky plot` で生成） |
| `index.html` | 伏角・方位の切替。色だけの表示と共通露出の表示。輝度の対数グラフ |
| `sky.nc` | 分光I/Q/U、太陽成分・背景・合計、XYZ、座標、単位、気圧・温度・オゾン、計算設定 |
| `colors.csv` | 絶対XYZ・xy・輝度・正規化プレビューsRGB・注意フラグ |
| `spectra.csv` | 波長区間ごとの放射輝度 W/(m² sr nm)。Q/Uは負値も許す |
| `bands.csv` | 0–5、5–15、15–35、35–60°の立体角重み付きXYZ平均。指定範囲に含まれる帯域のみ |
| `manifest.json` | 実行条件・版・ハッシュ・背景の出典ファイル |
| `validation.json` | 検証を実行した場合の数値感度。高度5°未満とそれ以上を分けて記録 |

NetCDFは次のように読み出せます。

```python
import xarray as xr
sky = xr.open_dataset('results/my_sky/sky.nc')
Y = sky.total_XYZ.sel(tristimulus='Y')
spectrum = sky.total_radiance.sel(depression=6, azimuth=0, elevation=15)
```

sRGBは各点の最大値で規格化した色比較です。暗順応時の知覚や発光体の制御値ではありません。ゼロ光量のxyは未定義（NetCDFではNaN、CSVでは空欄）、プレビューは黒とします。

## 背景光

夜光・星・街明かりを推測で青色に置き換える処理は行いません。必要なら、観測位置での分光放射輝度をCSVで与えます。

```text
wavelength_nm,radiance_W_m2_sr_nm
...
```

380–780 nmを完全に覆い、波長は昇順、放射輝度は非負にしてください。`--background measured_background.csv` で区間積分して太陽成分に加算します。これは**観測者位置で既に減衰を受けた、全天で一様・無偏光と仮定する背景**です。大気上端の放射や方位依存の街明かりをそのまま渡すことはできません。

## 精度の解釈

Mie計算は球形粒子と仮定した粒径・屈折率に対する物理計算です。既定値は特定の観測日の測定値ではなく、HGからMieへの変更だけで実空に近づいたことは証明できません。既存のimproved・overview・smokeは吸湿成長を含みません。モデル2.1から、指定したRHとκによる吸湿成長と霧層を選択的に計算できます。[霧と塵の説明](fog-and-dust.md) を参照してください。

`aerosol_ssa`、`aerosol_g`、`angstrom_exponent` はHGの場合のみ有効です。Mieではこれらを固定せず、`median_radius_nm`、`mode_width`、複素屈折率から計算します。逆にHGでは粒径と屈折率の設定を使いません。`scenarios/` に清澄・霞・成層圏層の仮想的な条件例を収録しています。

観測高度を変えても地表半径は海面基準のままです。高い観測者を配置する機能であり、山岳地形を再現する機能ではありません。

太陽円盤の有限サイズと直接光、雲、水平の大気変化、視線・散乱光の屈折、可視域のO2/H2Oの細い吸収線、NO2、暗順応は未実装です。深い薄明では夜空背景の欠落も支配的です。太陽伏角を単純に0.57°ずらすような一律の大気差補正はしていません。

`akutou-sky validate` は実際に計算した角度の一部を選び、物理条件を同じに保って刻みを変更します。水平・鉛直の細分化量はメモリ用の上限内で選び、実際の設定差を記録します。予算内で細分化できない検査はskippedとなります。相対輝度差5%・xy距離0.005を目安として記録しますが、測定誤差や誤差保証ではありません。未実施・不合格を収束済みと表示しません。

## 旧版との関係

直下の日本語CSV・PNG、`runs/*.npz`、`summary.json`、`provenance.json`、`科学的根拠と利用上の注意.md` および隣のPDFは **2026-10-06の旧版結果**です。v2の計算値は `results/` の中にまとめています。旧モデル・旧PDF生成スクリプトと互換アダプターは整理時に削除しました。現行モデルの実装は `src/akutou_sky/` に集約しています。

## 参照

- [SASKTRAN2 逐次散乱法・2D屈折の対応範囲](https://sasktran2.readthedocs.io/en/latest/components/source_terms/successive_orders.html)
- [SASKTRAN2 空気の屈折率・Ciddor式](https://sasktran2.readthedocs.io/en/latest/users_guide/refraction.html)
- [SASKTRAN2 Mie散乱データベース](https://sasktran2.readthedocs.io/en/latest/users_guide/mie.html)
- [CIE 1931 2°等色関数](https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer)
- [TSIS-1 Hybrid Solar Reference Spectrum v2](https://doi.org/10.1029/2022EA002637)

波長区間の中で輸送係数は中心波長の値で代表します。太陽スペクトルと等色関数を正しく積分しても、狭い吸収線を解像する分光モデルにはなりません。
