# akutou-sky — 朝焼け・夕焼け・薄明のモデリングライブラリ

球面大気の多重散乱から分光放射輝度、CIE XYZ・xy・輝度、表示用sRGBを計算するPythonライブラリです。SASKTRAN2を使い、Rayleigh散乱、MieまたはHGエアロゾル、オゾン吸収、太陽入射光の屈折を扱います。ライブラリ版は **0.3.0**、物理モデルの版は **2.2** です。2.1で微粒子の吸湿成長と独立した霧の水滴層を追加し、2.2で前方散乱の処理と数値分解能を改善しました。

**物理モデルは数値収束未達の研究用試算です。** 格子依存性が大きい条件と偏光の左右対称性の不合格が残っています。吸湿・霧を無効にした標準設定は従来の計算を維持します。[改善内容と検証](改善内容と検証.md) に従来モデルの検証を、[霧と塵の説明](docs/fog-and-dust.md) に追加機能を記録しています。

## 霧と塵を変える4種類の空

[霧と塵の比較ページ](results/atmosphere_explorer/index.html) は、都会・田舎・高湿度・低湿度を切り替え、塵4段階と霧4段階、朝・夕・薄明、見る方向、表示露出を調整できます。伏角は−6°〜＋6°、0.5°刻みの25点です。全64条件の計算結果を収めたオフラインHTMLで、物理的に同じ粒子なしの大気は再利用します（52種類の大気を計算）。色・輝度CSV、計算設定JSON、分光NetCDFも保存できます。条件は比較用の仮想例で、実測に基づく地域の代表値ではありません。

```sh
akutou-sky explore --output results/atmosphere_explorer --cache-dir model_cache --offline --resume
```

初めて作る保存先では `--resume` を省略できます。`--resume` はソース・設定・入力データのハッシュが一致する計算を再利用します。初回は数十分以上かかり、細分化した検証計算にも時間が必要です。より自由な値は次のように指定できます。

```sh
akutou-sky simulate --preset fog_refined --aod550 0.08 --fog-aod550 0.12 \
  --relative-humidity 0.85 --hygroscopicity-kappa 0.30 --d -6 -3 0 3 6 \
  --elevations 0 1 5 15 30 60 90 --azimuths 0 90 180 \
  --output results/custom_fog --cache-dir model_cache --offline
```

塵は乾燥状態の550 nm AOD、霧は水滴だけの550 nm AODです。相対湿度は0〜0.95の割合で与えます。吸湿効果は `hygroscopicity_kappa > 0` のとき有効です。`fog_refined` は低層格子と前方散乱処理を改良した設定です。旧設定の `fog_preview` も維持しています。いずれも精度の保証ではありません。

## インストール

Python 3.11–3.13に対応し、検証は3.13で行っています。このフォルダで実行します。

```sh
python -m pip install .             # 計算・色変換・NetCDF・HTML/CSV
python -m pip install '.[plot]'      # PNGの作図も使う場合
python -m pip install -e '.[plot]'   # ソースを編集しながら使う場合
```

既存のMac用環境は `sh setup_environment.sh` で依存パッケージと編集可能なライブラリをインストールします。`requirements.txt` は既存計算を検証したPython 3.13環境の固定版一覧です。ライブラリの直接の依存は `pyproject.toml` に定義しています。

## Pythonから使う

```python
from akutou_sky import ModelConfig, simulate, export_results

config = ModelConfig.from_preset("preview", aod550=0.05, ozone_du=330)
sky = simulate(
    depressions=[-2, 0, 6, 12],
    config=config,
    elevations=[0, 1, 5, 15, 30, 60, 90],
    azimuths=[0, 90, 180],
    cache_dir="model_cache",  # 既存データを再利用
    allow_download=False,
)
luminance = sky.total_XYZ.sel(tristimulus="Y")  # cd/m²
spectrum = sky.total_radiance.sel(depression=6, azimuth=0, elevation=15)
export_results(sky, "results/library_example")
```

`simulate()` はメモリ上の `xarray.Dataset` を返します。ログが必要なら `progress=print` を渡します。太陽成分のI/Q/Uのみが必要なら `compute_radiance()` を使います。

角度は幾何学的な角度、単位は度です。`depressions` は太陽中心の伏角（正なら地平線下、−18〜24°）、`elevations` は高度角（0〜90°）、`azimuths` は太陽からの相対方位（−180〜180°、0°が太陽側）です。朝焼けは伏角が減る順、夕焼けは増える順で与えます。入力の順番は維持されます。

`ModelConfig` は変更不能な設定オブジェクトです。`from_preset()`、`from_dict()`、`from_json()` で検証付き設定を作り、変更には `dataclasses.replace()`、保存には `config.to_json("scenario.json")` を使います。プリセット名は精度の保証ではありません。

## 保存と再利用

```python
from akutou_sky import load_dataset, save_dataset, export_results, plot_profiles

sky = load_dataset("results/improved/sky.nc")
save_dataset(sky, "results/copied.nc")
export_results(sky, "results/copied_viewer")
plot_profiles(sky, "results/copied_viewer/preview.png", azimuth=0)
```

`save_dataset()` は圧縮NetCDFを保存します。`export_results()` は `sky.nc`、再現条件の `manifest.json`、`colors.csv`、`spectra.csv`、`bands.csv`、`summary.json`、オフラインで動く `index.html` を生成します。既定では結果ファイルの上書きを拒否します。置き換える場合は `overwrite=True` を明示します。`plot_profiles()` は指定PNGを保存し、既存PNGも置き換えます。

公開API・返り値・単位は [APIリファレンス](docs/api.md)、実行例は [examples/](examples/) にあります。

## コマンドから使う

```sh
akutou-sky --help
akutou-sky simulate --preset preview --d 0 6 12 --elevations 0 5 30 90 \
  --output results/example --cache-dir model_cache --offline
akutou-sky export results/example/sky.nc --output results/example_viewer
akutou-sky plot results/example/sky.nc --output results/example/preview.png
akutou-sky validate --input results/example/sky.nc --checks spectral iterations \
  --cache-dir model_cache --offline
```

`python -m akutou_sky` も同じ入口です。設定の優先順位は「標準値 → プリセット → `--config` のJSON → 個別オプション」です。上書きには `--overwrite` が必要です。`validate` は入力と同じフォルダに `validation.json`、その `validation/` に比較計算を保存します。比較が基準外・未完了なら終了コード1、引数・入出力エラーなら2を返します。

`./sky.sh` は `akutou-sky simulate` と同じライブラリCLIを起動するショートカットです。旧Pythonスクリプトの互換入口と旧モデルのソースは削除しました。出力の再生成・作図・検証には上記の `export`・`plot`・`validate` を使います。

## キャッシュと副作用

`import akutou_sky` はダウンロード、保存先の作成、環境変数やMatplotlibの設定変更を行いません。計算時のキャッシュの優先順位は次のとおりです。

1. `cache_dir` / `--cache-dir`
2. `AKUTOU_SKY_CACHE_DIR`
3. 既存の `SASKTRAN2_DATABASE_ROOT`
4. OSのユーザーキャッシュ内の `akutou-sky/model-data`

初回の計算では不足するFASCODE・オゾン断面積・太陽スペクトルを取得します。`allow_download=False` / `--offline` なら不足をエラーにします。Mieデータは同じキャッシュ内に生成されます。複数プロセスの同時書込みは避け、計算は順番に実行してください。

`./sky.sh` はこのフォルダの `model_cache/` と `results/` を使います。CLIの既定出力は実行場所の `results/sky/` です。インストール先のパッケージ内には結果やキャッシュを作りません。

## 構成

```text
akutou_sky/
├── pyproject.toml             # インストール・依存・CLIの定義
├── src/akutou_sky/
│   ├── api.py / config.py     # 公開計算API・設定
│   ├── model.py / data.py     # 放射伝達・データキャッシュ
│   ├── colorimetry.py         # XYZ・xy・輝度・表示RGB・帯域平均
│   ├── io.py / provenance.py  # NetCDF・出力一式・再現条件
│   ├── export.py / plotting.py# CSV・HTML・任意依存のPNG
│   └── validation.py / cli.py # 数値感度の比較・コマンド
├── examples/ / docs/ / tests/
├── scenarios/                 # 大気条件のJSON例
├── results/ / model_cache/    # 既存結果とデータ
├── tools/verify_physics.py     # 現行モデルの物理チェック
└── sky.sh / setup_environment.sh
```

`results/improved/` と `results/overview/`、旧CSV・PNG・NPZは既存の計算結果です。元の詳しい説明は [モデル・既存結果の説明](docs/model-and-results.md) に保存しています。

## 開発と検証

```sh
python -m pip install -e '.[plot]'
python -m unittest discover -s tests -v
python -m pip wheel --no-deps . --wheel-dir dist
```

テストは積分・色変換・保存・公開APIと、保存済み小規模計算に対する実ソルバーの回帰を確認します。実計算の回帰は既存の `model_cache/` と `runs/v2_smoke.nc` を使い、ネットワーク取得は行いません。データのない配布先ではその検査をスキップします。物理モデルの収束不合格と、ライブラリ構成のテストの成功は別の判定です。

物理法則の診断を再実行する場合は `python tools/verify_physics.py` を使います。この診断は既知の偏光対称性の不合格を含み、全結果を記録したうえで終了コード1を返します。

伏角を細かく比較する新しいビューアーは **−6°〜＋6°、0.5°刻み（25点）** です。内部格子・前方散乱処理も改良しました。[数値感度の改善内容](docs/numerical-refinement.md)を参照してください。
