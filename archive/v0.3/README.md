# akutou-sky — 朝焼け・夕焼け・薄明のモデリングライブラリ

更新日：2026-10-09。ライブラリ **0.3.0**、物理モデル **2.2**、SASKTRAN2 **2026.10.0** を使用します。配布名は `akutou-sky`、Pythonのパッケージ名とこのフォルダ名は `akutou_sky` です。

球面大気の多重散乱から分光放射輝度、CIE XYZ・xy・輝度、表示用sRGBを計算するPythonライブラリです。Rayleigh散乱、MieまたはHGエアロゾル、オゾン吸収、太陽入射光の屈折を扱います。モデル2.1で微粒子の吸湿成長と独立した霧の水滴層を追加し、2.2で前方散乱処理、角度・空間分解能、一回散乱の独立した高度格子を改善しました。

**物理モデルは数値収束未達の研究用試算です。** これは、数値格子を細かくした際の差が、設定した比較基準まで十分に小さくなっていない条件があることを意味します。最新の代表条件の検査では最大10.394%の輝度差が残り、観測との校正も未実施です。最新の比較ページはスカラー計算で、別設定の偏光計算には既存の粗い格子の左右対称性の不合格も残っています。[改善内容と検証](改善内容と検証.md) と [科学的根拠と利用上の注意](科学的根拠と利用上の注意.md) に、根拠・結果・利用範囲を記載しています。

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

塵は乾燥状態の550 nm AOD、霧は水滴だけの550 nm AODです。AODは無次元の鉛直光学的厚さで、質量濃度・視程を直接指定するものではありません。相対湿度は0〜0.95の割合で与え、RHと `hygroscopicity_kappa` がともに正のとき吸湿成長を計算します。RHと霧の量は独立した入力です。

`explore` は `fog_refined` を使います。角度110点・水平41列・次数16のDelta-M近似と一回散乱の独立格子（低層50 m、中層250 m、上層500 m刻み）を採用しています。Python APIや `simulate` で同じ改善を使う場合も、`fog_refined` を明示してください。`ModelConfig()` / `standard`、軽量の `preview`、従来の `fog_preview` は別の数値設定です。プリセット名は精度の保証ではありません。[霧と塵の詳しい設定](docs/fog-and-dust.md) と [計算方法](docs/numerical-refinement.md) も参照できます。

## 保存済みのresultsの選び方

| フォルダ | 内容・用途 |
|---|---|
| [atmosphere_explorer](results/atmosphere_explorer/index.html) | 最新の霧・塵・湿度の比較。全64条件、伏角25点、一回散乱の独立格子を反映 |
| [atmosphere_explorer 2](<results/atmosphere_explorer 2/index.html>) | 一回散乱の独立格子を導入する前の中間結果。表示上の版番号は最新と同じだが、計算内容とソースのハッシュが異なる |
| [atmosphere_explorer/history/0.2.0](results/atmosphere_explorer/history/0.2.0/index.html) | 霧・塵の比較ページの初期版 |
| [improved](results/improved/index.html) | 以前の太陽方向の試算。5 nm幅・38角度点で、最新の霧・吸湿の出力ではない |
| [overview](results/overview/index.html) | 以前の3方位の試算。5 nm幅・26角度点。improvedとは方位以外の設定も異なる |
| [smoke](results/smoke/index.html) | 以前の動作確認用の小規模出力 |
| [numerical_refinement](results/numerical_refinement/README.md) | 数値手法を選ぶための診断計算・中間記録 |

通常の霧・塵の比較には `atmosphere_explorer` を使用します。メインフォルダ直下の日本語CSV・PNG、`summary.json`、`provenance.json` は初期計算の成果物です。出力ごとの設定は、そのフォルダのmanifestまたはNetCDFで確認してください。

## 最新の検証状況

保存済みの直近のライブラリテストは **48件すべて通過**しました。全64条件について、合計分光放射輝度の有限・非負、XYZの再計算、AODの積分、HTML用データ、ソースのハッシュの整合性を確認しています。一回散乱成分を更新したデータと、一から解いた全計算も、高湿度の90視線で一致しました。

数値感度検査は **11件中10件完了、5件基準内、5件基準未達、1件未実施**です。比較基準は最大輝度差5%以内かつ色度xy距離0.005以内。高湿度の鉛直格子差は10.394%、一回散乱格子差は輝度3.352%・xy距離0.005191で、どちらも基準未達です。Delta-M次数16→32は検証のメモリ負荷上限を超えるため未実施です。

これらは代表条件の数値設定間の差で、実空に対する誤差の上限や全64条件の収束証明ではありません。結果は [検証の集計](results/atmosphere_explorer/validation_summary.md)、[全記録JSON](results/atmosphere_explorer/numerical_sensitivity.json)、[ライブラリテストのログ](library_test_results.txt) にあります。

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

config = ModelConfig.from_preset(
    "fog_refined", aod550=0.08, fog_aod550=0.12,
    relative_humidity=0.85, hygroscopicity_kappa=0.30,
)
sky = simulate(
    depressions=[-6, -3, 0, 3, 6],
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

sky = load_dataset("results/atmosphere_explorer/humid/sky.nc")
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
│   ├── microphysics.py       # 吸湿成長・微粒子・霧
│   ├── delta_m.py            # 前方散乱の打切り・一回散乱補正
│   ├── colorimetry.py         # XYZ・xy・輝度・表示RGB・帯域平均
│   ├── io.py / provenance.py  # NetCDF・出力一式・再現条件
│   ├── export.py / plotting.py# CSV・HTML・任意依存のPNG
│   ├── explorer.py / templates/ # 4種類の霧・塵の比較ページ
│   └── validation.py / cli.py # 数値感度の比較・コマンド
├── examples/ / docs/ / tests/
├── scenarios/                 # 大気条件のJSON例
├── results/ / model_cache/    # 既存結果とデータ
├── tools/                    # 全64条件の検査・数値感度・物理診断
└── sky.sh / setup_environment.sh
```

以前の `improved`・`overview` 等の計算方法は [モデル・既存結果の説明](docs/model-and-results.md) に保存しています。この資料は旧設定を扱う説明で、現在の基準は本READMEと更新した日本語2文書です。

## 開発と検証

```sh
python -m pip install -e '.[plot]'
python -m unittest discover -s tests -v
python -m pip wheel --no-deps . --wheel-dir dist
```

テストは積分・色変換・保存・公開APIと、保存済み小規模計算に対する実ソルバーの回帰を確認します。実計算の回帰は既存の `model_cache/` と `runs/v2_smoke.nc` を使い、ネットワーク取得は行いません。データのない配布先ではその検査をスキップします。物理モデルの収束不合格と、ライブラリ構成のテストの成功は別の判定です。

最新の比較ページを検査する場合は、次を実行します。保存データの検査、代表条件の感度検査、検証文書とページ内表示の更新を行います。計算は順番に実行してください。

```sh
python tools/validate_explorer.py --directory results/atmosphere_explorer --cache-dir model_cache
```

検証ツールの既定の作業量上限は `--max-phase-terms 520000000` です。これはRAMのバイト数ではありません。保存されている今回の結果では未実施が1件あるため、このツールは終了コード1を返します。個々の基準未達はJSONと集計文書に記録します。

物理法則の診断を再実行する場合は `python tools/verify_physics.py` を使います。この診断は既知の偏光対称性の不合格を含み、全結果を記録したうえで終了コード1を返します。
