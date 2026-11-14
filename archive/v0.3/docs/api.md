# 公開API

`from akutou_sky import ...` で読み込む入口です。

| API | 用途 | 返り値 |
|---|---|---|
| `ModelConfig` | 物理条件・数値設定。変更不能なdataclass | 設定 |
| `simulate(depressions, *, ...)` | 放射伝達と色変換 | `xarray.Dataset` |
| `compute_radiance(depressions, *, ...)` | 太陽成分のI/Q/U | `Dataset` |
| `add_colorimetry(dataset, background=None)` | XYZ等を追加 | 新しい `Dataset` |
| `xyz_from_radiance(radiance, wavelengths, bounds=None)` | 配列からXYZ | 最後の軸がX/Y/Zの配列 |
| `display_values(xyz)` | 色度と表示色 | `(xy, preview_sRGB, gamut_clipped)` |
| `band_mean(elevations, xyz, lower, upper)` | 高度角帯の立体角重み付きXYZ平均 | XYZの配列 |
| `save_dataset(dataset, path, *, overwrite=False)` | 圧縮NetCDFを保存 | `Path` |
| `load_dataset(path)` | ファイルを閉じてメモリに読む | `Dataset` |
| `export_results(dataset, directory, *, overwrite=False)` | NetCDF・CSV・HTML・再現条件を保存 | 出力先 `Path` |
| `plot_profiles(source, out, *, azimuth=0)` | DatasetまたはNetCDFからPNG | PNGの `Path` |
| `compare_runs(reference, candidate)` | 同じ物理条件・同じ角度の比較 | `dict` |
| `validate_sensitivity(input, *, ...)` | 格子を変えて数値感度を実計算 | `dict` |
| `build_explorer(directory, *, cache_dir=None, allow_download=True, resume=False, progress=None)` | 霧・塵の64条件とオフライン比較HTMLを作成 | 出力先 `Path` |

## 計算引数

`simulate()` と `compute_radiance()` の共通キーワードは次のとおりです。

- `config=None`: `ModelConfig`。未指定は通常設定。
- `elevations=None`: 高度角の配列。既定は0〜5°未満を0.25°刻み、5〜90°を1°刻み。
- `azimuths=None`: 太陽からの相対方位の配列。既定は `[0]`。
- `cache_dir=None`: 入力とMieデータの保存先。
- `allow_download=True`: 不足する入力の取得を許すか。
- `progress=None`: 文字列1つを受け取るコールバック。既定はログなし。

`simulate()` にだけ `background=None` があります。観測者位置での一様・無偏光背景のCSVで、列名は `wavelength_nm,radiance_W_m2_sr_nm`、380〜780 nmを完全に覆う必要があります。背景未指定では太陽散乱光だけです。

角度は重複なしの有限の1次元配列、単位は幾何学的な度です。伏角は−18〜24°、高度角は0〜90°、方位は−180〜180°。伏角が正なら太陽中心は地平線下、方位0°は太陽側です。

## Datasetの主な変数

主な座標は `depression, azimuth, elevation, wavelength, stokes`。波長は区間中心、端は `wavelength_bounds` です。

| 変数 | 意味・単位 |
|---|---|
| `solar_stokes_radiance` | 太陽成分I/Q/U。W m⁻² sr⁻¹ nm⁻¹ |
| `solar_irradiance` | 区間平均太陽スペクトル。W m⁻² nm⁻¹ |
| `total_radiance` | 太陽成分I＋背景。W m⁻² sr⁻¹ nm⁻¹ |
| `solar_XYZ`, `total_XYZ` | CIE 1931 2°。Yはcd m⁻²、X/Zは同じ683係数で積分 |
| `solar_xy`, `total_xy` | 色度。光量ゼロではNaN |
| `*_preview_sRGB` | 視線ごとの最大線形RGBを1に規格化した表示色。輝度情報なし |
| `*_gamut_clipped` | sRGBの負成分をクリップしたか |
| `background_radiance` | 観測者位置での背景スペクトル |
| `pressure`, `temperature`, `ozone_number_density` | 鉛直物性。Pa、K、m⁻³ |
| `aerosol_extinction_550`, `fog_extinction_550` | 微粒子・霧の550 nm消散係数。m⁻¹ |
| `solar_degree_linear_polarization` | 偏光計算時のDoLP。I=0ではNaN |

`compute_radiance()` には色変換後の変数がありません。`export_results()` は必要なら色変換を加えて出力します。渡したDatasetの値は変更しません。

## 設定・検証・例外

プリセットは `preview, standard, overview, polarized, fog_preview, fog_refined` です。未知のキーや非物理的な設定は拒否します。`num_stokes=3` では既知の対称性不合格を `ExperimentalPolarizationWarning` で通知します。

`relative_humidity`（0〜0.95の割合）と `hygroscopicity_kappa`（0〜1.5）を指定すると、Mie微粒子の吸湿成長を有効にします。この場合 `aod550` は乾燥状態のAODです。属性 `aerosol_aod550_ambient`、`aerosol_growth_factor`、`aerosol_wet_median_radius_nm` に湿った状態を記録します。標準値はRH=0・κ=0なので従来結果を維持します。

`fog_aod550 > 0` は独立した霧の水滴層です。`fog_scale_height_m`、`fog_median_radius_um`、`fog_mode_width` で高さ・粒径を指定できます。細かい `boundary_step_m`・`source_boundary_step_m`、256以上の `legendre_moments` と `delta_m_scaling=True` が必要です。`fog_preview` と `fog_refined` にこれらの設定を用意しています。`boundary_step_m=0` は従来の低層格子を使用し、正なら0〜1 kmだけ細かくします。[物理的な仮定と制限](fog-and-dust.md) を確認してください。

`build_explorer()` は64条件・伏角−6°〜＋6°の0.5°刻み（25点）を保存します。塵が0のケースは湿度・乾燥粒子の設定に依存しないため、同じ霧量の解を再利用します（52種類の異なる大気）。再利用元はNetCDFとmanifestに記録します。`resume=True` の場合は同じコード・設定・入力・座標のキャッシュを確認して再利用します。出力には比較ページと、各条件のNetCDF・再現記録が含まれます。計算過程を見たい場合は `progress=print` を指定します。

`validate_sensitivity(input, *, checks=(...), output_dir=None, run_directory=None, cache_dir=None, allow_download=True, progress=None, max_phase_terms=None)` はDatasetまたはNetCDFのパスを受け取ります。既定検査は `spectral, vertical, horizontal, angular, iterations`、追加は `moments, delta_m, single_scatter`。`output_dir=None` なら結果ファイルを書かず、レポートを返します。保存先を与えると `validation.json` と `runs/` を作ります。`run_directory` で比較計算の保存先を変更できます。`max_phase_terms` は大気層数×散乱源層数×水平列数×角度点数×Stokes数×初期散乱の係数数を制限し、大きな位相関数テーブルの生成を防ぎます。RAM容量そのものではありません。指定すると鉛直・水平の細分化量を調整し、それでも上限を超える検査は未実施として記録します。

`all_requested_checks_complete` は検査が実施できたか、`passes_all_sampled_checks` は相対Y差5%・xy距離0.005の比較目安を満たしたかです。物理的な誤差保証や実測校正ではありません。Dataset入力にはファイルハッシュがないため、HTMLの検証表示には自動で紐付けられません。表示まで必要なら先にNetCDFを保存し、そのファイルを検証します。

不正な設定・角度・比較は `ValueError`、上書き拒否は `FileExistsError`、オフラインでの入力不足は `FileNotFoundError`、取得・入出力失敗は `OSError` 等を返します。v2のDatasetの座標・変数名・単位は維持しています。旧スクリプトの互換入口は削除し、公開Python APIとライブラリCLIに統一しました。

`delta_m_method="split"` は整合した打切り多重散乱と一回散乱補正を組み合わせます（スカラー、Delta-M有効時のみ）。`source_surface_m=1` は地表付近の散乱源点を追加し、`horizontal_sampling="near_horizon"` は±8°を重点的に細分化します。`horizon_quadrature=True` はバックエンドの地平線に合わせた角度積分、FalseはLebedev積分です。旧設定の既定値は維持しています。新しい設定と検証範囲は[数値感度の改善](numerical-refinement.md)を参照してください。

`single_scatter_step_m=500` は観測視線の一回散乱を独立した高度格子で計算します。0〜1 kmは設定値の1/10、1〜10 kmは1/2、それより上は設定値の刻みです。霧がある場合、最下層は霧のスケール高度の1/4以下にします。0は従来の格子を使います。正の値は `delta_m_method="split"` の場合だけ使えます。結果の `solar_single_scatter_radiance` と `single_scatter_altitude` に成分と格子を保存します。
