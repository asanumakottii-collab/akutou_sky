# 数値手法を選ぶための診断計算

このフォルダは改良中の比較記録です。利用する比較ページは [`../atmosphere_explorer/index.html`](../atmosphere_explorer/index.html) です。

`final-*.nc` は最終設定を選んだ際の少数視線・40 nm幅の試算です。公開ビューアーの20 nm幅・25伏角の計算とは異なります。`final-angular.json` は110→146角度点の差、`cache_isolation.json` は高湿度→粒子なし→高湿度の順に計算したときの再現性を記録しています。

`independent-single-*.nc` と `single-grid-*.nc` は一回散乱の格子依存性を調べた記録です。`independent-single-full-regression.nc` は最終設定で一から全放射伝達を解き直した90視線の検算で、成分を置換して更新した公開データとの相対差は最大2.37×10⁻¹⁰でした。[検算の記録](../atmosphere_explorer/component_replacement_regression.json)に比較範囲とファイルのSHA-256があります。

`source_before_single_grid/` は成分置換前のソースを保存したものです。対応する中間計算は `../atmosphere_explorer/history/0.3.0-coarse-single/` に保存しています。

それ以外は格子、Delta-M次数、角度積分法を比較した中間試算です。採用しなかった設定も含み、最終版の精度を示すものではありません。最終出力の検査は [`../atmosphere_explorer/validation_summary.md`](../atmosphere_explorer/validation_summary.md) を参照してください。
