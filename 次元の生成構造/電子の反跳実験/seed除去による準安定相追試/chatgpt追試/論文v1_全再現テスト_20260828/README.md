# 論文v1_全再現テスト_20260828

論文 v1（DOI 10.5281/zenodo.22112009）と note 記事 n07c3e4c97e3a が引用する 15 パッケージ（Zenodo 14 zip ＋ decompactification zip）を `original/`（比較基準）と `rerun/` に展開し、プログラムを**パス変更 4 箇所以外は無変更**で再実行して保存データ・図と突合した。結果は `再テスト結果_論文v1引用プログラム全再実行_20260828.md`。

- 主要数値主張（Floquet 乗数 9 桁、onset–残差則、成長率、厳密保存、等分配、閉包探索 exact/MITM）は再現。
- 軌道後半の量（秩序化 step、H∥/H⊥ 最終値、モジュライ位相、出発深さの桁数）は丸め誤差増幅で計算機依存。
- N=5 の 3+3+2+2 時間分離（記事図 7、2627/4923 step）等 101 ファイルは生成プログラム未同梱で再テスト不能。

## 再現

```bash
cd 論文v1_全再現テスト_20260828 && bash run_all.sh   # zips 展開 → 再実行（約 2 分）→ 突合
```
必要: python3 + numpy/pandas/matplotlib/scipy/Pillow、c++（clang 可、`compat/` で `bits/stdc++.h` を補う）。zip は `../*.zip` と `../../../../ゼロ閉塞の幾何・代数構造/…/A2a_N5_ab_probe_20260825/N5_gamma_continuum_test_bundle_20260825.zip` から `zips/` に集める（`run_all.sh` が行う）。
