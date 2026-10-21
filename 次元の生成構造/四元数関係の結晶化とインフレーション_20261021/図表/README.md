# 図版の再現確認と英語版の図（2026-10-21）

対象：`../quaternion_crystallization_inflation_ja_v1.md`（v1.1 も同じ図）の図 1〜9。

## 再現確認（日本語の図）

- 元のプログラム `../make_figures_quaternion_crystallization_20261021.py` をそのまま実行して `ja/` に再生成し、フォルダ直下の `fig01〜fig09_*.svg` とバイト単位で比較した：**9 枚すべて一致**（`ja_compare.txt`）。
- 図 5 の確認出力：`fig05: hull->Delaunay match: True 53 53`（`ja_run_log.txt`）。付録 A の記載どおり。
- PNG（cairosvg、幅 2 倍）も `ja/` に生成した。

## 英語版の図

- 図版生成プログラムに英語出力を加えた版 `../make_figures_quaternion_crystallization_20261021_v1.1.py` を作った（元のプログラムは上書きしていない）。
  - 実行：`python make_figures_quaternion_crystallization_20261021_v1.1.py [出力フォルダ] [ja|en]`
  - 図の文字だけを `L(日本語, 英語)` で切り替える。数値・座標・乱数・図形は v1 と同一。
  - **v1.1 の ja 出力は v1 とバイト単位で一致**（9 枚、`ja_v1.1_compare.txt`）。
- 英語版の図は `en/fig0N_*_en.{svg,png}`（9 枚）。図 5 は英語版でも `hull->Delaunay match: True 53 53`（`en_run_log.txt`）。
- 英語版の SVG に日本語の文字（全角の空白・読点を含む）が残っていないことを機械的に確認した。
- 目視で直した点：図 4 の左下の注記が右の枠に重ならないよう略記にした（内容は同じ）。図 8 の式の読点「、」を「, 」に、図 6 の全角空白を半角にした（英語版のみ）。

## 再現の手順

```
図表/run_all.sh
```

日本語の再生成と比較、v1.1 の ja 出力の比較、英語版の生成、`SHA256SUMS.txt` の更新を行う。依存：numpy、scipy（ConvexHull、Delaunay）、PNG 変換に cairosvg。

**環境メモ（このマシン）：** システムの Python では Homebrew の cairo が見えないため、`~/lib/libcairo.2.dylib` を `/opt/homebrew/lib/libcairo.2.dylib` へのシンボリックリンクとして置いた（2026-10-21）。cairosvg がない環境では SVG だけが生成される。
