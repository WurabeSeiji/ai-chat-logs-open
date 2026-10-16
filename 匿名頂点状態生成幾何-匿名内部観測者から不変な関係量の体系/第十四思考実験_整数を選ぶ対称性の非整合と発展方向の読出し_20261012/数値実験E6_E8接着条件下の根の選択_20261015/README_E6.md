# 第十四思考実験 数値実験 E6（2026-10-15）：E₈ の接着条件の下での根の選択

## 目的
有限位数根 R_{n,m} ＝ cos²(πm/n) は稠密にあり、閉包だけでは 23/124 を隣の根から選べない。E4 は四次元フルヴィッツ格子が「分母が 4 の倍数」の類を選ぶことを示したが、その類の中で 23/124 を選んでいたのはセルの数 137（半径 3 の球）だった。本実験は、第 24 章の E₈ の接着条件（二つのフルヴィッツ四元数 x, y で x − y ∈ (1＋i)H。この対の集合が E₈ 格子、verify_e8.py で 240 本確認）を力学に毎ステップ課したとき、137 を入力せずに、mod 4 の類の中で 124 が 92・104・116・136・148・120 より選ばれるかを調べる。著者判断（10-15）：接着条件を初期値だけに課すのは不自然なので、毎ステップ課す。

## 方法（E4 part A の忠実な複製＋最小の変更）
- 元プログラム：E2 `te14_E2_hurwitz_124_23.py`、E3 `te14_E3_ground_and_levels.py`、E4 `te14_E4_selection_pullin_recovery.py`・`te14_E4_mod4_check.py` を `control/` に複製（SHA256 一致、`control/SHA256SUMS_copied_sources.txt`）。
- 対照走行：`control/` で E4 をそのまま実行。`te14_E4_summary.json`・`te14_E4_mod4_summary.json` は保存済みとバイト一致、`te14_E4_rows.csv` は改行コード以外一致（`control/SHA256SUMS_control_outputs.txt`、所要 9.3 s＋1.6 s）。
- 変更点（`te14_E6_e8_glued_selection.py`）：状態を x ∈ H から対 (x, y) に替える。回転 u ＝ cos φ ＋ i sin φ（軸 Z4）を両方に左から掛け、丸めを「接着格子への最近点」にする。z ∈ (1＋i)H ⟺ N(z) が偶数なので、接着条件は 4N(x−y) ≡ 0 (mod 8) で判定。丸めは E2 の round_h を各成分に当て、接着条件を満たさないときだけ近傍（整数点 3⁴＋半整数点 3⁴）から距離の和が最小の接着する対を取る。読み出しは E3 と同じ巻き数 W/P を x と y について計算し、x の W/P が target に一致すればロック。
- 三つのモード：`same`（y ＝ x、接着なし。乱数列も E4 と同じで、x の軌道と part A の数が E4 と一致しなければならない内部対照）、`free`（y は独立な別の点、接着なし）、`glued`（接着あり）。
- 対象：mod 4 の類 {17/92, 19/104, 21/116, 23/124, 25/136, 27/148, 23/120}、gray 深さの図の窓の隣 {23/122, 23/123, 23/125, 24/127, 22/119}、E4 の対照 {23/126, 24/131, 18/97}。半径 32・64・96、各 30 点。乱数種 20261015（same は 20261014）。MAXIT 3,000,000。観測値（137、α）は入力しない。

## 内部対照の結果
`same` モードは E4 part A の 21 行（7 対象 × 3 半径）とロック数・母数ともすべて一致（`te14_E6_rows_same.csv` と `control/te14_E4_rows.csv` の part A）。所要 4.0 s。

## 結果
（走行後に追記）

## 再現
```
cd control && python te14_E4_selection_pullin_recovery.py && python te14_E4_mod4_check.py   # 対照走行
cd .. && python te14_E6_e8_glued_selection.py same    # 内部対照（E4 part A と一致すること）
python te14_E6_e8_glued_selection.py free             # y 独立・接着なし
python te14_E6_e8_glued_selection.py glued            # 接着あり
```
出力：`te14_E6_rows_<mode>.csv`（mode, target, R, lock_x, lock_y, n, sec）、`te14_E6_summary_<mode>.json`。
