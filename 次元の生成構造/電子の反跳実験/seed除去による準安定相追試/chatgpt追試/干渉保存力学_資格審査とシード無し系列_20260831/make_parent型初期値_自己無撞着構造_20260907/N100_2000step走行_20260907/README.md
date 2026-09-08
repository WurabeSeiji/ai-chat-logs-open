# N=100・den=100（Δτ=2π/100）・2000 step 走行パッケージ（2026-09-07 作成、未実行）

親: ../N100_親生成_20260907/parents_N100/parent_static_N00100_makeparent_20260905.npz
（make_parent 型・資格審査合格: 残差 9.8e-13、Z4 位相 8e-13 rad、P_even=P_odd=0.5）

- wrapper_run_N100_2000_v1.py — 物理正本（SHA照合）を宣言行7行のみ置換して exec
  （N=[100]・STEPS=2000・OFFSETS=(0,)・den124除去・出力先・保存名 states_2000・metadata N_range）。
  物理の式は無変更。所要見込み 約2.5〜5時間（eigh 4950² × 2000歩）。
- plot_inflation_N100_v1.py — インフレーション単図（semilogy、正本パネル様式）
- plot_complex_plane_final_N100_v1.py — step2000 複素平面単図（step0 図と同様式）
- plot3d_N100_2000_v1.py — 3D対話HTML（実験14テンプレート。フレームは4歩刻み501本、軌跡は毎歩保持、
  スライダーlabel実step化、x軸range 2000）
- run_all_N100_2000.sh — 一括再現（走行→図化3本）

予測（事前登録）: onset ≈ 9.5×100 ≈ 950 歩（時計約9.5回転）、増幅率 ≈ 0.5·2π/100/ln10 ≈ 0.014 log10/step 級、
飽和後は 1歩 2π/100=3.6° の剛体旋回・100歩周期再帰（N=40 の実測則の外挿）。
