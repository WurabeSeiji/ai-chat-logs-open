# N=64・den=64（Δτ=2π/64）・2000 step 走行パッケージ（2026-09-07 作成）

背景: N=100 は1歩≈35秒（eigh 4950²）で2000歩≈17時間と判明し中断（実行ログ参照は N100 フォルダ）。
木原指示で N=64（M=2016、1歩≈2〜3秒、2000歩≈1.5時間見込み）・進捗ログ10歩単位に変更。

- 親: make_static_parent_N64_v1.py（正本生成器と diff 2行: parents_N64 / [64]）→
  parents_N64/parent_static_N00064_makeparent_20260905.npz
  監査（audit_parent_N64_v1.py → parents_N64/audit_parent_N64_v1.json）:
  残差 1.99e-13、Z4位相 1.4e-13 rad、CV 2.95%、σ=62.844、|Σz²|=1.8e-15、P_even=P_odd=0.5
- wrapper_run_N64_2000_v1.py — 物理正本（SHA照合）を宣言行8箇所置換
  （N=[64]・STEPS=2000・den=Nのみ・出力先・保存名 states_2000・metadata・**10歩ごと進捗print**）。物理式は無変更。
- plot_inflation_N64_v1.py / plot_complex_plane_final_N64_v1.py / plot3d_N64_2000_v1.py — 図化3本
  （N100版と同一様式。3DはフレームStride=4で501本・スライダー実step表示）
- run_all_N64_2000.sh — 一括再現

事前登録の予測: onset ≈ 9.5×64 ≈ 610 歩、増幅率 ≈ 0.5·2π/64/ln10 ≈ 0.021 log10/step、
飽和後は 1歩 2π/64=5.625° の剛体旋回・64歩周期再帰。
