# 荷電二体系・解析近似基準解 v1

後続の匿名状態相互作用系とは独立な、leading-order Einstein–Maxwell adiabatic quasi-circular inspiral の解析近似 reference branch。

## ファイル

- `solve_charged_binary_analytic_reference_v1.py` — 閉形式 $t(r),\phi(r)$ から raw data を生成
- `charged_binary_analytic_reference_v1_raw.csv` — 40001 行の一次データ
- `charged_binary_analytic_reference_v1_metadata.json` — 条件・式・主要結果・SHA-256
- `validate_charged_binary_analytic_reference_v1.py` — 独立検算
- `validation_results.json` — 検算結果
- `plot_charged_binary_analytic_reference_v1.py` — raw data 専用図化
- `analytical_equations_and_initial_conditions_ja.md` — 解析式・導出・初期条件
- `results_analysis_ja.md` — 軌道・GW・EM 放射・1PN 診断の分析
- `REPRODUCIBILITY_ja.md` — 再現手順
- `figures/*.svg` — Google Drive 保存対象
- `figures/*.png` — ローカル目視確認用。Google Drive 保存対象外

## 位置づけ

一般 Einstein–Maxwell 二体問題の厳密解ではなく、また 1PN/2PN 精度の一般軌道でもない。

一般解析解のない後続の反復相互作用実験に対して、同一条件で比較可能な最初の独立解析基準を固定することが目的である。
