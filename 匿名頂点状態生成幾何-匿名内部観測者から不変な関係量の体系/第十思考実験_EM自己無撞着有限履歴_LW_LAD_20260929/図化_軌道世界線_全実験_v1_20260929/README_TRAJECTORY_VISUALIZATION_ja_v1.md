# 第十思考実験系列：軌道・世界線・波束軌跡 図化 v1

過去の Paper 6 / Paper 7 の軌道図規約（明示的な読出し座標、空間図は同一x/yスケール、開始/終了点、SVG+PNG同時出力）を継承した再現パッケージ。

## 図の意味
- figure01: retarded EM/LAD 参照実装の粒子位置。
- figure02: exact anchor / strict spin の複素位相軌道（空間軌道ではない）。
- figure03: KT 二中心の物理世界線 `x_phys=a(t)x_coord`。
- figure04: RWZ / RN coupled の累積放射読出し履歴。
- figure05: strict RN の波束強度重心。
- figure06: GW/EM 同時バランス最適初期状態の波束強度重心。
- figure07: 同じ dual-balance run の時空強度図。

波束重心は `r_c = sum_i r_i |psi_i|^2 / sum_i |psi_i|^2`。readout は transition へ一切フィードバックしない。

## 再現
`python plot_experiment_trajectory_suite_v1.py`

依存: numpy, scipy, matplotlib。元生成コード、入力JSON、過去のPaper 6/7図化コードも同梱。
