# Stage A2c N=5 方向基底の系譜追跡 報告書

## 実行状態

**STAGE_A2C_COMPLETE**

記述的分類: **ROTATING_OR_MIXED_LINEAGE**

この分類は方向系譜だけを対象とする。H1/H2/H0、三方向成立step、物理的追加次元の存在は判定していない。

## 使用原本とSHA-256

| file | sha256 | absolute_path |
|---|---|---|
| run_n_scaling_lowrank_v1.py | ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py |
| run_plane_flow_exact_v1.py | 9cf28ca8c0d2ad8fac2f0f6dae045248695247c5809c21ccb2069ef91a94ab76 | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_exact_v1.py |
| run_n300_dimension_saturation_v2.py | 229938a66631057426f187ed80b17de08cfcb9107cfe509c30f5bbdcca3a03e6 | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/code/run_n300_dimension_saturation_v2.py |
| run_paper7_5color_timeseries.py | fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503 | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py |
| run_paper7_transverse.py | ac1073bea329971de3ff4c2fd1588d926029a8502c21e8cc01f406acb86ad60b | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_transverse.py |
| run_plane_flow_approx_v1.py | a9d247a8070d849fe989e35e00320e470968354614424c9bdceda57132d9f0fa | /Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/run_plane_flow_approx_v1.py |

固定原本、依存原本、Stage A0成果物は実行前にSHA-256照合した。原本コードと既存成果物は編集・上書きしていない。

## Python環境

- Python: `3.9.6 (default, Jan  9 2026, 11:03:41)  [Clang 17.0.0 (clang-1700.6.4.2)]`
- NumPy: `2.0.2`
- OS: `macOS-26.3.1-arm64-arm-64bit`

## 軌道再現

- 対象: `N=5`, `float64`, Stage A0の `delta=1e-15`, `seed_index=0`
- 再実行範囲: step `0..5000`
- fのStage A0とのbitwise一致: `True`（`5001`点）
- f最大絶対誤差: `0.00000000000000000e+00`
- q1〜q4、rank_q、gram_rank、dominant eigenvalueのStage A0とのbitwise一致: `True`（`671`保存点）
- q最大絶対誤差: `0.00000000000000000e+00`
- existing crossing: `1167`
- 親残差: `2.13989836008277647e-13`

したがって、以下は新しい軌道探索ではなく、Stage A0軌道内部で既存関数が生成する基底の抽出である。

## D34(t)の既存定義

各stepで既存 `gram_reduce` と `dominant_plane` により `Bdom(t)` を得て、既存

`D34(t) = s4_new_dirs(B0, Bdom(t))`

をそのまま使用した。比較の正本は列ではなく

`P34(t) = D34(t) D34(t)^T`

である。個別direction 3/4表示は、既存5色コードと同様に固定other空間へ射影・QR後、既存 `align_2d` で連続化した。

`S4(t)` は既存 `run_paper7_transverse.py:s4_basis` で再構成し、`PS4(t)=S4(t)S4(t)^T` を保存した。

## D34_lateの比較用構成

step `1800..2500` の全 `701` 個の `P34(t)` を平均し、平均射影行列の上位2固有ベクトルを `D34_late` とした。これは比較用代表部分空間であり、新しい物理方向定義ではない。

- 平均射影行列固有値: `[0.9412859657007462, 0.8882112083926906, 0.08336492160628962, 0.04300645195943384, 0.027927932922577742, 0.01100822855346245, 0.004572243751708384, 0.0006230471130909696, 1.6331388345634138e-16, -1.7166477719221818e-18]`
- late window内のD34対D34_late overlap中央値/最小値: `9.464955e-01` / `5.557283e-01`
- late window内の最大主角中央値/最大値: `2.726340e-01` / `9.545436e-01` rad

## 横摂動方向の既存コードからの再構成

既存 `run_paper7_transverse.py` の固定値をそのまま用いた。

- `t0 = crossing + 3000 = 4167`
- 方向PRNG seed: `70005`
- 既存方向seed数: `3`
- epsilon: `[1e-08, 1e-10, 1e-12, 1e-14]`
- 規則: `eta.real` と `eta.imag` をそれぞれ既存 `S4(t0)^perp` へ射影し、複素ベクトル全体を正規化
- 保存正本: 原本が使った複素 `eta`
- 2次元比較表現: `span(eta.real, eta.imag)` のQR基底
- 任意方向の追加: `False`
- Stage A0の最初の保存step 4200における全12レコード文字列一致: `True`
- 基準軌道状態のbitwise一致: `True`

同じt0のS4とTperpの最大直交誤差は `2.766343e-16` であり、構成上の直交を数値誤差内で確認した。

## 急拡大前のD34の数値解像度

crossing前step 0〜1166のq解像度帯別点数:

- `Q0_LT_1E-8`: 697 step
- `Q1_1E-8_TO_1E-6`: 286 step
- `Q2_1E-6_TO_1E-4`: 93 step
- `Q3_1E-4_TO_1E-2`: 91 step
- `Q4_GE_1E-2`: 0 step

- crossing前で `min(q3,q4)/q1 >= 1e-6` の比較点数: `184`
- その比較点でのD34対D34_late overlap中央値: `0.149478938071618`
- その比較点でのD34対Tperp最大overlap中央値: `0.38521163690339677`

`min(q3,q4)/q1 < 1e-8` の早期基底も削除していないが、数値解像度不足帯として分離した。この帯の列方向を物理的に確定した方向とは扱わない。

## 急拡大中のD34部分空間回転

固定step 900〜1400について:

- 連続step間最大主角の最大値: `3.721679e-01` rad
- 連続step間projector distance最大値: `3.636360e-01`
- 最大主角の区間内総和（経路長の記述量）: `3.871064e+00` rad
- basis column swap: `0` 回
- sign flip: `4` 列

総和はイベント閾値ではなく、固定観察窓内の連続step回転量の単純和である。

## 急拡大後D34との主角・overlap

分類評価は `min(q3,q4)/q1 >= 1e-6` を初めて満たすstep `983` からstep 5000までとした。

- D34(t)対D34_late overlap中央値: `7.028877e-01`
- 最大主角中央値: `6.697581e-01` rad

## 横摂動方向との主角・overlap

評価範囲でのD34(t)対Tperp overlap中央値:

seed 0: `5.436309e-02`, seed 1: `6.638024e-02`, seed 2: `6.074581e-02`

固定集約値（3 seed中央値の最大）: `6.638024e-02`

D34_late対各Tperp:

| seed | epsilons_sharing_this_direction | overlap | theta_1_rad | theta_2_rad | projector_distance | same_t0_S4_vs_Tperp_orthogonality_error |
|---|---|---|---|---|---|---|
| 0 | 1e-08;1e-10;1e-12;1e-14 | 3.29554650868904009e-02 | 1.33824918237557799e+00 | 1.45741394486574150e+00 | 1.39071530869053794e+00 | 1.56213833072082129e-16 |
| 1 | 1e-08;1e-10;1e-12;1e-14 | 3.01241105974697568e-02 | 1.32831027632905307e+00 | 1.51986066247149143e+00 | 1.39274971865194019e+00 | 2.76634302133747398e-16 |
| 2 | 1e-08;1e-10;1e-12;1e-14 | 1.81712190614377461e-02 | 1.41178998538540212e+00 | 1.46442721019391509e+00 | 1.40130566325735062e+00 | 1.81817830009204883e-16 |

epsilonは同一seedの方向を変えないため、方向部分空間比較では3個のunique seed方向を用い、12個のepsilon条件は再構成検証表に保持した。

## 個別direction列の交換と部分空間系譜

`consecutive_subspace_continuity` は、各連続stepについて以下を別々に保存した。

- raw列割当からのcolumn swap
- 対応列のsign flip
- signed permutationから残るD34内部回転
- `align_2d` 後の個別列内積
- P34主角・overlap・projector distanceによるambient部分空間回転

したがって、direction 3/4の色交換や符号反転だけを、物理的2次元部分空間の選び直しとは解釈していない。

## f水準別の方向系譜

| level_label | status | step | q_resolution_band | D34_vs_late_overlap | D34_vs_late_max_angle_rad | D34_vs_Tperp_seed0_overlap | D34_vs_Tperp_seed1_overlap | D34_vs_Tperp_seed2_overlap |
|---|---|---|---|---|---|---|---|---|
| 1e-20 | found | 289 | Q0_LT_1E-8 | 2.17346408090011101e-01 | 1.45121555299987737e+00 | 1.24956985743780669e-01 | 2.51823041730216246e-01 | 2.67487557602620774e-01 |
| 1e-18 | found | 382 | Q0_LT_1E-8 | 2.18169112706752377e-01 | 1.44683858246647445e+00 | 1.26190357071148868e-01 | 2.50751531463366928e-01 | 2.63156746137974129e-01 |
| 1e-16 | found | 476 | Q0_LT_1E-8 | 2.17968282863131024e-01 | 1.44687284512004588e+00 | 1.26371344805239288e-01 | 2.50969732389757061e-01 | 2.63003941433167965e-01 |
| 1e-14 | found | 569 | Q0_LT_1E-8 | 2.17833387434074155e-01 | 1.44702185954333418e+00 | 1.26474470655342258e-01 | 2.50982737014402268e-01 | 2.62873031793794176e-01 |
| 1e-12 | found | 662 | Q1_1E-8_TO_1E-6 | 2.16400017943768730e-01 | 1.44493015691649829e+00 | 1.29068360089831885e-01 | 2.50696860126402521e-01 | 2.59539125752811239e-01 |
| 1e-10 | found | 755 | Q0_LT_1E-8 | 1.89688724982191081e-01 | 1.41019435467696708e+00 | 1.97280010321965127e-01 | 2.64415856417339057e-01 | 2.12136221142953801e-01 |
| 1e-08 | found | 849 | Q1_1E-8_TO_1E-6 | 1.53485539623484046e-01 | 1.36564696479516612e+00 | 3.68643083945642158e-01 | 3.53652508821606537e-01 | 1.95515292411112235e-01 |
| 1e-06 | found | 942 | Q1_1E-8_TO_1E-6 | 1.51462325949341725e-01 | 1.36384054211152472e+00 | 3.84471814054445304e-01 | 3.64423014385491517e-01 | 1.98545439342969310e-01 |
| 1e-04 | found | 1035 | Q2_1E-6_TO_1E-4 | 1.50734480135932297e-01 | 1.36338930358272115e+00 | 3.85351006862287626e-01 | 3.65827046775420217e-01 | 1.98538877156239163e-01 |
| 1e-03 | found | 1082 | Q3_1E-4_TO_1E-2 | 1.49065268807031026e-01 | 1.36266126246915764e+00 | 3.85495513445842874e-01 | 3.67812930405358340e-01 | 1.98082473577756824e-01 |
| 1e-02 | found | 1130 | Q3_1E-4_TO_1E-2 | 1.43269244290135683e-01 | 1.36589201375409863e+00 | 3.82294889723597764e-01 | 3.72090066889873894e-01 | 1.96007198034510527e-01 |
| 0.05 | found | 1167 | Q3_1E-4_TO_1E-2 | 1.24237587729040713e-01 | 1.46977426294231739e+00 | 3.00474759679236647e-01 | 3.36900245566113365e-01 | 1.91955695663771020e-01 |
| 0.1 | found | 1188 | Q3_1E-4_TO_1E-2 | 1.54615870032162817e-01 | 1.29681606907522662e+00 | 4.19703963874024366e-01 | 3.81886431173036223e-01 | 2.13517547054683066e-01 |
| 0.2 | found | 1229 | Q4_GE_1E-2 | 1.52149315169303911e-01 | 1.42293477513353483e+00 | 3.07439181043833165e-01 | 3.09464566892228787e-01 | 1.91708797618798582e-01 |
| 0.5 | found | 1342 | Q4_GE_1E-2 | 1.54783244652477336e-01 | 1.52479029463764215e+00 | 3.00371502935150247e-01 | 2.66791387541483627e-01 | 1.99063921334719979e-01 |
| 0.8 | found | 1718 | Q4_GE_1E-2 | 4.95861628614275851e-01 | 1.16873874586829718e+00 | 1.47670576912387941e-01 | 1.11468817936436052e-01 | 6.66234436296587046e-02 |

## q解像度帯別の方向系譜

| q_resolution_band | point_count | first_step | last_step | median_D34_vs_late_overlap | median_D34_vs_late_max_angle_rad | median_max_D34_vs_Tperp_overlap | maximum_consecutive_max_angle_rad | basis_column_swap_count | sign_flip_total |
|---|---|---|---|---|---|---|---|---|---|
| Q0_LT_1E-8 | 697 | 0 | 880 | 2.16172576479320006e-01 | 1.44709103691088581e+00 | 2.63807871776682357e-01 | 1.50666577289229942e+00 | 0 | 101 |
| Q1_1E-8_TO_1E-6 | 286 | 89 | 982 | 1.58253000739000449e-01 | 1.37081869236138809e+00 | 3.36140068811611292e-01 | 1.84481830265565361e-01 | 0 | 5 |
| Q2_1E-6_TO_1E-4 | 93 | 983 | 1075 | 1.50836055730379337e-01 | 1.36344506877913374e+00 | 3.85314136451019640e-01 | 4.86755942611448867e-04 | 0 | 1 |
| Q3_1E-4_TO_1E-2 | 118 | 1076 | 1193 | 1.45973931521193812e-01 | 1.36320223438681931e+00 | 3.84331905715193822e-01 | 3.72167928512057267e-01 | 0 | 1 |
| Q4_GE_1E-2 | 3807 | 1194 | 5000 | 7.83470473086408603e-01 | 5.63782477763240863e-01 | 6.69770025827969018e-02 | 1.87356569883079943e-02 | 0 | 25 |

これらのq帯は数値解像度別表示であり、物理的方向成立閾値ではない。

## 記述的分類

**ROTATING_OR_MIXED_LINEAGE**

固定判定量:

- evaluation start: `983`
- D34対late overlap中央値: `0.7028877098650121`
- D34対late最大主角中央値: `0.6697580791073114` rad
- D34対Tperp seed別overlap中央値: `[0.05436309086624651, 0.06638024213580548, 0.06074580864567081]`
- transverse集約: `0.06638024213580548`

## 三つの問いへの回答

1. 急拡大前に微小に見えた方向部分空間は、急拡大後に最初に成立する三方向の部分空間へ連続しているか。

   単純な連続とは分類されなかった。比較可能範囲は回転または混合した系譜として記述された。

2. その微小方向部分空間は、準安定後の追加方向萌芽の横摂動部分空間に近いか。

   既存Tperpへの単純な一致とも分類されなかった。

3. 急拡大途中に、方向部分空間の交換・回転・選び直しが起きているか。

   固定step 900〜1400で、連続step間最大主角の最大値は `3.721679e-01` rad、projector distance最大値は `3.636360e-01`。同区間の列swapは `0` 回、sign flip総数は `4`。列交換・符号反転と2次元部分空間回転は別表で分離した。

## データから直接言えること

- Stage A0軌道のfとq診断をbitwise一致で再実行できた。
- P34によるD34部分空間系譜、D34_late、既存3個のTperp方向との主角・overlap・射影距離を全stepで比較できた。
- t0のS4とTperpは既存生成規則どおり数値誤差内で直交していた。
- 個別列の交換・符号反転と、2次元部分空間自体のambient回転を分離できた。

## データだけでは言えないこと

- q解像度床未満の最早期基底が物理的に確定した方向であるか
- 単一の三方向成立step
- H1/H2/H0の判定
- 追加方向萌芽が自然軌道上に実在する時刻
- 高精度親や別Delta/Nにおける同じ系譜

## 必須表

- [direction_basis_snapshots](../processed/direction_basis_snapshots.md)
- [direction_projector_quality](../processed/direction_projector_quality.md)
- [consecutive_subspace_continuity](../processed/consecutive_subspace_continuity.md)
- [early_vs_late_direction_overlap](../processed/early_vs_late_direction_overlap.md)
- [early_vs_transverse_overlap](../processed/early_vs_transverse_overlap.md)
- [late_direction_vs_transverse_overlap](../processed/late_direction_vs_transverse_overlap.md)
- [lineage_by_f_level](../processed/lineage_by_f_level.md)
- [lineage_by_q_resolution_band](../processed/lineage_by_q_resolution_band.md)
- [transverse_direction_reconstruction](../processed/transverse_direction_reconstruction.md)
- [numerical_health](../processed/numerical_health.md)

## 必須図

- `figure01_f_and_lineage_sampling_coordinates.png`
- `figure02_q_ratios_and_resolution_bands.png`
- `figure03_direction3_direction4_occupation.png`
- `figure04_D34_vs_late_principal_angles.png`
- `figure05_D34_vs_late_overlap.png`
- `figure06_D34_vs_Tperp_principal_angles.png`
- `figure07_D34_vs_Tperp_overlap.png`
- `figure08_late_D34_vs_each_Tperp_overlap.png`
- `figure09_consecutive_D34_max_angle.png`
- `figure10_consecutive_D34_projector_distance.png`
- `figure11_lineage_by_f_level.png`
- `figure12_lineage_by_q_resolution_band.png`
- `figure13_crossing_neighborhood_zoom.png`
- `figure14_metastable_region_zoom.png`
- `figure15_column_exchange_vs_subspace_continuity.png`

## 最終停止

Stage A2cの報告書完成をもって停止する。高精度、Delta掃引、`N=40`、`N=300`、Stage B/Cへ進まない。
