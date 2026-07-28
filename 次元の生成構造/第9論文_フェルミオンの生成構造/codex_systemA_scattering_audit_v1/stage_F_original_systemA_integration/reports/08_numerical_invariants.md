# 08 数値不変量

全F-B実行の最大残差:

- `unitarity_residual`: `3.3306690738754696e-16`
- `coefficient_orthogonality_residual`: `1.2174184636272442e-16`
- `path_sum_residual_A`: `2.7703533911349609e-15`
- `path_sum_residual_B`: `2.4420569733063502e-15`
- `total_norm_conservation_residual`: `4.9960036108132044e-15`
- `demodulation_reconstruction_residual`: `4.3958451321218078e-15`
- `parity_projection_sum_residual`: `2.9976021664879238e-15`

- NaN/Inf: `0`
- theta範囲違反: `0`
- 数値不安定判定閾値: `1.0e-08`

経路干渉項は今回も直交eta由来で機械精度内に抑えられた。`B_to_A_transfer`は `spectral_similarity_to_initial_B; not path flux` であり、経路フラックスではない。
