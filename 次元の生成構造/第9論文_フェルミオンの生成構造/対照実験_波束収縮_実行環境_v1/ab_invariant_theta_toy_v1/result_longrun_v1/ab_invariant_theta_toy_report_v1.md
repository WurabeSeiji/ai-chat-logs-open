# AB 回転不変量 theta 自動生成トイモデル v1 実行報告

- 衝突回数: 256
- high_n: 63
- 外部 R/theta 入力: なし
- 更新: 実直交 2x2 AB 回転、個別正規化なし

| case | theta(0) | R(0) | max theta drift | max closure drift | verdict |
|---|---:|---:|---:|---:|---|
| fundamental_control | 0 | 0 | 0.000e+00 | 0.000e+00 | PASS |
| even_boson_control_B62 | 0 | 0 | 0.000e+00 | 0.000e+00 | PASS |
| odd_fermion_candidate_B63 | 0.761952071831 | 0.4765625 | 2.220e-16 | 1.289e-16 | PASS |

この判定は数値的不変性の確認であり、ボゾン・フェルミオン対応の物理的実証ではない。
また theta が一定になることは同じ U の反復を与えるが、有限 n での U^n=I を自動的には保証しない。
