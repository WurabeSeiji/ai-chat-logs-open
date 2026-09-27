# 第6論文補講: G=q0=c=1 初期化同値検証

- 判定: **PASS**
- 比較対象: 保存済み5ケースの `raw_first_macro_microsteps.csv` の micro=0 行
- 新規側: `G=q0=c=1`, `q0=1`, `M/q0=20/3`, 等質量 `m_A=m_B=M/2` から初期化
- 比較方法: CSV読込後の各診断列の完全一致

| case | result | P | N | C | D | Q |
|---|---|---:|---:|---:|---:|---:|
| n1_attractive | PASS | 50.0 | 0.25 | 1.09 | 0.36 | 1.0 |
| n1_repulsive | PASS | 50.0 | 0.25 | 0.91 | 0.0 | 1.0 |
| n3_attractive | PASS | 50.0 | 0.25 | 1.81 | 3.24 | 1.0 |
| n3_repulsive | PASS | 50.0 | 0.25 | 0.19 | 0.0 | 1.0 |
| n4_attractive | PASS | 50.0 | 0.25 | 2.44 | 5.76 | 1.0 |

全ケースで micro, q_index, P, N, C, D, k1p, k2p, k3p, k4p, dP, Q が保存済み初期rowと完全一致した。
