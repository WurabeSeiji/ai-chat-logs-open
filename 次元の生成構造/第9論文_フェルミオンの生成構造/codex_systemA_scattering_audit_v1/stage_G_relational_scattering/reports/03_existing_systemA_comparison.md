# 既存System A代表条件の3モード比較

C0の既存正規化再現は最大絶対誤差0、Stage FのC0/reversed_C1 raw系列との12回帰比較も最大絶対誤差0で通過した。主系列はStage Fの結果に従い `raw_update` とした。

| mode | kappa | Gamma range | R_eff range | min L gap@col | min N_eff gap@col | class |
|---|---:|---:|---:|---:|---:|---|
| C0 | 0.01 | 1.402e-15 | 0 | 0.00011479@110 | 0.70804@110 | relational_term_inactive |
| C0 | 0.1 | 1.402e-15 | 0 | 0.00011479@110 | 0.70804@110 | relational_term_inactive |
| C0 | 1 | 1.402e-15 | 0 | 0.00011479@110 | 0.70804@110 | relational_term_inactive |
| relational_C1 | 0.01 | 3.712e-15 | 0 | 7.4823e-05@125 | 0.46151@125 | constant_relation_reparameterization |
| relational_C1 | 0.1 | 2.914e-16 | 0 | 3.8972e-05@46 | 0.24038@46 | constant_relation_reparameterization |
| relational_C1 | 1 | 3.331e-16 | 1.11e-16 | 9.7353e-05@101 | 0.60048@101 | constant_relation_reparameterization |
| reversed_C1 | 0.01 | 3.306e-15 | 0 | 4.9746e-05@29 | 0.30684@29 | constant_relation_reparameterization |
| reversed_C1 | 0.1 | 1.551e-15 | 1.11e-16 | 2.4706e-05@44 | 0.15239@44 | constant_relation_reparameterization |
| reversed_C1 | 1 | 7.737e-16 | 2.22e-16 | 4.3848e-05@89 | 0.27046@89 | constant_relation_reparameterization |

`N_A=1,N_B=63`では `c_A=c_B=-1` が保存された。relational_C1の補正はreversed_C1より小さいが、衝突ごとに動く補正にはならなかった。
