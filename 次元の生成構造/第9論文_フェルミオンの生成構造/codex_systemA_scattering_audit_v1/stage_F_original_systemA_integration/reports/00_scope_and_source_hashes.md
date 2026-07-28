# 00 範囲と参照正本

Stage Fの全出力は本ディレクトリ内に限定した。既存System A/B、Stage E、既存CSV・図・報告書は読取りだけで、変更していない。

## 正本

- `次元の生成構造/第9論文_フェルミオンの生成構造/対照実験_波束収縮_実行環境_v1/20260713/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py`
  - SHA-256: `f815320f5632ae1b23ccade3a53b01e9110ad770da8407e2b10fa6065ef1695c`
- `次元の生成構造/第9論文_フェルミオンの生成構造/対照実験_波束収縮_実行環境_v1/20260715/run_system_A_localization_exchange_R_sweep_preliminary_v1.py`
  - SHA-256: `91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec`


事前監視対象は合計9ファイル。事後値は `reference_hashes_after.json` と `reference_hash_comparison.json` に保存する。

20260713正本は `N_A=1,N_B=63`、128衝突、R点 `(0,0.51,0.55,0.60,0.70,0.90,1)` の直接参照元である。31系列は別条件の20260715 custom packet `A=(1), B=(1,2), R=0.697177927` に属する。

## 事後照合

9参照すべてについて、path・size・mtime・SHA-256が事前値と一致した。 `all_references_unchanged=true`。
