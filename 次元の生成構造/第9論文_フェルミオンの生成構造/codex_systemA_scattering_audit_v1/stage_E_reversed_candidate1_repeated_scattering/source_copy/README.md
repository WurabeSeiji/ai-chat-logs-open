# System A実験コピーの由来

`run_system_A_localization_exchange_R_sweep_preliminary_v1_ORIGINAL_SNAPSHOT.py`
は、既存System A原本を実行せず、Stage E開始時にバイト単位で複製した
読み取り専用スナップショットである。

- 原本SHA-256:
  `91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec`
- スナップショットSHA-256:
  `91a1a19a5e11be80626b34630e353fccc59b0197782c2fcd5417b9e18a2766ec`

実際のStage E実験実装は
`../code/system_A_experimental_copy.py` に自己完結形で置いた。
原本・スナップショットのどちらもimportまたは実行していない。
