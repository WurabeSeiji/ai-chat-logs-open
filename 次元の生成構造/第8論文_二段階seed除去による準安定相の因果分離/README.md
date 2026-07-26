# 第8論文：二段階seed除去による準安定相の因果分離

第1予備実験（二段階seed除去対照）。第7論文の実験系を seed の有無以外そのまま再利用し、
初期seed(S0)と準安定横摂動seed(S1)を別々に除去した対照軌道でベースラインを確定する。

- 条件A: initial seed OFF / metastable seed OFF（Z0 = v）
- 条件B: initial seed ON / metastable seed OFF（Z0 = (v+δg)/‖·‖, δ=1e-15、以後自然発展のみ）
- 条件D: initial seed ON / metastable seed ON（第7論文の二段階seed基準）
- N = 5, 40, 300。共通最終step = 55000（第7論文と同一）。
- 独自解釈・追加解析（FFT/自己相関/Lyapunov/seed掃引/スケール不変性）は行わない。

構成: `instructions/` `code/` `config/` `raw/` `summary/` `figures/` `diagnostics/` `logs/` `reports/`
第7論文コードは read-only import で再利用（不変更）。SHA-256 は `config/source_file_hashes.json`。
