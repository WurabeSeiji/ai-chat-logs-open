# 第二思考実験 v1.2 — 実数配列版

日付: 2026-09-21

## 目的

旧版は入力一覧で「複素数を使わない」と宣言していた一方、T14〜T16、T20、T22 の証明と検証コードで複素数演算を計算上の近道として使用していた。v1.2 はこの不整合を除去し、二つの実数値 `a,b` と実行列・実配列演算だけで全導出と検証を再構成する。

## 中核の実写像

`X=(a,b)^T`、`J0=[[0,-1],[1,0]]` とし、

```text
B(X) = a I + b J0 = [[a,-b],[b,a]]
Phi(X) = B(X) X = [a^2-b^2, 2ab]^T
r = a^2+b^2
```

を使う。`J0^2=-I` のマイナスは実相互作用演算の行列に明示的に含まれ、複素単位を入力しない。

## v1.2 のファイル

- `thought_experiment_02_coulomb_readout_ja_v1.md` — 修正版論文本文（日本語）
- `thought_experiment_02_coulomb_readout_en_v1.md` — 修正版論文本文（英語）
- `thought_experiment_02_coulomb_readout_{ja,en}_v1.tex` / `.pdf` — 組版した論文（日英）
- `run_coulomb_readout_experiments_v1.py` — 数値検証 E1〜E7。実数配列のみ
- `coulomb_readout_experiments_results_v1.json` — v1.2 の数値結果
- `verify_coulomb_readout_identities_sympy_v1.py` — T13〜T22 の記号検算。実行列のみ
- `make_coulomb_readout_figures_v1.py` — 図1〜図5生成。幾何計算は実数配列のみ
- `fig01...fig05 ... _v1.svg/.png` — 日本語・英語の図
- `run_v1_stdout.txt`, `verify_v1_stdout.txt`, `figures_v1_stdout.txt` — 再実行ログ
- `verification_v1_real_array.txt` — 監査結果
- `SHA256SUMS_v1.txt` — 固定点ハッシュ

## 再現順序

```bash
python3 verify_coulomb_readout_identities_sympy_v1.py
python3 run_coulomb_readout_experiments_v1.py
python3 make_coulomb_readout_figures_v1.py
```

## 固定した結論

複素演算を除去しても T14〜T22 の恒等式と数値判定は維持された。数値検証の PASS/FAIL は旧版と一致し、既知の E4c だけが FAIL のままである。したがって主結果は複素数の仮定を必要としない。ただし、旧版の導出と実装は宣言と不整合であり、v1.2 で初めて「複素数を使わない」という方法上の要求と一致した。
