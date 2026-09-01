# denominator=124 recheck (2026-09-01)

目的: 元の `Δτ=2π/124` 条件を 500 step で再計算し、後続の分母変更実験との基準にする。

## 固定条件
- N=3..16
- STEPS=500
- 分母=124, `Δτ=2π/124`
- seed追加なし
- `H_ef=A_ef conj(z_e) z_f`
- frozen-H の Hermitian eigendecomposition による `exp(-i Δτ H)` 更新
- 観測 `Hperp/H` は初期状態の real/imag 2-plane を基準に計算

## 重要な監査上の注意
この再計算は、元の正本 `pass2_run.py` をそのまま実行したものではなく、同じ数式を再実装した監査コードである。したがって、正本再現性の最終確認には、元 `pass2_run.py` を500 stepへだけ短縮して state-by-state 比較する追加監査が必要。

N=6 の500 step値は正本CSVとほぼ同桁で整合したが、これをもって全系列のbitwise再現とみなしてはいけない。
