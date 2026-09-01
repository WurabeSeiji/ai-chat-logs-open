# denominator N±2 controls (2026-09-01)

目的: `Δτ=2π/N` が特別かどうかを、隣接分母 `N-2,N-1,N,N+1,N+2` で対照する。

## 固定条件
- N=3..16
- STEPS=500
- seed追加なし
- 初期状態は同じ hm 保存状態
- 相互作用 `H_ef=A_ef conj(z_e) z_f`
- frozen-H Hermitian eigendecomposition による `exp(-i Δτ H)` 更新
- 変更点は分母だけ

## 監査上の注意
この系列は、後に元の denominator=124 の正本との比較手順に混乱が生じたため、現時点では「未検証派生実験」として保存する。研究上の証拠として採用する前に、正本 `pass2_run.py` の短縮版で state-by-state 再現を確認すること。
