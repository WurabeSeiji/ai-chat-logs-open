# make_parent 段階1 振幅選択の再現解析
## 目的
`original_engine.py` を変更せず、`_make_parent_phase_only` の1反復を外部診断として再現し、振幅を決める `v=W y` のノルム、`JG` 固有値、シンプレクティック偏極量を記録する。
## 固定条件
- N=3..16
- 反復 400
- beta=0.5
- seed = 40260721+1000*N
- `v` の正規化なし
## 恒等式の検算
各反復で `JG y=-i sigma y`, `||y||=1`, `v=W y`, `G=W^T W` より `||v||^2 = sigma * chi`, `chi=i y^† J y`。CSV の `identity_error` で確認。
## 収束残差の注意
現行 `_eigenmode_residual` は `mu=v^†(iKv)` を使い、`v` 非正規化後も分母 `v^†v` がない。比較のため現行残差と正しい Rayleigh quotient を使うスケール不変残差を両方保存した。元コードは変更していない。
## N=3..6 exact 候補
- N=3: {'sigma2': '1/2', 'norm2': '1/5', 'r2': '1/15', 'Nr2': '1/5', 'sigma': 'sqrt(2)/2', 'chi': 'sqrt(2)/5'}
- N=4: {'sigma2': '2', 'norm2': '2/3', 'r2': '1/9', 'Nr2': '4/9', 'sigma': 'sqrt(2)', 'chi': 'sqrt(2)/3'}
- N=5: {'sigma2': '7/2', 'norm2': '49/57', 'r2': '49/570', 'Nr2': '49/114', 'sigma': 'sqrt(14)/2', 'chi': '7*sqrt(14)/57'}
- N=6: {'sigma2': '6', 'norm2': '6/5', 'r2': '2/25', 'Nr2': '12/25', 'sigma': 'sqrt(6)', 'chi': 'sqrt(6)/5'}
## 科学的注意
`exact_candidates_N3_N6.json` は倍精度値からの代数数同定であり、記号的証明ではない。今後の証明対象は位相写像上の不変軌道と `y^†Gy` の閉形式。
