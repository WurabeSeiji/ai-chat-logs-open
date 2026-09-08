# 高対称理論床・偶数系光子対精密解析 2026-09-08

高対称巡回 Fourier 相対平衡 `z_ij=r_d exp(i 2π(i+j)/N)` を定義式から再構成し、偶数 N=4..40 の全波pairについて `z_b=±i z_a` を監査する。

## 実行

```bash
./run_all.sh
```

## 入力仕様

外部NPZは使わない。理論床の定義式と縮約固有値問題をコード内で再構成する。比較基準文書は同研究系列の `theoretical_floor_derivation_N3_N40_20260906.md`。

## 主結果

- exact pair 出現: N=8,16,24,32,40 のみ
- 各Nで全M波をM/2 pairへ完全matching可能
- exact候補pair数: `N(4N-11)/2`
- 候補グラフは bipartite
- component size は 4/8/16 のみ
- `N=8k` 選択則は位相90°条件 `Δ(i+j)=N/4` と parity `i+j≡d (mod2)` の整合から説明できる
- unordered pair の ±i 符号は順序反転で反転するため、右旋/左旋の物理的個数には別の orientation が必要
