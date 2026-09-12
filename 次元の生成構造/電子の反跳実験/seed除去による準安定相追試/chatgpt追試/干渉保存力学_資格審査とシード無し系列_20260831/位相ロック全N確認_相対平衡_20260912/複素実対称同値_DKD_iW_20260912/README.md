# 90°骨格：複素自己無撞着問題 ≡ 実対称0-1行列Wの固有値問題（D^{-1}KD=iW）＋Perron–Frobenius

日付: 2026-09-12。読出しのみ・物理無変更・新規走行なし。論文1 §5 の核心の追加結果を全N=3-40で検証。

## 主張
D=diag(e^{iφ_e}) とすると、90°骨格（ψ∈(π/2)Z, sin2ψ=0）で sinψ·e^{iψ}=i sin²ψ ゆえ D^{-1}KD=iW（W_ef=A_ef sin²ψ_ef, 実対称0-1・非負）。
よって Kz=iωz（z=Dr）⇔ Wr=ωr。W既約なら Perron–Frobenius で正振幅主固有ベクトル r>0 が一意（ω=ρ(W)）、z=Dr が自己無撞着床。

## 実行
```
bash run_all.sh
```
正本 run_N3_N40_stage123_v1.py（SHA256 1abf2353…）から adjacency を import し、整数K床・make_parent床（N=3-40）で
‖D^{-1}KD-iW‖、Wr=ωrの一様性、ω=ρ(W)=σ(iK)、r>0、Perron重なり、W連結成分を算出。

## 結果（詳細は REPORT.md）
両床・全Nで D^{-1}KD=iW（≤2e-12）、ω=ρ(W)=σ一致、r>0、Perron重なり=1.000000、W連結（成分1）。
→ 90°骨格上の自己無撞着正振幅波は Perron–Frobenius により一意に構成できる（全N=3-40）。

## 依存・入力
numpy, scipy。入力は読み取りのみ（既存床 npz）。出力は results/。
