# Release Notes: 三次元正方配列の三角形ループ積における過冷却時間結晶背景、プローブ選別、および 4 段階・倍々エンタングルメント相転移（第 2 速報）

## DOI
- Concept DOI: https://doi.org/10.5281/zenodo.24618399
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.24618400
- Zenodo record: https://zenodo.org/records/24618400
- 前報（第 1 速報 Concept DOI）: https://doi.org/10.5281/zenodo.24614040

## 公開ファイル（Zenodo, 9点）
- `三次元正方配列の三角形ループ積における4段階相転移とプローブ選別_第2速報.md` / `triangle_loop_four_stage_transitions_probe_selection_rapid2_en.md`
- `三次元正方配列の三角形ループ積における4段階相転移とプローブ選別_第2速報.tex` / `triangle_loop_four_stage_transitions_probe_selection_rapid2_en.tex`
- `三次元正方配列の三角形ループ積における4段階相転移とプローブ選別_第2速報.pdf` / `triangle_loop_four_stage_transitions_probe_selection_rapid2_en.pdf`
- `build_tex_rapid2.py`
- `figures_and_programs_rapid2.zip`（全再現スクリプト 11 本、全結果 TXT/NPZ、日英 PNG 図 4 枚）
- `SHA256SUMS_rapid2`

## 主な内容
1. **512 マス状態の厳密解剖と自己修復**: 精度 7 桁で厳密直積 $H = O_8 \otimes M_8$（内側パウリ分解でスピン 2 的 $\sigma_x, \sigma_z \sim 10^{-16}$）、内部実効結合定数 $\rho$ は離散 7 値。背景 $+0.1$ 衝撃を 5 歩で $10^{-165}$ へ弾き、内部 1 マス穴あけを 1 歩で修復。
2. **時間結晶背景 $H$ による近傍プローブの 4 分岐と「衣」の形成**: $i, k$ 共有 $H^2 P$ 結合で、①大域複素波は 5 桁不変共存、②ヌル・カイラル波（$\sum P^2=0, I+\sigma_y$）は吸収されつつ $H$ 自身の $\|P\|$ を厳密 2 倍へ相転移、③スピン 2 波はエネルギーを 10 倍超へ吸収し周期 2 ロック、④局在凝縮体（単一サイト・$H_4$）は周囲に「衣（$E_{\mathrm{oth}}$）」をまとって安定束縛（$H_4$ は 1:1、単一サイト塊は大幅減速）。
3. **$\tau = 10,000$ フルスイープによる 4 段階・倍々エンタングルメント相転移と 2 つの解析定理**: Case 1（実数種）・Case 2（大域丸め誤差種：再加熱的バースト `512 → 3,072 → 1,560` と約 5,000 歩の第 3 層を経由）の両方で、$\text{rank } 1\;(1.0) \to \text{rank } 2\;(\frac{1}{\sqrt{2}}\times 2, \|P\|\times 2) \to \text{rank } 4\;(\frac{1}{2}\times 4, \|P\|\times 2, d_1=0)$ の静止固定点へ終着。第 3 層＝単一複素構造ロック、第 4 層＝3 重外積 $a_i b_j c_k$ × ベル状態ペアを解析証明。
4. **3 軸（虚数 3 次元）への自発的凝縮**: Case 2 の最終固定点（$\tau=10,000$）で実数軸 $A=0$ が `0.04%` へ消滅し、虚数 3 軸（$A=4,2,7$）が 84.84% を独占。

## 変更履歴
| 日付 | 内容 |
|---|---|
| 2026-11-28 | v1.0 公開。Zenodo 9 ファイル（DOI `10.5281/zenodo.24618400` / Concept `10.5281/zenodo.24618399`）、note・X・Facebook 原稿作成 |
