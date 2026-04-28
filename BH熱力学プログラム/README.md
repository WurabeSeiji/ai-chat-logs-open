# BH熱力学プログラム

中心投影と離散性という最小限の仮定から、4+1次元 Tangherlini ブラックホール熱力学を再構成する6篇の研究プログラム。

**著者**: 木原 範昭（Noriaki Kihara）／WF System Co., Ltd.／ORCID [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**開始日**: 2026-04-28
**ステータス**: 草稿フェーズ

---

## プログラム構成

| # | 主題 | 主要結果 | 状態 |
|---:|---|---|---|
| 1 | 俯瞰・仮説提示 | 仮説1（写像）+ 仮説2（離散ドリフト） | 全章試案完成（日英）|
| 2 | 中心投影の数学 | グノモン投影下の Schwarzschild–de Sitter 計量 | 全章試案完成（日英）|
| 3 | 離散構造の数学 | $\Delta(R) \sim (16\pi/3) R^3$, $c = 8/(3\pi)$ | 全章試案完成（日英）|
| 4 | エントロピー対応 | $\Delta/S_{BH} \to 32/(3\pi) \approx 3.40$ | 全章試案完成（日英）|
| 5 | 動的接続 | $R \propto M^{1/2}$, $T_H \propto M^{-1/2}$, $\tau \propto M^{3/2}$ | 全章試案完成（日英）|
| 6 | 統合と 4+1→3+1 削減 | 候補機構3つ（KK/ブレーン/ホログラフィー） | 全章試案完成（日英）|

## ディレクトリ

```
BH熱力学プログラム/
├── README.md                       # このファイル
├── papers/
│   ├── paper1_overview/            # 論文1：俯瞰・仮説提示
│   │   ├── outline.md
│   │   ├── draft_introduction.md
│   │   ├── draft_section3_hypotheses.md
│   │   └── draft_section5_critiques.md
│   ├── paper2_projection/          # 論文2：中心投影の数学
│   │   └── outline.md
│   ├── paper3_packing/             # 論文3：離散構造の数学
│   │   ├── outline.md
│   │   └── draft_sections_2_4.md
│   ├── paper4_entropy/             # 論文4：エントロピー対応
│   │   └── outline.md
│   ├── paper5_dynamics/            # 論文5：動的接続
│   │   └── outline.md
│   └── paper6_synthesis/           # 論文6：統合または残課題
│       └── outline.md
├── computations/                   # 数値・解析計算
│   ├── packing_numerical.py        # N(k) の対称性付き列挙
│   ├── asymptotic_fit.py           # Δ(R) の多項式フィット
│   ├── jacobi_verification.py      # Jacobi 接続の検証
│   └── results/
│       ├── N_table_k0_30.tsv
│       ├── N_table_k0_60.tsv
│       └── math_report.md          # 数学レポート（c の理論的決定）
└── notes/
    ├── strategy.md                 # 戦略原則・進捗管理
    └── correspondence_record.md    # 対応・通信記録
```

## 主要結果（確立済み）

### 数学的核心結果

充填条件
$$\sum_{i=1}^{4}(|x_i|+1/2)^2 \le (2k+1)^2$$
を満たす整数組の個数 $N(k)$ について：

1. **正確な数値**（[`N_table_k0_60.tsv`](computations/results/N_table_k0_60.tsv)）：
   - $N(0) = 1$, $N(1) = 137$, $N(10) = 818{,}145$, $N(60) = 1{,}028{,}515{,}513$
2. **漸近係数の理論値**：$c = 8/(3\pi) \approx 0.84883$
3. **副主要項**：$\Delta(R) = (16\pi/3) R^3 - 6\pi R^2 + O(R)$
4. **数値検証**：$k_{\min}=30$ の3次フィットで $c$ は理論値と $0.024\%$ 一致
5. **Jacobi 接続**：$N(k)$ は Jacobi の四平方表現数 $r_4(N)$ を用いて閉形式に近い表現が可能（境界補正の精密化は **後で検討**）

詳細：[`computations/results/math_report.md`](computations/results/math_report.md)

### 物理的対応

- 4+1 次元 Tangherlini Bekenstein–Hawking エントロピー：$S_{BH} \propto R^3$（$R = r_h$）
- $\Delta(R) \propto R^3$ と一致するスケーリング
- 漸近的係数比：$\Delta/S_{BH} \to 32/(3\pi) \approx 3.40$
- $R \propto M^{1/2}$, $T_H \propto M^{-1/2}$, 寿命 $\propto M^2$（4+1 Tangherlini と一致）
- 蒸発レート指数の不一致は残課題として明示

## 戦略原則

1. **自己完結性**: 各論文は外部論文への引用を極力ゼロ。
2. **自己引用ゼロ**: 過去の自分の論文への参照は一切しない。
3. **統一モデルの背景化**: $xyztRQ$ 統一モデルは前面に出さない。本プログラムは $xyztR$ 5軸幾何学に集中。
4. **独立した経路**: ホログラフィー原理、AdS/CFT 等を仮定しない。
5. **残課題の明示**: 完全一致を目指さず、不一致点は今後の課題として明示。
6. **シンプルな仮定**: 中心投影 + 離散性のみ。

戦略原則の詳細：[`notes/strategy.md`](notes/strategy.md)

## 進捗

完了：
- フォルダ構造作成
- $N(k)$ の数値計算（$k = 0..60$）
- 漸近係数 $c = 8/(3\pi)$ の理論的決定と数値検証
- 副主要項 $-6\pi R^2$ の理論導出
- Jacobi 接続の検証コード（直接列挙との一致確認）
- 数学レポート完成
- 論文1 全章試案（§1 序論、§2 既存研究レビュー、§3 仮説、§4 検証プログラム、§5 批判先回り、§6 結論）日英完成
- 論文2 全章試案（§1–§7、付録）日英完成
- 論文3 全章試案（§1–§8、§2–§4 詳細版）日英完成
- 論文4 全章試案（§1–§5）日英完成
- 論文5 全章試案（§1–§6）日英完成、蒸発レート指数の混同解明（§5.1 で 3+1 vs 4+1 規約混同を明示）

未完了：
- 各論文の参考文献の正確な書誌情報
- LaTeX 化、PDF 生成
- Zenodo 公開
- 一部の **後で検討** マーカー：Christoffel 記号の閉形式（論文2 付録 A）、Jacobi 接続の正確な係数（論文3 §6.3）、$G_5$/$\ell_P$ 規格化詳細（論文4・5）、$E(R)$ の振動構造解析（論文3 §7.2）、4+1→3+1 削減機構の定量化（論文6 §4）

## ライセンス

CC BY 4.0（予定）
