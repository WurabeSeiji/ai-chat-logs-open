# RELEASE NOTES — 自己無撞着な関係波閉鎖系の基礎（全10本シリーズ）

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems
**著者:** 木原 範昭（WF System Co., Ltd.）／ORCID 0009-0004-6753-4020
**ライセンス:** CC BY 4.0

$N$ 個の存在の全二体関係 $M=N(N-1)/2$ を複素関係波とし、外部シードを除去した閉鎖系写像 $z\to\exp(\Delta\tau K)z$ の力学を、総括1本（論文0）＋基礎小論文9本（論文1〜9）として統一形式で書き、個別に DOI を取得したシリーズ。外部参照は Concept DOI（常に最新版へ転送）、版の明記は Version DOI。

## 論文一覧・DOI

| # | 日本語題（短縮） | English title (short) | Version DOI | Concept DOI |
|---|---|---|---|---|
| 0 | 総括：自己無撞着な関係波閉鎖系の力学 | Dynamics of a Self-Consistent Relational-Wave Closed System (capstone) | 10.5281/zenodo.22729147 | 10.5281/zenodo.22729146 |
| 1 | 単一波の運動方程式 | The Equation of Motion of a Single Wave | 10.5281/zenodo.22728761 | 10.5281/zenodo.22728760 |
| 2 | 点火の資格＝不安定な相対平衡 | The Qualification for Ignition Is an Unstable Relative Equilibrium | 10.5281/zenodo.22728824 | 10.5281/zenodo.22728823 |
| 3 | 90度床の解析的厳密解 | Analytic Exact Solution of the 90-Degree High-Symmetry Floor | 10.5281/zenodo.22728996 | 10.5281/zenodo.22728995 |
| 4 | ゼロ閉塞因数分解 | Zero-Closure Factorization of the High-Symmetry Floor | 10.5281/zenodo.22729081 | 10.5281/zenodo.22729080 |
| 5 | 自発的対称性の破れと縮退真空 | Spontaneous Symmetry Breaking and Degenerate Vacua | 10.5281/zenodo.22729089 | 10.5281/zenodo.22729088 |
| 6 | インフレは床の線形不安定 | The Inflation Curve Is the Amplification of the Linearly Unstable Mode | 10.5281/zenodo.22729096 | 10.5281/zenodo.22729095 |
| 7 | ロック後の幾何 | Post-Lock Geometry: Double Rotation and the 4-Dimensional Span | 10.5281/zenodo.22729108 | 10.5281/zenodo.22729107 |
| 8 | 加法法則 a∞=√(H/M) | The Additivity Law of the Terminal Equal Amplitude | 10.5281/zenodo.22729116 | 10.5281/zenodo.22729115 |
| 9 | 無名性と情報の最小単位（検討論文） | Anonymity and the Minimal Unit of Information (exploratory) | 10.5281/zenodo.22729128 | 10.5281/zenodo.22729127 |

## 収録ファイル（各論文）

- 本文：`paperN_<slug>_ja.md` / `paperN_<slug>_en.md`（日本語正本は `論文N_*.md`）
- 組版：`paperN_<slug>_{ja,en}.tex` / `paperN_<slug>_{ja,en}.pdf`（lualatex・ltjsarticle、エラー0で検証）
- 再現：`paperN_<slug>_repro.zip`（プログラム＋参照図＋分析md/README/SHA256SUMS/run_all.sh＋小集計CSV。大容量の生データ〔.npy/.npz/大CSV/入れ子zip/対話HTML〕は同梱せず `_LARGE_FILES_SHA256.txt` に SHA-256 を記録）

## 親論文（下敷き・既公開）

- 自己無撞着な関係波閉鎖系におけるインフレーション的急拡大の機構（v2）: Version 10.5281/zenodo.22176949 / Concept 10.5281/zenodo.22112008
- 開始様式判別（増幅⟺不安定な相対平衡）: Concept 10.5281/zenodo.21798854

## 力学の正本

全論文で同一：`run_N3_N40_stage123_v1.py`（SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`）、物理式無変更。

## 変更履歴

- 2026-09-13 v1：全10本を Zenodo に初版公開（論文1→…→9→0 の順、外部参照は Concept DOI）。日英 md/tex/pdf＋再現zip をアップロード。Zenn 記事 `relational-wave-closed-system-foundations` を追加。
