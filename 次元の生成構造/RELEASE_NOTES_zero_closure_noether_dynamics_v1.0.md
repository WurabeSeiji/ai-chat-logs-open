# Release Notes: 離散零閉包からの Noether 保存則と関係位相動力学 —— 零閉包を厳密保存する離散自己写像、N→∞ 連続場方程式、局所ゲージ幾何、および標準模型一世代表現と chirality 選択

## DOI
- Concept DOI: https://doi.org/10.5281/zenodo.22040735
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.22040736
- Zenodo record: https://zenodo.org/record/22040736
- 公開日: 2026-08-21

## Zenn
- 記事: https://zenn.dev/noriaki_kihara/articles/zero-closure-noether-dynamics
  （元原稿: `articles/zero-closure-noether-dynamics.md`）

## note
- 日本語版: https://note.com/kiharanoriaki/n/n91202fa73800
  （元原稿: `次元の生成構造/電子の反跳実験/note_zero_closure_noether_dynamics_ja_v1.md`）
- 英語版: https://note.com/kiharanoriaki/n/n6db6d9f175d8
  （元原稿: `次元の生成構造/電子の反跳実験/note_zero_closure_noether_dynamics_en_v1.md`）

## Facebook
- 日本語版: （公開後に追記）
  （元原稿: `次元の生成構造/facebook_zero_closure_noether_dynamics_ja_v1.md`）

## 公開ファイル（Zenodo, 6点）
- `zero_closure_noether_dynamics_ja_public_v1.0.md` / `zero_closure_noether_dynamics_en_public_v1.0.md`
- `zero_closure_noether_dynamics_ja_public_v1.0.tex` / `zero_closure_noether_dynamics_en_public_v1.0.tex`
- `zero_closure_noether_dynamics_ja_public_v1.0.pdf`（63ページ） / `zero_closure_noether_dynamics_en_public_v1.0.pdf`（48ページ）

作業原本: `次元の生成構造/電子の反跳実験/zero_closure_noether_SM_chirality_theory_and_test_spec_v9.md`（思考実験ログ: 同フォルダ `CHATGPTネータ定理.md`）。公開版は v9 にヘッダー（DOI 等）を付加し、参考文献 [6][7][8] を検証可能な文献（PDG 2024 / Peskin–Schroeder / Wilson 1974）に差し替え、[9] Georgi–Glashow 1974 と [10] 前論文 DOI を追加したもの。本文は v9 と同一。

## 位置づけ

「次元の生成構造」シリーズ、閉包公理からの対称性導出 v1.0（Concept DOI 10.5281/zenodo.22028072）の続編。前論文で残った二課題──Noether 型保存則と、次状態を決める動力学──を扱う。

## 主な内容

- **動力学の中心主張**: 状態の書換えを動力学として暗黙に許さない。自己無撞着が基礎条件なので、動力学は許容状態空間 𝒵_N の自己写像 𝓕_N: 𝒵_N→𝒵_N として構成されねばならない
- **離散 current と離散作用**: J_ij = A² sin(φ_j−φ_i)、S_N = −A² Σ cos(φ_j−φ_i)、停留条件＝離散 continuity equation（有限 N で厳密な離散 Noether 保存則）
- **零閉包保存自己写像**: 零閉包 Σe^{2iφ}=0 の接空間射影 P_φ と有限 retraction R_φ により、各反復で Σ X_i² = 0 を厳密保存
- **s ≠ t**: 自己写像パラメータ s は構成・選択パラメータであり物理時間ではない。物理時間 t は固定点配置の Lorentz 読出し（虚軸 it）の内部に含まれる
- **連続極限**: N→∞ で ∂_μ J^μ = 0、S → (A²/2)∫ g^{μν}∂φ∂φ、□φ = 0（一般形は制約付き wave/Laplace 方程式）
- **局所ゲージ幾何**: 局所位相原点の無名性 → 辺 connection θ_ij → simplex 面 curvature Θ_ijk → 連続極限で D_μ, F_μν, Maxwell 作用、C^r 化で Yang–Mills 非可換項
- **標準模型一世代表現**: S(U(3)×U(2)) の trace-zero 条件 3y₃+2y₂=0 で hypercharge 比固定、V*⊕Λ²V = d^c⊕L⊕u^c⊕Q⊕e^c（15 左手 Weyl 成分、5̄⊕10 と同型だが SU(5) は仮定しない）、全 perturbative anomaly と SU(2) global anomaly の相殺を直接検算
- **chirality 選択**: 共役二 Weyl sector を過去論文の A/B 二チャネル選択系と同定、mirror-odd 相関 J=Im(B*²C) と非線形選択項から最小 normal form Ṡ_χ = λJ + gS_χ(1−S_χ²)
- **数値検証仕様（§21）**: 既存 A/B Fermi 型実験系を変更せず、観測量追加と mirror run で検証。**公理保存監査**（零閉包・有限回帰・simplex 閉包・ノルム・current 保存）を動力学採用の必須条件とする

## 残る課題

- §21 の数値検証（後続論文）
- 3⊕2 読出しの選択則の一般化、V* 選択の根拠、N 自体の選択則
- Higgs 動径モード・Yukawa・世代数・質量階層・量子補正（本論文の範囲外）
