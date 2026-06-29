# RELEASE NOTES — simplified_double_slit（二重スリット思考実験シリーズ）

位置揺らぎをもつ光源による二重スリット遠方場干渉の思考実験シリーズ。光源位置分布が観測量（縞シフト）へ写される**押し出し（形の保存）**を扱い、これを**局在光源**へ拡張する。すべて指定幾何（$L=10,W=5$）上の厳密計算。物理法則の導出・測定問題の解決は主張しない。外部引用は確立した教科書事項に限る。

---

## 論文1：位置揺らぎ光源による二重スリット干渉 ― 光源位置分布の縞シフト量分布への押し出し（形の保存）（v0.3）

A Thought Experiment on Double-Slit Interference from a Source with Positional Fluctuation: Push-forward of the Source-Position Distribution to the Fringe-Shift Distribution (Shape Preservation)

- **Concept DOI**: 10.5281/zenodo.21035808（外部参照用・最新版へ自動転送）
- **Version DOI (v0.3, 最新)**: 10.5281/zenodo.21035809
- **Zenodo**: https://zenodo.org/records/21035809
- **公開日**: 2026-06-29
- **ライセンス**: CC BY 4.0
- **位置づけ**: 思考実験（モデル計算）。単一波長点光源の位置揺らぎ $P(y)=\cos^2$ が、各試行で同形の縞を幾何学的経路差だけシフトさせ、反復試行の縞シフト分布が $P$ の押し出し（近軸線形ゆえ形を保存）になることを厳密計算で例示。観測モード (a)（縞シフトを読む）と (b)（強度積算＝可視度低下・vCZ）を峻別。

### 収録ファイル（6）

- `paper_doubleslit_position_readout_ja_v0_3.md` / `.tex` / `.pdf`（日本語）
- `paper_doubleslit_position_readout_en_v0_3.md` / `.tex` / `.pdf`（英語）
- 作図：`fig_setup_double_slit.py`, `fig_setup_source_uncertainty.py`, `fig_decomposition_static.py`, `fig_shift_histogram.py`

### 内容

- 厳密位相 $\Phi_k(s;y)=\frac{2\pi}{\lambda_0}[r_k(y)-y_{{\rm slit},k}s]$、$I=|e^{i\Phi_1}+e^{i\Phi_2}|^2$。$y$ 依存は位相オフセット $\Delta r(y)$ のみ＝波形は厳密に同形でシフト。
- 縞シフト $u(y)=2\pi\Delta r/\lambda_0\approx-2\pi Wy/(L\lambda_0)$（近軸線形）。押し出し $\rho(u)=P(y(u))|dy/du|$ が $\cos^2$ 形を保存、端の $\sim3\%$ ずれは非近軸 $\tfrac12\tan^2\theta$。
- 一意性 $|u|<180^\circ$（エイリアスなし）、arcsine（射影座標）との座標差の明示。

---

## 論文2：局在奇数倍音光源による二重スリット干渉 ― 形の保存は条件付きで脆く、単一波長 N=1 が頑健な特別な場合（v0.1）

A Thought Experiment on Double-Slit Interference from a Localized Odd-Harmonic Source: Shape Preservation Is Conditional and Fragile, and the Single-Wavelength N=1 Is the Robust Special Case

- **Concept DOI**: 10.5281/zenodo.21035830（外部参照用・最新版へ自動転送）
- **Version DOI (v0.1, 最新)**: 10.5281/zenodo.21035831
- **Zenodo**: https://zenodo.org/records/21035831
- **公開日**: 2026-06-29
- **ライセンス**: CC BY 4.0
- **自己参照**: 論文1（Concept DOI 10.5281/zenodo.21035808）を `isSupplementTo` で参照（唯一の自己参照）
- **位置づけ**: 思考実験。論文1 の押し出しを**局在光源**（奇数倍音孤立ピーク波 $S_N$）へ拡張。正味は**否定的結果**＝局在化は押し出しを改善せず、整列拘束と off-axis 脆さを加える。形の保存は条件付き・脆く、$N=1$ が頑健な特別な場合。

### 収録ファイル（6）

- `paper_localized_source_fluctuation_ja_v0_1.md` / `.tex` / `.pdf`（日本語）
- `paper_localized_source_fluctuation_en_v0_1.md` / `.tex` / `.pdf`（英語）
- 作図・計算：`fig_oddharm_interference.py`（パラメータ化：`--L --W --lam0 --N --halfdeg --dlam --ngrid`）、`make_paper2_localized_wave.py`
- 図：`fig_paper2_localized_wave_N17`（$S_N$）、`fig_oddharm_interference_L10_W5_lam1p0308_N17_pm180`（整列形保存干渉）、`fig_setup_source_uncertainty`（論文1 揺らぎ系）、`fig_oddharm_fluct_L10_W5_lam1_N1_dlam1`（N=1）、`..._lam1_N17_dlam1`（非整列散乱）、`..._lam1p0308_N17_dlam1`（整列）

### 内容

- 整列条件 $\sqrt{L^2+(W/2)^2}=(m/2)\lambda_0$（スリットがスパイクに乗る）。整列時 $I=4S_N(\theta)^2$（局在波の二乗＝中央に鋭い局在縞、周期 $180^\circ$ 縞列の中央）。
- 許容帯 $\sim1/(2N)$（$N$ 大で狭い）。$\lambda_0$ の約3%差で干渉成立／散乱が分かれる。中央整列の許容帯は $W/L$ 非依存、$W/L$ 大は off-axis 脆さ（$|dr_1/dy|_{y=0}\propto W$）。
- 揺らぎ下：$N=1$ は無条件で論文1 と機械精度一致（検証済み）。整列 $N\ge2$ は $\cos^2$ を継承するが off-axis 散乱で縁が下方へずれる（(i) 良性非線形性＋(ii) 新規散乱）。非整列は $y=0$ から散乱。
- 空間コヒーレンス／van Cittert–Zernike（モード b＝可視度低下）との区別を明示（本稿はモード a＝押し出し）。

### 改訂履歴

- **v0.1 (2026-06-29)**: 初版公開。AI 査読（claude.ai）2 ラウンド反映：(4.2) 長さ/位相分離、§3.3 の $W/L$ 帰属是正（中央許容帯は $W/L$ 非依存／$W$ 効果は off-axis）、ずれの (i) 非線形性＋(ii) 散乱分解、「単一縞」を周期列の中央に限定＋位相ロック注記、表の包絡基準／重み基準の脚注、否定的結果としての前面化（題・要旨・結論）。$N=1$ で論文1（fig_decomposition_static）と全行＋曲線 $2.2\times10^{-16}$ 一致を検証。

### Zenn 記事

- [double-slit-localized-source-fragility](https://zenn.dev/noriaki_kihara/articles/double-slit-localized-source-fragility)（論文1・論文2 をまとめて紹介）

### note 記事（一般向け、図2枚・正直な但し書き付き）

- 日本語: https://note.com/kiharanoriaki/n/n65be6bf06c9b （ドラフト：`note_doubleslit_ja.md`）
- English: https://note.com/kiharanoriaki/n/n701e9d57d7bb （ドラフト：`note_doubleslit_en.md`）

---

## 引用方針

- 論文2 の本文の自己参照は論文1 のみ（`isSupplementTo`）。外部引用は確立した教科書事項に限定：Born & Wolf『Principles of Optics』（二スリット遠方場・vCZ）、Mandel & Wolf『Optical Coherence and Quantum Optics』（vCZ・空間コヒーレンス）、高木『解析概論』（奇数倍音余弦和の閉形式）。
- 先行研究調査（2026-06-29）：本統合アプローチ（奇数倍音局在＋形保存干渉＋位置 push-forward の $\cos^2$ 読み出し）の先行論文は確認範囲で見当たらず。最も近い確立概念は空間コヒーレンス／vCZ（モード b＝可視度低下）で、本稿のモード a（押し出し）とは別観測量として §6.1 で区別。
