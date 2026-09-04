# 新プログラム N=40 自己対照（2026-09-04）

## 目的

新相互作用系プログラム `run_and_plot_N3_N40_mixedseed_20260903.py`（N=40 を含む唯一の新系プログラム）
をコピーし、N=40 について正本結果（2026-09-03 走行）と同一になるかを対照する。
新→旧への一因子変形実験（N=40 ライン）の出発点。

## ファイル

- `run_and_plot_N3_N40_mixedseed_20260903.py` — 正本の bit 同一コピー（diff で確認、無変更）
- `run_and_plot_N40_only_selfcontrol.py` — 上のコピーに対し変更は **2行のみ**（diff 全文を確認済み）:
  1. `OUT` を本フォルダ `results/` に変更（無変更だと正本 results/ を上書き破壊するため必須）
  2. `for N in range(3,41)` → `range(40,41)`（N=40 のみ実行。各 N は独立）

物理・入力・dtype・演算順序は一切変更なし。N=40 の入力は正本と同じ
`hm_mp_free_N3_N40_20260901/data/hm_N40/parent_v.npz` の `v`。

## 実行と合格ゲート（2026-09-04）

`python3 run_and_plot_N40_only_selfcontrol.py` → `done N 40` / `ALL DONE`（exit 0）

- **GATE PASS**: 6分母（38,39,40,41,42,124）の `hm_N40_den_*_states_500.npz` 全配列
  （Z 501×780・N・denominator・steps）が正本と **bit 一致**
- **GATE PASS**: timeseries CSV の N=40 行 3006 行、summary CSV の N=40 行 6 行が正本と完全一致

## 記録事項

- N=40 単独実行では matmul RuntimeWarning（divide by zero / overflow / invalid value）が
  N=40 処理中に表示された。正本の全系列走行では同種警告が N=8 で先に出て以降抑制されて
  いた（numpy は同一箇所の警告を1プロセス1回しか出さない）ため、N=40 でも同じ計算が
  起きていたことが今回可視化された。結果は bit 一致しており挙動差ではない。
- 図 `fig_Hperp_denominator_controls_with_124_N3_N40.png` は N=40 のパネルのみ描画される
  （legend 警告は空パネル由来、無害）。

## 追加読出し図（plot_complex_plane_N40_selfcontrol_v1.py、読み出しのみ）

Δτ=2π/40 の `results/hm_N40_den_40_states_500.npz` から。描画部は
`自発的分裂予備実験_v1/N40_state_readout_20260904/plot_complex_plane_N40_v1.py` と同一様式。

- `fig_N40_selfcontrol_complex_plane_step0.png` — step0（hm_mp_free 解析親）
- `fig_N40_selfcontrol_complex_plane_final.png` — 最大 step（τ=500）
- `fig_N40_selfcontrol_complex_plane_final_cluster_zoom.png` — 凝縮部（角クラスター）拡大

## 静的親差し替え実験（2026-09-04 実行・run_all_staticparent.sh）

初期データのみを 7月 make_parent 静的親
`自発的分裂予備実験_v1_N40対照実験系_20260904/largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz`
の `Z0` に差し替え（`run_N40_staticparent_v1.py`、自己対照版との差分は
OUT先・STATIC_PARENT 定数・z0 読込・図名・メタデータ名のみ。物理・分母系列・step数は不変）。
出力は `results_staticparent/`（既存 results/ は無傷）。

- **INPUT GATE PASS**: 全6分母 npz の Z[0] が静的 Z0 と bit 一致（check_staticparent_inputs_v1.py）
- 図: `results_staticparent/fig_Hperp_denominator_controls_with_124_N40_staticparent.png`
- 複素読出し図（plot_complex_plane_N40_staticparent_v1.py、Δτ=2π/40）:
  `fig_N40_staticparent_complex_plane_step0.png` / `_final.png` / `_final_cluster_zoom.png`

## 変形第1段: 振幅正規化の追加（2026-09-04 実行・run_all_staticparent_phaseonly.sh）

`run_N40_staticparent_phaseonly_v1.py` — 静的親版との差分は力学1行のみ:
`H=H_of(np.exp(1j*np.angle(z)),A)`（生成子を単位振幅化した波から構成。
旧 set_theta(np.angle(Z)) と同一の数学）。他は出力先・図名・メタデータ名のみ変更。
出力 `results_staticparent_phaseonly/`。

- **INPUT GATE PASS**（全6分母 Z[0] が静的 Z0 と bit 一致）
- 結果（振幅込み版との比較）: 振幅込みでは全分母が天井 ~1.16e-3 で安定だったのに対し、
  **位相のみでは全6分母が 0.05 を交差し（τ=198/148/49/24/18/4）、H⊥/H≈0.94〜0.9999 まで
  ほぼ完全に親平面から離脱**（den=38 のみ max 0.48・振動）。
- step1 ミスマッチは分母依存が強い: den=38 で 1.41e-9、39→124 で 8.5e-5→8.9e-3。
  den=38 は step100 の 2.1e-6 から step200 の 0.147 へ指数的成長区間を示す。
- **注意（解釈の限定）**: どの分母も f(1) は種スケール（3.6e-32）に留まらないため、
  これは δ 種のインフレーションではなくミスマッチ起動の不安定性。ただし
  「振幅正規化の有無だけで安定（~1e-3 天井）⇔ 強不安定（ほぼ完全離脱）が反転する」
  ことは単一因子の結果として確定。

## 変形第3段: 生成子を虚部のみに（2026-09-04 実行・run_all_staticparent_imK.sh）

`run_N40_staticparent_imK_v1.py` — 段2との差分は力学のみ:
`H = 1j*Im(H_of(exp(1j*angle(z))))` = i·K（K=sin(Δθ) 実反対称）。
exp(−iΔτ·iK)=exp(Δτ·K) の実直交回転。出力 `results_staticparent_imK/`。

- **INPUT GATE PASS**（全6分母 Z[0] bit一致）
- **f(1) ゲート合格**: step1 H⊥/H = 4.0〜9.1×10⁻²⁹（種スケール。段1/2 の 10⁻⁸〜10⁻²
  のミスマッチ注入が消滅）→ 親の相対平衡性が予測どおり回復
- **閉塞保存回復**: |Z·Z|/H は 1.2e-15 → 500 step 後も 2.8e-14（段1/2 は 10⁻³ まで成長していた）
- **緩和曲線（インフレーション）復活**: 全6分母で 10⁻²⁸ 付近から真っ直ぐな指数増幅。
  レートは分母依存（約4〜16 steps/decade、den=124 最速・den=39 最遅）で、
  飽和は 0.038〜0.10。7月正本（64 steps/decade・飽和0.166）との残差は時計（段4）の領分
- 図: `results_staticparent_imK/fig_inflation_N40_staticparent_imK.png`（単独大判）・
  `fig_Hperp_..._imK.png`（グリッド）・複素3枚 `fig_N40_staticparent_imK_complex_plane_*.png`

## 検証実験: 段2だけ除去（2026-09-05 実行・run_all_staticparent_ampimK.sh）

「段2（振幅正規化）は本当に必要か」の直接検証。段3との差分は力学1行のみ:
`H = 1j*Im(H_of(z))` = i·K_amp（K'_ef=|z_e||z_f|sin(θ_f−θ_e)、振幅込み・虚部のみ・実直交回転）。
出力 `results_staticparent_ampimK/`。事前予測: f(1) はミスマッチ水準へ跳び緩和曲線は出ない。

- **INPUT GATE PASS**（全6分母 Z[0] bit一致）
- **結果は予測どおり**: f(1)=1.8〜2.3×10⁻⁸（ミスマッチ注入。段3の 10⁻²⁹ と対照的）、
  全分母 0.05 未交差、~3.5×10⁻³ への緩慢な漂いのみ。**緩和曲線なし**
- **結論（2×2 が同一アーキテクチャ内で完結）**: 緩和曲線＝インフレーションには
  「振幅正規化（段2）」と「cos対称部の除去（段3）」の**両方が必要**。どちらか一方では
  親の平衡性が壊れ、種の指数増幅は観測不能

## 変形第4段: 段2+段3＋σ正規化時計＝段1の固定Δτを除去（2026-09-05 実行・run_all_staticparent_sigmaclock.sh）

`run_N40_staticparent_sigmaclock_v1.py` — 段3との差分は時計のみ:
位相因子を固定 exp(−i(2π/den)w) からスペクトル的 Cayley (1+iγw̃)/(1−iγw̃)
（w̃=w/σ_max、γ=tan(π/144)、σ_max は eigh の厳密値）に置換。旧写像の回転角を厳密再現。
den は力学に入らない。出力 `results_staticparent_sigmaclock/`。

- **INPUT GATE PASS**・**時計独立性確認**: 6分母の走行が完全同一（den 不使用の設計どおり）
- **f(1)=5.84×10⁻³⁰（種スケール）**、閉塞 2.8e-14 維持
- **緩和曲線が7月レートに収束**: 後半区間（step250→500）の傾き **65.8 steps/decade**
  （7月正本 64.0 steps/decade、差 ~3%）。f(500)=3.4e-20 で、外挿交差は τ~1700-2000
  （7月の crossing 2011 と整合）。STEPS=500 のため交差自体は未到達
- **結論（変形完了）**: 段2（振幅正規化）＋段3（cos除去）＋σ時計で、7月インフレーションの
  物理（種スケール平衡・指数レート・閉塞保存）が新アーキテクチャ上で定量再現される

## 次段（予定）

N=40 の初期データを、`次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904` の
make_parent（正本引数）で生成した**静的親ファイル**に差し替える。make_parent は動的には
呼ばず、v・g・Z0 を別々に保存し、Z0 が既存 `states_N00040_delta1e-15_seed0.npz` の Z0 と
bit 一致することを合格条件とする。
