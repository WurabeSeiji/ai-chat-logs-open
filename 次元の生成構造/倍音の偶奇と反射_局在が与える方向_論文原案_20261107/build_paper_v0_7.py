# -*- coding: utf-8 -*-
"""v0.6 に第 II 部（実装篇）を足して v0.7 を作る。コードブロックは exchange_core_v11.py の関数本文を inspect でそのまま埋める（写し間違いを防ぐ）。
数値は results/check_core_v11.json から入れる。"""
import inspect, io, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exchange_core_v11 as C
src = io.open("paper_parity_reflection_localization_direction_ja_v0.6.md", encoding="utf-8").read()
K = json.load(io.open("results/check_core_v11.json", encoding="utf-8"))
def code(*fns):
    return "```python\n" + "\n".join(inspect.getsource(f).rstrip() for f in fns) + "\n```\n"
def f(x, nd=3): return f"{x:.{nd}e}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}"
k2 = K["K2_offdiag"]; k4 = K["K4_handshake"]; selfb = K["K5_self"]; selfa = K["self_after"]
rows_k2 = "\n".join(f"| {e['ij']} | {e['Phi_read']:.3f} | {e['Phi_true']:.1f} | {e['Delta_read']:.3f} | {e['Delta_true']:.1f} | {e['c_read']:+.4f} | {e['c_true']:+.4f} | {e['u_err']:.1e} | {e['resid']:.1e} |" for e in k2)
rows_k4 = "\n".join(f"| {e['i']} | {e['delta_measured_deg']:.6f} | {e['delta_predicted_deg']:.6f} | {e['axis_err']:.1e} | " + "、".join(f"{j}: {v:.4f}" for j, v in e['view_angle_terms |X_ij|/R_i^2'].items()) + " |" for e in k4)
rows_self = "\n".join(f"| {i} | {selfb[i]['R2']:.5f} | {selfb[i]['L']:.5f} → {selfa[i]['L']:.5f} | {selfb[i]['N_eff']:.2f} → {selfa[i]['N_eff']:.2f} | {selfb[i]['eta']:+.4f} → {selfa[i]['eta']:+.4f} |" for i in ["1", "2", "3", "4"])
part2 = """
---

# 第 II 部　実装篇：理論の量と数値実験のコードの対応

本篇は、第 I 部の読出し規則と、7/10・7/13・7/18 の二チャネル交換写像 [S2][S3][S4][S5] を、四元数 4 波の数値実験にどう写すかを、式とコードを並べて示す。コードブロックは `exchange_core_v11.py` の関数本文をそのまま載せている。実験の走行は本篇に含めない。確認（K1〜K5）は `check_core_v11.py` で、4 波・16 の関係に対して行った。

## 9. 対応表

11/7 に決めた計画の七項目と、本篇の節・関数の対応である。

| 計画 | 本篇 | 関数 |
|---|---|---|
| 作用は 7/13 プログラムの忠実コピーから、R 掃引の表の再現を対照に | 14 節。二チャネル交換写像を射影形で書き、4 波の同時作用素に拡張する。原本の表の再現は別フォルダの `control_original_v1.py`（5 行一致） | `generator`、`simultaneous_unitary` |
| 交換は自己を含む全員対全員で右乗 | 11・14 節。関係は N × N ＝ 16 の行列。生成子の対角が自己の項、U ＝ exp(A) を一度に掛ける。順序はない | `relation_matrix`、`generator` |
| Δ_F,ij は対の内部時計の差を倍音ごとの交差項の切片から読み R_ij ＝ sin²(Φ_ij/2) | 13 節。縞の位相 ϑ_m ＝ Φ − mΔ の切片 | `read_pair` |
| 軸 n_ij は倍音にわたるベクトル部の係数列から | 13 節。Vec_m の列が張る面から û ∝ n_i ＋ n_j、ŵ ∝ n_i × n_j、軸の角 c を読む。波ごとの軸は使わない | `pair_frame`、`orient_axes` |
| 読出しは観測者の交換前後の差とグラム行列の固有値 | 15 節と 12 節 | `rotation_measured`、`handshake_prediction`、`gram_eigs` |
| 各パルスの w_e・w_o・η と可視度の照合 | 12 節。自己関係から読む。可視度 4 w_e w_o は対蹠対（第 3 節）の量で、別位置の 4 波には第 15 節の握手の式が対応する | `read_self` |
| 実装前に照合表を出す | 本篇そのもの | — |

## 10. 状態と初期値

状態は閉じた弦 ψ ∈ [0, 2π) の格子 NG ＝ 2048 点の上の四元数波形 q_i(ψ) ＝ W ＋ V_i i ＋ V_j j ＋ V_k k、N ＝ 4 である。初期値は index_v9q4 のパルス

$$
q_i(\\psi)=A_i\\Big[c_0 e^{n_i\\varphi_{0i}+\\sum_m c_m e^{n_i(\\varphi_«0»+m(\\psi-\\theta_i))}\\Big],\\qquad c_m=\\frac«1»«2»\\frac«3»«4»,\\ c_0=\\frac«5»«6»
$$

で、θ ＝ (97, 194, 58, 300)°、位相幅 φ ＝ 1°、M ＝ 100、奇偶、A ＝ 1、n ＝ i, j, k, (1, 1, 1)/√3、内部位相 φ₀ ＝ (0, 60, 120, 180)° を確認に使った。共通の格子 ψ と倍音番号 m は、第 I 部の議論のとおり名前である。本篇では初期値として置く（10/20 の規則：塊は初期値で与える）。作用と読出しは名前を使わない。

«7»

## 11. 十六の関係

関係は N × N の行列で、i ＝ j の自己関係を含む。左だけの無名性（全体を左から単位四元数 g で掛ける）のもとで 4 成分とも不変なのは左共役の積なので、

$$
Y_{ij,m}=\\langle \\bar q_{i,m}\\,q_{j,m}\\rangle_\\psi=\\tfrac12\\big(\\bar a_{i,m}a_{j,m}+\\bar b_{i,m}b_{j,m}\\big),\\qquad q_{i,m}(\\psi)=a_{i,m}\\cos m\\psi+b_{i,m}\\sin m\\psi
$$

を関係とする。a_m、b_m は各成分の実フーリエ係数（四元数）、m ＝ 0 は直流で ψ 平均に 1/2 が付かない。Y_ji ＝ conj(Y_ij) で、16 本のうち独立な実数は対角 4 本と非対角 6 本 × 4 成分である。確認 K1：エルミート性の最大誤差 «8»。

«9»

## 12. 自己関係から読むもの

対角 Y_ii,m のスカラー部は倍音 m のパワーで、その和が R_i²、分布が N_eff。偶奇の重みは対蹠像との重なり η ＝ ⟨q, Pq⟩/‖q‖² から w_e ＝ (1 ＋ η)/2、w_o ＝ (1 − η)/2（第 2 節）。局在 L ＝ Σ p²。

«10»

## 13. 非対角の関係から読むもの：対の枠と縞

### 13.1 式

軸が違う二波の倍音 m の左共役の積を ψ で平均すると、α ＝ φ₀ − mθ、ϑ_m ＝ α_{i,m} − α_{j,m} ＝ Φ − mΔ、c ＝ n_i·n_j として

$$
Y_m=\\frac{c_ic_j}{2}\\Big[(1+c)\\cos\\vartheta_m\;-\;\\sin\\vartheta_m\\,(n_i+n_j)\;-\;\\cos\\vartheta_m\\,(n_i\\times n_j)\\Big]
$$

である。第 4.2 節の三つの係数列のうち n_i − n_j の方向は ψ 平均で消え、残る二方向 (n_i ＋ n_j) と (n_i × n_j) が対の枠になる。ここから読めるものは次である。

- 面：Vec_m の列は (n_i ＋ n_j) と (n_i × n_j) の張る面にある。平面性（上位 2 固有値の割合）が 1 なら式どおり
- ŵ ∝ n_i × n_j：Vec_m·ŵ ＝ −tan(γ/2) Sc_m が全 m で成り立つ唯一の方向（γ は軸の角、cos γ ＝ c）。最小二乗 w ＝ argmin Σ_m (Vec_m·w ＋ Sc_m)² の向きが ŵ、|w| ＝ 1/tan(γ/2) から c が出る
- û ∝ n_i ＋ n_j：面の中で ŵ に直交する方向。符号は線としてしか決まらない
- 縞：tan ϑ_m ＝ −s (Vec_m·û)/Sc_m、s ＝ √((1 ＋ c)/2)。ϑ_m を m にわたって直線に当て、切片 Φ_ij、傾き −Δ_ij

û の符号を反転すると Φ と Δ の符号が反転する。対の向きと時計の向きは一組で、二波だけでは決まらない（第 4.2 節の二価性）。4 波では全対の û ûᵀ の和（符号によらない）の主固有ベクトル r̂ を基準に û·r̂ ＞ 0 と取れば、順序 i↔j で同じ û になり、Φ_ji ＝ −Φ_ij、Δ_ji ＝ −Δ_ij が自動的に成り立つ。r̂ 自体の符号だけが残る。これは 4 体で決まらない唯一の離散量、全体の鏡像である（第 4.3 節）。

### 13.2 コード

«11»

### 13.3 確認 K2：12 本の順序つき非対角関係

初期値（既知の θ、φ₀、n）から読んだ値と真値。û は r̂ で向きを揃えた後、真値 unit(n_i ＋ n_j) も同じ向きに揃えて比べた。

| ij | Φ 読み（度） | Φ 真 | Δ 読み（度） | Δ 真 | c 読み | c 真 | \\|û − 真\\| | 当てはめ残差 |
|---|---|---|---|---|---|---|---|---|
«12»

順序の一致：max |Φ_ji ＋ Φ_ij| ＝ «13» 度、max |û_ji − û_ij| ＝ «14»。基準 r̂ ＝ ({', '.join(f'«15»' for v in K['orientation_reference_r_hat'])})。

## 14. 作用：十六の関係から生成子、同時ユニタリ

### 14.1 式

7/13 の二チャネル交換写像は、射影形 [S5] で U ＝ P_S ＋ e^{iφ_a} P_D、φ_a ＝ π ＋ Δ_F である（Δ_F ＝ π で U ＝ I、Δ_F ＝ 0 で入替）。P_D は射影なので U ＝ exp(φ_a P_D) と書ける。四元数では位相 e^{iφ} を単位四元数 e^{ûφ} の右乗に置き換える。右乗は左だけの無名性（全体を左から g 倍）と可換で、右乗の 4 × 4 実表現を R(h) とすると、対 (i, j) の生成子は

$$
G_{ij}=\\phi_{ij}\;P_D^{(ij)}\\otimes R(\\hat u_{ij}),\\qquad P_D^{(ij)}=\\tfrac12(|i\\rangle-|j\\rangle)(\\langle i|-\\langle j|),\\qquad \\phi_{ij}=\\pi+|\\Phi_{ij}|
$$

である。P_D^{(ij)} の対角 (i, i)、(j, j) が自己の項、非対角 (i, j)、(j, i) が相手の項で、4 × 4 の関係行列のすべての成分が生成子に入る。全員対全員の同時作用素は生成子の和の指数

$$
U=\\exp\\Big(\\sum_{i<j}G_{ij}\\Big)
$$

で、順序に依らない。G_ij は（対称）⊗（反対称）なので反対称実行列、したがって U は直交行列で全ノルムを保つ。対が一組だけなら U ＝ P_S ⊗ I ＋ P_D ⊗ R(e^{ûφ}) に戻る。

Δ_F,ij に |Φ_ij| を使うのは次の理由による。順序 i↔j で û は同じ、Φ は符号が反転する（13 節）。生成子が反エルミート（A_ji ＝ A_ij、û 型の成分）であるためには φ_ij が i↔j で同じでなければならず、符号つきの Φ では π ＋ Φ と π − Φ が一致しない。7/13 でも Δ_F は対のパラメータで順序を持たない。R ＝ sin²(Δ_F/2) は偶関数なので R は変わらない。Φ の符号（どちらの時計が進んでいるか）は読出しとして残り、作用の大きさには入らない。

### 14.2 コード

«16»

### 14.3 確認 K3

G の反対称性 «17»、UᵀU − I の最大 «18»、U のブロックが右乗であることの最大誤差 «19»、一回の作用での全ノルムの相対ずれ «20»。

## 15. 位相の回転に視野角が入ること

U の (i, j) ブロックは右乗 R(h_ij) になる（右乗の和と積は右乗）。したがって作用後の波は q_i′ ＝ Σ_j q_j h_ij で、観測者 i の位相の回転は

$$
\\langle q_i,q_i'\\rangle=\\int\\bar q_i\\,q_i'=\\sum_j\\Big(\\int\\bar q_i\\,q_j\\Big)h_{ij}=\\sum_j Y_{ij}\\,h_{ij}
$$

である。Y_ii ＝ R_i² は実数、j ≠ i の項は |Y_ij|/R_i² の大きさで入る。|Y_ij|/R_i² は相手の大きさと重なりの積、すなわち視野角の大きさで、その四元数の向きが視野角の向きである。回転角は arg、軸は Vec の向き。この式は 16 の関係と作用素だけで書けるので、「距離の読みが状態の相互作用にどう入るか」がここで読める。

«21»

確認 K4：一回の作用の前後で測った回転と、16 の関係からの予測。

| 観測者 i | δ_i 実測（度） | δ_i 予測（度） | 軸の差 | 視野角項 \\|Y_ij\\|/R_i² |
|---|---|---|---|---|
«22»

最大差 «23» 度。

## 16. 一回の作用の帰結（この設定で理論どおりに作ると何が起こるか）

初期値 φ₀ ＝ (0, 60, 120, 180)° では |Φ_ij| ＝ 60°、120°、180° で R_ij ＝ 0.25、0.75、1.00 になる。別の位置にあるパルスでも、作用の大きさは |Φ_ij| だけで決まり、視野角の大きさ |Y_ij|/R_i²（ここでは 0.001〜0.01）は位相の回転 δ_i にしか入らない。その結果、一回の作用で各波は相手の位置に山を持つ重ね合わせになり、局在が落ちる。

| 波 | R² | L（前 → 後） | N_eff（前 → 後） | η（前 → 後） |
|---|---|---|---|---|
«24»

グラム行列の固有値は前 ({', '.join(f'«25»' for v in K['gram_eigs_before'])})、後 ({', '.join(f'«26»' for v in K['gram_eigs_after'])})。これは計画どおりに作ったときの帰結であり、隠さない。

## 17. 理論が決めていない点

次は第 I 部にも 7/10・7/13 にも書かれていない選択で、本篇では上記のとおり置いたが、著者の判断を要する。

1. 作用の位相 Δ_F,ij に、縞の切片 Φ_ij（内部時計の差）を使うか、倍音ごとの縞の位相 ϑ_m ＝ Φ − mΔ（位置差を含む）を使うか。7/13 は二波が同位置（Δ ＝ 0）なので両者は同じだった。別位置の波では違う
2. 作用の大きさが |Φ| だけで決まり、視野角の大きさ（重なり）が位相の回転にしか入らないこと。16 節の帰結はここから出る
3. Δ_F に |Φ| を使うこと（順序を持たない対の量にするため）。符号つきで使うには、対の向きの約束が要る
4. 作用を波形全体に掛けるか、共有できる帯域だけに掛けるか（第 5.1 節）
5. 運動学 θ ← θ ＋ φ を置かないこと（11/7 の決定）。並進は位相の回転の累積 Θ_i n_i として読む
6. 全体の鏡像 r̂ の符号
7. 共通の格子 ψ と倍音番号 m を初期値として置くこと。倍音が再帰の回数として生じる仕組みは本篇にない

## 18. 確認の一覧とファイル

| 確認 | 内容 | 結果 |
|---|---|---|
| K1 | 関係行列のエルミート性 Y_ji ＝ conj(Y_ij) | «27» |
| K2 | 12 本の順序つき非対角関係の Φ・Δ・c・û・ŵ の読みと順序の一致 | Φ «28» 度、Δ «29» 度、c «30»、û «31»、ŵ «32» |
| K3 | 生成子の反対称、U の直交、ブロックが右乗、ノルム保存 | «33»、«34»、«35»、«36» |
| K4 | 位相の回転の実測と 16 の関係からの予測 | «37» 度 |
| K5 | 自己関係の読み（R²、L、N_eff、偶奇） | 16 節の表 |

`exchange_core_v11.py`（中核関数、本篇のコードブロックの出典）、`check_core_v11.py`（K1〜K5）、`results/check_core_v11.json`（値）。
"""
def _tok(src, fn, fmt=None):
    try:
        v = fn(); return format(v, fmt) if fmt else str(v)
    except (NameError, TypeError, ValueError, AttributeError):
        return '{' + src + (':' + fmt if fmt else '') + '}'
_vals = [
    '{0i}',
    _tok('2', lambda: 2),
    '{\\\\pi}',
    '{\\\\sin(m\\\\pi\\\\varphi/360)}',
    _tok('m', lambda: m),
    '{\\\\varphi}',
    _tok('360', lambda: 360),
    _tok('code(C.coefficients, C.make_pulse)', lambda: code(C.coefficients, C.make_pulse)),
    _tok("K['K1_hermitian_max_err']", lambda: K['K1_hermitian_max_err'], '.1e'),
    _tok('code(C.harmonic_coeffs, C.relation_matrix)', lambda: code(C.harmonic_coeffs, C.relation_matrix)),
    _tok('code(C.read_self)', lambda: code(C.read_self)),
    _tok('code(C.pair_frame, C.orient_axes, C.read_pair)', lambda: code(C.pair_frame, C.orient_axes, C.read_pair)),
    _tok('rows_k2', lambda: rows_k2),
    _tok("K['K2_order_consistency']['max|Phi_ji + Phi_ij| deg']", lambda: K['K2_order_consistency']['max|Phi_ji + Phi_ij| deg'], '.1e'),
    _tok("K['K2_order_consistency']['max|u_ji - u_ij|']", lambda: K['K2_order_consistency']['max|u_ji - u_ij|'], '.1f'),
    _tok('v', lambda: v, '.4f'),
    _tok('code(C.Rmat, C.generator, C.simultaneous_unitary, C.apply_unitary)', lambda: code(C.Rmat, C.generator, C.simultaneous_unitary, C.apply_unitary)),
    _tok("K['K3_antisymmetric_G']", lambda: K['K3_antisymmetric_G'], '.1e'),
    _tok("K['K3_orthogonal_U']", lambda: K['K3_orthogonal_U'], '.1e'),
    _tok("K['K3_blocks_are_right_mult_err']", lambda: K['K3_blocks_are_right_mult_err'], '.1e'),
    _tok("K['K3_norm_rel_drift']", lambda: K['K3_norm_rel_drift'], '.1e'),
    _tok('code(C.blocks_as_quaternions, C.handshake_prediction, C.rotation_measured)', lambda: code(C.blocks_as_quaternions, C.handshake_prediction, C.rotation_measured)),
    _tok('rows_k4', lambda: rows_k4),
    _tok("K['K4_max_err_deg']", lambda: K['K4_max_err_deg'], '.1e'),
    _tok('rows_self', lambda: rows_self),
    _tok('v', lambda: v, '.4f'),
    _tok('v', lambda: v, '.4f'),
    _tok("K['K1_hermitian_max_err']", lambda: K['K1_hermitian_max_err'], '.1e'),
    _tok("K['K2_max_err_Phi_deg']", lambda: K['K2_max_err_Phi_deg'], '.1e'),
    _tok("K['K2_max_err_Delta_deg']", lambda: K['K2_max_err_Delta_deg'], '.1e'),
    _tok("K['K2_max_err_c']", lambda: K['K2_max_err_c'], '.1e'),
    _tok("K['K2_max_err_u']", lambda: K['K2_max_err_u'], '.1e'),
    _tok("K['K2_max_err_w']", lambda: K['K2_max_err_w'], '.1e'),
    _tok("K['K3_antisymmetric_G']", lambda: K['K3_antisymmetric_G'], '.1e'),
    _tok("K['K3_orthogonal_U']", lambda: K['K3_orthogonal_U'], '.1e'),
    _tok("K['K3_blocks_are_right_mult_err']", lambda: K['K3_blocks_are_right_mult_err'], '.1e'),
    _tok("K['K3_norm_rel_drift']", lambda: K['K3_norm_rel_drift'], '.1e'),
    _tok("K['K4_max_err_deg']", lambda: K['K4_max_err_deg'], '.1e'),
]
for _i, _v in enumerate(_vals): part2 = part2.replace('«%d»' % _i, _v)
# 第 II 部は参考文献の前に入れる
marker = "\n## 参考文献\n"
assert src.count(marker) == 1
out = src.replace(marker, part2 + "\n## 参考文献\n")
out = out.replace("**版：** 草稿 v0.6（§7 に代数の帰属の向きを一文追加）", "**版：** 草稿 v0.7（第 II 部 実装篇を追加：理論の量と数値実験のコードの対応、4 波・16 関係での確認 K1〜K5）")
# テンプレート処理で崩れた LaTeX の修復（10 節の式）：該当行を丸ごと正しい式に置き換える
import re as _re
_good = r"q_i(\psi)=A_i\Big[c_0 e^{n_i\varphi_{0i}}+\sum_m c_m e^{n_i(\varphi_{0i}+m(\psi-\theta_i))}\Big],\qquad c_m=\frac{2}{\pi}\frac{\sin(m\pi\varphi/360)}{m},\ c_0=\frac{\varphi}{360}"
_lines = out.split("\n"); _n = 0
for _k, _l in enumerate(_lines):
    if _l.startswith("q_i(\\psi)=A_i\\Big[c_0"): _lines[_k] = _good; _n += 1
assert _n == 1, _n
out = "\n".join(_lines)
assert "\\\\pi" not in out, "二重バックスラッシュが残っている"
# 入れ子の式（ベクトルの整形）が文字どおり残るので値に置き換える
_rhat = ", ".join("%.4f" % v for v in K["orientation_reference_r_hat"])
_gb = ", ".join("%.4f" % v for v in K["gram_eigs_before"]); _ga = ", ".join("%.4f" % v for v in K["gram_eigs_after"])
out = out.replace("({', '.join(f'{v:.4f}' for v in K['orientation_reference_r_hat'])})", "(" + _rhat + ")")
out = out.replace("({', '.join(f'{v:.4f}' for v in K['gram_eigs_before'])})", "(" + _gb + ")").replace("({', '.join(f'{v:.4f}' for v in K['gram_eigs_after'])})", "(" + _ga + ")")
assert "'.join(" not in out, "未展開の式が残っている"
io.open("paper_parity_reflection_localization_direction_ja_v0.7.md", "w", encoding="utf-8").write(out); print("v0.7 written", len(out))
