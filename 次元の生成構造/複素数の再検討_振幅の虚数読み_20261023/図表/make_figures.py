# -*- coding: utf-8 -*-
"""複素数の再検討：振幅はどこに住むか ― 論文用の静的な図（PNG・SVG・PDF）を作る。

図 1  指数平面 w ＝ a ＋ iθ と、i を掛ける 90° 回転で振幅 a と位相 θ が入れ替わること
図 2  振幅の置き場所：動径（|z| ＝ R と一致して区別できない）と、円周と直交する軸 Q（円柱 S¹ × ℝ_Q）
図 3  w ＝ log u ＝ φ ＋ iθ の二つの断面：|cosh(w/2)|²、|sinh(w/2)|²、その和 cosh φ ＝ γ（複素角 ω ＝ θ − iφ では cos・sin）
図 4  結合係数 |k| ＝ |tanh(w/2)| と、自己双対点 |k| ＝ 1（θ ＝ ±π/2、φ によらない）。分母 |1＋u|² ＝ 2e^φ(cosh φ ＋ cos θ)
図 5  四元数：振幅軸 m と位相軸 n の回転の交換子が第三の軸 m × n を生む（e^{ma}e^{nθ} − e^{nθ}e^{ma} ＝ 2 sin a sin θ (m × n)）
図 6  帰属の台帳：ℂ にすでにあるもの／ℍ が足すもの／𝕆 が足すもの

点はすべて式から直接計算する。乱数は使わない。
実行:  python make_figures.py [en]     （en で英語のラベル。計算と描画は同じで文字だけが違う）
出力:  figures/figN_*.png / .svg / .pdf（英語版は figN_*_en）
依存:  numpy, matplotlib
記法の宣言（検討ルール）：u ＝ e^w、w ＝ φ ＋ iθ（φ ＝ 対数振幅＝ラピディティ、θ ＝ 位相角）。k ＝ (1−u)/(1＋u) ＝ −tanh(w/2)、
  分母は合成 1 ＋ u（読み出しの最後に一度だけ割る）。四元数の積は左右を区別して書く。
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

EN = len(sys.argv) > 1 and sys.argv[1] == "en"
L = (lambda ja, en: en if EN else ja)
SUF = "_en" if EN else ""
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.unicode_minus"] = False
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

C_A, C_P, C_T, C_INK, C_MUTED, C_GRID = "#c4470f", "#2a56c6", "#0e7a6e", "#18212c", "#596577", "#c9d0d8"


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{name}{SUF}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", name + SUF)


# ---------------------------------------------------------------- 四元数の小道具（[w, x, y, z]）
def qmul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])


def qexp_pure(n, t):
    """純虚単位 n（3 ベクトル）と角 t で e^{n t} ＝ cos t ＋ n sin t"""
    return np.concatenate([[np.cos(t)], np.sin(t) * np.asarray(n, float)])


# ================================================================ 図 1
def fig1():
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    ax.set_aspect("equal")
    ax.axhline(0, color=C_GRID, lw=1); ax.axvline(0, color=C_GRID, lw=1)
    a, th = 1.6, 0.9
    w = complex(a, th)
    iw = 1j * w
    for z, col, lab in ((w, C_A, L("w ＝ a ＋ iθ", "w = a + iθ")), (iw, C_P, L("i·w ＝ −θ ＋ ia", "i·w = −θ + ia"))):
        ax.add_patch(FancyArrowPatch((0, 0), (z.real, z.imag), arrowstyle="-|>", mutation_scale=16, lw=2, color=col))
        ax.text(z.real + 0.08, z.imag + 0.08, lab, color=col, fontsize=11)
    # 90° の弧
    r = 0.9
    ang = np.linspace(np.angle(w), np.angle(iw), 60)
    ax.plot(r*np.cos(ang), r*np.sin(ang), color=C_T, lw=1.5)
    ax.text(r*np.cos(ang[30]) + 0.05, r*np.sin(ang[30]) + 0.05, "90°", color=C_T, fontsize=11)
    # 成分の目印
    ax.plot([a, a], [0, th], ls=":", color=C_A); ax.plot([0, a], [th, th], ls=":", color=C_A)
    ax.text(a, -0.22, L("a（対数振幅）", "a (log-amplitude)"), ha="center", color=C_A, fontsize=10)
    ax.text(a/2, th + 0.07, L("θ（位相）", "θ (phase)"), ha="center", va="bottom", color=C_A, fontsize=10)
    ax.plot([-th, -th], [0, a], ls=":", color=C_P); ax.plot([0, -th], [a, a], ls=":", color=C_P)
    ax.text(-th, -0.22, "−θ", ha="center", color=C_P, fontsize=10)
    ax.text(-th/2, a - 0.07, "a", ha="center", va="top", color=C_P, fontsize=10)
    ax.set_xlim(-1.6, 2.4); ax.set_ylim(-0.6, 2.4)
    ax.set_xlabel(L("実部（振幅の指数）", "real part (exponent of amplitude)"))
    ax.set_ylabel(L("虚部（位相）", "imaginary part (phase)"))
    ax.set_title(L("指数平面：i を掛けると振幅と位相が入れ替わる（双対）", "Exponent plane: multiplying by i swaps amplitude and phase (duality)"), fontsize=11)
    txt = L("同じ回転の現れ：\n数と位相の共役　時間と振動数\nWick 回転（減衰 ↔ 振動）\nLaplace の s ＝ σ ＋ iω",
            "The same rotation appears as:\nnumber–phase conjugacy, time–frequency,\nWick rotation (decay ↔ oscillation),\nLaplace s = σ + iω")
    ax.text(0.98, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=C_MUTED,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C_GRID))
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    save(fig, "fig1_exponent_plane_duality")


# ================================================================ 図 2
def fig2():
    fig = plt.figure(figsize=(11, 4.8))
    # 左：平面の円。振幅＝動径だと |z| ＝ R と一致する
    ax = fig.add_subplot(1, 2, 1); ax.set_aspect("equal")
    R, th = 1.0, 0.8
    ax.add_patch(Circle((0, 0), R, fill=False, color=C_INK, lw=1.5))
    z = R*np.exp(1j*th)
    ax.add_patch(FancyArrowPatch((0, 0), (z.real, z.imag), arrowstyle="-|>", mutation_scale=14, lw=2, color=C_A))
    ax.text(z.real*0.5 - 0.02, z.imag*0.5 + 0.1, L("|z| ＝ R ？", "|z| = R ?"), color=C_A, fontsize=11)
    ax.plot([0, R], [0, 0], color=C_T, lw=1.5); ax.text(R*0.5, -0.14, L("曲率半径 R", "curvature radius R"), color=C_T, ha="center", fontsize=10)
    arc = np.linspace(0, th, 40); ax.plot(0.35*np.cos(arc), 0.35*np.sin(arc), color=C_P, lw=1.2); ax.text(0.42, 0.12, "θ", color=C_P, fontsize=11)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(L("平面の円：振幅を動径に置くと\n曲率半径 R と同じ数になり、区別できない", "Circle in the plane: amplitude on the radius\ncoincides with R; the two cannot be told apart"), fontsize=10.5)
    # 右：円柱 S¹ × ℝ_Q。振幅は円周とも R とも直交する軸 Q に
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    t = np.linspace(0, 2*np.pi, 200)
    for q0 in (0.0,):
        ax.plot(R*np.cos(t), R*np.sin(t), q0*np.ones_like(t), color=C_INK, lw=1.5)
    for q0 in (-0.9, 0.9):
        ax.plot(R*np.cos(t), R*np.sin(t), q0*np.ones_like(t), color=C_GRID, lw=1)
    for tt in np.linspace(0, 2*np.pi, 12, endpoint=False):
        ax.plot([R*np.cos(tt)]*2, [R*np.sin(tt)]*2, [-0.9, 0.9], color=C_GRID, lw=0.6)
    A = 0.7
    p0 = (R*np.cos(th), R*np.sin(th), 0.0); p1 = (R*np.cos(th), R*np.sin(th), A)
    ax.plot([0, p0[0]], [0, p0[1]], [0, 0], color=C_T, lw=1.5)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], color=C_A, lw=2.5)
    ax.scatter(*p1, color=C_A, s=40)
    ax.plot([0, 0], [0, 0], [-1.0, 1.1], color=C_P, lw=1.2, ls="--")
    ax.text(0, 0, 1.18, L("Q 軸（読めない方向）", "Q axis (unreadable)"), color=C_P, fontsize=9)
    ax.text(p1[0]+0.05, p1[1], p1[2]+0.05, L("振幅 A", "amplitude A"), color=C_A, fontsize=10)
    ax.text(0.45*np.cos(th), 0.45*np.sin(th), -0.12, "R", color=C_T, fontsize=10)
    ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_zlim(-1.0, 1.2)
    ax.set_box_aspect((1, 1, 0.9)); ax.view_init(elev=22, azim=-60)
    ax.set_axis_off()
    ax.set_title(L("円柱 S¹ × ℝ_Q：振幅は円周とも R とも直交する\n読めない方向 Q に住む（Q 軸の由来）", "Cylinder S¹ × ℝ_Q: the amplitude lives on the unreadable\naxis Q orthogonal to both the circle and R (origin of the Q axis)"), fontsize=10.5)
    save(fig, "fig2_amplitude_radius_vs_Q_axis")


# ================================================================ 図 3
def complex_half_angle(theta, phi):
    """w ＝ log u ＝ φ ＋ iθ。|cosh(w/2)|², |sinh(w/2)|² を式どおり。
    複素角 ω ＝ θ − iφ（u ＝ e^{iω}）で書けば cosh(w/2) ＝ cos(ω/2)、sinh(w/2) ＝ i sin(ω/2)。"""
    w = phi + 1j*theta
    c = np.cosh(w/2); s = np.sinh(w/2)
    return np.abs(c)**2, np.abs(s)**2


def fig3():
    th = np.linspace(-np.pi, np.pi, 361); ph = np.linspace(-2.0, 2.0, 241)
    TH, PH = np.meshgrid(th, ph)
    c2, s2 = complex_half_angle(TH, PH)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, Z, ttl in ((axs[0], c2, "|cosh(w/2)|² ＝ |cos(ω/2)|²"), (axs[1], s2, "|sinh(w/2)|² ＝ |sin(ω/2)|²")):
        cs = ax.contourf(TH, PH, Z, levels=np.linspace(0, 4, 17), cmap="viridis", extend="max")
        ax.contour(TH, PH, Z, levels=[0.5, 1, 2, 3], colors="white", linewidths=0.7)
        ax.axhline(0, color=C_A, lw=1.5, ls="--")
        ax.set_xlabel(L("θ（位相）", "θ (phase)")); ax.set_ylabel(L("φ（ラピディティ＝対数振幅）", "φ (rapidity = log-amplitude)"))
        ax.set_title(ttl, fontsize=11)
        ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi]); ax.set_xticklabels(["−π", "−π/2", "0", "π/2", "π"])
        fig.colorbar(cs, ax=ax, shrink=0.85)
    axs[0].text(-np.pi+0.1, 0.08, L("φ ＝ 0 の断面：cos²(θ/2)（ボルン則の形）", "section φ = 0: cos²(θ/2) (the form of the Born rule)"), color=C_A, fontsize=9)
    # 右：和 ＝ cosh φ ＝ γ
    ax = axs[2]
    ph1 = np.linspace(-2, 2, 241)
    for th0, col in ((0.0, C_P), (np.pi/2, C_T), (np.pi, C_A)):
        c, s = complex_half_angle(th0, ph1)
        ax.plot(ph1, c + s, color=col, lw=2, label=L(f"θ ＝ {th0/np.pi:.1f}π の和", f"sum at θ = {th0/np.pi:.1f}π"))
    ax.plot(ph1, np.cosh(ph1), "k:", lw=1.5, label="cosh φ ＝ γ")
    ax.axhline(1, color=C_GRID, lw=1)
    ax.set_xlabel("φ"); ax.set_ylabel(L("|cosh(w/2)|² ＋ |sinh(w/2)|²", "|cosh(w/2)|² + |sinh(w/2)|²"))
    ax.set_title(L("和は θ によらず cosh φ。φ ＝ 0 でだけ 1", "The sum is cosh φ for every θ; equal to 1 only at φ = 0"), fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper center")
    fig.suptitle(L("w ＝ log u ＝ φ ＋ iθ（複素角 ω ＝ θ − iφ、u ＝ e^{iω}）の二つの断面：ボルン則の正規化は cosh φ ＝ γ で割ること",
                   "Two sections of w = log u = φ + iθ (complex angle ω = θ − iφ, u = e^{iω}): normalizing to the Born rule divides by cosh φ = γ"), fontsize=12)
    fig.tight_layout()
    save(fig, "fig3_complex_angle_sections_gamma")


# ================================================================ 図 4
def fig4():
    th = np.linspace(-np.pi, np.pi, 361); ph = np.linspace(-2.0, 2.0, 241)
    TH, PH = np.meshgrid(th, ph)
    u = np.exp(PH + 1j*TH)
    k = (1 - u) / (1 + u)                      # 分母 ＝ 合成 1 ＋ u（最後に一度だけ割る）
    ak = np.abs(k); ak = np.where(ak > 0, ak, np.nan)   # u ＝ 1（同位相）で k ＝ 0、log の極を外す
    den = np.abs(1 + u)**2
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axs[0]
    cs = ax.contourf(TH, PH, np.log10(ak), levels=np.linspace(-2, 2, 17), cmap="RdBu_r", extend="both")
    ax.contour(TH, PH, ak, levels=[1.0], colors=[C_T], linewidths=2.5)
    ax.axhline(0, color=C_A, lw=1.2, ls="--")
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi]); ax.set_xticklabels(["−π", "−π/2", "0", "π/2", "π"])
    ax.set_xlabel(L("θ（相対位相の角）", "θ (angle of relative phase)")); ax.set_ylabel(L("φ（振幅の不一致）", "φ (amplitude mismatch)"))
    ax.set_title(L("|k| ＝ |tanh(w/2)|（log₁₀）。緑：自己双対点 |k| ＝ 1", "|k| = |tanh(w/2)| (log₁₀). Green: self-dual locus |k| = 1"), fontsize=10.5)
    ax.text(-np.pi + 0.12, 1.85, L("緑の線 θ ＝ ±π/2：\n|k| ＝ 1 は φ によらない\n＝ 光円錐・結合の閾値", "green lines θ = ±π/2:\n|k| = 1 for every φ\n= light cone, binding threshold"), color=C_T, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_T, alpha=0.9))
    fig.colorbar(cs, ax=ax, shrink=0.85, label="log₁₀|k|")
    ax = axs[1]
    for ph0, col in ((0.0, C_P), (0.5, C_T), (1.0, C_A)):
        d = np.abs(1 + np.exp(ph0 + 1j*th))**2
        ax.plot(th, d, color=col, lw=2, label=f"|1＋u|²,  φ ＝ {ph0:.1f}")
        ax.plot(th, 2*np.exp(ph0)*(np.cosh(ph0) + np.cos(th)), "k:", lw=1)
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi]); ax.set_xticklabels(["−π", "−π/2", "0", "π/2", "π"])
    ax.set_xlabel("θ"); ax.set_ylabel("|1 + u|²")
    ax.set_title(L("分母 |1＋u|² ＝ 2e^φ(cosh φ ＋ cos θ)（点線）\n＝ 4e^φ |cosh(w/2)|²：ボルン則の正規化と同じ割り算", "Denominator |1+u|² = 2e^φ(cosh φ + cos θ) (dotted)\n= 4e^φ |cosh(w/2)|²: the same division as the Born normalization"), fontsize=10.5)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig4_k_selfdual_denominator")


# ================================================================ 図 5
def fig5():
    m = np.array([1.0, 0, 0]); n = np.array([0, 1.0, 0]); mxn = np.cross(m, n)
    a, th = np.pi/5, np.pi/3
    P = qmul(qexp_pure(m, a), qexp_pure(n, th))   # e^{ma} e^{nθ}
    Q = qmul(qexp_pure(n, th), qexp_pure(m, a))   # e^{nθ} e^{ma}
    comm = P - Q
    pred = np.concatenate([[0.0], 2*np.sin(a)*np.sin(th)*mxn])
    assert np.allclose(comm, pred), (comm, pred)
    fig = plt.figure(figsize=(6.8, 6.0)); ax = fig.add_subplot(111, projection="3d")
    def arrow(v, col, lab, lw=2.2):
        ax.quiver(0, 0, 0, *v, color=col, lw=lw, arrow_length_ratio=0.12)
        ax.text(*(np.asarray(v)*1.1), lab, color=col, fontsize=10.5)
    arrow(m, C_A, L("m（振幅の単位）", "m (unit of amplitude)"))
    arrow(n, C_P, L("n（位相の単位）", "n (unit of phase)"))
    arrow(mxn, C_T, L("m × n（第三の方向＝空間）", "m × n (third direction = space)"))
    # 合成回転の軸（順序で m×n 成分の符号が反転する）
    for q, col, lab in ((P, "#7a3fb0", L("e^{ma}e^{nθ} の軸", "axis of e^{ma}e^{nθ}")), (Q, "#b08a00", L("e^{nθ}e^{ma} の軸", "axis of e^{nθ}e^{ma}"))):
        ax_ = q[1:]/np.linalg.norm(q[1:])
        arrow(ax_*0.9, col, lab, lw=1.6)
    ax.set_xlim(-0.2, 1.2); ax.set_ylim(-0.2, 1.2); ax.set_zlim(-0.8, 1.2); ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=-55); ax.set_axis_off()
    ax.set_title(L("四元数：振幅軸 m と位相軸 n の回転は可換でなく、\n交換子 e^{ma}e^{nθ} − e^{nθ}e^{ma} ＝ 2 sin a sin θ (m × n) が第三の方向を生む",
                   "Quaternions: rotations about the amplitude axis m and the phase axis n do not commute;\nthe commutator e^{ma}e^{nθ} − e^{nθ}e^{ma} = 2 sin a sin θ (m × n) generates the third direction"), fontsize=10)
    ax.text2D(0.02, 0.02, L(f"a ＝ π/5, θ ＝ π/3　交換子 ＝ {comm[3]:.4f} (m×n)　予測 2 sin a sin θ ＝ {pred[3]:.4f}",
                             f"a = π/5, θ = π/3   commutator = {comm[3]:.4f} (m×n)   predicted 2 sin a sin θ = {pred[3]:.4f}"),
              transform=ax.transAxes, fontsize=8.5, color=C_MUTED)
    save(fig, "fig5_quaternion_commutator_third_axis")


# ================================================================ 図 6
def fig6():
    cols = [
        (L("ℂ にすでにあるもの", "Already in ℂ"), C_P, [
            L("振幅と位相の双対（指数平面の 90° 回転）", "amplitude–phase duality (90° rotation of the exponent plane)"),
            L("虚数振幅 → 不定計量・光円錐", "imaginary amplitude → indefinite norm, light cone"),
            L("複素角の二断面：cos²（ボルン則）と cosh φ ＝ γ", "two sections of the complex angle: cos² (Born) and cosh φ = γ"),
            L("自己双対点 |k| ＝ 1（α ＝ π/2）", "self-dual point |k| = 1 (α = π/2)"),
            L("Σz ＝ 0 と Σz² ＝ 0、U^n ＝ I、m/n の数え上げ", "Σz = 0 and Σz² = 0, U^n = I, counting m/n"),
            L("割り算の罠、指数部での線形性", "the division trap, linearity in the exponent"),
            L("振幅 ≠ 曲率半径（円柱 S¹ × ℝ で見える）", "amplitude ≠ curvature radius (visible on S¹ × ℝ)"),
        ]),
        (L("ℍ が足すもの", "Added by ℍ"), C_T, [
            L("振幅と位相が区別できる（m ⊥ n）", "amplitude and phase become distinguishable (m ⊥ n)"),
            L("非可換性 [m, n] ＝ 2m × n → 第三の方向＝空間", "non-commutativity [m, n] = 2m × n → third direction = space"),
            L("本当の二乗の符号 (1, 3)、向きを持つ光円錐", "signature (1, 3) of the true square; light cone with direction"),
            L("二つの定常波が位置を記憶する", "two standing waves remember a position"),
            L("凝縮体が作れる（複素数では方向がない）", "condensates become possible (ℂ has no direction)"),
            L("3 体の判別式 Δ₃ ＝ ¼|uv − vu|²", "3-body discriminant Δ₃ = ¼|uv − vu|²"),
        ]),
        (L("𝕆 が足すもの", "Added by 𝕆"), C_A, [
            L("非結合性：結合子は 3 個で立つ、閉塞まで課すと 4 個", "non-associativity: associator needs 3; with closure, 4"),
            L("R 型と Q 型の振幅が代数的に別物 → 質量と電荷の区別", "R-type and Q-type amplitudes become distinct → mass vs charge"),
            L("分解型 (4, 4) と |a| ＝ |b|：時空方程式", "split form (4, 4) and |a| = |b|: the spacetime equation"),
            L("逐次の結合の順序が結果を変える（履歴）", "order of sequential bonding matters (history)"),
            L("Dyson の集団が存在しない（前例なし）", "no Dyson ensemble (no precedent)"),
        ]),
    ]
    import textwrap
    fig, ax = plt.subplots(figsize=(13, 6.8)); ax.set_axis_off()
    x0 = 0.01; wcol = 0.32; width = 38 if EN else 24; lh = 0.048   # 1 行の高さ（軸座標）
    for i, (title, col, items) in enumerate(cols):
        x = x0 + i*(wcol + 0.01)
        ax.add_patch(FancyBboxPatch((x, 0.02), wcol, 0.96, boxstyle="round,pad=0.01", fc="white", ec=col, lw=2, transform=ax.transAxes))
        ax.text(x + wcol/2, 0.93, title, ha="center", va="center", fontsize=13, color=col, weight="bold", transform=ax.transAxes)
        y = 0.85
        for it in items:
            lines = textwrap.wrap(it, width=width) if EN else textwrap.wrap(it, width=width)
            txt = "• " + "\n  ".join(lines)
            ax.text(x + 0.015, y, txt, ha="left", va="top", fontsize=9.4, color=C_INK, transform=ax.transAxes, linespacing=1.25)
            y -= lh*len(lines) + 0.035
    fig.suptitle(L("帰属の台帳：段 L の性質と言えるのは、段 L − 1 の代数では述べられないことを示したときだけ",
                   "Attribution ledger: a property belongs to level L only if it cannot be stated with the algebra of level L − 1"), fontsize=11.5, y=0.995)
    save(fig, "fig6_attribution_ledger")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
