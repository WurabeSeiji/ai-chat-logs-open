#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2編（動力学導出論文）用の図 G1〜G7 を生成する（PNG/SVG・系列データJSONも保存）

G1 正単体定理: 真空スペクトルの完全縮退と導出値 1/√(2M)（v4 実測との一致）
G2 力の法則の検証: Dyn-1 式 vs 頂点実装（散布・機械精度一致）
G3 運動学の三分解: Ω ブロック構造の概念図
G4 セクター＝復調の裁定（本命・3パネル）:
   (a) 重力側は巻きに盲目（3ケースの分離角が重なる）
   (b) 巻き一致選択則（等巻き類 vs 不等巻き類の 4.3° 差）
   (c) 復調側は Δw の符号を区別（η モードスペクトルの差）
G5 束縛振動: 静止対の大振幅有界振動（平均分離 60°→約39°へ低下、衝突のみ）
G6 調和閉鎖と逆二乗: |ω_n|Δθ_n=Ω と α_n∝Δθ^{-2}（log-log 勾配 −2）
G7 離散 Gauss 定理: S・∂S・境界流束の概念図

決定的（プローブと同一シード・同一構成の再計算）。系列は fig_data_f2_v1.json に保存。
使い方: python3 make_f2_figures_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def _save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  {stem}.png/.svg")


# ---------------------------------------------------------------- 共有: 二体正本
_uni = _load("uni_f2fig", UNI / "unified_interaction_v1.py")
_cr0 = _load("cr0_f2fig", EXP / "run_cr0_control_no_theta_v2.py")
base = _uni.two_body_base
step = _uni.collision_step_exact
SP = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
NC, NE = int(SP.chi_grid_n), int(SP.eta_grid_n)
SLOPE, ICEPT, _ = _cr0.calibrate_shift(SP, NC, NE)
PACKET = tuple(range(1, 18))
ETA = 2.0 * np.pi * np.arange(NE) / NE


def make_ab(mA, mB, packet_b=PACKET):
    case = base.explicit_packet_case(
        mode=f"f2fig_{mA}_{mB}", packet_a=PACKET, packet_b=packet_b,
        packet_a_shift=_cr0.shift_for_deg(-30.0, SLOPE, ICEPT),
        packet_b_shift=_cr0.shift_for_deg(+30.0, SLOPE, ICEPT))
    a = base.make_case_state(SP, case, "A", hair_enabled=True)
    b = base.make_case_state(SP, case, "B", hair_enabled=True)
    a = (a.reshape(NC, NE) * np.exp(1j * mA * ETA)[None, :]).reshape(-1)
    b = (b.reshape(NC, NE) * np.exp(1j * mB * ETA)[None, :]).reshape(-1)
    return a, b


def sep_of(a, b):
    ta, _ = _cr0.circle_position(a, NC, NE)
    tb, _ = _cr0.circle_position(b, NC, NE)
    return abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))


def run_seps(mA, mB, T):
    a, b = make_ab(mA, mB)
    out = []
    for _ in range(T):
        a, b, _ = step(a, b, SP)
        out.append(sep_of(a, b))
    return np.array(out)


def eta_spec(psi):
    f = np.fft.fft(psi.reshape(NC, NE), axis=1)
    P = np.sum(np.abs(f) ** 2, axis=0)
    m = ((np.arange(NE) + NE // 2) % NE) - NE // 2
    o = np.argsort(m)
    return m[o], (P / P.sum())[o]


DATA = {}


# ---------------------------------------------------------------- G1
def g1():
    N = 16
    M = N * (N - 1) // 2
    lam = 1.0 / (2 * M)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax1.bar(np.arange(1, N), [np.sqrt(lam)] * (N - 1), color="tab:blue", width=0.7)
    ax1.axhline(0.0645, color="tab:red", ls="--", lw=1.4,
                label="v4 実測 0.0645（三系模型・後期）")
    ax1.set_title("真空スペクトルの完全縮退（導出）", fontsize=10)
    ax1.set_xlabel("主軸番号")
    ax1.set_ylabel("符号付き √λ")
    ax1.set_ylim(0, 0.08)
    ax1.legend(fontsize=8)
    ax1.annotate(r"$\sqrt{\lambda}=1/\sqrt{2M}=0.06455$", (7.5, 0.068),
                 ha="center", fontsize=10)

    ax2.axis("off")
    ax2.set_title("導出値と実測（v4 主張14）の対照", fontsize=10)
    rows = [("λ（非自明・15重縮退）", "1/(2M)=1/240", "—"),
            ("√λ（主軸の値）", "0.064550", "0.0645"),
            ("tr(B)", "1/N=0.062500", "0.0625000"),
            ("r_rms", "1/N=0.062500", "0.062500"),
            ("上位3占有率", "3/(N−1)=0.200", "0.2000"),
            ("虚方向・rank", "0本・15", "0本・15")]
    txt = "量                       導出              実測\n" + "-" * 46 + "\n"
    for r in rows:
        txt += f"{r[0]:<22s} {r[1]:<15s} {r[2]}\n"
    ax2.text(0.02, 0.9, txt, family="monospace", fontsize=9, va="top",
             transform=ax2.transAxes)
    fig.suptitle("G1  真空正単体定理——無名性＋規格化だけから全判定値が出る", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, "fig_g1_simplex_v1")


# ---------------------------------------------------------------- G2
def g2():
    rng = np.random.default_rng(7)
    n = 6
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    ia = np.array([e[0] for e in edges]); ib = np.array([e[1] for e in edges])
    m = len(edges)
    x = rng.normal(size=m) + 1j * rng.normal(size=m)
    R = rng.random(m)
    a2, z2 = np.abs(x) ** 2, x ** 2
    A = np.zeros(n); B = np.zeros(n, complex)
    AR = np.zeros(n); BR = np.zeros(n, complex)
    np.add.at(A, ia, a2); np.add.at(A, ib, a2)
    np.add.at(B, ia, z2); np.add.at(B, ib, z2)
    np.add.at(AR, ia, R * a2); np.add.at(AR, ib, R * a2)
    np.add.at(BR, ia, R * z2); np.add.at(BR, ib, R * z2)
    cA = A[ia] + A[ib] - 2 * a2; cB = B[ia] + B[ib] - 2 * z2
    cAR = AR[ia] + AR[ib] - 2 * R * a2; cBR = BR[ia] + BR[ib] - 2 * R * z2
    rate = 0.5j * (R * (cA * x - cB * np.conj(x)) + (cAR * x - cBR * np.conj(x)))
    dd2_num = 2 * np.real(np.conj(x) * rate)
    xb2 = np.conj(x) ** 2
    dd2_th = R * np.imag((B[ia] + B[ib]) * xb2) + np.imag((BR[ia] + BR[ib]) * xb2)

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot(dd2_th, dd2_num, "o", ms=7, color="tab:blue")
    lim = 1.1 * float(np.max(np.abs(dd2_th)))
    ax.plot([-lim, lim], [-lim, lim], "-", color="gray", lw=1)
    ax.set_xlabel(r"理論: $g_1[\,R\,\mathrm{Im}((B_u{+}B_v)\bar x^2)+\mathrm{Im}((B^R_u{+}B^R_v)\bar x^2)\,]$")
    ax.set_ylabel(r"実装: $2\,\mathrm{Re}(\bar x\,\dot x)$（頂点レート）")
    ax.set_title("G2  力の法則（Dyn-1）：式 vs 実装——最大差 7×10⁻¹⁵", fontsize=11)
    _save(fig, "fig_g2_force_law_v1")


# ---------------------------------------------------------------- G3
def g3():
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.axis("off")
    ax.set_title("G3  運動学の三分解——$\\Omega$ ブロックと体ゲージ成分の一対一対応",
                 fontsize=11)
    def box(x, y, w, h, text, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03",
                                    fc=fc, ec="black", lw=1.1,
                                    transform=ax.transAxes))
        ax.annotate(text, (x + w / 2, y + h / 2), ha="center", va="center",
                    fontsize=9.5, xycoords="axes fraction")
    box(0.05, 0.62, 0.40, 0.26,
        "内部ブロック $\\Omega_{ij}$ ($i,j\\leq 3$)\n可視枠の剛体回転\n= $\\dot O$（ゲージ・ホロノミー）",
        "#dce9f7")
    box(0.55, 0.62, 0.40, 0.26,
        "混合ブロック $\\Omega_{i,m}$ ($i\\leq 3<m$)\n上位空間の尾部への回転\n= $\\dot\\theta$（体の運動の正体）",
        "#ffe6cc")
    box(0.05, 0.20, 0.40, 0.26,
        "$\\dot\\lambda_{1..3}$（トレースレス）\n楕円体の形の変化\n= $\\dot{\\widehat G}$（加速度・潮汐）",
        "#e8f5dc")
    box(0.55, 0.20, 0.40, 0.26,
        "トレース部 $\\sum\\dot\\lambda = 0$（系 Dyn-1a）\n= $\\dot\\rho = 0$\nスケールは頂点では動かない",
        "#f0f0f0")
    ax.annotate("$\\Omega_{ij}=\\frac{v_j^{T}\\dot B v_i}{\\lambda_i-\\lambda_j}$"
                "   （$\\dot B$ の源は端点の $B_v$——Dyn-1）",
                (0.5, 0.06), ha="center", fontsize=11, xycoords="axes fraction")
    _save(fig, "fig_g3_decomposition_v1")


# ---------------------------------------------------------------- G4
def g4():
    T = 200
    s00 = run_seps(0, 0, T)
    s33 = run_seps(3, 3, T)
    s3m3 = run_seps(3, -3, T)
    seq = run_seps(1, 0, T)          # 等巻き w=(2,2)
    DATA["G4_sep_uneq"] = s00.tolist()
    DATA["G4_sep_eq"] = seq.tolist()

    # 復調スペクトル（Δw=±2 終端）
    a1, b1 = make_ab(3, 0)    # w=(4,2) Δw=+2
    a2, b2 = make_ab(3, 4)    # w=(4,6) Δw=−2
    for _ in range(T):
        a1, b1, _ = step(a1, b1, SP)
        a2, b2, _ = step(a2, b2, SP)
    m1, P1 = eta_spec(a1)
    m2, P2 = eta_spec(a2)
    DATA["G4_demod_modes"] = m1.tolist()
    DATA["G4_demod_P_dwp2"] = P1.tolist()
    DATA["G4_demod_P_dwm2"] = P2.tolist()

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    ax = axes[0]
    ax.plot(s00, lw=2.2, color="tab:blue", label="(0,0)")
    ax.plot(s33, lw=1.2, ls="--", color="tab:orange", label="(+3,+3)")
    ax.plot(s3m3, lw=0.9, ls=":", color="tab:green", label="(+3,−3)")
    ax.set_title("(a) 重力側読出しは巻きに盲目\n3ケースが機械精度で重なる", fontsize=9.5)
    ax.set_xlabel("步"); ax.set_ylabel("分離角 [deg]"); ax.legend(fontsize=8)

    ax = axes[1]
    diff = np.abs(seq - s00)
    ax.plot(diff, lw=1.2, color="tab:red")
    ax.set_title("(b) 巻き一致選択則\n|等巻き類 − 不等巻き類|（最大 4.3°・恒常的に有限）",
                 fontsize=9.5)
    ax.set_xlabel("步"); ax.set_ylabel("分離角の差 [deg]")

    ax = axes[2]
    w = 0.38
    ax.bar(m1 - w / 2, P1, width=w, color="tab:purple", label=r"$\Delta w=+2$")
    ax.bar(m2 + w / 2, P2, width=w, color="tab:brown", label=r"$\Delta w=-2$")
    ax.set_title("(c) 復調側は $\\Delta w$ の符号を区別\nA チャネル η モードスペクトル（終端）",
                 fontsize=9.5)
    ax.set_xlabel("η モード（符号付き巻き）"); ax.set_ylabel("パワー比"); ax.legend(fontsize=8)

    fig.suptitle("G4  セクター＝復調——同一の相互作用の、盲目な読出し（重力）と符号を見る読出し（電荷）",
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    _save(fig, "fig_g4_sectors_v1")


# ---------------------------------------------------------------- G5
def g5():
    T = 400
    seps = run_seps(0, 0, T)
    DATA["G5_collapse"] = seps.tolist()
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.plot(seps, lw=0.9, color="tab:blue", alpha=0.85)
    ax.axhline(60.0, color="gray", ls="--", lw=1.2, label="初期分離 60°（衝突なしなら不変）")
    mu = float(np.mean(seps))
    ax.axhline(mu, color="tab:red", ls="-", lw=1.6, label=f"平均分離 {mu:.1f}°（引力側へ低下）")
    ax.set_xlabel("步"); ax.set_ylabel("分離角 [deg]")
    ax.set_title("G5  束縛振動——静止対は衝突のみで大振幅の有界振動（平均低下・符号盲目）",
                 fontsize=10.5)
    ax.legend(fontsize=9, loc="lower right")
    _save(fig, "fig_g5_collapse_v1")


# ---------------------------------------------------------------- G6
def g6():
    omega1 = np.pi / 72.0
    ns = np.arange(1, 33)
    dth = 2 * np.pi / ns
    alpha = (ns * omega1) ** 2  # R=1
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    ax1.plot(ns, np.abs(ns * omega1) * dth, "o-", ms=4, color="tab:blue")
    ax1.axhline(2 * np.pi * omega1, color="tab:red", ls="--", lw=1.2,
                label=r"$\Omega_N=2\pi\omega_1=2\pi^2/72$")
    ax1.set_xlabel("倍音 n"); ax1.set_ylabel(r"$|\omega_n|\,\Delta\theta_n$")
    ax1.set_title("調和閉鎖の双対（恒等式・ずれ $10^{-17}$）", fontsize=10)
    ax1.legend(fontsize=8)
    ax2.loglog(dth, alpha, "o-", ms=4, color="tab:green")
    ax2.set_xlabel(r"$\Delta\theta_n$"); ax2.set_ylabel(r"$\alpha_n=R\omega_n^2$")
    ax2.set_title("逆二乗（log-log 勾配 厳密 −2）", fontsize=10)
    fig.suptitle("G6  N体調和閉鎖——[AB] 双対の N 体版（集団時計 ω₁=π/72）", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    _save(fig, "fig_g6_harmonic_closure_v1")


# ---------------------------------------------------------------- G7
def g7():
    rng = np.random.default_rng(5)
    n = 10
    pos = rng.normal(size=(n, 2)) * 1.2
    pos -= pos.mean(axis=0)
    S = [i for i in range(n) if np.linalg.norm(pos[i]) < 1.1]
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    ax.set_aspect("equal")
    for i in range(n):
        for j in range(i + 1, n):
            inS = (i in S) + (j in S)
            if inS == 2:
                c, lw, ls = "tab:blue", 1.2, "-"
            elif inS == 1:
                c, lw, ls = "tab:red", 1.8, "-"
            else:
                c, lw, ls = "lightgray", 0.7, ":"
            ax.plot(*np.array([pos[i], pos[j]]).T, color=c, lw=lw, ls=ls, zorder=1)
    for i in range(n):
        ax.plot(*pos[i], "o", ms=10, zorder=3,
                color=("tab:blue" if i in S else "gray"))
    circ = Circle((0, 0), 1.1, fill=False, color="black", ls="--", lw=1.4)
    ax.add_patch(circ)
    ax.annotate("領域 S", (-1.05, 1.12), fontsize=11)
    ax.annotate("境界 ∂S を跨ぐ辺 = 境界流束（赤）", (0.02, 0.02),
                xycoords="axes fraction", fontsize=9.5,
                bbox=dict(boxstyle="round", fc="white", ec="gray"))
    ax.set_title("G7  離散 Gauss 定理:  $\\sum_{v\\in S}B_v = 2\\sum_{e\\subseteq S}z_e^2 + \\sum_{e\\in\\partial S}z_e^2$",
                 fontsize=11)
    ax.axis("off")
    _save(fig, "fig_g7_discrete_gauss_v1")


# ----------------------------------------------------------------
if __name__ == "__main__":
    t0 = time.time()
    print("第2編 図 G1〜G7 を生成:")
    g1(); g2(); g3(); g4(); g5(); g6(); g7()
    (HERE / "fig_data_f2_v1.json").write_text(
        json.dumps(DATA, indent=1, ensure_ascii=False))
    print(f"系列データ: fig_data_f2_v1.json  ({time.time()-t0:.1f}s)")
