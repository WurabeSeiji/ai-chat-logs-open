#!/usr/bin/env python3
"""層 2：吸収体に温度 T を与える（STATUS §12 の三層のうち生成子が変わる層）。

固定一式 ../exchange_rel_20261005/exchange_cascade.py の準位 ETOT と自発放出の A 係数をそのまま使い、
電磁の各辺 i→f（ħω = E_i − E_f）に Planck の占有数 n̄(ω) = 1/(e^{ħω/kT} − 1) で
  下向き（自発＋誘導）  R_fi = A (1 + n̄)
  上向き（吸収）        R_if = A n̄ g_i/g_f、g = 2F + 1
を入れる（詳細釣り合い：平衡で p_i/p_f = (g_i/g_f) e^{−ħω/kT}）。重力波の辺には温度を掛けない（熱的な重力子の背景は置かない）。
二光子の辺も掛けない（ħω ≈ 5 eV、n̄ ≈ e^{−2×10⁴}）。
T = 0 で固定一式に戻る。T = 2.7255 K（宇宙背景放射、CODATA/Planck）で 1s の超微細が背景と平衡する。
出力：absorber_temperature.md、audit_temperature_T*.txt
"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
FIXED = HERE.parent / "exchange_rel_20261005"
sys.path.insert(0, str(FIXED))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None
import exchange_cascade as X  # noqa: E402

K_B_EV = 8.617333262e-5          # Boltzmann 定数 [eV/K]
T_CMB = 2.7255                   # 宇宙背景放射の温度 [K]（Fixsen 2009）
EV = X.EV
G_STAT = np.array([int(2 * X.jf_of(l, se, sp)[1] + 1) for (n, l, se, sp) in X.ST])   # 2F+1


def nbar(hw_eV, T):
    if T <= 0.0:
        return 0.0
    x = hw_eV / (K_B_EV * T)
    return 0.0 if x > 700.0 else 1.0 / np.expm1(x)


def generator_T(T):
    """生成行列 R(T) と、辺ごとの (src, dst, kind, A, n̄, 上向き率)。"""
    ed = X.edges(X.W_GW)
    R = np.zeros((X.DIM, X.DIM))
    rows = []
    for src, dst, kind, A in ed:
        hw = (X.ETOT[src] - X.ETOT[dst]) * EV
        n = nbar(hw, T) if kind in ("e1", "e2", "m1") else 0.0
        down = A * (1.0 + n)
        up = A * n * G_STAT[src] / G_STAT[dst]
        R[dst, src] += down
        R[src, src] -= down
        if up > 0.0:
            R[src, dst] += up
            R[dst, dst] -= up
        rows.append((src, dst, kind, A, n, up, hw))
    return R, rows


def run_T(T, t_max=1.0e17, n_times=300):
    """時間発展（Radau）。T > 0 では t_max = 1e17 s ≫ 最遅の緩和（1s 超微細、2×10¹² s）なので p(t_max) が定常分布 π。
    率が 6e8〜1e-15 s⁻¹ と 24 桁またぐため、R の固有分解や SVD で零空間を取るのは数値的に当てにならず、長時間積分で取る。
    定常の検査は残差 max|R π| で行う。"""
    R, rows = generator_T(T)
    p0 = np.zeros(X.DIM)
    p0[X.IDX[(5, 1, 1, -1)]] = 1.0
    times = np.concatenate([[0.0], np.logspace(-12, np.log10(t_max), n_times)])
    sol = solve_ivp(lambda t, p: np.dot(R, p), (0.0, t_max), p0, method="Radau", jac=lambda t, p: R,
                    t_eval=times, rtol=1e-11, atol=1e-30)
    hist = sol.y.T
    pi = hist[-1] if T > 0 else None
    resid = float(np.max(np.abs(np.dot(R, hist[-1]))))
    return dict(R=R, rows=rows, p0=p0, times=times, hist=hist, pi=pi, resid=resid)


def label(i):
    n, l, se, sp = X.ST[i]
    j, F = X.jf_of(l, se, sp)
    return "%d%s_{%d/2} F=%d" % (n, {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}[l], int(2 * j), int(F))


L = []
L.append("# 層 2：吸収体の温度。固定一式の準位と A 係数に誘導放出と吸収を加える")
L.append("")
L.append("初期 5p₃/₂ F=1、t_max = 1e17 s。下向き A(1+n̄)、上向き A n̄ g_i/g_f、n̄ = 1/(e^{ħω/kT} − 1)。重力波と二光子には温度を掛けない。")
L.append("")
i1, i0 = X.IDX[(1, 0, 1, 1)], X.IDX[(1, 0, 1, -1)]
hw21 = (X.ETOT[i1] - X.ETOT[i0]) * EV
for T in (0.0, T_CMB):
    res = run_T(T)
    R, rows, hist, pi = res["R"], res["rows"], res["hist"], res["pi"]
    pT = hist[-1]
    L.append("## T = %.4f K" % T)
    L.append("")
    if T > 0:
        n21 = nbar(hw21, T)
        L.append("| 量 | 値 |")
        L.append("|---|---|")
        L.append("| kT | %.4e eV |" % (K_B_EV * T))
        L.append("| 21 cm：ħω、ħω/kT、n̄ | %.4e eV、%.5f、%.2f |" % (hw21, hw21 / (K_B_EV * T), n21))
        L.append("| Lyman-α：ħω/kT、n̄ | %.3e、%.1e |" % (10.1988 / (K_B_EV * T), nbar(10.1988, T)))
        L.append("| 1s F=1→0 の率：自発 A、誘導 A n̄、吸収 0→1 A n̄ g₁/g₀ | %.3e、%.3e、%.3e s⁻¹ |" % (-R[i0, i0] * 0 + [r[3] for r in rows if r[0] == i1 and r[1] == i0][0], [r[3] * r[4] for r in rows if r[0] == i1 and r[1] == i0][0], [r[5] for r in rows if r[0] == i1 and r[1] == i0][0]))
        relax = [r[3] * (1 + r[4]) + r[5] for r in rows if r[0] == i1 and r[1] == i0][0]
        L.append("| 1s 超微細の緩和率 A(1+n̄) + A n̄ g₁/g₀、時定数 | %.3e s⁻¹、%.3e s = %.2e 年 |" % (relax, 1 / relax, 1 / relax / 3.15576e7))
        L.append("")
    L.append("### 終状態 p(t_max)（T > 0 では定常分布 π。残差 max|R p| = %.1e s⁻¹）" % res["resid"])
    L.append("")
    L.append("| 準位 | p(t_max) |")
    L.append("|---|---|")
    for i in range(X.DIM):
        if pT[i] > 1e-12:
            L.append("| %s | %.6f |" % (label(i), pT[i]))
    L.append("")
    if T > 0:
        x = hw21 / (K_B_EV * T)
        ratio_th = 3.0 * np.exp(-x)
        ratio = pT[i1] / pT[i0]
        T_spin = hw21 / (K_B_EV * np.log(3.0 * pT[i0] / pT[i1]))
        L.append("1s の F=1/F=0：p(t_max) の比 %.6f、理論 3e^{−ħω/kT} = %.6f、差 %.1e。読み出せる spin 温度 T_spin = (ħω/k)/ln(3p₀/p₁) = %.5f K（与えた T = %.4f K）。"
                 % (ratio, ratio_th, ratio - ratio_th, T_spin, T))
        # 定常での放出と吸収の釣り合い
        P_down = sum(r[3] * (1 + r[4]) * pi[r[0]] * r[6] for r in rows)
        P_up = sum(r[5] * pi[r[1]] * r[6] for r in rows)
        L.append("定常状態での放出パワー Σ A(1+n̄) π_i ħω = %.6e eV/s、吸収パワー Σ A n̄ (g_i/g_f) π_f ħω = %.6e eV/s、差 %.1e（詳細釣り合い）。"
                 % (P_down, P_up, P_down - P_up))
        drop = (res["p0"] @ X.ETOT - pi @ X.ETOT) * EV
        L.append("初期から定常までの正味の放出 p₀·E − π·E = %.9f eV（T = 0 では 13.054539268 eV。差 %.3e eV は 1s F=1 に残る分 p(F=1) × 5.87e-6）。" % (drop, 13.054539268 - drop))
        L.append("")
        L.append("5p からの多段（E1、二光子、E2）は T で変わらない（hν/kT ≥ 4×10⁴）。変わるのは 1s 超微細だけで、吸収状態 F=0 が定常分布 F=1 : F=0 = 3e^{−ħω/kT} : 1 に置き換わり、"
                 "緩和が自然放出だけの 1.1×10⁷ 年から %.1e 年に縮む。原子は背景の温度計になる。" % (1 / relax / 3.15576e7))
    else:
        L.append("T = 0 は固定一式そのもの（吸収状態 1s F=0、放出 13.054539268 eV）。")
    L.append("")
    lines = ["T_K=%.4f" % T, "final_p=%s" % [(X.ST[i], float(pT[i])) for i in range(X.DIM) if pT[i] > 1e-12],
             "steady_residual_max_abs_Rp=%.3e" % res["resid"], "sum_p=%r" % float(pT.sum())]
    (HERE / ("audit_temperature_T%.4f.txt" % T)).write_text("\n".join(lines) + "\n", encoding="utf-8")

out = "\n".join(L) + "\n"
(HERE / "absorber_temperature.md").write_text(out, encoding="utf-8")
print(out)
