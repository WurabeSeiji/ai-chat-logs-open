#!/usr/bin/env python3
"""プログラム構造文書（プログラム構造_ja_20261006.md）の付表を本体から生成する。手で転記しない。
出力：structure_tables.md
  表 A  50 準位：添字、(n,ℓ,j,F)、ORB_C・ORB_G・HFS・ETOT [eV]、全崩壊率 Γ_i、寿命、滞在時間 τ_i、終状態 p_∞
  表 B  辺の種別ごとの本数と、A の大きい順 25 本
  表 C  種別ごとの放出エネルギー（流量 A·τ·ΔE の和）
  表 D  導出定数
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib import font_manager

font_manager.fontManager.addfont = lambda path: None
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exchange_cascade as X  # noqa: E402

names = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}
res = X.run(X.W_GW)
R, tau, p_inf, ed = res["R"], res["tau"], res["p_inf"], res["ed"]
L = []
L.append("# 付表（structure_tables.py が exchange_cascade.py から生成）")
L.append("")
L.append("## 表 D　導出定数")
L.append("")
L.append("| 記号 | 値 | 式 |")
L.append("|---|---|---|")
mu = X.MU_P / (1.0 + X.MU_P)
L.append("| α | %.12e | 1/137.035999177（CODATA 2022） |" % X.ALPHA)
L.append("| ħ（q₀ = c = 1、e = 3） | %.9f | 9/α |" % X.HBAR)
L.append("| μ_red/m_e | %.12f | m_p/(m_e + m_p) |" % mu)
L.append("| α_G | %.6e | G m_e²/(ħc)（SI） |" % X.ALPHA_G)
L.append("| MU（m_e、G = q₀ = c = 1） | %.6e | 3√(α_G/α) |" % X.MU)
L.append("| g = G m_e m_p/e² | %.6e | MU_P·MU²/9 |" % (X.MU_P * X.MU**2 / 9.0))
L.append("| W_GW | %.6e | 4 g m_e m_p/(m_p − m_e)² |" % X.W_GW)
L.append("| a（μ の Bohr 半径、SI） | %.6e m | ħ/(μ c α) |" % (X.HBAR_SI / (X.ME_SI * mu * X.C_SI * X.ALPHA)))
L.append("| K（超微細、m_e c² 単位） | %.6e | g_p α⁴ μ³/(m_p/m_e) |" % (X.G_P * X.ALPHA**4 * mu**3 / X.MU_P))
L.append("| m_e c² | %.5f eV | CODATA 2022 |" % X.ME_EV)
L.append("")
L.append("## 表 A　50 準位（状態ベクトル p の成分）")
L.append("")
L.append("| i | (n,ℓ,s_e,s_p) | 準位 | ORB_C [eV] | ORB_G [eV] | HFS [eV] | ETOT [eV] | Γ_i [s⁻¹] | 寿命 [s] | τ_i [s] | p_∞ |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
for i, (n, l, se, sp) in enumerate(X.ST):
    j, F = X.jf_of(l, se, sp)
    G_i = -R[i, i]
    life = "%.4e" % (1.0 / G_i) if G_i > 0 else "∞"
    L.append("| %d | (%d,%d,%+d,%+d) | %d%s_{%d/2} F=%d | %.9f | %.4e | %+.4e | %.9f | %.4e | %s | %.4e | %.3f |"
             % (i, n, l, se, sp, n, names[l], int(2 * j), int(F), X.ORB_C[i] * X.EV, X.ORB_G[i] * X.EV, X.HFS[i] * X.EV,
                X.ETOT[i] * X.EV, G_i, life, tau[i], p_inf[i]))
L.append("")
L.append("## 表 B　辺（遷移）")
L.append("")
kinds = {}
for src, dst, k, A in ed:
    kinds[k] = kinds.get(k, 0) + 1
L.append("種別ごとの本数：" + "、".join("%s %d" % kv for kv in sorted(kinds.items())) + "（合計 %d）" % len(ed))
L.append("")
L.append("A の大きい順 25 本：")
L.append("")
L.append("| src → dst | 種別 | A [s⁻¹] | ΔE [eV] | 分岐比 A/Γ_src |")
L.append("|---|---|---|---|---|")
def lab(i):
    n, l, se, sp = X.ST[i]
    j, F = X.jf_of(l, se, sp)
    return "%d%s_{%d/2}F%d" % (n, names[l], int(2 * j), int(F))
for src, dst, k, A in sorted(ed, key=lambda e: -e[3])[:25]:
    L.append("| %s → %s | %s | %.4e | %.6f | %.4f |" % (lab(src), lab(dst), k, A, (X.ETOT[src] - X.ETOT[dst]) * X.EV, A / (-R[src, src])))
L.append("")
L.append("種別ごとの最小・最大 A：")
L.append("")
L.append("| 種別 | 最小 A [s⁻¹] | 最大 A [s⁻¹] |")
L.append("|---|---|---|")
for k in sorted(kinds):
    As = [A for s, d, kk, A in ed if kk == k]
    L.append("| %s | %.3e | %.3e |" % (k, min(As), max(As)))
L.append("")
L.append("## 表 C　種別ごとの放出エネルギー Σ A·τ_src·ΔE")
L.append("")
L.append("| 種別 | 放出 [eV] |")
L.append("|---|---|")
for k, v in sorted(res["flux_E"].items()):
    L.append("| %s | %.9e |" % (k, v * X.EV))
tot = sum(res["flux_E"].values()) * X.EV
drop = (res["p0"] @ X.ETOT - p_inf @ X.ETOT) * X.EV
L.append("| 合計 | %.9f |" % tot)
L.append("")
L.append("準位の落差 p₀·E − p_∞·E = %.9f eV、差 %.1e eV。" % (drop, tot - drop))
L.append("")
# ---- 表 E・F：時間方向の幅と空間方向の幅 ----
# 時間：準位 i の幅は Γ_i = −R_ii。遷移 i→f の線幅（角振動数の半値全幅）は Γ_i + Γ_f、可干渉時間 1/(Γ_i+Γ_f)、
#       可干渉長 c/(Γ_i+Γ_f)。待ち時間の指数分布と Lorentz 線形は Fourier 対。
# 空間：準位の広がり ⟨r⟩ = (a/2)[3n² − ℓ(ℓ+1)]、⟨r²⟩ = (a²n²/2)[5n² + 1 − 3ℓ(ℓ+1)]、Δr = √(⟨r²⟩ − ⟨r⟩²)、
#       運動量の広がり Δp = ħ/(na)（⟨p²⟩ = ħ²/(n²a)²）。遷移の双極子長 a|R̃⁽¹⁾|、波長 λ = 2πc/ω、多重極の展開径数 x = ωa/c。
#       反跳シフト (ħω)²/(2Mc²) は重心運動を凍結しているので模型には無い（参考値）。
a_si = X.HBAR_SI / (X.ME_SI * mu * X.C_SI * X.ALPHA)
M_SI = X.ME_SI * (1.0 + X.MU_P)
L.append("## 表 E　準位の空間的な幅")
L.append("")
L.append("| 準位 | Γ_i [s⁻¹] | ⟨r⟩ [pm] | Δr [pm] | Δp [ħ/a] | Δr·Δp/ħ |")
L.append("|---|---|---|---|---|---|")
seen = set()
for i, (n, l, se, sp) in enumerate(X.ST):
    if (n, l) in seen:
        continue
    seen.add((n, l))
    r_mean = 0.5 * (3 * n**2 - l * (l + 1))
    r2_mean = 0.5 * n**2 * (5 * n**2 + 1 - 3 * l * (l + 1))
    dr = np.sqrt(r2_mean - r_mean**2)
    dp = 1.0 / n
    L.append("| %d%s | %.3e | %.2f | %.2f | %.3f | %.3f |" % (n, names[l], -R[i, i], r_mean * a_si * 1e12, dr * a_si * 1e12, dp, dr * dp))
L.append("")
L.append("## 表 F　遷移の時間的・空間的な幅（A の大きい順 25 本と二光子）")
L.append("")
L.append("| src → dst | 種別 | ħω [eV] | λ [nm] | 線幅 (Γ_i+Γ_f)/2π [MHz] | 可干渉時間 [s] | 可干渉長 [m] | 双極子長 a\\|R̃⁽¹⁾\\| [pm] | x = ωa/c | 反跳シフト [MHz]（参考） |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
sel = sorted(ed, key=lambda e: -e[3])[:25] + [e for e in ed if e[2] == "2ph"]
for src, dst, k, A in sel:
    dE = (X.ETOT[src] - X.ETOT[dst]) * X.EV
    w = dE * X.E_SI / X.HBAR_SI
    lam = 2 * np.pi * X.C_SI / w
    G = -R[src, src] - R[dst, dst]
    n, l, _, _ = X.ST[src]; n2, l2, _, _ = X.ST[dst]
    dip = a_si * abs(X.radial(n, l, n2, l2, 1)) * 1e12 if abs(l - l2) == 1 else float("nan")
    rec = (dE * X.E_SI) ** 2 / (2 * M_SI * X.C_SI**2) / X.E_SI / (1e6 * 4.135667696e-15)
    L.append("| %s → %s | %s | %.6f | %.3f | %.4g | %.3e | %.3e | %s | %.3e | %.2f |"
             % (lab(src), lab(dst), k, dE, lam * 1e9, G / (2 * np.pi) / 1e6, 1.0 / G, X.C_SI / G,
                ("%.2f" % dip) if dip == dip else "—", w * a_si / X.C_SI, rec))
L.append("")
L.append("二光子の行の線幅は和の振動数の幅（Γ_2s + Γ_1s）。個々の光子の帯域は連続体の全幅 ≈ ΔE/ħ で、可干渉時間 ħ/ΔE ≈ %.1e s、長さ %.0f nm。"
         % (X.HBAR_SI / (10.2 * X.E_SI), X.C_SI * X.HBAR_SI / (10.2 * X.E_SI) * 1e9))
L.append("")
(HERE / "structure_tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("structure_tables.md: %d lines" % len(L))
