#!/usr/bin/env python3
"""連続状態 → 束縛準位の放射会合 H + H → H₂ + 光子（一重項 X¹Σg⁺ 上、E2・M1・GW）。
箱（R ≤ R_max）で規格化した離散化連続状態 ψ_n（E_n > 0、部分波 J）から束縛準位 (v',J') への率 A を、束縛–束縛と同じ式（system_h2.rates）で計算する。
  通過 1 回あたりの捕獲確率   P_J(E_n) = T_box Σ_f A(n→f)、T_box = 2πħ ρ_J(E_n)（箱の準位密度 ρ = dn/dE から。2πħρ は箱の往復時間）
  断面積                      σ(E) = (π/k²) Σ_J w_J (2J+1) P_J(E)、k² = 2μE
  率係数                      k(E) = σ(E) v、v = ħk/μ。熱平均 k(T) = ∫ k(E) f_MB(E) dE
w_J = 電子スピン一重項 1/4 × 核スピン（J 偶：パラ 1/4、J 奇：オルト 3/4）。三重項 b³Σu⁺（電子スピンの 3/4）は斥力で束縛準位が無く、この部品に無い。
P_J は R_max に依らない（ψ_n の井戸内振幅² ∝ 1/R_max、ρ ∝ R_max）。R_max = 50 → 100 で確認する。
出力：radiative_association_h2.md
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import system_h2 as S  # noqa: E402

CM = S.HARTREE_CM
K_B_CM = 0.695034800            # cm⁻¹/K
A0_CM = 5.29177210903e-9        # bohr in cm
E_MAX = 0.02                    # hartree（4390 cm⁻¹、0.54 eV）
J_MAX = 24                      # E = 3000 cm⁻¹ で k b ≈ 5 × 4 = 20 まで井戸に届くので、それより上まで
E_SEL_CM = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]
W_SPIN = {0: 0.25 * 0.25, 1: 0.25 * 0.75}   # J 偶／奇


def levels_J(J, emax, dr=0.025, rmax=50.0):
    E, F, R = S.dvr_levels(J, dr=dr, rmax=rmax, emax=emax)
    return [dict(v=(v if E[v] < 0 else None), J=J, E_h=float(E[v]), f=F[:, v]) for v in range(len(E))], R


def capture_table(J, bound, R, cont):
    """部分波 J の各連続状態について (E_n, ρ_n, ΣA_e2, ΣA_m1, ΣA_gw, 終準位ごとの A_e2)。"""
    Ec = np.array([c["E_h"] for c in cont])
    rho = np.empty(len(Ec))
    rho[1:-1] = 2.0 / (Ec[2:] - Ec[:-2]); rho[0] = 1.0 / (Ec[1] - Ec[0]); rho[-1] = 1.0 / (Ec[-1] - Ec[-2])
    finals = [b for b in bound if abs(b["J"] - J) in (0, 2)]
    rows = []
    for n, c in enumerate(cont):
        a_e2 = a_m1 = a_gw = 0.0
        per = {}
        for b in finals:
            e2, m1, gw = S.rates(c, b, R)
            a_e2 += e2; a_m1 += m1; a_gw += gw
            if e2 > 0:
                per[(b["v"], b["J"])] = e2
        rows.append((c["E_h"], rho[n], a_e2, a_m1, a_gw, per))
    return rows


def main():
    L = ["# H + H → H₂ + 光子：連続状態から束縛準位への放射会合（一重項 X¹Σg⁺、E2・M1・GW）", "",
         "箱 R ≤ 50 bohr（dr 0.025）で規格化した離散化連続状態から、束縛–束縛と同じ率の式で捕獲率を計算。P_J(E) = 2πħρ_J ΣA は通過 1 回あたりの捕獲確率、σ(E) = (π/k²)Σ_J w_J(2J+1)P_J、k(E) = σv。",
         "重み w_J = 電子一重項 1/4 × 核スピン（J 偶 1/4、J 奇 3/4）。三重項（電子スピンの 3/4）は斥力で捕獲が無い。E ≤ %.0f cm⁻¹、J ≤ %d。" % (E_MAX * CM, J_MAX), ""]
    # 束縛準位（J ≤ J_MAX+2）と連続状態（J ≤ J_MAX）
    bound, cont = {}, {}
    for J in range(J_MAX + 3):
        lv, R = levels_J(J, E_MAX if J <= J_MAX else 0.0)
        bound[J] = [d for d in lv if d["E_h"] < 0]
        cont[J] = [d for d in lv if d["E_h"] > 0]
    allb = [b for J in bound for b in bound[J]]
    tabs = {J: capture_table(J, allb, R, cont[J]) for J in range(J_MAX + 1)}

    # 1. 部分波ごとの P_J(E_n)（箱の準位で）
    L += ["## 1. 通過 1 回あたりの捕獲確率 P_J(E_n)（箱の離散化連続準位そのまま、J = 0・1 の最初の 12 点と J = 4 の最初の 6 点）", "",
          "| J | E_n [cm⁻¹] | E_n/k_B [K] | ρ_J [1/cm⁻¹] | T_box = 2πħρ [s] | ΣA_e2 [s⁻¹] | ΣA_m1 | ΣA_gw | P_J（E2＋M1） | 主な終準位 (v',J')：分岐 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for J, npts in ((0, 12), (1, 12), (4, 6)):
        for E_h, rho, a2, a1, ag, per in tabs[J][:npts]:
            Tbox = 2 * np.pi * rho * S.AU_TIME_S
            top = sorted(per.items(), key=lambda kv: -kv[1])[:3]
            L.append("| %d | %.3f | %.2f | %.3e | %.3e | %.3e | %.1e | %.1e | %.3e | %s |" % (
                J, E_h * CM, E_h * CM / K_B_CM, rho / CM, Tbox, a2, a1, ag, Tbox * (a2 + a1),
                "、".join("(%d,%d) %.2f" % (k[0], k[1], v / a2) for k, v in top)))
    r4 = tabs[4][0]
    L += ["", "J = 4 の最下の箱準位（E = %.3f cm⁻¹、P₄ = %.2e）は (14,4) の準束縛準位：本計算（非断熱なし）では解離限界のすぐ上に出る（全補正では 0.027 cm⁻¹ 下で束縛、validation_h2.md §4）。"
          "この共鳴が E ≈ 1 cm⁻¹ の σ を押し上げる（下の表の E = 1 cm⁻¹ の行）。非断熱補正を入れれば束縛準位に戻り、連続状態からは消える。" % (r4[0] * CM, 2 * np.pi * r4[1] * S.AU_TIME_S * (r4[2] + r4[3])), ""]

    # 2. 共通の E 格子へ補間（log–log）、σ(E)、k(E)
    def interp_P(J, E_h):
        Ec = np.array([r[0] for r in tabs[J]]); P = np.array([2 * np.pi * r[1] * S.AU_TIME_S * (r[2] + r[3]) for r in tabs[J]])
        lo, hi = Ec[0], Ec[-1]
        if E_h < lo:   # Wigner の閾値則 P ∝ E^{J+1/2} で下へ外挿
            return P[0] * (E_h / lo) ** (J + 0.5), True
        if E_h > hi:   # 最上の箱の準位より上（E_MAX との隙間 ≤ 1 準位間隔）は最上の値で代用
            return float(P[-1]), True
        return float(np.exp(np.interp(np.log(E_h), np.log(Ec), np.log(P)))), False

    L += ["## 2. 断面積 σ(E) と率係数 k(E) = σv（部分波の和、log–log 補間）", "",
          "| E [cm⁻¹] | E/k_B [K] | P₀ | P₁ | P₂ | P₃ | P₄ | σ [cm²] | k(E) [cm³ s⁻¹] | J=4 の寄与/全体 | J=%d の寄与/全体 |" % J_MAX, "|---|---|---|---|---|---|---|---|---|---|---|"]
    kE_grid = []
    for Ecm in E_SEL_CM:
        E_h = Ecm / CM
        k2 = 2 * S.MU_N * E_h
        tot = 0.0; last = 0.0; c4 = 0.0; Ps = {}
        for J in range(J_MAX + 1):
            P, ext = interp_P(J, E_h)
            Ps[J] = P
            c = W_SPIN[J % 2] * (2 * J + 1) * P
            tot += c
            if J == J_MAX:
                last = c
            if J == 4:
                c4 = c
        sigma_au = np.pi / k2 * tot
        v_au = np.sqrt(2 * E_h / S.MU_N)
        sigma_cm2 = sigma_au * A0_CM**2
        kE = sigma_au * v_au * A0_CM**3 / S.AU_TIME_S
        kE_grid.append((E_h, kE))
        L.append("| %g | %.1f | %.2e | %.2e | %.2e | %.2e | %.2e | %.3e | %.3e | %.2f | %.1e |" % (Ecm, Ecm / K_B_CM, Ps[0], Ps[1], Ps[2], Ps[3], Ps[4], sigma_cm2, kE, c4 / tot if tot > 0 else 0, last / tot if tot > 0 else 0))
    L.append("")

    # 3. 熱平均 k(T)
    Eg = np.logspace(np.log10(tabs[0][0][0]), np.log10(E_MAX), 400)
    kEf = []
    for E_h in Eg:
        k2 = 2 * S.MU_N * E_h
        tot = sum(W_SPIN[J % 2] * (2 * J + 1) * interp_P(J, E_h)[0] for J in range(J_MAX + 1))
        kEf.append(np.pi / k2 * tot * np.sqrt(2 * E_h / S.MU_N) * A0_CM**3 / S.AU_TIME_S)
    kEf = np.array(kEf)
    L += ["## 3. Maxwell–Boltzmann 熱平均 k(T) = ∫ k(E) (2/√π)(E/kT)^{1/2} e^{−E/kT} dE/kT（E ≤ %.0f cm⁻¹ で打ち切り）" % (E_MAX * CM), "",
          "| T [K] | k(T) [cm³ s⁻¹] | 打ち切りの重み（E > E_max の MB 確率） |", "|---|---|---|"]
    for T in (10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0):
        kT = K_B_CM * T / CM
        x = Eg / kT
        w = 2 / np.sqrt(np.pi) * np.sqrt(x) * np.exp(-x) / kT
        kT_avg = float(np.trapezoid(kEf * w, Eg))
        from scipy.special import gammaincc
        L.append("| %g | %.3e | %.1e |" % (T, kT_avg, gammaincc(1.5, E_MAX / kT)))
    L += ["", "T ≲ 30 K の k(T) には上の (14,4) 準束縛共鳴（E ≈ 0.7 cm⁻¹）が入っている。共鳴を除いた非共鳴の値は k(E) の 3〜10 cm⁻¹ の行（4〜6e-30 cm³ s⁻¹）が目安。T ≳ 100 K では共鳴の寄与は 1 割以下。", ""]

    # 4. 箱に依らないこと（J = 0、R_max 100）
    lv100, R100 = levels_J(0, E_MAX, rmax=100.0)
    b100 = [d for d in lv100 if d["E_h"] < 0]
    lv2, _ = levels_J(2, 0.0, rmax=100.0)
    c100 = [d for d in lv100 if d["E_h"] > 0]
    tab100 = capture_table(0, b100 + lv2, R100, c100)
    Ec50 = np.array([r[0] for r in tabs[0]]); P50 = np.array([2 * np.pi * r[1] * S.AU_TIME_S * r[2] for r in tabs[0]])
    Ec100 = np.array([r[0] for r in tab100]); P100 = np.array([2 * np.pi * r[1] * S.AU_TIME_S * r[2] for r in tab100])
    L += ["## 4. 検査", ""]
    for Ecm in (10.0, 100.0, 1000.0):
        E_h = Ecm / CM
        p50 = np.exp(np.interp(np.log(E_h), np.log(Ec50), np.log(P50))); p100 = np.exp(np.interp(np.log(E_h), np.log(Ec100), np.log(P100)))
        L.append("- J = 0、E = %g cm⁻¹：P₀ は R_max 50 で %.4e、100 で %.4e（相対差 %.1e）。箱の準位数は %d → %d。" % (Ecm, p50, p100, p100 / p50 - 1, len(tabs[0]), len(tab100)))
    # 往復時間の古典値との比較（J=0、E = 100 cm⁻¹）
    E_h = 100.0 / CM
    Rg = R[(S.v_el(R) + S.e_ad(R)) < E_h]
    vcl = np.sqrt(2 * (E_h - (S.v_el(Rg) + S.e_ad(Rg))) / S.MU_N)
    T_cl = 2 * float(np.trapezoid(1.0 / vcl, Rg)) * S.AU_TIME_S
    n0 = int(np.argmin(np.abs(Ec50 - E_h)))
    L.append("- 箱の往復時間 2πħρ（J = 0、E ≈ 100 cm⁻¹）：%.3e s、古典の往復 2∫dR/v（内側転回点〜R_max）：%.3e s。" % (2 * np.pi * tabs[0][n0][1] * S.AU_TIME_S, T_cl))
    L.append("- 事前の粗い見積もり（ħω = D_e の E2 率 α⁵ω⁵/15·Θ(R_e)² ≈ 2 s⁻¹ × 井戸の通過時間 5e-15 s ≈ 1e-14）は 5 桁過大だった。実際の捕獲は ħω ≈ 0.4〜1.1 eV で v' = 8〜11 へ入り（v' = 0 への ħω = 4.7 eV の行列要素は連続状態の振動との不整合で小さい）、ω⁵ が 1e-4 倍になる。")
    L.append("- 入れていないもの：三重項チャネル（斥力、捕獲なし）、電子励起状態への経路、非断熱・相対論・QED（準位と同じ）。σ は E2・M1 のみで GW は 1e-34 倍。")
    out = "\n".join(L) + "\n"
    (HERE / "radiative_association_h2.md").write_text(out, encoding="utf-8")
    print(out)
    # 図化用 JSON（上と同じ数値）
    import json
    dump = dict(E_MAX=E_MAX, J_MAX=J_MAX, W_SPIN={str(k): v for k, v in W_SPIN.items()},
                tabs={str(J): [dict(E_h=r[0], rho=r[1], a_e2=r[2], a_m1=r[3], a_gw=r[4],
                                    top=[[list(k), v] for k, v in sorted(r[5].items(), key=lambda kv: -kv[1])[:3]]) for r in tabs[J]] for J in tabs},
                sigma_rows=[], kT_rows=[], check_rmax=dict(E50=Ec50.tolist(), P50=P50.tolist(), E100=Ec100.tolist(), P100=P100.tolist()))
    for Ecm in E_SEL_CM:
        E_h = Ecm / CM
        k2 = 2 * S.MU_N * E_h
        Ps = {J: interp_P(J, E_h)[0] for J in range(J_MAX + 1)}
        tot = sum(W_SPIN[J % 2] * (2 * J + 1) * Ps[J] for J in Ps)
        dump["sigma_rows"].append(dict(E_cm=Ecm, P={str(J): Ps[J] for J in Ps}, sigma_cm2=np.pi / k2 * tot * A0_CM**2,
                                       k_cm3s=np.pi / k2 * tot * np.sqrt(2 * E_h / S.MU_N) * A0_CM**3 / S.AU_TIME_S))
    dump["kE_fine"] = dict(E_h=Eg.tolist(), k_cm3s=kEf.tolist())
    for T in (10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0):
        kT = K_B_CM * T / CM
        x = Eg / kT
        w = 2 / np.sqrt(np.pi) * np.sqrt(x) * np.exp(-x) / kT
        dump["kT_rows"].append(dict(T=T, k_cm3s=float(np.trapezoid(kEf * w, Eg))))
    (HERE / "radiative_association_h2.json").write_text(json.dumps(dump, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
