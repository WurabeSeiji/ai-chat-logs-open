#!/usr/bin/env python3
"""H₂ の走行：共通エンジン＋水素分子の部品。水素原子の一式と同じ読み出し（放出の合計、滞在時間、線の幅、層 2・1・3）。
使い方：python3 run_h2.py [v J]   （初期準位、既定 1 3 = 1-0 S(1) 線の上準位）
出力：audit_h2_v{v}J{J}.txt
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import engine as EN  # noqa: E402
import system_h2 as S  # noqa: E402

CM = S.HARTREE_CM
v0q, J0q = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) >= 3 else (1, 3)
states, channels, init, aux = S.build(init=(v0q, J0q))
idx = aux["idx"]
n = len(states)
p0 = np.zeros(n); p0[init] = 1.0
E_cm = np.array([s["E_eV"] / S.HARTREE_EV * CM for s in states])
L = ["# H₂（X¹Σg⁺ 振動回転）：共通エンジンによる走行、初期準位 (v,J) = (%d,%d) %s" % (v0q, J0q, states[init]["name"]), "",
     "準位：BO（Pachucki 2010）＋断熱（Pachucki–Komasa 2014）、sinc-DVR dr = %.3f、R ≤ %.0f bohr、μ_n = m_p/2。束縛準位 %d、辺 %d（e2 %d、m1 %d、gw %d）。" % (
         aux["dr"], aux["rmax"], n, len(channels), sum(c[2] == "e2" for c in channels), sum(c[2] == "m1" for c in channels), sum(c[2] == "gw" for c in channels)),
     "率：E2（Θ = Q_WSD/2）、M1（g(R)）、GW（4(α_Gp/α)|⟨q_m⟩|²/|⟨Θ⟩|² × E2）。検証は validation_h2.md。", ""]

# ---------- T = 0：厳密解 ----------
R0, rows = EN.generator(states, channels, 0.0)
tau, p_inf = EN.exact(R0, p0)
fl = EN.fluxes(states, channels, tau)
out_rate = -np.diag(R0)
E_init = states[init]["E_eV"]
E_final = float(np.dot(p_inf, [s["E_eV"] for s in states]))
L += ["## 1. T = 0：厳密解（τ = −R_TT⁻¹ p₀）", "",
      "初期準位の解離に対する深さ −E = %.4f cm⁻¹ = %.6f eV。" % (-E_cm[init], -E_init),
      "放出合計 Σ_kind = %.9f eV、E_init − E_final = %.9f eV、差 %.1e eV。" % (sum(fl.values()), E_init - E_final, sum(fl.values()) - (E_init - E_final)), "",
      "| 種別 | 放出エネルギー [eV] | 光子数（期待値） |", "|---|---|---|"]
nph = {}
for a, b, kind, k, A in channels:
    nph[kind] = nph.get(kind, 0.0) + A * tau[a]
for kind in sorted(fl):
    L.append("| %s | %.9e | %.6f |" % (kind, fl[kind], nph[kind]))
L += ["", "終状態 p_∞：" + "、".join("%s %.6f" % (states[i]["name"], p_inf[i]) for i in range(n) if p_inf[i] > 1e-12),
      "（オルトとパラは結ばれないので、初期準位と同じ核スピン系列の最下準位で止まる。）", "",
      "滞在時間 τ > 0 の準位（経路）。最右列は原子間重力 ⟨−Gm_H²/R⟩ の別勘定（点質量 m_H の Newton 項。電子の交差項 ~1e-3・1PN ~1e-10 は未実装。倍精度の床の 20 桁下なので E には加わらない）：", "",
      "| 準位 | −E [cm⁻¹] | τ [s] | 1/Γ（寿命）[s] | Γ [s⁻¹] | ⟨−Gm_H²/R⟩ [cm⁻¹] |", "|---|---|---|---|---|---|"]
for i in np.argsort(-tau):
    if tau[i] <= 0.0:
        continue
    L.append("| %s | %.4f | %.6e | %.6e | %.6e | %.4e |" % (states[i]["name"], -E_cm[i], tau[i], 1.0 / out_rate[i] if out_rate[i] > 0 else float("inf"), out_rate[i], states[i]["E_grav_cm"]))
L += ["", "放出線（光子数 A τ_src > 1e-9、上位 30 本）。幅：Γ = Γ_src + Γ_dst、Δt = 1/Γ、Δx = c/Γ。最右列は重力による線のずれ Δ(ħω)_grav = E_grav(src) − E_grav(dst)：", "",
      "| 線 | 種別 | ħω [eV] | λ [μm] | A [s⁻¹] | 光子数 | Γ [s⁻¹] | Δt [s] | Δx [m] | Γħ/ħω | Δ(ħω)_grav [eV] |", "|---|---|---|---|---|---|---|---|---|---|---|"]
lines = sorted(((A * tau[a], a, b, kind, A) for a, b, kind, k, A in channels if A * tau[a] > 1e-9), reverse=True)
for nq, a, b, kind, A in lines[:30]:
    hw = states[a]["E_eV"] - states[b]["E_eV"]
    G = out_rate[a] + out_rate[b]
    lam = 2 * np.pi * EN.HBAR_EVS * EN.C / hw * 1e6
    L.append("| (%d,%d)→(%d,%d) | %s | %.6f | %.4f | %.4e | %.6f | %.3e | %.3e | %.3e | %.1e | %+.3e |" % (
        *states[a]["label"], *states[b]["label"], kind, hw, lam, A, nq, G, 1.0 / G, EN.C / G, G * EN.HBAR_EVS / hw,
        states[a]["E_grav_eV"] - states[b]["E_grav_eV"]))
L += ["", "原子間重力の別勘定：⟨−Gm_H²/R⟩ は (0,0) で %.4e cm⁻¹、(%d,%d) で %.4e cm⁻¹。GW の放出 %.3e eV と合わせ、重力はこの系で (i) 準位の 1e-31 cm⁻¹ のずれと (ii) 1e-35 eV の放出としてだけ現れる。" % (
    states[idx[(0, 0)]]["E_grav_cm"], v0q, J0q, states[init]["E_grav_cm"], fl.get("gw", 0.0)), ""]

# ---------- 層 2：温度 ----------
T0 = 2.7255
R2, rows2 = EN.generator(states, channels, T0)
times, hist = EN.evolve(R2, p0, t_max=1e17, n_times=60)
pT = hist[-1]
g_i, g_f = idx[(0, 1 if J0q % 2 else 0)], idx[(0, 3 if J0q % 2 else 2)]
dE = states[g_f]["E_eV"] - states[g_i]["E_eV"]
boltz = (2 * states[g_f]["J"] + 1) / (2 * states[g_i]["J"] + 1) * np.exp(-dE / (EN.K_B_EV * T0))
rmax = max(rows2, key=lambda r: r[5])
L += ["## 2. 層 2：吸収体の温度 T = %.4f K" % T0, "",
      "最下準位からの最低励起 %s→%s：ħω = %.4f eV = %.1f kT。全辺で Planck 占有数が最大なのは近接準位対 (%d,%d)→(%d,%d)（ħω = %.4f cm⁻¹）の ⟨n̄⟩ = %.2f で、これは高励起準位間の偶然の近接。" % (
          states[g_i]["name"], states[g_f]["name"], dE, dE / (EN.K_B_EV * T0), *states[rmax[0]]["label"], *states[rmax[1]]["label"], rmax[7] / S.HARTREE_EV * CM, rmax[5]),
      "t = 1e17 s での占有：%s %.6f、%s %.1e。Boltzmann 比 g_f/g_i e^{−ħω/kT} = %.1e。" % (states[g_i]["name"], pT[g_i], states[g_f]["name"], pT[g_f], boltz),
      "水素原子の 21 cm（n̄ = 39.5、F=1:F=0 = 0.745:0.255）と違い、H₂ の最低遷移は kT の %.0f 倍で、2.7255 K の背景は何も励起しない。定常状態は核スピン系列の最下準位 1 つ。" % (dE / (EN.K_B_EV * T0)), ""]

# ---------- 層 1：方向・重心速度・反跳 ----------
rng = np.random.default_rng(11)
vdrift = np.array([0.0, 0.0, 370e3])
ev_all = []
n_mol = 200
for a in range(n_mol):
    ev_all += [e for e in EN.simulate_atom(states, channels, T0, vdrift, rng, init, n_events=12) if e["hw_lab"] > 0.0 and "2ph" not in e["kind"]]
Nmat = np.array([e["n_lab"] for e in ev_all]); dvec = np.array([e["delta"] for e in ev_all])
b0 = np.linalg.lstsq(Nmat, dvec, rcond=None)[0]
sol = least_squares(lambda bv: dvec - np.array([EN.doppler_delta(nn, bv * EN.C) for nn in Nmat]), x0=b0)
vfit = sol.x * EN.C
kicks = [np.linalg.norm(e["v_after"] - e["v_before"]) for e in ev_all]
rw = [e["rel_width"] for e in ev_all]
L += ["## 3. 層 1：吸収体の静止系に印、重心速度 v₀ = 370 km/s、反跳（量子跳躍モンテカルロ、T = %.4f K）" % T0, "",
      "%d 分子、%d 光子。1 光子あたりの反跳 |Δv| = %.3e〜%.3e m/s（M = 2(m_p+m_e)、ħω/Mc）。線の自然幅 Γħ/ħω = %.1e〜%.1e。" % (
          n_mol, len(ev_all), min(kicks), max(kicks), min(rw), max(rw)),
      "記録 (δ = ω_lab/ω′−1, n̂) だけから厳密な Doppler 模型で v⃗ を復元：v⃗_fit = (%.3f, %.3f, %.3f) m/s、|v⃗_fit − v⃗₀| = %.3f m/s（1 次の線形解からの初期値 |Δ| = %.1f m/s）。" % (
          *vfit, np.linalg.norm(vfit - vdrift), np.linalg.norm(b0 * EN.C - vdrift)),
      "t ≤ 1e17 s に起きた吸収・誘導放出の事象数：%d（最下準位からの励起は ⟨n̄⟩ ~ 1e-135）。" % sum(1 for e in ev_all if "abs" in e["kind"] or "stim" in e["kind"]), ""]

# ---------- 層 3：m 分解（初期準位から下へ到達できる同系列の準位に限る） ----------
sub = [i for i in range(n) if states[i]["E_eV"] <= E_init + 1e-12 and (states[i]["J"] - J0q) % 2 == 0]
remap = {i: k for k, i in enumerate(sub)}
sstates = [states[i] for i in sub]
schannels = [(remap[a], remap[b], kind, k, A) for a, b, kind, k, A in channels if a in remap and b in remap]
Rm, mb, edm, worst = EN.generator_m(sstates, schannels, 0.0)
midx = {s: k for k, s in enumerate(mb)}
p0m = np.zeros(len(mb)); p0m[midx[(remap[init], 0)]] = 1.0
taum, pinfm = EN.exact(Rm, p0m)
ph = EN.line_photons_m(mb, edm, taum)
L += ["## 4. 層 3：m 分解（Wigner–Eckart）、初期 m = 0", "",
      "対象：初期準位以下の同じ核スピン系列 %d 準位、m 副準位 %d、辺 %d 本、和則の最大相対差 %.1e。" % (len(sub), len(mb), len(edm), worst), "",
      "| 線 | 種別 | N₀ : N₊ : N₋（Δm = 0, +1, −1 の光子数）| 偏光度 P(90°) |", "|---|---|---|---|"]
for key in sorted(ph, key=lambda kk: -sum(ph[kk].values()))[:8]:
    Nq = ph[key]; tot = sum(Nq.values())
    if tot < 1e-6:
        continue
    L.append("| (%d,%d)→(%d,%d) | %s | %.4f : %.4f : %.4f | %+.4f |" % (*sstates[key[0]]["label"], *sstates[key[1]]["label"], key[2],
                                                                        Nq.get(0, 0) / tot, Nq.get(1, 0) / tot, Nq.get(-1, 0) / tot, EN.polarization_90(Nq)))
gJ = states[g_i]["J"]
pm_fin = {mm: pinfm[midx[(remap[g_i], mm)]] for mm in range(-gJ, gJ + 1)}
tot = sum(pm_fin.values())
align = sum(mm * mm * p for mm, p in pm_fin.items()) / tot - gJ * (gJ + 1) / 3.0
L += ["", "終状態 %s の m 分布：%s。整列 ⟨m²⟩ − J(J+1)/3 = %+.4f。" % (states[g_i]["name"], "、".join("m=%+d %.4f" % (mm, p / tot) for mm, p in sorted(pm_fin.items())), align),
      "T = 2.7255 K の背景は ⟨n̄⟩ ≈ 0 なので m 副準位間の緩和・整列（水素の 21 cm では双極子 5e-9）は起こらず、初期 m の記憶がそのまま残る。", ""]

out = "\n".join(L) + "\n"
(HERE / ("audit_h2_v%dJ%d.txt" % (v0q, J0q))).write_text(out, encoding="utf-8")
print(out)

# ---------- 図化用に同じ数値を JSON で保存（計算は上と同一、出力形式だけ） ----------
import json  # noqa: E402
dump = dict(
    init=[v0q, J0q], T=T0,
    states=[dict(label=list(s["label"]), E_eV=s["E_eV"], J=s["J"], I=s["I"], E_grav_eV=s["E_grav_eV"], E_grav_cm=s["E_grav_cm"]) for s in states],
    tau=tau.tolist(), p_inf=p_inf.tolist(), out_rate=out_rate.tolist(), fluxes=fl, photons=nph,
    lines=[dict(src=int(a), dst=int(b), kind=kind, A=A, n=A * tau[a], hw_eV=states[a]["E_eV"] - states[b]["E_eV"],
                Gamma=out_rate[a] + out_rate[b]) for a, b, kind, k, A in channels if A * tau[a] > 1e-12],
    layer2=dict(ground=int(g_i), first=int(g_f), p_ground=float(pT[g_i]), p_first=float(pT[g_f]), boltzmann=float(boltz),
                times=times.tolist(), hist_ground=hist[:, g_i].tolist(), hist_init=hist[:, init].tolist()),
    layer1=dict(v0=vdrift.tolist(), vfit=vfit.tolist(), v_lin=(b0 * EN.C).tolist(), n_mol=n_mol,
                events=[dict(n=[float(x) for x in e["n_lab"]], delta=float(e["delta"]), hw_eV=float(e["hw_lab"]), kind=e["kind"],
                             dv=float(np.linalg.norm(e["v_after"] - e["v_before"])), rel_width=float(e["rel_width"])) for e in ev_all]),
    layer3=dict(n_sub=len(sub), n_m=len(mb), lines={"%d,%d->%d,%d:%s" % (*sstates[k[0]]["label"], *sstates[k[1]]["label"], k[2]): {str(q): v for q, v in Nq.items()}
                                                   for k, Nq in ph.items() if sum(Nq.values()) > 1e-6},
                final_m={str(mm): p / tot for mm, p in pm_fin.items()}, align=align),
)
(HERE / ("audit_h2_v%dJ%d.json" % (v0q, J0q))).write_text(json.dumps(dump, indent=1), encoding="utf-8")
