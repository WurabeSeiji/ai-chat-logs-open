#!/usr/bin/env python3
"""回帰テスト：エンジン＋水素の部品が、固定一式と ../exchange_general_20261006/ の層 1〜3 の数値に戻ること。
出力：regression_hydrogen.md
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import engine as EN  # noqa: E402
import system_hydrogen as SH  # noqa: E402
X = SH.X
EV = X.EV

L = ["# 回帰テスト：エンジン＋水素部品 vs 固定一式・層 1〜3", ""]
states, channels, init = SH.build()
p0 = np.zeros(len(states)); p0[init] = 1.0

# 1. 固定一式（T = 0）
R, rows = EN.generator(states, channels, 0.0)
tau, p_inf = EN.exact(R, p0)
fl = EN.fluxes(states, channels, tau)
res_f = X.run(X.W_GW)
emit_f = sum(res_f["flux_E"].values()) * EV
emit_e = sum(fl.values())
L.append("## 1. T = 0：固定一式との比較")
L.append("")
L.append("| 量 | 固定一式 | エンジン | 差 |")
L.append("|---|---|---|---|")
L.append("| 放出合計 [eV] | %.9f | %.9f | %.1e |" % (emit_f, emit_e, emit_e - emit_f))
for k in sorted(fl):
    L.append("| 放出 %s [eV] | %.9e | %.9e | %.1e（相対） |" % (k, res_f["flux_E"][k] * EV, fl[k], fl[k] / (res_f["flux_E"][k] * EV) - 1))
L.append("| τ の最大相対差 | | | %.1e |" % max(abs(tau[i] / res_f["tau"][i] - 1) for i in range(len(states)) if res_f["tau"][i] > 0))
L.append("| p_∞ の最大差 | | | %.1e |" % np.max(np.abs(p_inf - res_f["p_inf"])))
L.append("| 生成行列 R の最大差 [s⁻¹] | | | %.1e |" % np.max(np.abs(R - res_f["R"])))
L.append("")

# 2. 層 2（T = 2.7255 K）
R2, _ = EN.generator(states, channels, 2.7255)
times, hist = EN.evolve(R2, p0)
pT = hist[-1]
i1, i0 = X.IDX[(1, 0, 1, 1)], X.IDX[(1, 0, 1, -1)]
L.append("## 2. 層 2（T = 2.7255 K）：absorber_temperature.md の値（F=1 0.745286、F=0 0.254714、比 2.925977）")
L.append("")
L.append("エンジン：F=1 %.6f、F=0 %.6f、比 %.6f、理論 3e^{−ħω/kT} = %.6f" % (pT[i1], pT[i0], pT[i1] / pT[i0], 3 * np.exp(-(states[i1]["E_eV"] - states[i0]["E_eV"]) / (EN.K_B_EV * 2.7255))))
L.append("")

# 3. 層 3（整列）
L.append("## 3. 層 3：absorber_alignment.md の値（等方 0、双極子のみ +5.074e-09／偏光 −3.806e-09、四重極込み +2.506e-08／−1.880e-08）")
L.append("")
Rm0, mb, edm, worst = EN.generator_m(states, channels, 0.0)
L.append("m 分解 %d 状態、辺 %d 本、和則の最大相対差 %.1e" % (len(mb), len(edm), worst))
midx = {s: k for k, s in enumerate(mb)}
p0m = np.zeros(len(mb)); p0m[midx[(init, 0)]] = 1.0
taum, pinfm = EN.exact(Rm0, p0m)
ph = EN.line_photons_m(mb, edm, taum)
key = (init, X.IDX[(1, 0, 1, -1)], "e1")
tot = sum(ph[key].values())
L.append("m_F = 0 の 5p₃/₂F1→1s F0：N₀:N₊:N₋ = %.4f:%.4f:%.4f、偏光度 %.4f（期待 +1）" % (ph[key][0] / tot, ph[key][1] / tot, ph[key][-1] / tot, EN.polarization_90(ph[key])))
for name, beta, a2 in (("等方", 0.0, 0.0), ("双極子", 370e3 / EN.C, 0.0), ("双極子＋四重極", 370e3 / EN.C, 4e-6)):
    Rm, mb, edm, _ = EN.generator_m(states, channels, 2.7255, beta, a2)
    _, histm = EN.evolve(Rm, p0m, n_times=2, rtol=1e-12)
    pm = {mm: histm[-1][midx[(i1, mm)]] for mm in (-1, 0, 1)}
    A2 = (pm[1] + pm[-1] - 2 * pm[0]) / sum(pm.values())
    P90 = (2 * pm[0] - pm[1] - pm[-1]) / (2 * pm[0] + pm[1] + pm[-1])
    L.append("- %s：整列 %+.3e、21 cm 偏光度 %+.3e" % (name, A2, P90))
L.append("")

# 4. 層 1（モンテカルロ、集団の速度読み出し）
L.append("## 4. 層 1：absorber_direction_recoil.md（集団 200 原子で v⃗ の差 0.09 m/s）")
L.append("")
rng = np.random.default_rng(7)
v0 = np.array([0.0, 0.0, 370e3])
ev_all = []
for a in range(200):
    ev_all += [e for e in EN.simulate_atom(states, channels, 2.7255, v0, rng, init, n_events=8) if e["hw_lab"] > 1.0 and "2ph" not in e["kind"]]
from scipy.optimize import least_squares  # noqa: E402
Nmat = np.array([e["n_lab"] for e in ev_all]); dvec = np.array([e["delta"] for e in ev_all])
b0 = np.linalg.lstsq(Nmat, dvec, rcond=None)[0]
w = np.array([e["rel_width"] for e in ev_all])
sol = least_squares(lambda bv: (dvec - np.array([EN.doppler_delta(n, bv * EN.C) for n in Nmat])) / w, x0=b0, loss="cauchy", f_scale=1.0)
vfit = sol.x * EN.C
L.append("エンジンの MC（別の乱数列）：%d 光子、推定 v⃗ = (%.2f, %.2f, %.2f) m/s、差 %.2f m/s（目安 0.5 m/s）" % (len(ev_all), *vfit, np.linalg.norm(vfit - v0)))
L.append("")
out = "\n".join(L) + "\n"
(HERE / "regression_hydrogen.md").write_text(out, encoding="utf-8")
print(out)
