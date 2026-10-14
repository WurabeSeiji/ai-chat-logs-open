# -*- coding: utf-8 -*-
"""図化：光円錐読出し（円周）"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["svg.fonttype"] = "none"
d = json.load(open("results_lightcone_readout_s1.json", encoding="utf-8"))
S = d["scan"]; lr = np.array([s["ln_ratio"] for s in S])
v = np.array([s["v"] for s in S]); ig = np.array([s["inv_gamma"] for s in S])
fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
ax = axs[0]; x = np.linspace(-3, 3, 400)
ax.plot(x, np.tanh(x), color="0.7", lw=3, label="tanh(ln(a/b))"); ax.plot(x, 1/np.cosh(x), color="0.85", lw=3, label="sech(ln(a/b))")
ax.plot(lr, v, "o", ms=3.5, color="#c0392b", label="v = k/R（場から）"); ax.plot(lr, ig, "s", ms=3.5, color="#1f77b4", label="1/γ = m/R（場から）")
ax.set_xlabel("ln(a/b)（± 回転項の重み比の対数）"); ax.set_title("速さと 1/γ"); ax.legend(fontsize=8); ax.grid(alpha=.3)
ax = axs[1]; t = np.linspace(-np.pi/2, np.pi/2, 200)
ax.plot(np.sin(t), np.cos(t), color="0.75", lw=2); ax.plot(v, ig, "o", ms=3.5, color="#6a51a3")
ax.plot([1, -1], [0, 0], "*", ms=12, color="#e67e22"); ax.text(0.98, 0.06, "光子的", ha="right", fontsize=9)
ax.plot([0], [1], "*", ms=12, color="#2ca02c"); ax.text(0.04, 1.02, "静止（a=b）", fontsize=9)
ax.set_aspect("equal"); ax.set_xlabel("v = k/R"); ax.set_ylabel("1/γ = m/R"); ax.set_title("v² + (1/γ)² = 1（光円錐）"); ax.grid(alpha=.3)
ax = axs[2]; A = d["velocity_addition"]; i = np.arange(len(A))
ax.plot(i, [a["v_galilean"] for a in A], "x", ms=9, color="0.5", label="単純和 v1+v2")
ax.plot(i, [a["v_relativistic_formula"] for a in A], "o", ms=11, mfc="none", color="#1f77b4", label="相対論の合成則")
ax.plot(i, [a["v_composed_field"] for a in A], ".", ms=10, color="#c0392b", label="重み比の積から場で読出し")
ax.axhline(1, color="0.7", ls=":"); ax.axhline(-1, color="0.7", ls=":")
ax.set_xticks(i); ax.set_xticklabels(["%.2f\n%+.2f" % (a["v1"], a["v2"]) for a in A], fontsize=8)
ax.set_xlabel("v1 / v2"); ax.set_title("速度合成"); ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.suptitle("円周：Σx_n² = R²（本当の二乗）を光円錐として読む", fontsize=12)
fig.tight_layout(); fig.savefig("fig_lightcone_readout_s1.png", dpi=160); fig.savefig("fig_lightcone_readout_s1.svg"); print("ok")
