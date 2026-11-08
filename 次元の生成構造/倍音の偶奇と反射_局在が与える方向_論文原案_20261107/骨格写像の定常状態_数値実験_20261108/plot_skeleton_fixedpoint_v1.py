#!/usr/bin/env python3
"""skeleton_fixedpoint_v1 の図化。results/ の JSON と NPZ を読み、figures/ に PNG と SVG を出力する"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent
R = json.loads((HERE / "results/skeleton_fixedpoint_v1.json").read_text(encoding="utf-8"))
D = np.load(HERE / "results/skeleton_fixedpoint_v1.npz")
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
for f in font_manager.findSystemFonts():
    if "NotoSansCJK" in f or "NotoSerifCJK" in f:
        font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "Noto Sans CJK SC", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
NG = R["params"]["NG"]; deg = np.arange(NG) * 360 / NG
eps = R["E3_half_size_rule"]["threshold"]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=160); fig.savefig(FIG / f"{name}.svg")
    plt.close(fig)


# 図1 反復の経過（実グラム固有値）と関係の変化
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for k, (name, lab) in enumerate((("complex", "複素（共通軸）"), ("quaternion", "四元数（4軸）"))):
    h = R["E1"][name]["history"]
    st = [0] + [x["step"] for x in h]
    ev = np.array([R["E1"][name]["real_gram_eig_initial"]] + [x["eig"] for x in h])
    for c in range(4):
        ax[k].plot(st, ev[:, c], "o-", ms=3)
    ax[k].axhline(eps, ls="--", c="gray", lw=0.8)
    ax[k].set_title(f"{lab}：実グラム固有値"); ax[k].set_xlabel("反復回数"); ax[k].set_ylim(-0.1, 4.1)
    ax[2].semilogy([x["step"] for x in h], [max(x["rel_change"], 1e-17) for x in h], "o-", ms=3, label=lab)
ax[2].set_title("関係 Y の1回あたり相対変化"); ax[2].set_xlabel("反復回数"); ax[2].legend()
save(fig, "fig1_iteration")

# 図2 密度：初期の各波と固定点の共通密度
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
for k, (name, lab) in enumerate((("complex", "複素（共通軸）"), ("quaternion", "四元数（4軸）"))):
    Q0 = D[f"E1_{name}_Q0"]; Qf = D[f"E1_{name}_Qfixed"]
    for i in range(4):
        ax[k].plot(deg, (Q0[i] ** 2).sum(-1), lw=0.9, alpha=0.6, label=f"初期 波{i+1}")
    ax[k].plot(deg, (Qf[0] ** 2).sum(-1), "k", lw=2.0, label="固定点（4波共通）")
    ax[k].set_title(f"{lab}：密度 |q_i(ψ)|²"); ax[k].legend(fontsize=8, ncol=5)
ax[1].set_xlabel("弦の位相 ψ（度）")
save(fig, "fig2_density")

# 図3 集団：共通密度の山の数と実グラム固有値の分布
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
w = 0.38
for k, (name, lab) in enumerate((("complex", "複素"), ("quaternion", "四元数"))):
    npk = D[f"E4_{name}_npeaks"]
    cs = np.arange(1, 7)
    ax[0].bar(cs + (k - 0.5) * w, [np.sum(npk == c) for c in cs], w, label=lab)
    ev = D[f"E4_{name}_eigs"].ravel()
    ax[1].hist(ev, bins=80, range=(-0.05, 4.05), alpha=0.6, label=lab)
ax[0].set_title("固定点の共通密度の山の数（200配置）"); ax[0].set_xlabel("山の数"); ax[0].legend()
ax[1].axvline(eps, ls="--", c="gray", lw=0.8)
ax[1].set_title("固定点の実グラム固有値の分布"); ax[1].set_xlabel("固有値"); ax[1].legend()
save(fig, "fig3_ensemble")
print("figures written to", FIG)
