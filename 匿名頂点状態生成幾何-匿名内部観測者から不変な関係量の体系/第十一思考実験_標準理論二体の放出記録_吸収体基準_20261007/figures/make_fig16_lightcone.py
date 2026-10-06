#!/usr/bin/env python3
"""図 16：吸収体は光円錐の null 面の写像（模式図、データなし）。
左：放出事象を頂点とする未来光円錐と過去光円錐。吸収体の静止系の時刻面で切れば円（球）、二体の静止系（速度 β）の時刻面で切れば
    r(θ) = t₀/(1 − β cos θ) の楕円。これが層 1 の Doppler の幾何。過去円錐の断面は最終散乱面（入射する背景 2.7255 K）。
右：二体の世界線から出た各事象の記録が、未来の無限遠（地平線）の断面 (u, n̂) に着地する。位相はヌル測地線に沿って不変。
出力：fig16_lightcone_absorber.png
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150, "savefig.dpi": 150})

fig = plt.figure(figsize=(12.0, 5.6))

# ---------------- left: 2+1 light cone ----------------
ax = fig.add_subplot(1, 2, 1, projection="3d")
th = np.linspace(0, 2 * np.pi, 160)
tmax = 3.0
tt = np.linspace(0, tmax, 13)
TH, TT = np.meshgrid(th, tt)
ax.plot_wireframe(TT * np.cos(TH), TT * np.sin(TH), TT, color="0.78", lw=0.4, rstride=1, cstride=16)
tp = np.linspace(0, 1.7, 8)
THp, TP = np.meshgrid(th, tp)
ax.plot_wireframe(TP * np.cos(THp), TP * np.sin(THp), -TP, color="0.86", lw=0.4, rstride=1, cstride=16)
t0, beta = 1.7, 0.33
ax.plot(t0 * np.cos(th), t0 * np.sin(th), t0 * np.ones_like(th), color="C0", lw=2.2)
r = t0 / (1 - beta * np.cos(th))
ax.plot(r * np.cos(th), r * np.sin(th), t0 + beta * r * np.cos(th), color="C3", lw=2.2)
tc = 1.35
ax.plot(tc * np.cos(th), tc * np.sin(th), -tc * np.ones_like(th), color="C1", lw=2.2)
for ang in (0.5, 2.5, 4.4):
    s = np.linspace(0, tmax, 2)
    ax.plot(s * np.cos(ang), s * np.sin(ang), s, color="k", lw=0.8)
    for sk in np.linspace(0.35, 2.75, 7):
        ax.plot([sk * np.cos(ang)], [sk * np.sin(ang)], [sk], marker="o", color="k", ms=2.2)
ax.scatter([0], [0], [0], color="k", s=28, zorder=6)
ax.text(0.15, -0.9, -0.05, "emission event", fontsize=8)
ax.text(-2.9, 2.2, 3.05, "future null cone", fontsize=8, color="0.35")
ax.text(1.4, -2.2, -1.75, "past null cone", fontsize=8, color="0.35")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("t")
ax.set_xlim(-3.1, 3.1); ax.set_ylim(-3.1, 3.1); ax.set_zlim(-1.8, 3.2)
ax.view_init(elev=24, azim=-48)
ax.set_box_aspect((1, 1, 0.85))
from matplotlib.lines import Line2D  # noqa: E402
handles = [Line2D([], [], color="C0", lw=2.2, label="cut in the absorber rest frame, t = t₀: a sphere"),
           Line2D([], [], color="C3", lw=2.2, label="cut in the two-body rest frame (β): r(θ) = t₀/(1−β cosθ), an ellipsoid"),
           Line2D([], [], color="C1", lw=2.2, label="past cone ∩ last scattering: incoming background 2.7255 K"),
           Line2D([], [], color="k", lw=0.8, marker="o", ms=2.2, label="null generators; phase marks keep their spacing (dφ = 0)")]
ax.legend(handles=handles, loc="upper center", fontsize=7.4, frameon=False, bbox_to_anchor=(0.5, 0.02), ncol=1)
ax.set_title("the absorber is a cut of the future null cone", fontsize=10, pad=2)

# ---------------- right: records landing on future null infinity ----------------
ax2 = fig.add_subplot(1, 2, 2)
ax2.set_xlim(-0.62, 1.55); ax2.set_ylim(-1.25, 1.18); ax2.set_aspect("equal"); ax2.axis("off")
ax2.plot([0, 0], [-1.0, 1.0], color="k", lw=2.2)
ax2.text(-0.03, 1.03, "two-body worldline (r = 0)", fontsize=8, ha="right", va="bottom")
ax2.plot([1.0, 0.0], [0.0, 1.0], color="C0", lw=2.6)
ax2.annotate("future null infinity / horizon = absorber\ncuts (u, n̂) ∈ ℝ_u × S²", xy=(0.22, 0.78), xytext=(0.46, 0.93), fontsize=7.6, color="C0",
             ha="left", va="center", arrowprops=dict(arrowstyle="-", color="C0", lw=0.7))
ax2.plot([1.0, 0.0], [0.0, -1.0], color="C1", lw=2.6)
ax2.annotate("past null infinity: incoming background\n(last scattering, 2.7255 K)", xy=(0.22, -0.78), xytext=(0.46, -0.93), fontsize=7.6, color="C1",
             ha="left", va="center", arrowprops=dict(arrowstyle="-", color="C1", lw=0.7))
events = ((-0.78, "event 1", "u₁"), (-0.36, "event 2", "u₂"), (0.08, "event 3", "u₃"))
for te, lab, u in events:
    s = (1 - te) / 2
    ax2.plot([0, s], [te, te + s], color="C0", lw=1.0)
    ax2.scatter([0], [te], color="k", s=18, zorder=5)
    ax2.scatter([s], [te + s], color="C0", s=24, zorder=5)
    ax2.text(-0.03, te, lab, fontsize=7.6, ha="right", va="center")
    ax2.text(s + 0.035, te + s - 0.01, "record (%s, n̂, ω, pol)" % u, fontsize=7.2, va="center", ha="left")
for ta in (-0.9, -0.56, -0.14):
    s = (1 + ta) / 2
    ax2.plot([s, 0], [ta - s, ta], color="C1", lw=0.8, ls="--")
ax2.annotate("", xy=(-0.12, 0.62), xytext=(-0.12, 0.26), arrowprops=dict(arrowstyle="-|>", lw=1.1, color="0.3"))
ax2.text(-0.16, 0.44, "between events:\nno readable change\n(stationary, phase locked)", fontsize=7.4, color="0.3", ha="right", va="center")
ax2.text(-0.6, -1.17, "Doppler = clock mismatch at the two ends of a null ray, not a change along it.\nEach record is the emission phase relation, carried unchanged.", fontsize=7.4, color="0.3", ha="left", va="bottom")
ax2.set_title("every emission event lands on its own cut of the absorber", fontsize=10)

fig.tight_layout()
fig.savefig(HERE / "fig16_lightcone_absorber.png")
print("saved", HERE / "fig16_lightcone_absorber.png")
