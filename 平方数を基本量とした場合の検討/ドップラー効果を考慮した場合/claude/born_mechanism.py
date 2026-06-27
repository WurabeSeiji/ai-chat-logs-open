#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 3 (mechanism): WHY the localized kernel reproduces the Born form,
and why EQUAL amplitude is essential.  Generates born_mechanism.png / .svg.

Top row  : equal-amplitude kernel  K = S_N = sum_{m} cos((2m+1)phi)      (= reproducing kernel)
Bottom   : tapered kernel          K~ = sum_{m} 1/(2m+1) cos((2m+1)phi)  (= NOT a reproducing kernel)

Columns:
 (1) the kernel samples psi_base at each phi0 (localized bumps under psi_base)
 (2) reconstruction (K*psi_base)/(pi/2)  vs target psi_base
 (3) Born form |K*psi_base|^2 / norm     vs target |psi_base|^2

Same psi_base as the main paper:  cos(phi) + 0.5 cos(3phi) - 0.3 cos(5phi)  (band up to mode 5).
Numerical fact reproduced here: equal-amp reconstruction is exact (~1e-16),
tapered-amp reconstruction is distorted (max dev ~0.53).
x-axis: phase in % of the half-wavelength interval (edge phi=+-pi/2 -> +-100%).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 9
HALF = np.pi/2
psi_base = lambda x: np.cos(x) + 0.5*np.cos(3*x) - 0.3*np.cos(5*x)

def kernel(x, amps):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    for m, a in enumerate(amps):
        out += a*np.cos((2*m+1)*x)
    return out

EQUAL   = np.ones((N+1)//2)
TAPERED = np.array([1.0/(2*m+1) for m in range((N+1)//2)])

def conv(phi0, amps, ng=200001):
    phi = np.linspace(-HALF, HALF, ng); b = psi_base(phi)
    return np.array([np.trapezoid(kernel(p0-phi, amps)*b, phi) for p0 in phi0])

# fine grids
phi  = np.linspace(-HALF, HALF, 1600)
xpct = phi/HALF*100
phi0 = np.linspace(-HALF*0.98, HALF*0.98, 41)
x0   = phi0/HALF*100
samples = np.array([-0.55, 0.0, 0.45])*HALF      # a few illustrative sampling points

fig, AX = plt.subplots(2, 3, figsize=(15.5, 8.6), sharex=True)

def row(r, amps, color, tag_ok):
    # (1) kernel samples psi_base
    ax = AX[r, 0]
    ax.plot(xpct, psi_base(phi), color="k", lw=2.0, label=r"$\psi_{\rm base}(\varphi)$")
    peak = kernel(np.array([0.0]), amps)[0]
    for s in samples:
        ax.axvline(s/HALF*100, color="#2ca02c", ls=":", lw=1.0, alpha=0.8)
        ax.plot([s/HALF*100], [psi_base(s)], "o", color="k", ms=5, zorder=6)
        bump = kernel(phi - s, amps)/peak*0.55       # scaled localized kernel
        ax.fill_between(xpct, 0, bump, color="#2ca02c", alpha=0.18, lw=0)
        ax.plot(xpct, bump, color="#2ca02c", lw=0.8, alpha=0.7)
    ax.set_ylim(-1.0, 1.9)
    ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
    ax.set_title(("(1) kernel $S_N$ samples $\\psi_{\\rm base}$ at each $\\varphi_0$" if r == 0
                  else "(1') tapered kernel $\\tilde S_N$, amps $1/(2m+1)$"), fontsize=10)

    # (2) reconstruction vs target
    rec = conv(phi0, amps)/(HALF)
    ax = AX[r, 1]
    ax.plot(xpct, psi_base(phi), color="0.6", lw=4.0, solid_capstyle="round",
            label=r"target $\psi_{\rm base}$")
    ax.plot(x0, rec, "o", color=color, ms=4.5, label=r"$(K*\psi_{\rm base})/(\pi/2)$")
    ax.set_ylim(-1.0, 1.9)
    ax.legend(fontsize=8.5, loc="lower center", framealpha=0.95)
    ax.set_title(("(2) reproduces faithfully" if r == 0
                  else "(2') reproduction is distorted"), fontsize=10)

    # (3) Born form
    born = rec**2
    tgt  = psi_base(phi)**2
    ax = AX[r, 2]
    ax.plot(xpct, tgt, color="0.6", lw=4.0, solid_capstyle="round",
            label=r"target $|\psi_{\rm base}|^2$")
    ax.plot(x0, born, "o", color=color, ms=4.5, label=r"$|K*\psi_{\rm base}|^2$ (norm.)")
    ax.set_ylim(-0.05, 0.85)
    ax.legend(fontsize=8.5, loc="lower center", framealpha=0.95)
    ax.set_title(("(3) Born form matches" if r == 0
                  else "(3') does NOT match $|\\psi_{\\rm base}|^2$"), fontsize=10)
    return np.max(np.abs(rec - psi_base(phi0)))

dev_eq  = row(0, EQUAL,   "#1f77b4", True)
dev_tap = row(1, TAPERED, "#d62728", False)

AX[0, 0].set_ylabel("EQUAL amplitude\n= reproducing kernel", fontsize=11,
                    color="#1f77b4", fontweight="bold")
AX[1, 0].set_ylabel("UNEQUAL amplitude\n= not reproducing", fontsize=11,
                    color="#d62728", fontweight="bold")
for c in range(3):
    AX[1, c].set_xlabel("phase  (% of half-wavelength)")

fig.suptitle("Why the localized kernel reproduces the Born form -- equal amplitude is essential\n"
             r"($N=9$,  $\psi_{\rm base}=\cos\varphi+0.5\cos3\varphi-0.3\cos5\varphi$)",
             fontsize=13, y=1.0)
fig.tight_layout()
fig.savefig("born_mechanism.png", dpi=150, bbox_inches="tight")
fig.savefig("born_mechanism.svg", bbox_inches="tight")
plt.close(fig)
print("saved born_mechanism.png/.svg")
print("equal-amp   max|recon - psi| = %.2e   (faithful)" % dev_eq)
print("tapered-amp max|recon - psi| = %.2e   (distorted)" % dev_tap)
