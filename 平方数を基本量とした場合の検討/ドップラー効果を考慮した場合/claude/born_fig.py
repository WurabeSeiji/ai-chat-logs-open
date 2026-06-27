#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1 for the Born-from-localization result (single localized kernel, two conditions).
Generates born_from_localization.png / .svg  (2x2):
 (a) projector+square reproduces base^2 = cos^2 EXACTLY for every N
 (b) general real band-limited base reproduced once N >= its bandwidth
 (c) envelope reading FAILS (1/sin^2 diverges, non-normalizable)
 (d) COMPLEX base: |conv|^2/norm = |psi|^2 = Re^2 + Im^2  (true Born modulus)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def S_N(x, N):
    x = np.asarray(x, float)
    s = np.sin(x)
    return np.where(np.abs(s) < 1e-12, (N+1)/2.0, np.sin((N+1)*x)/(2.0*s))

def conv(phi0_arr, N, base, ng=120001):
    phi = np.linspace(-np.pi/2, np.pi/2, ng); b = base(phi)
    return np.array([np.trapezoid(S_N(p0-phi, N)*b, phi) for p0 in phi0_arr])

def conv_c(phi0_arr, N, base, ng=120001):
    phi = np.linspace(-np.pi/2, np.pi/2, ng); b = base(phi)
    return np.array([np.trapezoid(S_N(p0-phi, N)*b, phi) for p0 in phi0_arr], dtype=complex)

cos   = lambda x: np.cos(x)
multi = lambda x: np.cos(x) + 0.5*np.cos(3*x) - 0.3*np.cos(5*x)
psi_c = lambda x: np.exp(1j*x) + (0.5-0.3j)*np.exp(3j*x) + (0.2+0.4j)*np.exp(-5j*x)

p = np.linspace(-np.pi/2*0.995, np.pi/2*0.995, 400)
xpct = p/np.pi*100   # phase as % of half-wavelength (paper's axis convention)

fig, AX = plt.subplots(2, 2, figsize=(13.5, 9.2))
ax = AX.ravel()

# (a) base = cos : projector+square is EXACT base^2 for every N
axa = ax[0]
axa.plot(xpct, np.cos(p)**2, color="k", lw=2.4, label="base² = cos²φ  (target, Born)")
for N, c in [(1,"#d62728"), (9,"#1f77b4"), (99,"#2ca02c")]:
    P = conv(p, N, cos)**2 / (np.pi/2)**2
    axa.plot(xpct[::12], P[::12], "o", ms=5, color=c, label=f"|S_N*cos|²/norm,  N={N}")
axa.set_title("(a) projector + square  →  base²  (EXACT, every N)", fontsize=10)
axa.set_xlabel("phase φ  (% of half-wavelength)"); axa.set_ylabel("normalized probability")
axa.legend(fontsize=7.5, loc="lower center"); axa.grid(alpha=0.25); axa.set_ylim(-0.05, 1.1)

# (b) general real band-limited base : reproduced once N >= highest mode (=5)
axb = ax[1]
axb.plot(xpct, multi(p)**2, color="k", lw=2.4, label="|psi_base|²  (real, modes 1,3,5)")
for N, c, ok in [(3,"#ff7f0e",False),(5,"#1f77b4",True),(31,"#2ca02c",True)]:
    P = conv(p, N, multi)**2 / (np.pi/2)**2
    lab = f"N={N}  " + ("(covers modes → exact)" if ok else "(misses mode 5)")
    axb.plot(xpct[::10], P[::10], "o" if ok else "x", ms=5, color=c, label=lab)
axb.set_title("(b) general real base: exact once N ≥ its bandwidth", fontsize=10)
axb.set_xlabel("phase φ  (% of half-wavelength)"); axb.set_ylabel("normalized probability")
axb.legend(fontsize=7.5, loc="upper center"); axb.grid(alpha=0.25)

# (c) envelope reading : 1/sin^2 diverges at center -> NOT a probability, != cos^2
axc = ax[2]
pe = np.linspace(-np.pi/2*0.97, np.pi/2*0.97, 1200); xe = pe/np.pi*100
N = 99
env = 0.5/((N+1)*np.sin(pe))**2
env_disp = env/np.max(env[np.abs(pe) > 0.15])   # display scaling only
axc.plot(xe, np.cos(pe)**2, color="k", lw=2.4, label="base² = cos²φ")
axc.plot(xe, np.clip(env_disp, 0, 1.5), color="#9467bd", lw=1.8,
         label="envelope ~ 1/sin²φ  (diverges at 0)")
axc.axvline(0, color="#9467bd", ls=":", lw=1)
axc.set_title("(c) envelope reading FAILS  (1/sin²φ, non-normalizable)", fontsize=10)
axc.set_xlabel("phase φ  (% of half-wavelength)"); axc.set_ylabel("display-scaled")
axc.set_ylim(0, 1.5); axc.legend(fontsize=8, loc="upper center"); axc.grid(alpha=0.25)

# (d) COMPLEX base : |conv|^2/norm = |psi|^2 = Re^2 + Im^2  (genuine modulus, !=Z^2)
axd = ax[3]
mod2 = np.abs(psi_c(p))**2
re2  = np.real(psi_c(p))**2
im2  = np.imag(psi_c(p))**2
axd.plot(xpct, mod2, color="k", lw=2.4, label="|psi_c|² = Re²+Im²  (complex, modes 1,3,−5)")
axd.fill_between(xpct, 0, re2, color="#1f77b4", alpha=0.20, label="Re² part")
axd.fill_between(xpct, re2, re2+im2, color="#d62728", alpha=0.20, label="Im² part")
for N, c, mk in [(3,"#ff7f0e","x"), (5,"#2ca02c","o"), (31,"#7f7f7f","o")]:
    P = np.abs(conv_c(p, N, psi_c))**2 / (np.pi/2)**2
    lab = f"|S_N*psi_c|²/norm, N={N}" + ("" if N>=5 else " (misses 5)")
    axd.plot(xpct[::13], P[::13], mk, ms=5, color=c, label=lab)
axd.set_title("(d) complex base: |Z|² = Re²+Im²  (true Born modulus, ≠ Z²)", fontsize=10)
axd.set_xlabel("phase φ  (% of half-wavelength)"); axd.set_ylabel("normalized probability")
axd.legend(fontsize=7.0, loc="upper center"); axd.grid(alpha=0.25)

fig.suptitle("Born distribution from a single localized odd-harmonic kernel  "
             "(phase-difference observation = convolution + squaring)", fontsize=12, y=1.0)
fig.tight_layout()
fig.savefig("born_from_localization.png", dpi=150, bbox_inches="tight")
fig.savefig("born_from_localization.svg", bbox_inches="tight")
plt.close(fig)
print("saved born_from_localization.png/.svg")
