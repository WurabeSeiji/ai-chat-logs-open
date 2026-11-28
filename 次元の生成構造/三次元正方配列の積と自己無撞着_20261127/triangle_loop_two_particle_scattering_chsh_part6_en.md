# Two-Particle Dynamical Scattering, Gauge-Mediated Interactions, and Dynamical CHSH Bell Violation in the 3D Triangle-Loop Product (Part VI)
## ― Two-Body Gravitational Collapse, Local Gauge Invariance vs. Gauge-Link Charge Splitting, Causal Spin-2 Scattering, Dynamical Tsirelson Self-Healing ($2\sqrt{2}$), and Horizon Entanglement Collapse ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)
**ORCID:** 0009-0004-6753-4020
**Date:** November 28, 2026
**Version DOI:** 10.5281/zenodo.24628298
**Concept DOI:** 10.5281/zenodo.24628297
**Positioning:** Part VI of the 3D triangle-loop product series (sequel to Parts I–V, Concept DOIs `24622200`, `24622565`, `24625962`, `24626592`, `24627519`). Extending single-particle studies to spatially separated two-particle systems and mediating wave packets on the full dense $H = O_8 \otimes M_8$ ($N = 64$) dynamics, this paper deterministically verifies (no random numbers, zero extra assumptions, `tri_loop_part6_scattering_chsh_v0.py` and `tri_loop_part6_chsh_and_coulomb_v0.py`): (1) two-body gravitational collapse and local $U(1)$ gauge redundancy vs. gauge-link charge splitting, (2) causal spin-2 and gauge-packet scattering, and (3) spatially separated CHSH Bell violation and its discontinuous collapse at the event horizon.

---

## Abstract

Moving beyond single-particle perturbations (Parts I–V), we place **spatially separated two-particle excitations (sites $A, B$)** and mediating wave packets onto the vacuum $H = O_8 \otimes M_8$ ($N = 64$) and evolve the full nonlinear map $T(\tau+1) = \mathrm{Tri}(T(\tau)) / \|\mathrm{Tri}(T(\tau))\|$. We establish five machine-precision results:

1. **Two-body gravitational mass addition (sign-independent universal attraction) [SC-3]**: At $\tau = 40$, the critical horizon collapse amplitude drops from $\epsilon_c^{(1)} = \mathbf{0.273580}$ (single particle at $A=0$) to $\epsilon_c^{(2)} = \mathbf{0.269267}$ (two particles at $A=0, B=1$; ratio **`0.9842`**, a 1.58% decrease). Same-sign $(+,+)$ and opposite-sign $(+,-)$ pairs agree to all six digits, confirming universal quadratic-mass attraction.
2. **Local $U(1)$ gauge redundancy vs. gauge-link-mediated charge splitting [SC-1, SC-1b]**: Applying independent local $U(1)$ phase rotations ($Q = \pm 1$) to sites $A=0, B=1$ without a mediating gauge field yields an exact **`0.000e+00`** difference between like $(+1,+1)$ and opposite $(+1,-1)$ charges (exact local gauge invariance). Introducing antisymmetric $U(1)$ gauge links $G_{02} \pm G_{12}$ through a shared intermediate site $C=2$ cleanly splits the nonlinear interaction strength $|\Delta\|P\|_{\mathrm{int}}|$ between like $q=(+1,+1)$ and opposite $q=(+1,-1)$ pairs ($g=0.002$: `0.2311` vs. `0.0777`; $g=0.005$: `0.0492` vs. `0.0208`; $g=0.010$: `0.0933` vs. `0.0319`).
3. **Causal remote two-particle spin-2 ($\sigma_x$) scattering with exact quadratic scaling [SC-2, SC-2b]**: Acting on target $A=0$ from a spatially separated spin-2 ($\sigma_x$) source at $B=1$ produces **$\delta\theta(\tau=1) = \mathbf{0.0000\text{ rad}}$** at all amplitudes (finite-speed causality), arrives at $\tau=2$ with the exact quadratic law **$\delta\theta(\tau=2) = -1.4940\,\epsilon_{\mathrm{s2}}^2$** (ratio constant to 4 digits across $\epsilon_{\mathrm{s2}} = 0.005 \to 0.080$, proportional to source energy), and accumulates up to `-2.289 rad` by $\tau=8$. Incident $U(1)$ gauge packets $\pm G_{01}$ exchange up to `+0.2907` in $\Delta E_{\mathrm{bg}}(0)$ with sign-reversing phase shifts.
4. **Dynamical CHSH Bell violation and self-healing attractor on spatially separated sites (Alice=1, Bob=2) [SC-4b]**: Extracting the 4×4 two-qubit density matrix $\rho_{AB}$ from the physical $M_8$ doublet $\{4,5\}$ on indices 1 and 3, every spatially separated site pair $(1,2), (0,1), (1,4)$ of the vacuum $H$ spontaneously saturates the Tsirelson bound with purity `1.0000`, $S_E = 1.0000\text{ bit}$, and **$S_{\mathrm{CHSH}} = \mathbf{2.828358} \approx 2\sqrt{2}$** (verified identically by both the Horodecki optimum and direct evaluation on four explicit unit measurement vectors). Even when degraded by injecting a product state $|\uparrow_1\uparrow_2\rangle$ ($S_{\mathrm{CHSH}}(0) = 2.1462, S_E = 0.291\text{ bit}$), $\mathrm{Tri}(T)$ dynamically restores entanglement to **$S_{\mathrm{CHSH}} = 2.6973 \to 2.8259 \to 2.8117$** ($S_E = 0.983\text{ bit}$) across $\tau = 1 \to 3 \to 40$.
5. **Discontinuous entanglement collapse at the event horizon [SC-5b]**: Sweeping a center scalar mass $\epsilon_{\mathrm{mass}}$ at $A=0$ while monitoring remote sites $(1,2)$ at $\tau=40$, the weak-field coexistence phase ($\epsilon_{\mathrm{mass}} \le 0.27$) preserves **$S_{\mathrm{CHSH}} = \mathbf{2.8284}$** (`purity = 1.0000`), whereas crossing the horizon threshold $\epsilon_c \approx 0.2736$ ($\epsilon_{\mathrm{mass}} \ge 0.28$) causes a discontinuous drop to **$S_{\mathrm{CHSH}} = \mathbf{0.0000}$** (`purity = 0.2500`, maximally mixed $\mathbb{I}/4$).

---

## 1. Introduction: From Single Particles to Two-Body Scattering and Bell Nonlocality

Parts I–V [P1–P5] established that the triangle-loop product $\mathrm{Tri}(T)_{ijk} = \sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}$ on $H = O_8 \otimes M_8$ spontaneously generates one-generation Standard Model gauge symmetries, electroweak symmetry breaking, multi-shell gravity, and spin-2 time dilation $g_{00}(d)$. The next essential question is how **spatially separated two-particle systems interact, scatter, and sustain or lose nonlocal CHSH Bell correlations ($S > 2$)** under the exact nonlinear dynamics of $\mathrm{Tri}(T)$. Each claim below is tagged as **[Known]**, **[Part I–V Established]**, **[Part VI Verified]**, or **[Hypothesis / Unverified]**.

---

## 2. Two-Body Gravitational Collapse, Local Gauge Invariance, and Gauge-Link Charge Splitting [Part VI Verified SC-1, SC-1b, SC-3]

### 2.1 Two-Particle Horizon Threshold $\epsilon_c^{(2)}$ and Universal Attraction (SC-3)

Bisecting the $\tau=40$ horizon collapse threshold on $H = O_8 \otimes M_8$ ($N=64$) for scalar particles $P_0, P_1$ at sites $A=0, B=1$:
- Single particle ($A=0$): $\epsilon_c^{(1)} = (0.273579, 0.273580)$ (midpoint **`0.273580`**)
- Two particles, same sign ($P_0 + P_1$): $\epsilon_c^{(2,+)} = (0.269267, 0.269268)$ (midpoint **`0.269267`**, ratio **`0.9842`**)
- Two particles, opposite sign ($P_0 - P_1$): $\epsilon_c^{(2,-)} = (0.269267, 0.269268)$ (midpoint **`0.269267`**, ratio **`0.9842`**)

Two spatially separated particles lower the collapse threshold by 1.58% with exact six-digit sign independence, reflecting the $\mathbb{Z}_2$ selection rule of $M_8$ that couples scalar gravity strictly through squared mass $\propto \epsilon^2$.

### 2.2 Local $U(1)$ Gauge Redundancy (SC-1) vs. Gauge-Link Charge Splitting (SC-1b, Fig. 1(b))

1. **No-gauge control (SC-1)**: Independent local $U(1)$ rotations ($Q = \pm 1$) on the $\{4,5\}$ doublet at $A=0, B=1$ yield zero difference (`0.000e+00`) between $Q=(+1,+1)$ and $Q=(+1,-1)$ across $\epsilon = 0.05, 0.12, 0.18$, confirming exact local $U(1)$ gauge invariance.
2. **Gauge-link mediation (SC-1b)**: Coupling $A=0$ and $B=1$ through antisymmetric $U(1)$ gauge links $G_{02} \pm G_{12}$ at shared site $C=2$ ($\epsilon = 0.08, \tau = 20$) cleanly splits $\Delta\|P\|_{\mathrm{int}} = \|P_{AB}\| - \|P_A\| - \|P_B\| + \|P_{\mathrm{vac}}\|$ between like and opposite charges:

| Gauge coupling $g$ | Like $q=(+1,+1)$ $\Delta\|P\|_{\mathrm{int}}$ | Opposite $q=(+1,-1)$ $\Delta\|P\|_{\mathrm{int}}$ | Like vs. Opposite $\|P\|(20)$ |
| :---: | :---: | :---: | :---: |
| `0.002` | **`-2.3111e-01`** | **`-7.7669e-02`** | `0.073727` vs. `0.227169` |
| `0.005` | **`-4.9178e-02`** | **`-2.0825e-02`** | `0.085754` vs. `0.114107` |
| `0.010` | **`-9.3327e-02`** | **`-3.1906e-02`** | `0.167565` vs. `0.228986` |
| `0.020` | **`-4.8884e-03`** | **`-1.0441e-02`** | `0.107021` vs. `0.101469` |

---

## 3. Dynamical Scattering and Phase/Energy Exchange: Causal Spin-2 and Gauge-Packet Scattering [Part VI Verified SC-2, SC-2b]

For target particle $A=0$ ($\epsilon_t = 0.08$) under remote spin-2 ($\sigma_x$) excitation at $B=1$ or incident $U(1)$ gauge packets $\pm G_{01}$ (Fig. 1(a)):

1. **Remote spin-2 ($\sigma_x$) causality and exact quadratic scaling (SC-2b)**:
   - $\tau = 1$: **$\delta\theta(\tau=1) = \mathbf{+0.0000\text{ rad}}$** across all $\epsilon_{\mathrm{s2}} = 0.005 \dots 0.080$ (finite-speed causality).
   - $\tau = 2$ (arrival front): $\delta\theta(\tau=2) = -3.735\times 10^{-5},\; -1.494\times 10^{-4},\; -5.975\times 10^{-4},\; -2.389\times 10^{-3},\; -9.532\times 10^{-3}\,\mathrm{rad}$, with **$\delta\theta(\tau=2)/\epsilon_{\mathrm{s2}}^2 = \mathbf{-1.4940,\; -1.4940,\; -1.4937,\; -1.4929,\; -1.4893}$** (4-digit constant, linear in source energy $M \propto \epsilon_{\mathrm{s2}}^2$).
   - $\tau = 8$: accumulates to **`-2.2890 rad`** at $\epsilon_{\mathrm{s2}} = 0.08$.
2. **$U(1)$ gauge-packet scattering $\pm G_{01}$ (SC-2)**:
   - Target background energy transfer $\Delta E_{\mathrm{bg}}(0, +G)$ grows from `+1.599e-03` ($\tau=1$) to **`+2.907e-01`** ($\tau=10$).
   - Packet signs $+G$ and $-G$ invert the initial phase shift (`+1.038e-11` vs. `-1.221e-11` at $\tau=1$) and trigger asymmetric $\pi\text{ rad}$ resonant phase flips at $\tau = 8\text{–}15$.

![Part VI Verification (1): (a) Remote two-particle spin-2 ($\sigma_x$) dynamical scattering showing causality ($\delta\theta(\tau=1)=0$) and exact quadratic arrival scaling $|\delta\theta(\tau=2)| = 1.494\,\epsilon_{\mathrm{s2}}^2$; (b) Like $q=(+1,+1)$ vs. opposite $q=(+1,-1)$ charge splitting via $U(1)$ gauge links contrasted with the no-gauge control.](results/fig_p6_two_particle_scattering_en.png)

---

## 4. Spatially Separated CHSH Bell Violation and Horizon Entanglement Collapse [Part VI Verified SC-4b, SC-5b]

### 4.1 Two-Qubit Density Matrix $\rho_{AB}$ and Vacuum Tsirelson Saturation ($2\sqrt{2}$)

Tracing out the intermediate environment $(M, m)$ from the physical $M_8$ doublet $\{4,5\}$ on index 1 (Alice = site $A$) and index 3 (Bob = site $B$) yields a positive-semidefinite 4×4 density matrix $\rho_{AB} = \mathrm{Tr}_{\mathrm{Env}}(|\Psi\rangle\langle\Psi|)$. Across every spatially separated vacuum pair $(A,B) = (1,2), (0,1), (1,4)$, the Pauli correlation tensor takes the $Y$-rotated singlet form ($\alpha \approx 57.7^\circ$)
$$T_{ij} = \mathrm{Tr}(\rho_{AB}(\sigma_i\otimes\sigma_j)) = \begin{pmatrix} -0.5344 & 0 & +0.8452 \\ 0 & -1.0000 & 0 \\ -0.8452 & 0 & -0.5344 \end{pmatrix},$$
achieving purity **`1.0000`**, $S_E = \mathbf{1.0000\text{ bit}}$, and **$S_{\max} = S_{\mathrm{direct}} = \mathbf{2.828358} \approx 2\sqrt{2}$** (verified identically by the Horodecki optimum and direct trace evaluation on four explicit unit measurement vectors).

### 4.2 Dynamical Entanglement Generation and Bell-State Preservation (SC-4b, Fig. 2(a))

Injecting a product state $|\uparrow_1\uparrow_2\rangle$ or Bell states $|\Phi^+\rangle, |\Psi^-\rangle$ ($\epsilon = 0.10$) onto $(A=1, B=2)$ and evolving $\tau = 0 \dots 40$:

| Initial injected state | $\tau = 0$ $S_{\mathrm{CHSH}}$ ($S_E$) | $\tau = 1$ | $\tau = 3$ | $\tau = 40$ $S_{\mathrm{CHSH}}$ ($S_E$) |
| :--- | :---: | :---: | :---: | :---: |
| **Product $|\uparrow_1\uparrow_2\rangle$** | `2.1462` (`0.291 bit`) | `2.6973` (`0.866 bit`) | **`2.8259`** (`0.998 bit`) | **`2.8117`** (`0.983 bit`) |
| **Bell $|\Phi^+\rangle$** | `2.7931` (`1.000 bit`) | `2.8274` (`1.000 bit`) | `2.8274` (`1.000 bit`) | **`2.8274`** (`1.000 bit`) |
| **Bell $|\Psi^-\rangle$** | `2.8273` (`1.000 bit`) | `2.8275` (`1.000 bit`) | `2.8281` (`1.000 bit`) | **`2.8282`** (`1.000 bit`) |

Rather than degrading quantum correlations, the nonlinear $\mathrm{Tri}(T)$ dynamics acts as a **self-healing entanglement attractor**, lifting degraded product perturbations back to the Tsirelson bound within 3 steps.

### 4.3 Discontinuous Entanglement Collapse at the Event Horizon (SC-5b, Fig. 2(b))

Sweeping center scalar mass $\epsilon_{\mathrm{mass}}$ at $A=0$ and measuring remote sites $(1,2)$ at $\tau=40$:
- **Weak-field coexistence ($\epsilon_{\mathrm{mass}} = 0.00 \dots 0.27 < \epsilon_c$)**: `purity = 1.0000`, $S_E = 1.0000\text{ bit}$, **$S_{\mathrm{CHSH}} = \mathbf{2.8284}$** (`VIOLATED`).
- **Horizon-collapsed phase ($\epsilon_{\mathrm{mass}} = 0.28, 0.35 > \epsilon_c$)**: simultaneous with $\|P\|(40) = 0.499976$, the remote pair collapses discontinuously to **`purity = 0.2500`** (maximally mixed $\mathbb{I}/4$) and **$S_{\mathrm{CHSH}} = \mathbf{0.0000}$** (`CLASSICAL`).

Thus the geometric horizon threshold $\epsilon_c \approx 0.2736$ coincides exactly with a **discontinuous entanglement-severing transition** ($2\sqrt{2} \to 0$).

![Part VI Verification (2): (a) Dynamical CHSH generation from a product state and Bell-state preservation on spatially separated sites (Alice=1, Bob=2); (b) Discontinuous CHSH collapse ($2\sqrt{2}\to 0$, purity $1.0\to 0.25$) at the event horizon threshold $\epsilon_c \approx 0.2736$.](results/fig_p6_chsh_and_horizon_transition_en.png)

---

## 5. Open Problems and Reproducibility

- **Open problems (deferred to Part VII+)**: (1) distance-dependent two-particle wave-packet collision cross-sections on $O_8^{\otimes K}\;(K\ge 3)$; (2) three-generation structure via the exceptional Jordan algebra $J_3(\mathbb{O})$, precision fermion mass ratios, and first-principles derivation of $\alpha$ and $\alpha_s$ (planned Part VII).
- **Reproducibility files**: `tri_loop_part6_scattering_chsh_v0.py` (SC-1–SC-5), `results/tri_loop_part6_scattering_chsh_v0.txt`, `tri_loop_part6_chsh_and_coulomb_v0.py` (SC-1b–SC-5b), `results/tri_loop_part6_chsh_and_coulomb_v0.txt`, `plot_part6_figures_v0.py`, and `results/fig_p6_*_{ja,en}.png`.

---

## References

### A. Self-Citations
- **[P1]** N. Kihara, *Spontaneous Emergence of Quantum, Gauge, Gravitational, and Fermionic Structures from the 3D Triangle-Loop Product (Integrated Part I)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24622200`.
- **[P2]** N. Kihara, *Unified Connection Hypothesis for Gravity, Electromagnetism, QCD, and Dirac Fermions in the 3D Triangle-Loop Product (Part II)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24622565`.
- **[P3]** N. Kihara, *Spacetime and Internal-Space Readout Connection Hypothesis in the 3D Triangle-Loop Product (Part III)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24625962`.
- **[P4]** N. Kihara, *Electroweak Symmetry Breaking and Maximal Parity Violation in the 3D Triangle-Loop Product (Part IV)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24626592`.
- **[P5]** N. Kihara, *Multi-Shell Large-Lattice Gravity, Causal Propagation, and Spin-2 Unlock of $g_{00}(d)$ in the 3D Triangle-Loop Product (Part V)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24627519`.

### B. External Primary Literature
- **[E1]** J. F. Clauser, M. A. Horne, A. Shimony, and R. A. Holt, *Proposed Experiment to Test Local Hidden-Variable Theories*, Phys. Rev. Lett. **23**, 880–884 (1969).
- **[E2]** R. Horodecki, P. Horodecki, and M. Horodecki, *Violating Bell inequality by mixed spin-1/2 states: necessary and sufficient condition*, Phys. Lett. A **200**, 340–344 (1995).
- **[E3]** B. S. Cirel'son (Tsirelson), *Quantum generalizations of Bell's inequality*, Lett. Math. Phys. **4**, 93–100 (1980).
