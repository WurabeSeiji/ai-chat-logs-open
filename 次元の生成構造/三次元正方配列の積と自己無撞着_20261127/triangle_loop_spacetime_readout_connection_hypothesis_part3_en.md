# Spacetime and Internal-Space Readout Connection Hypothesis in the 3D Triangle-Loop Product (Part III)
## — Logarithmic Transform $O = e^{\ln A + n\theta}$, Observer-Centered Past Light Cone, Resolution Gauge, and $(x, y, z, t, R, Q_1, Q_2, Q_3)$ Mapping —

**Author:** Noriaki Kihara (WF System Co., Ltd.)
**ORCID:** 0009-0004-6753-4020
**Date:** November 28, 2026
**Version DOI:** 10.5281/zenodo.24625963
**Concept DOI:** 10.5281/zenodo.24625962
**Status:** Readout connection hypothesis paper (**Part III**) following the Unified Full Paper Part I (Concept DOI: `10.5281/zenodo.24622200`) and Unified Connection Hypothesis Part II (Concept DOI: `10.5281/zenodo.24622565`). Motivated directly by two concrete gaps remaining after Part I and Part II—(i) the lack of an explicit readout map from the $\mathrm{Tri}(T)$ state space to $(x, y, z, t, R, Q_1, Q_2, Q_3)$ and (ii) the absence of a resolution-based definition of visibility—this paper directly connects the readout geometry developed in prior thought experiments [K9–K14] to the actual deterministic $\mathrm{Tri}(T)$ states.

---

## Abstract

Part I [P1] established without external axioms or random numbers that the single normalized real 3-tensor map $\mathrm{Tri}(T)_{ijk} = \sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}$ spontaneously generates Hurwitz crystal selection, four quantum theorems (Q1–Q4), spherical 4-class metric distortion with critical collapse, Hurwitz mass hierarchy, 28D $Spin(8)$ Triality spin-$1/2$, and CAR `1:3:3:1` fermion charges. Part II [P2] introduced the $O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$ separation hypothesis to connect gravity, electromagnetism, QCD, and Dirac matter.

Nevertheless, two essential gaps remained before closing the connection to general relativity and the Standard Model: (1) how and why the $\mathrm{Tri}(T)$ state space (phase space) maps onto physical spacetime and internal coordinates $(x, y, z, t, R, Q_1, Q_2, Q_3)$, and (2) how "visible vs. hidden" is defined by an observer's resolution rather than by inspecting all $64^3$ components from an external God's-eye view. In this Part III, we bridge our prior readout frameworks—central projection and split-octonion null balance [K9], two sections of the complex exponent [K10], resolution as a gauge [K11], light-cone projection axes [K12, K13], and cross-term visibility [K14]—directly onto the deterministic $\mathrm{Tri}(T)$ states (`tri_loop_spacetime_log_resolution_bridge_v0.py`, no random numbers) and establish five results:

1. **`(1+3)+(1+3)` Split-Octonion Null Balance & Projection-Axis Covariance**: On the Part I vacuum time crystal $H = O_8 \otimes M_8$, all seven Fano quaternion `(1+3)+(1+3)` splits strictly satisfy $|a|^2 = |b|^2 = 0.500000$ (max difference **`3.331e-16`**), and 1D ($\mathbb{R}$), 2D ($\mathbb{C}$), and 4D ($\mathbb{H}$) projection-axis choices preserve the exact ratio $d/8$ across all selections (max deviation `5.55e-17` to `2.22e-16`).
2. **`1 + 3 + 4` Split of Outer $O_8^{\mathrm{ext}}$ and Pure $SU(2)$**: Inside $G_2$ (14D), the subgroup preserving the `3` ($\{e_1,e_2,e_3\}$) and `4` ($\{e_4,e_5,e_6,e_7\}$) sectors is 6D $SO(4)$, containing a **3D pure $SU(2)$ subgroup that leaves 3-space strictly invariant (max error `8.327e-17`) while rotating the 4-space as spin-$1/2$ ($2\pi$ rotation $= -I$, error `5.778e-15`)**. Thus $O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$ (64 real DOF) accommodates gravity (`1+3`), electroweak $SU(2)_L$/Higgs 4-component doublet (outer `4`), and $SU(3)_C \times U(1)_Q$ (inner `2+6`) without overlap.
3. **Logarithmic Transform $O = e^{\ln A + n\theta}$ (Origin of Duality & Local Curvature Radius $R$)**: Because $\mathrm{Tri}$ operates in the multiplicative state space while spacetime and rapidity live in the additive exponent, $\ln A$ and $n\theta$ sit side by side, yielding Born's rule $\cos^2(\theta/2)$ and the Lorentz factor $\cosh\phi$ as two sections of one identity (error `2.220e-16`). Re-reading the Part I §6 mass sweep via $\Delta\phi = \frac{1}{2}\ln(E_{00}/E_{\mathrm{far}})$ and $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$ shows monotonic shrinking `1.0000 → 0.8328` for $\epsilon = 0 \to 0.270$, followed by a discontinuous horizon cut $R_{\mathrm{loc}}/R_0 \to 0$ at $\epsilon = 0.275$.
4. **Observer-Centered Past Light Cone & Tidal Splitting**: Interpreting $n$ not as shared Cartesian $(x,y,z)$ but as sightlines on the observer's own past light cone, the logarithmic past distance $D_{\log}(A_{\mathrm{obs}}\to B)$ is strictly isotropic in vacuum ($0.000000 \pm 10^{-16}$, equivalence principle) for any observer, whereas in the presence of mass P4 a distant observer ($A_{\mathrm{obs}}=1$) measures a clear **tidal splitting $\Delta D_{\mathrm{tidal}} = +0.028615$** between the radial mass direction (`-0.006094`) and the six tangential directions (`-0.034709 ± 5.27e-16`).
5. **Resolution Gauge Definition of Visibility**: (a) Defining visibility at $\tau = 75$ by $L_3 > \delta_{\mathrm{res}}$ and sweeping $\delta_{\mathrm{res}}$ across 18 orders of magnitude (`1e-2` to `1e-20`) leaves the visible set (64 triangles) and weighted view angle $\langle\rho\rangle_w = 0.7767719701$ invariant to 10 digits (resolution = gauge). (b) Under the $\mathbb{Z}_2^3$ Walsh–Hadamard low-pass window $K_{\mathrm{cut}} \in \{0,1,2,3\}$, the minimum resolution required to fully resolve Hurwitz particles $\mathbb{O}_8, \mathbb{H}_4, \mathbb{C}_2, \mathbb{R}_1$ (codimension $m=0,1,2,3$) strictly matches $K_{\mathrm{cut}} = m$.

Section 6 cleanly separates the remaining open problems (multi-shell $O_8^{\otimes K}$ continuum limit, $J_3(\mathbb{O})$ and three generations, dynamical electroweak breaking and parity violation, $\alpha \approx 1/137$, two-particle dynamical scattering, and dynamical CHSH violation).

---

## 1. Introduction: Three Remaining Dissatisfactions as the Motivation Chain

In a paper centered on physical readout and interpretation, stating conclusions first invites the impression of post-hoc fitting. The motivation for Part III is the direct chain of three concrete dissatisfactions remaining after Part I [P1] and Part II [P2]:

1. **Dissatisfaction 1 (Unclosed Spacetime Connection)**: Although Part I §6 derived 4-class metric distortion and critical collapse and Part II connected gravity, EM, QCD, and Dirac matter, our prior work established central projection $\sum z_n^2 = 0 \mapsto \sum z_n^2 = R^2$ [K9, K12] and the 8D sign-vector null equation
   $$x^2 + y^2 + z^2 + (it)^2 = R^2 + Q_1^2 + Q_2^2 + Q_3^2. \tag{1}$$
   Part I and Part II lacked the explicit bridge explaining how and why the $\mathrm{Tri}(T)$ state space maps onto $(x, y, z, t, R, Q_1, Q_2, Q_3)$.
2. **Dissatisfaction 2 (Missing Resolution Definition of Visibility)**: Part I §2.3 measured the view angle $\rho = L_3/L_1$ by inspecting all $64^3$ tensor entries externally. What was missing is the principle from [K11, K14]: observation is not "existence vs. non-existence" but "visible vs. hidden relative to an observer's resolution," where resolution is a gauge freedom and the view angle is its invariant.
3. **Dissatisfaction 3 (Misreading $n$ as Shared $xyz$ and the Role of `1–7` Octonion Sectors)**: Reading the unit imaginary vector $n$ in $O = e^{\ln A + n\theta}$ as shared Cartesian $(x,y,z)$ violates background independence; the apex of the light cone is always the observer, so $n$ labels sightlines on the observer's **past light cone (sector coordinates)**. Moreover, when reading `1–7` as `(1+3)+(1+3)` or `1+3+4` inside $O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$, the physical role of the extra `4` components must be identified.

As in Part II, every claim is explicitly tagged as **【Known】**, **【Part I/II Established】**, **【Part III Verified】**, or **【Hypothesis / Unverified】**.

---

## 2. `(1+3)+(1+3)` and `1+3+4` Readings of $O_8$ and the 64D Role Table

### 2.1 Split-Octonion Null Balance and 1D/2D/4D Projection-Axis Covariance 【Part III Verified P3-1】

Testing the actual Part I vacuum time crystal $H = O_8 \otimes M_8$ across all seven Fano quaternion subalgebras $H_4 \subset O_8$ and across 1D ($\mathbb{R}$), 2D ($\mathbb{C}$), and 4D ($\mathbb{H}$) projection axes [K9, K13] yields:

- **All 7 `(1+3)+(1+3)` splits**: $|a|^2 = |b|^2 = 0.500000$ with max difference **`3.331e-16`** (exact split-octonion null balance).
- **1D/2D/4D projection-axis covariance**: 1D ($\mathbb{R}$, $1+7$) = `0.125000` (dev `5.551e-17`), 2D ($\mathbb{C}$, $2+6$) = `0.250000` (dev `2.220e-16`), 4D ($\mathbb{H}$, $4+4$) = `0.500000` (dev `2.220e-16`).
- **With localized mass P4 ($\epsilon=0.16, \tau=40$)**: $|a|^2 = 0.499002$ vs. $|b|^2 = 0.499488$ (ratio `0.999026`, $\Delta\phi_{ab} = -0.000487$), reflecting the greater multiplicity of tangential edges (`+7.2%`) relative to the central site (`+32.5%`).

### 2.2 `1 + 3 + 4` Split of Outer $O_8^{\mathrm{ext}}$: Pure $SU(2)$ Subgroup 【Part III Verified P3-5】

A single $O_8$ (8 real DOF) lacks sufficient degrees of freedom for all Standard Model quantum numbers; Dixon and Furey's $\mathbb{C}\otimes\mathbb{H}\otimes\mathbb{O}$ likewise has $2\times 4\times 8 = 64$ real DOF [E1, E2]. Extracting via SVD the subgroup of $G_2$ (14D) on $O_8^{\mathrm{ext}}$ that preserves the `3` ($\{e_1,e_2,e_3\}$) and `4` ($\{e_4,e_5,e_6,e_7\}$) sectors yields:

- **Block-preserving subgroup**: strictly **6D $SO(4) \cong SU(2)_{\mathrm{space}} \times SU(2)_{\mathrm{pure4}}$**.
- **Pure 4-space subgroup**: **3D $SU(2)_{\mathrm{pure4}}$**, whose action on 3-space $\{e_1,e_2,e_3\}$ vanishes to **`8.327e-17`**.
- **Spin-$1/2$ doublet transformation**: under a $2\pi$ rotation by $SU(2)_{\mathrm{pure4}}$, the 4-space $\{e_4,e_5,e_6,e_7\}$ flips sign to $-I$ with error **`5.778e-15`**, matching the geometric carrier of **weak isospin $SU(2)_L$ / the 4-real-component Higgs doublet** 【Hypothesis / Unverified】.

Table 1 (and Fig. 2(b)) summarizes the complete, non-overlapping 64D assignment of $O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$.

| Sector | Split | Exact Symmetry | Physical DOF / Coordinates | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Outer $O_8^{\mathrm{ext}}$** | **`1`** ($e_0$) | Scalar ($\ln A$) | Time reference $t$, scale, local curvature radius $R$ | 【Part III Verified】+【Hypothesis】 |
| | **`3`** ($e_1,e_2,e_3$) | $SO(3)$ (3D) | Past light-cone sightlines $n$ ($x,y,z$ readout), gravity tides | 【Part I/III Verified】 |
| | **`4`** ($e_4,e_5,e_6,e_7$) | **Pure $SU(2)$ (3D)** | **Weak isospin $SU(2)_L$, Higgs 4-comp. doublet, spin-$1/2$** | 【Part III Verified】+【Hypothesis】 |
| **Inner $O_8^{\mathrm{int}}$** | **`2 + 6`** ($e_7$ fixed) | **$SU(3)_C$ (8D) $\times U(1)_Q$** | **Strong/EM forces, CAR `1:3:3:1` charges $(Q_1,Q_2,Q_3)$, clock $M_8$** | 【Part I/II Established】 |

---

## 3. Logarithmic Transform $O = e^{\ln A + n\theta}$: Origin of Duality and Local Curvature Radius $R$

### 3.1 Multiplicative Dynamics vs. Additive Spacetime 【Known】+【Part III Verified P3-2】

Any octonion $O \in \mathbb{O}$ admits the single-exponent form
$$O = A\,e^{n\theta} = e^{\ln A + n\theta} \tag{2}$$
(round-trip error `1.665e-16`). While $\mathrm{Tri}(T)$ evolves in the **multiplicative** state space, physical spacetime intervals, rapidities, and energies live in the **additive** exponent ($\ln A + n\theta$). Because the logarithm places scale $\phi = \ln A$ and phase $n\theta$ in the same linear space, the identity $|\cosh(w/2)|^2 \pm |\sinh(w/2)|^2 = (\cosh\phi, \cos\theta)$ (error **`2.220e-16`**) [K10] renders Born's rule ($\cos^2(\theta/2)$ at $\phi=0$) and the Lorentz factor ($\cosh\phi$ at $\theta=0$) dual sections related by a $90^\circ$ rotation. Normalization by Frobenius norm corresponds merely to a uniform origin shift in $\ln A$, and a compact bounded state space [K14] maps to an open horizon ($A \to 0 \Rightarrow \ln A \to -\infty$).

### 3.2 Local Curvature Radius $R_{\mathrm{loc}}$ Shrinking and Horizon Cut 【Part III Verified P3-2】

Re-reading the Part I §6 mass P4 sweep ($\tau=40$) via $\Delta\phi = \frac{1}{2}\ln(E_{00}/E_{\mathrm{far}})$ and $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$ [K12] yields Table 2 and Fig. 1(a).

| Mass amplitude $\epsilon$ | $E_{00}/E_{\mathrm{far}}$ | $\Delta\phi = \frac{1}{2}\ln(E_{00}/E_{\mathrm{far}})$ | **$R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$** | $\cosh(\Delta\phi)$ | $E_{\mathrm{far}}$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| `0.000` | `1.000000` | `0.000000` | **`1.000000`** | `1.000000` | `0.015625` |
| `0.050` | `1.047528` | `+0.023217` | **`0.977051`** | `1.000270` | `0.015495` |
| `0.100` | `1.174017` | `+0.080216` | **`0.922917`** | `1.003219` | `0.015162` |
| `0.160` | `1.324838` | `+0.140645` | **`0.868798`** | `1.009907` | `0.014789` |
| `0.220` | `1.400375` | `+0.168370` | **`0.845041`** | `1.014208` | `0.014686` |
| `0.270` | `1.441867` | `+0.182969` | **`0.832794`** | `1.016786` | `0.014492` |
| **`0.275`** | **$\infty$ (horizon cut)** | **$+\infty$** | **`0.000000`** | **$\infty$** | **`0.000e+00`** |

This directly connects the working hypothesis of [K12] ("gravity as shrinking of the local curvature radius $R$ and event horizon") with the Part I critical collapse of $\mathrm{Tri}(T)$.

---

## 4. Observer-Centered Past Light Cone and Tidal Splitting 【Part III Verified P3-3】

Treating $n$ as sightlines from an observer site $A_{\mathrm{obs}}$ on its own past light cone, we define the logarithmic past distance to site $B \neq A_{\mathrm{obs}}$ by
$$D_{\log}(A_{\mathrm{obs}} \to B) \equiv -\frac{1}{2}\ln\!\left(\frac{E_{A_{\mathrm{obs}}, B}}{E_{A_{\mathrm{obs}}, A_{\mathrm{obs}}}}\right). \tag{3}$$

- **Vacuum $H$ (equivalence principle)**: For both $A_{\mathrm{obs}}=0$ and $A_{\mathrm{obs}}=1$, all seven directions are strictly isotropic: $D_{\log} = 0.000000 \pm 8.32\times 10^{-16}$ (and $\pm 2.29\times 10^{-11}$ at $A=1$). Uniform curvature is invisible; the local readout is flat [K9].
- **With mass P4 ($A=0, \epsilon=0.16$)**:
  - Observer on the mass ($A_{\mathrm{obs}}=0$): all 7 directions remain strictly spherical ($D_{\log} = +0.134551 \pm 7.77\times 10^{-16}$).
  - **Distant observer ($A_{\mathrm{obs}}=1$)**: the radial mass direction ($B=0$) gives `D_log = -0.006094`, whereas the six tangential directions ($B=2\dots 7$) give `D_log = -0.034709 ± 5.27e-16`, producing a clear **tidal splitting $\Delta D_{\mathrm{tidal}} = +0.028615$** (Fig. 1(b)).

This directly verifies the prediction of [K9] that uniform gravity disappears in an observer-centered readout while directional tidal distortion survives.

![Part III Verification (1): (a) Local curvature radius shrinking $R_{\mathrm{loc}}/R_0 = e^{-\Delta\phi}$ and critical horizon cut via the logarithmic map; (b) Observer-centered past light cone showing vacuum isotropy (equivalence principle) vs. distant-observer tidal splitting $\Delta D_{\mathrm{tidal}} = +0.0286$.](results/fig_p3_log_curvature_and_tidal_en.png)

---

## 5. Resolution-Gauge Definition of Visibility 【Part III Verified P3-4】

To replace the external God's-eye view in Part I, we implement two resolution filters [K11, K14] directly on $\mathrm{Tri}(T)$ states (Fig. 2(a)):

### 5.1 18-Order Threshold Sweep: Resolution as Gauge, View Angle as Invariant

At $\tau = 75$ (representative particle component $i=4, j=0, k=4$), defining visibility by **$L_3 > \delta_{\mathrm{res}}$** and sweeping $\delta_{\mathrm{res}} \in \{10^{-2}, 10^{-5}, 10^{-10}, 10^{-20}\}$ leaves the visible set strictly at **64 triangles** and the weighted view angle at **$\langle\rho\rangle_w = 0.7767719701$** (10-digit invariant). Detached background triangles fall below resolution (`null` [K14]), confirming that resolution is a gauge and the view angle is a gauge invariant [K11].

### 5.2 $\mathbb{Z}_2^3$ Harmonic Low-Pass Window $K_{\mathrm{cut}}$ and Hurwitz Visibility Boundary

Applying the $\mathbb{Z}_2^3$ Walsh–Hadamard low-pass window $K_{\mathrm{cut}} \in \{0,1,2,3\}$ to the exact Part I §7 Hurwitz 4-level probes ($\epsilon=0.05, \tau=30$) yields Table 3.

| Hurwitz Particle (codimension $m$) | $K_{\mathrm{cut}}=0$ (DC) | $K_{\mathrm{cut}}\le 1$ | $K_{\mathrm{cut}}\le 2$ | $K_{\mathrm{cut}}\le 3$ (Full) | Full-visibility $k_{\max}$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **$\mathbb{O}_8$ ($m=0$: background)** | `1.000000` | `1.000000` | `1.000000` | `1.000000` | **`0`** |
| **$\mathbb{H}_4$ ($m=1$: P5)** | `0.900000` | `1.000000` | `1.000000` | `1.000000` | **`1`** |
| **$\mathbb{C}_2$ ($m=2$)** | `0.749250` | `0.916417` | `1.000000` | `1.000000` | **`2`** |
| **$\mathbb{R}_1$ ($m=3$: P4)** | `0.556705` | `0.746689` | `0.936672` | `1.000000` | **`3`** |

The codimension $m$ strictly equals the minimum resolution window $K_{\mathrm{cut}} = m$ needed to fully resolve each particle; for $K_{\mathrm{cut}} < m$, finer modes ($k > K_{\mathrm{cut}}$) remain hidden as internal rigid-crystal structure [K11].

![Part III Verification (2): (a) Hurwitz 4-level particle visibility boundary at $K_{\mathrm{cut}} = m$ and 18-order threshold invariance of the view angle; (b) Complete 64D readout mapping from $O_8^{\mathrm{ext}} \otimes O_8^{\mathrm{int}}$ to $(x, y, z, t, R, Q_1, Q_2, Q_3)$ and $SU(2)_L$.](results/fig_p3_resolution_gauge_and_mapping_en.png)

---

## 6. Clean Separation of Remaining Open Problems

The following essential problems remain **【Hypothesis / Unverified】** and are explicitly deferred to subsequent work:

1. **Multi-shell $O_8^{\otimes K}$ continuum limit**: Deriving $1/r$ Newtonian gravity, $1/r^2$ Coulomb force, light-cone propagation speed $c$, Einstein's field equations, and cumulative gravitational time dilation $g_{00}(r)$ by interpreting nested shells $O_8^{\otimes K}$ as concentric depths of the past light cone.
2. **Dynamical electroweak breaking and maximal parity violation**: Dynamically coupling the outer pure $SU(2)$ (4-space) with the inner $M_8$ chiral asymmetry (`-1.41%`), left-handed weak selection, and $W/Z$ mass generation.
3. **Three generations, precision mass ratios, and coupling constants**: Connecting the exceptional Jordan algebra $J_3(\mathbb{O})$ [E3] or $O_8^{\otimes K}$ nesting to three fermion generations, exact lepton/quark mass ratios, and first-principles derivation of $\alpha \approx 1/137.036$ and $\alpha_s$.
4. **Dynamical scattering and multi-body quantum correlations**: Two-particle attractive/repulsive scattering, electron–photon scattering, and dynamical CHSH Bell inequality violation.

---

## 7. Reproducibility and Published Files

- Verification script: `tri_loop_spacetime_log_resolution_bridge_v0.py` (deterministic, ~23 s)
- Numerical log: `results/tri_loop_spacetime_log_resolution_bridge_v0.txt`
- Figure script: `plot_part3_figures_v0.py`
- Figures: `results/fig_p3_log_curvature_and_tidal_{ja,en}.png`, `results/fig_p3_resolution_gauge_and_mapping_{ja,en}.png`

---

## References

### A. Self-Citations
- **[P1]** N. Kihara, *Spontaneous Emergence of Quantum, Gauge, Gravitational, and Fermionic Structures from the 3D Triangle-Loop Product (Unified Full Paper, Part I)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24622200`.
- **[P2]** N. Kihara, *Unified Connection Hypothesis for Gravity, Electromagnetism, QCD, and Dirac Matter in the 3D Triangle-Loop Product (Part II)*, Zenodo (2026-11-28), Concept DOI: `10.5281/zenodo.24622565`.
- **[K9]** N. Kihara, *Relational Readout of Quaternions and Inflation as Crystallization*, Zenodo (2026-10-21), Concept DOI: `10.5281/zenodo.23523237`.
- **[K10]** N. Kihara, *Did We Properly Understand Complex Numbers?*, Zenodo (2026-10-23), Concept DOI: `10.5281/zenodo.23585436`.
- **[K11]** N. Kihara, *Observer Resolution Is a Gauge*, Zenodo (2026-11-10), Concept DOI: `10.5281/zenodo.24097611`.
- **[K12]** N. Kihara, *Emergence of Spacetime and Elementary Particles from a Bounded Hilbert Space and Octonions*, Zenodo (2026-11-16), Concept DOI: `10.5281/zenodo.24267347`.
- **[K13]** N. Kihara, *Projection-Axis Selection and Covariance of Conservation Laws*, Zenodo (2026-11-16), Concept DOI: `10.5281/zenodo.24273723`.
- **[K14]** N. Kihara, *The Coin-and-Cup Thought Experiment*, Zenodo (2026-11-26), Concept DOI: `10.5281/zenodo.24562867`.

### B. External Primary Literature
- **[E1]** G. M. Dixon, *Division Algebras: Octonions, Quaternions, Complex Numbers and the Algebraic Design of Physics*, Kluwer Academic Publishers (1994).
- **[E2]** C. Furey, *$SU(3)_C \times SU(2)_L \times U(1)_Y (\times U(1)_X)$ as a symmetry of division algebraic ladder operators*, Eur. Phys. J. C **78**, 375 (2018), arXiv:1806.00612.
- **[E3]** J. C. Baez, *The Octonions*, Bull. Amer. Math. Soc. **39**(2), 145–205 (2002), arXiv:math/0105155.
- **[E4]** S. Coleman and J. Mandula, *All Possible Symmetries of the $S$ Matrix*, Phys. Rev. **159**(5), 1251–1256 (1967).
