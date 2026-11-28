# Emergent Quantum Structure from the 3D Triangle-Loop Product: Exact Discretization, Spontaneous Single Complex Structure, Born Rule, and Discrete Uncertainty Saturation (Fourth Rapid Communication)
## ― Analytical Derivation of Complex Hilbert Structure, Probability-Squared Law, and Minimum Resolution from a Real 3-Tensor Without Quantum Postulates ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)
**ORCID:** 0009-0004-6753-4020
**Date:** November 28, 2026
**Status:** Fourth Rapid Communication following Rapid Reports 1 (DOI: 10.5281/zenodo.24616931), 2 (DOI: 10.5281/zenodo.24618400), and 3 (DOI: 10.5281/zenodo.24619357). Establishes four analytical theorems with machine-precision numerical verification showing that the core pillars of quantum mechanics emerge directly from a purely real 3-tensor iteration without quantum postulates.

---

## Abstract

We investigate the background-free, purely real 3-tensor triangle-loop product $\mathrm{Tri}(T)_{ijk} = \sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}$ without assuming spacetime, an action, complex amplitudes, or probabilistic postulates. Synthesizing the $N = 64$ ($O_8 \otimes O_8$) results established in Rapid Communications 1–3, we prove and numerically verify four foundational theorems showing that quantum structure emerges intrinsically from the algebra itself:

1. **Theorem Q1 (Exact Discretization)**: The cubic phase map $\theta \mapsto 3\theta$ combined with the null/chiral balance condition $\sum z_n^2 = 0$ quantizes periodic phases to rational angles $2\pi\mathbb{Z}/(3^p - 1)$, restricts spatial crystals to Hurwitz dimensions ($1, 2, 4, 8$), and discretizes internal view-angle couplings $\rho = L_3/L_1$ to exactly 7 values.
2. **Theorem Q2 (Spontaneous Single Complex Structure $\mathbb{C}$)**: Although the system starts from a real array $T \in \mathbb{R}^{64\times 64\times 64}$, the inner sector $M_8$ of Layer 1 ($H = O_8 \otimes M_8$) suppresses symmetric anisotropic Pauli components ($\sigma_x, \sigma_z$) to $\sim 10^{-16}$, spontaneously selecting a pure complex field $\mathbb{C} \cong \{c_0 I + c_J(i\sigma_y)\}$. Furthermore, Layer 2 ($\mathbb{C}_{\mathrm{out}} \otimes_{\mathbb{R}} \mathbb{C}_{\mathrm{in}} \cong \mathbb{C}_+ \oplus \mathbb{C}_-$) is unstable and collapses into Layer 3, $U_{\text{Layer 3}} = \frac{1}{\sqrt{2}}(u\otimes v + (Ju)\otimes(-Jv))$, locking a single global imaginary unit $J$ while exactly doubling $\|P\|$ ($|\sigma_{1,2} - 1/\sqrt{2}| = 3.23 \times 10^{-7}$).
3. **Theorem Q3 ($O(N)$ Covariance, $L^2$ Uniqueness, and the Born Rule)**: The map $\mathrm{Tri}$ is strictly $O(64)$-covariant (relative error $8.84 \times 10^{-16}$). Among $L^p$ norms, only $p = 2$ is invariant under $O(64)$ rotations (relative change $2.22 \times 10^{-16}$, vs. $+351.8\%$ for $L^1$ and $-42.0\%$ for $L^4$). Together with the emergent complex structure $\mathbb{C}$ (Theorem Q2) and dimension $d \ge 3$, Gleason's theorem uniquely selects the Born rule $P(\Pi) = \mathrm{Tr}(\rho\Pi) = |\langle\phi|\psi\rangle|^2$.
4. **Theorem Q4 ($\mathbb{Z}_2^3$ Discrete Uncertainty Saturation and Minimum Resolution)**: On the outer crystal $O_8 \cong \mathbb{Z}_2^3$, any normalized state $f$ and its Walsh–Hadamard spatial harmonic transform $\hat{f} = H_8 f$ satisfy the exact Donoho–Stark bound $|\mathrm{supp}(f)|\cdot|\mathrm{supp}(\hat{f})| \ge 8$ and Shannon entropic uncertainty bound $H_{\mathrm{space}}(f) + H_{\mathrm{harmonic}}(\hat{f}) \ge \ln 8 \approx 2.079442$. All three bare localization prototypes—nonlocal ($8\times 1$), quaternion core $H_4$ ($4\times 2$), and 1-site core $\delta_{000}$ ($1\times 8$)—strictly saturate both bounds (excess $0.00 \times 10^0$), whereas $\tau = 30$ dressed states ($L_x$, P5, P4) shift strictly above the bound ($+0.051$ to $+0.313$) due to their background clouds. Minimum spatial resolution $\Delta x_{\min} = 1$ site (Planck length $\ell_P$ prototype), UV cutoff $k_{\max} = 3$, minimum step $\Delta\tau_{\min} = 1$ (Planck time $t_P$ prototype), and minimum phase-space cell $\ln 8$ (dimensionless $\hbar$ prototype) arise directly from the algebra.

---

## 1. Introduction and Scope

In Rapid Communications 1–3 [1–3], we showed that iterating the normalized triangle-loop product $U^{(\tau+1)} = \mathrm{Tri}(U^{(\tau)})/\|\mathrm{Tri}(U^{(\tau)})\|_F$ on $N = 64$ ($O_8 \otimes O_8$) selects Hurwitz crystals, undergoes a universal four-stage doubling-entanglement phase transition (rank $1 \to 2 \to 4$), exhibits $SO(3)$ spin selection ($1$ vs. $(1/3)^\tau$) and $G_2$ ($14$) $\supset SU(3)$ ($8$) gauge symmetries, separates nonlocal gauge waves from localized condensates via $\mathbb{Z}_2^3$ spatial harmonics, and supports two-particle coexistence.

A central foundational question remains: is this merely a classical nonlinear tensor map, or does it inherently generate the mathematical axioms of quantum mechanics? Here we demonstrate—via four analytical theorems and deterministic machine-precision verification (`tri_loop_n64_quantum_foundations_verify_v0.py`, zero randomness)—that discretization, a single complex Hilbert structure, the Born rule, and discrete uncertainty saturation with a UV cutoff emerge without external quantum postulates. Throughout, $\tau$ denotes the discrete update step, distinct from emergent physical time $t$.

---

## 2. Theorem Q1: Exact Discretization via $\sum z_n^2 = 0$ and Cubic Self-Consistency

**Theorem Q1 (Phase, Coupling, and Crystal Discretization)**:
Because $\mathrm{Tri}$ is homogeneous of degree 3, phase evolution acts as the tripling map $\theta_n \mapsto 3\theta_n \pmod{2\pi}$. Any self-consistent fixed point or period-$p$ cycle is therefore restricted to discrete rational angles:
$$(3^p - 1)\theta_n \equiv 0 \pmod{2\pi} \quad \Longrightarrow \quad \theta_n \in \frac{2\pi}{3^p - 1}\mathbb{Z}.$$
Under the null/chiral condition $\sum_n z_n^2 = \sum_n r_n^2 e^{2i\theta_n} = 0$ (the $I + \sigma_y$ sector in real representation [2]), arbitrary continuous phase distributions cannot balance the quadratic invariant, forcing discrete polygonal equilibria. Spatially, Hurwitz's theorem restricts self-consistent division-algebra crystals to dimensions $1, 2, 4, 8$, and the internal triangle view-angle couplings $\rho = L_3/L_1$ of Layer 1 ($H = O_8 \otimes M_8$) quantize into exactly 7 discrete values ($1.0000, 0.7646, 0.7463, 0.7268, 0.5558, 0.5424, 0.4147$) [2].

---

## 3. Theorem Q2: Spontaneous Emergence of a Single Complex Structure $\mathbb{C}$

A defining postulate of quantum mechanics is that state spaces are Hilbert spaces over $\mathbb{C}$ (not $\mathbb{R}$) and composite systems combine via $\otimes_{\mathbb{C}}$. Starting from a purely real tensor $T \in \mathbb{R}^{64\times 64\times 64}$, this structure emerges in two stages:

**Theorem Q2 (Single Complex Structure Lock)**:
1. **Local Complexification (Layer 1)**: Decomposing the $2\times 2$ real inner blocks ($b \in \{0,1\},\; a,c \in \{4,5\}$) of $M_8$ in the Pauli basis $B = c_0 I + c_J J + c_x \sigma_x + c_z \sigma_z$ ($J = -i\sigma_y,\; J^2 = -I$) yields $|c_x|, |c_z| \sim 10^{-16}$ (machine zero), leaving a pure complex field $\mathbb{C} \cong \{c_0 I + c_J J\}$ (Table 1).
2. **Global Single-Clock Lock (Layer 2 $\to$ Layer 3)**: A real tensor product of two independent complex structures (Layer 2: $\mathbb{C}_{\mathrm{out}} \otimes_{\mathbb{R}} \mathbb{C}_{\mathrm{in}} \cong \mathbb{C}_+ \oplus \mathbb{C}_-$) is unstable and collapses into Layer 3,
   $$U_{\text{Layer 3}} = \frac{1}{\sqrt{2}}\Bigl(u \otimes v + (Ju)\otimes(-Jv)\Bigr),$$
   locking $J_{\mathrm{out}} = -J_{\mathrm{in}}$ into a single global imaginary unit $J$ (a single clock) while strictly doubling $\|P\|$ (Case 1 in [2]: $0.240854 \to 0.481707$; verified here with $|\sigma_{1,2} - 1/\sqrt{2}| = 3.231 \times 10^{-7}$).

| Verification Item (Theorem Q2) | Measured Value | Interpretation |
| :--- | :--- | :--- |
| Layer 1 $M_8$ ($b=0$) complex coeffs $c_0, c_J$ | $+0.27343360,\; +0.49200681$ | Pure complex field $\mathbb{C}$ |
| Layer 1 $M_8$ ($b=0$) anisotropic $|c_x|, |c_z|$ | $1.943 \times 10^{-16},\; 6.661 \times 10^{-16}$ | Exact zero (machine precision) |
| Layer 1 $M_8$ ($b=1$) anisotropic $|c_x|, |c_z|$ | $5.551 \times 10^{-17},\; 5.274 \times 10^{-16}$ | Exact zero (machine precision) |
| Layer 3 single-complex lock $|\sigma_{1,2} - 1/\sqrt{2}|$ | $3.231 \times 10^{-7}$ ($\|P\| \times 2$) | $\mathbb{C}_{\mathrm{out}} \otimes_{\mathbb{R}} \mathbb{C}_{\mathrm{in}} \to \mathbb{C}_{\mathrm{single}}$ |

---

## 4. Theorem Q3: $O(N)$ Covariance, $L^2$ Uniqueness, and the Born Rule

**Theorem Q3 ($L^2$ Uniqueness and Born Rule $P = |\psi|^2$)**:
For any orthogonal transformation $R \in O(N)$ ($R^T R = I$), the index contraction of $\mathrm{Tri}$ satisfies strict covariance:
$$\mathrm{Tri}\bigl((R \otimes R \otimes R) U\bigr) = (R \otimes R \otimes R)\,\mathrm{Tri}(U).$$
Applying a deterministic Givens-product rotation $R_{64} \in SO(64)$ ($\|R_{64}^T R_{64} - I\| = 1.890 \times 10^{-15}$) to $H$ confirms a relative covariance error of **$8.843 \times 10^{-16}$** (Table 2, Fig. 1(b)).

Under $O(N)$ covariance, the **unique** basis-invariant $L^p$ norm is **$p = 2$** (Hilbert/Euclidean norm). Under $R_{64}$, the $L^1$ norm changes by $+351.8\%$ and $L^4$ by $-42.0\%$, whereas $L^2$ is invariant to **$2.220 \times 10^{-16}$**. Thus any basis-independent partitioning across orthogonal subspaces (channels, sites, harmonics, spins, or Bell pairs) is fixed to the squared-amplitude measure $|c_m|^2$. Combined with the emergent complex structure $\mathbb{C}$ (Theorem Q2) and dimension $d \ge 3$, **Gleason's theorem** [4] uniquely forces every non-contextual probability measure on orthogonal projections $\Pi$ to take the Born-rule form:
$$P(\Pi) = \mathrm{Tr}(\rho\,\Pi) = |\langle \phi | \psi \rangle|^2.$$
This directly explains why the spin-$s$ amplitude law $1/\sqrt{s+1}$ (squared weights $1:1/2:1/3$) [3], Bell-pair singular values $(1/2)^2 \times 4 = 1$ [2], and two-particle coexistence energies [3] all obey exact $L^2$ conservation.

| $O(64)$ Rotation Test (Theorem Q3) | Before Rotation | After Rotation | Relative Change |
| :--- | :---: | :---: | :---: |
| $\mathrm{Tri}$ covariance error $\|\mathrm{Tri}(R^3 H)-R^3\mathrm{Tri}(H)\|/\|\mathrm{Tri}(H)\|$ | — | — | **$8.843 \times 10^{-16}$** |
| $L^1$ norm ($p = 1$) | $21.531074$ | $97.269545$ | $+351.8\%$ (non-invariant) |
| **$L^2$ norm ($p = 2$: Hilbert / Born)** | **$1.000000$** | **$1.000000$** | **$2.220 \times 10^{-16}$ (exact invariant)** |
| $L^4$ norm ($p = 4$) | $0.227654$ | $0.132001$ | $-42.0\%$ (non-invariant) |

![Verification of Theorems Q2 and Q3. (a) Pauli decomposition of inner $M_8$ in Layer 1 and single-complex-structure locking in Layer 3. (b) Exact $O(64)$ covariance ($8.84 \times 10^{-16}$) and rotational invariance test of $L^p$ norms ($p=1,2,4$), uniquely singling out $p=2$.](results/fig_r4_complex_and_born_en.png)

---

## 5. Theorem Q4: $\mathbb{Z}_2^3$ Discrete Uncertainty Saturation and Minimum Resolution

**Theorem Q4 (Donoho–Stark and Entropic Uncertainty Bounds)**:
For any normalized distribution $f: \mathbb{Z}_2^3 \to \mathbb{R}$ ($\sum_x |f(x)|^2 = 1$) on the outer crystal $O_8 \cong \mathbb{Z}_2^3$ and its normalized Walsh–Hadamard spatial harmonic transform $\hat{f} = H_8 f$, two exact uncertainty inequalities hold [5]:
$$|\mathrm{supp}(f)| \cdot |\mathrm{supp}(\hat{f})| \;\ge\; 8, \qquad H_{\mathrm{space}}(f) + H_{\mathrm{harmonic}}(\hat{f}) \;\ge\; \ln 8 = 3\ln 2 \approx 2.079442.$$

As shown in Table 3 and Fig. 2, all three bare localization prototypes identified in Rapid Communication 3 [3] **strictly saturate both uncertainty bounds** (zero excess). By contrast, when evolved on the time-crystal background $H$ to $\tau = 30$, the dressed stable excitations ($L_x$, P5, P4) strictly respect the bound while shifting upward ($+0.0512$ to $+0.3130$) in direct proportion to their interaction clouds.

| State (on $\mathbb{Z}_2^3$) | $\|\mathrm{supp}(x)\| \times \|\mathrm{supp}(k)\|$ | $H_{\mathrm{space}}$ | $H_{\mathrm{harmonic}}$ | $H_x + H_k$ | Excess over $\ln 8$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **[Bare] Nonlocal wave ($H$, P3a core)** | **$8 \times 1 = 8$** | $2.079442$ | $0.000000$ | **$2.079442$** | **$0.000000$ (exact saturation)** |
| **[Bare] Quaternion core ($H_4$)** | **$4 \times 2 = 8$** | $1.386294$ | $0.693147$ | **$2.079442$** | **$0.000000$ (exact saturation)** |
| **[Bare] 1-site core ($\delta_{000}$)** | **$1 \times 8 = 8$** | $0.000000$ | $2.079442$ | **$2.079442$** | **$0.000000$ (exact saturation)** |
| [Dressed $\tau=30$] Background $H$ / P3a | $8 \times 1 = 8$ | $2.079442$ | $0.000000$ | $2.079442$ | $+0.000000$ (harm. $0.00\%$) |
| [Dressed $\tau=30$] Spin-1 gauge $L_x$ | — | $2.058650$ | $0.072019$ | $2.130669$ | $+0.051228$ (harm. $3.66\%$) |
| [Dressed $\tau=30$] $H_4$ + cloud P5 | — | $2.022809$ | $0.129825$ | $2.152634$ | $+0.073192$ (harm. $10.00\%$) |
| [Dressed $\tau=30$] 1-site + cloud P4 | — | $1.793664$ | $0.598750$ | $2.392414$ | $+0.312972$ (harm. $45.85\%$) |

**Dimensionless Planck-Scale Prototypes**:
The finite crystal enforces a minimum spatial resolution $\Delta x_{\min} = 1$ site (Planck length $\ell_P$ prototype), a maximum harmonic order $k_{\max} = 3$ (8-mode UV cutoff), a minimum update step $\Delta\tau_{\min} = 1$ (Planck time $t_P$ prototype), and a minimum phase-space cell $\ln 8$ (dimensionless $\hbar$ prototype), eliminating continuum UV divergences at the axiomatic level.

![Verification of Theorem Q4. (a) $\mathbb{Z}_2^3$ entropic uncertainty plane ($H_x + H_k \ge \ln 8$): bare cores strictly saturate the boundary, while $\tau=30$ dressed states shift above the bound by their cloud weights. (b) Donoho–Stark support-product saturation $xy = 8$ with minimum resolution $\Delta x_{\min}=1$ and UV cutoff $k_{\max}=3$.](results/fig_r4_uncertainty_and_cutoff_en.png)

---

## 6. Discussion, Symmetry Inventory, and Outlook

Theorems Q1–Q4 establish that the real 3-tensor triangle-loop product intrinsically generates discretization, a single complex Hilbert structure, the Born rule, and discrete uncertainty saturation with a UV cutoff. To maintain strict scientific precision, we separate established results from open problems:

- **Established (Rapid Reports 1–4)**: Hurwitz crystal selection ($O_8$), 4-stage doubling-entanglement cascade and Bell pairs, $SO(3)$ spin selection ($1$ vs. $(1/3)^\tau$) and $1/\sqrt{s+1}$ law, $G_2$ ($14$) $\supset SU(3)$ ($8$ gluons) and $U(1)$ phase symmetries, two-particle coexistence, and the four quantum structure theorems Q1–Q4.
- **Open Problems for Subsequent Reports**: (1) Quantitative relation between localized energy, local time dilation (first seen as P4 slowdown `d1: 0.904 -> 0.120` in [2]), signal-propagation speed limit $c$, and Lorentz/gravitational factors; (2) explicit chiral $SU(2)_L$ parity violation; and (3) fermionic spin-$1/2$ exchange antisymmetry via octonionic triality ($8_v, 8_s, 8_c$).

---

## 7. Reproducibility and Released Files

All numerical results and figures are deterministically reproducible without randomness:
- Verification script: `tri_loop_n64_quantum_foundations_verify_v0.py` (~15 s)
- Numerical log: `results/tri_loop_n64_quantum_foundations_verify_v0.txt`
- Figure generator: `plot_fourth_report_figures_v0.py`
- Figures: `results/fig_r4_complex_and_born_{ja,en}.png` and `results/fig_r4_uncertainty_and_cutoff_{ja,en}.png`
- TeX builder: `build_tex_rapid4.py`

## References

1. N. Kihara, *Self-Consistent Crystal Selection, Dimensional Melting, and Two-Scale Particle Condensation in the Triangle-Loop Product of 3D Cubic Arrays (Rapid Communication)*, Zenodo (2026), DOI: 10.5281/zenodo.24616931.
2. N. Kihara, *Four-Stage Doubling-Entanglement Phase Transitions, Single-Clock Complex Locking, Exact Triple-Outer-Product Bell Pairs, and Probe-Type Branching in the 3D Triangle-Loop Product (Second Rapid Communication)*, Zenodo (2026), DOI: 10.5281/zenodo.24618400.
3. N. Kihara, *SO(3) Spin Selection, G2/SU(3) Gauge Symmetries, Spatial Harmonic Localization, and Two-Particle Coexistence in the 3D Triangle-Loop Product (Third Rapid Communication)*, Zenodo (2026), DOI: 10.5281/zenodo.24619357.
4. A. M. Gleason, *Measures on the Closed Subspaces of a Hilbert Space*, J. Math. Mech. **6**, 885–893 (1957).
5. D. L. Donoho and P. B. Stark, *Uncertainty Principles and Signal Recovery*, SIAM J. Appl. Math. **49**, 906–931 (1989).
