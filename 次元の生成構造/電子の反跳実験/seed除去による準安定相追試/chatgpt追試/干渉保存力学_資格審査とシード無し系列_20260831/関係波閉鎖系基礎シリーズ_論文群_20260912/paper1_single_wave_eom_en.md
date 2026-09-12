# The Equation of Motion of a Single Wave in a Relational-Wave Closed System — Kuramoto-Type Phase Coupling on $L(K_N)$ and the Condition for a Self-Consistent Floor on the 90-Degree Phase-Difference Skeleton

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 1 of 10)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22728761
**Concept DOI:** 10.5281/zenodo.22728760
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22728761

---

## Abstract

Consider the closed-system map that places a complex relational wave $z_e\in\mathbb C$ on every edge of the complete graph $K_N$ and evolves it by a phase-only real antisymmetric generator (defined in the prior work [1]). The question of this paper arises from an empirical motivation: the floor produced by make_parent — which builds a self-consistent fixed point $KZ=i\sigma Z$ from a random seed — turned out, for all $N$, to be **naturally a 90-degree two-axis phase lattice** (Experiment 4-1, Study 28). We therefore do not ask "why 90 degrees"; taking 90 degrees as given, we ask **"can a highly symmetric floor that rotates self-consistently be constructed on the 90-degree phase-difference skeleton?"** First, decomposing the one-step map $z\to\exp(\Delta\tau K)z$ into amplitude and phase, each edge wave couples in Kuramoto form $\sin(\varphi_f-\varphi_e)$ only to edges sharing a vertex, and we derive

$$\frac{dr_e}{d\tau}=\frac12\sum_{f\sim e}r_f\sin\!\big(2(\varphi_f-\varphi_e)\big)\quad(\text{amplitude}),\qquad r_e\frac{d\varphi_e}{d\tau}=\sum_{f\sim e}r_f\sin^2(\varphi_f-\varphi_e)\quad(\text{phase}).$$

On the 90-degree skeleton (all phase differences are multiples of $90^\circ$) the second-harmonic amplitude drive $\sin(2\psi)$ **vanishes termwise**, giving $\dot r_e=0$ (amplitude freezing). The remaining self-consistency condition (uniform rotation $\dot\varphi_e=\omega$) reduces to the eigenvalue problem of the amplitude vector, $\omega r_e=\sum_{f\sim e}r_f\sin^2\psi_{ef}$ (i.e. $Wr=\omega r$). Moreover, with $D=\operatorname{diag}(e^{i\varphi_e})$, on the 90-degree skeleton the identity $\boxed{D^{-1}KD=iW}$ ($W_{ef}=A_{ef}\sin^2\psi_{ef}$ real symmetric, non-negative) holds, so the **complex self-consistency problem $Kz=i\omega z$ ($z=Dr$) is exactly equivalent to the eigenvalue problem $Wr=\omega r$ of a real symmetric non-negative integer matrix**. If $W$ is irreducible, Perron–Frobenius guarantees a unique positive-amplitude principal eigenvector ($\omega=\rho(W)$), and $z=Dr$ constructs the self-consistent wave on the 90-degree skeleton (verified irreducible for all $N=3$–$40$). Its algebraic exact solution is given by Paper 3 (the make_parent floor is an iterative approximation of the Perron eigenvector of the $W$ fixed by its own 90-degree labels, a different member from the integer-$K$ floor, with different labels and $W$). Placed side by side at the same $N$, the make_parent floor has amplitudes that differ across almost all edges ($699$ values at $N=40$), whereas the intentionally constructed high-symmetry floor has $\le 3$ amplitude values and a markedly higher symmetry (Fig. 2). Orbit reproduction of the exact map is $\le3.0\times10^{-15}$; finite-difference verification of the EOM is $\le1.4\times10^{-6}$ (amplitude) and $\le2.4\times10^{-7}$ (phase); amplitude freezing on the 90-degree skeleton is $\le1.1\times10^{-13}$ ($N=6,22,40$). Classification: a derived consequence of the canonical one_step, with machine-precision verification.

---

## 1. Introduction

This series studies the dynamics of a closed system in which the full pairwise relations $M=N(N-1)/2$ among $N$ entities are complex relational waves $z_e$, with the external seed removed. The **make_parent** procedure used initially builds a self-consistent fixed point ($KZ=i\sigma Z$) from a random seed by iteration, and its floor turned out, for all $N=3$–$40$, to be **naturally a 90-degree two-axis phase lattice** (relative phase difference $90^\circ$; Fig. 1, Experiment 4-1, Study 28). The phase is naturally $90^\circ$, but the amplitudes are edge-by-edge scattered ($699$ values at $N=40$, coefficient of variation of a few percent) and the centroid $\Sigma z\neq0$, so the symmetry is not high.

![Fig. 1: Complex plane of the random make_parent floor (self-consistent fixed point KZ=iσZ from a random seed) for all N=3-40, step 0. For all N every wave lies on two orthogonal axes (relative phase difference 90 degrees) = naturally a 90-degree two-axis lattice. Amplitudes are scattered edge by edge.](../make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/fig_complex_plane_step0_makeparent_N3_N40.png)

**Fig. 1** Complex plane of the random make_parent floor ($N=3$–$40$, step 0). The self-consistent fixed point built from a random seed lies, for all $N$, naturally on two orthogonal axes (relative phase difference $90^\circ$). This is the starting point of this paper — the empirical fact that "90 degrees is given"; the reason for it (90 degrees is forced by self-consistency plus termwise freezing) is given by Paper 0.

Taking this empirical fact as the starting point, the question of this paper is not "why 90 degrees" but — **taking 90 degrees as given (it turned out that way naturally)** — **"can a highly symmetric floor that rotates self-consistently be constructed within the 90-degree phase-difference skeleton?"** Indeed, placed side by side at the same $N$ (Fig. 2), against the random make_parent floor (amplitudes nearly $M$ distinct values) the intentionally constructed high-symmetry floors have few amplitude values (high-symmetry v2: $\le 2$ values, $D_N$-symmetric / integer-$K$ analytic floor: $\le 3$ values, small centroid, relative equilibrium), with markedly higher symmetry.

![Fig. 2: Symmetry comparison of floors (same N, phase all 90 degrees). Left = random make_parent floor (amplitudes nearly M values, centroid ≠ 0); middle = high-symmetry v2 floor (amplitudes ≤ 2 values but largest-class centroid); right = integer-K analytic floor (≤ 3 values, small centroid, relative equilibrium = Paper 3).](../位相ロック全N確認_相対平衡_20260912/床の対称性比較_makeparent_vs_高対称_20260912/fig_floor_symmetry_compare.png)

**Fig. 2** Symmetry comparison of floors ($N=6,12,17,40$). Against the 90-degree floor generated naturally by random make_parent (amplitudes nearly $M$ values, $699$ at $N=40$), the intentionally constructed high-symmetry floors have few amplitude values and higher symmetry (details in Papers 2, 3).

## 2. Dynamics (the map)

The one_step of the canonical program `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353…`) is, with $u=e^{i\arg z}$, $H_{ef}=A_{ef}\,u_f\bar u_e$, $H\leftarrow i\,\Im H$: $z_{n+1}=\exp(-i\Delta\tau H)z_n$. Rearranged,

$$K_{ef}=A_{ef}\sin(\varphi_f-\varphi_e),\qquad z_{n+1}=\exp(\Delta\tau K)\,z_n,\qquad \Delta\tau=\frac{2\pi}{N},$$

where $A$ is the adjacency matrix of the line graph $L(K_N)$ of $K_N$ (edges $e,f$ share a vertex $\Rightarrow A_{ef}=1$) and $\varphi_e=\arg z_e$. **$K$ is determined by phases only and contains no amplitude; it is a real antisymmetric matrix**, so $\exp(\Delta\tau K)$ is a real orthogonal matrix. Acting identically on $\Re z$ and $\Im z$, it exactly conserves $\|z\|^2$ and the square-closure $z^{\mathsf T}z=\sum_e z_e^2$.

## 3. Equation of motion of a single wave (derivation)

Writing the generator flow $dz/d\tau=Kz$ in components, for $f\sim e$ (edges sharing a vertex),

$$\frac{dz_e}{d\tau}=\sum_{f\sim e}\sin(\varphi_f-\varphi_e)\,z_f=\sum_{f\sim e}r_f\sin(\varphi_f-\varphi_e)\,e^{i\varphi_f}.$$

With $z_e=r_e e^{i\varphi_e}$, writing the LHS as $e^{i\varphi_e}(\dot r_e+ir_e\dot\varphi_e)$ and using $e^{i\varphi_f}=e^{i\varphi_e}e^{i\psi}$ ($\psi=\varphi_f-\varphi_e$) on the RHS to cancel $e^{i\varphi_e}$ gives $\text{RHS}=\sum r_f\sin\psi(\cos\psi+i\sin\psi)$. From the real and imaginary parts,

$$\boxed{\ \dot r_e=\frac12\sum_{f\sim e}r_f\sin(2\psi)\ }\qquad\boxed{\ r_e\dot\varphi_e=\sum_{f\sim e}r_f\sin^2\psi\ }$$

(using $\sin\psi\cos\psi=\tfrac12\sin2\psi$). Each edge wave is an oscillator coupled in Kuramoto form only to the adjacent edges sharing a vertex.

## 4. This equation explains the orbit

- **Phase (displacement modulo rotation and rigid rotation):** since $\sin^2\ge0$, every wave $\varphi_e$ advances monotonically = rotation. The common part is a rigid rotation; the **difference** of each wave's $(1/r_e)\sum r_f\sin^2\psi$ is the phase displacement with the rigid rotation removed.
- **Amplitude (growth/decay):** the drive is the second harmonic $\sin(2\psi)$ — a weighted sum of the sine of twice the phase difference with adjacent edges.
- **Floor (90-degree skeleton):** all phase differences are multiples of $90^\circ$ $\Rightarrow\sin(2\psi)=0\Rightarrow\dot r_e=0$ (termwise amplitude freezing). Further, in a configuration where the amplitudes satisfy $Wr=\omega r$ ($W_{ef}=A_{ef}\sin^2\psi_{ef}$), all waves rotate at the same rate $\dot\varphi_e=\omega$, forming a self-consistent relative equilibrium. The existence and uniqueness of this floor is treated in §5.

The subject of this paper is the **construction** of the self-consistent floor on the 90-degree skeleton; that the floor is unstable (the qualification for ignition) and the dynamics after leaving the floor (inflation, re-locking) are not treated here (Papers 2, 5, 6).

## 5. Amplitude freezing and the self-consistency condition on the 90-degree skeleton

90 degrees is a given skeleton (the make_parent self-consistent floor turned out that way naturally, §1). On top of that, this section shows that the 90-degree skeleton freezes the amplitudes termwise, and that the condition for a self-consistent floor becomes an eigenvalue problem of the amplitudes.

**Amplitude freezing (termwise):** the amplitude drive is the second harmonic $\sin(2\psi)$, whose zeros on $[0,360^\circ)$ are only $\psi=0,90,180,270^\circ$ (by enumeration). Hence if all phase differences are multiples of $90^\circ$, then $\sin(2\psi)=0$ for every pair, i.e. $\dot r_e=0$ holds **term by term** (termwise). This is a direct consequence of the generator being $\sin\psi$ (first harmonic), so that the amplitude drive is the second harmonic; it is why the 90-degree skeleton is the skeleton on which amplitudes "can be placed without being destroyed." (Even at general phases, $\dot r=0$ can occur if the per-edge sum cancels = cancellation type; this corresponds to the lock endpoint.)

**Theorem (real-symmetrization of the 90-degree skeleton).** Imposing a self-consistent single rotation ($\dot\varphi_e=\omega$ uniform) on the 90-degree skeleton, the phase equation becomes $\omega r_e=\sum_{f\sim e}r_f\sin^2\psi_{ef}$, i.e. $Wr=\omega r$ ($W_{ef}=A_{ef}\sin^2\psi_{ef}$). This is more than an eigenvalue problem merely "appearing" — it is an **exact equivalence**: for $\psi_{ef}\in\frac{\pi}{2}\mathbb Z$ and $D=\operatorname{diag}(e^{i\varphi_e})$, the identity $\sin\psi\,e^{i\psi}=\tfrac12\sin2\psi+i\sin^2\psi=i\sin^2\psi$ (⟺ $\sin2\psi=0$) gives
$$\boxed{\ D^{-1}KD=iW\ },\qquad W_{ef}=A_{ef}\sin^2\psi_{ef}\in\{0,1\}\ (\text{real symmetric, non-negative}).$$
Therefore $\boxed{\,Kz=i\omega z\ (z=Dr,\ r\ \text{real})\iff Wr=\omega r\,}$. The complex self-consistent eigen-wave problem on the 90-degree skeleton is exactly equivalent to the eigenvalue problem of the real symmetric non-negative integer matrix $W$.

**Corollary (existence and uniqueness of the positive-amplitude self-consistent wave).** If $W$ is irreducible (its associated graph is connected), then by Perron–Frobenius there exists, up to scale, a unique $r$ with $Wr=\rho(W)r$, $r_e>0$. Hence $\boxed{\,z=Dr\,}$ is a self-consistent single rotation $Kz=i\rho(W)z$ on the 90-degree skeleton. That is, the answer to this paper's question "can a self-consistent positive-amplitude wave be built on the 90-degree skeleton?" is, under irreducibility of $W$, **yes (unique up to scale)**.

**Verification (all $N=3$–$40$, follow-up).** For both the integer-$K$ floor and the make_parent floor, $\|D^{-1}KD-iW\|\le2\times10^{-12}$; the $\omega=\rho(W)=\sigma$ ($iK$ spectral radius) of $Wr=\omega r$ agrees; $r>0$; the overlap between $r$ and the $W$ principal eigenvector is $1.000000$; and $W$ is irreducible (connected) for all $N$. In summary, $\boxed{\text{90-degree skeleton}\overset{D^{-1}KD=iW}{\longrightarrow}\text{real non-negative symmetric eigenproblem}\overset{\rm PF}{\longrightarrow}\text{existence/uniqueness of positive-amplitude self-consistent wave}}$.

The **closed form** of that high-symmetry solution (the $N$-dependence of $\rho(W)$, the 2/3-value degeneracy of the amplitudes, and the agreement with $\sigma_{\max}$ of the integer $K$) is solved concretely by Paper 3. The make_parent floor approximates the Perron principal eigenvector of the $W$ ($W_{\rm mp}$) fixed by its own 90-degree phase labels obtained through iteration (the sum rule $(\sum_f r_f\sin^2\psi_{ef})/r_e=\sigma$ holds for all $N$, Study 31); it is a different member from the integer-$K$ floor, with different labels and different $W$, yet consistent in that each is the Perron solution of its own $W$ ($W_{\rm mp}\neq W_{\rm integerK}$, Paper 0).

That is, 90 degrees is not a derived necessity but a given skeleton, and what this paper establishes is that "on the 90-degree skeleton the amplitudes freeze termwise, and the existence of a self-consistent floor reduces to the eigenvalue problem $Wr=\omega r$ of the amplitudes." The qualification for ignition of the floor (instability) is treated in Paper 2, the exact solution in Paper 3, and the rapid expansion from it in Papers 5, 6.

## 6. Numerical verification ($N=6,22,40$)

| Test | Content | Max error |
|---|---|---|
| T1 | The exact map $\exp(\Delta\tau K(\varphi_n))z_n$ reproduces the canonical orbit $z_{n+1}$ ($K$ is the generator) | $\le3.0\times10^{-15}$ |
| T2 | The amplitude/phase components of the substep $(\exp(\varepsilon K)z-z)/\varepsilon$ agree with the EOM of §3 | amplitude $\le1.4\times10^{-6}$ / phase $\le2.4\times10^{-7}$ |
| T3 | On the 90-degree lattice of the floor, $\dot r_e\approx0$ (amplitude freezing) | $|\dot r|\le1.1\times10^{-13}$ |

(T2 is due to rounding of the $\varepsilon$-scale finite difference; the formula is exact.) Measured 90-degree skeleton vs. lock endpoint: the 90-degree skeleton (step 0) has individual terms $|\sin2\psi|\sim10^{-12}\text{–}10^{-14}$ (termwise, lattice deviation $0.0^\circ$); the lock endpoint has individual terms $\sim1.0$ but per-edge weighted sums $\sim10^{-7}\text{–}10^{-13}$ (cancellation, lattice deviation $\sim45^\circ$).

![Fig. 3: Verification of the single-wave EOM (N=6). Left = amplitude, right = phase; measured (vertical) and derived formula (horizontal) agree on y=x.](../位相ロック全N確認_相対平衡_20260912/単一波運動方程式_導出と検証_20260912/fig_single_wave_eom_verify.png)

**Fig. 3** Verification of the single-wave equation of motion ($N=6$).

![Fig. 4: Amplitude freezing (termwise) on the 90-degree skeleton and locking (cancellation). Left = zeros of sin(2ψ) = k·90 degrees; right = measured distributions of the 90-degree skeleton (termwise, individual terms ~0) and locking (cancellation, individual terms O(1)) for N=6.](../位相ロック全N確認_相対平衡_20260912/単一波運動方程式_導出と検証_20260912/fig_90deg_necessity.png)

**Fig. 4** The 90-degree skeleton (termwise amplitude freezing) and the lock endpoint (cancellation), measured ($N=6$).

## 7. Claim classification / open items

- Classification: a **derived consequence** of the canonical one_step (plus machine-precision verification). Verdict: retained.
- Assertion (main result of this paper): on the 90-degree skeleton $D^{-1}KD=iW$ holds, and the complex self-consistency problem $Kz=i\omega z$ is **exactly equivalent** to the eigenvalue problem $Wr=\omega r$ of the real symmetric non-negative integer matrix $W_{ef}=A_{ef}\sin^2\psi_{ef}$. If $W$ is irreducible, Perron–Frobenius makes the positive-amplitude principal eigenvector unique ($\omega=\rho(W)$) = **a self-consistent positive-amplitude floor can be constructed on the 90-degree skeleton** (verified irreducible at machine precision for all $N=3$–$40$). This is the mathematical answer to "can a self-consistent high-symmetry wave be built on the 90-degree skeleton?"
- Open: the overall scale $r^2$ of the floor (scale is privately held), the moduli of relative equilibria (the 90-degree floor is the most symmetric member, not the unique solution), the $N=3$ exception ($\mu/r^2=-3/2$, rank 1), and the general condition for irreducibility of $W$ at an arbitrary 90-degree label configuration. These are treated in Papers 2, 3, 5.

## Related work (structural coincidence, not a derivation source)

- Kuramoto model (Kuramoto 1975/1984): $K_{ef}=A_{ef}\sin(\varphi_f-\varphi_e)$ is a Kuramoto-type phase coupling on $L(K_N)$.

## References

1. Noriaki Kihara, "The Mechanism of Inflationary Rapid Expansion in a Self-Consistent Relational-Wave Closed System," Concept DOI 10.5281/zenodo.22112008. (Source of the definition of the dynamics of this paper and of $\tau\neq$ time.)
2. Y. Kuramoto, *Chemical Oscillations, Waves, and Turbulence*, Springer (1984); "Self-entrainment of a population of coupled non-linear oscillators," 1975.

## Reproduction

The verification scripts and data are packaged in full in `位相ロック全N確認_相対平衡_20260912/単一波運動方程式_導出と検証_20260912/` (`verify_single_wave_eom_20260912.py` = T1/T2/T3, `verify_90deg_necessity_20260912.py`, `make_figures_eom_90deg_20260912.py`, `REPORT.md`, `SHA256SUMS.txt`, `run_all.sh`). The **all-$N$ verification of the complex↔real-symmetric equivalence $D^{-1}KD=iW$ and Perron–Frobenius (§5)** is in `位相ロック全N確認_相対平衡_20260912/複素実対称同値_DKD_iW_20260912/` (`verify_DKD_iW_equivalence_20260912.py`, `results/DKD_iW_summary.csv`, README/REPORT/SHA/run_all; read-only). The symmetry comparison of Fig. 1 is in `同/床の対称性比較_makeparent_vs_高対称_20260912/`. The input is the existing den=N npz ($N=6$: 500 step, $N=22,40$: 2000 step of this series). The dynamical map is identical to the canonical `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`), with no change to the physics equations.
