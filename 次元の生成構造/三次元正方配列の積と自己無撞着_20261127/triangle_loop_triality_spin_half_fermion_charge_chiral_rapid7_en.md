# Octonionic Spin(8) Triality Spin-1/2, CAR Fermion Charge Spectrum (0, 1/3, 2/3, 1), and Time-Crystal Chiral Splitting in the 3D Triangle-Loop Product (Seventh Rapid Communication)
## ― Exact 28D Triality Sign Flip Under $2\pi$ Rotation, Pauli Exclusion and $SU(3)$ Casimirs in a 1:3:3:1 Hierarchy, and Inner-$M_8$ Chiral Asymmetry ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)
**ORCID:** 0009-0004-6753-4020
**Date:** November 28, 2026
**Status:** Seventh Rapid Communication following Rapid Reports 1–6 (DOIs: 10.5281/zenodo.24616931, 24618400, 24619357, 24620262, 24621139, 24621454). Establishes spin-$1/2$ sign flips via $Spin(8)$ triality, canonical anticommutation relations (CAR) with exact one-generation Standard Model electric charges ($0, 1/3, 2/3, 1$) and $SU(3)$ Casimirs, and dynamical chiral splitting driven by the inner time crystal $M_8$.

---

## Abstract

We analyze the internal algebraic and dynamical fermionic structures of the outer octonion crystal $O_8$ and Layer-1 time crystal $H = O_8 \otimes M_8$ selected by the background-free real 3-tensor triangle-loop product $\mathrm{Tri}(T)_{ijk} = \sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}$ at $N = 64$ (`tri_loop_n64_triality_charge_v0.py`, zero randomness). Three exact results are established to machine precision:

1. **28D $Spin(8)$ Triality and Spin-$1/2$ Sign Flip ($-I$ under $2\pi$ Rotation)**: The 3-index triality invariance system ($512\times 84$) on $O_8[i, j, k]$ has `rank = 56` and an exact invariant kernel of dimension **`28`** (max zero singular value `4.489e-15`, finite rotation error `8.934e-16`). Under a $2\pi$ rotation in the spatial $(e_1, e_2)$ plane, the vector representation returns to $+I$ (`7.622e-16`), whereas both left- and right-handed spinors flip sign strictly to **$R_s(2\pi) = -I$ (`1.584e-14`) and $R_c(2\pi) = -I$ (`1.793e-14`)**, returning to $+I$ only after $4\pi$.
2. **CAR Fermion Algebra (Pauli Exclusion) and One-Generation Charge Spectrum ($1:3:3:1$)**: From the left-multiplication Clifford algebra $\{L_a, L_b\} = -2\delta_{ab}I$ (error **`0.000e+00`**) and the three color complex planes $(e_1,e_6), (e_2,e_5), (e_3,e_4)$ induced by fixing $e_7$, three fermionic ladder operators $\alpha_1, \alpha_2, \alpha_3$ satisfy exact canonical anticommutation relations $\{\alpha_m, \alpha_n\} = 0$ (Pauli exclusion $\alpha_m^2 = 0$) and $\{\alpha_m, \alpha_n^\dagger\} = \delta_{mn}I$ (both error **`0.000e+00`**). The charge operator $Q = \frac{1}{3}\sum_{m=1}^3 \alpha_m^\dagger \alpha_m$ commutes strictly with all 8 $SU(3)$ gluon generators ($\|[\Lambda_A, Q]\| = \mathbf{0.000\times 10^0}$) and decomposes the 8 states into $Q = 0$ (singlet, $C_2=0$: $\nu$-type), $Q = +1/3$ (triplet, $C_2=16/3$: $\bar{d}$-type), $Q = +2/3$ (triplet, $C_2=16/3$: $u$-type), and $Q = +1$ (singlet, $C_2=0$: $e^+$-type)—an exact **`1 : 3 : 3 : 1`** binomial degeneracy matching the $\mathbb{Z}_2^3$ spatial harmonic law of Rapid Reports 3 and 6.
3. **Dynamical Chiral Splitting Driven by Inner $M_8$**: On $H = O_8 \otimes M_8$ ($\tau = 30$), spatial parity inversion ($O_8 \leftrightarrow O_8^T$) leaves $\|P\|$ invariant, whereas reversing the inner time-crystal rotation ($M_8 \leftrightarrow M_8^T$) splits $\|P\|$ from `0.046873` to `0.046212` ($-1.41\%$).

---

## 1. Introduction: From Bosons and Spacetime to Fermions and Electric Charge

Rapid Communications 1–6 [1] established Hurwitz crystal selection, four-stage entanglement cascades, $SO(3)$ spin selection, $G_2 \supset SU(3)$ gauge dimensions, four emergent quantum theorems, 4-class spherical metric distortion with critical collapse, and Hurwitz 4-stage mass hierarchies. Three foundational questions remained open: (i) the emergence of spin-$1/2$ and Pauli exclusion, (ii) the algebraic origin of fractional/integer electric charges ($0, 1/3, 2/3, 1$) and their link to the `1:3:3:1` $\mathbb{Z}_2^3$ spatial harmonics, and (iii) dynamical left/right chiral asymmetry. Here we resolve all three at machine precision.

---

## 2. 28D $Spin(8)$ Triality and Spin-$1/2$ Sign Flip

While identical rotations $(R, R, R)$ on all three indices of $O_8$ preserve only $G_2$ (14D), assigning distinct representations $(8_v, 8_s, 8_c)$ yields the triality invariance condition $\sum_p ((A_v)_{ip}O_8[p,j,k] + (A_s)_{jp}O_8[i,p,k] + (A_c)_{kp}O_8[i,j,p]) = 0$. SVD of this $512 \times 84$ system confirms `rank = 56` and an exact **`28`-dimensional** $Spin(8)$ kernel (Table 1, Fig. 1). Exponentiating the spatial $(e_1, e_2)$ generator demonstrates exact spin-$1/2$ behavior: $R_v(2\pi) = +I$, while $R_s(2\pi) = R_c(2\pi) = -I$.

| Verification Item ($Spin(8)$ Triality) | Residual Norm | Exact Property |
| :--- | :---: | :--- |
| Triality kernel dimension ($512\times 84$ SVD) | `28` (max zero SV `4.489e-15`) | Full 28D $Spin(8)$ invariance |
| Finite triality rotation ($\theta = 0.7$) error | `8.934e-16` | $(R_v \otimes R_s \otimes R_c)O_8 = O_8$ |
| $2\pi$ rotation: $\|R_v(2\pi) - (+I)\|$ | `7.622e-16` | Vector $8_v$ (period $2\pi$) |
| **$2\pi$ rotation: $\|R_s(2\pi) - (-I)\|,\; \|R_c(2\pi) - (-I)\|$** | **`1.584e-14`, `1.793e-14`** | **Spin-$1/2$ sign flip ($-I$)** |
| $4\pi$ rotation: $\|R_s(4\pi) - (+I)\|,\; \|R_c(4\pi) - (+I)\|$ | `3.168e-14`, `3.587e-14` | Spinors $8_s, 8_c$ (period $4\pi$) |

![(a) Exact 28D $Spin(8)$ triality invariance and machine-precision residuals. (b) Phase sign along spatial $(e_1, e_2)$ rotation $0 \le \theta \le 4\pi$ contrasting vector $8_v$ (period $2\pi$) against spinors $8_s, 8_c$ (period $4\pi$, exact $-I$ flip at $2\pi$).](results/fig_r7_triality_spin_half_en.png)

---

## 3. CAR Fermion Algebra and One-Generation Charge Spectrum (`1:3:3:1`)

Left multiplication matrices $L_a$ on $O_8$ satisfy $\{L_a, L_b\} = -2\delta_{ab}I$ (`0.000e+00`). Fixing $e_7$ defines three color complex planes $(e_1,e_6), (e_2,e_5), (e_3,e_4)$ and three Furey ladder operators $\alpha_m = \frac{1}{2}(-L_a + i\,\mathrm{sgn}\cdot L_b)$ [2] obeying exact CAR ($\{\alpha_m, \alpha_n\} = 0$ and $\{\alpha_m, \alpha_n^\dagger\} = \delta_{mn}I$, both error **`0.000e+00`**). The electric charge $Q = \frac{1}{3}\sum_{m=1}^3 \alpha_m^\dagger \alpha_m$ commutes strictly with all 8 $SU(3)$ generators ($\|[\Lambda_A, Q]\| = \mathbf{0.000\times 10^0}$) and partitions the 8 states into one Standard Model generation (Table 2, Fig. 2(a)):

| Number $N$ | Charge $Q = N/3$ | Multiplicity | $SU(3)$ Casimir $C_2$ | One-Generation Assignment |
| :---: | :---: | :---: | :---: | :--- |
| $0$ | **`0`** | **1** | `0.000000` | Color singlet ($\nu$-type) |
| $1$ | **`+1/3`** | **3** | `5.333333` ($= 16/3$) | Color triplet ($\bar{d}$-type quark) |
| $2$ | **`+2/3`** | **3** | `5.333333` ($= 16/3$) | Color triplet ($u$-type quark) |
| $3$ | **`+1`** | **1** | `0.000000` | Color singlet ($e^+$-type) |

This exact **`1 : 3 : 3 : 1`** multiplicity shares the identical 3-bit $\mathbb{Z}_2^3$ binomial structure $\binom{3}{k}$ discovered in the spatial harmonics of Rapid Reports 3 and 6 [1].

---

## 4. Dynamical Chiral Splitting Driven by Inner Time Crystal $M_8$

Background $H$ has parity decomposition $E_{\mathrm{sym}} = 0.700395$ and $E_{\mathrm{asym}} = 0.299605$. Evolving four chiral/parity/CP probes ($\epsilon = 0.10$) on $H$ to $\tau = 30$ shows that spatial reversal ($O_8 \leftrightarrow O_8^T$) leaves $\|P\|$ unchanged, whereas reversing the inner time-crystal rotation ($M_8 \leftrightarrow M_8^T$) produces a clean $-1.41\%$ splitting (Table 3, Fig. 2(b)):

| Probe (Chiral / Parity / CP Test) | `nz` | $E_{\mathrm{probe}}$ | Norm $\|P\|$ at $\tau = 30$ | Outcome |
| :--- | :---: | :---: | :---: | :--- |
| $O_8 \otimes M_8$ (Co-chiral R,R) | 1024 | `0.009901` | **`0.046873`** | Reference (co-rotating) |
| **$O_8 \otimes M_8^T$ (Anti-chiral R,L)** | 1024 | `0.009901` | **`0.046212`** | **$-1.41\%$ chiral split** |
| $O_8^T \otimes M_8$ (Parity L,R) | 1024 | `0.009901` | **`0.046873`** | Invariant under spatial flip |
| **$O_8^T \otimes M_8^T$ (CP-flip L,L)** | 1024 | `0.009901` | **`0.046212`** | **$-1.41\%$ chiral split** |

![(a) Exact CAR fermion algebra and one-generation charge/color spectrum (error `0.000e+00`, `1:3:3:1` degeneracy). (b) Chiral splitting driven solely by inner time-crystal rotation $M_8$ vs. $M_8^T$ ($\Delta\|P\| = -1.41\%$).](results/fig_r7_fermion_charge_and_chiral_en.png)

---

## 5. Discussion and Outlook

Rapid Communication 7 completes the algebraic identification of spin-$1/2$ ($Spin(8)$ triality), Pauli exclusion (CAR), one-generation electric charges ($0, 1/3, 2/3, 1$) with $SU(3)$ Casimirs, and inner time-crystal chiral splitting. Extending to three generations, explicit $SU(2)_L$ weak interactions, and quantitative electron-photon scattering (fine-structure constant) remain clear targets for future work beyond the unified full paper integrating Rapid Communications 1–7.

---

## 6. Reproducibility and Released Files

- Verification script: `tri_loop_n64_triality_charge_v0.py`
- Numerical log: `results/tri_loop_n64_triality_charge_v0.txt`
- Figure generator: `plot_seventh_report_figures_v0.py`
- Figures (Fig. 1): `results/fig_r7_triality_spin_half_{ja,en}.png`
- Figures (Fig. 2): `results/fig_r7_fermion_charge_and_chiral_{ja,en}.png`
- TeX builder: `build_tex_rapid7.py`

## References

1. N. Kihara, *Rapid Communications 1–6*, Zenodo (2026), DOIs: 10.5281/zenodo.24616931, 24618400, 24619357, 24620262, 24621139, 24621454.
2. C. Furey, *Standard model physics from an algebra?*, Ph.D. thesis, University of Waterloo (2015), arXiv:1611.09182.
3. J. C. Baez, *The Octonions*, Bull. Amer. Math. Soc. **39**, 145–205 (2002).
