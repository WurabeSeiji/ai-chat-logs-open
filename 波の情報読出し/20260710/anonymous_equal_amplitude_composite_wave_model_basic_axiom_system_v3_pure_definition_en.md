# Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v3

**Date:** 2026-07-12  
**Author:** Noriaki Kihara  
**Status:** Definition paper  
**Version DOI:** 10.5281/zenodo.21315736
**Concept DOI:** 10.5281/zenodo.21315735

---

# 1. First Principles

## Axiom 0: Component Anonymity

Classification: Axiom

Basic components are not assigned individual names, privileged axes, privileged signs, or privileged types.

The index `n` is a convenient label and is not an intrinsic name of the component.

Forbidden:

```text
making only one specific component the negative-sign axis
making only one specific component the time axis
making only one specific component real-valued
making only one specific component complex-valued
giving an intrinsic name to only one specific component
```

---

## Axiom 0+: Formulation Anonymity

Classification: Axiom

At the first-principle layer, theory names, formulation names, and readout-logic names are not assigned intrinsic names.

Equations, projections, or readout rules at the first-principle layer must not be changed on the basis of physical names, object names, interaction names, or existing theory classifications.

Forbidden:

```text
adopting a circular-motion equation because it is called gravity
adopting a hyperbolic equation because it is called Coulomb force
erasing signs because it is called mass
preserving signs because it is called charge
moving a component to a time axis because it is called energy
using a different closure equation because the target is a two-body system
using a different closure equation because the target is a three-body system
using a special readout rule because the target is an ABC system
choosing a right-hand-side representation to obtain a desired physical name
```

Allowed:

```text
different readout operations applied to the same closure condition
different symmetry breaking applied to the same closure condition
different gauge stability applied to the same closure condition
different sign preservation applied to the same closure condition
different sign erasure applied to the same closure condition
different phase-branch selection applied to the same closure condition
different closure target sets applied to the same closure condition
```

Order:

```text
place the same all-positive-sign zero closure
define the readout operation first
define symmetry breaking first
define sign preservation or sign erasure first
confirm gauge stability
assign working names after readout
```

Target set:

```math
Q(P)=\sum_n p_n^2
```

```math
Q(A)=0
```

```math
Q(A+B)=0
```

```math
Q(A+B+C)=0
```

Differences in target set are allowed.

Changing the definition of `Q` is not allowed.

---

## Axiom 0.5: Dimension-Wise Central Projection Readout

Classification: Auxiliary geometric axiom

No curvature correction is introduced for readout along a single axis.

Curvature-correction candidates are introduced only when two or more independent directions form an area, volume, closed path, or transport mismatch.

---

## Axiom 1: All-Positive-Sign Zero Closure

Classification: Axiom

```math
\sum_{n=1}^{N}x_n^2=0
```

Expanded form:

```math
x_1^2+x_2^2+\cdots+x_N^2=0
```

Coefficient signs:

```text
all terms +
```

Not adopted:

```math
\sum_n |x_n|^2=0
```

Term used for closure:

```math
x_n^2
```

Term not used for closure:

```math
x_n\bar{x}_n
```

Applicable states:

```text
stationary closed state
stable closed state
```

Quasi-stationary states:

```text
immediately after external influence
closure recovery process
metastable state
```

Temporary nonzero closure residual in a quasi-stationary state:

```math
\sum_n x_n^2\ne0
```

After stationarization:

```math
\sum_n x_n^2=0
```

---

## Axiom 2: Nontrivial Existence

Classification: Axiom

```math
\exists n,\quad x_n\neq0
```

---

# 2. Consequences from the First Principles

## Consequence 1: Insufficiency of Real Positive-Definite Components

Classification: Derived consequence

```math
x_n\in\mathbb{R}
\Rightarrow
x_n^2\ge0
```

```math
\sum_n x_n^2=0
\Rightarrow
x_n=0\quad\forall n
```

---

## Consequence 2: Internal Generation of Sign Reversal

Classification: Derived consequence

```math
x_1=A,\qquad x_2=iA
```

```math
x_1^2+x_2^2=A^2+(iA)^2=A^2-A^2=0
```

```math
i^2=-1
```

---

## Consequence 3: Necessity of Phase Algebra

Classification: Derived consequence

```math
x_n=r_n e^{i\phi_n}
```

```math
x_n^2=r_n^2e^{i2\phi_n}
```

Closure condition:

```math
\sum_{n=1}^{N}r_n^2e^{i2\phi_n}=0
```

Equal-amplitude condition:

```math
r_n=r
```

```math
\sum_{n=1}^{N}e^{i2\phi_n}=0
```

---

## Consequence 4: Type Anonymity

Classification: Derived consequence

```math
x_n\in\mathbb{C}\quad\forall n
```

or

```text
all components are assigned the same phase-algebra type
```

---

## Consequence 5: Ninety-Degree Phase Difference

Classification: Derived consequence

```math
A^2+(iA)^2=0
```

---

# 3. Anonymous Equal-Amplitude Composite Wave

## Axiom 3: Anonymous Equal-Amplitude Composite Wave

Classification: Axiom

```math
\Psi_N=\frac{A_0}{N}\sum_{k=1}^{N}e^{i\theta_k}
```

```math
|\psi_k|=\frac{A_0}{N}
```

When all phases are aligned:

```math
\max|\Psi_N|=A_0
```

---

## Axiom 4: Unreadability of Individual Phase

Classification: Axiom

```math
\theta_k\rightarrow\theta_k+\alpha
```

---

## Axiom 5: Phase-Difference Readout

Classification: Axiom

```math
\Delta_{ij}=\theta_i-\theta_j
```

Number of readable pair channels:

```math
\binom{N}{2}=\frac{N(N-1)}{2}
```

---

## Axiom 6: System Identification

Classification: Axiom

```math
(N,A_0,\{\Delta_{ij}\})
```

---

## Axiom 7: Composite Phase Projection

Classification: Axiom

```math
S_N=\sum_k e^{i\theta_k}
```

```math
\Psi_N=\frac{A_0}{N}S_N
```

```math
\Theta_N=\arg\Psi_N
```

---

## Axiom 8: No Requirement of Inversion Symmetry

Classification: Axiom

Inversion symmetry is not required.

---

# 4. External Readout

## Axiom 9: External-Axis Response

Classification: Axiom candidate

```math
M_\alpha=e^{i\alpha}
```

```math
Y(\alpha)=\frac{A_0}{N}\sum_{k=1}^{N}e^{i(\theta_k-\alpha)}
=e^{-i\alpha}\Psi_N
```

Continuous response:

```math
S_N(\alpha)=\operatorname{Re}\left[\frac{1}{N}\sum_{k=1}^{N}e^{i(\theta_k-\alpha)}\right]
```

Binary response:

```math
S_N(\alpha)=\operatorname{sgn}\operatorname{Re}\left[\frac{1}{N}\sum_{k=1}^{N}e^{i(\theta_k-\alpha)}\right]
```

---

## Axiom 10: External Readout Phase

Classification: Axiom candidate

```math
\Delta\phi_\mu=\phi_\mu^{(\Psi)}-\phi_\mu^{(M)}
```

```math
\mu=x,y,z,t
```

---

## Axiom 11: Position-Phase Readout

Classification: Axiom candidate

```math
x,y,z
```

---

## Axiom 12: Time-Phase Readout

Classification: Axiom candidate

```math
\phi_t\equiv\omega t_{\mathrm{read}}\pmod{2\pi}
```

---

## Axiom 13: Covering Phase

Classification: Axiom candidate

```text
2pi-periodic readout
4pi-periodic readout
winding-number-difference readout
```

---

# 5. Observational Phase Discreteness

## Axiom 14: Observational Phase Discreteness

Classification: Axiom candidate

```math
\sum_{\mathrm{cycle}}\Delta\phi=2\pi m,\qquad m\in\mathbb{Z}
```

```math
W=\frac{1}{2\pi}\oint d\phi,\qquad W\in\mathbb{Z}
```

Readout-cell form:

```math
\phi_{\mathrm{obs}}=m\epsilon_\phi,\qquad m\in\mathbb{Z}
```

Phase-area candidate:

```math
Q_\phi=\sum_{i<j}\sin^2\frac{\Delta_{ij}}{2}
```

---

# 6. Waveform and Localization

## Working Axiom A: State Waveform

Classification: Working hypothesis

```math
A(\theta)=\sum_n A_n\cos(n\theta+\phi_n)
```

```math
Z(\theta)=\sum_n r_n e^{i(n\theta+\phi_n)}
```

---

## Definition A1: Normalized Equal-Amplitude Odd-Harmonic Localized Wave

Classification: Definition

```math
N_h=2K-1
```

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

```math
Z_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}e^{i(2m+1)u}
```

```math
\Psi_{N_h}(u)=A S_{N_h}(u)
```

```math
\Psi_{N_h}(u)=A Z_{N_h}(u)
```

```math
\frac{A}{K}
```

---

## Theorem A1: Representative-Amplitude Preservation

Classification: Derived theorem

```math
S_{N_h}(0)=1
```

```math
\Psi_{N_h}(0)=A
```

---

## Theorem A2: Highest Odd-Harmonic Order and Localization Width

Classification: Derived theorem

```math
K=\frac{N_h+1}{2}
```

```math
\lambda_{N_h}=\frac{\lambda_0}{N_h}
```

```math
\Delta x_{N_h}\propto\frac{\lambda_0}{N_h}
```

```math
N_{h,\rm req}\sim\frac{\lambda_0}{\Delta x}
```

```math
N_{h,\rm req}^{(t)}\sim\frac{T_0}{\Delta t}
```

---

## Consequence A1: Identical Construction of Spatial and Temporal Localization

Classification: Derived consequence

```math
\Psi(\chi,\tau)=A S_{N_{h,\chi}}(\chi-\chi_0)S_{N_{h,\tau}}(\tau-\tau_0)e^{i\phi}
```

---

## Working Axiom B: Localization

Classification: Working hypothesis

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

```math
Z_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}e^{i(2m+1)u}
```

---

## Working Axiom C: Coupling Entrance

Classification: Working hypothesis

The entrance to interaction is low-order base-phase alignment.

---

## Working Axiom D: Observation Map

Classification: Working hypothesis

```math
s_\alpha=\mathrm{Loc}\left[\int A(\theta)M_\alpha(\theta)d\theta\right]
```

---

## Working Axiom E: Squared Correlation Strength

Classification: Working hypothesis

```math
P(\alpha)=\frac{|C_\alpha|^2}{\sum_\beta |C_\beta|^2}
```

---

## Working Axiom F: Shared Low-Order Base Synchronization

Classification: Working hypothesis

```math
\theta_A-\theta_B=\Delta
```

---

## Working Axiom G: Synchronized Demodulation Failure

Classification: Working hypothesis

```math
V(t)=\left|\sum_n w_n e^{i\epsilon_n(t)}\right|
```

```math
V(t)=\left|\sum_n w_n e^{-\frac12\sigma_n^2(t)}\right|
```

---

## Working Axiom H: Hierarchy

Classification: Working hypothesis

```text
low-order base mode
+ high-order localized mode
```

---

# 7. Readout, Memory, and Complex Representation

## I/Q Readout

Classification: Theoretical observation

```math
z=a+ib
```

```math
a=r\cos\theta,\qquad b=r\sin\theta
```

```math
u=r\cos\theta+r\sin\theta
=\sqrt2 r\cos\left(\theta-\frac{\pi}{4}\right)
```

---

## Harmonic-Information Readout

Classification: Theoretical observation

```math
Z(\theta)=\sum_n r_n e^{i(n\theta+\phi_n)}
```

```math
P(Z)=\operatorname{Re}(Z)+\operatorname{Im}(Z)
```

```math
P(Z)=\sum_n \sqrt2 r_n
\cos\left(n\theta+\phi_n-\frac{\pi}{4}\right)
```

---

## Preservation Side and Readout Side

Classification: Theoretical observation

```math
z\bar z=a^2+b^2
```

```math
z^2=a^2-b^2+2iab
```

---

## Square Registry

Classification: Working hypothesis

```math
\mathcal{A}_0=\sum_k A_k^2
```

```math
a_0=\sum_k |A_k|^2=\sum_k A_k\bar A_k
```

---

# 8. Statistical Classification

## Pair Phase Difference

Classification: Definition

```math
\Delta_{ij}=\theta_i-\theta_j
```

```math
\binom{N}{2}=\frac{N(N-1)}{2}
```

Pair composite:

```math
\psi_i+\psi_j
=2A\cos\frac{\Delta_{ij}}{2}e^{i(\theta_i+\theta_j)/2}
```

Pair interference strength:

```math
I_{ij}=4A^2\cos^2\frac{\Delta_{ij}}{2}
```

---

## Overlap Degree

Classification: Definition

```math
R_N=\left|\frac{A_0}{N}\sum_{k=1}^{N}e^{i\theta_k}\right|
```

```math
\rho_N=\frac{R_N}{A_0}
=\frac{1}{N}\left|\sum_{k=1}^{N}e^{i\theta_k}\right|
```

```math
0\le\rho_N\le1
```

---

## Fermionic Degree

Classification: Definition

```math
F_N=
\frac{2}{N(N-1)}
\sum_{i<j}\sin^2\frac{\Delta_{ij}}{2}
```

```math
B_N=1-F_N
```

```math
B_N=
\frac{2}{N(N-1)}
\sum_{i<j}\cos^2\frac{\Delta_{ij}}{2}
```

```math
F_N=\frac{N}{2(N-1)}(1-\rho_N^2)
```
