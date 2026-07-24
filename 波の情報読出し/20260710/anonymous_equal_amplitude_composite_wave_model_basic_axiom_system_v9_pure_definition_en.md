# Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v9

**Date:** 2026-07-24  
**Author:** Noriaki Kihara  
**Status:** Definition paper  
**Version DOI:** 10.5281/zenodo.21522310<br>
**Concept DOI:** 10.5281/zenodo.21315735

Changes in v9 (introduction of the anonymous partial-projection existence axiom): (1) Axiom 0.7 is added as a general principle stating that physically realized state space is the image of a nontrivial admissible partial projection selected from the candidate space of the foundational axioms. (2) The numbers and texts of the existing Axioms 0, 0.5, 0.6, and 1--17, as well as the working axioms, are unchanged. (3) No specific symmetry-group name, eight-component assumption, two-quartet structure, or derived number 240, 248, or 30 is introduced into the axioms. These remain a conditional identification under additional branch assumptions in the independent paper. Concept DOI: 10.5281/zenodo.21521899

---

# 1. First Principles

## Axiom 0: Anonymity

Classification: Axiom

This axiom consists of two clauses (Axiom 0 and Axiom 0+ of v7 and earlier, merged as clauses; no change to the text).

### Clause 1: Component Anonymity

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

### Clause 2: Formulation Anonymity

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

## Axiom 0.5: Scale Anonymity

Classification: Axiom

The relations, update rules, and readouts of the system must reference no absolute scale.

Every construction must be covariant under

```math
Z\to\lambda Z,\qquad\lambda\in\mathbb{C}^\ast
```

That is:

```text
relations (closure conditions) are homogeneous
the generator is zeroth-order homogeneous, K(λZ)=K(Z)
the state map is first-order homogeneous
readouts are ratio quantities only
```

No means exists to measure an absolute scale from inside the system.

This axiom is an unprovable formation rule of the same rank as Axiom 0 (anonymity). The two are clauses of a single idea: no names, no scale — only relations and ratios can be read.

The compact part `|λ|=1` of the covariance group `C*` is the common-phase invariance of the zero closure; the unreadability of phase (Axiom 5) and the unreadability of scale are the argument part and the modulus part of one and the same group.

The zero on the right-hand side of Axiom 1 is the unique choice compatible with this axiom within the family

```math
\sum_n x_n^2=c
```

The first introduction and the discussion are given in Paper 5. Concept DOI: 10.5281/zenodo.21486233

Theoretical background: writing `x_n=q_n+ip_n`, the imaginary part of Axiom 1 is the condition that the generator of scale transformations (the dilatation) `Σq_np_n` vanishes, so the identical elimination of the scale degree of freedom is built into the closure. The geometric organization including this observation is given in the expository note. Concept DOI: 10.5281/zenodo.21495305

---

## Axiom 0.6: Dimension-Wise Central Projection Readout

Classification: Auxiliary geometric axiom (Axiom 0.5 of v7 and earlier; renumbered upon placing scale anonymity at 0.5; no change to the text)

No curvature correction is introduced for readout along a single axis.

Curvature-correction candidates are introduced only when two or more independent directions form an area, volume, closed path, or transport mismatch.

---

## Axiom 0.7: Anonymous Partial-Projection Existence

Classification: meta-selection axiom for physical realization

An absolute-scale difference between states is not counted as a difference between candidates. By Axiom 0.5, set

```math
X\sim_{\mathrm{sc}}\lambda X,
\qquad
\lambda\in\mathbb C^\ast.
```

Let

```math
\mathscr S_{\mathrm{ax}}
```

be the set of nontrivial candidate systems satisfying the adopted mandatory axioms other than Axiom 0.7, modulo this scale equivalence.

For a branch in which working axioms or additional conditions are adopted, the candidate space of that branch is the subset satisfying those conditions.

A partial projection is a pair

```math
\mathcal P
=
\left(
\mathscr D_{\mathcal P},
p_{\mathcal P}
\right),
```

where

```math
\mathscr D_{\mathcal P}
\subseteq
\mathscr S_{\mathrm{ax}}
```

and

```math
p_{\mathcal P}:
\mathscr D_{\mathcal P}
\longrightarrow
\mathscr Y_{\mathcal P}.
```

Here “partial projection” is not restricted to an orthogonal projection in linear algebra. It may include selection retaining only part of the candidates as its domain, a quotient identifying unreadable differences, and a map sending readable relational quantities to a representation space.

A partial projection is nontrivial if at least one of the following holds:

```math
\mathscr D_{\mathcal P}
\subsetneq
\mathscr S_{\mathrm{ax}},
```

or

```math
\exists s_1\ne s_2
\quad\text{such that}\quad
p_{\mathcal P}(s_1)
=
p_{\mathcal P}(s_2).
```

A nontrivial partial projection is admissible when all of the following conditions hold:

```text
it is independent of individual component names and formulation names
it is independent of absolute scale
it is consistent with the foundational closure and the composition rule on the image
two states identified by the projection cannot be distinguished using only the adopted readouts
no physical name or symmetry-group name available only after projection is used in the projection condition
```

Write the set of admissible partial projections as

```math
\mathfrak P(\mathscr S_{\mathrm{ax}}).
```

The physically realized state space is the image of at least one nontrivial admissible partial projection:

```math
\boxed{
\exists\mathcal P_*
\in
\mathfrak P(\mathscr S_{\mathrm{ax}})
\quad\text{such that}\quad
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*.
}
```

This axiom does not assert:

```text
that the admissible partial projection is unique
that a functional or dynamics selecting the actual projection has already been obtained
that the foundational axioms alone uniquely force a specific symmetry group
```

Let the set of readable transition differences in the selected image be

```math
\Delta_{\mathcal P_*}.
```

The physical symmetry of the selected image is defined as the automorphism group preserving its readout relations, transition-difference set, and composition rule:

```math
G_{\mathrm{phys}}
:=
\operatorname{Aut}
\left(
\operatorname{Im}\mathcal P_*,
\Delta_{\mathcal P_*}
\right).
```

The order is therefore:

```text
foundational axioms
candidate space
anonymous admissible partial projection
physical partial image
subsequent identification of the symmetry of the image
```

No specific symmetry-group name is included in this axiom. The conditional identification of the `E8` lattice under the additional branch assumptions of a positive-definite eight-component readout, a two-quartet decomposition, local `D4` transition lattices, and a common-center neutrality projection is given in the independent paper. Concept DOI: 10.5281/zenodo.21521899

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

Theoretical background: the solution set of this closure is the complex isotropic cone, and nontrivial solutions require imaginary parts. Under the projectivization of Axiom 0.5, the state space becomes a projective quadric (a compact Kähler manifold), and compactness makes discrete spectra and finite-dimensional state spaces follow from the geometry. The expository identification by known mathematics is given in the expository note. Concept DOI: 10.5281/zenodo.21495305

---

## Axiom 2: Nontrivial Existence

Classification: Axiom

```math
\exists n,\quad x_n\neq0
```

---

## Axiom 3: Projection-Axis Degeneration and Imaginary-Direction Readout

Classification: Axiom

When an intrinsic coordinate system is constructed by a projection method, the projection-axis direction is not distinguishable as a direction for the residents of that intrinsic coordinate system.

However, if a value in the projection-axis direction remains inside the closure condition, that value is not erased.

When the compensating quantity in the projection-axis direction is read as `z`, it is treated inside the intrinsic coordinate system as

```math
iz
```

This component is unobservable as a first-order direction.

When squared inside the closure expression,

```math
(iz)^2=-z^2
```

it appears with a reversed sign.

Therefore, an unobservable compensating quantity in the projection-axis direction appears inside the intrinsic coordinate system as a sign-reversed compensating quantity through squared readout.

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

## Axiom 4: Anonymous Equal-Amplitude Composite Wave

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

## Axiom 5: Unreadability of Individual Phase

Classification: Axiom

```math
\theta_k\rightarrow\theta_k+\alpha
```

---

## Axiom 6: Phase-Difference Readout

Classification: Axiom

```math
\Delta_{ij}=\theta_i-\theta_j
```

Number of readable pair channels:

```math
\binom{N}{2}=\frac{N(N-1)}{2}
```

---

## Axiom 7: System Identification

Classification: Axiom

```math
(N,A_0,\{\Delta_{ij}\})
```

---

## Axiom 8: Composite Phase Projection

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

## Axiom 9: No Requirement of Inversion Symmetry

Classification: Axiom

Inversion symmetry is not required.

---

# 4. External Readout

## Axiom 10: External-Axis Response

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

## Axiom 11: External Readout Phase

Classification: Axiom candidate

```math
\Delta\phi_\mu=\phi_\mu^{(\Psi)}-\phi_\mu^{(M)}
```

```math
\mu=x,y,z,t
```

---

## Axiom 12: Position-Phase Readout

Classification: Axiom candidate

```math
x,y,z
```

---

## Axiom 13: Time-Phase Readout

Classification: Axiom candidate

```math
\phi_t\equiv\omega t_{\mathrm{read}}\pmod{2\pi}
```

---

## Axiom 14: Covering Phase

Classification: Axiom candidate

```text
2pi-periodic readout
4pi-periodic readout
winding-number-difference readout
```

---

# 5. Observational Phase Discreteness

## Axiom 15: Observational Phase Discreteness

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

---

# 9. Observation Selection and Curvature Projection

## Axiom 16: Unique Observation Selection by Complete Pairwise Relational Waves

Classification: observation-connection axiom

This axiom is an independent observation principle not derived from Axiom 0 or Axiom 1.

For a closed system with component count `N`, define the unordered complete pairwise relation set without proper names as

```math
\mathcal E_N
=\left\{\{i,j\}\mid 1\le i<j\le N\right\}
```

For each relation

```math
e=\{i,j\}\in\mathcal E_N
```

the corresponding

```math
X_e\in\mathbb C
```

is a physical relational wave that constitutes the system, not an observation-device channel.

The relational wave `X_e` is not a quantity derived from component waves; it is an independent physical state component. The correspondence map to component waves is not specified by this axiom.

Write the relational-wave state as

```math
X_{\mathcal E}
=\left(X_e\right)_{e\in\mathcal E_N}
```

A permutation of individual names acts as a permutation of relational waves and does not change the closure equation, the update rule, or the observation selection.

Define the quadratic form of the relational-wave set as

```math
q_{\mathcal E}(X)
=\sum_{e\in\mathcal E_N}X_e^2
```

The closure of the relational waves is written, using the compensation quantity of Axiom 3, as

```math
\sum_{e\in\mathcal E_N}X_e^2+(iR)^2=0
```

The transposed form

```math
q_{\mathcal E}(X)=R^2
```

may be used.

Let the number of relational waves be

```math
M=\left|\mathcal E_N\right|
```

and let a real-coefficient generator acting identically on all relational waves be

```math
K\in\mathbb R^{M\times M}
```

The relational wave `X_e` remains complex, and `K` acts as the same real-coefficient transformation on its real and imaginary parts.

Write the local linear update of relational waves as

```math
\dot X=KX
```

If this update preserves the quadratic form identically,

```math
\frac{d}{d\tau}q_{\mathcal E}(X)
=X^{\mathsf T}\left(K^{\mathsf T}+K\right)X
=0
```

Therefore, the closure-preserving generator of relational waves satisfies

```math
K^{\mathsf T}=-K
```

This antisymmetry is not introduced from any particular physical name or axis name.

A real antisymmetric generator decomposes, without additional proper names, privileged axes, or external references, into a direct sum of two-dimensional rotation planes and a null space.

Let the number of rotation planes be

```math
r=\frac12\operatorname{rank}K
```

and the null-space dimension be

```math
\dim\ker K=M-2r
```

The adoption conditions for intrinsically observable directions are as follows.

First, rotation planes become candidates for distinct directions only when they can be mutually distinguished by the relational quantity of independent rotation frequencies. Distinguishability is a necessary condition and does not by itself license adoption as spatial directions.

Second, when a normal direction is read from the plane spanned by two linearly independent relational directions, the normal may be adopted as a direction unique up to sign only when the normal-candidate space is one-dimensional. For a real `d`-dimensional spatial display, the dimension of the normal-candidate space is

```math
d-2
```

and this condition holds only for

```math
d=3
```

Third, the number of intrinsically observable directions is not counted from the relational-wave count `M`, the generator rank, or the null-space dimension. It is counted as the number of directions satisfying the uniqueness conditions above.

The minimal example is `M=3`: when

```math
\operatorname{rank}K=2,
\qquad
\dim\ker K=1
```

the two directions of one rotation plane and the one normal direction uniquely determined from that plane, three directions in total, may be adopted.

When `M=6` and three rotation planes are distinguished by independent rotation frequencies, three directions may be adopted.

When the number of distinguishable rotation planes exceeds the range of spatial display satisfying the second uniqueness condition, the excess rotation planes are not adopted as spatial directions and are retained as internal phase modes.

When multiple isomorphic rotation planes remain, or when a normal-candidate space or null space of dimension two or more remains and its internal directions cannot be distinguished by relational quantities alone, no single direction among them may be selected as an observable spatial direction. Such a residue is retained as an internal subspace without a unique basis.

The primitive readout of each relational wave is constructed from the phase difference of Axiom 6,

```math
\Delta_{ij}=\theta_i-\theta_j
```

and the bilinear relation obtained from the quadratic form of Axiom 1,

```math
q(X)=\sum_n X_n^2
```

```math
B(U,V)
=\frac{1}{2}\left(q(U+V)-q(U)-q(V)\right)
```

Readout for three or more bodies is constructed from the set of complete pairwise relational waves and the invariant decomposition of their closure-preserving generator.

No independent higher-order readout rule may be added on the grounds of physical names, object names, or dimension names.

The specific coupling strengths, phase-difference dependence, and update period of the relational-wave generator are not specified by this axiom.

They are defined as independent relational update rules consistent with Axiom 0, Axiom 0.5, Axiom 1, and this axiom.

---

## Working Axiom I: Phase-Difference Sine Relational Update Rule

Classification: working hypothesis

As a concrete form of the relational-wave generator of Axiom 16, define the endpoint-sharing adjacency

```math
A_{ef}
=
\begin{cases}
1,& e\ne f\ \text{and}\ e\cap f\ne\varnothing,\\
0,& \text{otherwise}
\end{cases}
```

and, from the initial phases `theta_e = arg X_e`, set

```math
\widetilde K_{ef}
=
A_{ef}\sin(\theta_f-\theta_e)
```

It satisfies

```math
\widetilde K^{\mathsf T}=-\widetilde K
```

and thus the antisymmetry condition of Axiom 16.

This coupling rule is not derived from Axiom 0, Axiom 0.5, Axiom 1, or Axiom 16. It is a working hypothesis placed as one concrete closure-preserving relational update rule. Note that this generator, as a function of phase differences only, is zeroth-order homogeneous and satisfies Axiom 0.5 without adjustment (Paper 5, Theorem 3. Concept DOI: 10.5281/zenodo.21486233).

---

## Consequence 16.1: Vertex Decomposition and Linear Rank Bound

Classification: derived consequence (under Working Axiom I)

Let `B` be the unsigned incidence matrix with individuals as rows and relations as columns, and from its `k`-th row build

```math
c_k=(B_{ke}\cos\theta_e)_{e},
\qquad
s_k=(B_{ke}\sin\theta_e)_{e}
```

Then

```math
\widetilde K=\sum_{k=1}^{N}\left(c_k s_k^{\mathsf T}-s_k c_k^{\mathsf T}\right)
```

holds exactly.

Each term is an antisymmetric matrix of rank at most 2, so for all `N` and all initial phases,

```math
\operatorname{rank}\widetilde K
\le
2\min\!\left(N,\left\lfloor\frac{M}{2}\right\rfloor\right)
```

The number of relational waves

```math
M=\binom{N}{2}
```

grows quadratically, while the number of rotation modes

```math
r=\frac12\operatorname{rank}\widetilde K
```

grows at most linearly.

Proof: Paper 3 of the Dimensional Generation Structure series. Concept DOI: 10.5281/zenodo.21465898

---

## Consequence 16.2: Generic-Position Rank Equality

Classification: derived consequence (computer-assisted proof, `3 <= N <= 12`)

For each `N` with `3 <= N <= 12`, for all initial phases except a set of Lebesgue measure zero,

```math
\operatorname{rank}\widetilde K
=
2\min\!\left(N,\left\lfloor\frac{M}{2}\right\rfloor\right)
```

The proof combines exact rational witness configurations via the tan half-angle parametrization with the fact that the zero set of a nontrivial real-analytic function has Lebesgue measure zero. Paper 3. Concept DOI: 10.5281/zenodo.21465898

The generic null-space dimension is then

```math
\dim\ker\widetilde K
=
M-2\min\!\left(N,\left\lfloor\frac{M}{2}\right\rfloor\right)
```

which for `N >= 5` equals

```math
\frac{N(N-5)}{2}
```

Equality for all `N` is a new hypothesis and is not included in this consequence.

---

## Consequence 16.3: Normal Uniqueness

Classification: derived consequence

When a normal direction is reconstructed from the plane spanned by two linearly independent relational directions, without additional background axes or selection rules, the normal-candidate space of a real `d`-dimensional display has dimension

```math
d-2
```

and the normal is unique up to sign only for

```math
d=3
```

For `d >= 4`, the orthogonal group

```math
O(d-2)
```

acting on the normal-candidate space while fixing the plane leaves a continuous selection freedom, and under anonymity (Axiom 0) no single direction can be selected. The projector and squared quantities are unique, but the individual linear directions inside are not.

This consequence is the basis of the second adoption condition of Axiom 16. Proof: Paper 3. Concept DOI: 10.5281/zenodo.21465898
---

## Consequence 16.4: Plane Decomposition of Fixed Generators and Kinematic Isomorphism

Classification: derived consequence (under a fixed generator)

A fixed closure-preserving generator decomposes into a direct sum of two-dimensional rotation planes and a null space. A Cayley-type orthogonal update acts on each rotation plane as a rotation by the constant angle

```math
\theta_j=2\arctan(\gamma\sigma_j)
```

and on the null space as the identity. The coordinate of each plane is therefore a constant-amplitude single-frequency rotation, kinematically isomorphic to the rotational kinematics of a single two-body relation.

The orientation of each plane is fixed by the generator itself through

```math
Kp_j=\sigma_j q_j
```

The residual gauge preserving the standard form is `SO(2)`; on planes of nonzero amplitude the signed phase progression, the frequency, and the amplitude are readable, while the absolute phase is not. This means that the unreadability of Axiom 5 is reproduced also at the level of the derived plane coordinates. If reflections are allowed (`O(2)`), the invariants are the absolute values of progression and frequency.

The null-space component is static and carries no phase progression, so progression-based readouts do not hold from the null space (the dynamical realization of Consequence 16.3).

Frequency nondegeneracy is an independent general-position assumption, not derivable from the generic rank equality (Consequence 16.2).

Proof: Paper 4 of the Dimensional Generation Structure series. Concept DOI: 10.5281/zenodo.21468959

---

## Consequence 16.5: Double Conservation Register

Classification: derived consequence (under a fixed generator)

Both the energy-like quantity and the squared closure decompose plane-wise and null-space-wise, with each term individually conserved under the orthogonal update of a fixed generator.

```math
H=\sum_{j}H_j+H_{\ker},
\qquad
H_j=\lvert p_j^\dagger X\rvert^2+\lvert q_j^\dagger X\rvert^2
```

```math
X^{\mathsf T}X=\sum_{j}Q_j+Q_{\ker}=R^2,
\qquad
Q_j=(p_j^{\mathsf T}X)^2+(q_j^{\mathsf T}X)^2
```

The fixed-generator system thus carries a double conservation register `(H_j, Q_j)` on each plane and on the null space.

Proof: Paper 4. Concept DOI: 10.5281/zenodo.21468959

---

## Consequence 16.6: No-Go for Spontaneous Splitting

Classification: derived consequence (under a fixed generator)

Under a fixed generator, neither the closure quantities `Q_j` nor the energy-like quantities `H_j` move between planes, or between planes and the null space.

Spontaneous splitting of the state accompanied by redistribution among registers therefore cannot occur under fixed-generator dynamics. Describing splitting requires a relational update rule beyond the fixed generator: sequentially reconstructed generators, nonlinear couplings, or inter-plane exchange.

This consequence provides the starting boundary for any additional definition of splitting dynamics consistent with Axiom 0, Axiom 0.5, Axiom 1, and Axiom 16.

Proof: Paper 4. Concept DOI: 10.5281/zenodo.21468959

---

## Axiom 17: Curvature Reaction Projection Centered on the Future Phase Position

Classification: dynamical projection axiom

This axiom is an independent dynamical and projection principle not derived from Axiom 0 or Axiom 1.

Assign an order index

```math
\tau\in\mathbb{N}
```

to the state sequence of a closed system.

In the state update

```math
X(\tau)\longrightarrow X(\tau+1)
```

the phase position on the `tau+1` side is defined as the future phase position of that state sequence.

When the phase progression on a rotation plane uniquely constructed without additional names by Axiom 16 has a curvature-radius readout `rho_c` and an angular-velocity readout `omega`, and can be expressed as a motion whose relational virtual rotation center is the future phase position, the center-directed compensation that closes this virtual rotation — that is, the reaction to the centrifugal force — is mapped and read out as an observable acceleration in the circumferential direction of progression.

The center-directed compensation is the internal representation, and the circumferential acceleration is the external readout. This axiom defines the projection connecting the two as a rule.

The magnitude of the acceleration is

```math
|\alpha_{\mathrm{read}}|
=\rho_c|\omega|^2
```

Let `t_hat_tau` be the unit vector of the circumferential direction of progression pointing from the current phase position toward the future phase position. The readout acceleration is

```math
\boldsymbol{\alpha}_{\mathrm{read}}
=\rho_c|\omega|^2\,\hat{\boldsymbol{t}}_\tau
```

Here `rho_c`, `omega`, and `t_hat_tau` are relational readout quantities without physical names in the first-principles layer. The circumferential direction is not an absolute axis of a background space; it is determined each time from the relation among the current phase position, the future phase position, and the closure radius.

This axiom may not be modified on the grounds of gravity, Coulomb force, mass, charge, or any other physical name.

The same projection rule applies to every closure mode satisfying the same conditions.

If the phase branch opposite to the future side, a reversal of the acceleration direction, or sign retention or sign elimination is adopted, it must not be selected from physical names; it must be defined in advance as an independent phase relation, branch selection, symmetry breaking, or readout operation.

This axiom does not identify the readout acceleration with standard gravity, the standard Coulomb force, or any known external force.

Physical names are given only after the acceleration map, distance exponent, sign, closure preservation, and gauge stability have been confirmed.
