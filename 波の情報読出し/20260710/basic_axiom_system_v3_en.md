# Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v3 Interpretation Note

**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Status:** Interpretation and working note for Basic Axiom System v3

This document cites the definition paper "Anonymous Equal-Amplitude Composite-Wave Model: Basic Axiom System v3" and records the present interpretation, development, and working explanations of that axiom system.

The canonical axiom definition is placed in the following definition paper:

```text
anonymous_equal_amplitude_composite_wave_model_basic_axiom_system_v3_pure_definition_en.md
```

The interpretations, explanations, working hypotheses, and readout examples in this document may be extended or revised.

However, a change in interpretation is not a rejection of the axioms.

If an error, omission, or required change is found in the axioms themselves, the canonical axiom system must be explicitly revised rather than being reinterpreted inside this note.

Thus this document is not the axiom system itself. It is the current interpretation layer for Basic Axiom System v3.

In v3, anonymity is extended beyond basic components, physical quantities, and observables to theory names, formulation names, and readout-logic names. In other words, the names gravity, Coulomb force, mass, charge, energy, two-body, three-body, and similar labels must not be used as reasons to change equations, readout rules, or projection rules at the first-principle layer. Physical names are labels assigned after readout results have been obtained from the same axioms.

---

# 1. First Principles

## Axiom 0: Anonymity

Classification: Axiom

Basic components are not assigned individual names, privileged axes, privileged signs, or privileged types. The index `n` is a convenient label and is not an intrinsic name of the component.

Forbidden:

```text
making only one specific component the negative-sign axis
making only one specific component the time axis
making only one specific component real-valued
making only one specific component complex-valued
giving an intrinsic name to only one specific component
```

Therefore,

```math
x_1^2+x_2^2-x_3^2=0
```

is not adopted at the first-principle layer, because `x_3` alone becomes the negative-sign axis.

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

At the first-principle layer, only the same closure condition,

```math
\sum_n x_n^2=0
```

and readout operations, symmetry breaking, gauge stability, sign preservation, sign erasure, phase-branch selection, and closure target sets defined in advance for that same closure condition are allowed.

The following order is forbidden:

```text
call it gravity
therefore adopt an R-type readout
therefore erase signs
```

Likewise, the following order is forbidden:

```text
call it Coulomb force
therefore adopt a Q-type readout
therefore preserve signs
```

The correct order is:

```text
place the same all-positive-sign zero closure
define the readout operation first
define symmetry breaking first
define sign preservation or sign erasure first
confirm gauge stability
assign working names after readout
```

Here `R`, `Q`, and `E` are not first-principle names. They are working names assigned to readout results obtained from the same closure condition.

In particular, an already defined `R_read` must not be changed into another definition because one wishes to read gravity.

What is allowed is to apply the same `R_read` to a single closed wave, a pair composite closed wave, or a whole-system closed wave, and to read differences in the target set as:

```text
R_A
R_{AB}
R_{ABC}
```

Likewise, the distinction between a two-body and a three-body system is treated as a difference in closure target set, not as a physical-name-based change in equation.

For example, the following are allowed because the target sets differ:

```math
Q(A)=0
```

```math
Q(A+B)=0
```

```math
Q(A+B+C)=0
```

However, each `Q` must have the same form:

```math
Q(P)=\sum_n p_n^2
```

Different readouts may be obtained from different target sets. Changing the definition of `Q` itself to obtain a desired physical name is not allowed.

Therefore, formulation anonymity extends Axiom 0, component anonymity, to the level of theoretical construction.

This axiom requires:

```text
closure equation before physical name
readout operation before classification name
symmetry breaking before theory name
gauge stability before interpretation name
```

Physical names, classification names, theory names, and interpretation names are not first principles. They are names assigned after readout.

## Axiom 1: All-Positive-Sign Zero Closure

Classification: Axiom

```math
\sum_{n=1}^{N}x_n^2=0
```

Expanded form:

```math
x_1^2+x_2^2+\cdots+x_N^2=0
```

All coefficient signs are the same positive sign. This is not the conjugate norm.

```math
\sum_n |x_n|^2=0
```

is not the intended closure. The closure uses `x_n^2`, not `x_n\bar{x}_n`.

Axiom 1 is the closure condition for a state in which the left-hand side of the closure sum has become stationary or stable.

It does not deny the existence of a quasi-stationary or metastable state between an external influence, interaction, observation, or readout re-embedding that changes the closure sum and the later return to stationarity or stability.

Therefore, in a quasi-stationary state, the following temporary nonzero closure residual is allowed:

```math
\sum_n x_n^2\neq0
```

However, in a state treated as stationary closure or stable closure,

```math
\sum_n x_n^2=0
```

holds.

Thus, later papers may treat closure recovery processes, metastable transitions, or closure residuals immediately after observation without rejecting Axiom 1.

Axiom 1 itself defines the closure condition after stationarization, not the quasi-stationary process.

## Axiom 2: Nontrivial Existence

Classification: Axiom

```math
\exists n,\quad x_n\neq0
```

The closed system is not restricted to the all-zero solution.

---

# 2. Consequences from the First Principles

## Consequence 1: Insufficiency of Real Positive-Definite Components

Classification: Derived consequence

```math
x_n\in\mathbb{R}
\Rightarrow
x_n^2\ge0
```

Therefore,

```math
\sum_n x_n^2=0
\Rightarrow
x_n=0\quad\forall n
```

To satisfy Axioms 0, 1, and 2 simultaneously, real positive-definite components alone are insufficient.

## Consequence 2: Internal Generation of Sign Reversal

Classification: Derived consequence

A negative sign as an external coefficient violates Axiom 0. The sign reversal required for closure must arise from the square of a component.

Minimal two-component solution:

```math
x_1=A,\qquad x_2=iA
```

```math
x_1^2+x_2^2=A^2+(iA)^2=A^2-A^2=0
```

```math
i^2=-1
```

The negative sign is not a privileged coefficient. It is a squared phase.

## Consequence 3: Necessity of Phase Algebra

Classification: Derived consequence

Nontrivial closure requires a common phase algebra that internally generates sign reversal by squaring.

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

Under the equal-amplitude condition `r_n=r`,

```math
\sum_{n=1}^{N}e^{i2\phi_n}=0
```

The complex phase type is the minimal implementation of this requirement.

## Consequence 4: Type Anonymity

Classification: Derived consequence

Making only one component complex-valued violates Axiom 0. Therefore,

```math
x_n\in\mathbb{C}\quad\forall n
```

or all components must be given the same phase-algebra type.

## Consequence 5: Ninety-Degree Phase Difference

Classification: Derived consequence

```math
A^2+(iA)^2=0
```

A ninety-degree phase difference is the minimal two-component structure of all-positive-sign zero closure.

---

# 3. Anonymous Equal-Amplitude Composite Waves

## Axiom 3: Anonymous Equal-Amplitude Composite Wave

Classification: Axiom

```math
\Psi_N=\frac{A_0}{N}\sum_{k=1}^{N}e^{i\theta_k}
```

```math
|\psi_k|=\frac{A_0}{N}
```

Components are not individually identified by amplitude. Since anonymous basic components are not distinguished by amplitude differences, all components are taken to have equal amplitude in the basic implementation of the composite wave.

Normalization condition:

When all phases are aligned,

```math
\max|\Psi_N|=A_0
```

Thus the normalization is

```math
\Psi_N=\frac{A_0}{N}\sum_{k=1}^{N}e^{i\theta_k}
```

What is preserved is the maximum composite amplitude readable from outside.

## Axiom 4: Unreadability of Individual Absolute Phase

Classification: Axiom

```math
\theta_k\rightarrow\theta_k+\alpha
```

does not change the internal structure. Absolute phase is not an intrinsic observable.

## Axiom 5: Phase-Difference Readout

Classification: Axiom

```math
\Delta_{ij}=\theta_i-\theta_j
```

The number of readable pair channels is

```math
\binom{N}{2}=\frac{N(N-1)}{2}
```

This is not the number of degrees of freedom, but the number of phase-difference relations.

## Axiom 6: System Identification

Classification: Axiom

Two systems with the same

```math
(N,A_0,\{\Delta_{ij}\})
```

are identified as the same system.

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

`Θ_N` is an external projected phase defined when `Ψ_N != 0`. Even when `Ψ_N=0`, the internal phase-difference network exists.

## Axiom 8: Reflection Symmetry Is Not Required

Classification: Axiom

A composite wave is not required to be reflection symmetric.

## Consequence 6: Reflection Symmetry Is a Special Condition

Classification: Derived consequence from Axiom 8

Reflection symmetry is a special condition.

```math
\Psi(-x)=\Psi(x)
```

In general, reflection-asymmetric configurations are allowed.

```math
\Psi(-x)\neq\Psi(x)
```

A reflection-symmetric wave is a special configuration corresponding to equal weighting of forward and reverse phase progression.

---

# 4. External Readout

## Axiom 9: External-Axis Response

Classification: Candidate axiom

Let the external axis be

```math
M_\alpha=e^{i\alpha}
```

Then

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

Spin-like information is defined by correlation readout with an external axis.

## Axiom 10: External Readout Phase

Classification: Candidate axiom

```math
\Delta\phi_\mu=\phi_\mu^{(\Psi)}-\phi_\mu^{(M)}
```

```math
\mu=x,y,z,t
```

Position-like, time-like, and spin-like readouts are relative phases with respect to an external reference.

## Axiom 11: Position-Phase Readout

Classification: Candidate axiom

```math
x,y,z
```

are not intrinsic labels of the composite wave. They are reconstructed from relative phase with an external spatial readout wave.

The position notation used here is a convenient reinterpretation and does not prove connection to physical position quantities.

## Axiom 12: Time-Phase Readout

Classification: Candidate axiom

```math
\phi_t\equiv\omega t_{\mathrm{read}}\pmod{2\pi}
```

Time is a relative-phase readout with respect to an external frequency reference. The parameter `t_read` here is not an a priori absolute time, but a convenient parameter reconstructed from time-phase readout.

## Axiom 13: Covering Phase

Classification: Candidate axiom

When the readout system preserves winding number, it distinguishes:

```text
2π-period readout
4π-period readout
winding-difference readout
```

---

# 5. Observed Phase Discreteness

## Axiom 14: Observed Phase Discreteness

Classification: Candidate axiom

A phase fixed as an observed value is discretized by a loop condition or winding number.

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

Candidate phase-area quantity:

```math
Q_\phi=\sum_{i<j}\sin^2\frac{\Delta_{ij}}{2}
```

## Readout Example: Bounded Closure

Classification: Implementation example

In a closed readout system,

```math
kL=2\pi m
```

```math
k=\frac{2\pi m}{L}
```

The observable phase change is integerized as a winding number. The first-principle closure form is Axiom 1.

---

# 6. Waveform, Localization, and Observation

## Working Axiom A: State Waveform

Classification: Working hypothesis

```math
A(\theta)=\sum_n A_n\cos(n\theta+\phi_n)
```

```math
Z(\theta)=\sum_n r_n e^{i(n\theta+\phi_n)}
```

A state is a hierarchical waveform.

## Definition A1: Normalized Equal-Amplitude Odd-Harmonic Localized Wave

Classification: Definition

Let the highest odd-harmonic order be

```math
N_h=2K-1
```

Define the normalized equal-amplitude odd-harmonic localized wave by

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

In complex notation,

```math
Z_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}e^{i(2m+1)u}
```

A localized wave with representative amplitude `A` is

```math
\Psi_{N_h}(u)=A S_{N_h}(u)
```

or

```math
\Psi_{N_h}(u)=A Z_{N_h}(u)
```

Each odd-harmonic component has amplitude

```math
\frac{A}{K}
```

and is normalized by the inverse of the number `K` of components. This normalization keeps the representative amplitude at in-phase alignment equal to `A`, even as `N_h` increases.

## Theorem A1: Preservation of Representative Amplitude

Classification: Derived theorem

Consider

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

At the central phase `u=0`,

```math
\cos((2m+1)0)=1
```

and hence

```math
S_{N_h}(0)=\frac{1}{K}\sum_{m=0}^{K-1}1=1
```

For

```math
\Psi_{N_h}(u)=A S_{N_h}(u)
```

we obtain

```math
\Psi_{N_h}(0)=A
```

Therefore, even as the highest odd-harmonic order `N_h` increases, the normalization factor `1/K` preserves the composite peak representative amplitude as `A`.

## Theorem A2: Highest Odd-Harmonic Order and Localization Width

Classification: Derived theorem

Let the base wavelength of the parent system be

```math
\lambda_0
```

Here `N_h` is not the number of odd-harmonic components, but the highest odd-harmonic multiplier relative to the parent base wave. The number of odd-harmonic components is

```math
K=\frac{N_h+1}{2}
```

The wavelength corresponding to the highest harmonic is conceptually

```math
\lambda_{N_h}=\frac{\lambda_0}{N_h}
```

Thus the localization width obtained by normalized equal-amplitude odd-harmonic synthesis is estimated as

```math
\Delta x_{N_h}\sim\frac{\lambda_0}{N_h}
```

or, using a half-wavelength criterion,

```math
\Delta x_{N_h}\sim\frac{\lambda_0}{2N_h}
```

The exact coefficient depends on whether one uses main-lobe width, half-maximum width, or first-zero width. This document uses

```math
\Delta x_{N_h}\propto\frac{\lambda_0}{N_h}
```

as a design estimate.

Therefore, the highest odd-harmonic order required to obtain a local width `Δx` from a base wavelength `λ_0` is approximately

```math
N_{h,\rm req}\sim\frac{\lambda_0}{\Delta x}
```

In implementation, an odd order at least as large as `N_{h,\rm req}` is used.

For the temporal direction, if the parent base period is `T_0` and the temporal localization width is `Δt`, then

```math
N_{h,\rm req}^{(t)}\sim\frac{T_0}{\Delta t}
```

## Consequence A1: Common Construction of Spatial and Temporal Localization

Classification: Derived consequence

The normalized equal-amplitude odd-harmonic localized wave can be applied not only to spatial phase but also to temporal phase.

Let spatial phase be `χ` and temporal phase be `τ`. A local wave localized in both spatial and temporal directions may be written as

```math
\Psi(\chi,\tau)=A S_{N_{h,\chi}}(\chi-\chi_0)S_{N_{h,\tau}}(\tau-\tau_0)e^{i\phi}
```

where `χ_0` is position phase, `τ_0` is time phase, `N_{h,χ}` is the highest odd-harmonic order in the spatial direction, and `N_{h,τ}` is the highest odd-harmonic order in the temporal direction.

This form represents spatial and temporal localization by the odd-harmonic structure of the waveform itself, without introducing an externally invisible temporal window function.

## Note A1: Quantity Preserved

Classification: Note

In this axiom system, what is preserved by the normalized equal-amplitude odd-harmonic localized wave is the representative amplitude `A` and the square register governed by Axiom 1, all-positive-sign zero closure.

Normalizing each odd-harmonic component to amplitude

```math
\frac{A}{K}
```

keeps the in-phase composite representative amplitude at

```math
A
```

Thus, increasing `N_h` does not make the in-phase representative amplitude diverge.

This document places the basis of conservation in

```math
\sum_n x_n^2=0
```

It is not a conjugate-norm conservation law

```math
\sum_n |x_n|^2
```

and it is not a conservation law with externally imposed weights on higher harmonic orders.

Localization by higher harmonics is treated as refinement of phase structure while preserving the representative amplitude `A`. Increasing the highest odd-harmonic order is therefore not an operation that externally increases the conserved quantity. It is an operation that refines phase structure with the same representative amplitude inside the closure condition of Axiom 1.

## Working Axiom B: Localization

Classification: Working hypothesis

Localization is represented as phase alignment of higher odd-harmonic components with respect to a low-order base wave. This document does not introduce localization as an external window function.

Localization is constructed by

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

or by its complex representation

```math
Z_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}e^{i(2m+1)u}
```

The larger `N_h` is, the narrower the localization width becomes.

If the base wavelength of the parent system is `λ_0` and the required localization width is `Δx`, the required highest odd-harmonic order is estimated as

```math
N_{h,\rm req}\sim\frac{\lambda_0}{\Delta x}
```

Similarly, for a temporal base period `T_0` and temporal localization width `Δt`,

```math
N_{h,\rm req}^{(t)}\sim\frac{T_0}{\Delta t}
```

Thus, when the parent base wavelength or period is very large, extremely large highest odd-harmonic order is required to construct local spatial or temporal localization.

## Working Axiom C: Coupling Entrance

Classification: Working hypothesis

The entrance to interaction is alignment of low-order base phases.

## Working Axiom D: Observation Map

Classification: Working hypothesis

```math
s_\alpha=\mathrm{Loc}\left[\int A(\theta)M_\alpha(\theta)d\theta\right]
```

`Loc` is a map that produces a stable localized output from a correlation amplitude. Observation here is an internal definition of this axiom system and does not claim connection to standard-theory observation.

## Working Axiom E: Correlation-Strength Squared

Classification: Working hypothesis

```math
P(\alpha)=\frac{|C_\alpha|^2}{\sum_\beta |C_\beta|^2}
```

This formula is a normalization of correlation readout strength within this axiom system.

## Working Axiom F: Low-Order Base Synchronization Sharing

Classification: Working hypothesis

```math
\theta_A-\theta_B=\Delta
```

Low-order base synchronization sharing is a state in which two waveforms share a low-order base phase difference.

## Working Axiom G: Loss of Synchronous Demodulability

Classification: Working hypothesis

```math
V(t)=\left|\sum_n w_n e^{i\epsilon_n(t)}\right|
```

Phase-fluctuation average:

```math
V(t)=\left|\sum_n w_n e^{-\frac12\sigma_n^2(t)}\right|
```

Loss of synchronous demodulability is a state in which phase fluctuations prevent a correlation readout including higher orders from producing a stable output.

## Working Axiom H: Hierarchy

Classification: Working hypothesis

A closed wave system has the hierarchy

```text
low-order base mode
+ higher-order localized mode
```

---

# 7. Readout, Memory, and Complex Representation

## I/Q Readout

Classification: Theoretical insight

```math
z=a+ib
```

```math
a=r\cos\theta,\qquad b=r\sin\theta
```

This is read as a pair of ninety-degree phase-shifted wave sources of the same physical quantity.

```math
u=r\cos\theta+r\sin\theta
=\sqrt2 r\cos\left(\theta-\frac{\pi}{4}\right)
```

Complex representation is the minimal memory of a present value and a ninety-degree delayed component.

## Harmonic Information Readout

Classification: Theoretical insight

```math
Z(\theta)=\sum_n r_n e^{i(n\theta+\phi_n)}
```

Readout map:

```math
P(Z)=\operatorname{Re}(Z)+\operatorname{Im}(Z)
```

```math
P(Z)=\sum_n \sqrt2 r_n
\cos\left(n\theta+\phi_n-\frac{\pi}{4}\right)
```

The order is preserved, and the readout phase is shifted by forty-five degrees.

## Conservation Side and Readout Side

Classification: Theoretical insight

```math
z\bar z=a^2+b^2
```

is a conjugate norm.

```math
z^2=a^2-b^2+2iab
```

preserves interference phase.

Axiom 1 uses the `z^2` type, not `z\bar z`.

## Square Register

Classification: Working hypothesis

Internal register:

```math
\mathcal{A}_0=\sum_k A_k^2
```

External standard intensity:

```math
a_0=\sum_k |A_k|^2=\sum_k A_k\bar A_k
```

`\mathcal{A}_0` preserves interference phase. `a_0` removes interference phase through conjugation.

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

Pair composition:

```math
\psi_i+\psi_j
=2A\cos\frac{\Delta_{ij}}{2}e^{i(\theta_i+\theta_j)/2}
```

Pair interference intensity:

```math
I_{ij}=4A^2\cos^2\frac{\Delta_{ij}}{2}
```

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

## Fermion Degree

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

Relation to overlap degree:

```math
F_N=\frac{N}{2(N-1)}(1-\rho_N^2)
```

Therefore,

```math
\rho_N=0\Rightarrow F_N=\frac{N}{2(N-1)}
```

`B_N` is the pair in-phase degree.

Classification:

```text
F_N=0: boson-like
N=2, F_N=1: purely fermion-like
0<F_N<1: elmion-like
```

## Boson-Like Composite Wave

Classification: Definition

```math
\theta_1=\theta_2=\cdots=\theta_N
```

```math
F_N=0,\qquad B_N=1,\qquad \rho_N=1,\qquad R_N=A_0
```

## Purely Fermion-Like Composite Wave

Classification: Definition

```math
N=2,\qquad \Delta_{12}=\pi
```

```math
e^{i\theta_1}+e^{i\theta_2}=0
```

```math
F_2=1,\qquad \rho_2=0,\qquad R_2=0
```

For three or more components, the condition that all pairs be opposite in phase cannot hold.

```math
\Delta_{12}=\pi,\quad \Delta_{23}=\pi
\Rightarrow
\Delta_{13}=2\pi\equiv0
```

## Elmion-Like Composite Wave

Classification: Definition

The main condition is

```math
0<F_N<1
```

A subset with nonzero external composite amplitude is

```math
0<\rho_N<1
```

```math
\rho_N=0\not\Rightarrow F_N=1
```

Zero external composite amplitude and pure two-body opposite-phase character are not the same.
