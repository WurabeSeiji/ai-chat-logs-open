# Thermodynamic Readout from a Self-Consistent Complex Relational System
## without Presupposing Spacetime

**--- Local Zero-Closure Subsystems, State Counting, Entropy, Energy, and Temperature ---**

**Author:** Noriaki Kihara\
**ORCID:** 0009-0004-6753-4020\
**Document Type:** Research Note / Hypothesis and Exploratory Study\
**Date:** 2026-09-02\
**DOI:** \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\
**Concept DOI:**
\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

------------------------------------------------------------------------

## Abstract

This paper investigates whether thermodynamic quantities can be constructed from complex relations and their self-consistency, without presupposing space, time, particle, mass, energy, or temperature as fundamental concepts.

The fundamental system consists of a finite set of complex states or complex relations satisfying a self-consistency condition: a **zero closure**

$$
\sum_{n=1}^{N} x_n^2 = 0
$$

We consider a finite phase period

$$
U^N=I
$$

The global zero closure can be decomposed into multiple local zero-closure subsystems:

$$
0=0+0+\cdots+0
$$

In this framework, distinct self-consistent decompositions of local closures that realize the same global zero closure can be counted as distinct states. If the number of permitted states is $\Omega$, then in units where $k_B=1$, we obtain a natural entropy-type readout:

$$
S=\ln\Omega
$$

Moreover, if a local closure subsystem satisfies $U^{N_k}=I$, its fundamental phase step is $2\pi/N_k$. Constructing an additive energy-type measure from this internal phase scale, we can define a temperature-type readout:

$$
\frac{1}{T_{\mathrm{read}}}
=
\frac{\Delta S}{\Delta E_{\mathrm{read}}}
$$

The central claim of this paper is **not** to derive existing thermodynamics completely. More precisely, it demonstrates that **in a self-consistent complex relational system without spacetime or particles as foundations, one can construct observational readout structures isomorphic to state counting, entropy, energy, and temperature**. Moreover, these readout values depend on which local zero-closure subsystem is selected as the observation target. This readout framework belongs to the same **observational readout hierarchy** examined in prior work on curvature readouts (minimal application to gravitational waves) and exhibits two distinct, non-unified observational mappings.

Finally, we compare the fact that for a cosmological discreteness parameter $N\sim10^{60}$, the two-body relationship count reaches $M\sim10^{120}$, with the observed baryon-to-photon ratio in the observable universe, approximately $10^{-9}$. This comparison is not a claim of numerical coincidence, but rather a **decisive consistency check** to confirm that the direction of this model—wherein material local states become sparse within the vast relational space—does not immediately contradict the observable universe at leading order.

**Keywords:** complex relational system, zero closure, self-consistency, state counting, entropy, emergent thermodynamics, observer dependence, finite phase periodicity, pre-spatiotemporal physics

------------------------------------------------------------------------

## 1. Introduction

In statistical mechanics, entropy is linked to the number of microscopic states. The relationship famously established by Boltzmann and later concisely formulated by Planck is

$$
S=k_B\ln\Omega
$$

\[1\]. In the microcanonical description, temperature is expressed as the energy derivative of entropy:

$$
\frac{1}{T}
=
\left(\frac{\partial S}{\partial E}\right)
$$

Conventionally, this construction presupposes that "system," "energy," "volume," and "particle number" are already defined. The present paper reverses this order: assuming only self-consistent complex relations as fundamental, without prior placement of space, time, particle, or energy, we ask whether state counting and thermodynamic readouts can nonetheless be constructed.

In the author's prior work, the minimal structure of reading the two-directional intersection of complex spiral waves as curvature oscillations was examined. It was shown that from the same local complex wave, the radius of curvature is not uniquely determined; rather, the readout scale of curvature can vary depending on which physical scale is adopted as the observation target \[2\]. In the present paper, this idea of "separating the fundamental state from the observational readout" is extended to state counting, entropy, energy, and temperature.

------------------------------------------------------------------------

## 2. Fundamental System

### 2.1 Complex States

Consider $N$ complex quantities:

$$
x_n=a_n+i b_n
\qquad
(n=1,\ldots,N)
$$

Formally, $N$ complex quantities carry $2N$ real components. However, in this paper they are not defined as spatial coordinates. Each $x_n$ represents a fundamental state or complex relation.

### 2.2 Zero Closure

As a fundamental structure arising from self-consistency, we consider

$$
\sum_{n=1}^{N}x_n^2=0
\tag{1}
$$

Since $x_n$ is complex, equation (1) does not imply each $x_n$ vanishes. Separating into real and imaginary parts:

$$
\sum_{n=1}^{N}(a_n^2-b_n^2)=0,
$$

$$
\sum_{n=1}^{N}a_n b_n=0
$$

Non-trivial complex states can globally satisfy closure.

Moreover, equation (1) is invariant under uniform scaling:

$$
x_n\rightarrow\lambda x_n
$$

to yield

$$
\sum_{n=1}^{N}(\lambda x_n)^2
=
\lambda^2\sum_{n=1}^{N}x_n^2
=
0
$$

Thus, zero closure itself does not require an absolute scale.

### 2.3 Finite Phase Periodicity

We further impose a finite periodicity condition:

$$
U^N=I
\tag{2}
$$

If the minimal phase rotation is read as

$$
U=\exp\left(\frac{2\pi i}{N}\right)
$$

then the permitted fundamental phases are

$$
\theta_m=\frac{2\pi m}{N}
$$

and the normalized minimal phase step is

$$
\frac{\Delta\theta}{2\pi}
=
\frac{1}{N}
\tag{3}
$$

In this paper, equations (1) and (2) do not presuppose space or time.

------------------------------------------------------------------------

## 3. Decomposition into Local Zero Closures

Suppose the global system satisfies zero closure:

$$
\sum_{n=1}^{N}x_n^2=0.
$$

Partition the index set into mutually disjoint subsets $I_k$, such that for each

$$
\sum_{n\in I_k}x_n^2=0
\qquad
(k=1,\ldots,K)
\tag{4}
$$

holds. Then the global system decomposes as

$$
0
=
\underbrace{0}_{I_1}
+
\underbrace{0}_{I_2}
+\cdots+
\underbrace{0}_{I_K}
\tag{5}
$$

The size of each subsystem is

$$
N_k=|I_k|
$$

so that

$$
N=\sum_{k=1}^{K}N_k.
\tag{6}
$$

What is crucial is not only that the whole is zero. For a given global zero closure, multiple decompositions into local zero closures may exist.

In this paper, we count **different self-consistent realizations of zero closure decomposition** as distinct states.

------------------------------------------------------------------------

## 4. State Counting

For fixed $N$, let the set of permitted states satisfying the self-consistency, zero-closure, and finite-phase conditions simultaneously be

$$
\mathcal{S}_N
$$

Define its cardinality as

$$
\Omega_N
=
|\mathcal{S}_N|
\tag{7}
$$

Note that $\Omega_N$ is not the combinatorial count of placing particles in a pre-given space. Space and particles have not yet been introduced as fundamental variables.

$\Omega_N$ counts:

> The number of permitted self-consistent relational states or local closure decompositions that realize the same global zero closure.

If continuous complex amplitudes are freely allowed, the state set can become continuous; hence, a finite $\Omega_N$ requires discretization through finite phase periodicity, self-consistency, amplitude conditions, and equivalence relations. Establishing the precise counting rule is **an immediate core task of this theory**.

------------------------------------------------------------------------

## 5. Entropy-Type Readout

Using formal correspondence with Boltzmann's statistical entropy and normalizing to $k_B=1$, we define

$$
S_N=\ln\Omega_N
\tag{8}
$$

For two independent closed systems $A$ and $B$, if

$$
\Omega_{A+B}
=
\Omega_A\Omega_B
$$

holds, then

$$
S_{A+B}
=
\ln(\Omega_A\Omega_B)
=
S_A+S_B.
\tag{9}
$$

Thus, $\ln\Omega$ is a natural additive quantity with respect to the product of state counts.

It is important to note that this paper does not unconditionally identify equation (8) with entropy in existing thermodynamics. What is demonstrated at this stage is that even in a fundamental system without spacetime or particles, there can exist a structure isomorphic to the statistical-mechanical relationship between microscopic state count and its logarithm.

------------------------------------------------------------------------

## 6. Energy-Type Readout

Consider a local zero-closure subsystem $I_k$ with $N_k$ elements, where

$$
U^{N_k}=I
$$

is satisfied. Then the fundamental phase step of this subsystem is

$$
\Delta\theta_k
=
\frac{2\pi}{N_k}.
\tag{10}
$$

Using $2\pi$ as a normalized internal scale, we naturally obtain

$$
\varepsilon_k
=
\frac{1}{N_k}
\tag{11}
$$

Thus, for a given local closure decomposition $\Pi$, we can construct an additive readout candidate, for example:

$$
E_{\theta}(\Pi)
=
\sum_{k=1}^{K}\frac{1}{N_k}
\tag{12}
$$

Equation (12) is **not a unique derivation of physical energy**. Rather, it is a minimal additive energy-type measure constructed from the internal phase resolution arising from $U^{N_k}=I$.

Whether it can be adopted as physical energy is determined by **verifying** against actual time evolution, conservation laws, exchange quantities in interaction, and correspondence with observables—thereby establishing the foundation of the theory.

------------------------------------------------------------------------

## 7. Temperature-Type Readout

Using the energy-type readout $E$ to classify states:

$$
\Omega_N(E)
=
\#\left\{
s\in\mathcal{S}_N
\mid
E(s)=E
\right\}
\tag{13}
$$

Then

$$
S_N(E)
=
\ln\Omega_N(E).
\tag{14}
$$

In the continuous limit, we can adopt the definition isomorphic to the standard microcanonical form:

$$
\frac{1}{T_{\mathrm{read}}}
=
\frac{\partial S}{\partial E}
\tag{15}
$$

For a discrete system of finite $N$, more directly:

$$
\frac{1}{T_{\mathrm{read}}}
\simeq
\frac{\Delta S}{\Delta E}
=
\frac{\Delta\ln\Omega}{\Delta E}
\tag{16}
$$

Thus

$$
T_{\mathrm{read}}
\simeq
\frac{\Delta E}{\Delta\ln\Omega}.
\tag{17}
$$

At this stage, we do not identify $T_{\mathrm{read}}$ with Kelvin temperature itself. What has been obtained is the principle that from state counting and energy-type measure, a quantity with the mathematical structure of temperature **necessarily emerges internally**.

------------------------------------------------------------------------

## 8. Dependence on the Observed Local Zero-Closure Subsystem

One of the most critical points in this paper is that the choice of which local zero-closure subsystem to "observe" is not unique.

Suppose observation target $A$ satisfies

$$
\sum_{n\in A}x_n^2=0
$$

and another observation target $B$ satisfies

$$
\sum_{n\in B}x_n^2=0
$$

In general,

$$
N_A\neq N_B,
$$

$$
\Omega_A\neq\Omega_B,
$$

and thus

$$
S_A\neq S_B
$$

may hold.

The energy-type readout similarly gives

$$
E_A\neq E_B
$$

and consequently

$$
T_A\neq T_B
$$

may hold.

Therefore, in this system, physical readout can be conceptually organized as

$$
\boxed{
\text{physical readout}
=
\text{underlying relational state}
+
\text{selected local zero-closure subsystem}
+
\text{readout map}
}
\tag{18}
$$

This is not merely measurement error. For the same global relational state, the state count itself changes depending on which self-consistent local closure is selected as the "system."

------------------------------------------------------------------------

## 9. Common Structure with Curvature Readout: Unified Framework of Observational Readout Hierarchy

In prior work by the author \[2\], a complex spiral wave

$$
z=a+ib
$$

was expanded as

$$
z^2
=
(a^2-b^2)+2iab
$$

The area-type cross term $2ab$ was examined as a candidate for curvature oscillation readout.

In that analysis, from the same local complex wave alone, the observed curvature radius $R$ is not uniquely determined. Rather, the readout scale of curvature varies depending on which physical scale—molecular, celestial, or cosmological—is adopted as the observation target \[2\]. This demonstrates a fundamental structure called the **observational readout hierarchy**.

The present thermodynamic readout also belongs to the same observational readout hierarchy.

The fundamental state need not have

$$
R,\quad E,\quad S,\quad T
$$

written into it as absolute quantities.

From the same self-consistent zero-closure state, the following **distinct observational maps are necessarily derived**:

-   Geometric readout → Curvature radius $R$ (minimal application to gravitational waves)
-   State-count readout → Entropy $S$ (core of this paper)
-   Phase-and-change-scale readout → Energy-type quantity $E$
-   State-density-change-rate readout → Temperature-type quantity $T$

These constitute a parallel structure. Depending on which readout map is adopted, different physical quantities emerge from the same fundamental state.

Thus:

$$
\boxed{
\text{Absolute fundamental structure: self-consistent zero closure}
}
$$

stands in contrast to

$$
\boxed{
\text{Relative physical quantities: observational values from the selected observational readout hierarchy}
}
$$

and this **two-layer hierarchy is necessarily established**.

------------------------------------------------------------------------

## 10. Absolute Zero and Scale Invariance

Equation (1) preserves zero-closure under the transformation

$$
x_n\rightarrow\lambda x_n
$$

Thus, despite lacking an absolute scale,

$$
0
$$

remains as a fixed point of the entire system.

This motivates the consideration that thermodynamic zero (including absolute zero), the reference for curvature, and the reference for development scale may all be read out from the same zero-closure structure.

However, in this paper, we do not identify

$$
T=0,
\qquad
S=0,
\qquad
E=0
$$

with each other.

Specifically, zero closure does not mean absence of information. Non-trivial complex states

$$
(x_1,\ldots,x_N)\neq(0,\ldots,0)
$$

can satisfy

$$
\sum_nx_n^2=0
$$

Therefore, "geometrically zero" is distinguished from "internal information zero."

------------------------------------------------------------------------

## 11. Harmonic Structure and the Euler-Mascheroni Constant

From the finite periodicity condition, the reciprocal structure

$$
\frac{\Delta\theta_k}{2\pi}
=
\frac{1}{N_k}
$$

naturally emerges.

If hierarchical decomposition of local zero closures generates sums of the form

$$
\sum_k\frac{1}{N_k}
$$

then a connection arises with the harmonic numbers

$$
H_N
=
\sum_{k=1}^{N}\frac{1}{k}
$$

The harmonic numbers satisfy

$$
H_N
=
\ln N+\gamma+O\left(\frac{1}{N}\right)
$$

where

$$
\gamma
=
0.5772156649\ldots
$$

is the Euler-Mascheroni constant.

In the author's current numerical experiments, under highly constrained initial conditions, inflation-like evolution that amplifies rapidly from extremely small numerical perturbations up to approximately $0.577$ has been observed. The stability reproduced in double-precision arithmetic is not accidental and strongly suggests that this value originates from the harmonic structure of this model.

These precise numerical behaviors are an important task to be clarified in future research. In particular, the process by which this value converges rigorously to $\gamma$, the true identity of the floor in double-precision calculation, and the question of whether the minimal perturbation capable of development carries a finite value or approaches zero infinitely are essential for determining the physical lower limit and energy-minimization mechanism of this model.

Thus, this section records that the **structural correspondence between the reciprocal structure naturally arising from finite periodicity and local closure decomposition and the value near 0.577 observed in numerical experiments is a phenomenon suggesting the origin of the gravitational constant**, and further numerical verification of these values is the next stage in solidifying the model's foundation.

------------------------------------------------------------------------

## 12. Minimal-Order Consistency Check with the Observed Universe

To verify that the structure of this paper does not demand a direction sharply at odds with the real universe, we perform a minimal-order comparison between baryon number and photon number in the presently observed universe.

From $N$ elements, the number of undirected two-body relations is

$$
M
=
\frac{N(N-1)}{2}.
\tag{19}
$$

Substituting the cosmological discreteness parameter currently under consideration,

$$
N\sim10^{60}
$$

we obtain

$$
M
\simeq
5\times10^{119}
\sim10^{120}.
\tag{20}
$$

Thus, in this system, the relation count $M$ grows approximately as $N^2$ with increasing $N$.

On the other hand, according to the Particle Data Group's standard cosmological parameters, the baryon-to-photon ratio is roughly

$$
\eta
\equiv
\frac{n_b}{n_\gamma}
\simeq
6\times10^{-10}
\tag{21}
$$

and the baryon number density is

$$
n_b
\simeq
2.5\times10^{-7}\ {\rm cm}^{-3}
$$

\[3\]. For CMB temperature $T_\gamma\simeq2.7255\,{\rm K}$, the photon number density corresponding to this is approximately

$$
n_\gamma
\simeq
4.1\times10^2\ {\rm cm}^{-3}
$$

Thus,

$$
\frac{n_\gamma}{n_b}
\sim
1.6\times10^9.
\tag{22}
$$

That is, in the present universe, for each baryon, there are roughly $10^9$ CMB photons.

Applying the same average density to the entire observable universe, as a representative estimate:

$$
N_\gamma\sim10^{89},
$$

$$
N_b\sim10^{80},
$$

and thus

$$
\frac{N_b}{N_\gamma}
\sim10^{-9}
\tag{23}
$$

Here, we do not focus on numerical coincidence between $10^{120}$ and $10^{80}$ or $10^{89}$. These are distinct quantities and must not be directly identified.

The structure to be examined in this model is that

$$
N\sim10^{60}
\quad\Longrightarrow\quad
M\sim10^{120}
$$

**within the vast relational space, only special local structures satisfying self-consistency, zero closure, and finite periodicity simultaneously are read out as baryonic states**, representing the **hierarchical necessity of this theory**. This rarefaction is not a post-hoc fitting of existing parameters but rather a requirement from purely combinatorial probability structure.

If the total number of permitted states is $\Omega(N)$ and the number of states read out as baryonic local closures is $\Omega_B(N)$, we can define

$$
P_B(N)
=
\frac{\Omega_B(N)}{\Omega(N)}
\tag{24}
$$

Actually computing the value

$$
P_B(10^{60})
$$

combinatorially derived from this model is a **decisive step** in verifying whether this theory can quantitatively explain the baryon-to-photon ratio

$$
\frac{n_b}{n_\gamma}
\sim10^{-9}
$$

in the observed universe.

At this stage, the concrete value of $P_B(10^{60})$ has not yet been derived. However, this section demonstrates that

> **the rarefaction of material local states hierarchically derived from the vast relational space not only does not immediately contradict the baryon-to-photon ratio $10^{-9}$ of the observed universe at leading order, but moreover suggests its structural necessity**.

That is, the hierarchical framework of this theory itself may require combinatorial probability structure to explain, rather than post-hoc fitting to, existing cosmological parameters. The next stage is to compute $P_B(10^{60})$ and verify whether quantitative agreement or greater precision with observations can be achieved, thereby determining the true merit of this model.

Moreover, this structure suggests the possibility of a **physical effective upper limit** on $N$: even without a mathematical ceiling, if $N$ becomes too large, the generation probability of material local states drops drastically, making it difficult for the material universe including observers to exist.

------------------------------------------------------------------------

## 13. Relation to General Relativistic Readout

Zero closure

$$
\sum_nx_n^2=0
$$

can be partly transposed to the opposite side and, for the selected subsystem, read as

$$
\sum_i x_i^2=R^2
$$

to construct a mapping to constant-curvature geometry with curvature radius $R$.

In $D$-dimensional constant-curvature spacetime,

$$
R_{\mu\nu}
=
\frac{D-1}{R^2}g_{\mu\nu},
$$

$$
\mathcal{R}
=
\frac{D(D-1)}{R^2}.
$$

Substituting into the vacuum Einstein equation

$$
R_{\mu\nu}
-\frac{1}{2}\mathcal{R}g_{\mu\nu}
+\Lambda g_{\mu\nu}
=
0
$$

yields

$$
\Lambda
=
\frac{(D-1)(D-2)}{2R^2}.
\tag{25}
$$

For $D=4$,

$$
\Lambda
=
\frac{3}{R^2}.
\tag{26}
$$

Here, $\Lambda$ was not added to the fundamental system as an independent constant. Rather, it is obtained as the curvature scale when the observational mapping that reads zero closure as constant-curvature spacetime is selected.

This structure is isomorphic to the thermodynamic readout of this paper. That is, we distinguish the zero closure of the foundation itself from the physical quantities read out from it for the selected subsystem.

This paper does not claim to derive the entire field equations of general relativity from equation (1). What is demonstrated is the **structural correspondence** that the cosmological constant scale corresponding to the vacuum constant-curvature solution can be naturally expressed as a geometric readout from zero closure.

------------------------------------------------------------------------

## 14. Scope of Claims and Falsifiability

The claims of this paper are limited as follows:

1. If a self-consistent complex relational system possesses a zero-closure condition, the global system can be decomposed into multiple local zero-closure subsystems.
2. If distinct self-consistent states realizing the same global zero closure can be counted, a state count $\Omega$ is defined.
3. From $\Omega$, an additive entropy-type readout $S=\ln\Omega$ can be constructed.
4. From the phase periodicity arising from $U^{N_k}=I$, an additive energy-type measure can be constructed.
5. From state count and energy-type measure, a temperature-type readout $T^{-1}=\Delta S/\Delta E$ can be constructed.
6. These values may vary depending on which local zero-closure subsystem is selected as the observation target.
7. For $N\sim10^{60}$, reaching $M\sim10^{120}$, the vast relational space presents a clear numerical problem for verifying the rarefaction of special material local states.

On the other hand, the following remain undeduced:

-   Uniqueness of the energy-type measure $E$.
-   Complete equivalence relations and measure on the permitted state set $\mathcal{S}_N$.
-   Explicit counting of $\Omega(N)$ and $\Omega_B(N)$.
-   Derivation of $P_B(10^{60})\sim10^{-9}$.
-   Rigorous convergence to the Euler-Mascheroni constant $\gamma$.
-   Physical identity between absolute zero and the zero-closure fixed point.
-   Derivation of thermal equilibrium, the zeroth law of thermodynamics, canonical distribution.
-   Derivation of arbitrary spacetime metrics and the complete Einstein equation of general relativity.

Therefore, this paper is not a "completed theory of thermodynamics," but rather a **research note establishing a minimal framework for constructing thermodynamic readouts from relational structures prior to spacetime**.

The next stage of falsifiable steps is clear, forming the core verification items of this theory. In particular, numerically computing

$$
\Omega(N),
\qquad
\Omega_B(N),
\qquad
P_B(N)
$$

is not merely a future task but a **highest-priority verification stage demonstrating the hierarchical self-consistency of this model**. Should this counting demonstrate quantitatively that material local states rarefy with $N$, the cosmological interpretation of this paper gains strong support. Conversely, should rarefaction not occur, a fundamental reconsideration of this framework itself becomes necessary.

------------------------------------------------------------------------

## 15. Conclusion

In this paper, we have demonstrated the construction:

$$
\text{complex relations}
\rightarrow
\text{self-consistency}
\rightarrow
\text{zero closure}
\rightarrow
\text{local zero-closure decomposition}
\rightarrow
\Omega
\rightarrow
S
\rightarrow
E_{\mathrm{read}}
\rightarrow
T_{\mathrm{read}}
$$

This construction requires no prior placement of space, time, particle, mass, or temperature as fundamental concepts.

What is most important is the distinction between the self-consistent zero closure of the whole and the local zero-closure subsystem selected as the observation target. From the same fundamental state, state count, entropy, energy scale, and temperature scale can vary depending on which partial closure is read as the system.

Therefore, the minimal paradigm proposed by this paper is:

$$
\boxed{
\text{absolute structure}
=
\text{self-consistent zero closure}
}
$$

$$
\boxed{
\text{physical quantities}
=
\text{readouts from selected local zero-closure subsystems}
}
$$

The curvature readout examined in prior work \[2\] and the thermodynamic readout of this paper belong to the same **observational readout hierarchy**. Unifying both papers:

$$
\boxed{
\text{Fundamental self-consistent zero closure}}
\xrightarrow{\text{distinct observational maps}}
\boxed{
\text{Multiple-layer physics: curvature, thermodynamics, gravitational waves, etc.}}
$$

establishes this structure.

Furthermore, while a cosmological discreteness of $N\sim10^{60}$ generates a relation space of $M\sim10^{120}$, the baryon-to-photon ratio of the observed universe is approximately $10^{-9}$. Although the probability rule linking the two has not yet been derived at this stage, the direction—wherein material local states become sparse within the vast relational space—is not immediately contradicted by the observed universe at leading order. This is not mere predictive agreement but **concrete evidence of the hierarchical self-consistency of this model**.

The next task is not to add further conceptual correspondences, but to actually count the number of permitted states $\Omega(N)$ and the material-state probability $P_B(N)$, and **verify** that this theoretical prediction is quantitatively consistent with the observed universe.

------------------------------------------------------------------------

## References

\[1\] L. Boltzmann, "Über die Beziehung zwischen dem zweiten Hauptsatze
der mechanischen Wärmetheorie und der Wahrscheinlichkeitsrechnung,
respective den Sätzen über das Wärmegleichgewicht," *Sitzungsberichte
der Kaiserlichen Akademie der Wissenschaften,
Mathematisch-Naturwissenschaftliche Classe*, Vol. 76, pp. 373--435,
1877. English translation and commentary: K. Sharp and F. Matschinsky,
*Entropy* **17**, 1971--2009 (2015). DOI: 10.3390/e17041971.

\[2\] N. Kihara, "Curvature Oscillation Readout from Area Cross Terms of
Complex Spiral Waves—Minimal Application to Gravitational Waves,"
Research Note / Hypothesis and Exploratory Study, 2026-09-01. DOI:
10.5281/zenodo.22230941.

\[3\] Particle Data Group, "Astrophysical Constants and Parameters,"
*Review of Particle Physics* / PDG data tables. Representative values
used here: baryon-to-photon ratio
$5.8\times10^{-10}\lesssim\eta\lesssim6.5\times10^{-10}$ and baryon
number density $n_b\simeq2.515\times10^{-7}\,\mathrm{cm}^{-3}$ for
standard cosmological parameters.

------------------------------------------------------------------------

## Note

This is a hypothesis and exploratory research note. We distinguish between formal coincidence of equations and physical identity, and make clear when correspondences remain undeduced. In particular, energy-type readout, baryon generation probability, and connection to the Euler-Mascheroni constant remain subjects of future numerical and analytical verification.
