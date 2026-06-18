# On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit
## — A Limited Verification Report from a Minimal Axiom System (Revised)

**Author:** Noriaki Kihara (ORCID: 0009-0004-6753-4020)
**Date:** 18 June 2026 (revised)
**DOI (this version, v4):** 10.5281/zenodo.20746811 (v1: 10.5281/zenodo.20740842; v2: 10.5281/zenodo.20742386; v3: 10.5281/zenodo.20743142)
**Concept DOI (all versions):** 10.5281/zenodo.20740841
**Zenodo:** https://zenodo.org/records/20746811
**Revision history (2026-06-18):** v2 = wording fix in §6.3 (what is preserved is the product of observables, not the product of fluctuations) + added a list of related supplements (3). v3 = expanded the supplements list to 6 (added complex structure, g/ω separation, discreteness). v4 = §7 reframed to canonicalization (privileged-scale separation, fully invertible, no full symmetry required; removed full-symmetry-dependent wording), added the forward / inverse map formalization, and expanded the supplements list to 7 (added the map appendix (ν,δ)↔(N,Δ)).

---

## Note to the reader (please read first)

This report does not claim the "discovery" of a new general fact. Because it starts from a deliberately limited set of premises, it is easily misread on a first pass. Most misreadings arise from **blurring the system-interior (the unknowable, anonymous internal viewpoint) and the system-exterior (the observer's viewpoint)**. Moreover, in this system the question "which side is real?" is ill-posed (see the note in §0.5); the objective anchors are placed not on "reality" but on fixed points and conserved quantities. To prevent misreading, each statement is, in principle, prefixed with a tag indicating whether it is an interior or an exterior statement. Please read §0.5 (interior and exterior) and §1 (conventions) first.

---

## Abstract

This report begins with a naive question about the conjugate relation of two quantities, νλ=k (a quantity read off by an observer). This rigid relation alone generates no structure. As the only crevice through which structure can enter, we posit that the integer-valued ν carries an intrinsic fluctuation, and that the system possesses high anonymity and high symmetry. In a limited regime requiring extremely high symmetry and anonymity, we report a **result within a limited scope**: a squared quantity representing scale can be mapped, through a logarithmic connection, onto a flat additive structure (a flat hyperplane).

Introducing the fundamental quantity as the logarithm of a squared quantity, N=log₂ν, the conjugate relation νλ=k becomes an additive identity ΣNᵢ=K, which is a flat hyperplane in any dimension. By canonicalization (separating the conserved quantity as a privileged scale C and taking the normalized shape N̂=log₂(ν/C)), the map ν↔N is fully invertible: both the individual values and the fluctuation δ are recovered exactly (full symmetry is not required). What matters is that, while this logarithmic map factors scale (g) out as a privileged, separable coordinate, the fluctuation Δ (the product / area / ω) remains invariant. That is, the map **separates scale (g) from quantization (ω)**. The exact formalization of the forward / inverse maps is given in the appendix "The maps (ν, δ) ↔ (N, Δ)." Here R / C is the scale (magnitude) of the system, and is neither a curvature tensor nor the radius of a sphere.

This report treats the magnitude of the fluctuation, ±1/2, as a **posit, not a derivation** (reasons in §6). Whether the result remains a special case of known frameworks (Weyl transformation, conformal flatness, etc.) cannot be determined here, and an exhaustive search for counterexamples has not been carried out. The purpose of this report is confined to recording this limited structure.

---

## §0 Motivation (author's note)

The starting point of this report is a naive question about the conjugate relation νλ=k.

[Exterior (observer)] Writing νλ=k gives λ=k/ν: as ν grows, λ shrinks. This is a rigid see-saw relation read off from the outside.

[Exterior (observer)] But this rigid relation alone yields no structure from the two quantities ν and λ. There is merely a single proportionality.

Where, then, can structure enter? This report loosens two points.

[Discrete structure] First, ν is discrete (integer-valued) and has no continuous fractional part. It is not that "a continuous value really exists but appears integral due to the limits of observational precision"; rather, **the fractional part does not exist** (this "does not exist" is a statement about the system-exterior, i.e., the ontic level, and is distinct from the interior indistinguishability discussed later). Yet the integer value is not fixed: it carries an intrinsic fluctuation. This is not observational error but a fluctuation intrinsic to the system. (Integrality itself is not a claim about reality, but a consequence of discreteness; its derivation is in §3.)

[Interior] Second, in νλ=k there is no name distinguishing ν from λ. From within the system, which is which cannot be told. The system thus possesses high anonymity.

This report deliberately builds the system under this high anonymity and high symmetry. Ordinarily, imposing such strong anonymity and symmetry constrains the degrees of freedom and makes structure hard to arise. The aim of this report is to see what nonetheless stands up.

[Exterior (observer) / Interior] One last point. νλ=k is, after all, a quantity read off from the observer's viewpoint, exterior to the system. Apart from it, might there be an invariant (a conserved quantity) or a fixed point (a stable point) intrinsic to the interior? This is the question running through the report.

---

## §0.5 Author's note: interior and exterior (the unknowable and invariants)

Throughout, this report strictly distinguishes the system-interior from the system-exterior. Blurring this causes most of the statements below to be misread.

- [Interior] The interior is the viewpoint of an anonymous element (one with no name, color, or position). What cannot be distinguished within the interior is treated as equivalent to *null* (nothing) within the interior. From the interior, exterior quantities cannot, in principle, be observed directly (except when an exterior quantity manifests inside as an effect; see Note 1).
- [Exterior] The exterior is the viewpoint of one who observes from outside. The observer obtains only a limited set of observable quantities, and does not necessarily reach the interior invariants (fixed points, conserved quantities). Hence the observer cannot assert what is (or is not) inside the system.

**[Important note: reality cannot be asked.]** This system rests on a dual / conjugate relation (νλ=k). Hence the question "which of ν and λ is real?" is ill-posed — declaring one real makes the reality of the other ambiguous. Nor do we deny interior reality. What carries meaning is not "which is real" but what is held invariant under the duality — the **fixed points (k=1, ν=1) and conserved quantities (the product δ₁δ₂=1/2, the sum of squares Σν²)**. In what follows, all objective anchors are stated as fixed points / conserved quantities, and no claim is made about the locus of reality.

Under this distinction, we derive the most basic counting from the interior viewpoint. The only thing registrable from the interior is "whether there is something other than oneself, or not."

- [Interior] The case 0: there is no subject to judge a state. Hence 0 cannot be distinguished from null. This only states interior indistinguishability and does not exclude 0 as an integer-structure entity.
- [Interior] The case 1: the sole element has no "other than itself." Being anonymous, "I exist" cannot even stand without a contrast against an other. Hence, in the interior, 1 cannot be distinguished from 0. This does not exclude 1 as an integer-structure entity.
- [Interior] The case 2: for the first time there is an "other than oneself." Through it, "I exist" stands. Two is the minimum at which anything can be known from the interior.
- [Interior] The case 2-or-more: the other carries no color, name, or position, and so cannot be counted. Hence, in the interior, 2 and 2-or-more cannot be distinguished. This does not exclude 2-or-more as an integer-structure entity.

[Interior] Therefore the interior viewpoint collapses to a binary — "no other" (0 and 1 indistinguishable) and "an other present" (2 and 2-or-more indistinguishable).

[Conserved quantity] On the other hand, as a structure of conserved quantities / fixed points, the integer sequence 0,1,2,3,… stands. The interior collapse to a binary is a collapse by unknowability, and does not negate this integer structure.

[Interior / Exterior] In one line: what is unknowable in the interior is, in the interior, equivalent to null. But it cannot be negated from the exterior either — what is unknowable is not even an object of negation.

**Note 1 (an example of an exterior quantity manifesting inside as an effect):** spatial curvature in a multidimensional spacetime is such an example. Curvature itself is an exterior quantity, and the interior cannot perceive its direction; but its effect — for instance, the sum of the interior angles of a triangle deviating from 180° — can be perceived from the interior. In the present system, however, this example does not yet clearly arise, and curvature lies outside the scope of this report (§7, §8).

---

## §1 Conventions (how to read interior/exterior)

We now enter the formal description. For readability, two conventions are set.

1. **[Tag]** Each statement is, in principle, prefixed with a tag indicating whether it belongs to the interior (the unknowable, anonymous internal viewpoint) or the exterior (the observer's viewpoint). Anchors held invariant under the duality are marked as [Fixed point], [Conserved quantity], or [Discrete structure]; no claim is made about the locus of reality.
2. **[Verb discipline]** When stating interior indistinguishability, we do not write "X does not exist" but "X cannot be distinguished in the interior." "Does not exist" is allowed only for objective statements about discrete structure / fixed points / conserved quantities, and is not used for interior indistinguishability.

**[Declaration]** Among the axioms and derivations below, all statements about 0, 1, 2 are statements about interior distinguishability, and do not negate the integer structure (discreteness / fixed points / conserved quantities).

---

## §2 Axioms (assumptions)

The system rests on the following axioms. No external premises or references to prior literature are used. Physical interpretation is kept separate; we restrict ourselves to statements about structure.

**Axiom 1 (existence of difference):** There is difference. Only "there is another quantity" can be said. What differs, which differs from which, and how many there are cannot be said (anonymity). Speaking of difference itself presupposes difference, so this is self-referential.

**Axiom 2 (the limit of anonymity = 2):** [Interior] The minimum at which difference stands is two. By anonymity, the maximum is also two. In the interior it closes on 0, 1, 2, with 2 as the ceiling (by the derivation in §0.5: 0 and 1 cannot be distinguished, and 2 and 2-or-more cannot be distinguished). Treating 3-or-more *as* 3-or-more is counting = labeling = the genesis of dimension, and is a breaking of anonymity. This is the limit of interior distinguishability, and does not exclude 3-or-more as an integer-structure entity.

**Axiom 3 (fluctuation):** Each quantity has an intrinsic fluctuation δ. δ is an attribute of that quantity and is a scalar. The magnitude ±1/2 of the fluctuation is a posit, not a derivation (§6).

**Axiom 4 (conservation of the product: the edge of being / non-being):** [Interior] The product of the two fluctuations is conserved: δ₁δ₂=1/2. This is the invariant of the edge of the binary of existence ("being / non-being") (product, hyperbolic SO(1,1), area), and corresponds to the imaginary part ω of the complex inner product introduced later. (Which one is taken as the imaginary part is, by the principle of anonymity, itself a matter of choice and is arbitrary.)

**Axiom 5 (conservation of the sum of squares: the binding of difference):** [Interior] The system is bounded and symmetric (anonymity invariant under relabeling). This makes an invariant binding the two on equal footing stand: Σν²=R² (sum of squares, SO(2) symmetry, scale R = distance g). Here SO(2) refers only to the symmetry preserving the sum of squares, and does not assume a circle / sphere as a geometric image. This is the invariant of the existence of difference ("there is another quantity"), and corresponds to the real part g of the complex inner product introduced later. (Which one is taken as the real part is, by the principle of anonymity, also arbitrary.)

**Lemma (closure is not an axiom but a consequence of anonymity):** [Interior] An end (boundary) is a privileged point with a neighbor only on one side, and violates anonymity. Hence a system preserving anonymity cannot have ends and closes — moreover, it closes periodically (rotationally symmetric). Therefore "closed system" is not an independent axiom but is derived as a consequence of anonymity (Axioms 1, 2). The closure assumption placed as an independent axiom in the old version is removed in this revision.

**Note (mutual underivability of Axioms 4 and 5):** [Interior] The product (Axiom 4, ω) and the sum of squares (Axiom 5, g) cannot be derived from each other. This is not a technical matter: "whether it is" and "whether it differs" are two root distinctions that stand independently under anonymity. If one followed from the other, the limit of anonymity would be wider than two.

---

## §3 Derivation of discreteness (discreteness is not an axiom but a derived quantity)

The discreteness of the system is not posited as an axiom. It is derived in two stages.

**Stage one (fundamental discreteness: from Axiom 2).** [Interior] In the interior, 1 does not stand alone (Axiom 2). That is, the unit "1" for subdividing a continuum, lacking contrast, cannot be distinguished from 0. A continuum (continuous cardinality) stands only when "a unit that stands alone can be subdivided infinitely," but that unit does not stand. Hence what can be constructed in the interior is only the repetition of binary presence/absence judgments = doubling / halving (×2, ÷2), and counting is discrete (integer). The base being 2, integrality, and the half-bit introduced later all descend from this "1 does not stand alone."

**Stage two (discreteness of the mode number: from boundedness).** [Discrete structure] By the lemma of §2, the boundedness of Axiom 5 appears as a **periodically closed system** (a system with no ends that returns upon going around once). That it is not a string [0,L] with ends is consistent with the lemma (ends violate anonymity). For a standing wave to stand on a periodically closed system (circumference L), the periodic boundary condition must be satisfied, and the admissible wavelengths are quantized to integer fractions of the circumference (λ_n=L/n, n=1,2,3,…). The spectrum is numbered by the integer n. The genesis of this integer mode number n is counting = the genesis of dimension, and is the mechanism that breaks through the ceiling (2) of Axiom 2 toward 3-or-more.

That is, discreteness is the name for what Axiom 5 does to Axioms 1–4, and is not an independent principle.

---

## §4 Introduction of the fundamental quantity N

[Conserved quantity] The conserved quantity is not the linear observables ν, λ but the squared quantity Σν²=R² (Axiom 5). The linear quantities ν, λ are observables that swing reciprocally under the conjugate relation νλ=k, and behind their product an invariant squared quantity is conserved.

We introduce the fundamental quantity as the logarithm of a squared quantity. The base is taken to be 2 by Axiom 2 (the interior binarity). Attaching a factor 1/2 gives an integer sequence:

$$N \equiv \tfrac{1}{2}\log_2(\nu^2) = \log_2\nu .$$

[Fixed point] When ν=2ⁿ, N=n (integer). N=0 corresponds to ν=1 (the self-dual fixed point). [Interior] In the interior, ν=1 cannot be distinguished from ν=0 (Axiom 2), but this does not negate 1 as a fixed-point / integer-structure entity. N represents the whole integer set ℤ, symmetric in sign about 0. The sign of N corresponds to which side of the conjugate pair one views from.

[Interior / Exterior] The imaginary unit i is merely a symbol folding two real components with a 90° phase difference into complex notation; the axioms of this system close over the reals (complex numbers are not introduced). As noted later (§7 note), the complex structure arises as a consequence of this real two-ness.

---

## §5 The additive identity of the conjugate relation (the g side)

[Exterior (observer)] Taking the base-2 logarithm of the conjugate relation νλ=k turns the product into a sum:

$$\nu\lambda=k \;\xrightarrow{\ \log_2\ }\; N_\nu + N_\lambda = K \quad(K\equiv\log_2 k).$$

This is an identity, and the conserved quantity is K (the sum of the N of the conjugate pair). Whichever of the three terms N_ν, N_λ, K is chosen on the right-hand side,

$$N_\lambda = K - N_\nu,\qquad N_\nu = K - N_\lambda,\qquad N_\nu + N_\lambda = K$$

are preserved as identities, and the only thing produced by the interchange is a relabeling — a subtraction of logarithms (i.e., division). [Interior] Hence, without privileging any quantity, anonymity is maintained.

[Exterior (observer)] A change of scale (a change of k) corresponds to a translation of the N axis, and scale invariance appears as translational symmetry. At k=1, K=0 and N_λ=−N_ν (sign-reversal symmetry about 0).

**Note (the identity generates no structure):** This additive identity is the very definition of the logarithm, and by itself generates no new structure. What generates structure is the fluctuation Δ (the ω side) of the next section. This identity is merely a rephrasing of the g side (distance, sum of squares).

---

## §6 The fluctuation Δ and uncertainty (the ω side, the side that generates structure)

What generates structure is the fluctuation Δ. To place this precisely, we first kill two misreadings.

### 6.1 Denial of hidden variables / observational error

[Discrete structure] ν is discrete (integer-valued) and has no continuous fractional part. It is not that "a continuous value ν=2.0±δ really exists and appears integral due to the limits of observational precision"; rather, the system is discrete, so ν is integral. This is not a claim about reality but a consequence of discrete structure (see the note in §0.5).

ν=2 means that ν fluctuates while remaining an integer, with its central value (median) being 2. For example, the value can take the neighboring integers 1, 2, 3, with center 2. The value itself jumps among integers. This is a fluctuation intrinsic to the system, not an observational error riding on a fixed true value.

[Interior] The interior quantity itself is an integer. Hence it cannot take a fractional value such as 0.5. If fluctuation appears in the interior, it appears only as an integer-step shift, i.e., ±1.

### 6.2 ±1/2 is a posit, not a derivation

**[Negation]** One is tempted here to say "the standard deviation is ±1/2." But this report does not derive it. There is, in this system, no basis on which to even construct the concept of a standard deviation.

**[Reason]** Defining a standard deviation requires an operation of observing the same object repeatedly and aggregating its scatter. But —

- [Interior] By anonymity, observations cannot be individuated from one another (this is the first, this is the second).
- Since no time evolution is defined, there is no order in which to arrange observations.
- It is not assumed that the true value is preserved the same between observations.

Hence the very phrase "observe multiple times" is self-contradictory in this system. That single phrase smuggles in many hidden premises — individuation, temporal order, persistence of value.

**[Posit]** Therefore ±1/2 is not a computed value but a posit (an assumption that it is presumably so). Whatever method one uses to derive 1/2, one necessarily smuggles in one of the hidden premises (individuation, temporal order, persistence). It is not that the standard deviation cannot be computed exactly; it is that the very basis for constructing the concept of a standard deviation is absent. This report places ±1/2 as an assumption, not a derivation.

### 6.3 Conservation law of the fluctuation (sum-zero is the definition, constant-product is derived)

**Note (the fluctuations of §6.1 and §6.3 are different things):** The "the interior integer quantity fluctuates in integer steps (±1)" stated in §6.1 and the half-bit fluctuation Δ=±1/2 placed in N-space in this §6.3 are **separate posits; one is not derived from the other**. No relation (identity, derivation, or approximation) is claimed between them. §6.1 is a statement about the fluctuation of the interior integer quantity; §6.3 is a statement about the fluctuation posited on the side of the fundamental quantity N. This report treats them as independent assumptions, and leaves their reconciliation or unification as future work. Joining what is not joined inevitably smuggles in a hidden premise (the same discipline as ±1/2 in §6.2), so this report does not join them.

[Interior] We place the fluctuation on the side of the fundamental quantity N (positing the fluctuation in N-space makes it scale-dependent in the observable space, which is natural as a fluctuation of integer counting). Each fundamental quantity is given a fluctuation Δ=±1/2 (a half-bit; a posit; scale-invariant). The conjugate pair satisfies the conservation of the sum N_a+N_b=K, and the fluctuations satisfy

$$\boxed{\ \Delta_a + \Delta_b = 0\ }$$

(which is taken as the definition).

[Exterior (observer)] By the inverse map a=2^{N_a+Δ_a},

$$k = ab = 2^{(N_a+N_b)+(\Delta_a+\Delta_b)} = 2^{K+0} = 2^K \quad(\text{constant})$$

is obtained, and the product of the observables ab is held constant under fluctuations (the lower bound of minimal uncertainty). That is, the uncertainty relation in the observable space (constant product) is derived from the sum-zero of the fluctuations in N-space. Note that what is preserved is the product of the observables ab, not the product of the fluctuations δ₁δ₂ (the latter depends on the allocation and is not constant; see the appendix "From the Sum-Zero of Fluctuations to Uncertainty"). The minimal unit of uncertainty, 1/2, corresponds to the posited half-bit ±1/2. Note that a single fluctuation ±1/2 and the product of two fluctuations δ₁δ₂=1/2 (Axiom 4) are not the same quantity; they are placed in correspondence as quantities belonging to the same ω (area) side.

---

## §7 Main result: the logarithmic map of scale and the g/ω separation (a limited regime)

**[Note: on "scale" here.]** The object of flattening in this section is R in Σν²=R², i.e., the scale (magnitude) of the system. This is neither a curvature tensor nor the radius of a sphere. Σν²=R² is a conservation law stating that the sum of squares equals R², and does not require a circle or sphere geometry. This report does not claim "flattening of curvature" but states "linearization (flattening) of scale by the logarithm." Extension to a curvature tensor lies outside the scope of this report (§8, §9). This report does not claim that the continuum limit of this system extends to a curvature tensor; rather, it regards it as likely to hold only in the discrete system.

### 7.1 Complementarity of linearization in the two spaces

[Exterior (observer)] The conserved quantity of Axiom 5 (an interior invariant) appears to the observer as the sum of squares Σν²=R² (g, scale R), and this is the conservation law representing scale. The conjugate relation νλ=k is a product (nonlinear). In the N-space via the logarithm, this relation is reversed. νλ=k is additive (N_ν+N_λ=K, §5), i.e., a flat hyperplane, whereas the sum of squares Σν² is in general a sum of exponentials Σ2^{2N_n} and is not linearized in N.

### 7.2 Flattening by canonicalization (no full symmetry required; fully invertible)

[Exterior (observer)] The conjugate (product) relation becomes, under the logarithm, the additive identity: ΣNᵢ=K for n variables (§5). This is a **flat hyperplane** in n-dimensional space, and it **holds in any dimension and does not require full symmetry (that all ν_n be equal)**.

**Canonicalization:** separate the conserved quantity as a privileged scale C=k^{1/n} (=2^{K/n}) and take the normalized shape N̂ᵢ=log₂(νᵢ/C) (ΣN̂ᵢ=0). The map ν↔N is then **fully invertible**: both the individual νᵢ and the fluctuation δ are recovered exactly. The privileged scale C is the sole coordinate passed through unchanged; every other quantity is mapped by the same logarithm. The exact formalization of the forward map N̂ᵢ=log₂(νᵢ/C) and inverse map νᵢ=C·2^{N̂ᵢ}, the conserved quantities, and the handling of the privileged scale, are given in the appendix "The maps (ν, δ) ↔ (N, Δ): recovery of values and conserved quantities."

Flattening means **factoring scale (g) out as a privileged, separable coordinate, leaving the rest as a flat additive structure**. It discards no information (fully invertible) and requires no full symmetry. Furthermore, since the system is scale-invariant, scale dependence is flattened as a straight line along the translation.

### 7.3 The heart of the map: separation of scale (g) and quantization (ω)

This is the heart of the report. The flattening does not "throw away the scale information"; it factors scale out as a privileged coordinate.

[Exterior (observer)] Scale (g) — appearing to the observer as the sum of squares Σν² (distance / metric, the Born-norm side), and in canonicalization as the privileged scale C — can be **separated and factored out as a flat coordinate (gauge)** by the logarithmic map. The g side carries scale but not quantization.

[Interior] On the other hand, the fluctuation Δ = δ₁δ₂=1/2 (Axiom 4) is area / the symplectic form ω (the imaginary part of the complex inner product), the side that carries quantization. Note that δ is mapped to Δ by the same map as ν (appendix "The maps (ν, δ) ↔ (N, Δ)"), but **Δ=±1/2 is scale-invariant**, so after the privileged scale C is separated it persists on the flat surface as an invariant. This is what "ω remains invariant" means.

That is, the logarithmic map **separates scale (g, gauge, separable) from quantization (ω, area, invariant)**. Scale can be factored out as a privileged coordinate (discard it and a scale-invariant structure remains); quantization resides invariantly in that structure. **This separation requires no full symmetry; it holds in any dimension, on a fully invertible map.** The physics (quantization) was, from the start, not in g (scale) but resided in ω (the fluctuation Δ). Here "separation" means, operationally, that the map carries both ν and δ, but separates g out as a privileged scale while ω (Δ) remains invariant because Δ is scale-invariant; it is not a claim that it is a canonical transformation strictly preserving the symplectic form (that rigorization is future work).

### 7.4 Structure of the map

[Exterior (observer)] The map has the following structure (the self-dual fixed point k=1 may be taken as reference; full symmetry is not required):

- the scale-bearing system in the observable space (ν-space)
- is separated, via the logarithmic map, into the privileged scale C and the normalized shape N̂ (a flat hyperplane ΣN̂=0),
- and the individual values are **exactly restored** by the inverse map νᵢ=C·2^{N̂ᵢ}=2^{Nᵢ} (fully invertible).
- And the fluctuation Δ (the ω side) is preserved invariant before and after the map.

The exact formalization of the forward map N̂ᵢ=log₂(νᵢ/C) and inverse map νᵢ=C·2^{N̂ᵢ}, with the conserved quantities and the handling of the privileged scale, is given in the appendix "The maps (ν, δ) ↔ (N, Δ): recovery of values and conserved quantities."

**Note (the complex structure as a consequence of real two-ness):** [Interior / Exterior] Reinterpreting ν₁, ν₂ as a radius ν and a phase θ gives ν₁=ν cosθ, ν₂=ν sinθ, and the two are geometrically isomorphic to the polar representation of z=ν e^{iθ}. The sum of squares ν²=ν₁²+ν₂²=|z|² corresponds to the real part g, the area of the fluctuation δ₁δ₂=1/2 to the imaginary part ω, and the cos↔sin rotation to the complex structure J. Complex numbers are not placed in the axioms; the complex structure arises as a consequence of real two-ness (the "2" of Axiom 2) and the sum of squares (Axiom 5).

---

## §8 Limitations

The scope of this result is strictly limited.

1. **Scope of the result.** The flattening of the additive structure ΣNᵢ=K and the full invertibility of the map (canonicalization: privileged-scale separation, forward / inverse maps) hold in any dimension and for any k (no full symmetry required). What lies outside the scope is the extension to a curvature tensor (item 5 below) and the linearization of the sum of squares Σν² itself (in N-space it is in general a sum of exponentials Σ2^{2N} and is not linearized, §7.1).

2. **±1/2 is a posit and has not been derived.** As in §6.2, this system has no basis on which to construct the concept of a standard deviation. Every computation deriving ±1/2 smuggles in a hidden premise. Note that the main result of this report (the g/ω separation of §7.3) presupposes this posit ±1/2. If the posit fails, the conclusion changes.

3. **No exhaustive search for counterexamples has been carried out.** It has not been shown that the result has no counterexamples.

4. **The relation to known frameworks is undetermined.** Whether the result remains a special case of known frameworks (Weyl transformation, conformal flatness, conformal cosmology, etc.) or exceeds them cannot be determined at present. In particular, how the g/ω separation of §7.3 differs from or agrees with known Kähler geometry / conformal flatness has not been rigorously compared.

5. **This report treats only the scale R, not a curvature tensor.** The R of Σν²=R² is the scale (magnitude) of the system, not a curvature tensor or the radius of a sphere. A rigorous proof including the differential structure (the curvature tensor) is not given, and extension to curvature / general curved spacetimes lies outside the scope of this report.

6. **No physical identification is made.** This system is an algebraic / geometric (informational) structure, and no physical interpretation (time, space, energy, particles, etc.) is attached to the observables / fundamental quantities. Isomorphism with a physical theory is a separate task.

7. **No relation or unification of the two fluctuations is given.** The integer-step fluctuation (±1) of §6.1 and the half-bit fluctuation Δ=±1/2 of §6.3 are separate posits, and no relation (identity, derivation, approximation) is claimed between them. The derivation of placing the fluctuation in N-space from the interior integer fluctuation is not carried out in this report. Unification is left as future work.

---

## §9 Positioning relative to known frameworks (for reference)

Computing physical quantities in a curved spacetime is generally difficult. One existing approach is the Weyl transformation $g_{\mu\nu}\to\Omega^2 g_{\mu\nu}$; but the standard action is not invariant under Weyl rescaling, and the curvature scalar transforms as

$$R \to e^{-2\sigma}\left(R - 6\,\Box\sigma - 6\,g^{\mu\nu}\partial_\mu\sigma\,\partial_\nu\sigma\right),$$

so the curvature information remains as derivative terms of the rescaling function σ. That is, even after flattening, the curvature does not entirely vanish. Conformal flatness always holds locally in two dimensions but not for general metrics in higher dimensions.

From a starting point different from these known frameworks — only a few axioms — this report presents the structure that, in a limited regime, a quantity representing scale (not a curvature tensor) is mapped onto a flat structure leaving no derivative terms, and quantization (ω) is preserved invariant under the map. What this report treats differs in scope from the Weyl transformation, which rescales a curvature tensor: it is a logarithmization of scale. It is neither a rehash of Weyl nor a claim to treat curvature. Extension to curvature is left as future work. No generalization is claimed. The possibility that it is a special case of a known framework is not excluded either.

---

## Supplements

This report has appendices that concretize points compressed in the main text and the "arising structure." All are `isSupplementTo` this report.

1. **A Magic-Coin Toy Model** — visualizes the limit of anonymity (0, 1, 2) and the denial of hidden variables of §0.5 as a tabletop magic. Concept DOI: [10.5281/zenodo.20741264](https://doi.org/10.5281/zenodo.20741264)
2. **Multidimensional Structure from the Logarithmic Fundamental Quantity** — §3, §5, §7. Exhibits the (n−1)-dimensional hyperplane defined by the additive identity ΣNᵢ=K (the minimal example of the "arising structure"). Concept DOI: [10.5281/zenodo.20741712](https://doi.org/10.5281/zenodo.20741712)
3. **From the Sum-Zero of Fluctuations to Uncertainty** — §6.3. Makes explicit the direction in which the constant product of observables is derived from the sum-zero in N-space (what is preserved is the product of observables, not the product of fluctuations). Concept DOI: [10.5281/zenodo.20742277](https://doi.org/10.5281/zenodo.20742277)
4. **The Complex Structure is a Consequence, Not a Posit** — §4, §7.4. Without positing the imaginary unit i, the complex structure follows from two real components and the rotation (the 90-degree rotation J, J²=−1). Concept DOI: [10.5281/zenodo.20742984](https://doi.org/10.5281/zenodo.20742984)
5. **The Logarithmic Map Separates Scale and Quantization** — §7.3 (the core of this report). Flattening is not a loss of information but a separation of g (scale, flattened) and ω (quantization, invariant). Concept DOI: [10.5281/zenodo.20742988](https://doi.org/10.5281/zenodo.20742988)
6. **Discreteness Arises Not from an Outer Bound but from Nodes** — §3. Discreteness arises intrinsically from the nodes generated by Axiom 2 (binarity), not from boundedness (avoiding the circularity of borrowing the string's boundary condition). Concept DOI: [10.5281/zenodo.20742990](https://doi.org/10.5281/zenodo.20742990)
7. **The maps (ν, δ) ↔ (N, Δ): recovery of values and conserved quantities** — §5, §7. Rigorous formalization of the forward map N̂=log₂(ν/C) and inverse map ν=C·2^{N̂}, the conserved quantities, the privileged scale C, and full invertibility (both the individual values and δ are recovered exactly). DOI: [10.5281/zenodo.20746818](https://doi.org/10.5281/zenodo.20746818)

---

## §10 Conclusion

Starting from only a few axioms (difference; the limit of anonymity = 2; fluctuation; conservation of the product; conservation of the sum of squares), strictly distinguishing the system-interior (the unknowable) from the system-exterior (observation), and placing the objective anchors not on reality but on fixed points / conserved quantities, this report records, as a limited structure, that when the fundamental quantity is introduced as the logarithm of a squared quantity N=log₂ν, in the extremely high-symmetry, anonymity limit, the sum of squares representing scale (Σν²=R², g) can be mapped via the logarithmic connection onto a flat hyperplane (affine structure), and that this map separates scale (g) from quantization (ω, the fluctuation Δ), preserving quantization as an invariant on the flat surface. Here R is the scale of the system, not a curvature tensor.

This is not claimed as the discovery of a general fact; it is recorded as a result within a limited scope, with the posited nature of ±1/2 and the incompleteness of the counterexample search explicitly stated. Whether this structure is a special case of a known framework, and whether a generalization is possible, are left to future study.

---

*This report is a preliminary, limited record, and contains no definitive claim regarding the generality, rigor, or novelty of the result.*
