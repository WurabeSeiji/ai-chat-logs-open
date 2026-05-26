# Thought Experiments through the First Three Chapters of "Foundations of Quantum Theory"

## ——From the Identification Wall in Measurement to the Algebra of Observables

### Abstract

Against the background of the framework of complex Hilbert spaces, observables, and measurement developed in Chapters 1 to 3 of Akira Shimizu's *New Edition: Foundations of Quantum Theory* [1], this paper organizes five stepwise thought experiments to unify the treatment of measurement precision, the uncertainty relation, wave packets, quantum correlations, and the algebra of observables. The starting point is the "identification wall" of classical measurement theory, and the destination is the algebraic picture in which "physical quantities satisfying an uncertainty relation are different projections extracted from the same wave packet." The present paper does not contradict the mathematical system of standard quantum theory; rather, it offers one reading that supplements its conceptual outlook.

---

## Thought Experiment I: The Identification Wall in Measurement

### Setting

Suppose a physical quantity *A* takes real values. We measure the distance *L* between two points *A*₁ and *A*₂ using a measure based on the concept of distance, which returns only discrete values with a minimum spacing of Δ.

At each measurement, one of the two adjacent lattice points {*k*Δ, (*k*+1)Δ} bracketing the true value is returned with probability depending on the position *r* of the true value within the cell:

P(up) = *r*/Δ,   P(down) = 1 − *r*/Δ

(position-dependent dithering).

### Behavior of the N-fold Average

When *A*₁ = *A*₂, the single-shot distance *D* = ε₂ − ε₁ takes values {−Δ, 0, +Δ} with probabilities {1/4, 1/2, 1/4}.

- E[*D*] = 0
- Var[*D*] = Δ²/2
- Standard deviation of the N-fold average: σ_N = Δ/√(2N)

As N → ∞, σ_N → 0, and L̂_N → *L* with probability 1.

### Condition for Distinguishability

For *A*₁ and *A*₂ to be recognized as two distinct objects, the two must belong to different Δ cells. If they fall within the same cell, the measure cannot distinguish "one point or two," and the quantity "distance between two points" cannot be constituted.

Therefore the condition for the problem of measuring *L* to be well-posed is

*L* ≥ Δ.

The region *L* < Δ is not "buried in measurement error"; it is a region in which the very concept of "two points" is not defined.

### Conclusion

> When a measure based on distance is used to measure distance, measurement at scales below the measure's spacing Δ is **undefined** prior to any invocation of the uncertainty principle. Δ functions not as a resolution but as the threshold of individuation.
>
> Within the range *L* ≥ Δ, arbitrary precision is attainable by repeated measurement regardless of the absolute value of Δ.

---

## Thought Experiment II: The Locus of Fluctuation

### Setting

Consider the reverse situation. The measure has infinite precision and returns its input value unchanged. The physical quantity *A* itself, however, fluctuates about its true value *A'* according to

*A*_i = *A'* + ξ_i,   ξ_i ~ Uniform[−Δ/2, +Δ/2].

We assume *L* > Δ.

### N-fold Average

- E[L̂_N] = *A'*
- Var[*L*_i] = Δ²/12
- σ_N = Δ/√(12N) → 0

### Comparison with Thought Experiment I

Instantiating both cases with true value 7.7 and spacing Δ = 1, each yields the same empirical distribution: out of 100 measurements, roughly 30 return 7 and 70 return 8, with the average converging to 7.7.

| | Fluctuation on the measure side | Fluctuation on the quantity side |
|---|---|---|
| Single-shot output space | Discrete {*k*Δ, (*k*+1)Δ} | Continuous |
| Probability distribution | Position-dependent Bernoulli | Additive uniform |
| N → ∞ convergence target | True value | True value |
| Distinguishability from observed data | Impossible | |

### Two-sided Fluctuation

When fluctuations ξ (width δ) on the quantity side and η (width σ) on the measure side coexist independently:

- Observed value: *L*_i = *A'* + ξ_i + η_i
- σ_N = √(Var[ξ] + Var[η]) / √N → 0

The 1/√N decay persists, but separating δ and σ from a single observed series requires additional assumptions external to the observation (independent calibration of the measure, time-series decomposition, prior knowledge of distribution shapes, etc.).

### Conclusion

> Defining the observed fluctuation magnitude as Δ_obs, one cannot determine from observed data alone whether it originates from the quantization width of the measure, from the fluctuation of the quantity itself, or from a combination of both.
>
> Rather than "convergence to the objective true value," the precise statement is "convergence to the effective true value defined by the observer," and the 1/√N decay applies rigorously in this sense.

---

## Thought Experiment III: Wavelength Representation of Momentum and the Uncertainty Relation

### Setting

We extend the object of Thought Experiment II to a quantum-theoretic wave packet. For position *x*, consider the center value *x*₁ and the spread δ*x* of the wave packet. Here "center value" is merely a label of the packet distribution; we do not take the position that there exists an objective true value external to observation.

Similarly for momentum *p*, consider the center value *p*₁ and the spread δ*p*.

### Rewriting via Wavelength

From the de Broglie relation [2]:

$$p = \frac{h}{\lambda} = \hbar k, \quad k = \frac{2\pi}{\lambda}$$

the momentum fluctuation is

$$\delta p = \frac{h}{\lambda^2}\delta\lambda = \hbar \, \delta k.$$

Rewriting the uncertainty relation Δ*x*Δ*p* ≥ ℏ/2 [7] in terms of the wavenumber:

$$\Delta x \cdot \Delta k \geq \frac{1}{2}.$$

### Observation

In this form, Planck's constant ℏ disappears from both sides. What remains is a mathematical inequality of Fourier transforms concerning the product of position and wavenumber.

Planck's constant ℏ functions as a dimensional conversion factor connecting the unit of position (m) and the unit of momentum (kg·m/s).

Qualitatively:

- Pinning position precisely ⟺ narrowing the position wave packet ⟺ broadening the wavelength distribution ⟺ making momentum vague
- Pinning momentum precisely ⟺ narrowing the wavelength to one ⟺ spreading over all of space ⟺ making position vague

| Packet form | Spatial spread | Wavelength distribution |
|---|---|---|
| Delta function | δ*x* = 0 | All wavelengths |
| Rectangular pulse | Finite | sinc function |
| Pure sine wave | All of space | Single wavelength |

### Conclusion

> In the wavenumber representation, the position–momentum uncertainty relation becomes Δ*x*Δ*k* ≥ 1/2, described as a lower bound on the product of conjugate quantities in Fourier transforms. Planck's constant appears as a dimensional conversion factor.

---

## Thought Experiment IV: Quantum Correlations as a Composite Wave Packet

### Setting

Consider a state localized in two spatial regions *A* and *B* with correlations mediated by a shared conserved quantity (momentum, energy, etc.). In standard descriptions, this system is treated as a "two-particle entangled state" [3, 4, 5].

We attempt to describe this instead as "a state in which a single composite wave packet is localized in two regions in space."

### String Analogy

When a long, taut string undergoing baseline vibration has one end suddenly fixed, a soliton-like deformation appears at the opposite end. This is a result of the entire string obeying conservation laws as a single system; no information has propagated from one end to the other.

### Geometry of the Composite Wave Packet

- The composite wave packet occupies one area element in phase space
- In configuration space it has localized peaks in two regions *A* and *B*
- A local interaction occurs on the *A* side, contracting Δ*x*_A
- Due to area conservation (the lower bound of the Robertson inequality [6]), the conjugate Δ*p*_A broadens
- To satisfy the conserved quantities of the entire composite packet, the spread of the wave packet on the *B* side also undergoes geometric deformation

In this description, no information transfer from *A* to *B* is introduced. The observed correlation is a geometric consequence of a local update to one wave packet being reflected throughout via the conservation laws.

### The True Nature of the Conserved Quantity

What is conserved in the composite system is not individual wavelengths or positions, but the product of the spreads of conjugate physical quantities (an area element in phase space):

$$\Delta x \cdot \Delta p \geq \frac{\hbar}{2}, \quad \Delta A \cdot \Delta B \geq \frac{1}{2}\left|\langle [A,B] \rangle\right|.$$

The left side has the dimension of area; the inequality means a lower bound on the area.

### Conclusion

> Quantum correlations are describable as "a state in which one composite wave packet is localized across multiple regions in space." The observed correlations are geometric consequences of conservation laws acting on the entire packet, and this description does not require information transfer between spatially separated systems.

---

## Thought Experiment V: The Algebra of Observables

### Setting

We list several pairs of physical quantities that satisfy the uncertainty relation Δ*A*Δ*B* ≥ |⟨[*A*,*B*]⟩|/2:

- Position *x* ↔ momentum *p*
- Spin *S*_x ↔ *S*_y ↔ *S*_z
- Polarization (H/V ↔ D/A ↔ R/L)
- Angular momentum *L*_x ↔ *L*_y ↔ *L*_z
- Energy ↔ time

All of these share a common structure.

### Common Structure

| Pair of quantities | Shared structure | Role of distinct operators |
|---|---|---|
| *x* ↔ *p* | Area element in phase space | Direct coordinate ↔ Fourier coordinate |
| *S*_x ↔ *S*_y ↔ *S*_z | Area element in spin phase space | Projection onto distinct axes |
| Polarization basis | Area element of polarization plane | Projection onto distinct bases |
| *L*_x ↔ *L*_y ↔ *L*_z | Angular-momentum phase space | Projection onto distinct axes |
| *E* ↔ *t* | Energy–time phase space | Distinct representations |

### Unified Picture

In a complex Hilbert space:

| Formal object | Geometric interpretation |
|---|---|
| Eigenvector \|ψ⟩ | Wave packet (one area element in phase space) |
| Operator Â | Projection extracting a real value from the area element |
| Eigenvalue *a* | Result of the projection |
| Non-commutativity [Â, B̂] ≠ 0 | Two projections viewing the same area element from different directions |

### Conclusion

> Physical quantities satisfying an uncertainty relation are described as distinct projections of the same wave packet (eigenvector). The product of non-commuting projections is bounded below by the area element of the wave packet (Robertson inequality).
>
> A physical quantity is a projection of a wave packet. The wave packet exists as an area element in phase space, and the area is conserved as the product of projections.

---

## Summary

Through the five thought experiments, the following structure emerged:

1. **The identification wall** (I) — Δ is not a resolution but a threshold of individuation. *L* ≥ Δ is a prerequisite for distance measurement.
2. **Indistinguishability of the locus of fluctuation** (II) — Δ on the measure side and Δ on the quantity side cannot be distinguished from observed data.
3. **Uncertainty relation in the wavenumber representation** (III) — Δ*x*Δ*k* ≥ 1/2 is a mathematical fact of Fourier transforms. ℏ appears as a unit conversion factor.
4. **Quantum correlations as a composite wave packet** (IV) — Described as the spatial division of a single wave packet, the correlations between spatially separated systems can be explained without assuming information transfer.
5. **The algebra of observables** (V) — Physical quantities satisfying an uncertainty relation are distinct projections of the same wave packet. The area element in phase space is the conserved quantity.

These are consistent with the framework of complex Hilbert spaces, observables, and measurement introduced in Chapters 1 to 3 of Shimizu's *New Edition: Foundations of Quantum Theory*, and present one reading that provides a conceptual outlook on quantum theory. The present paper does not modify the mathematical predictions of standard quantum theory; it confines itself to offering a geometric interpretation of each concept.

---

## Acknowledgment

The present thought experiments were progressively organized through dialogue with AI. Verbatim records including the trial-and-error process are publicly available in the ai-chat-logs-open repository (`新版量子論の基礎/思考実験(6)*`, `思考実験(7)*`).

## References

[1] Shimizu, A. *New Edition: Foundations of Quantum Theory — For an Easy Understanding of Its Essence* (in Japanese), Saiensu-sha, 2003.

[2] de Broglie, L. (1924) *Recherches sur la théorie des Quanta*, doctoral thesis, Faculté des Sciences de Paris.

[3] Einstein, A., Podolsky, B., Rosen, N. (1935) "Can Quantum-Mechanical Description of Physical Reality Be Considered Complete?" *Physical Review* **47**, 777–780.

[4] Schrödinger, E. (1935) "Discussion of Probability Relations Between Separated Systems," *Mathematical Proceedings of the Cambridge Philosophical Society* **31**, 555–563.

[5] Bell, J. S. (1964) "On the Einstein-Podolsky-Rosen Paradox," *Physics Physique Физика* **1**, 195–200.

[6] Robertson, H. P. (1929) "The Uncertainty Principle," *Physical Review* **34**, 163–164.

[7] Heisenberg, W. (1927) "Über den anschaulichen Inhalt der quantentheoretischen Kinematik und Mechanik," *Zeitschrift für Physik* **43**, 172–198.
