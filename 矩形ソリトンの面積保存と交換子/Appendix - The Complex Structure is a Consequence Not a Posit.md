# Appendix: The Complex Structure is a Consequence, Not a Posit
## — How the Complex Plane Arises from Two Real Components

**Author:** Noriaki Kihara (ORCID: 0009-0004-6753-4020)
**Date:** 18 June 2026
**DOI (this version):** 10.5281/zenodo.20742985
**Concept DOI (all versions):** 10.5281/zenodo.20742984
**Zenodo:** https://zenodo.org/records/20742985
**Related paper:** "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit," §4, §7.4 (Concept DOI: [10.5281/zenodo.20740841](https://doi.org/10.5281/zenodo.20740841))

---

## Purpose

This appendix develops independently a point compressed into a few lines in §4 and the note in §7.4 of the related paper — namely that, **without positing the imaginary unit i in the axioms, the complex structure arises as a consequence of two real components and one rotation**. It adds no new claim to the parent paper; its aim is to make explicit "where i comes from," which is easily misread.

Ordinarily, quantum theory gives the complex Hilbert space as the basic framework and posits the imaginary unit i at the outset. In this system i is not placed in the axioms. What is placed in the axioms is only two independent real conserved quantities — the sum of squares (the binding of difference, g) and the product (the edge of being/non-being, ω) — and the rotation connecting them. The complex structure (the complex plane, i, the complex inner product) arises as a consequence of this real two-ness.

This appendix is not a rigorous proof but an observation making explicit the structure of the parent paper. No physical identification is made. It does not claim to "have explained why quantum theory is described by complex numbers"; it only shows a path by which the complex structure formally appears from the minimal axioms of this system.

---

## 1. Starting point: two real components, no i

The axioms of this system close over the reals. Complex numbers are not introduced. The basic objects are two independent real conserved quantities (Axioms 4 and 5 of the parent paper).

Write the two real components as ν₁, ν₂. These are anonymous, and which is called first or second is arbitrary (the anonymity of the parent paper). The two form a pair of equal real axes.

The imaginary unit i does not exist at this stage. What is shown below is that on top of the two real components ν₁, ν₂, two conserved quantities — the sum of squares and the area — stand, and that they form the same algebraic structure as the real and imaginary parts of a complex number.

---

## 2. The two conserved quantities: the sum of squares (g) and the area (ω)

For the two real components ν₁, ν₂, consider two independent conserved quantities.

**Sum of squares (the binding of difference, g):**

$$g = \nu_1^2 + \nu_2^2.$$

This binds the two on equal footing and is invariant under relabeling (ν₁↔ν₂). It is the invariant of SO(2) — the rotation preserving the sum of squares (parent paper, Axiom 5).

**Area (the edge of being/non-being, ω):** the product of fluctuations

$$\omega = \delta_1\delta_2 = \tfrac{1}{2} \quad(\text{area / symplectic form; parent paper, Axiom 4}).$$

This is an area (a quantity spanned by the two axes), related to SO(1,1) — the hyperbolic, product-preserving transformation — and an invariant preserved under rotation (parent paper, Axiom 4). Note that the product of the components themselves, ν₁ν₂, is a different bilinear quantity that depends on the rotation θ, and is not ω (δ₁δ₂) itself (§3).

By the note of the parent paper (mutual underivability of Axioms 4 and 5), g and ω cannot be derived from each other. "Whether it is" (product, ω) and "whether it differs" (sum of squares, g) are two root distinctions standing independently under anonymity.

---

## 3. The rotation: introducing cos and sin

Represent the two real components by a radius ν and an angle θ. This is merely writing the two components of a plane in polar form, not a new assumption:

$$\nu_1 = \nu\cos\theta, \qquad \nu_2 = \nu\sin\theta.$$

cos and sin are parameters of the rotation (SO(2)) preserving the sum of squares. Here, by cos²θ+sin²θ=1,

$$g = \nu_1^2 + \nu_2^2 = \nu^2(\cos^2\theta + \sin^2\theta) = \nu^2,$$

so the sum of squares g equals the square of the radius ν², and is invariant under the rotation θ (determined by the radius alone). On the other hand, the product of the individual components is

$$\nu_1\nu_2 = \nu^2\cos\theta\sin\theta = \tfrac{1}{2}\nu^2\sin 2\theta,$$

which varies with θ. This is a component-level bilinear quantity and is not ω itself. The area ω (δ₁δ₂=1/2, §2) corresponding to the imaginary part of the complex structure is a symplectic invariant preserved under rotation.

---

## 4. The arising of the complex structure

Now bind the two real components (ν₁, ν₂) into a single object. Write the bound object as z, and call its two components the "real part" and "imaginary part":

$$z = \nu_1 + i\,\nu_2 = \nu(\cos\theta + i\sin\theta) = \nu\,e^{i\theta}.$$

The i that appears here is not posited. **i is the sign of the 90-degree rotation J connecting the two real axes** (the cos↔sin of §3). The two-ness (two independent axes) gives the "plane" on which i acts, and the 90-degree rotation J on that plane corresponds to i. That **J²=−1** — a 90-degree rotation done twice is a 180-degree rotation, i.e. multiplication by −1 — becomes the rule i²=−1. i is not a new number introduced from outside but a symbol writing the rotation on the two axes algebraically. (Two independent axes alone are merely the real plane ℝ²; only by adding the rotation J does it become the complex plane ℂ.)

Under this binding, three correspondences stand:

- **Sum of squares g = ν₁²+ν₂² = |z|²** — the square of the absolute value of the complex number (real part).
- **Area ω = δ₁δ₂** — a quantity belonging to the imaginary side of the complex number (symplectic form).
- **Rotation cos↔sin** — the complex structure J (90-degree rotation, corresponding to multiplication by i).

That is, the structure written entirely in reals — two real components, the sum of squares, the area, the rotation — corresponds one-to-one with the complex number z=νe^{iθ} and its algebra (absolute value, argument, multiplication by i). The complex plane is an image that arises on top of the real two-ness (Figure 1).

![Figure 1: The complex structure arising from two real components. The point z rides on the real axes ν₁, ν₂; the radius ν=√g is invariant under rotation (dashed circle), while the individual components vary with the angle θ. The pair (ν₁, ν₂) corresponds to the complex number z=νe^{iθ}, and i appears as the sign of the 90-degree rotation J (J²=−1) connecting the two axes.](fig3_complex_from_two_real.svg)

**Correspondence with the complex inner product.** Binding the two conserved quantities into a single complex quantity gives the form

$$\langle\,\cdot\,,\cdot\,\rangle = g + i\,\omega,$$

where the real part g (sum of squares, distance, metric) and the imaginary part ω (area, symplectic form) correspond to the real and imaginary parts of the complex inner product. The g/ω separation of §7 of the parent paper states that the real and imaginary parts of this complex inner product follow separate fates under the logarithmic map (g is flattened, ω is invariant).

---

## 5. Summary of the structure

| Structure written in reals | Correspondence in complex representation |
|---|---|
| Two real components ν₁, ν₂ | real and imaginary parts of the complex number z |
| The 90-degree rotation J connecting the two axes | i (J²=−1, i.e. i²=−1) |
| Sum of squares g = ν₁²+ν₂² | |z|² (square of the absolute value, real side) |
| Area ω = δ₁δ₂ | imaginary side (symplectic form) |
| Rotation cos↔sin (SO(2)) | the complex structure J (multiplication by i, argument θ) |
| Binding of the two conserved quantities | complex inner product ⟨·,·⟩ = g + iω |

Complex numbers are not placed in the axioms. What is placed is only the two independent real conserved quantities (the sum of squares g, the area ω) and the rotation. i appears, at the stage of binding, as the symbol of the 90-degree rotation. The complex plane / complex inner product is a consequence of this real two-ness (the two-ness given by the parent paper's "limit of anonymity = 2") and the rotation.

Why base 2, why two — that is by Axiom 2 of the parent paper (the limit of anonymity = 2). That the minimum and maximum at which difference stands in the interior is two is the ground for "two real components," and also the ground for the complex structure being two-dimensional (one real part, one imaginary part). The complex number being two-dimensional is a reflection of the limit of anonymity being 2.

---

## 6. Limitations (reading notes)

This appendix makes explicit the structure of §4 and §7.4 of the parent paper and adds no new claim.

1. **This appendix does not claim to have explained "why quantum theory is described by complex numbers."** It only shows one path by which the complex structure formally appears from the minimal axioms. Isomorphism with the complexity of physical theory is a separate task.

2. **i appears as the symbol of the 90-degree rotation J (J²=−1).** This is not a claim that "i has been given physical reality," but a consequence of the algebraic notation that binds the rotation on two real components into a single complex quantity. Two independent axes alone are merely the real plane ℝ²; only by adding the rotation J does it become the complex plane ℂ. Whether one writes the real two-ness and the rotation in complex notation or in two real components plus a rotation is a choice of representation; the content of this system closes over the reals.

3. **The mutual underivability of the sum of squares g and the area ω is a premise.** As in the note of the parent paper, g and ω are two independent root distinctions, and one is not derived from the other. The complex structure appears on top of these two standing independently.

4. **No physical identification is made.** The physical names of the two real components ν₁, ν₂ are anonymous labels, and this appendix makes no identification with a particular physical quantity.

5. **The ground for two-ness depends on Axiom 2 of the parent paper.** "Why two real components (why the complex number is two-dimensional)" is attributed to the limit of anonymity = 2 (Axiom 2). The derivation of that axiom itself is outside the scope of the parent paper and this appendix.

---

## Reference

- Kihara, N. "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit — A Limited Verification Report from a Minimal Axiom System," §4, §7.4 (the complex structure as a consequence of real two-ness), Axioms 4 and 5. Concept DOI: [10.5281/zenodo.20740841](https://doi.org/10.5281/zenodo.20740841).
