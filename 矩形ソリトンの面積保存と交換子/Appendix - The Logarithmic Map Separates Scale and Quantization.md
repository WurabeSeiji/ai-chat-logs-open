# Appendix: The Logarithmic Map Separates Scale and Quantization
## — g is Factored Out as a Privileged Scale, ω Remains Scale-Invariant

**Author:** Noriaki Kihara (ORCID: 0009-0004-6753-4020)
**Date:** 18 June 2026
**DOI (this version, v2):** 10.5281/zenodo.20746814 (v1: 10.5281/zenodo.20742989)
**Concept DOI (all versions):** 10.5281/zenodo.20742988
**Zenodo:** https://zenodo.org/records/20746814
**Related paper:** "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit," §7.3 (Concept DOI: [10.5281/zenodo.20740841](https://doi.org/10.5281/zenodo.20740841))

---

## Purpose

This appendix develops independently the core of the related paper, §7.3 ("The heart of the map: separation of scale (g) and quantization (ω)"). It adds no new claim; its aim is to make explicit the most easily misread point — that **flattening by the logarithmic map is not "throwing away information" but separating scale (g) and quantization (ω)**.

Ordinarily, flattening in a symmetric limit is read as "structure vanishes / becomes trivial." But the logarithmic map maps both ν and δ by the same map; it **factors the scale (g) out as a privileged, separable coordinate**, while ω (Δ) remains invariant because Δ is scale-invariant. This is the separation. That the physics (quantization) does not vanish is because quantization resided from the start not in g (scale) but in ω (the fluctuation, area — a scale-invariant half-bit). This holds in any dimension and requires no full symmetry.

This appendix is not a rigorous proof but an explication of the structure of §7.3 of the parent paper. The magnitude of the fluctuation ±1/2 is a posit (parent paper §6.2). No physical identification is made. "Separation" here is meant only operationally, and includes no claim of a canonical transformation strictly preserving the symplectic form (its rigorization is future work).

---

## 1. The two quantities: g and ω

This system has two independent real conserved quantities (parent paper, Axioms 4 and 5; and the appendix "The Complex Structure is a Consequence, Not a Posit"). The two correspond to the real and imaginary parts of the complex inner product ⟨·,·⟩=g+iω.

**g (real part, distance, scale):** the sum of squares

$$g = \sum \nu_n^2 = R^2.$$

This appears to the observer as the scale R (the size, distance, metric of the system), the side of normalization (the Born norm). It is the invariant of SO(2) — the rotation preserving the sum of squares.

**ω (imaginary part, area, quantization):** the product of fluctuations

$$\omega = \delta_1\delta_2 = \tfrac{1}{2}.$$

This is an area / symplectic form, the side carrying quantization. It is related to SO(1,1) — the hyperbolic, product-preserving transformation.

By the note of the parent paper, g and ω are independent quantities, mutually underivable. One carries scale (a continuous magnitude), the other quantization (the discrete minimal unit).

---

## 2. The logarithmic map separates g

Introduce the fundamental quantity as the logarithm of a squared quantity, N=log₂ν (parent paper §4). This map **separates the scale (g) side as a privileged coordinate**.

The product (conjugate) relation νλ=k becomes, under the logarithm, the additive identity ΣNᵢ=K. **Canonicalization:** separate the conserved quantity as a privileged scale C=k^{1/n} and take the normalized shape N̂ᵢ=log₂(νᵢ/C) (ΣN̂ᵢ=0). This **holds in any dimension and requires no full symmetry** (the exact formalization of the forward / inverse maps is in the appendix "The maps (ν, δ) ↔ (N, Δ)"). The source of scale is the square (power), and the logarithm linearizes the power.

**g can be factored out as a gauge.** The scale C is the sole coordinate passed through unchanged; separating it out (or discarding it) leaves a scale-invariant shape. No structure carrying quantization remains on the g side — because g carried scale and never carried quantization.

---

## 3. The fluctuation δ is mapped the same way, but Δ is scale-invariant

The fluctuation δ, too, is mapped by the **same map** func=log₂ as the central value ν: δ → Δ=log₂(the fluctuation factor) (appendix "The maps (ν, δ) ↔ (N, Δ)," §4). That is, **the logarithmic map does act on δ** — this is a canonicalization of the same form as ν→N̂ and must not be ignored.

However, the resulting Δ=±1/2 (the half-bit) is **scale-invariant**: independent of the privileged scale C, it is the same ±1/2 at any scale. Hence, after factoring out C as a gauge, Δ remains as an invariant. **This is what the invariance of ω really is** — not "the logarithmic map does not act on ω," but "δ is mapped the same way to Δ, yet Δ is scale-invariant."

Concretely, the half-bit ±1/2 is independent of the scale C and is unchanged when C is separated. The degree of freedom of which member of the conjugate pair the half-bit is allocated to (Δ_a+Δ_b=0, sum-zero; appendix "From the Sum-Zero of Fluctuations to Uncertainty") is also preserved independently of the separation of C.

**Δ (ω) remains as the scale-invariant minimal unit of quantization.** ω = area δ₁δ₂=1/2 is this scale-invariant half-bit unit, and remains as the minimal unit of quantization even after the scale is separated (gauged away).

---

## 4. The structure of the separation

From the above, the logarithmic map separates the two quantities.

| Quantity | What it carries | Under the logarithmic map |
|---|---|---|
| g (sum of squares, distance, real part) | scale (continuous magnitude) | flattened (removable as a gauge) |
| ω (product, area, imaginary part) | quantization (discrete minimal unit) | remains invariant |

**The operational meaning of "separation":** "to separate" means that the map maps both ν and δ by the same map, while it factors the scale g (real part) out as a privileged scale C, and ω (Δ, imaginary part) remains as an invariant because Δ is scale-invariant. g and ω are the real and imaginary parts of the complex inner product g+iω; the map carries both, but the scale (real part) is separable as a gauge, while ω (imaginary part) remains scale-invariant (Figure 1).

![Figure 1: g/ω separation by the logarithmic map. The scale g in ν-space is factored out as a privileged scale C (gauge) under the logarithm. On the other hand, the quantization ω (the half-bit ±1/2) remains invariant under the map. This is not a loss of information but a separation of g and ω (no full symmetry required).](fig4_g_omega_separation.svg)

**The substance of "the physics does not vanish even when flattened":** factoring the scale (g) out as a privileged coordinate, quantization (ω) nonetheless remains because quantization resided from the start not in g (scale) but in ω (the fluctuation, area). The map separates g but does not fold ω. Hence, even when scale is factored out, quantization remains.

This separation is a consequence of quantization (ω) being a conserved quantity independent of, and distinct from, scale (g). Scale is pushed into the flat surface and gauged away, and quantization is purified as an invariant on that flat surface — this is the work of the logarithmic map.

---

## 5. Summary of the structure

For the complex inner product g+iω, the logarithmic map N=log₂ν:

- **separates / factors out the real part g (scale) as a privileged coordinate (a gauge; no full symmetry required);**
- **maps the imaginary part ω (product, quantization, Δ) by the same map, yet because Δ is scale-invariant it remains as an invariant on the flat surface.**

Flattening is not a loss of information but a separation of two quantities. Scale (g) carries a continuous magnitude and is flattenable; quantization (ω) carries the discrete minimal unit and is invariant under the map. That the two are independent (the parent paper's note: the mutual underivability of Axioms 4 and 5) is the ground for this separation.

This is the heart of the parent paper — even when flattened in a high-symmetry limit, quantization does not vanish, because quantization resides not in the side that is flattened (g, scale) but in the side that remains invariant (ω, area).

---

## 6. Limitations (reading notes)

This appendix makes explicit the structure of §7.3 of the parent paper and adds no new claim.

1. **"Separation" is meant only operationally.** It means the map flattens g and keeps ω invariant. It includes no claim of a canonical transformation strictly preserving the symplectic form. That rigorization (a proof that the map is a canonical transformation) is future work.

2. **±1/2 is a posit.** As in §6.2 of the parent paper, the minimal unit of ω, ±1/2, is a posit, not a derivation, and the separation of this appendix presupposes this posit.

3. **The difference from known frameworks is uncompared.** Whether the g/ω separation of this appendix differs from or agrees with known Kähler geometry (the structure in which a metric g and a symplectic form ω are tied by a complex structure) or conformal flatness has not been rigorously compared. As in the limitation §8.4 of the parent paper, whether this result remains a special case of known frameworks cannot be determined.

4. **Holds in any dimension and for any k.** As in the revised §8.1 of the parent paper, the flattening of the additive structure and the full invertibility of the map (canonicalization: privileged-scale separation) require no full symmetry. The g/ω separation of this appendix likewise holds in any dimension. What lies outside the scope is the extension to a curvature tensor and the linearization of the sum of squares Σν² itself.

5. **No physical identification is made.** The physical names of g (scale) and ω (quantization) are anonymous labels, and this appendix makes no identification with a particular physical quantity.

---

## Reference

- Kihara, N. "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit — A Limited Verification Report from a Minimal Axiom System," §7.2, §7.3, §8.1, §8.4, Axioms 4 and 5. Concept DOI: [10.5281/zenodo.20740841](https://doi.org/10.5281/zenodo.20740841).
