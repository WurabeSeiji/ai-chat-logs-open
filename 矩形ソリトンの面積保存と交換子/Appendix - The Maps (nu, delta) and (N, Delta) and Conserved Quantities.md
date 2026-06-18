# Appendix: The forward / inverse maps (ν, δ) ↔ (N, Δ) — canonicalization by the conserved quantity, a unified map, the privileged scale, and full invertibility

**Author:** Noriaki Kihara (ORCID: 0009-0004-6753-4020)
**Date:** 18 June 2026
**DOI (this version, v1):** 10.5281/zenodo.20746818 / **Concept DOI (all versions):** 10.5281/zenodo.20746817 / **Zenodo:** https://zenodo.org/records/20746818
**Related paper:** "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in a High-Symmetry Limit," §4, §5, §6.3, §7 (Concept DOI: 10.5281/zenodo.20740841)

---

## Purpose

This appendix gives a **rigorous formulation, in a single unified form, of the forward and inverse maps** between the observable side $(\nu,\delta)_n$ and the fundamental-quantity side $(N,\Delta)_n$. The key is to **canonicalize (normalize) by the conserved quantity** on the right-hand side of the identity. As a result:

1. Every anonymous quantity (the central value $\nu$, the fluctuation $\delta$, the conserved value on the right-hand side and its fluctuation) is mapped by the **same map** $\mathrm{func}(\,\text{value}/\text{normalization}\,)$.
2. The **scale (the conserved quantity) used for normalization is the only privileged (non-anonymous) coordinate that is passed unchanged.**
3. The inverse map also has the same form, $\mathrm{func}^{-1}(\text{mapped value})\times(\text{normalization})$, so **both the individual values and the fluctuations are recovered exactly.**

This is a discrete system; the continuum is not treated. $1/2$ is a posit (parent paper §6.2). Interpretive matters such as the conjugacy of the roles of value and uncertainty are not treated (out of scope).

---

## 1. Domain (discrete; exact range)

Following §4 of the parent paper.

- Number of axes $n\in\mathbb{Z}_{\ge1}$.
- Fundamental quantity $N_i \in \mathbb{Z}$ (integers; parent paper §4: "when $\nu=2^n$, $N=n$, and $N$ ranges over all integers $\mathbb{Z}$").
- Central value $\nu_i = 2^{N_i}$ (**powers of 2**, $\nu_i\in\{2^m:m\in\mathbb{Z}\}$).
- The fluctuation $\delta_i$ is expressed as a multiplicative factor $\varphi_i=2^{\Delta_i}>0$ ($\Delta_i\in\{\pm\tfrac12\}$, the half-bit, a posit). That is, the value including the fluctuation is $\nu_i\,\varphi_i = 2^{N_i+\Delta_i}$.

Identity (conserved quantity): as in §5 of the parent paper,
$$
\prod_i \nu_i = k \quad\Longleftrightarrow\quad \sum_i N_i = K\ (=\log_2 k).
$$
$K$ (or $k$) is the **conserved quantity on the right-hand side.**

---

## 2. The principle of the unified map

Take the core of the map to be
$$
\mathrm{func} = \log_2,\qquad \mathrm{func}^{-1} = 2^{(\cdot)}.
$$
Define the **normalization (the privileged scale)** as the per-axis value of the conserved quantity,
$$
C \;=\; k^{1/n} \;=\; 2^{K/n}
$$
(so $\log_2 C = K/n$).

- **Forward map (same form for every anonymous quantity):** normalize the quantity $x$ by a reference $r$,
$$
\boxed{\ X \;=\; \mathrm{func}\!\Big(\frac{x}{r}\Big) \;=\; \log_2\frac{x}{r}.\ }
$$
- **Inverse map (same form throughout):**
$$
\boxed{\ x \;=\; r\cdot \mathrm{func}^{-1}(X) \;=\; r\cdot 2^{X} \;=\; \mathrm{func}^{-1}\big(X+\log_2 r\big).\ }
$$
(The right-most form folds the normalization inside $\mathrm{func}^{-1}$.)

- **The privileged scale $C$ is the only coordinate passed unchanged** (it is the basis of normalization and is not anonymous). Every other quantity is mapped by the same $\mathrm{func}$ as a ratio against $C$ (or a lower-level reference).

---

## 3. The map of the central value ($\nu \leftrightarrow N$, normalized form)

With reference $r=C$, normalize and map the central value:
$$
\hat N_i \;=\; \mathrm{func}\!\Big(\frac{\nu_i}{C}\Big) \;=\; \log_2\frac{\nu_i}{C} \;=\; N_i - \frac{K}{n}.
$$
The effect of normalization:
$$
\sum_i \hat N_i \;=\; \sum_i N_i - n\cdot\frac{K}{n} \;=\; K-K \;=\; 0,
\qquad \prod_i \frac{\nu_i}{C} \;=\; \frac{k}{C^n} \;=\; \frac{k}{k} \;=\; 1.
$$
That is, **canonicalization makes the conserved quantity on the right-hand side equal to $1$ (i.e., $0$ in the logarithm).** $\hat N_i$ is the **shape (direction)** of the central value, carrying the individuality of each axis.

Inverse map (same form):
$$
\nu_i \;=\; C\cdot \mathrm{func}^{-1}(\hat N_i) \;=\; C\cdot 2^{\hat N_i} \;=\; 2^{\hat N_i + K/n} \;=\; 2^{N_i}.
$$
The pair $(C,\hat N)$ is in **bijection** with $(N_i)$ (since $N_i=\hat N_i+K/n$).

---

## 4. The map of the fluctuation ($\delta \leftrightarrow \Delta$, same form)

Viewing the fluctuation as a multiplicative factor $\varphi_i=2^{\Delta_i}$ (reference $r=1$), map it by the same $\mathrm{func}$:
$$
\Delta_i \;=\; \mathrm{func}(\varphi_i) \;=\; \log_2 \varphi_i \;\in\;\{\pm\tfrac12\},
\qquad
\varphi_i \;=\; \mathrm{func}^{-1}(\Delta_i) \;=\; 2^{\Delta_i}.
$$
That is, **the fluctuation obeys exactly the same map rule as the central value** (by anonymity, ν and δ share the same forward map, and N and Δ share the same inverse map). On a conjugate pair we impose sum-zero $\sum_i \Delta_i=0$ (parent paper §6.3).

The full quantity including the fluctuation is $\nu_i\varphi_i = C\cdot 2^{\hat N_i+\Delta_i}$; conversely $\hat N_i+\Delta_i = \log_2(\nu_i\varphi_i/C)$.

---

## 5. The right-hand side (the conserved quantity) and the privileged scale

The quantity on the right-hand side of the identity fits in the same frame:

- **The conserved value $k$ → $1$ after normalization (i.e., $0$ in the logarithm).** The value itself is **passed unchanged** as the privileged scale $C=k^{1/n}$.
- **The fluctuation of the conserved value → $0$.** By the sum-zero of the conjugate pair $\sum\Delta_i=0$, the product including fluctuations is still $\prod(\nu_i\varphi_i)=k\cdot 2^{\sum\Delta_i}=k$, invariant. Hence the right-hand fluctuation maps $0\to0$.

Therefore the right-hand side (the conserved quantity) is **the basis of normalization = the privileged scale**: (i) its value is canonicalized to $1$, (ii) its fluctuation is $0$, (iii) only the scale $C$ is passed unchanged. This is what "the one extra ($+1$) dimension that is passed unchanged is the only non-anonymous coordinate" really means.

---

## 6. Full invertibility (both individual values and fluctuations are recovered exactly)

The data retained are $(C,\ \hat N,\ \Delta)$ = (privileged scale, central shape, fluctuation). By the inverse map (the same form as §2–§4),
$$
\nu_i = C\cdot 2^{\hat N_i} = 2^{N_i},\qquad \varphi_i = 2^{\Delta_i},
$$
so **the individual central values $\nu_i$ and the fluctuations $\Delta_i$ ($\delta_i$) are recovered exactly.** $(N_i)\leftrightarrow(C,\hat N_i)$ is a bijection, and so is $\Delta_i\leftrightarrow\varphi_i$. **This canonicalized map is fully invertible.**

Even if a non-integer (e.g., $C=2\sqrt2$) appears, it is **absorbed by the privileged scale $C$**, while the structure (the shape $\hat N$, the fluctuation $\Delta$) stays discrete and clean, remaining within the domain.

---

## 7. The canonicalization that discards the scale (gauge-fixing of g)

The scale $C$ is a privileged coordinate. If it is retained, the map is fully invertible as in §6. If instead it is **discarded and one works with $(\hat N,\Delta)$ only, what remains is a scale-invariant shape + fluctuation = the structure (the $\omega$ side).** This matches the g/ω separation of §7.3 of the parent paper (the scale $g$ can be removed as a gauge, while the quantization $\omega$ remains invariant). In either case the map preserves the structure $(\hat N,\Delta)$ (retain $C$ and the individual values are fully recovered; discard $C$ and the scale-invariant structure remains).

---

## 8. Numerical check (an actual sequence of operations)

Input $n=2$, $\nu=(2,4)$, fluctuation $\Delta=(+\tfrac12,-\tfrac12)$ (sum-zero).

1. **Conserved quantity / scale:** $k=\prod\nu_i=8$, $K=\sum N_i=1+2=3$, $C=k^{1/2}=2\sqrt2\approx2.828$ (the privileged scale).
2. **Forward map (central value, normalized):** $\hat N_i=\log_2(\nu_i/C)$:
 $\hat N_1=\log_2(2/2\sqrt2)=\log_2(2^{-1/2})=-\tfrac12$,
 $\hat N_2=\log_2(4/2\sqrt2)=\log_2(2^{+1/2})=+\tfrac12$. → $\hat N=(-\tfrac12,+\tfrac12)$, $\sum\hat N=0$.
3. **Forward map (fluctuation):** $\varphi_i=2^{\Delta_i}=(2^{+1/2},2^{-1/2})=(\sqrt2,\,1/\sqrt2)$.
4. **Inverse map (central value):** $\nu'_i=C\cdot2^{\hat N_i}=2\sqrt2\cdot(2^{-1/2},2^{+1/2})=(2,4)$ ← **recovered exactly.**
5. **Inverse map (fluctuation):** $\Delta'_i=\log_2\varphi_i=(+\tfrac12,-\tfrac12)=\Delta$ ← **invariant.**

Result: $\nu'=(2,4)=\nu$, $\Delta'=\Delta$. **Both the individual values and the fluctuations are conserved exactly.** The non-integer $2\sqrt2$ is contained in the privileged scale $C$, while the shape $(-\tfrac12,+\tfrac12)$ and the fluctuation $(\pm\tfrac12)$ stay clean.

---

## 9. Summary

| Quantity | Reference (normalization) | Forward map $\mathrm{func}(x/r)$ | Inverse map $r\cdot\mathrm{func}^{-1}(X)$ |
|---|---|---|---|
| Central value $\nu_i$ | Scale $C=k^{1/n}$ | $\hat N_i=\log_2(\nu_i/C)$ ($\sum\hat N=0$) | $\nu_i=C\cdot2^{\hat N_i}$ |
| Fluctuation $\varphi_i$ | $1$ | $\Delta_i=\log_2\varphi_i$ ($\sum\Delta=0$) | $\varphi_i=2^{\Delta_i}$ |
| Right-hand conserved value $k$ | $k$ | $1\to\log_2 1=0$ | passed unchanged as scale $C$ (privileged) |
| Right-hand fluctuation | — | $0$ (invariant by sum-zero) | $0$ |

- **Same map:** every anonymous quantity goes through $\mathrm{func}(x/r)=\log_2(x/r)$, with inverse $r\cdot 2^{X}$.
- **Only the privileged scale $C$** is passed unchanged — the single non-anonymous coordinate (the $+1$ dimension).
- **Fully invertible:** the individual $\nu_i$ and $\Delta_i$ are recovered exactly from $(C,\hat N,\Delta)$; non-integers are absorbed by the scale.
- **Gauge-fixing** = the canonicalization that discards the scale $C$ (retaining only the scale-invariant structure $(\hat N,\Delta)$; the g/ω separation, §7.3).

---

## 10. Limitations (notes for reading)

1. **$1/2$ is a posit** (parent paper §6.2).
2. **A discrete system:** $N\in\mathbb{Z}$, $\nu=2^N$. All non-integers are contained in the privileged scale $C$, and the structure $(\hat N,\Delta)$ is discrete.
3. **The fluctuation is treated as a multiplicative factor $\varphi=2^\Delta$** (so that it rides on the same $\mathrm{func}$ as the central value). The signed addition $\pm\tfrac12$ corresponds to $\Delta=\log_2\varphi$.
4. **The scale $C$ is a privileged (non-anonymous) coordinate.** The canonicalization that discards it corresponds to scale-invariance (gauge-fixing of g) and is consistent with §7.3.
5. **No physical identification is made.**

---

## References

- Kihara, N., "On the Logarithmic Map of Scale and the Separation of Scale and Quantization in the High-Symmetry Limit — A Limited Verification Report from a Minimal Axiom System (Revised Edition)," §4, §5, §6.3, §7, Axioms 2, 4, 5. Concept DOI: 10.5281/zenodo.20740841.
