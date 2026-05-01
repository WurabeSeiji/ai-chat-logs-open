# BH Thermodynamics Programme Paper 7 Supplement: Geometric Observations on the Second-Order Correction of the α Identity

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd.
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Version**: v1
**Date**: May 1, 2026
**License**: CC BY 4.0
**Concept DOI**: [10.5281/zenodo.19933729](https://doi.org/10.5281/zenodo.19933729)
**v1 Version DOI**: [10.5281/zenodo.19933730](https://doi.org/10.5281/zenodo.19933730)
**Zenodo page**: https://zenodo.org/records/19933730
**Note article (Japanese)**: <https://note.com/kiharanoriaki/n/ne1dc24c07455>
**Note article (English)**: <https://note.com/kiharanoriaki/n/n268fa839f6c4>

---

## Character of This Document

**This is an observation paper, not a proof paper.**

For the α identity $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ in Paper 7 [BH7], we provide one geometric interpretive framework concerning the coefficient $\pi^2/2$ and the residual at the level of 8.7 ppb against observation.

This document does NOT claim:

- That $\pi^2/2$ admits no interpretation other than the one presented here
- That the residual of 8.7 ppb is fully derived under our interpretation
- That other interpretive possibilities are excluded
- That this document "proves," "strengthens," or "extends" the claims of Paper 7

What this document presents is **one reading**. The possibility that the same equations and numbers admit other readings is not denied.

The claims of Paper 7 do not depend on this document. The numerical fact observed in Paper 7 (that $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ holds to 8.7 ppb precision) does not change regardless of whether one accepts or rejects this document's interpretation.

## Abstract

In Paper 7 [BH7], a geometric identity using unit cube-packing of a 4-dimensional ball was reported:

$$\alpha^{-1} = 137 + \frac{\pi^2}{2}\alpha,$$

agreeing with the observed value $\alpha^{-1} = 137.035999177$ at relative precision $8.7 \times 10^{-8}$ (8.7 ppb). This document offers **one geometric interpretation** of this structural relation. Specifically:

1. $\pi^2/2$ admits a reading as the "normalized boundary measure $A_{S^3}(R)/(4R^3)$" of the 4D ball boundary $S^3$, rather than as the "unit ball volume $V_4(1)$" (§2).
2. Under this reading, the volume deficit $G(R)$ of cube-covering of a 4D ball admits the asymptotic expansion (§3):
   $$G(R) = \frac{16\pi}{3}R^3 - \left(\frac{3\pi^2}{4} + \frac{9\pi}{2}\right)R^2 + O(R).$$
3. The second-coefficient is analytically **negative** in this expansion (§4).
4. Numerical experiments (R = 50–300) are consistent with this expansion (§5).
5. The back-calculated correction $\delta \approx 0.001619$ from the observed residual is consistent in sign and magnitude with this geometric structure (§6).

This document does not provide a complete derivation of $\delta$. Item 5 remains an observation; the normalization rule mapping the $R^2$ deficit term to the coefficient term in the α equation is left as an open question (§7).

## §1 Introduction

### 1.1 Context of Paper 7

Paper 7 [BH7] considers unit hypercube cells centered on integer lattice points in the 4D Euclidean lattice, and establishes the self-consistent identity:

$$\alpha^{-1} = 137 + \frac{\pi^2}{2}\alpha \tag{1.1}$$

based on the count $N(1) = 137$ of cells fully contained in a 4D ball of radius $R = 3$, and the unit 4-ball volume $V_4(1) = \pi^2/2$. The positive root yields $\alpha = 7.29735194 \times 10^{-3}$, agreeing with the CODATA 2022 value $\alpha = 7.2973525693 \times 10^{-3}$ at relative precision $8.7 \times 10^{-8}$.

### 1.2 The Residual Question

Paper 7 left the origin of the residual $6.3 \times 10^{-10}$ as an open problem, identifying it as a candidate for an $\alpha^3$-order geometric correction with coefficient $c_3 \approx 1.6 \times 10^{-3}$.

### 1.3 Anticipated Concerns

The structure of Paper 7 may be subject to such concerns as:

- **Concern A**: "$\pi^2/2$ as the unit ball volume may have been chosen post-hoc to match observation."
- **Concern B**: "Without a geometric origin for the residual, the coefficient $\pi^2/2$ cannot be excluded as numerical fitting."

This document offers a **partial response** to these concerns. It is not a complete refutation but an observation that one geometric interpretation is consistent with the structure of Paper 7.

### 1.4 Document Structure

§2 presents an alternative reading of $\pi^2/2$ (as the boundary measure normalization core). §3 derives the asymptotic expansion of the volume deficit $G(R)$. §4 shows that the second coefficient is analytically negative. §5 confirms consistency with numerical experiments. §6 records the relation to the observed residual as observation. §7 explicitly acknowledges the limitations and open questions. §8 organizes implications for Paper 7.

---

## §2 An Alternative Geometric Reading of $\pi^2/2$

### 2.1 Numerical Equivalence of Two Interpretations

In Paper 7, $\pi^2/2$ was introduced as the volume of the unit 4-ball:

$$V_4(R) = \frac{\pi^2}{2}R^4, \quad V_4(1) = \frac{\pi^2}{2}. \tag{2.1}$$

However, $\pi^2/2$ also appears naturally from a different route, via the boundary measure of the 4-ball:

$$A_{S^3}(R) = \frac{d}{dR}V_4(R) = 2\pi^2 R^3. \tag{2.2}$$

Normalizing this by $4R^3$:

$$\frac{A_{S^3}(R)}{4R^3} = \frac{2\pi^2 R^3}{4R^3} = \frac{\pi^2}{2}. \tag{2.3}$$

Numerically, $V_4(1) = A_{S^3}(R)/(4R^3) = \pi^2/2$; the two are equal.

### 2.2 The Two Readings Differ in Meaning

Although the numerical value is the same, the two readings differ in their **physical/geometric meaning**:

| Reading | Meaning | Scale dependence |
|---|---|---|
| $V_4(1)$ | Volume of the unit 4-ball | Choice of "unit" present |
| $A_{S^3}(R)/(4R^3)$ | Boundary measure of 4-ball normalized by $R^3$ | **Independent of R** (same value for any R) |

The second reading does not specify any particular $R$. It yields the same value $\pi^2/2$ for any $R$. This carries the property of a **scale-invariant geometric core**.

### 2.3 The Reading Proposed by This Document

This document proposes that $\pi^2/2$ in Paper 7 (1.1) admits the reading

> **as the "normalization core of the boundary measure $A_{S^3}(R)/(4R^3)$ of the 4D ball" rather than as the "volume $V_4(1)$ of the unit 4D ball,"**

as **one possible interpretation**.

Caveats:

- This is not the "more correct" interpretation (since the numerical values are identical, mathematics cannot select between them).
- It simply provides the observation that $\pi^2/2$ appears as a scale-invariant quantity, free of any arbitrary unit choice.
- Re-reading the structure of Paper 7 under this lens leads naturally to the developments of §3 onward.

---

## §3 Asymptotic Expansion of $G(R)$

To analyze the structure of Paper 7 under this interpretive framework, we quantify the geometric mismatch between the 4D ball and the cubic lattice.

### 3.1 Definition of Volume Deficit

On the 4D Euclidean integer lattice $\mathbb{Z}^4$, place a unit hypercube cell of side 1 centered at each lattice point. Let $N_{\text{full}}(R)$ denote the number of cells fully contained in the 4D ball $B_R$ of radius $R$.

Define the volume deficit $G(R)$ as:

$$G(R) := V_4(R) - N_{\text{full}}(R) = \frac{\pi^2}{2}R^4 - N_{\text{full}}(R). \tag{3.1}$$

This is "the residual when ball volume minus fully-contained cell volume is taken," corresponding to the contribution of cells partially cut by the spherical surface.

### 3.2 Hypercube Support Function

The outward support function of a unit hypercube with respect to the unit normal $u = (u_1, u_2, u_3, u_4) \in S^3$ on the spherical surface is defined by:

$$\tau(u) = \frac{1}{2}\sum_{i=1}^{4} |u_i|. \tag{3.2}$$

This represents the "outward protrusion of the hypercube" in direction $u$.

### 3.3 Sphere Integrals

To compute the asymptotic expansion coefficients, we introduce two basic integrals on $S^3$.

#### Lemma 3.1: $\int_{S^3} \tau(u) \, d\Omega = 16\pi/3$

**Proof**: From the symmetry of $\tau(u) = (1/2)\sum_{i=1}^4 |u_i|$,

$$\int_{S^3} \tau(u) \, d\Omega = 2 \int_{S^3} |u_1| \, d\Omega.$$

Parameterize $S^3$ with $u_1 = \cos\theta$, $\theta \in [0, \pi]$. The volume element of $S^3$ is $\sin^2\theta \, d\theta \, d\Omega_{S^2}$, and the surface area of $S^2$ is $4\pi$. Thus:

$$\int_{S^3} |u_1| \, d\Omega = 4\pi \int_0^\pi |\cos\theta| \sin^2\theta \, d\theta = 8\pi \int_0^{\pi/2} \cos\theta \sin^2\theta \, d\theta = \frac{8\pi}{3}.$$

Therefore $\int_{S^3} \tau \, d\Omega = 2 \cdot 8\pi/3 = 16\pi/3$. $\square$

#### Lemma 3.2: $\int_{S^3} \tau(u)^2 \, d\Omega = \pi^2/2 + 3\pi$

**Proof**: $\tau^2 = (1/4)(\sum |u_i|)^2 = (1/4)[\sum u_i^2 + 2\sum_{i<j} |u_i||u_j|]$.

First term: $\sum_{i=1}^4 u_i^2 = 1$ (definition of unit sphere). Surface area of $S^3$ is $2\pi^2$, so

$$\int_{S^3} \sum_i u_i^2 \, d\Omega = 2\pi^2.$$

Second term: by symmetry, all 6 pairs $i \neq j$ have equal $\int |u_i||u_j| \, d\Omega$. Parameterize $S^3$ with $u_1 = \cos\theta_1$, $u_2 = \sin\theta_1 \cos\theta_2$, $u_3 = \sin\theta_1 \sin\theta_2 \cos\theta_3$, $u_4 = \sin\theta_1 \sin\theta_2 \sin\theta_3$ (volume element $\sin^2\theta_1 \sin\theta_2 \, d\theta_1 d\theta_2 d\theta_3$):

$$\int_{S^3} |u_1||u_2| \, d\Omega = 2\pi \cdot \int_0^\pi |\cos\theta_1| \sin^3\theta_1 \, d\theta_1 \cdot \int_0^\pi |\cos\theta_2| \sin\theta_2 \, d\theta_2.$$

$\int_0^\pi |\cos\theta_1| \sin^3\theta_1 \, d\theta_1 = 2 \int_0^{\pi/2} \cos\theta_1 \sin^3\theta_1 \, d\theta_1 = 1/2$.

$\int_0^\pi |\cos\theta_2| \sin\theta_2 \, d\theta_2 = 2 \int_0^{\pi/2} \cos\theta_2 \sin\theta_2 \, d\theta_2 = 1$.

So $\int_{S^3} |u_1||u_2| \, d\Omega = 2\pi \cdot 1/2 \cdot 1 = \pi$.

Therefore $\sum_{i<j} \int |u_i||u_j| \, d\Omega = 6\pi$, and

$$\int_{S^3} \tau^2 \, d\Omega = \frac{1}{4}[2\pi^2 + 2 \cdot 6\pi] = \frac{\pi^2}{2} + 3\pi. \quad \square$$

### 3.4 Asymptotic Expansion of $G(R)$

We adopt the approximation that, on the spherical boundary $\partial B_R$ in direction $u$, the hypercube cell is shaved inward by $\tau(u)$. This makes the effective radius in direction $u$ equal to $R - \tau(u)$.

$$G(R) \approx \frac{1}{4}\int_{S^3}\left[R^4 - (R - \tau(u))^4\right] d\Omega. \tag{3.3}$$

Here the factor $1/4$ is the normalization derived from $V_4(R) = \pi^2 R^4/2$ and $A_{S^3}(R) = 2\pi^2 R^3$ (the relation $dV/dR = A_{S^3}/4 \cdot$ spherical mean).

Binomial expansion:

$$R^4 - (R - \tau)^4 = 4R^3 \tau - 6R^2 \tau^2 + 4R \tau^3 - \tau^4.$$

Taking $\int_{S^3}$ and multiplying by $1/4$:

$$G(R) = R^3 \int_{S^3}\tau \, d\Omega - \frac{3}{2}R^2 \int_{S^3}\tau^2 \, d\Omega + R \int_{S^3}\tau^3 \, d\Omega - \frac{1}{4}\int_{S^3}\tau^4 \, d\Omega.$$

Substituting Lemmas 3.1 and 3.2:

$$\boxed{G(R) = \frac{16\pi}{3}R^3 - \left(\frac{3\pi^2}{4} + \frac{9\pi}{2}\right)R^2 + O(R).} \tag{3.4}$$

### 3.5 Meaning of the Main Term

The first term $(16\pi/3)R^3$ is the surface area of the 4-ball boundary $S^3$, which is $2\pi^2 R^3$, multiplied by the average boundary thickness of the hypercube. This admits the standard interpretation as "the main contribution from the continuous ball's boundary layer."

---

## §4 Observation on the Sign of the Second Term

### 4.1 Sign Analysis

The coefficient of the second term in (3.4) is:

$$-\left(\frac{3\pi^2}{4} + \frac{9\pi}{2}\right) \approx -7.402 - 14.137 = -21.539. \tag{4.1}$$

This is a negative value. In general, $\int_{S^3}\tau^2 \, d\Omega > 0$, and since the coefficient $-3/2$ stands in front, **this term is necessarily negative under this interpretive framework**.

### 4.2 Interpretation

Adopting the approximation model "the spherical boundary is shaved inward by the hypercube structure":

- First term ($R^3$): main contribution from the continuous ball's boundary
- Second term ($R^2$): second-order effect of the "shaving"

The negative sign of this second term means, **under this approximation model**, that estimating the volume deficit using only the main term overestimates relative to the continuous ball approximation.

Caveats:

- This may be read as "a fundamental property of the mismatch between ball and hypercube"
- However, this does not claim "the $\pi^2/2$ in Paper 7 admits no explanation other than this interpretation"
- The possibility that similar asymptotic expansions can be obtained from other approximation models or other geometric constructions is not excluded.

---

## §5 Numerical Verification

### 5.1 Numerical Experiment Setup

For integer $R \in \{50, 100, 200, 300\}$, compute:

1. $N_{\text{full}}(R)$: number of integer centers $n \in \mathbb{Z}^4$ satisfying $\sum_{i=1}^4 (n_i + \epsilon_i)^2 \le R^2 \, (\forall \epsilon_i \in \{\pm 1/2\})$.
2. $V_4(R) = \pi^2 R^4/2$.
3. $G(R) = V_4(R) - N_{\text{full}}(R)$.
4. $G(R)/R^3$ (value normalized to the main term).
5. $[G(R) - (16\pi/3)R^3]/R^2$ (second-coefficient after removing main term).

The implementation is in Python (numpy + integer enumeration). Code is shown in Appendix A.

### 5.2 Results

| $R$ | Filling rate $N_{\text{full}}/V_4$ | $G(R)/R^3$ | $[G(R) - (16\pi/3)R^3]/R^2$ |
|---:|---:|---:|---:|
| 50 | 0.9339 | 16.3130 | $-22.11$ |
| 100 | 0.9665 | 16.5529 | $-20.22$ |
| 200 | 0.9831 | 16.6595 | $-19.13$ |
| 300 | 0.9887 | 16.6838 | $-21.42$ |

Theoretical values:
- $G(R)/R^3 \to 16\pi/3 \approx 16.755$
- Second coefficient $\to -(3\pi^2/4 + 9\pi/2) \approx -21.539$

### 5.3 Observation

- $G(R)/R^3$ asymptotically approaches $16\pi/3$ as $R$ increases.
- The second coefficient fluctuates in the range $-20$ to $-22$, distributed around the theoretical value $-21.539$.
- The fluctuation amplitude is consistent with finite-size effects of the discrete lattice.

The numerical experiment is consistent with the asymptotic expansion (3.4). This is an observed fact. However, consistency does not show "this interpretation is correct"; it only shows "this interpretation is not numerically inconsistent."

---

## §6 Relation to the Observed Residual

### 6.1 Back-Calculation from Observed Value

Introduce a correction term $\delta$ in the equation (1.1) of Paper 7:

$$\alpha^{-1} = 137 + \left(\frac{\pi^2}{2} - \delta\right)\alpha. \tag{6.1}$$

Back-calculating from CODATA 2022 value $\alpha^{-1} = 137.035999177$:

$$\delta = \frac{\pi^2}{2} - \alpha^{-1}(\alpha^{-1} - 137) \approx 0.001619. \tag{6.2}$$

### 6.2 Consistency in Magnitude and Sign

The second term in (3.4) is an $R^2$-order geometric correction. The $\delta$ in (6.2), on the other hand, is a dimensionless coefficient that multiplies the coefficient term in the α equation directly.

To establish a direct correspondence between the two, a normalization rule that converts the $R^2$ deficit term into "the coefficient term in the α equation" is needed (discussed in §7). This rule is not formulated in this document.

However, the following observations can be recorded:

- The coefficient $-21.539$ of the second term in (3.4) is **negative**.
- The $\delta$ in (6.2) is also a quantity **subtracted** from $\pi^2/2$.
- Their signs are consistent.
- In magnitude, $\delta/(\pi^2/2) \approx 3.3 \times 10^{-4}$, falling within the natural order-of-magnitude range for discrete corrections of hypercube lattices ($\sim 10^{-3}$).

### 6.3 Limited Claims

This document claims:

- The **sign and magnitude** of $\delta$ are not inconsistent with the geometric structures of §3–§4.

This document does NOT claim:

- That the specific value $\delta = 0.001619$ has been derived from geometry.
- A complete numerical derivation.
- That $\delta$ cannot be explained under interpretive frameworks other than this one.

---

## §7 Limitations and Open Questions of This Document

### 7.1 Unformalized Mappings

A clear gap exists in this document:

**Gap G1**: A normalization rule converting the $R^2$ deficit term into the coefficient term $\delta$ of the α equation is not formulated.

The second term in (3.4) is an $R^2$-order volume correction, while $\delta$ in §6 is a dimensionless coefficient in the α equation. To connect the two, a mapping rule is needed between:

- Volume corrections at the specific scale $R = 3$ (where Paper 7 yields 137), and
- The self-consistency condition of the α equation.

This document does not provide this rule.

### 7.2 Higher-Order Corrections

The third and subsequent terms ($O(R)$, $O(1)$) of (3.4) are not computed in this document. The complete asymptotic expansion takes the form:

$$G(R) = \frac{16\pi}{3}R^3 - \left(\frac{3\pi^2}{4} + \frac{9\pi}{2}\right)R^2 + c_1 R + c_0 + o(1),$$

where $c_1, c_0$ can also be analytically computed; however, evaluation of these contributions is beyond the scope of this document.

### 7.3 Possibility of Other Interpretations

This document presents only **one geometric interpretation**. The possibility that the same numerical value emerges naturally from other geometric structures (e.g., 4D cross-polytope packing, 4D root system lattices, etc.) is not excluded.

In particular, in the framework of structural correspondence with Wilson standard lattice gauge theory shown in Paper 8 [BH8], similar correction terms may be obtained via different routes.

---

## §8 Implications for Paper 7

### 8.1 Positioning as Reinforcement

The contribution of this document is limited. This document does NOT:

- **Modify** the numerical fact of 8.7 ppb precision in Paper 7
- **Justify** the choice of $\pi^2/2$ as a coefficient (it merely presents an alternative reading of the same numerical value)
- **Solve** the residual problem (no complete derivation of $\delta$ is provided)

However, this document does provide:

- A geometric viewpoint that $\pi^2/2$ admits reading as a scale-invariant boundary measure core
- The observation that, under this viewpoint, the sign and magnitude of the residual are interpretable
- A partial response to the concerns about Paper 7 (Concerns A and B in §1.3)

### 8.2 Effect on Defense Lines

This document does **NOT depend on** the claims of Paper 7:

- The 8.7 ppb observation of Paper 7 does not change regardless of acceptance or rejection of this document's interpretation.
- Even if this document is partially rejected, Paper 7 remains unaffected.
- This document is a reinforcement attempt and does not modify the logical foundation of Paper 7.

### 8.3 Problems Posed for Subsequent Research

Problems posed by this document:

1. Resolution of Gap G1 (§7.1): formulation of the rule converting the $R^2$ deficit term to the α equation
2. Analytical computation of higher-order terms ($c_1, c_0$)
3. Derivation of an isomorphic correction term from other geometric structures (4D cross-polytope, root system lattice, etc.)
4. Re-interpretation in the framework of Wilson structural correspondence [BH8]

---

## §9 Conclusion

This document provided one geometric interpretation of $\pi^2/2$ in the α identity $\alpha^{-1} = 137 + (\pi^2/2)\alpha$ of Paper 7 [BH7]:

> **Reading $\pi^2/2$ as the normalization core $A_{S^3}(R)/(4R^3)$ of the 4D ball boundary $S^3$ measure allows the structure of Paper 7 to admit a scale-invariant interpretation.**

Under this interpretation, the volume deficit $G(R)$ of the cube-covering of the 4D ball admits the expansion:

$$G(R) = \frac{16\pi}{3}R^3 - \left(\frac{3\pi^2}{4} + \frac{9\pi}{2}\right)R^2 + O(R),$$

where the second coefficient is necessarily negative under our approximation model. This is consistent in sign and magnitude with the back-calculated $\delta \approx 0.001619$ from observation.

However, this document:

- Is an observation paper, not a proof paper
- Does not provide a complete derivation of $\delta$
- Does not exclude other interpretations
- Does not modify the claims of Paper 7

The contribution of this document is limited to: presenting one geometric viewpoint on the structure of Paper 7, and observing that this viewpoint is not numerically inconsistent.

When future researchers, perhaps a hundred years from now, provide a more complete geometric explanation of the structure of Paper 7, this document may be referenced as an intermediate observation. That is sufficient.

---

## Appendix A: Numerical Experiment Implementation

```python
import numpy as np

def count_full_cells_4d(R):
    """Count unit hypercube cells fully contained in 4D ball of radius R"""
    R_sq = R * R
    count = 0
    R_int = int(R) + 1
    for n1 in range(-R_int, R_int + 1):
        for n2 in range(-R_int, R_int + 1):
            for n3 in range(-R_int, R_int + 1):
                for n4 in range(-R_int, R_int + 1):
                    # Check if the farthest of the 8 vertices is inside the ball
                    max_d_sq = 0
                    for s1 in [-0.5, 0.5]:
                        for s2 in [-0.5, 0.5]:
                            for s3 in [-0.5, 0.5]:
                                for s4 in [-0.5, 0.5]:
                                    d_sq = (n1 + s1)**2 + (n2 + s2)**2 + \
                                           (n3 + s3)**2 + (n4 + s4)**2
                                    if d_sq > max_d_sq:
                                        max_d_sq = d_sq
                    if max_d_sq <= R_sq:
                        count += 1
    return count

def G_R(R):
    """Compute volume deficit G(R)"""
    V_4 = (np.pi**2 / 2) * R**4
    N_full = count_full_cells_4d(R)
    return V_4 - N_full

# Example
for R in [50, 100, 200, 300]:
    g = G_R(R)
    print(f"R={R}: G(R)/R^3 = {g/R**3:.4f}, "
          f"[G(R) - 16π/3·R³]/R² = {(g - (16*np.pi/3)*R**3)/R**2:.2f}")
```

Note: The above is conceptual implementation. For $R = 300$, the total number of loops becomes $4 \cdot 600^4 \approx 5 \times 10^{11}$, which is impractical. Implementation must optimize the search region near the inner edge of the ball, or use Monte Carlo estimation.

---

## Appendix B: Details of Sphere Integrals

The detailed derivations of the sphere integrals used in Lemmas 3.1 and 3.2 are completed in the main text. Recap:

### B.1 $\int_{S^3} |u_1| \, d\Omega = 8\pi/3$

The volume element of $S^3$ can be parameterized as $u_1 = \cos\theta_1$ with $\sin^2\theta_1 \cdot d\theta_1 \cdot 4\pi$ (surface area of $S^2$).

$$\int_{S^3} |u_1| \, d\Omega = 4\pi \int_0^\pi |\cos\theta| \sin^2\theta \, d\theta = 8\pi \int_0^{\pi/2} \cos\theta \sin^2\theta \, d\theta = \frac{8\pi}{3}.$$

### B.2 $\int_{S^3} |u_1 u_2| \, d\Omega = \pi$

With $u_1 = \cos\theta_1$, $u_2 = \sin\theta_1 \cos\theta_2$:

$$\int |u_1 u_2| \, d\Omega = 2\pi \int_0^\pi |\cos\theta_1| \sin^3\theta_1 \, d\theta_1 \cdot \int_0^\pi |\cos\theta_2| \sin\theta_2 \, d\theta_2 = 2\pi \cdot \frac{1}{2} \cdot 1 = \pi.$$

### B.3 $\int_{S^3} u_i^2 \, d\Omega = \pi^2/2$

By symmetry, $\int u_i^2 = (1/4) \int \sum_j u_j^2 = (1/4) \cdot \text{Vol}(S^3) = (1/4) \cdot 2\pi^2 = \pi^2/2$.

---

## References

- [BH7] Kihara, N., *A Geometric Identity for the Fine-Structure Constant: From the 4D Unit Ball Volume and its Cube-Packing Deficit*, Concept DOI: [10.5281/zenodo.19869266](https://doi.org/10.5281/zenodo.19869266) (2026).
- [BH8] Kihara, N., *Chain Complex Structure on the 4D Hypercubic Lattice: Structural Correspondence between Kihara Cube-Packing and Wilson Lattice Gauge Theory*, Concept DOI: [10.5281/zenodo.19880467](https://doi.org/10.5281/zenodo.19880467) (2026).
- [Wilson74] Wilson, K. G., *Confinement of quarks*, Phys. Rev. D 10, 2445 (1974).
- [CODATA22] CODATA Recommended Values of the Fundamental Physical Constants: 2022.

---

Author: Noriaki Kihara
WF System Co., Ltd. / ORCID: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
Concept DOI: [10.5281/zenodo.19933729](https://doi.org/10.5281/zenodo.19933729)
v1 Version DOI: [10.5281/zenodo.19933730](https://doi.org/10.5281/zenodo.19933730)
License: CC BY 4.0
