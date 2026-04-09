# Relativity of Observation in Multiple Subjective Spaces: Geometric Consequences of the Symmetries of Central Projection

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (B.Eng., Osaka University, School of Engineering Science)  
**Date:** April 2026  
**Category:** Research note (geometric investigation)  
**Preceding papers:**  
[1] Geometric Formulation of 4-Dimensional Space via Central Projection. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780).  
[2] Geometric Symmetries of Central Projection: Mathematical Foundations for the Multi-Axis Model. DOI: [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932).

---

## Abstract

In the preceding two papers [1, 2], we constructed a map from an $n$-dimensional tangent hyperplane to the hypersphere $S^n(R)$ via central projection and proved four geometric symmetries of this construction. In the present paper, as a natural consequence of Symmetry II (axis equivalence) and Symmetry IV (convertibility of subjective coordinate systems) from [2], we show that **multiple subjective spaces with different projection center axes can be constructed from the same ambient space, and the observables in each subjective space depend on the choice of the internal observer's standpoint**.

Specifically:

1. From the ambient space $\mathbb{R}^{n+1}$, different subjective spaces $H_A^+, H_B^+, \ldots$ are constructed by choosing different projection center axes
2. The internal observer in each subjective space directly observes only quantities constructed from the metric
3. Internal observers in different subjective spaces share the same ambient space while obtaining different observables
4. By Symmetry IV of [2], there exists a regular composition map between different subjective spaces via the ambient space

This is a purely geometric consequence derived directly from Symmetries II and IV of central projection and contains no physical interpretation. This paper introduces no new assumptions and argues solely from the results of [1, 2].

---

## 1. Introduction

Paper [1] defined a central projection from the $n$-dimensional tangent hyperplane $\Pi_R$ to the hypersphere $S^n(R) \subset \mathbb{R}^{n+1}$ and derived the induced metric and the Einstein tensor. Paper [2] proved four geometric symmetries of this central projection:

- **Symmetry I** (discrete stability): The map is regular whether the coordinates take zero, positive integer, or negative integer values
- **Symmetry II** (axis equivalence): The structure of the map is identical regardless of which axis of the ambient space $\mathbb{R}^{n+1}$ is chosen as the projection center
- **Symmetry III** (geodesic deviation): The radius of curvature $R$ generates geodesic deviation $|\xi|^2$
- **Symmetry IV** (convertibility of subjective coordinate systems): Two subjective coordinate systems centered on different axes are mutually convertible

In this paper, we describe the structure that follows naturally from Symmetries II and IV in particular. Specifically, we show that multiple subjective spaces constructed from the ambient space $\mathbb{R}^{n+1}$ with different projection center axes can coexist simultaneously, and that the internal observer in each subjective space obtains different observables.

All claims in this paper are purely geometric propositions and contain no physical interpretation. The scope of claims in this paper is explicitly delimited in Section 7.

---

## 2. Preliminaries: Multiple Subjective Spaces

### 2.1 Central Projection and Notation for Subjective Spaces

Following the definition in [2] Section 2.1, we consider the central projection $\Phi_A$ with the basis vector $e_A$ ($A \in \{1, \ldots, n+1\}$) in the ambient space $\mathbb{R}^{n+1}$ as the projection center axis:

$$\Phi_A: \Pi_R^{(A)} \to S^n(R), \quad Y^\mu = \frac{Rx^\mu}{l_A} \quad (\mu \neq A), \qquad Y^A = \frac{R^2}{l_A} \tag{2.1}$$

where $l_A = \sqrt{R^2 + \sum_{\mu \neq A}(x^\mu)^2}$. The image of $\Phi_A$ is the open hemisphere

$$H_A^+ = \{Y \in S^n(R) \mid Y^A > 0\} \tag{2.2}$$

which we call the **subjective space with projection center axis $A$**.

### 2.2 Simultaneous Existence of Multiple Subjective Spaces

By Symmetry II (Theorem 4.1, axis equivalence) of [2], the central projection $\Phi_A$ has the same algebraic structure regardless of which axis $A$ of the ambient space $\mathbb{R}^{n+1}$ is chosen as the projection center, and the pullback metric, curvature tensor, and Einstein tensor are all independent of the choice of $A$.

Therefore, $n+1$ subjective spaces $\{H_A^+\}_{A=1}^{n+1}$ can be simultaneously constructed from the ambient space $\mathbb{R}^{n+1}$. These $n+1$ subjective spaces share the same ambient space while each having a different axis as its projection center.

**Definition 2.1** (family of subjective spaces)  
The family of subjective spaces constructed from the ambient space $\mathbb{R}^{n+1}$ is denoted by

$$\mathcal{H} = \{H_A^+ \mid A \in \{1, 2, \ldots, n+1\}\} \tag{2.3}$$

**Proposition 2.1** (structural identity of subjective spaces)  
For any $A, B \in \{1, \ldots, n+1\}$, the subjective spaces $H_A^+$ and $H_B^+$ possess structurally identical geometric properties:

**(i)** The pullback metrics take the same algebraic form in their respective subjective coordinates  
**(ii)** The Riemann curvature tensor, Ricci tensor, Einstein tensor, and scalar curvature all take identical values  
**(iii)** The sectional curvature $K = 1/R^2$ is common to all

**Proof**  
This follows directly from Theorem 4.1 (ii) of [2]. By axis equivalence, the pullback metric takes the form

$$g_{\mu\nu}^{(A)} = \frac{R^2}{l_A^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{l_A^2}\right) \quad (\mu, \nu \neq A) \tag{2.4}$$

which is an algebraic structure independent of the choice of $A$. All curvature tensors are determined by the metric and, as a property of a space of constant curvature ([1] Theorem 4.1), are independent of the choice of axis. $\square$

### 2.3 Internal Observers and Subjective Spaces

Following the glossary in [2] Section 2, we define an **internal observer** as a hypothetical observer who performs measurements using only the metric $g_{\mu\nu}^{(A)}$ within the subjective space $H_A^+$. Hereafter, we refer to the internal observer of the subjective space $H_A^+$ as **observer $A$**.

Observer $A$ describes observations using the subjective coordinates $(x^\mu)_{\mu \neq A}$ of the subjective space $H_A^+$. By Theorem 5.2 (i) of [2], the radius of curvature $R$ is not directly measurable by observer $A$ as a quantity of dimension $[L]$ and contributes to the metric only through $R^2$. The information accessible to observer $A$ is limited to scalar invariants constructed from the metric.

---

## 3. Observation within a Subjective Space

### 3.1 Information Obtained by Observer $A$

We organize the information that the internal observer $A$ of the subjective space $H_A^+$ can obtain within the framework of [1, 2].

**Proposition 3.1** (information obtained by observer $A$)  
Observer $A$ obtains the following information:

**(i)** Their own position and motion in the subjective coordinates $(x^\mu)_{\mu \neq A}$

**(ii)** Scalar invariants constructed from the metric $g_{\mu\nu}^{(A)}$ (geodesic deviation $|\xi|^2$, scalar curvature $n(n-1)/R^2$, etc.)

Information that observer $A$ does **not** obtain:

**(iii)** The value of the projection center axis $A$ (not measurable as a quantity of dimension $[L]$ by Theorem 5.2 (i) of [2])

**(iv)** Absolute position in the ambient space $\mathbb{R}^{n+1}$ (Proposition 6.5 of [2]: structural loss of positional information)

**(v)** The existence or structure of other subjective spaces $H_B^+$ ($B \neq A$) outside the subjective space $H_A^+$ (observer $A$ has access only to the interior of their own subjective space $H_A^+$)

**Proof**  
(i)--(ii) follow directly from the results of [1, 2]. (iii) follows from Theorem 5.2 (i) of [2]. (iv) follows from Proposition 6.5 of [2]. (v) follows from the definition of an internal observer (glossary of [2] Section 2): an internal observer is a hypothetical observer who performs measurements using only the metric of their own subjective space and does not have access to other subjective spaces. $\square$

### 3.2 Role of Axes in the Subjective Space of Observer $A$

Observer $A$ performs observations on each axis of their subjective coordinates $(x^\mu)_{\mu \neq A}$ within the subjective space $H_A^+$. These $n$ axes appear to observer $A$ as "the axes of the space in which I reside."

An important point must be stated here. Observer $A$ does not know that their subjective space was "constructed with axis $A$ as the projection center." For observer $A$, the subjective coordinates $(x^\mu)_{\mu \neq A}$ are simply "the coordinates of my space," and information about which axis was the projection center axis in the ambient space cannot be obtained, by Proposition 3.1 (iv).

For observer $A$, the subjective space $H_A^+$ appears as a self-contained world. Neither the existence of the ambient space $\mathbb{R}^{n+1}$ nor the existence of other subjective spaces $H_B^+$ is directly accessible information for observer $A$.

---

## 4. Comparison of Different Subjective Spaces

### 4.1 Comparison of Observer $A$ and Observer $B$

We compare the internal observers of two different subjective spaces $H_A^+$ and $H_B^+$ ($A \neq B$): observer $A$ and observer $B$. Both reside in subjective spaces constructed from the same ambient space $\mathbb{R}^{n+1}$, but the quantities each observes play different roles.

**Proposition 4.1** (comparison of observer $A$ and observer $B$)  
The following comparison holds for observer $A$ and observer $B$:

**(i)** Observer $A$ observes the subjective coordinates $(x^\mu)_{\mu \neq A}$ and observer $B$ observes the subjective coordinates $(x'^\nu)_{\nu \neq B}$. The subjective coordinates of both observers are related by the composition map $T_{A \to B} = \Phi_B^{-1} \circ \Phi_A$ via the ambient space $\mathbb{R}^{n+1}$, by Theorem 6.1 of [2].

**(ii)** For observer $A$, axis $A$ is "outside the subjective space" (the projection center axis) and is not directly accessible. On the other hand, axis $A$ is one of observer $B$'s subjective coordinates (since $A \neq B$, axis $A$ appears as the coordinate $x'^A$ within observer $B$'s subjective space $H_B^+$).

**(iii)** Similarly, axis $B$ is "outside the subjective space" for observer $B$ and is not directly accessible, but it appears as one of observer $A$'s subjective coordinates $x^B$.

**Proof**  
(i) By Theorem 6.1 of [2], the composition $T_{A \to B} = \Phi_B^{-1} \circ \Phi_A$ of $\Phi_A$ and $\Phi_B$ is a $C^\infty$ diffeomorphism. The subjective coordinates of observer $A$ and those of observer $B$ are related through this composition map.

(ii)--(iii) By Section 2.1 of [2], the coordinates obtained by the inverse map $\Phi_A^{-1}$ are $(x^\mu)_{\mu \neq A}$, which do not include axis $A$ itself. Axis $A$ is a basis vector of the ambient space $\mathbb{R}^{n+1}$ and does not appear in observer $A$'s subjective space $H_A^+$. On the other hand, in observer $B$'s subjective space $H_B^+$, as long as $A \neq B$, axis $A$ of the ambient space appears as the subjective coordinate $x'^A$. By formula (6.5) of Proposition 6.1 of [2], the coordinate $x'^A$ in observer $B$'s subjective coordinates is related to $x^B$ in observer $A$'s subjective coordinates by $x'^A = R^2/x^B$. $\square$

**Remark** (on the notation $x'^A$)  
Here $x'^A$ denotes the coordinate component in the $A$-direction within observer $B$'s subjective coordinate system (cf. Paper [2], Eq. (6.5)). Note that this refers to the component representing the $A$-direction in $B$'s coordinate system, not axis $A$ itself.

### 4.2 Exchange of Roles of Observables

We organize the structure derived from Proposition 4.1. For axes $A$ and $B$ of the ambient space:

- **From observer $A$'s perspective**: Axis $A$ is an external axis that is not directly accessible as the "projection center axis," and axis $B$ is the coordinate $x^B$ within their own subjective space
- **From observer $B$'s perspective**: Axis $B$ is an external axis that is not directly accessible as the "projection center axis," and axis $A$ is the coordinate $x'^A$ within their own subjective space

That is, two axes $A$ and $B$ that are equivalent in the ambient space have their roles exchanged by the choice of the observer's standpoint. While Symmetry II (axis equivalence) of [2] guarantees the equivalence of axes in the ambient space, the choice of the observer's standpoint renders the roles of the axes asymmetric.

**Remark 4.1** (axis equivalence and asymmetry of roles)  
Symmetry II of [2] asserts the equivalence of axes in the ambient space. Proposition 4.1 of this section shows that this equivalence is converted into an observer-dependent asymmetry by the construction of subjective spaces (central projection). While all axes are equivalent in the ambient space, choosing a particular axis as the projection center makes that axis "external" to the observer, giving it a different role from the other axes.

**Remark 4.2** (arbitrariness of the observer's standpoint)  
By Symmetry II of [2], the geometric structure is identical regardless of which axis is chosen as the projection center. Therefore, the question "which observer is in the correct standpoint?" has no geometric meaning. The standpoints of observer $A$ and observer $B$ are completely equivalent, and the difference in their observations arises from the choice of observer itself.

---

## 5. Relativity of Observation

### 5.1 Relationship between Choice of Subjective Space and Observables

We generalize the comparison of the preceding section. For the multiple subjective spaces $\{H_A^+\}_{A=1}^{n+1}$ constructed from the ambient space $\mathbb{R}^{n+1}$, the internal observer in each subjective space obtains different observables while sharing the same ambient space.

**Theorem 5.1** (relativity of observation)  
For the family of subjective spaces $\mathcal{H} = \{H_A^+\}_{A=1}^{n+1}$ constructed from the ambient space $\mathbb{R}^{n+1}$, the following holds:

**(i)** Each subjective space $H_A^+$ possesses structurally identical geometric properties (Proposition 2.1)

**(ii)** The internal observer $A$ in each subjective space $H_A^+$ obtains quantities constructed from the metric $g_{\mu\nu}^{(A)}$ (Proposition 3.1)

**(iii)** Internal observers $A$ and $B$ in different subjective spaces $H_A^+$ and $H_B^+$ share the same ambient space while obtaining observables with different roles in their respective subjective coordinates (Proposition 4.1)

**(iv)** By Symmetry II (axis equivalence) of [2], the standpoint of every observer is geometrically equivalent, and there exists no geometric basis for designating a particular observer's standpoint as "correct"

**Proof**  
(i)--(iii) follow directly from the results of the preceding sections. (iv) follows from Theorem 4.1 of [2]: the structure of the map, metric, and curvature tensor is identical for any choice of axis $A$, and which axis is chosen as the projection center is geometrically indistinguishable. $\square$

### 5.2 Structure of the Relativity of Observation

What Theorem 5.1 demonstrates is the relativity of observables within the framework of central projection. "Observables" are determined in the following two stages:

**First stage**: The geometric structure in the ambient space $\mathbb{R}^{n+1}$. By Symmetry II of [2], all axes are equivalent, and this structure is universally preserved.

**Second stage**: The choice of the projection center axis (= the choice of the observer's standpoint). Through this choice, one of the equivalent axes in the ambient space is externalized as the "projection center axis," and the remaining axes appear as coordinates within the observer's subjective space.

Different observers have different choices at the second stage and therefore obtain different observables. However, the ambient space at the first stage is common, and by Theorem 6.1 of [2] (convertibility of subjective coordinate systems), the observables of different observers are mutually convertible via composition maps through the ambient space.

**Remark 5.1** (character of the relativity)  
The "relativity of observation" described in this paper is a relativity that depends on the choice of the observer's subjective space (= the choice of the projection center axis). This is a different kind of relativity from the physical relativity that depends on the choice of inertial or accelerated motion of an internal observer (Galilean relativity, special relativity, general relativity, etc.), and is a more fundamental relativity concerning the choice of the subjective space itself.

---

## 6. Conclusion

In this paper, as natural consequences of Symmetry II (axis equivalence) and Symmetry IV (convertibility of subjective coordinate systems) of [2], we have shown the following:

1. From the ambient space $\mathbb{R}^{n+1}$, $n+1$ subjective spaces $\{H_A^+\}$ can be simultaneously constructed by the choice of the projection center axis (Proposition 2.1)

2. The internal observer in each subjective space obtains only quantities constructed from the metric of their own subjective space and cannot directly access other subjective spaces or the structure of the ambient space (Proposition 3.1)

3. Internal observers in different subjective spaces share the same ambient space while obtaining observables with different roles in their respective subjective coordinates (Proposition 4.1)

4. By Symmetry II of [2], the standpoint of every observer is geometrically equivalent, and the difference in observables arises from the choice of the observer's subjective space (Theorem 5.1)

These results are natural consequences of Symmetries II and IV of [2] and introduce no new assumptions whatsoever. All claims in this paper are purely geometric propositions and are argued solely from the results of [1, 2].

---

## 7. Scope and Reservations

To clarify the scope of claims in this paper, the following points are reserved.

### 7.1 On Physical Interpretation

This paper deals with purely geometric propositions. Central projection, subjective spaces, symmetries---all of these are defined as objects of differential geometry, and correspondences with specific concepts from physics (spacetime, time, gravitational field, electromagnetic field, charge, energy, etc.) are not claimed in this paper. It is possible for readers to interpret the results of this paper physically, but such interpretations exceed the scope of this paper and require separate consideration.

In particular, the "relativity of observation" described in this paper is a claim about geometric structure and does not assert correspondence with any specific theory in physics (general relativity, quantum mechanics, cosmology, etc.).

### 7.2 On Inhabitants of Subjective Spaces

This paper uses expressions such as "internal observer of the subjective space," "observer $A$," and "observer $B$." These are figurative expressions referring to internal observers in subjective spaces constructed by central projection and do not assert the existence of physically real observers. This paper treats multiple subjective spaces as mathematical structures and compares the geometric properties of each.

The expression "observer" follows the definition of "internal observer" in the glossary of [2] Section 2 (a hypothetical observer who performs measurements using only the metric within the subjective space). This definition is mathematical and does not engage with philosophical concepts such as physical consciousness, perception, or subjectivity.

### 7.3 On the General Solution of the Transformation Formula

By Theorem 6.1 (convertibility of subjective coordinate systems) of [2], between the subjective coordinates of two different subjective spaces $H_A^+$ and $H_B^+$, there structurally exists a composition map via points on the sphere $S^n(R)$ (God's-eye coordinates):

$$T_{A \to B} = \Phi_B^{-1} \circ \Phi_A : \Pi_R^{(A)} \to \Pi_R^{(B)}$$

This is regular as a $C^\infty$ diffeomorphism (within the range of subjective coordinates corresponding to the common domain $U_{AB} = H_A^+ \cap H_B^+$). The structure treated in this paper---"internal observers in multiple subjective spaces observe the same ambient space in different forms"---is premised on the existence of this composition map.

However, this paper does not formulate the explicit transformation formulas of this composition map as a general solution. In Proposition 6.1 of [2], explicit formulas for the transformation between two axes are given, but the construction of a closed-form general solution for arbitrary dimensions and arbitrary combinations of axes is outside the scope of this paper and remains a topic for future research.

The claims of this paper remain at the conceptual level based on the existence of the composition map (structural proof). The concrete construction of transformation formulas is in principle derivable from the central projection formulation of [1, 2] for individual combinations of axes, but the presentation of a closed-form general solution requires separate research.

### 7.4 On Anthropic Expressions

This paper uses expressions that could be interpreted anthropically or theologically, such as "observer $A$," "observer $B$," "the observer's standpoint," and "the observer's choice." Specifically, internal observers of central projection are described in personified terms, and the information obtained by internal observers of different subjective spaces is compared.

These expressions are figurative expressions introduced to aid the reader's understanding of the geometric structure of central projection and do not assert any of the following:

- The existence of God, or the reality of God's perspective
- The reality of human beings or other entities serving as internal observers of subjective spaces
- The physical existence of observers in a specific subjective space
- The philosophical or conscious nature of the act of observation itself
- The anthropic principle (neither the weak nor the strong form)
- Correspondence with parallel universe or multiverse theories
- The measurement problem in quantum mechanics or the physics of consciousness

All claims of this paper concern the properties of geometric objects in the mathematical construction called central projection. The expressions "observer," "observation," and "standpoint" are figurative expressions introduced to aid the reader's understanding, and it is not the intention of this paper to draw philosophical, theological, or ontological conclusions from these expressions.

This paper is neither a foundational paper in physics nor a paper in philosophy; it is a paper dealing with the structure of symmetries of central projection in differential geometry.

---

## 8. Summary of Results

| Item | Content | Basis |
| --- | --- | --- |
| **Simultaneous existence of multiple subjective spaces** | $n+1$ subjective spaces $\{H_A^+\}$ can be constructed from the ambient space $\mathbb{R}^{n+1}$ | Proposition 2.1, Theorem 4.1 of [2] |
| **Information obtained by internal observers** | Observer $A$ obtains quantities constructed from the metric of the subjective space $H_A^+$ and does not access the projection center axis $A$ or other subjective spaces | Proposition 3.1, Theorem 5.2 (i) of [2] |
| **Comparison of observers** | Observers $A$ and $B$ share the same ambient space while obtaining different observables. The roles of axes in the ambient space are exchanged by the observer's standpoint | Proposition 4.1, Theorem 6.1 of [2] |
| **Relativity of observation** | Observables depend on the choice of the observer's subjective space (choice of projection center axis). By Symmetry II of [2], the standpoint of every observer is geometrically equivalent | Theorem 5.1, Symmetry II of [2] |
| **Existence of transformation (general solution outside scope)** | A $C^\infty$ regular transformation via the ambient space exists between different subjective spaces, but this paper does not present a closed-form general solution | Section 7.3, Theorem 6.1 of [2] |

---

## References

[1] Kihara, N. (2026). Geometric Formulation of 4-Dimensional Space via Central Projection. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780).

[2] Kihara, N. (2026). Geometric Symmetries of Central Projection: Mathematical Foundations for the Multi-Axis Model. DOI: [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932).

[3] do Carmo, M. P., *Riemannian Geometry*, Birkhaeuser, Boston, 1992, Chapter 4, 8.

[4] Wald, R. M., *General Relativity*, University of Chicago Press, Chicago, 1984, Appendix C, D.

[5] Wolf, J. A., *Spaces of Constant Curvature*, 6th ed., AMS Chelsea Publishing, Providence, 2011, Chapter 2.

---

*This paper is limited to the proof of purely geometric propositions. No physical, philosophical, or ontological interpretation is given. The scope of claims in this paper is explicitly delimited in Section 7.*
