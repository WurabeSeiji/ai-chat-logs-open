---
title: "Composition of Central Projection and the Closed Form of the Composite Curvature Radius"
subtitle: "An Algebraic Formulation of High-Dimensional Reduction via One Central Projection and Commutative Cuts on the Sphere"
author: "Noriaki Kihara"
affiliation: "WF System Co., Ltd."
orcid: "0009-0004-6753-4020"
license: "CC BY 4.0"
date: "2026-05"
---

# Abstract

Dimensional reduction by central projection from an $n$-dimensional Euclidean space to a $d$-dimensional target space is often described as "the composition of successive central projections along multiple axes." This paper points out that this description is generally inaccurate and provides the correct formulation.

The main results are as follows:

- **First stage**: The central projection from the Euclidean background $\mathbb{R}^n$ onto the $(n-1)$-dimensional sphere $S^{n-1}(r_1)$ of radius $r_1$ occurs **only once**. This is the dimension-reducing map in the true sense.
- **Second stage onward**: "Axis-wise operations" applied to a point already on the sphere $S^{n-1}(r_1)$ are not central projections, but **simultaneity-section cuts on the sphere**. These operations only change the curvature radius and have no direct relation with the Euclidean background coordinates.
- **Commutativity**: Cut operations along distinct axes on the sphere are completely **commutative**.
- **Closed form of the composite curvature radius**: After successive cuts along the axis set $S = \{i_1, \ldots, i_k\} \subseteq \{1, \ldots, n\}$, the final spherical radius is

$$
r_{\rm final}^{2} \;=\; r_{1}^{2} \;-\; \sum_{j=1}^{k}\bigl(x_{i_j}^{*}\bigr)^{2}
$$

where $x_{i_j}^{*}$ is **the value of the axis $x_{i_j}$ component on the sphere immediately after the first central projection**.

- **Algebraic structure**: The set of axis-cut operations $\{\sigma_i\}$ forms an Abelian semigroup, and the composite $\sigma_S$ depends only on the axis subset $S$.

This provides a unified description of dimensional reduction from any $n$-dimensional space to a $d < n$ dimensional space as a closed operation on the combinatorics of axis selection.

Public information:

- DOI (Concept, auto-redirects to the latest version): [10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)
- DOI (v1): [10.5281/zenodo.20060729](https://doi.org/10.5281/zenodo.20060729)
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN = 6 files

---

# 1. Introduction

## 1.1 Motivation

Central projection (also known as radial projection or gnomonic projection) is a classical geometric operation that radially projects points of a Euclidean space onto a sphere or a tangent plane [1].

To realize dimensional reduction from an $n$-dimensional Euclidean space to a $d$-dimensional target space (with $d < n$), one selects $k = n - d$ axes and applies central projection successively, as is sometimes described [3]. For example, such descriptions take the following form:

> "Within an $n$-dimensional space, select axes $x_{i_1}, x_{i_2}, \ldots, x_{i_k}$ and apply central projection successively along each axis direction to reduce to an $(n-k)$-dimensional space."

However, this description has an ambiguity that allows multiple interpretations of the operation. Specifically, the fact that the first operation and the second-and-subsequent operations are **essentially different kinds of maps** is easily overlooked.

## 1.2 Purpose of this paper

The purpose of this paper is threefold:

1. To clarify the essential difference between the first and the second-and-subsequent operations in the composition of central projections.
2. To show that the correct second-and-subsequent operation is "a simultaneity-section cut on the sphere," and that this cut is commutative.
3. To derive the closed form of the composite curvature radius and establish the algebraic structure of the axis-cut operations (Abelian semigroup).

## 1.3 What this paper does not address

This paper is purely a paper on geometric algebraic structure. The following are not addressed:

- The concrete interpretation of the axis $x_i$ (its specific meaning in spacetime, observers, gauges, relativity, etc.)
- The concrete meaning of the axis selection $S$
- The physical necessity of specific dimensions $(n, d)$
- Generalization to non-Euclidean backgrounds (pseudo-Riemannian metrics, curved manifolds, etc.)

These are problems to be examined separately, taking the algebraic results of this paper as given.

---

# 2. Formulation of Central Projection and Cut Operations

## 2.1 Central projection $\pi$: from Euclidean background to sphere

We define the **central projection** as the map that sends a point $P = (x_1, \ldots, x_n)$ of $n$-dimensional Euclidean space $\mathbb{R}^n$ (with $P \neq 0$) radially to the $(n-1)$-dimensional sphere of radius $r_1$ centered at the origin $O$:

$$
S^{n-1}(r_1) \;:=\; \bigl\{X \in \mathbb{R}^{n} \;\big|\; \lVert X\rVert = r_1\bigr\}
$$

via the formula:

$$
\pi(P) \;=\; \frac{r_1}{\lVert P\rVert} P
$$

This is a well-defined map from $\mathbb{R}^n \setminus \{O\}$ to $S^{n-1}(r_1)$.

The central projection $\pi$ has the following properties:

- **(P1) Constraint on the image**: $\lVert \pi(P) \rVert = r_1$ (for any $P$).
- **(P2) Direction preservation**: $\pi(P)$ has the same radial direction as $P$.
- **(P3) Dimension of the image**: The image $\pi(P)$ is a point on the $(n-1)$-dimensional manifold $S^{n-1}(r_1)$.

## 2.2 Axis components of points on the sphere

The point after central projection $P' = \pi(P) \in S^{n-1}(r_1)$ has $n$ coordinate components in $\mathbb{R}^n$:

$$
P' \;=\; (x_1', x_2', \ldots, x_n'), \qquad x_i' \;=\; \frac{r_1}{\lVert P\rVert} x_i
$$

We call the component $x_i'$ corresponding to the axis $x_i$ the **axis $x_i$ component on the sphere**, and denote it $x_i^*$ in what follows:

$$
x_i^{*} \;:=\; \frac{r_1}{\lVert P\rVert} x_i
$$

This quantity is uniquely determined once $P$ is fixed.

## 2.3 Simultaneity-section cuts $\sigma_i$ on the sphere

Consider the operation of fixing the axis $x_i$ component of a point $X$ on the sphere $S^{n-1}(r)$ to $c$. That is, take the intersection of the hyperplane $\Sigma_i(c) := \{X \in \mathbb{R}^n \mid x_i = c\}$ with the sphere $S^{n-1}(r)$:

$$
S^{n-1}(r) \cap \Sigma_i(c) \;=\; \bigl\{X \;\big|\; \lVert X\rVert = r,\; x_i = c\bigr\}
$$

When $|c| \leq r$, this intersection is isomorphic to the $(n-2)$-dimensional sphere of radius

$$
r' \;=\; \sqrt{r^{2} - c^{2}}.
$$

We identify this with the sphere $S^{n-2}(r')$ within the $(n-1)$-dimensional space obtained by removing the axis $x_i$.

We call this operation the **cut along the axis $x_i$** and denote it $\sigma_i$:

$$
\sigma_{i}: S^{n-1}(r) \;\longrightarrow\; S^{n-2}(r'), \qquad r' = \sqrt{r^{2} - c^{2}}
$$

Here $c$ is a fixed value, and in this paper we adopt **$x_i^{*}$, the axis $x_i$ component on the sphere immediately after the first central projection**, as $c$.

## 2.4 Essential difference between central projection and cut

| Property | Central projection $\pi$ | Cut $\sigma_i$ |
|---|---|---|
| Input | Point of Euclidean background $\mathbb{R}^n$ | Point on the sphere $S^{n-1}(r)$ |
| Output | Point on the sphere $S^{n-1}(r_1)$ | Point on the sphere $S^{n-2}(r')$ |
| Dimension reduction | $n \to (n-1)$ manifold | $(n-1) \to (n-2)$ manifold |
| Radius | $r_1$ (given) | $r' = \sqrt{r^2 - c^2}$ (determined) |
| Reference to Euclidean background | Direct | None (intrinsic to the sphere) |

A central projection is "a map from the Euclidean background to a curved surface," whereas a cut is "a geometric operation on a curved surface." They are sometimes both referred to as "dimension-reduction operations," but they are **essentially different kinds of maps**.

---

# 3. Commutativity of Cut Operations

## 3.1 Main theorem: commutativity

**Theorem 1 (Commutativity of cuts).** For two distinct axes $x_i, x_j$ ($i \neq j$) on the sphere $S^{n-1}(r)$, the cut operations $\sigma_i, \sigma_j$ commute:

$$
\sigma_j \circ \sigma_i \;=\; \sigma_i \circ \sigma_j
$$

**Proof.** The action of $\sigma_i$ on a point $X = (x_1, \ldots, x_n) \in S^{n-1}(r)$ is:

- Fix the axis $x_i$ component to $c_i$ (with $c_i = x_i^{*}$).
- The remaining point in the $(n-1)$-dimensional space is $X$ with the axis $x_i$ removed.
- The remaining point lies on the sphere $S^{n-2}(r')$ where $r' = \sqrt{r^2 - c_i^2}$.

The axis $x_j$ component $c_j$ **does not change** under the cut $\sigma_i$ (since the operation removes only the axis $x_i$, the values along other axes are preserved).

Order 1: $\sigma_j \circ \sigma_i$:

- Cut $\sigma_i$ changes the radius from $r$ to $\sqrt{r^2 - c_i^2}$ and removes the axis $x_i$.
- Cut $\sigma_j$ changes the radius from $\sqrt{r^2 - c_i^2}$ to $\sqrt{r^2 - c_i^2 - c_j^2}$ and removes the axis $x_j$.
- Final radius: $\sqrt{r^2 - c_i^2 - c_j^2}$.

Order 2: $\sigma_i \circ \sigma_j$:

- Cut $\sigma_j$ changes the radius from $r$ to $\sqrt{r^2 - c_j^2}$ and removes the axis $x_j$.
- Cut $\sigma_i$ changes the radius from $\sqrt{r^2 - c_j^2}$ to $\sqrt{r^2 - c_j^2 - c_i^2}$ and removes the axis $x_i$.
- Final radius: $\sqrt{r^2 - c_j^2 - c_i^2} = \sqrt{r^2 - c_i^2 - c_j^2}$.

Both yield the same radius, and the remaining coordinate components are the same (the residue after removing both axes $x_i$ and $x_j$). Hence the maps coincide. $\Box$

## 3.2 Generalization: commutativity of $k$-axis cuts

**Corollary 1 (Commutativity of $k$-axis cuts).** For $k$ pairwise distinct axes $x_{i_1}, \ldots, x_{i_k}$, the cuts $\sigma_{i_1}, \ldots, \sigma_{i_k}$ yield the same map under any permutation:

$$
\sigma_{i_{\tau(1)}} \circ \sigma_{i_{\tau(2)}} \circ \cdots \circ \sigma_{i_{\tau(k)}} \;=\; \sigma_{i_1} \circ \sigma_{i_2} \circ \cdots \circ \sigma_{i_k}
$$

where $\tau$ is any permutation of $\{1, \ldots, k\}$.

**Proof.** By Theorem 1, any adjacent pair commutes, and the symmetric group $S_k$ is generated by adjacent transpositions. $\Box$

## 3.3 Composite cut by a set of axes

**Definition 2 (Composite cut).** For an axis set $S = \{i_1, \ldots, i_k\} \subseteq \{1, \ldots, n\}$, define the composite cut $\sigma_S$ by:

$$
\sigma_S \;:=\; \sigma_{i_1} \circ \sigma_{i_2} \circ \cdots \circ \sigma_{i_k}
$$

By Corollary 1, this is independent of the order and is determined solely by $S$.

---

# 4. Closed Form of the Composite Curvature Radius

## 4.1 Main theorem: the composite curvature radius

**Theorem 2 (Composite curvature radius).** After being projected onto the sphere $S^{n-1}(r_1)$ of radius $r_1$ by the first central projection, applying the composite cut $\sigma_S$ for the axis set $S = \{i_1, \ldots, i_k\}$ yields the final sphere of radius:

$$
\boxed{\;r_{\rm final}^{2} \;=\; r_{1}^{2} \;-\; \sum_{j=1}^{k}\bigl(x_{i_j}^{*}\bigr)^{2}\;}
$$

where $x_{i_j}^{*}$ is the axis $x_{i_j}$ component on the sphere immediately after the first central projection.

**Proof.** In the proof of Theorem 1 (or Corollary 1), each cut decreases $r^2$ by $(x_{i_j}^{*})^2$. Repeating this $k$ times yields:

$$
r_{\rm final}^{2} \;=\; r_{1}^{2} - (x_{i_1}^{*})^{2} - (x_{i_2}^{*})^{2} - \cdots - (x_{i_k}^{*})^{2}
$$

Since the result is independent of order, the closed form holds. $\Box$

## 4.2 Pythagorean structure

The closed form for the composite curvature radius takes the form of **a sum of $(x_{i_j}^{*})^2$** over the cut axes. This is nothing other than the higher-dimensional version of the Pythagorean theorem.

Specifically, for a point $P = (x_1, \ldots, x_n)$ on the sphere $S^{n-1}(r_1)$ centered at the origin,

$$
r_1^{2} \;=\; \sum_{i=1}^{n} x_i^{2}
$$

When we apply cuts along the axis set $S$, the axis $x_{i_j}$ components $(x_{i_j}^{*})$ are "subtracted" from the squared sphere radius, and the remaining $(n-k)$ axis components form the new sphere's squared radius:

$$
r_{\rm final}^{2} \;=\; \sum_{i \notin S} x_i^{2} \;=\; r_1^{2} - \sum_{i \in S} (x_i^{*})^{2}
$$

This means that "the coordinate components of a point on the sphere contribute independently as an orthogonal decomposition by the axis set."

## 4.3 Invariance of remaining coordinates

**Corollary 2 (Invariance of remaining coordinates).** Even after applying the composite cut $\sigma_S$ for the axis set $S$, the coordinate components of axes not in $S$ do not change:

$$
\bigl(\sigma_S(P)\bigr)_i \;=\; x_i \qquad \text{for all } i \notin S
$$

**Proof.** Each cut $\sigma_i$ is the operation that removes the axis $x_i$ component, leaving the components of other axes untouched. Hence, even after applying the entire composite $\sigma_S$, the components of axes not in $S$ remain invariant. $\Box$

## 4.4 Invertibility of the composite cut

**Corollary 3 (Invertibility).** The composite cut $\sigma_S$ admits an inverse map: from the set of cut-axis values $\{x_i^{*}\}_{i \in S}$ and the image point $\sigma_S(P)$ on the post-cut sphere $S^{n-|S|-1}(r_{\rm final})$, the original point $P$ on the pre-cut sphere $S^{n-1}(r_1)$ is uniquely recovered.

**Proof.** The image point $\sigma_S(P)$ has $(n - |S|)$ coordinate components, each unchanged from before the cut by Corollary 2. The cut-axis values $x_i^{*}$ ($i \in S$) are also preserved. Hence all $n$ coordinate components of the original point $P$ on $S^{n-1}(r_1)$ are recoverable. The radius is uniquely determined by Theorem 2 as $r_1^2 = r_{\rm final}^2 + \sum (x_i^{*})^2$. $\Box$

**Note**: The corollary holds even in the degenerate case $k = n - 1$ (where the final sphere is 0-dimensional, i.e., two points on a 1-dimensional number line).

---

# 5. Algebraic Structure: Abelian Semigroup

## 5.1 Composition by set operations

**Theorem 3 (Composition by set operations).** For axis sets $S, T \subseteq \{1, \ldots, n\}$ with $S \cap T = \emptyset$:

$$
\sigma_{S \cup T} \;=\; \sigma_S \circ \sigma_T
$$

**Proof.** By definition, both sides are the composition of cuts along all axes in $S \cup T$. By Corollary 1, the result is independent of order, so both sides are equal. $\Box$

## 5.2 Semigroup structure

The set of axis-cut operations $\{\sigma_S \mid S \subseteq \{1, \ldots, n\}\}$ has the following structure:

- **Operation**: composition $\circ$
- **Commutativity**: $\sigma_S \circ \sigma_T = \sigma_T \circ \sigma_S$ when $S \cap T = \emptyset$
- **Associativity**: $(\sigma_S \circ \sigma_T) \circ \sigma_U = \sigma_S \circ (\sigma_T \circ \sigma_U)$ for pairwise disjoint $S, T, U$
- **Identity**: $\sigma_\emptyset = \mathrm{id}$ (identity map)

This is the structure of an **Abelian semigroup** [2] isomorphic to subset union (when disjoint). Inverse elements do not exist (an inverse operation that increases dimension is not defined), so this is a semigroup, not a group.

## 5.3 Integration of the overall structure

The total reduction from an $n$-dimensional Euclidean space to a $d$-dimensional final space is:

1. **First stage (central projection)**: $\pi: \mathbb{R}^n \to S^{n-1}(r_1)$, exactly once.
2. **Second stage (composite cut)**: $\sigma_S: S^{n-1}(r_1) \to S^{d-1}(r_{\rm final})$, with $|S| = n - d$ axes.

The composite total map is:

$$
\Phi_{S} \;=\; \sigma_S \circ \pi \;:\; \mathbb{R}^{n} \setminus \{O\} \;\longrightarrow\; S^{d-1}(r_{\rm final})
$$

The final radius is computable in closed form as $r_{\rm final}^2 = r_1^2 - \sum_{i \in S}(x_i^{*})^2$.

---

# 6. Numerical Verification

## 6.1 Example: 7-dimensional $\to$ 5-dimensional

For the 7-dimensional Euclidean point $P = (1, 2, 3, 4, 10, 5, 3)$, we apply central projection with $r_1 = 10$, then a composite cut along 2 axes ($x_5, x_6$).

**Central projection**:

$\lVert P \rVert = \sqrt{1 + 4 + 9 + 16 + 100 + 25 + 9} = \sqrt{164} \approx 12.806$

After central projection: $P' = (10/\sqrt{164}) \cdot P \approx (0.781, 1.562, 2.343, 3.123, 7.809, 3.904, 2.343)$

Axis $x_5$ component: $x_5^{*} \approx 3.904$ (the original value 5 multiplied by $r_1 / \lVert P \rVert$)

Axis $x_6$ component: $x_6^{*} \approx 2.343$ (the original value 3 multiplied by $r_1 / \lVert P \rVert$)

**Composite cut**:

Order A ($x_5 \to x_6$):

- Radius after the $x_5$ cut: $\sqrt{10^2 - 3.904^2} = \sqrt{100 - 15.244} \approx \sqrt{84.756} \approx 9.206$
- Radius after the $x_6$ cut: $\sqrt{9.206^2 - 2.343^2} = \sqrt{84.756 - 5.488} \approx \sqrt{79.268} \approx 8.903$

Order B ($x_6 \to x_5$):

- Radius after the $x_6$ cut: $\sqrt{10^2 - 2.343^2} = \sqrt{100 - 5.488} \approx \sqrt{94.512} \approx 9.722$
- Radius after the $x_5$ cut: $\sqrt{9.722^2 - 3.904^2} = \sqrt{94.512 - 15.244} \approx \sqrt{79.268} \approx 8.903$

Both orders yield the same final radius $r_{\rm final} \approx 8.903$.

**Verification by closed form**:

$$
r_{\rm final}^{2} \;=\; r_{1}^{2} - (x_5^{*})^{2} - (x_6^{*})^{2} \;=\; 100 - 15.244 - 5.488 \;=\; 79.268
$$

$$
r_{\rm final} \;=\; \sqrt{79.268} \;\approx\; 8.903
$$

Order A, Order B, and the closed form all agree.

## 6.2 Verification of invariance of remaining coordinates

The 5-dimensional components of the final point after the composite cut:

$$
(0.781, 1.562, 2.343, 3.123, 7.809)
$$

These are the first 5 components of $P'$ immediately after central projection, unchanged by the cut operations. This verifies Corollary 2.

## 6.3 Verification of invertibility

Using the 5-dimensional components $(0.781, 1.562, 2.343, 3.123, 7.809)$ and the cut-axis values $x_5^{*} = 3.904$, $x_6^{*} = 2.343$, we recover the point on the pre-cut sphere $S^6(r_1 = 10)$:

$$
P' \;=\; (0.781, 1.562, 2.343, 3.123, 7.809, 3.904, 2.343)
$$

This is the original $P'$, and $\lVert P' \rVert = 10 = r_1$ is also recovered. This verifies Corollary 3.

---

# 7. Compatibility with Existing Descriptions

## 7.1 Interpretation of "successive central projections along multiple axes"

When the reduction from $n$-dimensional to $d$-dimensional is described as "successive central projections along multiple axes," the formulation in this paper interprets it as follows:

- "First central projection": The actual central projection $\pi$ from the Euclidean background to the sphere.
- "Second-and-subsequent central projections": Cut operations $\sigma_i$ on the sphere.

It is possible to refer to both as "central projection" without distinction, but in that case:

- First: $\mathbb{R}^n \to S^{n-1}(r_1)$, true dimension reduction.
- Second-and-subsequent: $S^m(r) \to S^{m-1}(r')$, dimension reduction between spheres.

One should note that this two-stage structure is implicitly assumed.

## 7.2 Note on order dependence

If "projection along axis $x_i$" in the Euclidean background is naively interpreted as "dividing all coordinates by $x_i$," then composition along multiple axes can become order-dependent. Specifically, treating a point on the sphere (after projection along axis $x_i$) again as a point of the Euclidean background and trying to project along axis $x_j$ leads to ambiguity in the meaning of the operation.

In the framework of this paper, after the first central projection, points are treated as "points on the sphere," and second-and-subsequent operations apply the cut $\sigma_j$ on the sphere. With this distinction, composition becomes completely commutative.

---

# 8. Conclusion

This paper has established:

- Dimensional reduction by central projection essentially consists of **one true central projection plus $k$ cut operations on the sphere**.
- The cut operations on the sphere are completely **commutative**, and a closed composition operation $\sigma_S$ for axis sets is defined.
- The composite curvature radius is given in closed form

$$
r_{\rm final}^{2} \;=\; r_{1}^{2} - \sum_{i \in S}(x_i^{*})^{2}
$$

with a Pythagorean orthogonal decomposition structure.

- The set of cut operations forms an Abelian semigroup.
- Remaining coordinates are invariant under cut operations.
- The composite cut is **invertible**: from the cut-axis values and the image point, the original point on the sphere can be uniquely recovered.

These algebraic structures function as a foundational language in various geometric and physical frameworks dealing with reduction from high to low dimensions.

---

# 9. Open Problems

This paper is purely a paper on geometric algebraic structure. The following problems lie outside its scope and are to be examined separately:

- **Concrete interpretation**: The specific meaning of axes $x_i$, the necessity of particular $(n, d)$, and relations to observers.
- **Algebraic structure of the inverse operation**: Invertibility (Corollary 3) exists at the level of maps; whether to extend this into a group structure on the algebraic structure of cut operations, or remain at the semigroup level, is a matter of choice.
- **Non-Euclidean backgrounds**: Analogous structures for pseudo-Riemannian metrics or curved background manifolds.
- **Group-theoretic extensions**: Group-theoretic structures of transformations that change the axis selection (continuous choices of axes), and correspondences with classical transformation groups.
- **Continuous cut values**: This paper fixes the cut value $c$ to the spherical component $x_i^{*}$ immediately after the first central projection; an algebraic structure for cuts with more general $c$ is left open.

---

# References

[1] Snyder, J.P. (1987). *Map Projections—A Working Manual*. U.S. Geological Survey Professional Paper 1395, U.S. Government Printing Office, Washington, D.C. (Modern standard reference for map projections including gnomonic projection. Cited in §1.1, §2.1 as the source of the classical definition of central projection.)

[2] Howie, J.M. (1995). *Fundamentals of Semigroup Theory*. London Mathematical Society Monographs, New Series, Vol. 12, Oxford University Press. (Standard reference for semigroup theory. Cited in §5 as the source of the Abelian semigroup structure.)

[3] Kihara, N. (2026). *A Geometric Formulation of 4-Dimensional Space via Central Projection* (in Japanese). Zenodo. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780). (Foundational paper of the central projection framework. Cited in §1.1 as one of the sources of the informal "successive central projections along multiple axes" description that this paper aims to correct. The algebraic results of the present paper hold independently of the application contexts in the cited paper.)

---

# Author Information

- **Author**: Noriaki Kihara
- **Affiliation**: WF System Co., Ltd. / Graduated from the School of Engineering Science, Osaka University (1983)
- **ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
- **License**: CC BY 4.0

---

**Note**: This paper provides an independent algebraic foundation for dimensional reduction by central projection. Its algebraic results do not depend on specific dimensions, axis selections, or physical interpretations, and provide a universal foundation for any application that requires them. Specific applications are developed in separate papers.
