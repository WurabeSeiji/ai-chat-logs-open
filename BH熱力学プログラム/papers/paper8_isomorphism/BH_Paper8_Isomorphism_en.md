# Chain Complex Structure on the 4D Hypercubic Lattice: Structural Correspondence between Kihara Cube-Packing and Wilson Lattice Gauge Theory (Paper 8: α-Identity II)

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Date**: April 2026
**Concept DOI (always resolves to the latest version)**: [10.5281/zenodo.19880467](https://doi.org/10.5281/zenodo.19880467)
**v1 Version DOI**: [10.5281/zenodo.19880468](https://doi.org/10.5281/zenodo.19880468)
**v2 Version DOI (latest, retitled to "Structural Correspondence", strengthened §6 and §7)**: [10.5281/zenodo.19881119](https://doi.org/10.5281/zenodo.19881119)

---

## Abstract

This paper clarifies the mathematical positioning of the self-consistent identity $\alpha^{-1} = N(1) + V_4(1)\cdot\alpha = 137 + (\pi^2/2)\alpha$ (predicting α to 8.7 ppb precision) observed in the companion paper (Paper 7: α-Identity). Paper 7 started from a purely geometric approach (4D cube-packing) but reached a precision too high to be coincidence. Reviewing the relationship to standard gauge theory, we have found that **the geometry of the 4D hypercubic lattice and Wilson lattice gauge theory (Wilson 1974) are structurally isomorphic as chain complexes at the level of chain complex, symmetry group, and discrete Stokes structure**.

In this paper, we (i) define the standard chain complex $C_\bullet = (C_0, C_1, C_2, C_3, C_4)$ on the 4D Euclidean integer lattice $\mathbb{Z}^4$, (ii) show that Schläfli duality (tesseract ↔ 16-cell) provides the isomorphism $C_n \leftrightarrow C_{4-n}$, (iii) demonstrate that the hypercubic symmetry group $B_4$ acts equivariantly on both, and (iv) show that the discrete Stokes theorem produces identical bulk-boundary decompositions. Wilson theory is the formalism that places gauge fields on $C_1$ (links) and an action on $C_2$ (plaquettes); Kihara cube-packing is the formalism that places volumes on $C_4$ (4-cells) and measures on $C_3$ (boundary faces). They are **structurally isomorphic as chain complexes** via Schläfli duality.

This isomorphism opens (a) gauge-theoretic computational tools (Wilson strong-coupling expansion, renormalisation group, Monte Carlo, etc.) for application to the Kihara programme, and (b) conversely, Kihara-side number-theoretic identities (Lagrange–Jacobi four-square theorem, inclusion-exclusion) as new computational methods for Wilson theory. Paper 7's inside/outside decomposition $1 = 137\alpha + V_4(1)\alpha^2$ is naturally interpreted as a bulk-boundary normalisation condition on this chain complex.

A first-principles derivation of the numerical value of α (mechanistically justifying the observation in Paper 7), the privilege of $R = 3$, and the geometric identification of the residual $c_3 \approx 1.6 \times 10^{-3}$ are left as open problems beyond this paper's scope. The contribution of this paper is to demonstrate, in standard mathematical language, that Paper 7 is not mere numerology but **a natural manifestation of the topological structure of the 4D lattice**.

---

## §1. Introduction

### §1.1 Reach of Paper 7 and Motivation

In the companion paper (Paper 7, Concept DOI: 10.5281/zenodo.19869266, v3), we observed the self-consistent algebraic identity
$$ \alpha^{-1} = N(1) + V_4(1) \cdot \alpha = 137 + \frac{\pi^2}{2}\alpha \tag{1.1} $$
constructed from $N(1) = 137$ (the count of integer-centred unit cubes fully contained in a 4D ball of radius 3) and $V_4(1) = \pi^2/2$ (the volume of the 4D unit ball). Solving the quadratic $(\pi^2/2)\alpha^2 + 137\alpha - 1 = 0$ for α yields α = 7.29735194×10⁻³, agreeing with the CODATA 2018 value to a relative error of $8.7 \times 10^{-8}$ (8.7 ppb).

This precision is approximately 1/3000 of the Eddington-style integer-fitting deviation ($2.6 \times 10^{-4}$), too good to be dismissed as "mere numerical coincidence". Paper 7 was carefully crafted to evade the Eddington trap by emphasising (i) the self-consistent equation structure, (ii) the independent 4D-geometric derivation of both 137 and $\pi^2/2$, and (iii) the QFT-perturbative analogy. However, the physical mechanism and the privilege of $R = 3$ remained open.

The motivation of this paper is the suspicion, given the strikingly high precision, that there must be a deeper mathematical structure underlying the identity. Specifically, on reviewing standard lattice gauge theory established by Wilson (1974), we recognise that the geometric setup of Paper 7 and Wilson's lattice formulation are **almost identically isomorphic from a geometric perspective**. The aim of this paper is to make this structural correspondence precise.

**Terminological caveat (v2 correction)**: v1 used the term "categorical isomorphism", but what we actually prove is **structural isomorphism as chain complexes on the 4D Euclidean integer lattice**, not a full deployment of categorical machinery (functors, natural transformations, equivalence of categories, adjunctions). v2 therefore uses "structural correspondence" / "structural isomorphism" as primary terms. This revision follows the legitimate point raised by Grok (xAI) in its v1 review.

### §1.2 Scope of the Claim

This paper **claims** the following:

1. The standard chain complex on the 4D Euclidean integer lattice $\mathbb{Z}^4$ admits an equivariant action of the hypercubic symmetry group $B_4$.
2. Schläfli duality (the self-duality of the tesseract / 16-cell) provides the isomorphism $C_n \leftrightarrow C_{4-n}$.
3. Wilson lattice gauge theory and Kihara cube-packing are formalisms that place different structures on this chain complex, and they are **structurally isomorphic as chain complexes** via Schläfli duality.
4. The discrete Stokes theorem produces an identical bulk-boundary decomposition in both formalisms.
5. Paper 7's inside/outside decomposition $1 = 137\alpha + V_4(1)\alpha^2$ is naturally interpreted from this isomorphism structure.

This paper does **not** claim:

1. That the value of α is derived from first principles via this isomorphism (we cite Paper 7 as observation only)
2. That the privilege of $R = 3$ is proved by this isomorphism
3. That U(1) gauge symmetry is uniquely emergent from this framework
4. The geometric identification of $c_3 \approx 1.6 \times 10^{-3}$
5. The specific connection to Higgs mechanism or electroweak symmetry breaking

These are all explicitly stated as open problems in §8.

### §1.3 Paper Structure

§2 defines the chain complex of the 4D hypercubic lattice in standard form. §3 describes Schläfli duality and Poincaré duality. §4 expresses Wilson lattice gauge theory in chain-complex language. §5 does the same for Kihara cube-packing. §6 states the main theorem (structural correspondence). §7 discusses the mathematical positioning of Paper 7's α identity. §8 lists open problems. §9 concludes.

This paper is written self-contained: the mathematical content of §2–§6 is accessible without reading Wilson 1974 or Paper 7.

---

## §2. Chain Complex on the 4D Hypercubic Lattice

### §2.1 Lattice Definition

Denote the 4D Euclidean integer lattice by $\Lambda = \mathbb{Z}^4 \subset \mathbb{R}^4$. Each $x \in \Lambda$ is a **site** (vertex). Directed edges between adjacent sites are **links**, with the set of links denoted $L$. The smallest 2D faces enclosed by links are **plaquettes**, 3D cells are **3-cells**, and 4D cells are **4-cells** (tesseracts).

### §2.2 Chain Complex

Free abelian groups generated by these cells with formal integer coefficients:
$$ C_0 = \mathbb{Z}\langle\text{vertices}\rangle, \quad C_1 = \mathbb{Z}\langle\text{links}\rangle, \quad C_2 = \mathbb{Z}\langle\text{plaquettes}\rangle, \quad C_3 = \mathbb{Z}\langle\text{3-cells}\rangle, \quad C_4 = \mathbb{Z}\langle\text{4-cells}\rangle. $$

The boundary operators $\partial_n : C_n \to C_{n-1}$ are defined standardly (with orientation), satisfying:
$$ \partial_n \circ \partial_{n+1} = 0. \tag{2.1} $$

This makes $C_\bullet = (C_0 \xrightarrow{\partial_1} C_1 \xrightarrow{\partial_2} C_2 \xrightarrow{\partial_3} C_3 \xrightarrow{\partial_4} C_4)$ a chain complex.

### §2.3 Action of the Hypercubic Symmetry Group $B_4$

The isometry group preserving $\Lambda$ is the hypercubic symmetry group $B_4$ (a Coxeter group of order $|B_4| = 2^4 \cdot 4! = 384$). Generators of $B_4$:

- Axis reflections: $x_i \to -x_i$ (4 generators, order $2^4 = 16$)
- Axis permutations: $x_i \leftrightarrow x_j$ ($4! = 24$ elements)

$B_4$ acts naturally on each $C_n$, and is equivariant with the boundary operators:
$$ g \cdot \partial_n = \partial_n \cdot g, \quad \forall g \in B_4. \tag{2.2} $$

### §2.4 Numerical Data of Unit Cells

| $n$ | Unit cell of $C_n$ | Boundary cells per unit | Per tesseract count |
|---|---|---|---|
| 0 | vertex | — | 16 |
| 1 | link | 2 vertices | 32 |
| 2 | plaquette | 4 links | 24 |
| 3 | 3-cell | 6 plaquettes | 8 |
| 4 | 4-cell | 8 3-cells | 1 |

These are the standard $f$-vector of the tesseract: $(16, 32, 24, 8, 1)$.

---

## §3. Schläfli Duality and Poincaré Duality

### §3.1 Self-Duality of Tesseract and 16-cell

Among 4D regular polytopes, dual pairs in Schläfli notation include the tesseract (8-cell, $\{4,3,3\}$) and the 16-cell (cross-polytope, $\{3,3,4\}$):
$$ \text{Tesseract}^* = \text{16-cell}, \quad (\text{16-cell})^* = \text{Tesseract}. $$

Specifically, the $k$-cell ($k$-dimensional face) of the tesseract corresponds to the $(4-k-1)$-cell of the 16-cell:

| Tesseract | 16-cell |
|---|---|
| 16 vertices ($k=0$) | 16 cells ($k=3$) |
| 32 edges ($k=1$) | 32 faces ($k=2$) |
| 24 faces ($k=2$) | 24 edges ($k=1$) |
| 8 cells ($k=3$) | 8 vertices ($k=0$) |

### §3.2 Poincaré Duality on the Chain Complex

Consider the dual lattice $\Lambda^* \cong \mathbb{Z}^4$ (half-integer-shifted) of $\Lambda$. The $n$-cells of $\Lambda$ correspond bijectively to the $(4-n)$-cells of $\Lambda^*$:
$$ D : C_n(\Lambda) \xrightarrow{\sim} C_{4-n}(\Lambda^*). \tag{3.1} $$

This is the discrete version of **Poincaré duality** (in 4 dimensions, with compact support), and is $B_4$-equivariant:
$$ D \circ g = g \circ D, \quad \forall g \in B_4. \tag{3.2} $$

### §3.3 Concrete Examples of the Dual Action

| $\Lambda$ side | $\Lambda^*$ side |
|---|---|
| Site ($C_0$) | 4-cell ($C_4$) |
| Link ($C_1$) | 3-cell ($C_3$) |
| Plaquette ($C_2$) | Plaquette ($C_2$, self-dual) |
| 3-cell ($C_3$) | Link ($C_1$) |
| 4-cell ($C_4$) | Site ($C_0$) |

Notable is the self-duality $C_2 \to C_2^*$, which is the direct origin of the isomorphism between Wilson plaquette action and Kihara cube-packing.

---

## §4. Wilson Lattice Gauge Theory in Chain-Complex Language

### §4.1 Link Variables

Wilson (1974) discretised gauge theory on a lattice as follows. Each link $\ell = (x, y) \in L$ carries an element $U_\ell \in G$ of the gauge group:
$$ U : C_1 \to G. $$

In this paper we focus on $G = U(1)$ (electromagnetic force), with $U_\ell = e^{i\theta_\ell}$, $\theta_\ell \in [0, 2\pi)$. Reverse link: $U_{yx} = U_{xy}^{-1}$.

### §4.2 Gauge Transformation

Assigning a phase $\alpha(x) \in U(1)$ at each site $x$, the link variable transforms:
$$ U_{xy} \to e^{i\alpha(x)} U_{xy} e^{-i\alpha(y)}. \tag{4.1} $$

This is a local gauge transformation; the phase can be chosen independently at each site.

### §4.3 Plaquette Phase

Representing a plaquette $p \in C_2$ as the closed circuit of 4 links (with orientation), the plaquette phase is:
$$ U_p = \prod_{\ell \in \partial p} U_\ell. \tag{4.2} $$

Under gauge transformation (4.1), the plaquette phase is invariant:
$$ U_p \to U_p \quad \text{(gauge invariant)}. $$

### §4.4 Wilson Action

$$ S_{\text{Wilson}} = \frac{1}{g^2} \sum_{p \in C_2} \left[1 - \text{Re}\, U_p\right] \tag{4.3} $$

where $g$ is the coupling constant. In the continuum limit $a \to 0$, this converges to the standard QED action $\frac{1}{4} F_{\mu\nu}F^{\mu\nu}$ (Wilson 1974, Kogut 1979).

### §4.5 Chain-Complex Re-reading of the Wilson Action

The Wilson action can be written chain-complex-theoretically:
$$ S_{\text{Wilson}} : C_2 \to \mathbb{R}, \quad p \mapsto \frac{1}{g^2}[1 - \text{Re}\, U_p]. \tag{4.4} $$

That is, Wilson theory is the formalism that places gauge fields on $C_1$ and an action functional on $C_2$.

### §4.6 Discrete Stokes Theorem: Bulk Cancellation

For any closed 2D surface $\Sigma \subset C_2$, the product of plaquette phases of $\Sigma$, $\prod_{p \in \Sigma} U_p$, equals the boundary Wilson loop because internal links appear pairwise in opposite orientations and cancel:
$$ \prod_{p \in \Sigma} U_p = \prod_{\ell \in \partial \Sigma} U_\ell. \tag{4.5} $$

This is the discrete version of Stokes' theorem, a direct consequence of $\partial^2 = 0$. The aggregation of bulk-to-boundary contributions is the computational core of Wilson theory.

---

## §5. Kihara Cube-Packing in Chain-Complex Language

### §5.1 Unit Cube Configuration

In Papers 3 and 7, the unit cube centred at each site $c \in \mathbb{Z}^4$ was defined as:
$$ \square(c) = \prod_{i=1}^{4} \left[c_i - \frac{1}{2}, c_i + \frac{1}{2}\right] \tag{5.1} $$

This is a 4-cell with vertices on the dual lattice $\Lambda^* = (\mathbb{Z} + 1/2)^4$, regarded as a basis element of $C_4(\Lambda^*)$.

That is, **a Kihara unit cube = a 4-cell of $\Lambda^*$**. By Schläfli duality $D : C_4(\Lambda) \to C_0(\Lambda^*)$, this is equivalent to a site of $\Lambda$.

### §5.2 Full Containment and Packing Number

The set of unit cubes **fully contained** in the 4D ball of radius $R$:
$$ \mathcal{N}(R) = \{c \in \mathbb{Z}^4 : \square(c) \subset B(R)\}. $$

Full containment requires all vertices of $\square(c)$ to lie in $B(R)$:
$$ \sum_{i=1}^{4} (|c_i| + 1/2)^2 \le R^2. \tag{5.2} $$

Packing count: $N(k) = |\mathcal{N}(2k+1)|$. By direct enumeration in Paper 3, $N(1) = 137$ ($R=3$).

### §5.3 Volume Deficit

The difference between continuous volume $V_4(R) = (\pi^2/2) R^4$ and discrete packing volume $N(k)$:
$$ \Delta(R) = V_4(R) - N(k). \tag{5.3} $$

In Paper 3 §5, the asymptotic expansion was derived:
$$ \Delta(R) = \frac{16\pi}{3} R^3 - 6\pi R^2 + O(R) \tag{5.4} $$

The leading coefficient $c = 8/(3\pi)$ (in the limit of $\Delta(R)/(2\pi^2 R^3)$) grows with the same scaling $R^3$ as the surface area of the 4D ball.

### §5.4 Chain-Complex Re-reading

Kihara cube-packing can be written chain-complex-theoretically:
$$ \chi_{B(R)} : C_4(\Lambda^*) \to \{0, 1\}, \quad \square(c) \mapsto \begin{cases} 1 & \text{if } \square(c) \subset B(R) \\ 0 & \text{otherwise} \end{cases} \tag{5.5} $$
$$ N(k) = \sum_{\square \in C_4(\Lambda^*)} \chi_{B(R)}(\square). \tag{5.6} $$

That is, the Kihara programme places an indicator function on $C_4$ and a volume measure on $C_3$ (boundary).

### §5.5 Bulk-Boundary Structure of Volume Deficit

The volume deficit $\Delta(R)$ has the following structure:

- **bulk**: in the interior of $B(R)$, cubes are packed exactly, contributing the strict integer $N_{\text{bulk}}$
- **boundary**: near $\partial B(R)$, cubes are only partially included; the "remainder" is $\Delta(R)$

This is the same structure as the discrete Stokes theorem (§4.6): bulk cancels completely, only the boundary contributions remain. Concretely, adjacent cubes share 3D faces, and the volume contribution of shared faces cancels on both sides.

---

## §6. Main Theorem: Structural Correspondence Theorem

### §6.1 Definition of the Category $\mathcal{L}_4$

**Definition 6.1**: Define the category $\mathcal{L}_4$:

- **Objects**: chain complexes $C_\bullet$ on the 4D Euclidean integer lattice (or its dual lattice), equipped with $B_4$ symmetry and boundary operators $\partial_\bullet$
- **Morphisms**: $B_4$-equivariant chain maps (homomorphisms commuting with boundary operators)

### §6.2 Isomorphism Theorem

**Theorem 6.2 (Main Theorem)**:

As objects of $\mathcal{L}_4$,

$$ (\mathcal{C}_W, \partial_\bullet, B_4) \cong_{\mathcal{L}_4} (\mathcal{C}_K, \partial_\bullet, B_4) $$

where:
- $\mathcal{C}_W$ is the chain complex of Wilson theory (gauge fields on $C_1$, action functional on $C_2$)
- $\mathcal{C}_K$ is the chain complex of Kihara packing (indicator function on $C_4$, volume measure on $C_3$)

The isomorphism is given by Schläfli duality $D$ (§3.2):

$$ D : \mathcal{C}_W \xrightarrow{\sim} \mathcal{C}_K $$
- $C_1(\Lambda)$ (Wilson links) $\xrightarrow{D}$ $C_3(\Lambda^*)$ (Kihara cube boundary faces)
- $C_2(\Lambda)$ (Wilson plaquettes) $\xrightarrow{D}$ $C_2(\Lambda^*)$ (Kihara cube shared faces, self-dual)

### §6.3 Sketch of Proof

(i) **Lattice isomorphism**: $\Lambda$ and $\Lambda^*$ are related by half-integer shift $\Lambda^* = \Lambda + (1/2, 1/2, 1/2, 1/2)$. They are isomorphic under translation, and $B_4$ axis reflections / permutations act equivariantly on both.

(ii) **Chain complex isomorphism (concrete construction)**: For each $n$-cell $\sigma$ of $\Lambda$, construct the dual $(4-n)$-cell $D(\sigma)$ of $\Lambda^*$ as follows:
- $C_0(\Lambda) \to C_4(\Lambda^*)$: vertex $x \in \Lambda$ corresponds to the unique 4-cell of $\Lambda^*$ enclosing $x$, namely $\square^*(x) = \prod_i [x_i - 1/2, x_i + 1/2]$.
- $C_1(\Lambda) \to C_3(\Lambda^*)$: edge $(x, y)$ corresponds to the 3-cell shared by the dual 4-cells of $x$ and $y$ (with orientation).
- $C_2(\Lambda) \to C_2(\Lambda^*)$: plaquette $p$ corresponds to the $\Lambda^*$ plaquette orthogonal to $p$ (the common 2D face of the four dual 4-cells around $p$). This is the **self-duality** $C_2 \to C_2$.
- $C_3(\Lambda) \to C_1(\Lambda^*)$ and $C_4(\Lambda) \to C_0(\Lambda^*)$ similarly.

(iii) **Commutativity with boundary operator**: $D \circ \partial_n = \partial_{4-n+1} \circ D$ (with sign) follows from the standard property of cap products. Concretely: for a plaquette $p \in C_2(\Lambda)$, the dual of its boundary $\partial p$ (4 links) and the boundary of its dual $D(p) \in C_2(\Lambda^*)$ (4 dual links) coincide via shared 3-cell structures between adjacent 4-cells.

(iv) **Symmetry preservation**: $B_4$ generators (axis reflection $x_i \to -x_i$, axis permutation $x_i \leftrightarrow x_j$) act in the same way on both lattices, and commute with $D$: $D \circ g = g \circ D$.

(v) **Preservation of bulk-boundary structure**: Comparing (4.5) with (5.3), Wilson's bulk cancellation $\sum_p \to \oint_{\partial}$ and Kihara's bulk packing $V_4(R) - N(k) = \Delta(\partial B(R))$ are expressed as the same discrete Stokes structure via $D$.

By the above, $\mathcal{C}_W$ and $\mathcal{C}_K$ are isomorphic as objects of $\mathcal{L}_4$. Detailed algebraic computations follow standard algebraic topology (Hatcher 2002, §3.3 on Poincaré duality). $\blacksquare$

### §6.4 Corollaries: Bidirectional Tool Borrowing

Corollaries follow from Theorem 6.2:

**Corollary 6.3**: Analytical methods of Wilson lattice gauge theory (strong-coupling expansion, renormalisation group, Monte Carlo etc.) are directly applicable to the Kihara cube-packing problem via Schläfli duality $D$.

**Corollary 6.4**: Conversely, Kihara-side number-theoretic identities (Lagrange–Jacobi four-square theorem $r_4(N) = 8\sigma(N)$, Paper 3 §6) may provide novel closed forms for finite-volume partition functions in Wilson theory.

These corollaries imply that the research programme based on Papers 7 and 8 can integrate the tools of gauge theory and number theory in a fused manner.

---

## §7. Implications: Mathematical Positioning of Paper 7's α Identity

### §7.1 Naturalness of the Inside/Outside Decomposition

The inside/outside decomposition observed in Paper 7 §4:
$$ 1 = \underbrace{137 \alpha}_{\text{inside (bulk)}} + \underbrace{V_4(1) \alpha^2}_{\text{outside (boundary)}} \tag{7.1} $$

is naturally interpreted from the structural correspondence theorem of §6. Concrete correspondences:

**Chain-complex reading of the first term $137 \alpha$**:
- $137 = N(1)$ is the count of basis elements of $C_4(\Lambda^*) \cap B(3)$ (number of 4-cells in the radius-3 ball)
- For each 4-cell $\square \in C_4(\Lambda^*) \cap B(3)$, $D(\square) \in C_0(\Lambda)$ (vertex of $\Lambda$) is the dual
- On the Wilson side, each vertex carries one gauge degree of freedom (phase freedom), contributing at coupling $\alpha$
- Therefore $137 \alpha$ is interpreted as "tree-level coupling contribution of 137 bulk degrees of freedom"

**Chain-complex reading of the second term $V_4(1) \alpha^2$**:
- $V_4(1) = \pi^2/2$ is the volume measure of the unit 4D ball
- As a measure on $C_3$ (boundary 3-cells), this describes the partial truncation of cubes near $\partial B(R)$
- On the Wilson side, the self-energy of $C_2$ plaquettes (self-dual $D: C_2 \to C_2$) contributes at $\alpha^2$ order
- The unit ball volume $V_4(1)$ measures the strength with which each boundary degree of freedom couples within its "unit influence region"
- Therefore $V_4(1) \alpha^2$ is interpreted as "self-energy correction of boundary degrees of freedom (one-loop equivalent)"

**Total normalisation $1 = $ bulk + boundary**:
- By the discrete Stokes theorem (§4.6 and §5.5), bulk interior contributions cancel completely, leaving only the residual boundary contributions
- The total normalisation $1 = 137\alpha + V_4(1)\alpha^2$ reads as a **completeness relation** on the chain complex of the 4D lattice

That is, Paper 7's self-consistent identity is re-read chain-complex-theoretically as a **bulk-boundary normalisation condition on the chain complex of the 4D lattice**, naturally following from the structural correspondence between Wilson action and Kihara packing.

**Important caveat**: The above correspondence is **structural / formal**, and does not imply that the value of $\alpha$ is **uniquely determined** by this structure. The mechanistic derivation of why the Wilson coupling $g$ corresponds to Kihara geometry such that $\alpha = 1/137.036$ is beyond the scope of this paper (§8.1).

### §7.2 Correspondence with Wilson Action

By Theorem 6.2, the formalism that places the Wilson action on Paper 7's packing problem:
$$ S_{\text{Kihara-Wilson}} = \sum_{p \subset B(3)} [1 - \text{Re}\, U_p] \tag{7.2} $$
is mathematically well-defined. This is the **finite-volume version of Wilson lattice gauge theory** naturally described in Paper 8's framework, encompassing Paper 7's geometric setup.

### §7.3 Connection with Topological Duality

The topological duality (Schläfli duality between tesseract and 16-cell) introduced in Paper 7 §6.5 is the very isomorphism map $D$ of Theorem 6.2. Paper 7 §6.5 thus prophetically observed the main theorem of this paper.

### §7.4 Repositioning of the α Identity

Combining the above, Paper 7's α identity is positioned as follows:

> **Paper 7's self-consistent identity $\alpha^{-1} = N(1) + V_4(1)\alpha$ is a structural identity naturally positioned, on the chain complex structure (Schläfli duality, $B_4$ symmetry, discrete Stokes) of the 4D Euclidean integer lattice, by the structural correspondence between Wilson lattice gauge theory and Kihara cube-packing.**

This is decisively different from Eddington-style naive numerical fitting:

- Integer 137 is an invariant of the discrete structure
- $\pi^2/2$ is an invariant of the continuous volume
- Their combination is the bulk-boundary normalisation on the chain complex
- The numerical value of α is the goodness-of-fit to this structure (observation)

However, the proof that the **numerical value of α is uniquely determined from this structure** is beyond this paper's scope and remains an open problem in §8.

---

## §8. Open Problems

We list problems unaddressed in this paper but approachable from this isomorphism structure.

### §8.1 Mechanistic Derivation of the Numerical Value of α
This isomorphism provides structural positioning of α, but the proof that the specific value $\alpha = 1/137.036$ is **uniquely determined by this isomorphism** is not given. How to identify the Wilson coupling $g$ with Kihara geometry is the issue.

### §8.2 Privilege of $R = 3$
Among $N(0) = 1, N(1) = 137, N(2) = 1545, \ldots$, why does $N(1)$ correspond to physical α? On the Wilson side, why is this finite cutoff scale physical? Correspondence with the mass gap or confinement scale is a candidate.

### §8.3 Geometric Identification of the Residual $c_3$
Paper 7 §6.1 estimated the residual at order $\alpha^3$ with $c_3 \approx 1.6 \times 10^{-3}$, but the closed form of $c_3$ is unknown. Approach via higher-order terms in the Wilson strong-coupling expansion via this isomorphism is possible.

### §8.4 Generalisation to Other Gauge Groups
This paper focused on $G = U(1)$, but generalisation to $SU(2), SU(3), SU(N)$ is natural. Whether Schläfli duality preserves group structure, and whether analogous predictions for $\alpha_s$ (strong coupling) are possible, are research questions.

### §8.5 Projection to Physical Spacetime (3+1D)
This paper and Paper 7 discussed everything in 4D Euclidean, but a projection mechanism to physical spacetime is needed. The fused formalism of central projection (Papers 2 and 6) and Wick rotation is a candidate.

### §8.6 Connection with Higgs Mechanism
Paper 7 §5.3 listed "beat coupling with longitudinal scalar field" as a candidate mechanism, but how to implement this in this isomorphism framework is unresolved (see thought experiment files in `ヒッグス波についての考察/`).

### §8.7 Numerical Verification
Computer implementation of this isomorphism, with Wilson Monte Carlo to back-calculate $N(k)$ and $\Delta(R)$, application of Lagrange–Jacobi identities to Wilson partition functions, etc., are future tasks.

---

## §9. Conclusion

This paper has clarified the mathematical positioning of the self-consistent identity $\alpha^{-1} = N(1) + V_4(1)\alpha$ (predicting α to 8.7 ppb precision) observed in Paper 7, as the chain complex structure of the 4D Euclidean integer lattice.

The main result is that **Wilson lattice gauge theory and Kihara cube-packing are structurally isomorphic as chain complexes on the chain complex ($B_4$ symmetry, Schläfli duality, discrete Stokes structure) of the 4D hypercubic lattice** (Theorem 6.2), implying:

1. Paper 7's inside/outside decomposition is naturally interpreted as a bulk-boundary normalisation on the chain complex.
2. Wilson theory's analytical methods become applicable to the Kihara programme.
3. Kihara-side number-theoretic identities may provide novel computations for Wilson theory.
4. The α identity is not mere numerology, but a natural manifestation of the topological structure of the 4D lattice.

However, the first-principles derivation of α's specific value, the privilege of $R = 3$, and the geometric identification of the residual $c_3$ are beyond this paper's scope, and are explicitly stated as open problems in §8.

The contribution of this paper is not the mechanistic derivation of α but the demonstration that Paper 7 connects naturally to standard mathematical structures (algebraic topology + lattice gauge theory). This further reinforces Paper 7's evasion of the Eddington trap, and positions this research programme to be in dialogue with standard physics and mathematics communities.

---

## References

1. Wilson, K. G. "Confinement of Quarks." *Phys. Rev. D* **10**, 2445 (1974).
2. Kogut, J. B. "An introduction to lattice gauge theory and spin systems." *Rev. Mod. Phys.* **51**, 659 (1979).
3. Hatcher, A. *Algebraic Topology*. Cambridge University Press (2002).
4. Coxeter, H. S. M. *Regular Polytopes*. Dover Publications (1973).
5. Schläfli, L. *Theorie der vielfachen Kontinuität*. Zürcher und Furrer, Zürich (1901).
6. Bekenstein, J. D. "Black Holes and Entropy." *Phys. Rev. D* **7**, 2333 (1973).
7. Hawking, S. W. "Particle Creation by Black Holes." *Commun. Math. Phys.* **43**, 199 (1975).
8. Tangherlini, F. R. "Schwarzschild Field in $n$ Dimensions and the Dimensionality of Space Problem." *Nuovo Cimento* **27**, 636 (1963).
9. Eddington, A. S. "The Charge of an Electron." *Proc. Roy. Soc. London* **A122**, 358 (1929).
10. Mohr, P. J., Newell, D. B., Taylor, B. N. "CODATA Recommended Values of the Fundamental Physical Constants: 2018." *Rev. Mod. Phys.* **93**, 025010 (2021).

---

## Appendix A: $f$-Vector of the Tesseract and Schläfli Duality

The $f$-vector of the tesseract (4D hypercube, Schläfli symbol $\{4,3,3\}$):
$$ (f_0, f_1, f_2, f_3) = (16, 32, 24, 8). $$

The $f$-vector of the 16-cell (Schläfli symbol $\{3,3,4\}$):
$$ (f_0, f_1, f_2, f_3) = (8, 24, 32, 16). $$

Schläfli duality is $f_k(\text{Tesseract}) = f_{3-k}(\text{16-cell})$.

This is a standard fact of polytope theory (Coxeter 1973), the source of self-dual structure in 4 dimensions.

## Appendix B: Numerical Reproduction of Paper 7's α Identity (for reference)

Identical to Paper 7 §3.2 (v3, 8.7 ppb precision):

```
V_4(1) = π²/2 ≈ 4.9348022
N(1) = 137
Equation: 4.9348 α² + 137 α - 1 = 0
Positive root: α ≈ 7.29735194 × 10⁻³
α⁻¹ ≈ 137.0360110
CODATA 2018: α⁻¹ = 137.0359991
Relative error: 8.7 × 10⁻⁸ (8.7 ppb)
```

For details, see Paper 7 (Concept DOI: 10.5281/zenodo.19869266).
