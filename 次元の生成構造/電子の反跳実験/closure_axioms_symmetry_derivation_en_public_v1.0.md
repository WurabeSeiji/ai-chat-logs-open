# Symmetry Generation from Zero Closure, Finite Order, and Self-Consistent Geometry — The Single External Parameter $N$ and the Remaining Tasks of Generalization and Dynamics

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** August 20, 2026<br>
**Version DOI:** 10.5281/zenodo.22028073<br>
**Concept DOI:** 10.5281/zenodo.22028072<br>
**Series:** "Generative Structure of Dimensions" series — Symmetry Derivation from Closure Axioms, public version v1.0<br>
**License:** CC BY 4.0

> **Subject**: This paper organizes how far geometry, symmetry, statistics, and readout structure can be derived from an axiomatic framework that, apart from the parts not yet derived or generalized, contains no explicit continuous tuning parameters and uses only a single discrete integer $N$ as its external parameter. While the five complex degrees of freedom are derived rigorously, the self-consistent selection rule among readout sectors including the $3+2$ decomposition, the generalization of some structures, and the derivation of complete dynamics remain the principal tasks ahead.



## Citation Policy and the Classification of Self-Hypotheses

This paper verifies correspondence with known mathematics and known physics primarily through external literature. However, among the constructions used here, those already published as the author's own numerical experiments and constructive derivations, and directly required by the derivation chain of this paper, are adopted as self-references — the following three papers only.

1. Noriaki Kihara, **"The Periodic Table of Waves v2 — Particle Classification by Winding-Number Address and Observation Clock, and the Unification of Mass, Lifetime, and Splitting via the Clock Field $\omega(x)$"**, Zenodo, Concept DOI: 10.5281/zenodo.21830706.
2. Noriaki Kihara, **"Two-Layer Separation of Waves and Fields — Unification of Gauge Fields and Gravitational Fields via a Universal Field-Readout Function"**, Zenodo, DOI: 10.5281/zenodo.21832257.
3. Noriaki Kihara, **"Zero Closure Was Four-Dimensional — 'Central Projection' Survives Even in the Complex World"**, Zenodo, Version DOI: 10.5281/zenodo.21902806, Concept DOI: 10.5281/zenodo.21902805.

Results depending on these three papers are never conflated with theorems established by external theory; they are explicitly marked as

$$
\boxed{\text{previously derived results by self-hypothesis and self-construction}}
$$

throughout.

In particular, what this paper imports from "The Periodic Table of Waves v2" is the $Z_2$ decomposition into odd and even harmonics, the sign under half-period translation, the reflectance determined by the odd-harmonic ratio, the clock covering degree, the antipodal two-point structure, and the Bose/Fermi/mixed statistical classification. These are not placed in the "underived" category here; they are classified as **derived, constructed, and numerically confirmed as self-hypotheses**.

Likewise, the separation of the wave-existence layer from the observational readout layer, and the placement of clock, curvature, and gravity on the readout side, imported from "Two-Layer Separation of Waves and Fields", are treated as **previously derived results as self-hypotheses**.

By this citation policy, the paper clearly separates "known results from external literature", "rigorous derivations within this paper", "previously derived self-hypothesis results originating in the three papers above", and "unresolved generalization tasks".


## Reconstructing Curvature, Spacetime Signature, Quantization, and Internal Symmetry from a Zero Closure Containing Unobservable Complex Axes

**Version:** public v1.0  
**Date:** 2026-08-20  
**Citation policy:** Self-citation is kept to the minimum. Among the author's published papers, only the three primary sources whose results are treated here as "already derived, constructed, and numerically confirmed" are self-cited; external literature is used for all other mathematical and mathematical-physics background, known theorems, and comparison with prior work.

---

## Abstract

The central feature of this paper is that, within the range of what has been derived, no continuous free parameters tuned from outside are introduced, and the only externally specified quantity is a discrete integer $N$. Extending some of the already constructed and numerically confirmed results into general theorems, and deriving complete dynamics including local interactions in closed form from the axiom system, are explicitly left as the principal tasks beyond this paper.

This paper examines whether the many symmetries and geometric structures introduced separately in modern theoretical physics can be reconstructed in a unified way from a small number of closure principles.

The central axiom is, throughout,

$$
\boxed{
\mathrm{A1}:\quad
\sum_{n=1}^{M} x_n^2=0,
\qquad x_n\in\mathbb C
}
$$

In this paper, the curvature radius $R$, the time axis $t$, and the internal axis $Q$ are not added outside this equation. They are all read as components of the complex $x_n$ that carry a sign different from the observable real directions.

For example, for the four components

$$
(x,y,z,iR)
$$

A1 gives

$$
x^2+y^2+z^2+(iR)^2=0
$$

that is,

$$
\boxed{
x^2+y^2+z^2=R^2
}
$$

Therefore $\sum x_n^2=R^2$ is neither an axiom separate from A1 nor a different level set. It is **the real-form display of A1 with one component read as the unobservable imaginary direction $iR$**.

Similarly, for

$$
(x,y,z,it,iR,iQ)
$$

we get

$$
x^2+y^2+z^2-t^2-R^2-Q^2=0.
$$

Thus the indefinite signature of spacetime, the curvature radius, and the internal axis all emerge from the same quadratic form as real-form displays of complex axes, without breaking the zero closure.

This is the central point of the paper.  
Curvature need not be added as a separate field. The real-form display of an unobservable complex axis appears on the observable side as a curvature radius. By the same principle, the negative sign of the time axis is not introduced from an external Lorentz metric: reading the unobservable axis $t$ as $it$ produces $(it)^2=-t^2$. Hence $R$, $t$, and $Q$ are unobservable complex axes of the same kind, and the origin of the negative signs is one and the same.

Furthermore,

$$
\boxed{
\mathrm{A2}:\quad U^N=I
}
$$

provides finite cyclicity, discrete phases, cyclotomic eigenvalues, and Born-type squared weights;

$$
\boxed{
\mathrm{A3}:\quad
\text{all pairwise distances of }N\text{ vertices are consistent as a simplex}
}
$$

converts the relational set into geometry; and

$$
\boxed{
\mathrm{A4}:\quad X=\mathcal F(X)
}
$$

provides self-consistent fixed points, stabilizers, and symmetry selection.

From these few principles, this paper classifies in the forward direction the orthogonal groups, cyclic groups, permutation groups, symplectic structure, Lorentz-type signature, constant-curvature geometry, conformal structure, simplicial chain complexes, fixed-point symmetries, and further the automorphism group of the five-degree-of-freedom structure obtained from this axiom system, and examines whether the internal symmetry groups used in the Standard Model appear as a result.

---

# 1. The Starting Point Is a Single Zero Closure

The first axiom of this paper is

$$
\boxed{
\sum_n x_n^2=0
}
$$

and nothing else.

What matters is that

$$
x_n\in\mathbb C
$$

When this equation is written as

$$
\sum_{\rm visible} x_a^2
=
R^2
$$

the right-hand side $R^2$ has not been introduced as a new conserved quantity.

The original equation contains a component

$$
x_R=iR
$$

and since

$$
x_R^2=(iR)^2=-R^2
$$

we have merely transposed

$$
\sum_{\rm visible}x_a^2-R^2=0
$$

Therefore

$$
\boxed{
\sum_{\rm visible}x_a^2=R^2
}
$$

is the decomposition display of

$$
\boxed{
\sum_nx_n^2=0
}
$$

into observable and unobservable components.

---


## Anonymity Does Not Mean Arbitrariness — Selection of Symmetry Sectors by Strong Constraints

In this construction, no physical labels such as "space", "time", "curvature", "charge", or "color" are assigned in advance to the fundamental relational components. The same zero closure can therefore be read as multiple contractions, refinements, and decompositions depending on the observation map.

For example, collecting

$$
r^2=x^2+y^2+z^2
$$

we may read

$$
r^2=t^2+R^2+Q^2
$$

Contracting further as

$$
R'^2=t^2+R^2+Q^2
$$

gives the three-dimensional ellipsoidal/spherical readout

$$
r^2=R'^2,
\qquad
x^2+y^2+z^2=R'^2
$$

On the other hand, the Lorentz-type decomposition

$$
x^2+y^2+z^2-t^2=R^2+Q^2
$$

is also possible, and refining the internal readout as

$$
Q^2=Q_1^2+Q_2^2+Q_3^2
$$

we may read

$$
x^2+y^2+z^2-t^2
=
R^2+Q_1^2+Q_2^2+Q_3^2
$$

However,

$$
\boxed{\text{anonymity}\neq\text{arbitrariness}}
$$

The existence of multiple readout candidates does not mean that physical interpretations can be freely invented to fit any observation. The states, decompositions, and readouts permitted in this construction must simultaneously satisfy the strong constraints imposed independently in this paper: zero closure, finite-order recurrence, simplex consistency, harmonic structure, and self-consistent fixed points.

Conceptually, the allowed set is expressed as

$$
\boxed{
\mathcal S_{\rm allowed}
=
\mathcal S_{\rm closure}
\cap
\mathcal S_{\rm recurrence}
\cap
\mathcal S_{\rm simplex}
\cap
\mathcal S_{\rm harmonic}
\cap
\mathcal S_{\rm self-consistent}
}
$$

Therefore, even if the space of input candidates is large, the allowed solution set that actually remains is restricted by this intersection.

The characteristic of this model is not the production of many structures through free interpretation. On the contrary, its direction is

$$
\boxed{
\text{extremely few inputs}
+
\text{extremely strong constraints}
\Longrightarrow
\text{dense geometry, algebra, and symmetry}
}
$$

In particular, the essential discrete parameter externally scannable in this paper is $N$, and

$$
M=\frac{N(N-1)}2
$$

is determined by $N$. Names such as space, time, curvature, and internal quantity are not additional parameters; they are observational namings of the symmetry sectors that survive the strong constraints.

Accordingly, "selecting a symmetry" does not mean choosing an arbitrary group from outside. It means that among the candidate decompositions, the sectors that simultaneously satisfy the axiom system and the self-consistency condition and close stably are selected, and the observer reads their stabilizer / automorphism as a symmetry.

$$
\boxed{
\text{anonymous relational system}
\longrightarrow
\text{selection of allowed sectors by strong constraints}
\longrightarrow
\text{stable geometry and symmetry}
\longrightarrow
\text{readout as physical quantities}
}
$$

This distinction bears directly on the predictive power of the paper. Even though multiple readout candidates exist, candidates that fail the constraints are rejected, so this is not a construction that "can explain anything". Rather, the fact that many mutually independent structures must be satisfied simultaneously from the same few axioms strengthens the falsifiability of this construction.


# 2. A1: Complex Zero Closure

Let

$$
x=a+ib,
\qquad a,b\in\mathbb R^M
$$

Then

$$
0=x^Tx
=
a^Ta-b^Tb+2i\,a^Tb.
$$

Therefore

$$
\boxed{
a^Ta=b^Tb
}
$$

and

$$
\boxed{
a^Tb=0
}
$$

hold simultaneously.

## Theorem 1: Isometric Orthogonal Two-Planes

A nonzero solution of A1 requires, within the real relational space,

$$
\boxed{
\|a\|=\|b\|,
\qquad
a\perp b
}
$$

That is, A1 is not a mere "zero sum": it has isometric orthogonal two-planes built in.

The global phase

$$
x\mapsto e^{i\theta}x
$$

is represented as a rotation within this plane.

---

# 3. Observable and Unobservable Axes

Split the components of A1 into components $r_a$ read as real directions on the observation side and components $s_\alpha$ carrying $i$ on the unobservable side.

$$
x=
(r_1,\ldots,r_p,
is_1,\ldots,is_q).
$$

Then A1 becomes

$$
\sum_{a=1}^{p}r_a^2
-
\sum_{\alpha=1}^{q}s_\alpha^2
=0.
$$

Hence, from the complex zero closure, an indefinite quadratic form of signature $(p,q)$ appears automatically in the real display. The signature here is not a metric signature imposed from outside, but the count of **the number of axes read as real, $p$, and the number of axes read as unobservable imaginary, $q$**.

$$
\boxed{
Q_{p,q}
=
r^2-s^2
}
$$

In this sense the negative sign is not the result of imposing a Minkowski metric from outside; it arises from the real-form display of complex axes,

$$
\boxed{
(i s_\alpha)^2=-s_\alpha^2
}
$$

---

# 4. The $R$ Axis Is Not a New Constant but a Complex Component

As the minimal example, consider

$$
x=(x,y,z,iR)
$$

A1 is

$$
x^2+y^2+z^2-R^2=0.
$$

Therefore

$$
\boxed{
x^2+y^2+z^2=R^2
}.
$$

Here $R$ is not a curvature radius added to the right-hand side.

It is a component inside the first axiom,

$$
\boxed{
x_R=iR
}
$$

If the observer does not read the $R$ direction directly and reads only $(x,y,z)$, the visible three directions satisfy the spherical condition of radius $R$.

Therefore **the coordinate value of a hidden complex axis is read in visible space as a curvature radius**.

This is the central geometric interpretation of the paper.

---

# 5. Local Curvature Is a Readout of the $R$ Axis

If a local section of visible space satisfies

$$
x^2+y^2+z^2=R^2
$$

the visible quadric is a sphere and the sectional curvature is

$$
\boxed{
K=\frac{1}{R^2}
}
$$

In general, on a $d$-sphere $S^d(R)$,

$$
R_{abcd}
=
\frac1{R^2}
(g_{ac}g_{bd}-g_{ad}g_{bc}),
$$

$$
R_{ab}
=
\frac{d-1}{R^2}g_{ab},
$$

$$
\boxed{
\mathcal R
=
\frac{d(d-1)}{R^2}
}.
$$

Therefore, if $R$ differs from place to place on the visible side,

$$
R=R(q)
$$

is **a local coordinate of the unobservable axis** and can simultaneously be read on the visible side as a local curvature radius.

What matters is that

$$
\boxed{
R(q)\text{ has not been added as an external gravitational field}
}
$$

$R$ is an internal component of A1.

---

# 6. The Time Axis Works the Same Way

For

$$
x=(x,y,z,it)
$$

A1 is

$$
x^2+y^2+z^2-t^2=0.
$$

This is a Lorentz-type null cone.

Therefore the negative sign of time also appears as the real-form display of an unobservable complex axis,

$$
\boxed{it}
$$

There is no need to assume a separate Minkowski metric here.

More generally, for

$$
x=(x,y,z,it,iR,iQ)
$$

we get

$$
\boxed{
x^2+y^2+z^2-t^2-R^2-Q^2=0
}.
$$

Therefore

- visible spatial axes,
- the time axis,
- the curvature axis $R$,
- the internal axis $Q$,

need not be posited as separate fundamental principles. In particular, the negative signs of $t$, $R$, and $Q$ all come from the same complex squaring rule

$$
\boxed{(iu)^2=-u^2}
$$

All are real-form displays of different complex directions of the same

$$
\boxed{
\sum_nx_n^2=0
}
$$


### Theorem 2: Lorentz Signature Arises from the Square of Imaginary Axes

Let the observable axes be $r_1,\ldots,r_p$ and the unobservable axes $u_1,\ldots,u_q$, and write the complex state as

$$
x=(r_1,\ldots,r_p,iu_1,\ldots,iu_q)
$$

A1 is

$$
\sum_{a=1}^{p}r_a^2+\sum_{\alpha=1}^{q}(iu_\alpha)^2=0
$$

that is,

$$
\boxed{
\sum_{a=1}^{p}r_a^2-\sum_{\alpha=1}^{q}u_\alpha^2=0
}
$$

Therefore the signature $(p,q)$ appearing in the real display is not an independently introduced metric signature but is determined by **the numbers of visible real axes and unobservable imaginary axes**.

In particular, for $p=3,q=1$,

$$
\boxed{
x^2+y^2+z^2-t^2=0
}
$$

and the Lorentz-type signature is a direct consequence of

$$
\boxed{(it)^2=-t^2}
$$

$\square$

In this sense the sign origins of the $R$ axis and the $t$ axis are exactly identical.

---

# 7. Conservation Groups Emerging from A1

The complex linear conservation group of the quadratic form

$$
Q(x)=x^Tx
$$

is

$$
\boxed{
O(M,\mathbb C)
}.
$$

When the signature is chosen as $(p,q)$ in the real display, the real conservation group is

$$
\boxed{
O(p,q)
}.
$$

The decomposition appearing in this six-component readout of the construction is not $3+1$ but

$$
\boxed{
3+3=(x,y,z)+(t,R,Q)
}
$$

That is, three axes on the visible side and three on the unobservable side.

For the $(x,y,z,t)$ partial readout, the automorphism group of its real quadratic form appears:

$$
\boxed{
O(3,1)
}
$$

The connected, orthochronous component is

$$
SO^+(3,1),
$$

and its double cover is

$$
\boxed{
Spin^+(3,1)\cong SL(2,\mathbb C)
}.
$$

---

# 8. Null Cone, Curvature, and Spacetime Are Different Readouts of the Same A1

A1 itself is always

$$
\sum_nx_n^2=0.
$$

But when the observer cannot read some axes directly, it appears as

$$
\sum_{\rm visible}r_a^2
=
\sum_{\rm hidden}s_\alpha^2
$$

If the right-hand side has just one axis, it can be read as a curvature radius $R$.

If the right-hand side includes the time axis, it can be read as a Lorentz-type null structure.

Including several unobservable axes gives

$$
r^2=t^2+R^2+Q^2+\cdots
$$

Therefore

$$
\boxed{
\text{zero closure}
\to
\text{signature}
\to
\text{curvature}
\to
\text{spacetime}
}
$$

is not a sequence of different axioms but differences of decomposition, projection, and readout of the same A1.

---

# 9. Normalization Erodes This Geometry

A1 is a homogeneous condition:

$$
Q(\lambda x)=\lambda^2Q(x).
$$

Therefore A1 does not kill the scale direction.

On the other hand, imposing

$$
x^\dagger x=1
$$

as a fundamental axiom fixes the radius of the Hermitian norm.

If the magnitudes of the unobservable complex axes are read on the visible side as

- curvature radius,
- clock scale,
- internal scale,

then global normalization risks freezing those degrees of freedom.

Here one must not confuse

$$
x^Tx
$$

with

$$
x^\dagger x
$$

$$
\boxed{
x^Tx=0
}
$$

is the geometric zero closure.

$$
\boxed{
x^\dagger x=1
}
$$

is a section choice of the Hermitian norm.

---

# 10. Complex Structure and Symplectic Structure

In the real display

$$
x=(a,b)
$$

the complex structure is

$$
J=
\begin{pmatrix}
0&-I\\
I&0
\end{pmatrix},
\qquad
J^2=-I.
$$

For a compatible metric $g$,

$$
\omega(u,v)=g(Ju,v)
$$

is an antisymmetric two-form.

Therefore

$$
\boxed{
\text{complex}
\to
\text{orthogonal}
\to
\text{symplectic}
}
$$

are not mutually independent extra assumptions; under the compatibility condition they are different aspects of the same complex geometry.

---

# 11. A2: $U^N=I$ Is Self-Closure of the Action

If

$$
\boxed{
U^N=I
}
$$

then every eigenvalue $\lambda$ satisfies

$$
\lambda^N=1
$$

Therefore

$$
\boxed{
\lambda_m=e^{2\pi i m/N}
}.
$$

If the operator has exact order $N$, then

$$
\boxed{
\langle U\rangle\cong\mathbb Z_N
}.
$$

Furthermore, the cyclotomic field

$$
\mathbb Q(\zeta_N)
$$

and

$$
\operatorname{Gal}
(\mathbb Q(\zeta_N)/\mathbb Q)
\cong
(\mathbb Z/N\mathbb Z)^\times
$$

are naturally attached.

A2 therefore provides simultaneously

- finite recurrence,
- phase discretization,
- cyclic groups,
- cyclotomic structure.

---

# 12. Born-Type $\cos^2/\sin^2$ Is a Readout of Finite Recurrence

Imposing

$$
U^N=I
$$

on the two-channel rotation

$$
U(\theta)
=
\begin{pmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{pmatrix}
$$

gives

$$
\theta=\frac{2\pi m}{N}.
$$

The squared projections onto the axes are

$$
\boxed{
\cos^2\theta,\qquad\sin^2\theta
}.
$$

Therefore Born-type squared weights can be generated, at least in two-channel systems, as

$$
\boxed{
\text{finite-order phase recurrence}
\to
\text{squared readout}
}
$$

Born-type squared weights are not a starting axiom.

---

# 13. A3: The Same $N$ Becomes the Vertex Count

Let the $N$ of A2 simultaneously be the number of geometric vertices.

Then the total number of pairwise relations is

$$
\boxed{
M=\frac{N(N-1)}2
}.
$$

Take each relation to be a scalar distance $d_{ij}$ and assume all distances are consistent as a simplex.

This establishes

$$
\boxed{
\text{number of phase divisions}
=
\text{number of vertices}
}
$$

A2 and A3 carry no independent integers; they are constrained to the same $N$.

This is a strong condition directly linking quantization and geometry.

---

# 14. Simplex Closure Creates Space

The pairwise distances of $N$ vertices are the edge data of

$$
K_N
$$

If the distances are consistent, the embedding dimension and shape are determined by the Gram matrix or the Cayley--Menger determinant.

Therefore

$$
\boxed{
\text{relational set}
\to
\text{distance geometry}
}
$$

holds.

Absolute positions are unnecessary.

Global

- translations,
- rotations,
- reflections,

do not change the distance data.

Hence they become redundancies in the relational description.

---

# 15. The Simplex Also Generates a Chain Complex

The boundary operator of oriented simplices satisfies

$$
\boxed{
\partial^2=0
}
$$

Dualizing,

$$
\boxed{
d^2=0
}.
$$

Therefore the simplex is not merely a "solid": it generates

- homology,
- cohomology,
- discrete differential forms,
- closed / exact structure.

This is the natural entry to the formal structure required by gauge geometry.

---

# 16. A4: Self-Consistency Selects Which Geometry Is Realized

Require

$$
\boxed{
X=\mathcal F(X)
}
$$

If $\mathcal F$ is equivariant under a group $G$, the conservation group of a fixed point $X_*$ is

$$
\boxed{
G_{X_*}
=
\{g\in G\mid gX_*=X_*\}
}.
$$

Therefore a symmetry selection

$$
\boxed{
G\to H
}
$$

occurs naturally.

The eigenvalues of the linearization

$$
\delta X'
=
D\mathcal F(X_*)\delta X
$$

select stable, neutral, and amplified directions.

Hence the phenomenon in which only a few directions become macroscopic can be described as the spectral selection of self-consistent fixed points.

---

# 17. A2 and A4 Are Different Levels of Self-Closure

$$
U^N=I
$$

is the finite recurrence of the whole operator.

$$
X=\mathcal F(X)
$$

is the fixed-point closure of states and geometry.

Extending to periodic fixed points

$$
\mathcal F^N(X)=X
$$

the two become different realizations of the same "finite self-recurrence".

Therefore, taking

$$
\boxed{
\text{self-closure}
}
$$

as the superordinate principle, we can classify as one family:

- zero closure of quadratic forms,
- finite closure of operators,
- simplex closure of distances,
- fixed-point closure of generating maps.

---

# 18. If $3+1$ Emerges, Lorentz Is Not an Additional Axiom

Suppose the complex-axis decomposition of A1 and the geometric selection of A3/A4 stably select

$$
\boxed{
3\text{ visible real axes}
+
1\text{ unobservable imaginary axis}
}
$$

Writing the unobservable axis as $it$, its square is automatically $-t^2$.

Then A1, in its real display, gives

$$
x^2+y^2+z^2+(it)^2=0
$$

that is,

$$
\boxed{x^2+y^2+z^2-t^2=0}
$$

The conservation group of this indefinite quadratic form is

$$
O(3,1)
$$

Therefore the unresolved point is neither the negative sign of Lorentz nor whether three visible directions arise. The negative sign is automatic from the square of $it$.

Moreover, three-dimensionality itself is not treated as a mere accident of numerical experiments. From the complex zero closure of the first axiom arises the equal-norm orthogonal structure of real and imaginary parts, and three directions are established by the orthogonal two-planes and their independent normal direction. In addition, from the consistency of phase closure and simplex geometry, the visible section is represented as a three-dimensional ellipsoidal structure with three mutually distinguishable principal axes.

Accordingly, the geometric path

$$
\boxed{
\text{zero closure}
\Longrightarrow
\text{equal-norm orthogonal structure}
\Longrightarrow
\text{orthogonal two-planes + normal}
\Longrightarrow
\text{3-dimensional ellipsoidal structure}
}
$$

is classified as the already-derived part of three-dimensionality.

In parallel, the author's published paper "Zero Closure Was Four-Dimensional" numerically confirms that, of at most $N-1$ nontrivial principal axes, the top three directions $A,B,C$ amplify selectively, the middle directions remain almost stationary, and some lower directions migrate to imaginary directions. Therefore

$$
\boxed{
\text{phase and closure geometry}
\Longrightarrow
\text{3-dimensional ellipsoidal structure}
}
$$

and

$$
\boxed{
\text{evolution of a self-consistent coherent parent}
\Longrightarrow
\text{spectral concentration and selective amplification onto 3 principal axes}
}
$$

are distinct grounds, yet point to the same three-dimensional visible structure.

Hence we do not classify "the principle generating three directions" as unresolved. The geometric establishment of three directions and their dominance and saturation are placed in the domain already derived and confirmed in the self-papers.

What remains is the **elevation of this coincidence to a general theorem of universality**. That is, for any sufficiently large $N$, a wide range of initial conditions, seeds, and admissible self-consistent maps, to determine analytically the necessary and sufficient conditions under which

$$
\boxed{
\text{3-dimensional ellipsoidal structure}
\Longleftrightarrow
\text{selective dominance of 3 principal axes}
}
$$

holds. In other words, the general analytic mechanism by which a self-consistent coherent parent, a fermionic seed with odd harmonics, and selective, inflation-like amplification preserving the zero closure choose three anisotropic principal axes out of many relational directions has not yet been determined.

Therefore the future task is not "to create three directions" but

$$
\boxed{
\text{to generalize the already reproduced spontaneous generation of three directions and to identify its selection principle analytically}
}
$$

---

# 19. Because of the $R$ Axis, Local Curvature Is Already Contained in A1

Conventional organizations tend to treat "the generation of local curvature" as a separate problem.

But if A1 contains the unobservable component

$$
iR
$$

then in the visible three directions

$$
x^2+y^2+z^2=R^2
$$

Therefore

$$
\boxed{
R^{-2}
}
$$

is precisely the curvature of the visible sphere.

If locally

$$
R=R(q)
$$

the curvature radius of the visible section observed at each point varies.

Therefore "the place where gravity lives" is not the addition of an external force to the state; it already exists inside A1 as

$$
\boxed{
\text{the local readout of an unobservable complex axis}
}
$$

General Riemannian dynamics remains a derivation task, but **local curvature as a geometric quantity need not be added as a new degree of freedom**.

---

# 20. Degrees of Freedom of the Complex Zero Closure, Readout Resolution, and Two Routes to the Standard Model Gauge Group


> **Note on layers and selection.** The $\mathbb C^6$, $3+3$, and $3+2$ below do not denote a unique ontological number of axes fixed in the fundamental relational system. The essential discrete parameter given externally on the fundamental side is $N$, and the number of relational components is fixed as $M=N(N-1)/2$. $\mathbb C^6$ is one effective representation reading the zero closure contracted into six registers, and the complex-dimension count $6-1=5$ after adopting that representation is rigorous. Meanwhile, refinements like $Q^2=Q_1^2+Q_2^2+Q_3^2$ and re-contractions like $r^2=t^2+R^2+Q^2$ are other symmetric readouts of the same closure. The problem is therefore not "to derive the true number of axes uniquely from the axioms" but to classify the readout sectors, and their symmetries, that close stably while satisfying the strong constraints simultaneously.

The first axiom is read from the start as a single complex zero closure on complex space

$$
\boxed{\sum_{n=1}^{6}X_n^2=0,\qquad X_n\in\mathbb C}
$$

Accordingly, the degree-of-freedom count here is not an argument subtracting one real constraint from six real axes. Correctly,

$$
\boxed{\dim_{\mathbb C}\mathbb C^6-1=5}
$$

In the real display, the single complex constraint corresponds to two real conditions, leaving 10 real dimensions out of 12 — again five complex dimensions.

Reading, as one observation map,

$$
\boxed{x^2+y^2+z^2=t^2+R^2+Q^2}
$$

the components $(x,y,z,it,iR,iQ)$ are ontologically all the same kind of complex axis, and the distinctions visible/unobservable and time/curvature/internal are attached on the readout side. In this display, $t$ can be placed as a dependent readout that is not read directly,

$$
\boxed{t^2=x^2+y^2+z^2-R^2-Q^2}
$$

In the sector where this readout adopts $t$ as a real clock quantity,

$$
\boxed{
x^2+y^2+z^2\ge R^2+Q^2
}
$$

is required, and on the equality surface

$$
\boxed{
x^2+y^2+z^2=R^2+Q^2
}
$$

we have $t=0$. This is a boundary arising automatically from the readout in question. At this stage the paper does not identify it with the event horizon of general relativity; it is classified as a **horizon-type boundary candidate**.


Therefore $3+2$ is not a number introduced from outside to match the Standard Model; it is one decomposition that reads, via an observation map, the five complex degrees of freedom obtained by imposing one complex zero closure on six complex axes.

## 20.1 Coarse-Grained and Refined Readouts

What matters is not to identify $Q$ with the number of fundamental components itself. $Q$ is a register that reads the unobservable-side closure quantity in contracted form; at a finer resolution it can be expanded as

$$
\boxed{Q^2=Q_1^2+Q_2^2+Q_3^2}
$$

Therefore

$$
x^2+y^2+z^2-t^2=R^2+Q^2
$$

and

$$
\boxed{x^2+y^2+z^2-t^2=R^2+Q_1^2+Q_2^2+Q_3^2}
$$

are not two competing ontologies; they can be treated as a coarse-graining/refinement relation reading the same zero closure at different internal resolutions.

With this distinction, one must not confuse the number of relational components counted by A3's

$$
M=\frac{N(N-1)}2
$$

with the number of readout axes $(r,t,R,Q_i,\ldots)$. Expanding readout registers is not, by itself, an operation that changes $M$.

## 20.2 The Coarse-Grained Route to the Standard Model Gauge Group

If the five complex degrees of freedom preserve a self-consistent $3\oplus2$ decomposition

$$
V=V_3\oplus V_2,
\qquad
\dim_{\mathbb C}V_3=3,
\qquad
\dim_{\mathbb C}V_2=2
$$

and each subspace is granted a Hermitian structure, the conservation group is

$$
U(3)\times U(2)
$$

If, further, the conservation condition removing the physically redundant global phase is realized as removal of the total determinant phase, we obtain

$$
\boxed{
S(U(3)\times U(2))
\cong
\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6}
}
$$

Note that it is not necessary to identify $V_3$ directly with visible space $(x,y,z)$ and let color $SU(3)$ act there. The central principle of this paper is the separation of the "complex zero closure of the existence layer" from the "observational readout", and the target on which internal symmetries act must be identified separately as a readout register.

Also, since in A1 the global phase

$$
X\mapsto e^{i\theta}X
$$

is treated as an unobservable common phase, the $\det=1$ condition is not a new independent axiom: it can very likely be derived as the operation removing this global-phase redundancy from the physical automorphism group. Rigorous treatment including the identification of the discrete center remains, but the residual task is not "assuming det=1 from outside"; it is "determining how far the global-phase quotient of A1 uniquely fixes the global form of $S(U(3)\times U(2))$".

## 20.3 The Refined Route Expanding the Internal Triplet

In the refined readout

$$
R^2+Q_1^2+Q_2^2+Q_3^2
$$

the internal triplet $(Q_1,Q_2,Q_3)$ appears independently of the three visible spatial axes. Preserving the complex Hermitian structure of this triplet gives an internal $U(3)$, and separating the global phase gives an entry to $SU(3)$.

In the readout collecting $(R,Q_1,Q_2,Q_3)$ as four complex components, known group theory provides the alternative route

$$
U(4)\supset SU(4),
\qquad
SU(4)\supset SU(3)\times U(1)
$$

Meanwhile, the four components on the $(x,y,z,t)$ side admit a Lorentz readout and a Euclidean-type readout; the double cover of the latter exhibits

$$
Spin(4)\cong SU(2)\times SU(2)
$$

Hence the refined readout also has a known mathematical connection to the Pati--Salam type

$$
SU(4)\times SU(2)\times SU(2)
$$

However, at the present stage this paper does not classify this as the physical gauge group uniquely selected by A1--A4. What matters is that the coarse-grained $3+2$ route and the internal-triplet expansion route arise from different readout resolutions of the same complex zero closure.

## 20.4 The $M=6$ Complex Zero Quadratic Form and Known Geometry

The projectivized zero set of the $M=6$ A1 can be treated as the nondegenerate complex quadric hypersurface $Q^4$ in the complex projective space $\mathbb{CP}^5$. Quadrics of this kind connect with Grassmann/Plücker geometry and with the accidental isomorphisms of six-dimensional orthogonal groups. Therefore, by choosing different real forms and observational readouts of the same complex zero closure, there exist external mathematical routes connecting to $SO(3,3)$, to the $SO(4,2)$ appearing in conformal geometry, and to $Spin(6)\cong SU(4)$.

This is not a claim to have proved the Standard Model directly from A1. But it shows that the central proposition of this paper,

$$
\boxed{\text{one complex zero geometry}\longrightarrow\text{different readouts}\longrightarrow\text{different physical symmetries}}
$$

connects naturally to the known mathematics of complex quadric hypersurfaces and their real forms. It is important as an external mathematical cross-check route independent of the coarse-grained $3+2$ route.

## 20.5 Derivation Stages at Present

**Directly obtained from A1**

$$
\boxed{\mathbb C^6\cap\{\sum X_n^2=0\}\Longrightarrow\dim_{\mathbb C}=5}
$$

and, as its observational readout,

$$
\boxed{x^2+y^2+z^2=t^2+R^2+Q^2}.
$$

**Already constructed and numerically confirmed, including the self-papers**

The dominance of the three visible principal axes, the readout decomposition including $R/Q$, the refinement of internal structure, and the particle-type classification by harmonic structure.

**Group-theoretically rigorous conditional part**

$$
3\oplus2+\text{preservation of the Hermitian decomposition}+\text{removal of the global phase}
\Longrightarrow
S(U(3)\times U(2))
$$

and

$$
S(U(3)\times U(2))
\cong
\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6}.
$$

**The remaining generalization and dynamics tasks** are to generalize which $3+2$ readout or internal triplet the rigorously obtained five complex degrees of freedom select self-consistently — not to simply rebuild readouts already constructed. They are to determine which readout resolution and which stabilizer are universally selected by self-consistent dynamics, how local gauge connections arise, and whether chirality, hypercharge, and anomaly cancellation close within the same derivation chain.

# 21. Revised Derivation Map


The top-level principle for reading the derivation map is

$$
\boxed{
N
\Longrightarrow
M=\frac{N(N-1)}2
\Longrightarrow
\text{anonymous relational system}
\Longrightarrow
\text{sector selection by strong constraints}
\Longrightarrow
\text{readout of symmetry and geometry}
}
$$

Accordingly, $3+3$, $3+2$, the internal $Q$ triplet, etc. are not mutually competing "fundamental axis counts"; they are classified as readouts of different resolutions and symmetry sectors of the same strongly constrained relational system. The admissibility of multiple readouts is not, in itself, the addition of free parameters.

| Structure | Origin | Status |
|---|---|---|
| Anonymity / strong-constraint sector selection | $N\to M=N(N-1)/2$; intersection of zero closure, finite recurrence, simplex, harmonics, self-consistency | Multiplicity of readouts treated as sector selection, not arbitrariness |
| Complex zero quadratic form | A1 | Axiom |
| Isometric orthogonal two-planes | A1 | Rigorous |
| Complex phase rotation | A1 | Rigorous |
| $O(M,\mathbb C)$ | Conservation group of A1 | Rigorous |
| Indefinite signature $O(p,q)$ | Real display of unobservable imaginary axes $iu$ | Rigorous |
| Null cone | Real display of A1 | Rigorous |
| $x^2+y^2+z^2=R^2$ | Including $iR$ in A1 | Rigorous |
| Curvature radius $R$ | Visible-sphere readout | Rigorous |
| $K=1/R^2$ | Spherical geometry | Rigorous |
| Negative sign of $t$ | Reading the unobservable axis as $it$, so $(it)^2=-t^2$ | Rigorous |
| Signature $(3,3)$ of the six-component readout | Real display of $(x,y,z,it,iR,iQ)$ | Rigorous |
| $O(3,3)$ | Conservation group of the $(3,3)$ quadratic form of the six-component readout | Rigorous |
| Lorentz $O(3,1)$ | Spacetime partial readout of $(x,y,z,t)$ | Conditionally rigorous |
| Distinguishability of 3 visible principal axes | Anisotropic ellipsoid $a\neq b\neq c$ | Rigorous within the range where self-consistent principal-axis solutions hold |
| Spin$(3,1)$ | Lorentz lift | Conditionally rigorous |
| Conformal structure | Null cone / scale | Strong known connection |
| Symplectic $J$ | Real display of the complex structure | Rigorous entry |
| $\mathbb Z_N$ | A2 | Rigorous |
| Cyclotomic eigenvalues | A2 | Rigorous |
| Galois symmetry | Cyclotomic field | Rigorous |
| Born-type squared weights | A2 + two-channel projection | Rigorous |
| $K_N$ | A3 | Rigorous |
| $S_N$ | Anonymous vertices | Rigorous |
| Simplex distance geometry | A3 | Rigorous |
| $\partial^2=0$ | Simplex chain complex | Rigorous |
| Cohomology | A3 | Rigorous |
| Stabilizer | A4 | Rigorous |
| $G\to H$ | A4 | Conditionally rigorous |
| 3-dimensional ellipsoidal structure | A1: equal-norm orthogonal structure → orthogonal two-planes + normal | Already derived |
| Selective dominance of 3 principal axes | A4 / `make_parent` self-consistent spectrum | Numerically confirmed in the self-papers; universality and necessary-and-sufficient conditions remain to be generalized |
| Complex six-axis zero closure | A1 on $\mathbb C^6$: one complex constraint | $\dim_{\mathbb C}=5$ derived rigorously |
| Coarse-grained readout | $x^2+y^2+z^2=t^2+R^2+Q^2$ | Derived as an observational display of the same complex zero closure |
| Refined internal readout | $Q^2=Q_1^2+Q_2^2+Q_3^2$ | Consistent with the readout conventions of the self-papers; nested with the coarse-grained readout |
| Independent degrees of freedom | Six complex dof $-$ one complex zero closure | Rigorous as $6-1=5$ complex dof; $3+2$ is one observational decomposition of it |
| $U(3)\times U(2)$ | Preservation of the Hermitian $3\oplus2$ decomposition of the 5 complex dof | Conditionally rigorous; the selection rule as a self-consistent stabilizer is the generalization task |
| $S(U(3)\times U(2))$ | The above + removal of the global-phase redundancy of A1 | Conditionally rigorous; rigorous treatment of the global form including the discrete center remains |
| Standard Model gauge Lie algebra | Lie algebra of $S(U(3)\times U(2))$ | $\mathfrak{su}(3)\oplus\mathfrak{su}(2)\oplus\mathfrak{u}(1)$: conditionally rigorous |
| SM global gauge group | $S(U(3)\times U(2))$ | $[SU(3)\times SU(2)\times U(1)]/\mathbb Z_6$: conditionally rigorous |
| Internal $Q$ triplet | Refined readout of $(Q_1,Q_2,Q_3)$ | Entry to an internal $SU(3)$ candidate independent of visible $(x,y,z)$ |
| $SU(4)$ route | Four-complex-component readout of $(R,Q_1,Q_2,Q_3)$ | Known connections to $Spin(6)\cong SU(4)$ and $SU(4)\supset SU(3)\times U(1)$ |
| Pati--Salam-type route | Internal 4 components + $Spin(4)\cong SU(2)\times SU(2)$ | Second, externally mathematical route; physical selection not yet generalized |
| Complex quadric $Q^4\subset\mathbb{CP}^5$ | Projective zero set of the $M=6$ A1 | Strong connection to known mathematics; unified cross-check route across real forms and readouts |
| Bose/Fermi/mixed sectors | Odd/even harmonic ratio, $Z_2$ parity, single/double cover | Derived, constructed, and numerically confirmed in the self-papers; mass identification and general spin-statistics correspondence remain |
| General Riemannian dynamics | $R(q)$ + connection | Not yet derived |
| Chirality/hypercharge/anomaly | SM representations | Not yet derived |

---

### 21.1 The Double Derivation Branch Added in This Revision

This revision retains the coarse-grained route

$$
\mathbb C^6\xrightarrow{\sum X_n^2=0}\dim_{\mathbb C}=5\xrightarrow{3\oplus2}S(U(3)\times U(2))
$$

while adding a second route that refines the internal readout to

$$
Q^2=Q_1^2+Q_2^2+Q_3^2
$$

This makes it unnecessary to let $SU(3)$ act directly on the three visible spatial axes; it can be read as a symmetry acting on the internal triplet.

Furthermore, the $M=6$ A1 has a known mathematical connection to the complex quadric hypersurface $Q^4\subset\mathbb{CP}^5$, providing an independent route through $Spin(6)\cong SU(4)$. The object of future verification is therefore not "constructing $SU(3)\times SU(2)\times U(1)$ from one accidental $3+2$ decomposition", but **to which stabilizer the multiple readout routes over the same complex zero geometry converge**.

---

# 22. The Problems That Truly Remain

By reading the first axiom correctly, the following cease to be "problems requiring new structure":

- curvature radius,
- the geometric location of local curvature,
- the negative sign of the time axis,
- the null cone,
- indefinite quadratic forms,
- the entry to the Lorentz group,
- the scale degree of freedom.

In particular, $t$, $R$, and $Q$ are unobservable complex axes of the same kind, and their sign differences arise solely from

$$
\boxed{(iu)^2=-u^2}
$$

The remaining central problems narrow down as follows.


## Statistical Structure Taken as Derived Self-Hypotheses

The construction obtained in "The Periodic Table of Waves v2" is adopted in this paper as the derivation chain

$$
\text{odd/even harmonics}
\longrightarrow
T_{\lambda_0/2}\psi=\pm\psi
\longrightarrow
Z_2\text{ parity}
\longrightarrow
\text{single/double cover}
$$

Furthermore, with the odd- and even-harmonic powers $P_{\rm odd},P_{\rm even}$,

$$
\boxed{
r=\frac{P_{\rm odd}}{P_{\rm odd}+P_{\rm even}}
}
$$

gives the exchange/reflection rate. This continuous quantity distinguishes the pure-even endpoint, the pure-odd endpoint, and the mixed sector in between.

Accordingly, in this paper,

$$
\boxed{
\begin{array}{ll}
P_{\rm odd}=0 &\rightarrow \text{Bose-type endpoint},\\[2mm]
P_{\rm even}=0 &\rightarrow \text{Fermi-type endpoint},\\[2mm]
P_{\rm odd}P_{\rm even}>0 &\rightarrow \text{mixed (Ermion-type) sector}
\end{array}
}
$$

is not a classification newly assumed from the external spin-statistics theorem, but **a derived result, as a self-hypothesis, obtained from the author's published model**.

Therefore the future task is not to "generate" the Bose/Fermi structure. What remains is to compare and formalize, under more general conditions, the correspondence between this harmonic-parity-derived statistical structure and the spin-statistics structure of standard quantum field theory.


## Problem 1

For the spontaneous generation of the three anisotropic visible principal axes already reproduced with `make_parent`, identify the generating principle analytically and generalize it beyond the specific numerical construction. In particular, determine which of the following are necessary and sufficient for three-direction selection: a self-consistent coherent parent, a fermionic seed with odd harmonics, zero closure, and selective amplification.

## Problem 2

Among the multiple readout resolutions admitted by the complex zero closure, generalize which Hermitian decomposition and stabilizer are universally selected by self-consistent dynamics. This is not the problem of deriving the $3+2$ degrees of freedom anew; it is the problem of determining the selection rule among the existing routes: the coarse-grained $3\oplus2$, the internal $Q$ triplet, and the $SU(4)$-type refinement.

## Problem 3

Can general Riemann curvature and geodesic dynamics be derived from the self-consistent variation of the unobservable axis

$$
R(q)
$$

## Problem 4

Can local gauge connections / curvature be derived from simplex cochains and phase data?

## Problem 5

The Bose/Fermi/mixed classification itself is already derived, constructed, and numerically confirmed in the self-papers. What remains is to generalize the classification by odd/even harmonic ratio and $Z_2$ parity to arbitrary sufficiently large $N$ and wide interaction conditions, and to close the correspondence with the mass identification and mass hierarchy of real particles and with the spin-statistics structure of standard quantum field theory.

## Problem 6

Can the chirality, hypercharge, and anomaly cancellation of the Standard Model be reproduced from stabilizer representations?

---

# 23. Falsifiability

The strong form of this construction is falsified by any of the following.

1. The self-consistent fixed points of A1--A4 fail to reproduce the observed three-principal-axis structure or the admissible readout sectors.
2. Identifying the $N$ of A2 with the vertex count of A3 eliminates all solutions.
3. Local variation of the $iR$ axis is inconsistent with the geodesic / curvature readout.
4. At large $N$ the three-principal-axis selection disappears, reducing to a mere finite-size effect.
5. The automorphism group of the five-degree-of-freedom structure fails to reproduce the internal symmetries and representation content of the Standard Model.

---

# 24. Conclusion

The central principle of this paper is, throughout,

$$
\boxed{
\sum_nx_n^2=0
}
$$

The curvature radius $R$ is not a quantity added to the right-hand side.

When the complex component

$$
\boxed{
x_R=iR
}
$$

is displayed in real form, it merely looks like

$$
\sum_{\rm visible}x_a^2=R^2
$$

Similarly, the time axis, as the complex direction

$$
it
$$

produces

$$
r^2-t^2=0
$$

Here $t$, $R$, and $Q$ are essentially unobservable complex axes of the same kind, and which axis is called time, curvature, or internal quantity is determined by which observation map reads it and how. The difference of names itself therefore need not be derived uniquely from the axioms.

Hence

$$
\boxed{
\text{curvature}
,\quad
\text{time signature}
,\quad
\text{null structure}
}
$$

are not additions external to the first axiom.

They all emerge from

$$
\boxed{
\text{the observable/unobservable axis decomposition of the complex zero closure}
}
$$

Furthermore, when the six directions are read as

$$
(x,y,z,it,iR,iQ)
$$

the first axiom gives

$$
\boxed{
x^2+y^2+z^2=t^2+R^2+Q^2
}
$$

In this six-component readout the geometry is thus expressed not as $3+1$ but as $3+3$, and the $3+1$ of observed spacetime is a partial readout of it. When the visible side becomes an anisotropic ellipsoid, its three principal axes appear autonomously as mutually distinguishable eigendirections.

The second axiom

$$
U^N=I
$$

provides the finite self-closure of the action;

the third condition, simplex closure, converts phase vertices into distance geometry;

and the fourth condition, self-consistency, selects the fixed points and symmetries realized among them.

The few principles therefore compress to

$$
\boxed{
\text{zero closure}
+
\text{finite recurrence}
+
\text{distance closure}
+
\text{self-consistency}
}
$$

Under this compression, many symmetries of theoretical physics can be reclassified not as independent axioms but as

$$
\boxed{
\text{automorphisms of a self-closing complex relational geometry}
}
$$

The most important task is not to add further symmetry groups. It is

$$
\boxed{
\text{to classify completely the fixed points and stabilizers of A1--A4}
}
$$

---

# References

External literature is primary; only the three primary sources of results derived and confirmed as self-hypotheses are self-cited, minimally.

1. E. W. Weisstein, "Sphere," *MathWorld*. Sphere geometry and curvature references.  
   https://mathworld.wolfram.com/Sphere.html

2. E. W. Weisstein, "Cayley-Menger Determinant," *MathWorld*.  
   https://mathworld.wolfram.com/Cayley-MengerDeterminant.html

3. M. Cvetič and L. Lin, "The Global Gauge Group Structure of F-theory Compactification with U(1)s," arXiv:1706.08521 (2017).  
   https://arxiv.org/abs/1706.08521

4. N. Raghuram, W. Taylor and A. P. Turner, "General F-theory models with tuned $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ symmetry," arXiv:1912.10991 (2019).  
   https://arxiv.org/abs/1912.10991

5. F. Klinker, "An explicit description of $SL(2,\mathbb C)$ in terms of $SO^+(3,1)$ and vice versa," arXiv:1712.02168 (2017).  
   https://arxiv.org/abs/1712.02168

6. J. T. Wheeler, "Weyl geometry," arXiv:1801.03178 (2018).  
   https://arxiv.org/abs/1801.03178

7. O. Macia and Y. Nagatomo, "Holomorphic isometric embeddings of complex Grassmannians into quadrics: The general case," *Kyoto Journal of Mathematics* 66(1), 67–85 (2026).  
   https://doi.org/10.1215/21562261-2024-0033

8. Standard differential geometry result: a $d$-sphere of radius $R$ has constant sectional curvature $1/R^2$, Ricci curvature $(d-1)g/R^2$, and scalar curvature $d(d-1)/R^2$.

9. Standard Lorentz geometry result: the real quadratic form of signature $(3,1)$ is preserved by $O(3,1)$; its proper orthochronous spin double cover is $SL(2,\mathbb C)$.

10. Standard algebraic topology result: the simplicial boundary operator satisfies $\partial^2=0$.

---

## Claim Strength

### Rigorous

- $x^Tx=0\Rightarrow \|a\|=\|b\|,\ a\perp b$
- That A1 containing the complex axis $iR$ gives $r^2=R^2$ on the visible side
- That A1 containing the complex axis $it$ gives the Lorentz-type negative sign on the visible side
- $K=1/R^2$ for visible spherical sections
- $\mathbb Z_N$ and cyclotomic eigenvalues from A2
- $K_N,S_N,\partial^2=0$ from A3
- Stabilizers from A4 + equivariance

### Conditionally Rigorous

- Lorentz / Spin / conformal in the $(x,y,z,t)$ partial readout of $3+3$
- $S(U(3)\times U(2))$ and the Standard Model global gauge group after the $3+2$ selection
- The construction reading $R(q)$ as a local visible curvature radius

## The Only Externally Specified Parameter Is the Discrete Integer $N$

An important feature of this construction is that, apart from the parts not yet derived or generalized, no explicit continuous tuning parameters for fitting phenomena are introduced.

The only structural parameter specified from outside is

$$
\boxed{N\in\mathbb N}
$$

This $N$ is not an ordinary coupling constant or fitting parameter. In this construction the same integer $N$ simultaneously binds

$$
\boxed{
U^N=I,\qquad
|V|=N,\qquad
M=\frac{N(N-1)}2
}
$$

specifying the degree of finite-order recurrence, the number of discrete divisions of the phase circle, the number of vertices of the closure geometry, and the total number of pairwise relations.

Therefore $N$ is not a value tuned continuously to observations; it is a **discrete structural degree** specifying the finiteness, closure degree, and relational-space scale of the system itself.

Meanwhile, the quantities appearing in this paper and the three self-referenced papers,

$$
R,\quad t,\quad Q,\quad
r,\quad
P_{\rm odd},\quad P_{\rm even},
$$

the principal-axis lengths, the clock field, the mass readout, the reflectance, the covering degree, and the Bose/Fermi/mixed classification, are not entered as independent free parameters; they are treated as quantities arising from the zero closure, finite order, simplex relations, harmonic structure, self-consistency, and observation maps.

At the present stage, therefore, the model is characterized as a

$$
\boxed{
\text{one-discrete-parameter structural model}
}
$$

that is, **a structural model whose only externally specified degree of freedom is the integer $N$**.

This is not merely a matter of "having few parameters". It means that no individual masses, coupling constants, curvature scales, statistical mixing rates, or anisotropies are added as separate tuning knobs to obtain the symmetries, geometry, statistics, and readout structures the construction sets out to explain.

Furthermore, if it is shown in the future that $N$ itself is selected uniquely or discretely by the self-consistency or closure conditions, then

$$
\boxed{
\text{external free input}=0
}
$$

becomes possible. This, however, is not treated as derived at present.

### The Remaining Generalization Task Concerning $N$

In the present theory $N$ is the only externally specified integer. One of the remaining fundamental problems is therefore whether

$$
\boxed{
\text{why this universe / observation sector selects a particular }N
}
$$

can be derived from self-consistency, finite-order closure, or the structure of the large-$N$ limit.

This is not a task of adding a new continuous parameter to the derived structure; it is **the task of finding the selection rule for the single remaining discrete input $N$**.


### Derived, Constructed, and Numerically Confirmed in the Self-Papers

- **The 3-dimensional ellipsoidal structure and the dominance of three principal axes**  
  Three-dimensionality is classified as already derived as a 3-dimensional ellipsoidal structure, from the equal-norm orthogonal structure arising from the zero closure, the orthogonal two-planes with independent normal, and the consistency of phase closure and simplex geometry. Furthermore, "Zero Closure Was Four-Dimensional" numerically confirms that of at most $N-1$ nontrivial principal axes the top three directions $A,B,C$ amplify selectively, the middle directions remain nearly stationary, and some lower directions migrate to imaginary directions. The establishment of three directions, the anisotropy, and the dominance and saturation in three directions are themselves not "underived".  
  **Remaining task:** to elevate the coincidence of geometric three-dimensionality and dynamical three-axis dominance to a general theorem for arbitrary sufficiently large $N$ and wide initial conditions.



The following are not items newly proved within this paper from A1--A4 alone, but are already derived, constructed, and numerically confirmed as self-hypotheses in the author's three published papers. They are therefore not classified as "underived".

- **Spontaneous generation of three directions**  
  The construction in which, from a self-consistent coherent parent via `make_parent`, a fermionic seed with odd harmonics, zero closure, and selective, inflation-like amplification, three anisotropic visible principal axes come to dominate among many relational directions, has been realized.  
  **Remaining task:** the analytic elucidation of why three directions are selected, the extraction of necessary and sufficient conditions, and generalization to arbitrary sufficiently large $N$.

- **Bose/Fermi/mixed (Ermion) classification**  
  The $Z_2$ parity under half-period translation of odd and even harmonics, the reflectance determined by the odd-harmonic ratio
  $$
  r=\frac{P_{\rm odd}}{P_{\rm odd}+P_{\rm even}},
  $$
  and the Bose/Fermi/mixed classification through single/double covers and the antipodal two-point structure are derived, constructed, and numerically confirmed in "The Periodic Table of Waves v2".  
  **Remaining task:** complete mass identification with real particle species, a closed derivation of the mass hierarchy, and the establishment of a general correspondence with the spin-statistics structure of standard quantum field theory.

- **The two-layer structure of mass, clock field, and gravity readout**  
  The construction separating the wave-existence layer from the observational readout layer, reading mass as the value of the clock field and gravity as its gradient, and its separation from gauge-type readout, are derived and numerically confirmed as self-hypotheses in "Two-Layer Separation of Waves and Fields".  
  **Remaining task:** fixing the mass invariant, generalization to the real mass hierarchy, extension of the gravity readout to general variable-curvature systems, and full dynamicization.

- **The readout-level distinction of $t,R,Q$**  
  $t,R,Q$ are fundamentally unobservable complex axes of the same kind; what is read as state updating is named time $t$, what is read as the curvature scale of visible geometry is named $R$, and what is read as an internal quantity is named $Q$. This is an identification on observation maps, not a "problem of deriving three kinds of axes from the axioms".  
  **Remaining task:** not a proof of principled role differentiation, but refining, as needed, the correspondence between each observation map and measured quantities.

### Not Yet Derived or Generalized, Even Including This Paper and the Self-Papers

Only what truly remains outside the closed derivation chain is placed here.

- **General classification of readout sectors from the fundamental $M$-component system**


  Apart from the three-direction generation by the existing `make_parent` construction, to obtain the selection of a specific number of axes as a general theorem from A1--A4 alone.

- **Complete dynamics of general variable curvature**  
  Not the local curvature structure itself, which reads $R$ as an unobservable direction of the complex square closure, but closed, complete dynamics for arbitrarily varying curvature fields.

- **Complete derivation of local gauge dynamics**  
  Beyond the gauge-type readouts and selection rules already obtained, to close local gauge connections and their complete dynamics generally from the axiom system.

- **Complete derivation of chirality / hypercharge / anomaly cancellation**  
  To close uniquely the correspondence candidates with the derived winding-number, parity, and internal-quantity structures, up to the chirality, hypercharge, and anomaly cancellation conditions of the Standard Model.




- **The selection rule for $N$**  
  To derive whether the single remaining externally specified structural parameter $N$ is selected autonomously from self-consistency, finite-order closure, the large-$N$ limit, or similar. If this closes, the model has zero external free input.

## Self-References

- Noriaki Kihara, "The Periodic Table of Waves v2 — Particle Classification by Winding-Number Address and Observation Clock, and the Unification of Mass, Lifetime, and Splitting via the Clock Field $\omega(x)$", Zenodo, Concept DOI: 10.5281/zenodo.21830706. **Used in this paper as the source of previously derived results as self-hypotheses.**
- Noriaki Kihara, "Two-Layer Separation of Waves and Fields — Unification of Gauge Fields and Gravitational Fields via a Universal Field-Readout Function", Zenodo, DOI: 10.5281/zenodo.21832257. **Used in this paper as the source of previously derived results as self-hypotheses.**
- Noriaki Kihara, "Zero Closure Was Four-Dimensional — 'Central Projection' Survives Even in the Complex World", Zenodo, Version DOI: 10.5281/zenodo.21902806, Concept DOI: 10.5281/zenodo.21902805. **Used in this paper as the source of the previously derived results, as self-hypotheses, on three-direction spectral concentration, anisotropic principal-axis generation, and inflation-like expansion.**
