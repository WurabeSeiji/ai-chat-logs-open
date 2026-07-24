# An Axiomatic Framework for Physical Symmetry Selection by Anonymous Partial Projection
## Conditional Identification of $E_8$ in an Eight-Component Two-Quartet Branch

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** July 24, 2026<br>
**Version DOI:** [10.5281/zenodo.21521900](https://doi.org/10.5281/zenodo.21521900)<br>
**Concept DOI:** [10.5281/zenodo.21521899](https://doi.org/10.5281/zenodo.21521899)<br>
**Status:** Independent paper, English version v1 / Mathematical physics and foundations of physics

---

## Abstract

When the solution set of a foundational axiom system contains more than one closed system, the symmetry observed in physics need not be identified with the symmetry of the entire solution set. This paper proposes an **anonymous partial-projection selection principle**, in which physical state space is defined as a selected partial-projection image of the candidate space allowed by the foundational axioms.

Let $\mathscr S$ be the candidate state space, and let $\mathfrak P(\mathscr S)$ be the set of admissible partial projections satisfying anonymity, scale anonymity, closure consistency, and readout consistency. The central axiom of this paper is

$$
\boxed{
\exists\mathcal P_*\in\mathfrak P(\mathscr S)
\quad\text{such that}\quad
\mathscr S_{\mathrm{phys}}
{}=
\operatorname{Im}\mathcal P_*.
}
$$

This axiom does not specify which symmetry group is to be selected. Physical symmetry is identified only after projection, as the automorphism group of the relational structure carried by the selected image.

As a nontrivial application, consider a branch in which transition differences in the selected image lie in an eight-dimensional positive-definite readout space, split into two four-component systems, with each local transition lattice of type $D_4$. Identify the two $D_4$ discriminant groups with a single common center and retain only the neutral sector under its diagonal action. The permitted glue classes are then restricted to

$$
(0,0),\qquad(v,v),\qquad(s,s),\qquad(c,c).
$$

The lattice obtained by this diagonal gluing is positive definite, eight-dimensional, even, and unimodular; it is therefore isomorphic to the $E_8$ lattice. Its minimal vectors number

$$
48+3\times64=240,
$$

and the dimension of the associated complex simple Lie algebra is

$$
8+240=248.
$$

Furthermore, a Coxeter element of $E_8$ has order 30.

The paper does not claim that the foundational axioms uniquely force $E_8$. The physical reasons for selecting eight components, the two-quartet $D_4$ structure, and the common-center projection remain selection problems. What is established is the axiomatic order: anonymous partial projection is formulated first, and $E_8$ appears only as the mathematical classification of one explicitly conditioned branch.

---

## 0. Position and Scope

### 0.1 What this paper addresses

This paper does three things.

1. It separates the admissible solution space of the foundational axioms from the physically realized state space.
2. It formulates physical realization as the image of an anonymous admissible partial projection.
3. It proves that, when a common-center neutrality projection is imposed on an eight-component two-quartet branch, the transition lattice of the selected image is $E_8$.

### 0.2 What this paper does not claim

This paper does not claim:

1. that the foundational axioms alone uniquely derive $E_8$;
2. that every physical partial-projection image has $E_8$ symmetry;
3. that a dynamics or selection functional determining the actual projection $\mathcal P_*$ has already been obtained;
4. that the eight readout quantities have already been mapped to particles, interactions, or quantities of the Standard Model;
5. that a finite-recurrence operator is physically identical to a Coxeter element of $E_8$; or
6. that $E_8$ is the symmetry of spacetime as a whole or of the entire state space.

### 0.3 Central claim

The central claim is an ordering principle:

$$
\boxed{
\begin{aligned}
&\text{foundational axioms}\\
&\Longrightarrow
\text{candidate space containing multiple admissible closed systems}\\
&\Longrightarrow
\text{set of anonymous admissible partial projections}\\
&\Longrightarrow
\text{selection of one physical partial image}\\
&\Longrightarrow
\text{mathematical identification of the symmetry of that image}.
\end{aligned}
}
$$

$E_8$ is not inserted as an axiom at the beginning of this chain. It is an identification obtained at the end of one particular branch.

---

## 1. Connection to the Foundational Axioms

### 1.1 Inherited formation rules

This paper inherits the following formation rules from the *Basic Axiom System for an Anonymous Equal-Amplitude Composite-Wave Model, v8* [1].

**Anonymity.** No elementary component is assigned an individual name, privileged axis, privileged sign, or privileged type. A physical name or the name of an existing theory may not be used as a reason to alter an equation or readout rule at the first-principles level.

**Scale anonymity.** The states $X$ and $\lambda X$, for $\lambda\in\mathbb C^\ast$, are equivalent as relational structures after removal of absolute scale.

**All-positive-sign zero closure:**

$$
Q(X):=\sum_{j=1}^{N}x_j^2=0.
$$

**Nontrivial existence:**

$$
X\neq0.
$$

### 1.2 What is added here

The inherited formation rules restrict candidate states, but they do not uniquely determine which candidates are physically readable. This paper fills that logical gap by formulating the **existence and admissibility of a partial projection**.

No particular symmetry group is added as an axiom. In particular,

$$
\text{“physical symmetry is }E_8\text{”}
$$

is not postulated. Such a postulate would violate anonymity by using the desired name of the resulting theory to select the readout rule.

---

## 2. Admissible Solution Space

### 2.1 Projectivized candidate states

Define scale equivalence by

$$
X\sim_{\mathrm{sc}}\lambda X,
\qquad
\lambda\in\mathbb C^\ast.
$$

The projectivized set of nontrivial zero-closure candidates is

$$
\mathscr S_0
:=
\left\{
[X]
\ \middle|\
X\neq0,\ Q(X)=0
\right\},
$$

where

$$
[X]
:=
\left\{
\lambda X
\mid
\lambda\in\mathbb C^\ast
\right\}.
$$

### 2.2 Candidate systems with evolution

Let $U$ be a complex-linear evolution operator on the state space. Linearity makes $U$ well defined on scale-equivalence classes. Compatibility with the same closed system is expressed by

$$
Q(UY)=Q(Y)
\qquad
(\forall Y).
$$

Define the candidate systems with evolution by

$$
\mathscr S
:=
\left\{
([X],U)
\ \middle|\
[X]\in\mathscr S_0,\ Q\circ U=Q
\right\}.
$$

### 2.3 Finite-recurrence branch

Candidate systems with finite recurrence form

$$
\mathscr S_{\mathrm{fin}}
:=
\left\{
([X],U,n)
\ \middle|\
([X],U)\in\mathscr S,\quad
U^n=I,\quad
n\in\mathbb N
\right\}.
$$

Finite recurrence is not required for the general partial-projection principle. It is an additional condition selecting the subfamily that may later be compared with a Coxeter period of $E_8$.

---

## 3. Anonymous Partial Projection

### 3.1 Definitions

**Definition 3.1 (partial projection).** A partial projection of the candidate space $\mathscr S$ is a pair

$$
\mathcal P
=
\left(
\mathscr D_{\mathcal P},
p_{\mathcal P}
\right),
$$

where $\mathscr D_{\mathcal P}\subseteq\mathscr S$ and

$$
p_{\mathcal P}:
\mathscr D_{\mathcal P}
\longrightarrow
\mathscr Y_{\mathcal P}.
$$

Here “partial projection” is not restricted to an orthogonal projection in linear algebra. It may include:

1. selecting only part of the candidate systems as its domain;
2. quotienting differences that are unreadable; and
3. mapping readable relational quantities into another representation space.

**Definition 3.2 (nontrivial partial projection).** A partial projection is nontrivial if at least one of the following holds:

$$
\mathscr D_{\mathcal P}
\subsetneq
\mathscr S,
$$

or

$$
\exists s_1\neq s_2
\quad\text{such that}\quad
p_{\mathcal P}(s_1)
=
p_{\mathcal P}(s_2).
$$

Thus either some candidates are removed, or some differences between candidates are identified as unreadable. The identity map alone is not counted as a nontrivial realization of the principle.

### 3.2 Transition-difference lattice

Assume that the selected image contains discrete transition differences. Let $\Delta_{\mathcal P}\subseteq V_{\mathcal P}$ be their set, and define its integer span by

$$
\Lambda_{\mathcal P}
:=
\operatorname{span}_{\mathbb Z}
\Delta_{\mathcal P}.
$$

$\Lambda_{\mathcal P}$ stores readable state differences or transition differences, not individual names of states.

### 3.3 Admissibility conditions

**Definition 3.3 (admissible partial projection).** A partial projection is admissible if it is nontrivial and satisfies the following conditions.

#### A. Component anonymity

For any permutation or anonymity-preserving transformation $g$, the projection result is independent of the relabeling of individual components:

$$
p_{\mathcal P}(g\cdot[X])
\sim_{\mathscr Y}
p_{\mathcal P}([X]).
$$

#### B. Scale anonymity

$$
p_{\mathcal P}([\lambda X])
=
p_{\mathcal P}([X]),
\qquad
\lambda\in\mathbb C^\ast.
$$

#### C. Closure consistency

The projection domain consists of candidates satisfying the foundational closure, and the projected transition differences are closed under the composition rule defined on the image.

#### D. Readout consistency

Two states identified by the projection cannot be distinguished using only the adopted readout quantities.

#### E. Prohibition of prior physical naming

Neither a physical name nor a symmetry-group name that becomes available only after projection may be used in the definition of the projection.

Write the set of all admissible partial projections as

$$
\mathfrak P(\mathscr S).
$$

---

## 4. Partial-Projection Existence Axiom

### 4.1 The axiom

**Axiom PP (anonymous partial-projection existence axiom).** For the candidate space $\mathscr S$ of the foundational axioms, the physically realized state space is the image of at least one nontrivial admissible partial projection:

$$
\boxed{
\exists\mathcal P_*
\in
\mathfrak P(\mathscr S)
\quad\text{such that}\quad
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*.
}
$$

### 4.2 What the axiom does not assert

Axiom PP does not assert that

$$
\mathcal P_*
\text{ is unique}.
$$

Nor does it specify a selection functional such as

$$
\mathcal P_*
=
\operatorname*{arg\,max}_{\mathcal P\in\mathfrak P(\mathscr S)}
\mathcal F(\mathcal P).
$$

It is therefore an existence axiom saying that physical realization is a partial-projection image. It is not yet a selection dynamics determining which image is realized.

### 4.3 Physical symmetry

**Definition 4.1 (physical symmetry of a projection image).** Define the physical symmetry group to be the automorphism group preserving the readout relations, transition-difference set, and composition rule of the selected image:

$$
G_{\mathrm{phys}}
:=
\operatorname{Aut}
\left(
\operatorname{Im}\mathcal P_*,
\Delta_{\mathcal P_*}
\right).
$$

For a description by a discrete transition lattice, write

$$
G_{\mathrm{lattice}}
:=
\operatorname{Aut}
\left(
\Lambda_{\mathcal P_*},
(\,\cdot\,,\,\cdot\,)
\right),
$$

where a lattice automorphism must preserve both the lattice and its inner product.

Physical symmetry is not a name assigned in advance to the full candidate space. It is identified after projection as the group preserving the relations that survive in the selected image.

### 4.4 Nonunique branches

Axiom PP does not exclude two admissible projections

$$
\mathcal P_1,\mathcal P_2
\in
\mathfrak P(\mathscr S)
$$

for which

$$
\operatorname{Im}\mathcal P_1
\not\cong
\operatorname{Im}\mathcal P_2.
$$

Different stability conditions, readout conditions, or boundary conditions may therefore realize different symmetry branches.

---

## 5. Eight-Component Two-Quartet Branch

We now impose additional conditions on one branch allowed by Axiom PP. They are not consequences of PP; they are explicit branch assumptions defining the $E_8$ branch.

### 5.1 Branch assumption Q1: positive-definite eight-component readout

Assume that the transition differences of the selected image lie in an eight-dimensional real readout space with a positive-definite inner product,

$$
V_{\mathrm{read}}
\cong
\mathbb R^8,
$$

and that

$$
\operatorname{rank}
\Lambda_{\mathcal P_*}
=8.
$$

The complex bilinear form inherited from Axiom 1,

$$
Q(X)=\sum_jx_j^2,
$$

and the positive-definite readout inner product

$$
(\alpha,\beta)
$$

introduced here are different structures. Positive definiteness is not derived directly from $Q(X)=0$.

### 5.2 Branch assumption Q2: two-quartet decomposition

Assume an orthogonal decomposition of the readout space into two four-component subspaces:

$$
V_{\mathrm{read}}
=
V_A\oplus V_B,
\qquad
\dim V_A
=
\dim V_B
=4.
$$

At this stage no component is named space, time, mass, charge, or any other physical quantity.

### 5.3 Branch assumption Q3: local $D_4$ transition lattices

Assume that the local transition lattice of each four-component subsystem is isomorphic to the lattice defined by

$$
D_4
:=
\left\{
z\in\mathbb Z^4
\ \middle|\
\sum_{i=1}^{4}z_i
\equiv0
\pmod2
\right\}.
$$

The symbol $D_4$ is used as its mathematical classification name. The content of the branch assumption is the displayed coordinate condition; the lattice is not selected because a symmetry name was desired.

Before gluing, the local lattice is therefore

$$
L_0
:=
D_4^{(A)}
\oplus
D_4^{(B)}.
$$

This local $D_4$ structure has not been derived from the foundational axioms or Axiom PP alone.

### 5.4 Discriminant classes of $D_4$

Let $D_4^*$ be the dual lattice of $D_4$. Its discriminant group is

$$
\mathcal A
:=
D_4^*/D_4
\cong
(\mathbb Z/2\mathbb Z)^2,
$$

with four classes

$$
\mathcal A
=
\{0,v,s,c\}.
$$

Using a standard orthonormal basis $e_1,\ldots,e_4$, representatives may be chosen as

$$
v=e_1+D_4,
$$

$$
s=
\frac12
(e_1+e_2+e_3+e_4)
+D_4,
$$

and

$$
c=
\frac12
(-e_1+e_2+e_3+e_4)
+D_4.
$$

The discriminant quadratic form

$$
q:
\mathcal A
\longrightarrow
\mathbb Q/2\mathbb Z,
\qquad
q(a)
=
(a,a)
\bmod
2\mathbb Z
$$

satisfies

$$
q(v)
=
q(s)
=
q(c)
=1
\pmod{2\mathbb Z}.
$$

---

## 6. Common-Center Neutrality Projection

### 6.1 Branch assumption Q4: common-center action

Assume that the discriminant groups of the two four-component systems are not retained as independent label groups. Instead, identify them with one common center label group

$$
\mathcal A
\cong
(\mathbb Z/2\mathbb Z)^2,
$$

on which the same element $z\in\mathcal A$ acts. Only the finite Abelian group structure of $\mathcal A$ and its discriminant pairing are used here.

The nondegenerate pairing induced by the discriminant form identifies $\mathcal A$ with its character group $\widehat{\mathcal A}$. Let $\chi_a$ denote the central character corresponding to $a\in\mathcal A$. Every character takes real values $\pm1$. Write the two discriminant classes as

$$
(a,b)
\in
\mathcal A\oplus\mathcal A.
$$

### 6.2 Branch assumption Q5: selection of the neutral sector

Assume that only components invariant under the common-center action remain in the physical readout image. Define the neutrality projector by

$$
\Pi_{\mathrm{com}}(a,b)
:=
\frac{1}{|\mathcal A|}
\sum_{z\in\mathcal A}
\chi_a(z)\chi_b(z).
$$

The orthogonality relation for characters of a finite Abelian group gives

$$
\Pi_{\mathrm{com}}(a,b)
=
\delta_{a+b,\,0}.
$$

### 6.3 Diagonal-class selection theorem

**Theorem 6.1 (diagonal-class selection by a common center).** Under branch assumptions Q4 and Q5, the only discriminant classes passing the common-center neutrality projection are

$$
\boxed{
(0,0),
\qquad
(v,v),
\qquad
(s,s),
\qquad
(c,c).
}
$$

**Proof.** Neutrality requires

$$
a+b=0.
$$

Every element of $\mathcal A\cong(\mathbb Z/2\mathbb Z)^2$ has order one or two, so

$$
-a=a.
$$

Hence

$$
b=-a=a.
$$

Substitution of $\mathcal A=\{0,v,s,c\}$ leaves precisely the four displayed classes. $\square$

### 6.4 Diagonal glue group

Define the admissible glue group by

$$
H_{\mathrm{diag}}
:=
\left\{
(0,0),
(v,v),
(s,s),
(c,c)
\right\}.
$$

A graph gluing obtained by simultaneously relabeling $v,s,c$ on one side by triality can be written in the same diagonal form after the identification of the two discriminant groups is relabeled. Thus “equal classes” is the standard form after both discriminant groups have been identified with the same physical center labels.

---

## 7. Conditional Identification of the $E_8$ Lattice

### 7.1 Glued lattice

Define the lattice after diagonal gluing by

$$
\Lambda
:=
\bigcup_{h\in H_{\mathrm{diag}}}
\left(
L_0+h
\right).
$$

The general theory of lattice gluing and discriminant forms follows Nikulin [2]. The index, determinant, and parity needed here are also calculated directly below.

### 7.2 Index and determinant

Since the determinant of $D_4$ is 4,

$$
\det L_0
=
\det(D_4\oplus D_4)
=
4^2
=16.
$$

The group $H_{\mathrm{diag}}$ has order 4, so

$$
[\Lambda:L_0]
=4.
$$

The determinant formula for a finite-index overlattice gives

$$
\det\Lambda
=
\frac{\det L_0}
{[\Lambda:L_0]^2}
=
\frac{16}{4^2}
=1.
$$

Thus $\Lambda$ is unimodular.

### 7.3 Evenness

Let $q\oplus q$ be the direct-sum discriminant form. For every nonzero diagonal class,

$$
(q\oplus q)(a,a)
=
q(a)+q(a)
=
1+1
=0
\pmod{2\mathbb Z},
$$

where

$$
a\in\{v,s,c\}.
$$

Therefore $H_{\mathrm{diag}}$ is isotropic with respect to the discriminant form, and the glued lattice $\Lambda$ is even.

### 7.4 Main theorem

**Theorem 7.1 ($E_8$ identification in the eight-component two-quartet branch).** Under branch assumptions Q1--Q5, the transition-difference lattice $\Lambda$ obtained by the common-center neutrality projection is isomorphic to the $E_8$ lattice:

$$
\boxed{
\Lambda
\cong
\Lambda_{E_8}.
}
$$

**Proof.** By Q1 and Q2, $\Lambda$ is positive definite and eight-dimensional. By Q3, the lattice before gluing is $D_4\oplus D_4$. By Theorem 6.1, the glue group is $H_{\mathrm{diag}}$. Section 7.2 gives $\det\Lambda=1$, and Section 7.3 shows that $\Lambda$ is even. Thus $\Lambda$ is a positive-definite, rank-eight, even unimodular lattice. In dimension eight, this lattice is unique up to isomorphism and is the $E_8$ lattice [3]. Hence $\Lambda\cong\Lambda_{E_8}$. $\square$

### 7.5 Logical status

Theorem 7.1 is not the unconditional statement

$$
Q(X)=0
\Longrightarrow
E_8.
$$

Its exact content is

$$
\boxed{
\begin{array}{c}
\text{Axiom PP}\\
{}+\text{Q1: positive-definite rank eight}\\
{}+\text{Q2: two quartets}\\
{}+\text{Q3: local }D_4\\
{}+\text{Q4--Q5: common-center neutrality projection}
\end{array}
\Longrightarrow
\Lambda_{E_8}.
}
$$

$E_8$ is therefore not a name inserted into the axioms. It is the classification of a projection image satisfying explicit conditions.

---

## 8. The 240 Roots and Dimension 248

### 8.1 Minimal transition vectors

Let the norm-two vectors of the $E_8$ lattice be

$$
\Phi
:=
\left\{
\alpha\in\Lambda
\mid
(\alpha,\alpha)=2
\right\}.
$$

### 8.2 Vectors inside $D_4\oplus D_4$

The roots of $D_4$ are

$$
\pm e_i\pm e_j,
\qquad
1\le i<j\le4,
$$

and their number is

$$
4\binom42
=24.
$$

Hence the number of norm-two vectors inside

$$
D_4\oplus D_4
$$

is

$$
24+24=48.
$$

### 8.3 Nonzero glue classes

Each of the discriminant classes $v,s,c$ has eight minimal representatives of norm 1. Each diagonal class

$$
(v,v),\qquad(s,s),\qquad(c,c)
$$

therefore contributes

$$
8\times8=64
$$

norm-two vectors. Since there are three nonzero diagonal classes,

$$
\boxed{
|\Phi|
=
48+3\times64
=240.
}
$$

### 8.4 Closure under reflections

For each $\alpha\in\Phi$, define

$$
s_\alpha(x)
:=
x
-\frac{2(x,\alpha)}
{(\alpha,\alpha)}
\alpha
=
x-(x,\alpha)\alpha.
$$

Because $\Lambda$ is integral and $(\alpha,\alpha)=2$,

$$
s_\alpha(\Lambda)
=
\Lambda.
$$

The group generated by these reflections is the Weyl group $W(E_8)$ of the $E_8$ root system.

### 8.5 Lie algebra

Using an eight-dimensional Cartan subalgebra $\mathfrak h$ and a one-dimensional root space $\mathfrak g_\alpha$ for every root,

$$
\mathfrak e_8(\mathbb C)
=
\mathfrak h
\oplus
\bigoplus_{\alpha\in\Phi}
\mathfrak g_\alpha.
$$

Consequently,

$$
\boxed{
\dim\mathfrak e_8
=
8+240
=248.
}
$$

The reflection symmetry of the discrete transition lattice is $W(E_8)$, while the complex simple Lie algebra containing the root spaces is $\mathfrak e_8(\mathbb C)$. If the compact real form corresponding to the positive-definite inner product is chosen, the associated simply connected Lie group is compact $E_8$.

---

## 9. Connection to Finite Recurrence

### 9.1 Coxeter element

The Coxeter number of the $E_8$ root system is 30, and a Coxeter element $C$ satisfies [4]

$$
\boxed{
C^{30}=I.
}
$$

Its eigenvalues on the Cartan space are

$$
\operatorname{Spec}(C)
=
\left\{
\exp
\left(
\frac{2\pi i m}{30}
\right)
\ \middle|\
m\in
\{1,7,11,13,17,19,23,29\}
\right\}.
$$

Under the Coxeter action, the 240 roots split into eight orbits of length 30 [4]:

$$
\Phi
=
\bigsqcup_{a=1}^{8}
\mathcal O_a,
\qquad
|\mathcal O_a|
=30.
$$

Thus

$$
\boxed{
248
=
8+240
=
8+8\times30
=
8(30+1).
}
$$

### 9.2 Logical relation to a finite-recurrence axiom

An $E_8$ projection image contains an order-30 finite-recurrence operator $C$. It does not follow from Theorem 7.1, however, that a general physical evolution operator $U$ satisfying

$$
U^n=I
$$

is identical to $C$:

$$
U=C.
$$

Identifying

$$
U_{\mathrm{phys}}
\longleftrightarrow
C
$$

requires an additional correspondence showing that physical evolution cycles through the complete root system along Coxeter orbits.

The result established here is only

$$
\boxed{
\Lambda_{\mathcal P_*}
\cong
\Lambda_{E_8}
\Longrightarrow
\exists C\in\operatorname{Aut}(\Lambda_{\mathcal P_*})
\text{ such that }
C^{30}=I.
}
$$

---

## 10. Consistency with Anonymity

### 10.1 No symmetry-group name in the projection condition

Axiom PP and the admissibility conditions contain no reference to the name $E_8$. The branch assumptions are stated as:

1. eight components;
2. two four-component systems;
3. an even-integer condition on local transition differences; and
4. a common-center neutrality projection.

The resulting lattice is then identified as $E_8$ by the known classification of lattices. The logical order is therefore

$$
\boxed{
\text{projection conditions}
\longrightarrow
\text{lattice structure}
\longrightarrow
\text{symmetry-group name}.
}
$$

### 10.2 Physical names are also assigned after projection

A later interpretation of the eight components as, for example,

$$
(x,y,z,R)
\oplus
(t,Q_1,Q_2,Q_3)
$$

is not part of the theorem. Such a correspondence is a physical interpretation that can be made only after readout operations, conserved quantities, exchange rules, and experimental correspondences have been specified for the components.

### 10.3 Branches other than $E_8$

If any of Q1--Q5 fails, this paper does not require $E_8$.

- If the readout rank is not eight, Theorem 7.1 does not apply.
- If the local lattice is not of type $D_4$, a different gluing problem results.
- Without diagonal gluing, $D_4\oplus D_4$ remains a decomposable lattice of determinant 16.
- Without positive definiteness, uniqueness of the positive-definite $E_8$ lattice does not apply.

The framework therefore permits other partial-projection images and other symmetry branches from the outset.

---

## 11. The Selection Problem

### 11.1 The unresolved core

The unresolved part of the framework is not the lattice calculation after gluing. It is why an actual physical system selects a projection satisfying Q1--Q5.

What remains necessary is

$$
\boxed{
\text{an anonymous selection rule for }\mathcal P_*.
}
$$

### 11.2 Representation by a selection functional

If the selection rule can be expressed as a functional, introduce

$$
\mathcal F:
\mathfrak P(\mathscr S)
\longrightarrow
\mathbb R
$$

and write

$$
\mathcal P_*
\in
\operatorname*{arg\,ext}_{\mathcal P\in\mathfrak P(\mathscr S)}
\mathcal F(\mathcal P).
$$

Here $\operatorname*{arg\,ext}$ denotes a maximum, minimum, or stationarity condition. No specific form of $\mathcal F$ is assumed in this paper.

### 11.3 Requirements on the selection rule

At minimum, the selection rule must:

1. not refer to individual component names;
2. not refer to an absolute scale;
3. preserve closure;
4. not insert a symmetry-group name directly into its objective; and
5. make the selection of eight components, local $D_4$ structure, and common-center action comparable under one criterion.

### 11.4 Status of the common-center condition

The common-center neutrality projection is a sufficient condition for diagonal gluing. However,

$$
\text{why the same center acts on both quartets}
$$

and

$$
\text{why only the neutral sector is readable}
$$

must be derived from a selection dynamics or an observational construction.

The common-center condition is therefore currently a **candidate selection axiom** for the $E_8$ branch, not a derived consequence of the foundational axioms.

---

## 12. Discriminating Tests

The framework calls for the following calculations.

### 12.1 Readout rank

Test whether the rank of physically independent transition differences is

$$
\operatorname{rank}\Lambda_{\mathcal P_*}
=8.
$$

### 12.2 Local lattice

Test whether the minimal transition vectors of each four-component system close under the $D_4$ root condition

$$
\pm e_i\pm e_j.
$$

### 12.3 Correlation of discriminant classes

Test whether simultaneous occurrence of the discriminant classes in the two four-component systems satisfies

$$
(a,b)
\text{ occurs}
\quad\Longleftrightarrow\quad
a=b.
$$

If off-diagonal classes such as

$$
(v,s),\qquad(v,c),\qquad(s,c)
$$

remain stable, the common-center neutrality-projection hypothesis is rejected.

### 12.4 Coxeter period

Test whether all selected transition differences split into eight orbits of length 30:

$$
240
\stackrel{?}{=}
8\times30.
$$

Even if this test succeeds, identifying physical time evolution $U$ with the Coxeter element $C$ still requires an operator-level correspondence.

---

## 13. Form of Incorporation into the Basic Axiom System

The general principle to be incorporated into the basic axiom system is Axiom PP, not $E_8$ itself.

### 13.1 Proposed minimal axiom

> **Anonymous partial-projection existence axiom.** For the admissible solution space of the foundational axioms, the physically realized state space is the image of at least one nontrivial admissible partial projection satisfying anonymity, scale anonymity, closure consistency, and readout consistency.

$$
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*,
\qquad
\mathcal P_*
\in
\mathfrak P(\mathscr S).
$$

### 13.2 Items not to be incorporated directly

The following are not inserted directly into the foundational axioms:

1. the symmetry-group name $E_8$;
2. the derived numbers 240, 248, and 30; and
3. a declaration that the eight-component two-quartet branch is the unique realization.

They remain theorems and consequences of this paper under additional branch conditions.

### 13.3 Consistency with a single source of truth

The published basic axiom system v8 [1] remains unchanged. After DOI publication of the present paper, only the general Axiom PP is to be added to the next version based on v8. The order is

$$
\begin{aligned}
&\text{fix definitions and theorems in an independent paper}\\
&\longrightarrow\\
&\text{incorporate only the general axiom into the next basic axiom system}.
\end{aligned}
$$

---

## 14. Conclusion

This paper has presented an axiomatic framework in which physical symmetry is not assigned in advance to the entire solution space of the foundational axioms, but is identified after projection as an automorphism of an admissible partial-projection image.

The central axiom is

$$
\boxed{
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*,
\qquad
\mathcal P_*
\in
\mathfrak P(\mathscr S).
}
$$

It does not specify a symmetry group. It states only that a partial-projection image consistent with anonymity, scale anonymity, closure, and readout gives the physical realization.

For one particular branch, positive-definite rank-eight readout, a two-quartet decomposition, local $D_4$ transition lattices, and a common-center neutrality projection were assumed. The only surviving classes are then

$$
\boxed{
(0,0),
\qquad
(v,v),
\qquad
(s,s),
\qquad
(c,c),
}
$$

and the glued lattice satisfies

$$
\boxed{
\Lambda
\cong
\Lambda_{E_8}.
}
$$

It follows that

$$
\boxed{
|\Phi|
=240,
\qquad
\dim\mathfrak e_8
=248,
\qquad
C^{30}
=I.
}
$$

The logical position of $E_8$ in this paper is therefore

$$
\boxed{
\begin{aligned}
&\text{$E_8$ was not selected by name;}\\
&\text{it is the classification of the image selected}\\
&\text{by an anonymous partial projection.}
\end{aligned}
}
$$

The central remaining problem is to derive an anonymous selection functional or stability condition that selects the actual $\mathcal P_*$ from the set of admissible partial projections.

---

## References

[1] N. Kihara, “Basic Axiom System for an Anonymous Equal-Amplitude Composite-Wave Model, v8,” Zenodo, Version DOI: [10.5281/zenodo.21495422](https://doi.org/10.5281/zenodo.21495422), Concept DOI: [10.5281/zenodo.21315735](https://doi.org/10.5281/zenodo.21315735), 2026.

[2] V. V. Nikulin, “Integral symmetric bilinear forms and some of their applications,” *Mathematics of the USSR-Izvestiya*, vol. 14, no. 1, pp. 103–167, 1980. DOI: [10.1070/IM1980v014n01ABEH001060](https://doi.org/10.1070/IM1980v014n01ABEH001060).

[3] L. J. Mordell, “The definite quadratic forms in eight variables with determinant unity,” *Journal de Mathématiques Pures et Appliquées*, 9e série, vol. 17, pp. 41–46, 1938. [NUMDAM](https://www.numdam.org/item/JMPA_1938_9_17_1-4_41_0/).

[4] D. A. Richter, “Triacontagonal coordinates for the $E_8$ root system,” arXiv:0704.3091, 2007. [arXiv](https://arxiv.org/abs/0704.3091).

---

## Appendix A: Audit of Claim Status

| Item | Status | Treatment in this paper |
|---|---|---|
| Anonymity, scale anonymity, and zero closure | Existing axioms | Inherited from the basic axiom system v8 [1] |
| Physical realization is an admissible partial-projection image | New Axiom PP | Central axiom of this paper |
| Multiple admissible partial projections may exist | Structure allowed by PP | Uniqueness is not assumed |
| Rule selecting the actual $\mathcal P_*$ | Not derived | A selection functional or dynamics is required |
| Positive-definite rank-eight readout | Branch assumption Q1 | Not derived from the foundational axioms |
| Two-quartet decomposition | Branch assumption Q2 | Not derived from the foundational axioms |
| Local $D_4$ transition lattice | Branch assumption Q3 | Not derived from the foundational axioms |
| Common-center action | Branch assumption Q4 | Candidate selection axiom |
| Retention of the neutral sector only | Branch assumption Q5 | Candidate selection axiom |
| Only diagonal classes survive | Theorem 6.1 | Derived from Q4--Q5 |
| Even unimodular rank-eight lattice | Derived consequence | Derived from the glue index and discriminant form |
| $\Lambda\cong\Lambda_{E_8}$ | Conditional Theorem 7.1 | Holds under Q1--Q5 |
| 240 roots | Derived consequence | Direct count by glue classes |
| Dimension 248 | Derived consequence | Rank 8 plus 240 roots |
| $C^{30}=I$ | Known mathematical consequence after $E_8$ identification | Richter [4] |
| Identification of physical $U$ with Coxeter element $C$ | Connection problem | Requires an operator correspondence |
| Mapping of eight components to $xyzR,tQ_1Q_2Q_3$ | Physical interpretation, not derived | Not part of the theorem |

---

## Appendix B: Branch Audit

| Condition removed or changed | Result |
|---|---|
| Remove Axiom PP | No principle distinguishes candidate space from physically realized space |
| Remove Q1 | The positive-definite rank-eight classification theorem cannot be applied |
| Remove Q2 | The construction cannot be expressed as a $D_4\oplus D_4$ gluing |
| Remove Q3 | The discriminant classes $\{0,v,s,c\}$ and the starting point of diagonal gluing are lost |
| Remove Q4 | The two quartet discriminant classes may be selected independently |
| Remove Q5 | There is no projection principle excluding off-diagonal classes |
| Remove diagonal gluing | $D_4\oplus D_4$ retains determinant 16 and is not self-dual |
| Replace positive definiteness by indefinite signature | Uniqueness of the positive-definite $E_8$ lattice no longer applies |
| Do not assume physical $U=C$ | An order-30 action exists inside $E_8$, but its identity with physical time evolution remains undetermined |

---

## Appendix C: Minimal Derivation Chain

$$
\boxed{
\begin{aligned}
&\text{anonymity, scale anonymity, and zero closure}\\
&\Longrightarrow
\mathscr S\\
&\xrightarrow{\text{Axiom PP}}
\operatorname{Im}\mathcal P_*\\
&\xrightarrow{\text{Q1--Q3}}
D_4\oplus D_4\\
&\xrightarrow{\text{Q4--Q5}}
H_{\mathrm{diag}}
=
\{(0,0),(v,v),(s,s),(c,c)\}\\
&\Longrightarrow
\text{positive-definite rank-eight even unimodular lattice}\\
&\Longrightarrow
\Lambda_{E_8}\\
&\Longrightarrow
240\text{ roots}\\
&\Longrightarrow
\dim\mathfrak e_8
=248\\
&\Longrightarrow
\exists C,\quad C^{30}=I.
\end{aligned}
}
$$
