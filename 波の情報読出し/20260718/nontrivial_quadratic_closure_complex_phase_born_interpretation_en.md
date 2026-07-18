# Complex Phase Structure Emerging from Nontrivial Quadratic Closure

## A Foundational Connection to Finite-Order Exchange Systems and Discrete Born-Type Weights

**Version:** English complete manuscript v1  
**Date:** July 18, 2026  
**Author:** Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**Version DOI:** 10.5281/zenodo.21422506  
**Concept DOI:** 10.5281/zenodo.21422505  
**Position in the series:** Additional interpretive paper deriving from the quadratic-closure axiom in the Wave Information Readout series

---

## Abstract

This paper aims to provide a more foundational interpretation of the previously reported result that finite-order closure of an exchange-symmetric two-channel unitary operator yields the discrete Born-type weights

$$
R_{n,m}=\cos^2\left(\frac{\pi m}{n}\right),
\qquad
1-R_{n,m}=\sin^2\left(\frac{\pi m}{n}\right).
$$

“Closure” in this paper does not mean a thermodynamically closed system isolated from the outside world. It means **algebraic quadratic closure**: the conjugation-free algebraic quadratic form $Q_0(x)=\sum_kx_k^2$ closes nontrivially to zero.

The preceding paper treated exchange symmetry, unitarity, and the finite-order condition within standard complex linear algebra. Here, instead, we take as a minimal axiom over a commutative field $K$ containing the real field $\mathbb R$ the nontrivial quadratic-closure condition

$$
\boxed{
\sum_{k=1}^{N}x_k^2=0,
\qquad
(x_1,\ldots,x_N)\neq(0,\ldots,0),
\qquad
x_k\in K,
\qquad
\mathbb R\subseteq K
}.
$$

If every component were real, each term would be nonnegative and the equation would have no nontrivial solution. The solution vector $\mathbf{x}:=(x_1,\ldots,x_N)$ must therefore lie outside $\mathbb R^N$; equivalently, at least one component is nonreal. In the minimal two-component system

$$
x^2+y^2=0,
\qquad x\neq0,
$$

we have

$$
\left(\frac{y}{x}\right)^2=-1.
$$

Writing an element whose square is $-1$ as $i$ gives

$$
y=\pm ix.
$$

The imaginary unit is not introduced externally to describe phase; it is required as the minimal extension direction that permits a nontrivial sum of squares to close to zero.

We interpret this nonreal direction as a quarter turn. With continuity, reversibility, and preservation of a positive-definite readout quantity added as further conditions, it is extended to a general phase rotation and connected to independently imposed exchange symmetry, the finite-order condition

$$
\zeta^n=1,
$$

and the exchange-symmetric channel projections

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2}.
$$

For the finite-order root $\zeta=e^{-2\pi im/n}$, taking squared weights as a positive-definite readout recovers the previously reported discrete Born-type sequence.

The rigorous derivation in this paper extends through the result that a nontrivial zero sum of squares cannot be formed from real numbers alone and that the minimal two-component system requires a square-root-of-minus-one direction. The subsequent progression to complex phase, an exchange-symmetric unitary representation, and Born-type weights uses additional conditions: continuity, linearity, reversibility, exchange symmetry, preservation of a positive-definite readout, and squared readout. This paper is therefore not a derivation of the complete Born rule. It classifies as derivation, additional construction, and unresolved connection the relation between the nonreal direction required by nontrivial quadratic closure and the finite-order exchange theorem reported previously.

**Keywords:** quadratic closure, imaginary number, complex phase, exchange symmetry, finite-order recurrence, Born rule, squared weight, closed system, unitary operator

---

# 1. Background and Objective

## 1.1 The Problem

In standard quantum mechanics, states are represented as complex vectors, time evolution is described by unitary operators, and observation probabilities are given by squared amplitude magnitudes:

$$
P(A)=|\langle A|\psi\rangle|^2.
$$

For a two-state system,

$$
|\psi\rangle
=
\cos\phi\,|A\rangle
+
e^{i\chi}\sin\phi\,|B\rangle,
$$

we have

$$
P(A)=\cos^2\phi,
\qquad
P(B)=\sin^2\phi.
$$

If complex numbers, conjugation, a positive-definite inner product, unitarity, and squared magnitudes are all assumed together before the $\cos^2$ form is obtained, however, the origin of that square law may not have been explained. Once complex amplitudes and conjugate products have been adopted, two-state projections are already constrained geometrically to have a $\cos^2/\sin^2$ form.

The preceding paper, “Emergence of Discrete Born-Type Weights in Iterated Two-Channel Exchange Systems,” did not take the Born rule as a numerical search target. It discovered discrete Born-type weights associated with finite-order recurrence in a model family addressing wave-packet localization transfer, metastable two-state dynamics, weak readout, and strong-observation selection. Its central result was that, up to an overall phase, an exchange-symmetric two-channel unitary operator can be written

$$
U=P_s+\zeta P_a,
\qquad |\zeta|=1,
$$

and that when the finite-order condition

$$
U^n=I
$$

selects

$$
\zeta=e^{-2\pi im/n},
$$

the channel weight is

$$
\left|\frac{1+\zeta}{2}\right|^2
=
\cos^2\left(\frac{\pi m}{n}\right).
$$

That result explains which phases and which squared weights are selected by exact recurrence. A more anterior question nevertheless remains: why should complex phase, conjugation, and unitarity be used at all?

## 1.2 Question Addressed in This Paper

This paper asks:

> Without taking complex numbers, conjugation, and unitarity as the first axioms, how far can one construct the connection from a nontrivial quadratic-closure condition to a square-root-of-minus-one direction, phase rotation, exchange-symmetric finite-order structure, and discrete Born-type weights?

We distinguish four layers:

$$
\boxed{
\begin{gathered}
\text{nontrivial quadratic closure}
\Rightarrow\text{nonreal direction}
\quad[\text{derived}]\\
\text{nonreal direction}
+\{\text{continuity, reversibility, }Q_+\text{ preservation}\}\\
\longrightarrow\text{phase rotation and exchange representation}
\quad[\text{additional construction}]\\
\text{exchange symmetry}+\text{finite order}+\text{squared readout}\\
\Rightarrow\text{discrete Born-type weights}
\quad[\text{conditional derivation}]\\
Q_0=0\text{ and the dynamical connection to the exchange kernel}\\
[\text{connection problem C2}]
\end{gathered}
}.
$$

## 1.3 Division of Roles Between the Two Papers

This paper does not replace the preceding finite-order exchange theorem.

The role of the preceding paper is to derive discrete Born-type weights rigorously within standard complex linear algebra from

$$
[U,X]=0,
\qquad
U^\dagger U=I,
\qquad
U^n=I.
$$

The role of the present paper is to show why a nonreal direction and phase description are required by the preceding nontrivial closure condition

$$
\sum_kx_k^2=0,
$$

and thereby to give a foundational interpretation of the mathematical structure of the preceding paper.

The relation between the papers is thus

$$
\boxed{
\begin{gathered}
\text{this paper: derives the necessity of a nonreal direction}\\
\xrightarrow[\text{C2 remains unresolved}]{\text{additional construction}}\\
\text{preceding paper: derives Born-type weights from finite-order exchange}
\end{gathered}
}.
$$

The preceding paper is permanently available at Concept DOI https://doi.org/10.5281/zenodo.21422470.

---

# 2. Axiom of Nontrivial Quadratic Closure

## 2.1 Axiom

Let $K$ be a commutative field containing the real field $\mathbb R$, and let there be finitely many components over it,

$$
x_k\in K,
\qquad
\mathbb R\subseteq K
\qquad(k=1,\ldots,N).
$$

We impose

$$
\boxed{
\sum_{k=1}^{N}x_k^2=0
}.
$$

To exclude the trivial solution alone, we also require

$$
\boxed{
(x_1,\ldots,x_N)\neq(0,\ldots,0)
}.
$$

We call this the **nontrivial quadratic-closure condition**.

Axiom A1 in the preceding work corresponds to the sector $K=\mathbb C$. Here we return the starting point to $\mathbb R\subseteq K$ in order to audit whether the conclusion arose merely because the complex field had been fixed in advance. After showing that two-component nontrivial closure generates a square root of $-1$ within $K$, we recover the minimal subfield $\mathbb R(i)\cong\mathbb C$.

It is essential that the condition is not

$$
\sum_k|x_k|^2=0,
$$

but

$$
\sum_kx_k^2=0.
$$

The former is a positive-definite norm using complex conjugation and has no nontrivial solution. The latter is a conjugation-free algebraic sum of squares and may have nontrivial solutions when the number system is extended beyond the reals.

## 2.2 Triviality over the Reals

Suppose every component is real:

$$
x_k\in\mathbb R.
$$

Then

$$
x_k^2\ge0.
$$

For

$$
\sum_{k=1}^{N}x_k^2=0
$$

to hold, every term must therefore vanish:

$$
x_1=x_2=\cdots=x_N=0.
$$

**Proposition 2.1 (absence of nontrivial real solutions).**  
Let $N\ge1$. For $x_k\in\mathbb R$, if

$$
\sum_{k=1}^{N}x_k^2=0,
$$

then

$$
x_k=0
\qquad(k=1,\ldots,N).
$$

Thus, if nontrivial quadratic closure is required, the solution vector lies outside the real vector space $\mathbb R^N$:

$$
\boxed{
\sum_kx_k^2=0,
\quad
(x_1,\ldots,x_N)\neq(0,\ldots,0)
\quad\Longrightarrow\quad
\mathbf{x}:=(x_1,\ldots,x_N)\in K^N\setminus\mathbb R^N
}.
$$

In components, this is equivalent to

$$
\boxed{
\exists j:\ x_j\notin\mathbb R
}.
$$

## 2.3 Meaning of the Proposition

The proposition does not require the entire extension field $K$ to equal the complex field. $K$ may be an extension larger than $\mathbb C$.

The element $i:=y/x$ obtained from two-component nontrivial closure nevertheless satisfies $i^2=-1$, so $K$ contains the subfield

$$
\mathbb R(i).
$$

As a field over $\mathbb R$, this is isomorphic to $\mathbb C$. The minimal commutative field extension obtained by adjoining a square root of $-1$ is therefore the complex field, while the solution vector of nontrivial quadratic closure lies outside $\mathbb R^N$.

In this paper, then, an imaginary number is not introduced arbitrarily as a convenient notation for a wave function. It appears as

> an extension direction required to close to zero a nontrivial sum of squares that cannot close using positive real squares alone.

---

# 3. The Imaginary Direction in the Minimal Two-Component System

## 3.1 Two-Component Closure

Consider the minimal nontrivial system

$$
x^2+y^2=0.
$$

If $x\neq0$, then

$$
y^2=-x^2,
$$

and

$$
\left(\frac{y}{x}\right)^2=-1.
$$

Writing an element whose square is $-1$ as $i$,

$$
i^2=-1,
$$

we obtain

$$
\boxed{
y=\pm ix
}.
$$

**Proposition 3.1 (minimal two-component closure).**  
In a commutative field extension satisfying $x^2+y^2=0$ with $x\neq0$, the ratio $y/x$ is a square root of $-1$.

No complex conjugation occurs in this proposition.

If $K$ contains $\mathbb R$, the minimal subfield generated by $i:=y/x$ is

$$
\mathbb R(i)\cong\mathbb C.
$$

Thus, in two-component closure the complex direction appears not as an externally supplied phase assumption, but as the minimal field extension generated by a nontrivial solution.

## 3.2 Interpretation as a Quarter Turn

In the complex plane, multiplication of a quantity $x$ along the real axis by $i$ corresponds to a 90-degree rotation:

$$
x\longmapsto ix.
$$

Therefore,

$$
y=\pm ix
$$

means that the two components are not merely positive or negative values on the same line; they point in directions displaced by one quarter of a cycle.

We call this relation a **quarter-period closure phase difference**:

$$
\Delta\phi=\pm\frac{\pi}{2}.
$$

The phase difference was not assumed first to construct quadratic closure. Rather, the nontriviality of quadratic closure requires a square-root-of-minus-one direction, whose geometric representation is a quarter turn.

## 3.3 Sign and Opposite Directions

The two solutions

$$
y=+ix,
\qquad
y=-ix
$$

have opposite directions of phase rotation:

$$
+\frac{\pi}{2},
\qquad
-\frac{\pi}{2}.
$$

This opposition may be related to the structure later represented by complex conjugation.

We do not, however, claim to have derived the uniqueness of the conjugation operation from this fact. What is rigorously obtained is that a square root of $-1$ has two signs, providing quarter-period closure rotations in opposite directions.

---

# 4. Nonreal Directions in a Multicomponent System

## 4.1 General Closure Equation

Consider generally

$$
\sum_{k=1}^{N}x_k^2=0.
$$

Suppose at least one component $x_j$ is nonzero. Since the square of a real component is nonnegative, closure cannot be achieved using real components alone.

For example, for $a,b\in\mathbb R$, let

$$
x_1=a,
\qquad
x_2=b,
\qquad
x_3=i\sqrt{a^2+b^2}.
$$

Then

$$
a^2+b^2+
\left(i\sqrt{a^2+b^2}\right)^2
=0.
$$

This example shows that many registered real squares can be closed by the square of a single imaginary direction.

## 4.2 Necessary and Sufficient Conditions

The solution vector of nontrivial closure lies outside $\mathbb R^N$. In fixed closure coordinates, this is equivalent to at least one component being nonreal.

The following stronger statements, however, do not hold in general:

- that there is exactly one nonreal component,
- that a nonreal component is purely imaginary,
- that each component separates into either a real or purely imaginary component, or
- that the complex field is the only possible extension.

The weaker proposition sufficient for this paper is

$$
\boxed{
\text{every solution of nontrivial quadratic closure lies outside }\mathbb R^N
}.
$$

This motivates the introduction of nonreal degrees of freedom representable as phase or rotation.

## 4.3 Connection to Unnamedness

If no external name or absolute coordinate is assigned to each component $x_k$, and only the closure condition of the entire system is taken as fundamental, then distinctions between components appear not as absolute values but as relative allocations and phase relations of the registered squares.

From this viewpoint,

$$
\sum_kx_k^2=0
$$

does not fix individual components as independent entities; it specifies that all components constitute a single closed relational system.

---

# 5. Distinguishing Quadratic Closure from Positive-Definite Readout

## 5.1 Two Kinds of Quadratic Form

We clearly distinguish the following two forms.

### Quadratic form defining closure

$$
Q_0(x):=\sum_kx_k^2.
$$

It contains no conjugation.

The nontrivial closure condition is

$$
Q_0(x)=0,
\qquad (x_1,\ldots,x_N)\neq(0,\ldots,0).
$$

### Quadratic form for observation or magnitude

For complex components, define

$$
Q_+(x):=\sum_kx_k^*x_k
=\sum_k|x_k|^2.
$$

This form is positive definite:

$$
Q_+(x)\ge0.
$$

## 5.2 Do Not Confuse the Closure Equation with the Norm Equation

Nontrivial quadratic closure means

$$
Q_0(x)=0,
$$

not

$$
Q_+(x)=0.
$$

For a nontrivial state, one ordinarily has

$$
Q_0(x)=0,
\qquad
Q_+(x)>0.
$$

For the minimal two-component example

$$
(x,y)=(a,ia),
$$

we have

$$
Q_0=a^2+(ia)^2=0,
$$

but

$$
Q_+=|a|^2+|ia|^2=2|a|^2>0.
$$

This distinction means that a zero sum of squares need not be read as a state in which all physical quantities have vanished. What vanishes is the directional algebraic closure register; the positive-definite readout may remain nonzero.

## 5.3 Role of Conjugation

In the construction presented here, conjugation is not a prerequisite for the imaginary direction to exist.

The imaginary direction is already required by the nontriviality of

$$
x^2+y^2=0.
$$

Conjugation is introduced afterward as an operation that reverses the nonreal direction,

$$
i\longmapsto-i,
$$

and converts

$$
z^*z
$$

into a real, positive-definite readout.

Separating the derived part from the added readout construction, the logical order is therefore

$$
\boxed{
\text{nontrivial quadratic closure}
\Rightarrow
\text{necessity of a nonreal direction}
\qquad[\text{derived}]
}
$$

and

$$
\boxed{
\begin{gathered}
\text{nonreal direction}+\text{definition of the opposing reversal}\\
+\text{definition of positive-definite readout}
\Rightarrow Q_+\\
[\text{additional construction}]
\end{gathered}
}.
$$

This paper does not prove that the standard complex conjugation is uniquely fixed by the nontrivial closure condition alone.

---

# 6. Generalization to Phase Rotation

## 6.1 Additional Construction from a Quarter-Period Closure Phase to Continuous Phase

Minimal two-component closure requires the quarter-turn relation

$$
1,
\qquad
i.
$$

If continuity, associativity, reversibility, and preservation of a positive-definite readout are added as conditions and this rotation is extended to an iterable continuous group, one obtains phase rotations by

$$
e^{i\phi}.
$$

They satisfy

$$
e^{i\phi}e^{i\psi}=e^{i(\phi+\psi)}.
$$

The continuous phase group $U(1)$ is not a consequence of the nontrivial closure condition alone. We additionally require

1. continuity of the rotation,
2. associativity of composition,
3. existence of an inverse, and
4. preservation of the positive-definite readout,

and adopt the unit-circle phase as the resulting one-parameter preservation group.

## 6.2 Real Two-Dimensional Representation

Without using complex numbers, the same transformation can be written as the real two-dimensional rotation

$$
R(\phi)
=
\begin{pmatrix}
\cos\phi&-\sin\phi\\
\sin\phi&\cos\phi
\end{pmatrix}.
$$

Then

$$
R(\phi)^TR(\phi)=I,
$$

so the real two-dimensional norm is preserved.

Complex phase can therefore be regarded as a compressed representation of real two-dimensional rotation.

The essential element in this paper is not the imaginary symbol itself, but

> extending the nonreal direction required by quadratic closure into a reversible, norm-preserving rotational degree of freedom.

## 6.3 Finite-Order Phase

The condition that the rotation return to its initial state after finitely many iterations is

$$
R(\phi)^n=I,
$$

or, in complex notation,

$$
\zeta^n=1.
$$

Thus,

$$
\zeta=e^{-2\pi im/n}.
$$

This integer-ratio phase gives a discrete closed orbit within the continuous phase space.

---

# 7. Connection to an Exchange-Symmetric Two-Channel System

## 7.1 Two-Channel Exchange

Let the exchange operator for the A/B channels be

$$
X=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix}.
$$

The symmetric and antisymmetric projectors are

$$
P_s=\frac{I+X}{2},
\qquad
P_a=\frac{I-X}{2}.
$$

Imposing

$$
[U,X]=0
$$

on an exchange-symmetric preserving map $U$ prevents it from mixing the symmetric and antisymmetric subspaces.

## 7.2 Phase-Fixed Representation

If a unitary representation is adopted as a linear reversible map preserving the positive-definite readout, then, up to an overall phase,

$$
\boxed{
U=P_s+\zeta P_a,
\qquad |\zeta|=1
}
$$

can be written.

Here,

- the symmetric component is fixed, and
- only the antisymmetric component undergoes phase rotation.

This is a **constructive hypothesis** placing the nonreal direction required by nontrivial quadratic closure in the antisymmetric degree of freedom of the two-channel exchange system.

### Connection Problem C2: Dynamical Nonpreservation of the Quadratic-Closure Surface

We now state explicitly the connection problem that remains between nontrivial quadratic closure and exchange dynamics. Write the two-channel state as

$$
x=(A,B)^T,
\qquad
Q_0(x)=x^Tx=A^2+B^2.
$$

For the exchange kernel

$$
U=P_s+\zeta P_a,
$$

the identities $P_sP_a=0$, $P_s^T=P_s$, and $P_a^T=P_a$ give

$$
U^TU=P_s+\zeta^2P_a.
$$

Therefore,

$$
\boxed{
Q_0(Ux)
=
x^T(P_s+\zeta^2P_a)x
}.
$$

For the nontrivial point on the closure surface

$$
x_0=(1,i)^T,
\qquad
Q_0(x_0)=1+i^2=0,
$$

one finds

$$
\boxed{
Q_0(Ux_0)=i(1-\zeta^2)
}.
$$

Thus, if $\zeta^2\neq1$, $Ux_0$ leaves the closure surface. Moreover, this exchange kernel satisfies

$$
Q_0(Ux)=Q_0(x)
$$

for every state only when $\zeta^2=1$. The previously reported exchange kernel with general finite-order roots is therefore not a self-map of the nontrivial quadratic-closure surface $Q_0=0$.

This paper does not claim that the preceding exchange kernel directly preserves $Q_0=0$. We designate the construction of an exchange map preserving the closure surface, a projection rule back onto that surface, or an enlarged state space that preserves the closure condition as **connection problem C2**. This is an open problem concerning a dynamical connection between the first axiom and exchange dynamics; it is independent of the validity of the finite-order exchange theorem itself.

## 7.3 Reconstruction in the Channel Basis

Returning to the A/B basis gives

$$
U
=
\frac12
\begin{pmatrix}
1+\zeta&1-\zeta\\
1-\zeta&1+\zeta
\end{pmatrix},
$$

and hence

$$
\boxed{
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2}
}.
$$

This sum and difference are the interference between the fixed symmetric phase $1$ and the rotating antisymmetric phase $\zeta$.

---

# 8. Finite-Order Closure and Discrete Born-Type Weights

## 8.1 Finite-Order Condition

After phase fixing, the condition for the operator to recur exactly after $n$ iterations is

$$
U^n=I.
$$

This is equivalent to

$$
\zeta^n=1.
$$

If $n$ is the primitive order and $\gcd(m,n)=1$, then

$$
\zeta=e^{-2\pi im/n}.
$$

## 8.2 Squared Weights

For an A-basis input, the A/B channel amplitudes are

$$
r_{n,m}=\frac{1+e^{-2\pi im/n}}{2},
$$

$$
t_{n,m}=\frac{1-e^{-2\pi im/n}}{2}.
$$

Taking squared magnitudes as a positive-definite readout gives

$$
|r_{n,m}|^2
=
\frac{1+\cos(2\pi m/n)}{2}
=
\cos^2\left(\frac{\pi m}{n}\right),
$$

$$
|t_{n,m}|^2
=
\frac{1-\cos(2\pi m/n)}{2}
=
\sin^2\left(\frac{\pi m}{n}\right).
$$

Therefore,

$$
\boxed{
R_{n,m}
=
\cos^2\left(\frac{\pi m}{n}\right)
}
$$

and

$$
\boxed{
1-R_{n,m}
=
\sin^2\left(\frac{\pi m}{n}\right)
}
$$

are obtained.

**Theorem 8.1 (finite-order exchange weights under additional conditions).**  
Let $U$ be a two-channel linear map satisfying the following conditions.

1. $U$ is an invertible map preserving the positive-definite readout $Q_+$.
2. $U$ commutes with the exchange operator $X$.
3. The overall phase is fixed so that the eigenvalue of the symmetric sector is $1$.
4. $U^n=I$, and the primitive order of the antisymmetric sector is $n$.
5. The A/B channel weights are defined by a positive-definite squared readout of the amplitudes.

Then, for coprime integers $m,n$,

$$
U=P_s+e^{-2\pi im/n}P_a,
$$

and the A/B channel weights for an A-basis input are

$$
\boxed{
R_{n,m}=\cos^2\left(\frac{\pi m}{n}\right),
\qquad
1-R_{n,m}=\sin^2\left(\frac{\pi m}{n}\right)
}.
$$

*Proof.* Conditions 1–3 give $U=P_s+\zeta P_a$ with $|\zeta|=1$. Condition 4 gives $\zeta^n=1$, and hence $\zeta=e^{-2\pi im/n}$. Returning to the A/B basis gives the amplitudes $(1\pm\zeta)/2$, and the positive-definite squared readout in condition 5 gives the displayed result. $\square$

This theorem restates the finite-order exchange theorem of the preceding paper so that the present paper can be read independently. Nontrivial quadratic closure does not derive conditions 1–5 themselves; it supplies the nonreal direction required for the complex-phase representation. The dynamical connection between the nontrivial closure surface and $U$ is C2 and is not included among the assumptions of this theorem.

## 8.3 Foundational Interpretation

The preceding paper obtained the formulas above as a finite-order theorem for an exchange-symmetric unitary kernel.

Here we interpret their more anterior stage through nontrivial closure,

$$
\boxed{
\sum_kx_k^2=0
}.
$$

The logical sequence is therefore not a single chain of unconditional implications, but is classified into four layers:

$$
\boxed{
\begin{gathered}
\text{nontrivial quadratic closure}
\Rightarrow\text{a square-root-of-minus-one direction}
\quad[\text{derived}]\\
\text{nonreal direction}
+\{\text{continuity, reversibility, }Q_+\text{ preservation}\}\\
\longrightarrow U(1)\text{ phase rotation}
\quad[\text{additional construction}]\\
U(1)\text{ phase}+\text{exchange symmetry}+\text{finite order}\\
+\text{squared readout}
\Rightarrow\text{discrete Born-type weights}
\quad[\text{conditional derivation}]\\
Q_0=0\text{ and the dynamical connection to the exchange kernel }U\\
[\text{connection problem C2}]
\end{gathered}
}.
$$

---

# 9. Relation to the Born Rule

## 9.1 What Has Been Explained

Taken together, this paper and the preceding paper explain the following.

1. A nontrivial zero sum of squares cannot be formed from real numbers alone. [Derived]
2. Minimal two-component closure requires a direction whose square is $-1$. [Derived]
3. This direction can be represented as a quarter turn. [Standard geometric interpretation]
4. Continuity, reversibility, and $Q_+$ preservation extend it to phase rotation. [Additional construction]
5. The phase is placed in the antisymmetric sector of an exchange-symmetric two-channel system. [Constructive hypothesis]
6. The finite-order condition discretizes the phase to an integer ratio. [Conditional derivation]
7. The squared weights of the interference amplitudes returned to the A/B basis are $\cos^2/\sin^2$. [Conditional derivation]

In this sense, the present paper offers

$$
\text{nontrivial quadratic closure}
$$

as a candidate root preceding the Born-type weights.

## 9.2 What Has Not Yet Been Explained

This paper does not prove

1. the complete Born rule of standard quantum mechanics,
2. why a positive-definite readout can be identified with an experimental frequency,
3. how single-trial selection into A or B occurs,
4. convergence of repeated-trial frequencies to $|r|^2$ and $|t|^2$,
5. that the complex field is the unique possible extension,
6. that standard complex conjugation is uniquely derived from nontrivial quadratic closure alone,
7. that the unitary group is uniquely derived without additional conditions,
8. projection measures in arbitrary dimension, or
9. that the preceding exchange kernel preserves the nontrivial quadratic-closure surface $Q_0=0$.

The conclusion is therefore not

> The Born rule has been derived completely.

More precisely, it is

> Nontrivial quadratic closure necessitates complex-phase structure, and when connected to an exchange-symmetric finite-order system, it yields discrete Born-type squared weights.

## 9.3 Significance for Circularity

A conventional circular explanation has the form

$$
\text{complex amplitude}
\rightarrow
\text{conjugate product}
\rightarrow
\cos^2.
$$

The construction here is

$$
\begin{aligned}
\text{nontrivial quadratic closure}
&\rightarrow\text{imaginary direction}
\rightarrow\text{phase rotation}\\
&\rightarrow\text{finite closure}
\rightarrow\text{channel amplitude}
\rightarrow\text{squared readout}.
\end{aligned}
$$

The remaining inputs in this construction are

- the positive-definite readout $Q_+$,
- adoption of that readout as the channel weight, and
- an additional physical process connecting that weight to experimental frequency.

This paper therefore does not derive the origin of the squared norm itself or the Born frequency law.

What it derives, or determines conditionally, is

- that nontrivial quadratic closure requires a nonreal direction,
- that the finite-order condition selects integer-ratio phases, and
- that projection of those phases onto an exchange-symmetric A/B basis gives discrete $\cos^2/\sin^2$ weights.

Thus, the inputs removed are the external provision of a nonreal direction and the external selection of the phase angle that closes. The inputs remaining are positive-definite squared readout and its frequency interpretation. This distinction precisely delimits the noncircular part of the derivation.

---

# 10. Meaning of the Closed-System Assumption

## 10.1 Definition of “Closure”

Closure in this paper does not immediately mean a physical container completely isolated in space from the outside world.

Its minimal meaning is that the relation

$$
\sum_kx_k^2=0
$$

is complete within the set of components under consideration.

Thus, closure means

> algebraic closure in which no undefined component outside the system is needed to evaluate the sum of the registered quantities under consideration.

## 10.2 Distinction from Physical Closure

Algebraic closure is different from thermodynamic, causal, or experimental closure.

This paper does not claim that a real quantum experiment is completely isolated from the outside world.

It considers instead the possibility that an enlarged effective system containing, over finite time and finite degrees of freedom,

- the target channels,
- the measurement apparatus, and
- the environment

may admit an internal closure relation.

## 10.3 Why This Is a Separate Paper

The finite-order exchange theorem in the preceding paper holds using standard unitary-operator theory alone.

Nontrivial quadratic closure in this paper, by contrast, contains foundational assumptions about the number system and the meaning of physical closure.

Combining the two in one paper would conflate the validity of the finite-order theorem with the physical validity of the closure axiom.

This paper therefore leaves the preceding work unchanged and presents the foundational argument independently as an additional interpretation.

---

# 11. Falsifiable Tests

## 11.1 Reconstruction in a Real Representation

Reconstruct the route from

$$
\sum_kx_k^2=0
$$

to the finite-order weights of the exchange system using only real two-dimensional rotations and no complex-number notation.

This would confirm that the result does not depend on symbolic manipulations specific to complex notation.

## 11.2 Comparison of Preserved Forms

Compare the following three conditions:

1. the conjugation-free quadratic closure $Q_0(x)=\sum x_k^2$,
2. the positive-definite norm $Q_+(x)=\sum|x_k|^2$, and
3. the indefinite-metric form $Q_G(x)=x^TGx$.

Determine which form generates the same finite-order sequence as the preceding exchange kernel under the fewest assumptions.

## 11.3 Uniqueness of the Conjugation Map

For a map that reverses a nonreal direction and gives a positive-definite readout, examine whether standard complex conjugation is determined uniquely when one imposes

- additivity,
- multiplicative order,
- involutivity, and
- fixation of real numbers.

## 11.4 Derivation of the Preservation Group

Classify the linear transformation group preserving the nontrivial closure set

$$
\mathcal C
=
\left\{x\neq0\mid\sum_kx_k^2=0\right\}.
$$

Then clarify under which conditions the subgroup that also preserves an added positive-definite readout becomes an orthogonal or unitary group.

## 11.5 Connection to the Frequency Law

For the finite-order exchange weights

$$
R_{n,m}
=
\cos^2\left(\frac{\pi m}{n}\right),
$$

perform repeated trials with observational backreaction and phase fluctuations, and determine the conditions under which

$$
f_A\to R_{n,m},
\qquad
f_B\to1-R_{n,m}.
$$

---

# 12. Discussion

## 12.1 Imaginary Numbers as Directions Required for Closure

The simplest consequence of this paper is that requiring a nontrivial solution of

$$
\sum_kx_k^2=0
$$

is impossible using real numbers alone.

For two components, one necessarily obtains

$$
y=\pm ix.
$$

An imaginary number can therefore be interpreted not as an arbitrary symbol appended to the real numbers, but as the minimal opposing direction required to realize quadratic closure nontrivially.

## 12.2 Conjugation Belongs to Readout, Not to Closure

The nontrivial closure equation requires no conjugation.

The square-root-of-minus-one direction is already required by

$$
x^2+y^2=0.
$$

Conjugation is needed later, when that nonreal direction is reversed to read the positive-definite quantity

$$
z^*z.
$$

Thus,

- the necessity of imaginary numbers,
- the role of conjugation, and
- the Born frequency law

are distinct problems.

## 12.3 Roles of the Two Kinds of Closure

Two kinds of closure appear in this paper, but they do not have the same role.

The first is the algebraic quadratic closure

$$
\sum_kx_k^2=0.
$$

The second is the finite-order phase closure

$$
\zeta^n=1.
$$

The first closure rigorously requires a nonreal direction. The second selects integer-ratio phases from the phase rotations supplied by the additional construction.

Returning that phase to the exchange-symmetric A/B basis gives the interference amplitudes

$$
\frac{1\pm\zeta}{2},
$$

whose positive-definite squared readouts are

$$
\cos^2\left(\frac{\pi m}{n}\right),
\qquad
\sin^2\left(\frac{\pi m}{n}\right).
$$

The condition that directly derives the discrete Born-type weights is therefore

$$
\boxed{
\text{phase closure}
+
\text{exchange projection}
+
\text{positive-definite squared readout}
\Rightarrow
\text{discrete Born-type weights}
}.
$$

Nontrivial quadratic closure provides a foundation for the nonreal direction of the complex-phase structure entering this derivation. It does not automatically follow, however, that the second exchange dynamics preserves the first closure surface; that dynamical connection remains problem C2.

## 12.4 Scope of Derivation and Additional Conditions

The central result derived directly here is that nontrivial quadratic closure requires a direction beyond the real number system, even if complex numbers are not postulated initially as an axiom of quantum theory.

We do not claim to have uniquely derived all the mathematics of standard quantum mechanics from nontrivial closure.

In particular, further conditions are still required between

- the complex field,
- standard conjugation,
- a positive-definite inner product,
- the unitary group, and
- the Born frequency law.

This paper does not conceal those conditions; it separates the derived part from the interpretive connection.

---

# 13. Conclusion

This paper began from the nontrivial quadratic-closure condition

$$
\boxed{
\sum_{k=1}^{N}x_k^2=0,
\qquad
(x_1,\ldots,x_N)\neq(0,\ldots,0),
\qquad
x_k\in K,
\quad
\mathbb R\subseteq K
}.
$$

If every component is real, each squared term is nonnegative, so the equation has only the trivial solution. The solution vector of nontrivial closure therefore lies outside $\mathbb R^N$. In fixed closure coordinates this is equivalent to at least one component being nonreal.

In the minimal two-component system,

$$
x^2+y^2=0,
\qquad x\neq0,
$$

we obtain

$$
\boxed{
y=\pm ix
}.
$$

No complex conjugation is used in this result. A square-root-of-minus-one direction is required to form a nontrivial zero sum of squares.

If continuity, reversibility, and preservation of a positive-definite readout are added to this nonreal direction to extend it to a phase rotation, and the resulting phase is placed in the antisymmetric mode of an exchange-symmetric two-channel system, one obtains

$$
U=P_s+\zeta P_a.
$$

The finite-order condition

$$
\zeta^n=1
$$

selects the discrete phase

$$
\zeta=e^{-2\pi im/n}.
$$

The interference amplitudes in the A/B basis are

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2},
$$

and their positive-definite squared readouts are

$$
\boxed{
|r|^2
=
\cos^2\left(\frac{\pi m}{n}\right)
}
$$

and

$$
\boxed{
|t|^2
=
\sin^2\left(\frac{\pi m}{n}\right)
}.
$$

Separating derivation, additional construction, and the unresolved connection, the relation established here is

$$
\boxed{
\begin{gathered}
\text{nontrivial quadratic closure}
\Rightarrow\text{nonreal direction}
\quad[\text{derived}]\\
\text{nonreal direction}
+\{\text{continuity, reversibility, }Q_+\text{ preservation}\}\\
\longrightarrow\text{phase rotation}
\quad[\text{additional construction}]\\
\text{phase rotation}+\text{exchange symmetry}+\text{finite order}\\
+\text{squared readout}
\Rightarrow\text{discrete Born-type weights}
\quad[\text{conditional derivation}]\\
Q_0=0\text{ and its dynamical connection to the exchange kernel}\\
[\text{connection problem C2}]
\end{gathered}
}.
$$

This is not a derivation of the complete Born rule. What has been derived rigorously is that nontrivial quadratic closure cannot be achieved using real numbers alone, that the minimal two-component system requires a square-root-of-minus-one direction, and that combination with the preceding finite-order exchange theorem recovers discrete Born-type weights.

Future problems include unique derivations of conjugation, positive-definite readout, unitarity, and a frequency law, as well as construction of an exchange map that preserves $Q_0=0$ or a projection rule onto the closure surface. The central significance of this paper is that complex phase is positioned not as a tool specific to quantum theory assumed from the beginning, but as a necessary condition for nontrivial quadratic closure.

---

# References

## Self-Citations

1. Noriaki Kihara, “Discovery of Finite-Order Resonances in Iterated Exchange Scattering: Identification of the Origin of Peaks Near the Fine-Structure Values 137 and 128 and a Reproducible Wave-Packet Mathematical Model,” Version DOI: 10.5281/zenodo.21421367, Concept DOI: 10.5281/zenodo.21421366, 2026.
2. Noriaki Kihara, “First-Principles Derivation of the Minimal Nontrivial 137 Cells in a Closed Two-Wave Exchange System: Two Principal Candidates, a Conditional Dual Attractor, and a Correspondence Hypothesis with the Fine-Structure Constant,” v3.0, 2026. Repository version: [`20260717_complete_paper_v3.md`](../20260717/20260717_complete_paper_v3.md). [Japanese].
3. Noriaki Kihara, “Emergence of Discrete Born-Type Weights in Iterated Two-Channel Exchange Systems: A Finite-Order Recurrence Law from Wave-Packet Localization Transfer, Metastable Two-State Dynamics, and Observation Selection,” Version DOI: 10.5281/zenodo.21422471, Concept DOI: 10.5281/zenodo.21422470, 2026.
4. Noriaki Kihara, “Experimental Specification v1 for Low-Localization and Harmonic Transfer Readout in Fermion-Like Collisions Using an Exchange-Interference Scattering Matrix,” 2026. [Japanese].
5. Noriaki Kihara, “Preliminary Experimental Summary v1 of the Acceleration Basis and Localization Exchange in Fermion-Like Collisions Using an Exchange-Interference Scattering Matrix,” Zenodo Concept DOI: 10.5281/zenodo.21333766, 2026. [Japanese].
6. Noriaki Kihara, “Preliminary Experimental Summary v1 of C Weak Readout and D Strong-Observation Selection at a White-Cat, Black-Cat, and Gray-Cat Metastable Interface,” Zenodo Concept DOI: 10.5281/zenodo.21353208, 2026. [Japanese].

## External References

7. Max Born, “Zur Quantenmechanik der Stoßvorgänge,” *Zeitschrift für Physik* **37**, 863–867 (1926). DOI: 10.1007/BF01397477.
8. Andrew M. Gleason, “Measures on the Closed Subspaces of a Hilbert Space,” *Journal of Mathematics and Mechanics* **6**, 885–893 (1957). DOI: 10.1512/iumj.1957.6.56050.

---

# Appendix A. Classification of Claims

| Claim | Classification |
|---|---|
| Requiring a nontrivial solution of $\sum_kx_k^2=0$ | First axiom of this paper |
| No nontrivial solution of $\sum_kx_k^2=0$ exists using only real components | Derived proposition |
| A nontrivial solution vector belongs to $K^N\setminus\mathbb R^N$ | Derived proposition |
| In fixed closure coordinates, at least one component is nonreal | Equivalent component representation of the preceding proposition |
| In the two-component system $x^2+y^2=0$, $y/x$ is a square root of $-1$ | Derived proposition |
| The minimal subfield generated by two-component nontrivial closure is $\mathbb R(i)\cong\mathbb C$ | Derived consequence |
| If a square root of $-1$ is denoted by $i$, then $y=\pm ix$ | Derived consequence |
| Reading $\pm i$ as a quarter turn | Standard geometric interpretation |
| Extending the nonreal direction to the continuous phase group $U(1)$ | Construction adding continuity, reversibility, and preservation |
| Standard complex conjugation is uniquely derived from nontrivial quadratic closure | Not derived |
| The unitary group is uniquely derived from nontrivial quadratic closure | Not derived |
| The exchange-symmetric finite-order kernel gives discrete Born-type weights | Derived in the preceding paper |
| The preceding exchange kernel preserves the nontrivial quadratic-closure surface $Q_0=0$ | False in general; connection problem C2 |
| Nontrivial quadratic closure is the physical origin of complex-phase structure | Foundational interpretive hypothesis |
| Discrete Born-type weights give experimental frequencies | Untested |
| The complete Born rule has been derived | Not derived |

---

# Appendix B. Minimal Derivation Sequence

## B.1 Nontrivial Closure

$$
\sum_{k=1}^{N}x_k^2=0,
\qquad
(x_1,\ldots,x_N)\neq(0,\ldots,0),
\qquad
x_k\in K,
\qquad
\mathbb R\subseteq K.
$$

## B.2 Impossibility Using Real Numbers Alone

$$
x_k\in\mathbb R
\quad\Longrightarrow\quad
x_k^2\ge0.
$$

Therefore,

$$
\sum_kx_k^2=0
\quad\Longrightarrow\quad
x_k=0\ \forall k.
$$

For a nontrivial solution, first in vector notation,

$$
\mathbf{x}\in K^N\setminus\mathbb R^N.
$$

The equivalent component representation in fixed closure coordinates is

$$
\exists j:\ x_j\notin\mathbb R.
$$

## B.3 Minimal Two-Component System

From

$$
x^2+y^2=0,
\qquad x\neq0,
$$

we have

$$
\left(\frac{y}{x}\right)^2=-1.
$$

If a square root of $-1$ is denoted by $i$,

$$
y=\pm ix.
$$

## B.4 Additional Construction and Phase Closure

The following is not a direct consequence of nontrivial quadratic closure alone. We add continuity, reversibility, preservation of a positive-definite readout, and exchange symmetry, representing the nonreal direction as the antisymmetric-sector phase $\zeta$. Imposing finite order on this additional construction gives

$$
\zeta^n=1
\quad\Longrightarrow\quad
\zeta=e^{-2\pi im/n}.
$$

## B.5 Exchange Projection

$$
U=P_s+\zeta P_a
=
\frac12
\begin{pmatrix}
1+\zeta&1-\zeta\\
1-\zeta&1+\zeta
\end{pmatrix}.
$$

## B.6 Discrete Born-Type Weights

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2},
$$

and therefore

$$
\boxed{
|r|^2
=
\cos^2\left(\frac{\pi m}{n}\right)
}
$$

$$
\boxed{
|t|^2
=
\sin^2\left(\frac{\pi m}{n}\right)
}.
$$
