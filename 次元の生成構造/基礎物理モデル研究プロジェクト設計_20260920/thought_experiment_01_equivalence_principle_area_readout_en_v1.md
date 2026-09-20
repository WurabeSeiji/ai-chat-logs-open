# The First Thought Experiment: Tracing the Equivalence Principle Down to Two Nameless Degrees of Freedom
## — Deriving the Sign of the Sum of Squares, the Curvature Radius, and Three Systems (Inertial, Uniformly Accelerated, Rotating) from the Conservation of an Area Readout —

**Author:** Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**Version DOI:** 10.5281/zenodo.22857953  
**Concept DOI:** 10.5281/zenodo.22857952  
**Version:** v1.0  
**Date:** 2026-09-20

---

## Abstract

Einstein's equivalence principle is traced down to a minimal system that presupposes no background spacetime, no metric, and no complex numbers. The inputs are two principles, one axiom candidate, and three design assumptions (plus, in Sec. 6 only, a fourth design assumption: a fixed resolution $\varepsilon$). The principles are two nameless values $(a,b)\sim(b,a)$, as the minimal expression that distinguishes state from relation, and the readout principle that a state alone is unreadable and only products are readable. The axiom candidate is that a conserved readout exists. The design assumptions are real values, linearity, and bilinearity.

The following is then derived.

1. From the principle that a state alone is unreadable, the readout of two configurations is uniquely determined, up to a constant multiple, to be the area form $\omega(X,Y)=a_Xb_Y-b_Xa_Y$. Both the Euclidean and the Minkowski inner products assign a value to a single configuration, and so violate the principle.
2. Conservation of the signed readout is equivalent to $\det S=1$ (since only the magnitude $\lvert q\rvert$ is readable, the orientation-reversing branch $\det S=-1$ also remains; this paper treats the branch $\det S=1$). The law then becomes an undirected second-order rule $X_{k+1}+X_{k-1}=\tau X_k$, written with a single real number $\tau=\operatorname{tr}S$. In this form the velocity (first difference) changes sign with the choice of direction, and the first readable dynamical quantity is the centripetal acceleration proportional to the configuration, $\Delta^2X_k=-\kappa X_k$ ($\kappa=2-\tau$). The operational content of the equivalence principle — "only the acceleration is readable, and its origin is not" — appears not as a postulate but as the structure of readout.
3. The sum of squares is not an input; it is induced as the area spanned by a configuration and its interaction image, $Q_S(X)=\omega(X,SX)$. Its signature is set by $\tau$: $a^2+b^2$ for $|\tau|<2$, $a^2-b^2=a^2+(ib)^2$ for $|\tau|>2$, degenerate at $|\tau|=2$. $\sum z_n^2=C^2$ holds in all three classes, and whether $z_2$ is real, imaginary, or nilpotent is decided by the class of the interaction.
4. The right-hand side $C$ of $\sum z_n^2=C^2$ equals exactly, at any step size, the curvature radius of the orbit measured in the metric induced by the interaction itself: $|\Delta^2X|/|\Delta X|^2=1/|C|$.
5. The three classes are matched with the inertial, the uniformly accelerated (Born's hyperbolic motion, $C=c^2/g$, the 1907 clock effect), and the rotating systems. The two classical scenes of the equivalence principle are the two non-degenerate classes of one structure.
6. Zero closure has two substantive forms: the null cone (configurations whose direction is unchanged by the interaction), and the orbit sum $\sum_k z_k^2=0$ of a closed orbit ($n\ge3$; the orbit is isotropic with respect to its own induced metric).
7. The closure number $n$ need not be a postulate. Under a finite resolution $\varepsilon$, $n$ takes only continued-fraction denominators of the rotation number of the law, and satisfies $n<2\pi C/\varepsilon$. $n$ is not a degree of freedom independent of $C$; it is the curvature radius measured by the resolution.

Numerical verification used 15 items whose failure conditions were declared in advance, plus 2 items added after the failure of E4b (E4b' is a post-hoc diagnostic; E6 is a prediction built from the cause analysis and declared before execution). One preregistered item (E4b) actually failed; the cause and the treatment are given in Sec. 7. Finally, Sec. 8 returns to the minimal skeleton of the design document (the zeroth paper) what this thought experiment can say about it (for example, that no closed structure arises from a linear, real, first-order skeleton).

---

## 0. Position of This Paper

This is the first thought experiment under the design document [0]. Its subject is to start from the equivalence principle and read out a nameless state–relation system together with the sum of squares, the curvature radius, and zero closure. The construction makes every input explicit and shows how far each input carries.

To the six categories of Sec. 14 of the design document we add "Derived," and every statement in the text carries one of the following tags.

| Tag | Meaning |
|---|---|
| [Principle] | philosophical principle |
| [Axiom candidate] | a candidate to be placed as an axiom |
| [Goal image] | correspondence with the goal images of Sec. 7 of the design document |
| [Design assumption] | an assumption of preliminary design (may be removed; scheduled to be removed) |
| [Derived] | proved from the statements tagged above (proofs in Appendix A) |
| [Numerical] | confirmed by numerical experiment |
| [Comparison] | correspondence with known physics; not a derivation |

**List of inputs (nothing else is used)**

| Symbol | Tag | Content |
|---|---|---|
| P1 | [Principle] | To distinguish state from relation, two values are needed. Two nameless values as the minimal such expression. $(a,b)\sim(b,a)$. Which is the state and which is the relation cannot be decided |
| P2 | [Principle] | A state alone is unreadable. Only product-type operations are readable, and readout is interaction (design document Sec. 5) |
| A1 | [Axiom candidate] | A conserved readout exists (design document Sec. 7.1). This paper uses it as conservation of the signed $q$ ($\det S=1$). If only conservation of the readable quantity $\lvert q\rvert$ is required, the branch $\det S=-1$ remains (Sec. 3.1). Its relation to the condition "the law does not distinguish a direction" is T3' |
| D1 | [Design assumption] | Each of the two values is real. That they must be real is not derived (Sec. 10). A configuration is written $X=(a,b)\in\mathbb R^2$ |
| D2 | [Design assumption] | The interaction is an invertible linear map $X\mapsto SX$ |
| D3 | [Design assumption] | "Product-type operation" is understood as a bilinear form $B(X,Y)$ of two configurations |
| D4 | [Design assumption] | (Sec. 6 only) the resolution $\varepsilon$ is placed as a fixed value |

Not used: a metric, an inner product, complex numbers, the closure condition $S^n=I$, a time coordinate, the balance equation $GM/R^2=\omega^2R$.

On D2. Since absolute magnitudes are unreadable and only ratios are readable, the law should commute with uniform dilation of configurations. Linearity is adopted as the minimal such candidate. Homogeneous but nonlinear candidates are deferred to the next paper.

---

## 1. The Operational Content of the Equivalence Principle

The equivalence principle has two classical scenes. One is the 1907 equivalence between a uniformly accelerating box and a uniform gravitational field [1]. The other is the equivalence between the centrifugal acceleration of a rotating system and gravity [2,5].

Combining the operational content common to both scenes with the relativity principle (velocity is unreadable) into a single statement:

$$
\boxed{\text{From inside, velocity is unreadable. The first readable quantity is the acceleration, and its origin (inertia or gravity) is unreadable.}}
$$

The first half is the relativity principle; the second half is the content of the equivalence principle. What this paper derives is the readout structure corresponding to this statement (Sec. 3.3), not the equivalence principle itself. The question of this paper is: what is the minimal system in which this statement is meaningful, and what is forced there?

---

## 2. The Minimal System and Its Readout

### 2.1 Two nameless degrees of freedom [Principle P1]

To express a state and a relation as distinct, one value is not enough; at least two values are needed. Being a minimal system, we take two and write $(a,b)$. The two can be distinguished, but there is no internal ground for deciding which is the state and which is the relation. What is being counted is the number of values: a complex number $a+ib$ is not one value but two. Conversely, this paper does not deny reading the two values $(a,b)$ as $a+ib$; whether that reading holds is decided by the class of the interaction (Sec. 4.1). The swap is written

$$
P=\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$

### 2.2 The readout is uniquely determined [Derived T1]

By P2, a readout that assigns a value to a single configuration is not allowed. $B(X,X)$ is a quantity determined by the single configuration $X$ alone, so if it is nonzero, a value has been assigned to a single configuration. Under D3, P2 can therefore be written as

$$
B(X,X)=0\quad(\forall X)
$$

Then $B$ is antisymmetric, and in two degrees of freedom it is unique up to a constant multiple:

$$
\boxed{\;\omega(X,Y)=a_Xb_Y-b_Xa_Y\;}
$$

Four remarks.

(i) The symmetric candidates — $a_Xa_Y+b_Xb_Y$ (Euclidean), $a_Xa_Y-b_Xb_Y$ (Minkowski), $a_Xb_Y+b_Xa_Y$ — all assign $B(X,X)\neq0$ to a single configuration. Besides violating P2, choosing among them is a choice of geometry. Placing $a^2+b^2$ as the first readout amounts to that choice.

(ii) $\omega(PX,PY)=-\omega(X,Y)$. Renaming flips the sign, so under P1 the sign of $\omega$ is unreadable. What is readable is $|\omega|$ and ratios of two $\omega$'s.

(iii) Viewing each degree of freedom as a real sequence $a_k,b_k$, the quantity $\omega(X_k,X_{k+1})=a_kb_{k+1}-a_{k+1}b_k$ is the Wronskian of the two sequences. It is a relational quantity between the two nameless degrees of freedom, and the relational quantity with oneself is zero.

(iv) With a single value, a bilinear form satisfying $B(X,X)=0$ is identically zero, and there is nothing to read. Only with two values does a readout exist — and it is then unique. "Two" is minimal from the side of P1 and from the side of P2 alike.

---

## 3. Conservation and the Law

### 3.1 Conservation is equivalent to unimodularity [Axiom candidate A1] [Derived T2]

By P2, readout is interaction. Define the readout of a configuration and its interaction image:

$$
q_k:=\omega(X_k,X_{k+1}),\qquad X_{k+1}=SX_k
$$

Since $\omega(SX,SY)=\det S\cdot\omega(X,Y)$, we have $q_{k+1}=\det S\cdot q_k$. If $S$ is not a scalar multiple of the identity, there exist initial configurations with $q_0\neq0$, so

$$
\boxed{\;\text{for }S\neq\lambda I\text{,}\quad\text{the signed }q_k\text{ is conserved for all initial configurations}\iff\det S=1\;}
$$

For $S=\lambda I$, $q\equiv0$ and there is no readout at all (remark (iv) of Sec. 4.1). What is conserved is the area spanned by a configuration and the next. In the continuum limit this corresponds to constant areal velocity (Kepler's second law) — but no metric has been used here.

One distinction is needed at this point. As remark (ii) of Sec. 2.2 says, only $\lvert q\rvert$ is readable. If only conservation of $\lvert q\rvert$ is required, the condition is $\lvert\det S\rvert=1$, and the orientation-reversing branch $\det S=-1$ remains. On that branch the sign of $q$ flips at every step. $S^2$ has $\det S^2=1$ and $\operatorname{tr}S^2=\tau^2+2\ge2$, so viewed every two steps it falls into the hyperbolic class of Sec. 4 (the identity map if $\tau=0$). The induced form $\omega(X,SX)$ is always indefinite, and its sign alternates at every step. The minimal example is the swap $P$ itself, with $\omega(X,PX)=a^2-b^2$. From here on this paper treats the branch $\det S=1$ on which the signed $q$ is conserved. Whether A1 is to be read as conservation of $q$ or of $\lvert q\rvert$ is left undecided here (Sec. 10).

### 3.2 The law is written with a single real number [Derived T3]

For $\det S=1$, the Cayley–Hamilton theorem gives $S+S^{-1}=\tau I$ ($\tau=\operatorname{tr}S$). Therefore

$$
\boxed{\;X_{k+1}+X_{k-1}=\tau X_k\;}
$$

This form has three properties. $\tau I$ commutes with $P$, so the form is unchanged by renaming. It is unchanged under $k\to-k$. The parameter is a single real number. Conversely, given two independent configurations $X_0,X_1$ and $\tau$, the map $S$ is uniquely determined. That is, **the law is the single number $\tau$, and $S$ is the law together with initial data**.

The first-order form $X_{k+1}=SX_k$ cannot, in general, be defined on nameless configurations. Maps with $PSP=S$ are restricted to the form $\gamma I+\beta P$, whose eigenvalues are real. Rotations satisfy $PSP=S^{-1}$: renaming becomes a reversal of the direction of iteration, exchanging "next" and "previous." Writing a rotation in first-order form is incompatible with namelessness; the second-order undirected law has no such problem.

What about the reverse direction? [Derived T3'] For a general invertible $S$ without assuming $\det S=1$, writing $d:=\det S$, Cayley–Hamilton gives $X_{k+1}+d\,X_{k-1}=\tau X_k$. The sequence read backwards, $Y_k=X_{-k}$, obeys $Y_{k+1}+d^{-1}Y_{k-1}=(\tau/d)\,Y_k$. The two laws coincide iff $d^2=1$ and $\tau(d-1)=0$, so the law keeps its form under $k\to-k$ exactly in two cases: $d=1$, and $d=-1$ with $\tau=0$ (a reflection with $S^2=I$; the swap $P$ itself is an example). Meanwhile the readout obeys $q_{k+1}=d\,q_k$, so on orbits with $q_0\neq0$ the signed $q$ is conserved iff $d=1$, and the readable quantity $\lvert q\rvert$ is conserved iff $\lvert d\rvert=1$. In summary,

$$
\boxed{\;q\ \text{conserved}\ (d=1)\ \Longrightarrow\ \text{the law does not distinguish a direction}\ \Longrightarrow\ \lvert q\rvert\ \text{conserved}\ (\lvert d\rvert=1)\;}
$$

Neither converse holds. $S=P$ does not distinguish a direction, yet the sign of $q$ flips at every step. $S=\begin{pmatrix}1&1\\1&0\end{pmatrix}$ ($d=-1$, $\tau=1$) conserves $\lvert q\rvert$, yet the reversed law changes as $\tau\to-\tau$. If $\lvert d\rvert\neq1$, then on orbits with $q_0\neq0$ the quantity $\lvert q_k\rvert=\lvert d\rvert^k\lvert q_0\rvert$ changes monotonically and an arrow given by the law becomes readable. Hence "the law does not distinguish a direction" is weaker than conservation of the signed quantity and stronger than conservation of the readable quantity; it is not a paraphrase of A1. Note that this is a statement at the level of the law: even when the law is undirected, an arrow may be readable on individual orbits (remark (iii) of Sec. 4.1).

### 3.3 What is readable is the acceleration [Derived T4]

On orbits with $q\neq0$, the law itself is readable as a ratio of areas:

$$
\tau=\frac{\omega(X_{k-1},X_{k+1})}{\omega(X_{k-1},X_k)}
$$

Numerator and denominator flip sign together under renaming, so the ratio is readable; it is also scale-free.

Consider differences. The first difference $\Delta X_k=X_{k+1}-X_k$ changes sign with the choice of direction of $k$. The second difference does not:

$$
\boxed{\;\Delta^2X_k:=X_{k+1}-2X_k+X_{k-1}=-\kappa X_k,\qquad\kappa:=2-\tau\;}
$$

As readable statements: $\omega(X_k,\Delta^2X_k)=0$ (the acceleration points to the center), and $\kappa=-\omega(X_{k-1},\Delta^2X_k)/\omega(X_{k-1},X_k)$.

[Comparison] This is the operational content of Sec. 1. The first readable dynamical quantity is the acceleration. There is only one of it: calling it the bending of the orbit (the inertial side) or the central field $-\kappa X$ (the gravitational side) is the same readout. Since $\kappa$ does not depend on the orbit (the same for every $C$ and every phase), the universality of free fall also holds automatically by linearity. In this system the equivalence principle is not a balance of two quantities but the fact that there is only one readout.

---

## 4. The Three Classes

### 4.1 The induced sum of squares and its sign [Derived T5]

Write the conserved readout $q$ as a function of the configuration. With $\Omega=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ and $N:=(S-S^{-1})/2$,

$$
Q_S(X):=\omega(X,SX)=X^{\mathsf T}GX,\qquad G=\Omega N
$$

$$
\boxed{\;\det G=1-\frac{\tau^2}{4},\qquad N^2=-\Bigl(1-\frac{\tau^2}{4}\Bigr)I\;}
$$

The sum of squares is not an input; the interaction induces it. The signature is decided by $\tau$ alone, and the normal forms combine into a single formula:

$$
\boxed{\;a^2+s\,b^2=C^2,\qquad s=\operatorname{sgn}(4-\tau^2)\in\{+1,\,0,\,-1\}\;}
$$

Setting $z_1=a$ and $z_2$ equal to $b$, $\epsilon b$ ($\epsilon^2=0$), or $ib$ according to $s=+1,\,0,\,-1$, this is $\sum z_n^2=C^2$.

| Class | Condition | $N^2$ | Emerging unit | Normal form | Form of $\sum z_n^2=C^2$ | Orbit |
|---|---|---|---|---|---|---|
| Elliptic | $\lvert\tau\rvert<2$, $\tau=2\cos\theta$ | $-\sin^2\theta\,I$ | $i$ ($J=N/\sin\theta$) | $a^2+b^2$ | $z_2=b$ | bounded; closes or is quasi-periodic |
| Parabolic | $\lvert\tau\rvert=2$ | $0$ | nilpotent $\epsilon$ ($\epsilon^2=0$) | $a^2$ (degenerate) | $z_2=\epsilon b$ | $a$ constant, $b$ grows linearly in $k$ |
| Hyperbolic | $\lvert\tau\rvert>2$, $\lvert\tau\rvert=2\cosh H$ | $+\sinh^2H\,I$ | $j$ ($j^2=+1$) | $a^2-b^2$ or $2ab$ | $z_2=ib$ | grows exponentially |

Four remarks.

(i) The imaginary structure appears not as a value of the state but on the side of the interaction, as $N^2\propto-I$. It is, however, one of three cases: the class of the interaction decides which of $i$, $\epsilon$, $j$ appears. $\sum z_n^2=C^2$ holds in all three classes.

(ii) $a^2+b^2$ is a normal form, not an input. In the elliptic class with $PSP=S^{-1}$, the general description gives $Q=A(a^2+b^2)+2Bab$ ($A>\lvert B\rvert$). The cross term is removed by a change of description compatible with namelessness (commuting with $P$): $a'=\alpha a+\beta b,\ b'=\beta a+\alpha b$. In the general description the naive $a^2+b^2$ is not conserved (E1d of Sec. 7; Fig. 1(d)).

(iii) The two hyperbolic normal forms are two descriptions of the same thing. If $PSP=S$ (strictly nameless), $S$ is a boost $\begin{pmatrix}\cosh H&\sinh H\\\sinh H&\cosh H\end{pmatrix}$ with invariant $a^2-b^2$. If $PSP=S^{-1}$, then $a\mapsto e^Ha,\ b\mapsto e^{-H}b$ with invariant $ab$ (light-cone coordinates). The swap flips the sign of $a^2-b^2$, so which is timelike and which is spacelike cannot be decided at this stage — another appearance of P1's "which is the state and which is the relation cannot be decided." Only in the strictly nameless class can the direction of iteration be read internally, as the direction in which the swap-invariant combination $u=(a+b)/\sqrt2$ grows. The elliptic class has no direction.

(iv) On the branch $\det S=1$, closed orbits with $n\le2$ are unreadable. If $S^n=I$ with $n\le2$ then $S=\pm I$, so $N=0$ and $Q_S\equiv0$: a configuration and its image are parallel and span no area. Readable closed orbits on this branch are restricted to $n\ge3$. The reflections of the branch $\det S=-1$ ($S^2=I$, e.g. $P$) have period 2 with $\lvert q\rvert\neq0$, but are not treated in this paper (Sec. 3.1). ($\tau\le-2$ is the mirror image under $S\to-S$; apart from alternating signs it is the same.)

From here on we normalize $\hat G:=G/\sqrt{\lvert\det G\rvert}$, $C^2:=\lvert X^{\mathsf T}\hat GX\rvert$, and write $\lvert v\rvert_{\hat G}:=\sqrt{\lvert v^{\mathsf T}\hat Gv\rvert}$ (in the hyperbolic class $\hat G$ is indefinite, so this is a notation, not a norm). This normalization matches the area element of the metric with $\lvert\omega\rvert$, and is restricted to the non-degenerate classes ($\lvert\tau\rvert\neq2$). In the parabolic class $\det G=0$ and $C$ is determined only up to a constant multiple.

### 4.2 $C$ is the curvature radius [Derived T6]

From $S^{\mathsf T}GS=G$ one obtains $(S-I)^{\mathsf T}G(S-I)=\kappa G$. In the language of areas, $\omega(\Delta X_k,\Delta X_{k+1})=\kappa\,\omega(X_k,X_{k+1})$. Combined with T4, on every orbit of a non-degenerate class ($\lvert\tau\rvert\neq2$) with $C\neq0$,

$$
\lvert\Delta X\rvert_{\hat G}^2=\lvert\kappa\rvert C^2,\qquad\lvert\Delta^2X\rvert_{\hat G}=\lvert\kappa\rvert\,\lvert C\rvert
$$

$$
\boxed{\;\frac{\lvert\Delta^2X\rvert_{\hat G}}{\lvert\Delta X\rvert_{\hat G}^{2}}=\frac{1}{\lvert C\rvert}\;}
$$

The left-hand side — the acceleration divided by the square of the step length — is what this paper defines as the discrete curvature. Therefore the right-hand side $C$ of $\sum z_n^2=C^2$ equals exactly, independently of the step size, the curvature radius of the orbit measured in the metric induced by the interaction itself. The same formula holds in the Euclidean and in the Lorentzian type.

In this form, "$C$ is the curvature radius" is a theorem, not a naming. No radius needs to be imported from outside; both the metric and the curvature are read from the interaction. The absolute value of $C$ of a single orbit, however, is unreadable; what is readable is the ratio $C_1/C_2$ of two orbits, and $C/\varepsilon$ of Sec. 6.

### 4.3 The three classes and three systems [Comparison C1–C3]

**C1: parabolic class and inertial system.** $\kappa=0$ and $X_k=X_0+kV$; the acceleration is zero. In the normal form $a$ is constant and $b$ accumulates linearly in $k$. [Goal image] The behavior $L_t(n)\sim n$ of Sec. 7.6 of the design document (a direction increasing monotonically as accumulated causal processing) is the behavior of this class.

**C2: hyperbolic class and the uniformly accelerated system.** In the strictly nameless class $S$ is a boost and the orbit is the hyperbola $a^2-b^2=C^2$. Reading $(a,b)\leftrightarrow(x,ct)$, this is Born's hyperbolic motion [3], the worldline of proper acceleration $g$, with

$$
\boxed{\;C=\frac{c^2}{g}\;}
$$

T6's "curvature $=1/C$" becomes "proper acceleration $=c^2/C$." The family of orbits of the same $S$ with different $C$ is the uniformly accelerated rigid frame (the Rindler frame [4]). By T6 the proper length per step is $\sqrt{\lvert\kappa\rvert}\,C$, so the ratio of clock rates at heights $C_1$ and $C_2=C_1+h$ is

$$
\frac{C_2}{C_1}=1+\frac{gh}{c^2}
$$

This agrees with the gravitational clock effect Einstein derived from the equivalence principle in 1907 [1] (E4c of Sec. 7 is the numerical check of this ratio).

**C3: elliptic class and the rotating system.** The family of orbits of the same $S$ with different $C$ is a rigid rotation; the step length is proportional to $C$ and the centripetal acceleration per step is $\kappa C$. The centrifugal scene of the equivalence principle corresponds to this.

Thus the two classical scenes of the equivalence principle are the two non-degenerate classes of one structure. In both, "acceleration $\div$ speed$^2=1/C$" holds.

**Scope of the comparison.** What agrees: the shape of the orbits, the signature of the induced metric, the curvature–acceleration relation, and the clock-rate ratio — four items. What is not claimed to agree: 3+1 dimensions, the light cone as physical causal structure, the inverse-square law, and mass. In particular, the balance equation $GM/R^2=\omega^2R$ cannot be imposed here. In a linear law $\kappa$ does not depend on $C$, and the magnitude of the central field is linear in $C$: $\lvert\Delta^2X\rvert_{\hat G}=\lvert\kappa\rvert\,C$. In Newtonian gravity this corresponds not to the field of a point mass but to the field inside a uniform background ($g\propto r$). The inverse-square law requires a relation between $\kappa$ and $C$ — that is, nonlinearity or a family of systems — which is not available at this stage.

---

## 5. The Content of Zero Closure

Transposing $a^2+b^2=C^2$ into $a^2+b^2+(iC)^2=0$ is a rewriting of T5 and has no new content. The triple $(a,b,iC)$ consists of real, real, and purely imaginary constants that cannot be interchanged, so it is not even nameless. The forms in which zero has content are the following two.

### 5.1 The null cone consists of eigenconfigurations [Derived T7]

$$
Q_S(X)=0\iff SX\parallel X
$$

A configuration with zero readout area is one whose direction is unchanged by the interaction. None exist in the elliptic class; every nonzero configuration is readable. In the hyperbolic class there are two directions; in the boost normal form, $a=\pm b$, i.e. $a^2+(ib)^2=0$.

### 5.2 In the hyperbolic class the null cone is a projective attractor [Derived T8]

Scales are unreadable, so only the direction of the configuration is readable. Let $v_\pm$ be the two eigendirections of the hyperbolic class (eigenvalue magnitudes $e^{\pm H}$) and decompose $X_k=A_kv_++B_kv_-$. The ratio of the two components is readable as a ratio of areas:

$$
\rho_k:=\left\lvert\frac{B_k}{A_k}\right\rvert=\left\lvert\frac{\omega(v_+,X_k)}{\omega(X_k,v_-)}\right\rvert=\rho_0\,e^{-2Hk}
$$

Being a ratio of areas, this is scale-free and uses no externally supplied length (the normalization of $v_\pm$ only affects $\rho_0$ by a constant). The direction of the configuration converges exponentially to the null direction $v_+$. If $\rho$ below $\varepsilon_{\rm rel}$ cannot be distinguished from zero, then after roughly

$$
k^{*}\approx\frac{1}{2H}\ln\frac{\rho_0}{\varepsilon_{\rm rel}}
$$

steps the configuration becomes indistinguishable from a zero-closure configuration ($\rho=0$). The conserved quantity remains conserved but becomes unreadable. (This statement was obtained from the diagnosis of the failure of E4b. E6 of Sec. 7 is the numerical verification that the quantity normalized by the Euclidean length of the description, $\lvert q\rvert/(\lvert X_k\rvert\lvert X_{k+1}\rvert)$, decays at the same rate $-2H$.)

### 5.3 Orbit sums of closed orbits, and reading the metric [Derived T9]

For a closed orbit of order $n\ge3$,

$$
\boxed{\;\sum_{k=0}^{n-1}X_kX_k^{\mathsf T}=\frac{nC^2}{2}\,\hat G^{-1}\;}
$$

In the normal description, writing $z_k=a_k+ib_k$,

$$
\sum_kz_k=0,\qquad\boxed{\sum_kz_k^2=0},\qquad\sum_kz_k^{p}=0\quad(n\nmid p)
$$

The $n$ values $z_k$ are carried into one another cyclically by the interaction, so they are mutually nameless and of equal amplitude. The minimum is $n=3$; for $n=2$ it does not hold.

There are two meanings. First, the zero closure of the sum of squares is the normal-description expression of the statement that "the orbit is isotropic with respect to its own induced metric." In components, $\sum a_k^2=\sum b_k^2$ and $\sum a_kb_k=0$: this is the legitimate form of "no cross term." Second, the formula gives a procedure for reading the metric from the orbit alone.

The comparison with the nameless complex zero closure of the existing series is not carried out in this paper, in accordance with the appendix of the design document, which forbids using previously published results as grounds. T9 can serve as the minimal point of contact when that comparison is made.

---

## 6. The Closure Number Need Not Be a Postulate

There is no need to postulate the closure condition $S^n=I$. For a generic elliptic law, $\alpha:=\theta/2\pi$ is irrational and the orbit never closes exactly.

[Design assumption D4] Place a fixed resolution $\varepsilon$ measured in the induced metric. Sec. 8.3 of the design document asks that the resolution be determined by the interacting structures themselves; the fixed value is a stand-in for that, declared as an assumption scheduled to be removed.

[Derived T10] The first return to a configuration indistinguishable from the start is the smallest $n$ with

$$
2C\,\Bigl\lvert\sin\frac{n\theta}{2}\Bigr\rvert<\varepsilon\iff\|n\alpha\|<\delta:=\frac1\pi\arcsin\frac{\varepsilon}{2C}
$$

where $\|x\|$ is the distance from $x$ to the nearest integer. We call this the effective closure number; $S^n=I$ does not hold exactly. This $n$ is a time at which the record minimum of $\|k\alpha\|$ is updated, so by the theorem on best approximations of the second kind [7] it is one of the continued-fraction denominators $q_j$ of $\alpha$. Moreover $n<1/\delta\approx2\pi C/\varepsilon$.

$$
\boxed{\;n=n\bigl(\alpha,\;C/\varepsilon\bigr)\in\{q_j\},\qquad n<\frac{2\pi C}{\varepsilon}\;}
$$

With the law $\alpha$ and the resolution $\varepsilon$ fixed, $n$ and $C$ are not independent degrees of freedom. $n$ is the curvature radius measured by the resolution, quantized to the continued-fraction denominators of the rotation number of the law. The continuous quantity $C$ is unreadable on its own, but with a resolution it becomes readable as the integer $n$. [Goal image] This is the minimal working example of the policy of Sec. 9 of the design document — obtaining discreteness from finite resolution.

[Numerical] (Fig. 2) For the golden-ratio law, $n$ takes only Fibonacci numbers $3,5,8,13,\dots$, and $n/(2\pi C/\varepsilon)$ stays in the theoretical band $[1/\sqrt5,\ \varphi/\sqrt5]\approx[0.447,\ 0.724]$ (observed $0.4472$–$0.7225$). For $\alpha=\pi-3$, $n$ stays at $113$ from $C/\varepsilon\approx18$ up to $\approx5.3\times10^{3}$ — the $355/113$ resonance. For laws whose rotation number is anomalously close to a rational, plateaus appear in which raising the resolution does not increase the closure number.

---

## 7. Numerical Verification

The numerical verification was carried out in a form whose outcome was not predetermined. The commitments were:

- The verified forms are derived from $S$; $S$ is not constructed to preserve the forms being verified. $S$ is a random element of $SL(2,\mathbb R)$.
- Failure conditions are declared before execution in `FAIL_CONDITIONS` at the head of the script.
- No rounding, clipping, or normalization inside the state update (design document Sec. 12).
- The identities are also checked by exact rational arithmetic. The random seed is fixed (20260920).

| Item | Failure condition (preregistered) | Result | Verdict |
|---|---|---|---|
| E1a | exact rational arithmetic ever gives $q_k\neq q_0$ | 0 violations in 300 maps $\times$ 12 steps | PASS |
| E1b | scale-normalized drift of $q$ exceeds $10^{-9}$ | max $2.5\times10^{-13}$ in 3000 maps (1802 elliptic, 1198 hyperbolic) $\times$ 40 steps | PASS |
| E1c | signature of $G$ disagrees with the classification by $\tau$ | 0 of 3000 | PASS |
| E1d | (control) the naive $a^2+b^2$ is conserved in a general description | median fluctuation 122%; above 1% in every elliptic sample | PASS |
| E2a | residual of the undirected law exceeds $10^{-9}$ | 0 violations exactly; max $1.3\times10^{-14}$ in floating point | PASS |
| E2b | for rotations ($n=3..12$), $\lvert PSP-S^{-1}\rvert>10^{-12}$ or $\lvert PSP-S\rvert<10^{-3}$ | former max $3.8\times10^{-16}$; latter min $1.41$ | PASS |
| E2c | a permutation-equivariant real linear map $\gamma I+\beta\mathbf 1\mathbf 1^{\mathsf T}$ ($N=2..8$) has a complex eigenvalue | max imaginary part $3.4\times10^{-15}$ in 2100 maps | PASS |
| E3a | a first-return number is not a continued-fraction denominator | 0 violations in 4 laws $\times$ 241 resolutions | PASS |
| E3b | $n<1/\delta$ is violated | 0 violations | PASS |
| E4a | the two identities of T6 fail in exact arithmetic | 0 violations | PASS |
| **E4b** | in floating point, (discrete curvature)$\times\lvert C\rvert$ deviates from 1 by more than $10^{-7}$ | **max $1.9\times10^{-4}$** | **FAIL** |
| E4c | for two orbits of the same $S$, the ratio of proper lengths deviates from $\lvert C_1/C_2\rvert$ by more than $10^{-7}$ | max $2.4\times10^{-14}$ | PASS |
| E5a | the metric readout error for $n\ge3$ exceeds $10^{-10}$ | max $2.6\times10^{-11}$ for $n=3..12$, 200 trials each (Fig. 3) | PASS |
| E5b | a power sum with $n\nmid p$ exceeds $10^{-10}$ | max $3.9\times10^{-15}$ | PASS |
| E5c | at $n=2$ the self-readout of area is not zero | $\lvert G\rvert=0$, $\operatorname{rank}\sum X_kX_k^{\mathsf T}=1$ | PASS |

**On the failure of E4b.** The tolerance was not loosened afterwards to rewrite it as PASS. The diagnosis: only the hyperbolic class failed; the elliptic maximum error was $4.9\times10^{-14}$. The worst case had $\tau=61.1$ with $\lvert X_k\rvert^2\lVert\hat G\rVert/C^2=2.4\times10^{10}$ — catastrophic cancellation when an indefinite form is evaluated on exponentially grown configurations. The identity itself holds in the exact computation E4a. What failed was therefore not the theorem but the design of the verification.

From this failure, two things were done.

| Item | Character | Failure condition | Result | Verdict |
|---|---|---|---|---|
| E4b' | **post-hoc diagnostic**, not preregistered | evaluation at three points $(X_{-1},X_0,X_1)$ without intervening growth deviates by more than $10^{-9}$ | max $3.4\times10^{-14}$ | PASS |
| E6 | **prediction** built from the cause of E4b and declared before execution | in the hyperbolic class the slope of $\ln\bigl(\lvert q\rvert/\lvert X_k\rvert\lvert X_{k+1}\rvert\bigr)$ deviates from $-2H$ by more than 1% | max $0.05\%$ over 300 laws | PASS |

E6 is the verification of T8's decay rate $-2H$. A failed verification taught the statement that in the hyperbolic class the conserved quantity becomes exponentially unreadable.

**Re-execution in a different environment.** The numbers in the table above and the bundled JSON and figures come from a single run (Linux x86_64, Python 3.12.3, NumPy 2.4.4, OpenBLAS 0.3.31). Re-running the script in a different environment (macOS arm64, Python 3.9.6, NumPy 2.0.2, Accelerate) reproduced all 17 verdicts (only E4b FAIL). The exact rational computations (E1a, the exact part of E2a, E4a) and the integer-valued results (E3a, E3b, the list of closure numbers) agreed exactly; the floating-point verification values changed only in their trailing digits (e.g., the E1b maximum $2.5\times10^{-13}\to$ $3.0\times10^{-13}$; the E4b maximum $1.9\times10^{-4}\to$ $1.4\times10^{-4}$). Floating-point verification values depend on the linear-algebra implementation, so what is reproduced is the verdicts, not the numbers themselves. E4b failed in both environments, consistent with the cancellation diagnosis.

![Fig. 1: three classes and two readouts (labels in Japanese)](fig01_area_readout_three_classes_ja_v1.svg)

**Fig. 1** (a) An elliptic orbit becomes an ellipse in a general description. (b) Hyperbolic class: the family of orbits of the same $S$ with different $C$, and the null cone. (c) Parabolic class. (d) Along the orbit of (a), the naive $a^2+b^2$ versus the area readout $\omega(X_k,X_{k+1})$. (Figure labels are in Japanese.)

![Fig. 2: closure numbers under finite resolution (labels in Japanese)](fig02_area_readout_closure_staircase_ja_v1.svg)

**Fig. 2** For four laws ($\alpha=(\sqrt5-1)/2,\ \sqrt2-1,\ e-2,\ \pi-3$), the first-return number $n$ against $\rho=C/\varepsilon$. The dashed line is the bound $1/\delta$. (Figure labels are in Japanese.)

![Fig. 3: reading the metric from second moments (labels in Japanese)](fig03_area_readout_metric_from_orbit_ja_v1.svg)

**Fig. 3** Relative error between $\sum_kX_kX_k^{\mathsf T}$ of closed orbits and $(nC^2/2)\hat G^{-1}$. (Figure labels are in Japanese.)

---

## 8. Returning to the Design Document (the Zeroth Paper)

### 8.1 Status of the five feasibility conditions (design document Sec. 11)

| Condition | What this paper can say |
|---|---|
| 1. Can a global conserved quantity be maintained? | Yes: the area readout $q$. Conservation of the signed quantity is equivalent to $\det S=1$ (T2), and then the law does not distinguish a direction (T3'). On this branch there is thus no arrow at the level of the law; arrows can only come from individual orbits (remark (iii) of Sec. 4.1). The quantity is on the side of relations, not a per-state quantity $\sum C(z_i)$ |
| 2. Can interaction including self-interaction be defined? | Yes, as the self-interaction of the two-degree-of-freedom system. But what can be defined on nameless configurations is the second-order undirected law (T3) |
| 3. Can the conservation structure be maintained as the number of states changes? | Not treated in this paper (fixed at two degrees of freedom); Sec. 8.3 gives only a lead |
| 4. Does it stay finite? | The elliptic class is bounded. The hyperbolic class is unbounded, and by T8 the conserved quantity becomes unreadable |
| 5. Are relational directions generated without external axes? | Partially: the signature of the metric and the two null directions of the hyperbolic class are generated by the interaction |

### 8.2 On the skeleton of Sec. 6 of the design document [Derived T11]

(In this and the next subsection $N$ denotes the number of states — distinct from the matrix $N=(S-S^{-1})/2$ of Sec. 4.1.)

Implement the skeleton $z_{i,n+1}=\sum_jf(z_{i,n},z_{j,n})$ of Sec. 6 of the design document with real, linear $f(x,y)=px+qy$. Summing including the self-term,

$$
z_i'=Np\,z_i+q\sum_jz_j,\qquad S=Np\,I+q\,\mathbf 1\mathbf 1^{\mathsf T}
$$

The eigenvalues are $Np$ (multiplicity $N-1$) and $N(p+q)$ — all real. For $N=2$,

$$
S=\begin{pmatrix}2p+q&q\\q&2p+q\end{pmatrix}
$$

Therefore:

- From a linear, real, first-order skeleton, the elliptic class (closed orbits, phase, oscillation) does not arise for any $N$.
- Imposing A1 at $N=2$ gives $\det S=4p(p+q)=1$, and $S$ is a boost ($\cosh H=2p+q$, $\sinh H=q$). Exponential spreading (design document Sec. 7.5) and the invariant $a^2-b^2$ come out of this skeleton with no additional assumptions.
- There are three routes to closed structure (the condensation and particle picture of Secs. 7.7–7.8 of the design document): (i) allow complex coefficients — this inputs $i$; (ii) make the law a second-order undirected law — then $i$ is not an input but appears as the $J$ of T5; (iii) go nonlinear.

This is an instance of what Sec. 13 of the design document calls "outcomes about the relations between requirements."

### 8.3 Lifting the undirected law to $N$ states [Derived T12] [Goal image]

Lifting the second-order undirected law to $N$ states, the permutation-equivariant linear case is

$$
X_{k+1}+X_{k-1}=(\gamma I+\beta\,\mathbf 1\mathbf 1^{\mathsf T})X_k
$$

The uniform mode has $\tau_u=\gamma+N\beta$; the $N-1$ relative modes have $\tau_r=\gamma$; each falls independently into the three classes of T5. From this:

- If $\lvert\gamma\rvert<2$ and $\beta>0$, the relative modes stay bounded while only the uniform mode moves into the hyperbolic class for $N>N^{*}=(2-\gamma)/\beta$. What moves the class boundary is the number of states $N$.
- In a linear system with fixed $N$ and law, $\tau$ is constant and no class transition occurs. The phase transition of Sec. 7.7 of the design document requires either a change of $N$ or nonlinearity.
- The pair quantities $W_{ij}(k)=x_{i,k}x_{j,k+1}-x_{i,k+1}x_{j,k}$ satisfy $W_{ij}(k)-W_{ij}(k-1)=\beta\bigl(\sum_mx_{m,k}\bigr)(x_{i,k}-x_{j,k})$. Since the law is second order, a state is a pair of consecutive times $(X_{k-1},X_k)$. If $\sum_mx_m=0$ at those two times, the law gives $\sum_mx_{m,k+1}=\tau_u\cdot0-0=0$, so $\sum_mx_m=0$ persists forever. On this zero closure of the linear sum, all $W_{ij}$ are conserved. There the all-to-all sum interaction becomes ineffective and each state follows the same law $x_{k+1}+x_{k-1}=\gamma x_k$; hence for $\lvert\gamma\rvert<2$ the dynamics on the zero closure is bounded independently of $N$.

The spectral and conservation statements above are derivations. Reading them against Secs. 7.2, 7.5, 7.7 of the design document is a correspondence with the goal images, to be tested in the next thought experiment.

---

## 9. Relation to Known Mathematics, and What Is New

As mathematics this is classical. The group of linear maps preserving $\omega$ is $SL(2,\mathbb R)=Sp(2,\mathbb R)$; viewing $(a,b)$ as a canonically conjugate pair $(q,p)$, the elliptic class is the harmonic oscillator, the hyperbolic class the inverted oscillator, the parabolic class the free particle [8]. The correspondence of the three classes with $i$, $\epsilon$, $j$ and with the three plane geometries — Euclidean, Galilean, Minkowskian — is also known [6]. Hyperbolic motion and the Rindler frame are likewise known [3,4].

The closest existing framework is centro-affine geometry: the geometry of plane curves under linear transformations fixing the origin, in which the determinant of the position vector and the tangent vector is the basic invariant. Choosing the parameter with $\det(\gamma,\gamma_\sigma)=1$ and differentiating gives $\det(\gamma,\gamma_{\sigma\sigma})=0$, i.e. $\gamma_{\sigma\sigma}=-\mu\gamma$; the curves with constant $\mu$ are the origin-centered ellipses, hyperbolas, and straight lines [9,10]. In the discrete version, two discrete curvatures are defined as ratios of determinants of edge vectors $t_k=r_{k+1}-r_k$: $\kappa_k=[t_k,t_{k+1}]/[t_{k-1},t_k]$ and $\bar\kappa_k=[t_{k-1},t_{k+1}]/[t_{k-1},t_k]$, with $t_{k+1}=-\kappa_kt_{k-1}+\bar\kappa_kt_k$. A constant-curvature curve closes iff $\kappa=1$, $\bar\kappa=2\cos\theta$, $p\theta=2l\pi$ — an affinely regular polygon [11]. The undirected law of this paper is the case $\kappa_k\equiv1$, $\bar\kappa_k\equiv\tau$, and T4's formula for $\tau$ has the same form as the formula for $\bar\kappa$. Thus the form of the equations of T3 and T4, and the fact that closed orbits are affinely regular polygons, are known.

This correspondence gives P1 a reading: the nameless pair of state and relation corresponds, in the language of mechanics, to a canonically conjugate pair. "Which is the state and which is the relation cannot be decided" is "which is position and which is momentum cannot be decided," and the fact that the swap acts as time reversal matches as well.

What this paper takes to be new:

1. That the area form follows uniquely from the single principle "a state alone is unreadable" (T1). The reason for not choosing a metric becomes a principle, not a preference.
2. The sorting-out of namelessness versus the direction of iteration: that the only law definable on nameless configurations is the second-order rule of a single real $\tau$, and that the first readable quantity there is the acceleration (T3, T4). The form of the equations is known as the constant-curvature case of discrete centro-affine geometry [11]; what is new is the derivation side — that it is forced by P1, P2, A1.
3. That the right-hand side of $\sum z_n^2=C^2$ equals exactly the curvature radius in the induced metric, at any step size and in both signatures (T6); together with the comparison reading the two scenes of the equivalence principle as two classes (C1–C3).
4. The two contents of zero closure (T7, T9) and the projective convergence in the hyperbolic class (T8).
5. That the closure number is determined by $C/\varepsilon$ and the continued fraction of the law (T10).
6. The consequences for the skeleton of the design document (T11, T12).

---

## 10. Limitations

- A1 remains an axiom candidate; why a readout should be conserved is not derived. The condition "the law does not distinguish a direction" is weaker than conservation of the signed $q$ and stronger than conservation of the readable $\lvert q\rvert$ (T3'), so it is not a paraphrase of A1. Moreover, this paper treated only the branch $\det S=1$. The orientation-reversing branch $\det S=-1$ (minimal example: the swap $P$ itself) was only checked to conserve $\lvert q\rvert$ and to fall into the hyperbolic class when viewed every two steps; nothing further was examined. Whether A1 should be read as conservation of $q$ or of $\lvert q\rvert$ is also left undecided.
- D2 (linearity) is an assumption. As long as the law is linear, neither class transitions nor the inverse-square law arise.
- That there are two values follows from P1 (Sec. 2.1). What is not derived is that the two values must be real (D1). T1–T4 and the identities of T5, T6 hold without reality (E1a and E4a are exact rational computations). Where reality is essential is the distinction of the three classes: the elliptic and hyperbolic classes are distinct because $-1$ is not a square among the values. If $a,b$ themselves were complex, then $b\to ib$ would identify $a^2+b^2$ with $a^2-b^2$ and the sign $s$ would lose its meaning (and $a^2+b^2=0$ would acquire the nontrivial solutions $a=\pm ib$). Hence deriving reality is the same problem as deriving why the signature distinction exists.
- The comparisons C1–C3 reach only 1+1-dimensional kinematics. How relative readouts compose when two systems are coupled (the analogue of the velocity-composition law) is the next verification, which could break the comparison. Also, the physical validity of the correspondence of C2, $(a,b)\leftrightarrow(x,ct)$, remains a point that cannot be checked by verifying formulas.
- The fixed resolution D4 is a stand-in. How T10's $n(C/\varepsilon)$ changes when the resolution is determined by the interaction is unexamined.
- Change of the number of states (design document Sec. 7.2) is not treated.

**Candidates for the next thought experiment.** (i) Does the minimal nonlinearization, with $\kappa$ depending on the configuration, produce a transition from the hyperbolic to the elliptic class? (ii) Does the composition law of relative readouts of two coupled two-degree-of-freedom systems remain consistent with C2? (iii) Can the resolution be defined endogenously from ratios of $\omega$?

---

## 11. Working Record

Recorded in accordance with Sec. 14 of the design document: "preserve why the computation was performed."

- **Record of verification.** The failure of E4b was kept without rewriting. E4b' is a post-hoc diagnostic; E6 is a prediction declared before execution. The distinction is stated explicitly in the tables.
- **Reproducibility.** The scripts use only relative paths, and output file names match the text. The SVG figures are always generated by a dependency-free simple renderer. PNG figures are generated only when matplotlib is available; Japanese fonts are searched among several candidates, with fallback to English labels if none is found. The symbolic verification of the identities is kept in a separate script.

---

## References

### This series

[0] N. Kihara, **Designing a Research Project to Explore the Emergence of Spacetime, Particles, and Interactions from States and Relations — A Framework of Preliminary Design and Feasibility Verification for Foundational-Physics Model Exploration under Uncertainty**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22851944.

### External

[1] Einstein, A., **Über das Relativitätsprinzip und die aus demselben gezogenen Folgerungen**, *Jahrbuch der Radioaktivität und Elektronik*, **4** (1907), 411–462.

[2] Einstein, A., **Die Grundlage der allgemeinen Relativitätstheorie**, *Annalen der Physik*, **354**(7) (1916), 769–822. DOI: 10.1002/andp.19163540702.

[3] Born, M., **Die Theorie des starren Elektrons in der Kinematik des Relativitätsprinzips**, *Annalen der Physik*, **30** (1909), 1–56.

[4] Rindler, W., **Kruskal Space and the Uniformly Accelerated Frame**, *American Journal of Physics*, **34** (1966), 1174–1178.

[5] Norton, J. D., **What Was Einstein's Principle of Equivalence?**, *Studies in History and Philosophy of Science Part A*, **16**(3) (1985), 203–246. DOI: 10.1016/0039-3681(85)90002-0.

[6] Yaglom, I. M., **A Simple Non-Euclidean Geometry and Its Physical Basis**, Springer, 1979.

[7] Khinchin, A. Ya., **Continued Fractions**, University of Chicago Press, 1964. (On best approximations of the second kind and convergents.)

[8] Arnold, V. I., **Mathematical Methods of Classical Mechanics**, 2nd ed., Springer, 1989.

[9] Olver, P. J., **Moving Frames and Differential Invariants in Centro-Affine Geometry**, *Lobachevskii Journal of Mathematics*, **31**(2) (2010), 77–89. DOI: 10.1134/S1995080210020010.

[10] Qu, C. Z., Yang, Y., **An invariant second-order curve flow in centro-affine geometry**, *Journal of Geometry and Physics*, **174** (2022), 104447.

[11] Yang, Y., Yu, Y., **Moving frame and integrable system of the discrete centroaffine curves in $\mathbb R^3$**, arXiv:1601.06530v2 (2016). (The title of v1 was: The Frenet–Serret formulas of a discrete centroaffine curve.)

---

## Appendix A: Proofs

**A.1 (T1)** If $B(X,X)=0$ for all $X$, then $0=B(X+Y,X+Y)=B(X,Y)+B(Y,X)$, so $B$ is antisymmetric. Antisymmetric bilinear forms on $\mathbb R^2$ form a one-dimensional space: $B=c\,\omega$.

**A.2 (T2)** Since $\omega(X,Y)=\det[X\ Y]$, we get $\omega(SX,SY)=\det(S[X\ Y])=\det S\cdot\omega(X,Y)$, hence $q_{k+1}=\omega(SX_k,SX_{k+1})=\det S\cdot q_k$. If $S$ is not a scalar multiple of the identity, then for any $X_0$ that is not an eigenvector, $q_0=\omega(X_0,SX_0)\neq0$; therefore conservation of the signed $q$ for all initial configurations is equivalent to $\det S=1$ (for conservation of $\lvert q\rvert$: $\lvert\det S\rvert=1$).

**A.3 (T3)** By Cayley–Hamilton for $2\times2$ matrices, $S^2-\tau S+(\det S)I=0$. For $\det S=1$, multiplying by $S^{-1}$ gives $S+S^{-1}=\tau I$; acting on $X_k$ yields the undirected law. Conversely, for independent $X_0,X_1$, defining $S$ by $SX_0=X_1$, $SX_1=\tau X_1-X_0$ gives, in this basis, the matrix $\begin{pmatrix}0&-1\\1&\tau\end{pmatrix}$ with $\det=1$, $\operatorname{tr}=\tau$.

**A.3' (T3')** Set $d:=\det S$. For general $S$, $S^2-\tau S+dI=0$; acting on $X_{k-1}$ gives $X_{k+1}-\tau X_k+d\,X_{k-1}=0$. Replacing $k$ by $-k$, dividing by $d$ and rearranging, the law obeyed by $Y_k=X_{-k}$ is $Y_{k+1}-(\tau/d)Y_k+d^{-1}Y_{k-1}=0$. The residual of $Y$ in the original law is $\frac{d-1}{d}\bigl[-\tau Y_k+(d+1)Y_{k-1}\bigr]$, which vanishes for all configurations exactly when $d=1$, or $d=-1$ and $\tau=0$; in the latter case $S^2=\tau S-dI=I$. For the readout, A.2 gives $q_{k+1}=d\,q_k$. On the branch $\det S=-1$: $\operatorname{tr}S^2=\tau^2-2d=\tau^2+2$; the matrix of the induced form is $\Omega(S-dS^{-1})/2$ with determinant $d-\tau^2/4<0$; and $\omega(SX,S^2X)=d\,\omega(X,SX)$.

**A.4 (T4)** Substituting $X_{k+1}=\tau X_k-X_{k-1}$ into $\omega(X_{k-1},\cdot)$ gives $\omega(X_{k-1},X_{k+1})=\tau\,\omega(X_{k-1},X_k)$. The second difference is a rewriting of the undirected law.

**A.5 (T5)** $\det S=1$ is equivalent to $S^{\mathsf T}\Omega S=\Omega$, hence $S^{\mathsf T}\Omega=\Omega S^{-1}$. The symmetric part of $\Omega S$ is $(\Omega S-S^{\mathsf T}\Omega)/2=\Omega(S-S^{-1})/2=G$. From $S-S^{-1}=2S-\tau I$, $\det(S-S^{-1})=4\chi_S(\tau/2)=4-\tau^2$, so $\det G=1-\tau^2/4$. From $S^2+S^{-2}=(\tau^2-2)I$, $N^2=(S^2-2I+S^{-2})/4=(\tau^2/4-1)I$. Conservation follows from $Q_S(SX)=\omega(SX,S\,SX)=\omega(X,SX)$.

**A.6 (T6)** From $S^{\mathsf T}GS=G$ and $S^{\mathsf T}G=GS^{-1}$: $(S-I)^{\mathsf T}G(S-I)=2G-G(S^{-1}+S)=(2-\tau)G=\kappa G$. Hence $Q_S(\Delta X)=\kappa Q_S(X)$, and from $\Delta^2X=-\kappa X$, $Q_S(\Delta^2X)=\kappa^2Q_S(X)$. Normalize by $\hat G$, take absolute square roots, and form the ratio.

**A.7 (T7)** $\omega(X,SX)=\det[X\ SX]=0$ iff $X$ and $SX$ are linearly dependent.

**A.8 (T8)** In the hyperbolic class $S$ has real eigenvalues $\lambda_\pm$ ($\lvert\lambda_\pm\rvert=e^{\pm H}$) with eigendirections $v_\pm$. Decomposing $X_k=A_kv_++B_kv_-$: $A_k=\lambda_+^kA_0$, $B_k=\lambda_-^kB_0$. Since $\omega(X_k,v_-)=A_k\,\omega(v_+,v_-)$ and $\omega(v_+,X_k)=B_k\,\omega(v_+,v_-)$, for $A_0\neq0$: $\rho_k=\lvert B_k/A_k\rvert=\rho_0e^{-2Hk}$. Also $\lvert X_k\rvert\sim\lvert A_0\rvert e^{Hk}\lvert v_+\rvert$, so $\lvert q\rvert/(\lvert X_k\rvert\lvert X_{k+1}\rvert)$ decays at the same rate $e^{-2Hk}$.

**A.9 (T9)** In the normal description let $X_k=C(\cos\varphi_k,\sin\varphi_k)$, $\varphi_k=\varphi_0+2\pi mk/n$, $\gcd(m,n)=1$. For $n\ge3$, $n\nmid2m$, so $\sum_ke^{2i\varphi_k}=0$ and $\sum X_kX_k^{\mathsf T}=(nC^2/2)I$. In a general description $X=TY$ ($\hat G=T^{-\mathsf T}T^{-1}$): $\sum X_kX_k^{\mathsf T}=(nC^2/2)TT^{\mathsf T}=(nC^2/2)\hat G^{-1}$. The power sums are geometric series.

**A.10 (T10)** In the induced metric, $\lvert X_n-X_0\rvert=2C\lvert\sin(n\theta/2)\rvert$. At the smallest $n$ satisfying the condition, all earlier $\|k\alpha\|$ are at least $\delta$, so $n$ updates the record minimum of $\|k\alpha\|$. Record-update times are restricted to denominators of convergents [7]. For the previous update time $q_{j-1}$: $\delta\le\|q_{j-1}\alpha\|<1/q_j$, hence $n=q_j<1/\delta$.

**A.11 (T11)** The permutation representation $\mathbb R^N$ splits into the trivial and an $(N-1)$-dimensional irreducible representation, so linear maps commuting with permutations are spanned by $I$ and $\mathbf 1\mathbf 1^{\mathsf T}$. The eigenvalues of $\gamma I+\beta\mathbf 1\mathbf 1^{\mathsf T}$ are $\gamma+N\beta$ along $\mathbf 1$ and $\gamma$ on its orthogonal complement.

**A.12 (T12)** In the eigendecomposition each mode takes the form of T3. Substituting the law into $W_{ij}(k)-W_{ij}(k-1)=x_{i,k}(x_{j,k+1}+x_{j,k-1})-x_{j,k}(x_{i,k+1}+x_{i,k-1})$, the $\gamma$-terms cancel, leaving $\beta\bigl(\sum_mx_{m,k}\bigr)(x_{i,k}-x_{j,k})$. For the sum, $\sum_mx_{m,k+1}=\tau_u\sum_mx_{m,k}-\sum_mx_{m,k-1}$, so if the sum vanishes at two consecutive times it vanishes forever.

---

## Appendix B: Files and Reproduction

| File | Content |
|---|---|
| `run_area_readout_experiments_v1.py` | all verifications E1–E6 and figure generation; `FAIL_CONDITIONS` at the head declares the failure conditions |
| `area_readout_experiments_results_v1.json` | execution results (verdicts, values, list of closure numbers) |
| `verify_area_readout_identities_sympy_v1.py` | symbolic check of the identities of Appendix A (T2, T3, T3', T5, T6, T11, and the normal-form readouts) |
| `fig01_area_readout_three_classes_ja_v1.svg` | Fig. 1 |
| `fig02_area_readout_closure_staircase_ja_v1.svg` | Fig. 2 |
| `fig03_area_readout_metric_from_orbit_ja_v1.svg` | Fig. 3 |

Run with `python3 run_area_readout_experiments_v1.py`; dependencies are numpy and mpmath. Output is written to the same folder as the script. If matplotlib is available, PNG versions of the figures are also generated under the same names. Symbolic verification: `python3 verify_area_readout_identities_sympy_v1.py` (dependency: sympy).
