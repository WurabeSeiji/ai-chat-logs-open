# Noether Conservation Laws and Relational Phase Dynamics from Discrete Zero Closure — A Zero-Closure-Preserving Discrete Self-Map, the $N\to\infty$ Continuum Field Equations, Local Gauge Geometry, and the Standard-Model One-Generation Representation with Chirality Selection

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** 21 August 2026<br>
**Version DOI:** 10.5281/zenodo.22040736<br>
**Concept DOI:** 10.5281/zenodo.22040735<br>
**Position in the series:** "Generative Structure of Dimensions" series — sequel to *Symmetry Generation from Closure Axioms* (dynamics, conservation laws, Standard-Model representation), public version v1.0<br>
**Preceding paper:** Symmetry Generation from Zero Closure, Finite Order, and Self-Consistent Geometry v1.0 (Concept DOI 10.5281/zenodo.22028072)<br>
**License:** CC BY 4.0

> **Subject.** This paper addresses the two problems left open by the preceding paper: Noether-type conservation laws, and the connection to a dynamics that determines the next state. The central claim on dynamics is that state rewriting is not admitted implicitly; a dynamics must be constructed as a self-map of the admissible state space. We construct a discrete self-map that preserves zero closure exactly at every finite iteration and obtain the standard continuum field equations in its $N\to\infty$ limit. We further give the theoretical connection to local gauge geometry, the Standard-Model one-generation representation, hypercharge, anomaly cancellation, and chirality selection, together with a numerical verification specification.

## Abstract

In the preceding paper, "Must the symmetries of physics really be given from the start?", many symmetry structures were derived and organised without presupposing a background spacetime or a ready-made symmetry group, from the few conditions

$$
\sum_{a=1}^{M}X_a^2=0,
\qquad
U^N=I,
$$

together with simplex closure and self-consistency. What remained were the local conservation laws corresponding to Noether's theorem, and the connection to a dynamics that determines the next state.

In this paper, under equal amplitude

$$
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.}
$$

we obtain from the relational phase differences

$$
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
$$

an oriented discrete current. Constructing further the discrete action

$$
S_N[\phi]
=
-A^2\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i)
$$

we find

$$
\frac{\partial S_N}{\partial\phi_i}
=
-\sum_{j\sim i}J_{ij}
$$

so that the stationarity condition is

$$
\sum_{j\sim i}J_{ij}=0
$$



For the dynamics, instead of using a naive phase update as it stands, we construct a constrained self-map that preserves the zero closure

$$
C(\phi)
=
\sum_i e^{2i\phi_i}
=
0
$$

exactly. Zero closure is equivalent to the two real conditions

$$
C_R(\phi)=\sum_i\cos2\phi_i=0,
\qquad
C_I(\phi)=\sum_i\sin2\phi_i=0
$$

. Accordingly, projecting the unconstrained relational phase force

$$
F_i(\phi)
=
\sum_{j\sim i}\sin(\phi_j-\phi_i)
$$

onto the tangent space of the zero-closure manifold and setting

$$
\dot\phi
=
P_\phi F(\phi)
$$

preserves zero closure exactly along a continuous map parameter. For finite iterations we use a retraction that brings this tangent vector back onto the zero-closure manifold, and define

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}
\!\left(
\eta P_{\phi^{(n)}}F(\phi^{(n)})
\right)
}
$$

. With this,

$$
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N
$$

holds at every iteration.

Raising the readout resolution $N$ and letting the neighbour spacing $h_N\to0$, the discrete current goes over to

$$
J_N^\mu
\longrightarrow
A^2\partial^\mu\phi
$$

and the discrete continuity equation converges to

$$
\boxed{
\partial_\mu J^\mu=0
}
$$



Since the discrete action itself goes over in the continuum limit to

$$
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq
$$

the Euler–Lagrange equation becomes

$$
\boxed{
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi
\right)
=
0
}
$$

. If $A$ and $g^{\mu\nu}$ are constant,

$$
\boxed{
\Box_g\phi=0
}
$$



Here we state explicitly the central claim of this paper concerning dynamics. **Rewriting a state is not, by itself, implicitly admitted as dynamics.** Since self-consistency is a foundational condition of this axiom system, the map to the next state must be a self-map that does not destroy the admissible state space. In particular, for the zero closure treated explicitly in this paper, we construct, while updating only the phases $\phi$,

$$
\boxed{
\mathcal F_N:\mathcal Z_N\to\mathcal Z_N
}
$$

and show that at every finite iteration

$$
\sum_iX_i^2=0
$$

is preserved exactly. Hence **the dynamics of the discrete system can be realised as an admissible phase map that keeps zero closure, rather than by breaking the axiom once and re-imposing the constraint afterwards.** This is not a premise of the derivation of dynamics in this paper but a result to be verified. For the other previously derived conditions, including finite recurrence, simplex closure, and self-consistency, the subsequent numerical verification treats "preserved after the update" as an independent audit item, and any update rule that fails to preserve them is rejected as a dynamics of this axiom system.

Therefore, from the discrete relational phase system, this paper gives as a single derivation chain a self-map that preserves zero closure exactly, Noether-type local conservation laws, and the continuum limit to the standard partial differential field equations.

Furthermore, connecting the five complex degrees of freedom and $S(U(3)\times U(2))$ of the preceding paper to the local gauge dynamics of this paper, the trace-zero condition on the $U(1)$ generator fixes the hypercharge ratio, and from the five degrees of freedom and the simplex two-body relations

$$
V^*\oplus\Lambda^2V
$$

appears. This representation decomposes into the 15 left-handed Weyl components

$$
d^c,\ L,\ u^c,\ Q,\ e^c
$$

of one Standard-Model generation, and the cancellation of all perturbative gauge anomalies and of the $SU(2)$ global anomaly can be checked directly. As for the selection of the conjugate two Weyl sectors, the last remaining step in the derivation chain of the Standard-Model internal representation, we identify the A/B two-state selection system of our earlier papers with the conjugate Weyl sectors, and from the mirror-odd internal correlation
$$
J=\operatorname{Im}(B^{*2}C)
$$
and the existing nonlinear selection term we obtain the minimal normal form
$$
\dot S_\chi=\lambda J+gS_\chi(1-S_\chi^2)
$$
. With this, the theoretical connection targeted by this paper — gauge group, local gauge geometry, one-generation internal representation, hypercharge, anomaly cancellation, and chirality selection — closes. What remains is the numerical verification of whether this identification and selection rule hold on the existing numerical dynamics.

---

## 1. Motivation

The preceding paper [1] derived many symmetries from a few closure conditions.

Two problems, however, remained.

First,

$$
\boxed{
\text{the symmetries were derived, but}
\text{the Noether conservation laws were not}
}
$$



Second,

$$
\boxed{
\text{no connection was made from the self-consistent structure}
\text{to a dynamics that determines the next state}
}
$$



The present work examined these two points.

The second problem carries a strict condition peculiar to this axiom system. In ordinary field theory one may first posit an update rule from a state $\Phi$ to the next state $\Phi'$ and treat that update as "time evolution". In this axiom system, however, **self-consistency itself is a foundational condition. One therefore cannot assume state rewriting unconditionally.** If the updated state leaves zero closure, finite recurrence, simplex closure, or any other admissibility condition, then that update rule is not a dynamics of this theory.

Hence the question of dynamics to be asked in this paper is not simply whether one can give

$$
\phi^{(n)}\mapsto\phi^{(n+1)}
$$

but whether one can construct

$$
\boxed{
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N
}
$$

— a **self-map inside the admissible state space**.

As a result of the investigation, the Noether-type conservation law appeared as the $N\to\infty$ continuum limit of the finite-$N$ discrete current conservation law. Moreover, the same relational phase current turned out to be the generator of the self-map that determines the next phases. Further, by using the tangent-space projection of the zero-closure manifold and a finite retraction, one can update the phases $\phi$ while preserving zero closure exactly at every finite step. **"Being able to produce a next state" and "being able to produce a next state without breaking the axioms" are different problems; what this paper closes is the latter.**

However, for the naive update

$$
\phi_i^{(n+1)}
=
\phi_i^{(n)}
+
\eta F_i
$$

exact preservation of zero closure is not automatic.

This paper therefore introduces a tangent-space projection and a finite retraction that preserve zero closure, and closes the dynamics as one satisfying

$$
\boxed{
\mathcal F:
\mathcal Z_N\to\mathcal Z_N
}
$$



---

## 2. Basic States and Zero Closure

### 2.1 Equal-Amplitude States

$$
\boxed{
X_i=Ae^{i\phi_i},
\qquad
A=\mathrm{const.}
}
\tag{1}
$$



### 2.2 Zero Closure

The first axiom is

$$
\sum_iX_i^2=0.
$$

Substituting (1),

$$
A^2\sum_i e^{2i\phi_i}=0.
$$

Since $A\neq0$,

$$
\boxed{
C(\phi)
:=
\sum_i e^{2i\phi_i}
=
0
}
\tag{2}
$$



Separating real and imaginary parts,

$$
\boxed{
C_R(\phi)
=
\sum_i\cos2\phi_i
=
0
}
\tag{3}
$$

$$
\boxed{
C_I(\phi)
=
\sum_i\sin2\phi_i
=
0
}
\tag{4}
$$



The admissible phase set is therefore

$$
\boxed{
\mathcal Z_N
=
\left\{
\phi
\in
(\mathbb R/2\pi\mathbb Z)^M
\ \middle|\
C_R(\phi)=0,
\ C_I(\phi)=0,
\ U^N=I
\right\}.
}
\tag{5}
$$

---

## 3. Relational Phase Current

For adjacent states $i,j$,

$$
X_i^*X_j
=
A^2e^{i(\phi_j-\phi_i)}.
\tag{6}
$$

We define the oriented edge current as

$$
\boxed{
J_{ij}
=
\operatorname{Im}(X_i^*X_j)
=
A^2\sin(\phi_j-\phi_i)
}
\tag{7}
$$



Under edge reversal,

$$
\boxed{
J_{ji}=-J_{ij}.
}
\tag{8}
$$

Moreover, $J_{ij}$ is invariant under

$$
\phi_i\mapsto\phi_i+\alpha
$$



---

## 4. Discrete Action and Discrete Noether Conservation Law

From the real part on each edge we obtain

$$
s_{ij}
=
-\operatorname{Re}(X_i^*X_j)
=
-A^2\cos(\phi_j-\phi_i)
\tag{9}
$$



The total action is

$$
\boxed{
S_N[\phi]
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i).
}
\tag{10}
$$

Under the global phase transformation

$$
\phi_i\mapsto\phi_i+\alpha
$$

the phase differences do not change, so

$$
\boxed{
S_N[\phi+\alpha]=S_N[\phi].
}
\tag{11}
$$

Varying with respect to the vertex phase $\phi_i$,

$$
\frac{\partial S_N}{\partial\phi_i}
=
-A^2
\sum_{j\sim i}
\sin(\phi_j-\phi_i).
$$

Hence

$$
\boxed{
\frac{\partial S_N}{\partial\phi_i}
=
-\sum_{j\sim i}J_{ij}.
}
\tag{12}
$$

The stationarity condition

$$
\frac{\partial S_N}{\partial\phi_i}=0
$$

becomes

$$
\boxed{
\sum_{j\sim i}J_{ij}=0
}
\tag{13}
$$



Defining the discrete divergence as

$$
\boxed{
(\operatorname{div}_dJ)_i
=
\sum_{j\sim i}J_{ij}
}
\tag{14}
$$

we have

$$
\boxed{
\operatorname{div}_dJ=0.
}
\tag{15}
$$

---

## 5. The Unconstrained Minimal Relational Phase Force

We define the negative gradient of the action as

$$
F_i(\phi)
:=
-\frac{1}{A^2}
\frac{\partial S_N}{\partial\phi_i}
$$



By (12),

$$
\boxed{
F_i(\phi)
=
\sum_{j\sim i}
\sin(\phi_j-\phi_i).
}
\tag{16}
$$

This is

$$
\boxed{
\text{the present relational phase differences determine the next phase change}
}
$$

— the minimal local generator.

However,

$$
\phi^{(n+1)}
=
\phi^{(n)}
+
\eta F(\phi^{(n)})
$$

alone does not, in general, guarantee

$$
C(\phi^{(n+1)})=0
$$



This is closed in the next section.

---

## 6. Tangent Space of the Zero-Closure Manifold

The two real conditions of zero closure are (3) and (4).

Their gradients are

$$
\frac{\partial C_R}{\partial\phi_i}
=
-2\sin2\phi_i,
$$

$$
\frac{\partial C_I}{\partial\phi_i}
=
2\cos2\phi_i.
$$

Hence

$$
\boxed{
g_R(\phi)
=
(-2\sin2\phi_1,\ldots,-2\sin2\phi_M)^T
}
\tag{17}
$$

$$
\boxed{
g_I(\phi)
=
(2\cos2\phi_1,\ldots,2\cos2\phi_M)^T.
}
\tag{18}
$$

A phase-change vector $v\in\mathbb R^M$ is tangent to the zero-closure manifold if and only if

$$
g_R^Tv=0,
\qquad
g_I^Tv=0.
\tag{19}
$$

This directly means

$$
\frac{dC_R}{ds}=0,
\qquad
\frac{dC_I}{ds}=0
$$



---

## 7. Projecting the Unconstrained Force onto the Tangent Space

Define the $M\times2$ matrix

$$
\boxed{
G(\phi)
=
\begin{pmatrix}
|&|\\
g_R&g_I\\
|&|
\end{pmatrix}
}
\tag{20}
$$



At regular points where $g_R,g_I$ are independent,

$$
G^TG
$$

is invertible.

The orthogonal projection onto the tangent space is

$$
\boxed{
P_\phi
=
I
-
G(G^TG)^{-1}G^T.
}
\tag{21}
$$

Indeed,

$$
G^TP_\phi
=
G^T
-
G^TG(G^TG)^{-1}G^T
=
0.
$$

Hence

$$
\boxed{
G^TP_\phi=0.
}
\tag{22}
$$

From the unconstrained force $F$ we define the tangential force that does not break zero closure as

$$
\boxed{
F_T(\phi)
=
P_\phi F(\phi)
}
\tag{23}
$$



Then

$$
g_R^TF_T=0,
\qquad
g_I^TF_T=0.
\tag{24}
$$

---

## 8. Exact Preservation of Zero Closure along a Continuous Self-Map Parameter

Physical time is not yet introduced here.

Using a bare self-map parameter $s$, set

$$
\boxed{
\frac{d\phi}{ds}
=
F_T(\phi)
=
P_\phi F(\phi)
}
\tag{25}
$$



Then

$$
\frac{dC_R}{ds}
=
\nabla C_R^T
\frac{d\phi}{ds}
=
g_R^TP_\phi F
=
0.
$$

Similarly,

$$
\frac{dC_I}{ds}
=
g_I^TP_\phi F
=
0.
$$

Hence

$$
\boxed{
\frac{dC_R}{ds}
=
\frac{dC_I}{ds}
=
0.
}
\tag{26}
$$

If the initial state satisfies zero closure,

$$
C_R(\phi(0))=C_I(\phi(0))=0
$$

so for every $s$,

$$
\boxed{
C_R(\phi(s))
=
C_I(\phi(s))
=
0.
}
\tag{27}
$$

that is,

$$
\boxed{
\sum_iX_i(s)^2=0
}
\tag{28}
$$

is preserved exactly.

This closes the preservation of zero closure along a continuous self-map parameter.

---

## 9. Exact Preservation of Zero Closure under Finite Iteration

The Euler update

$$
\phi+\eta F_T
$$

is correct in the tangent direction but, for finite $\eta$, leaves the manifold by $O(\eta^2)$.

For finite iterations we therefore use a retraction.

Let the tentative tangential update be

$$
\widetilde\phi
=
\phi^{(n)}
+
\eta
P_{\phi^{(n)}}F(\phi^{(n)})
\tag{29}
$$



Next, using the two normal directions $g_R,g_I$, set

$$
\phi^{(n+1)}
=
\widetilde\phi
+
\lambda_Rg_R(\widetilde\phi)
+
\lambda_Ig_I(\widetilde\phi)
\tag{30}
$$



Determine $\lambda_R,\lambda_I$ from the two equations

$$
C_R(\phi^{(n+1)})=0,
$$

$$
C_I(\phi^{(n+1)})=0
\tag{31}
$$



At regular points, by the implicit function theorem, for sufficiently small $\eta$ there exists a locally unique

$$
(\lambda_R,\lambda_I)
$$



Writing this operation as

$$
R_\phi(v)
$$

we have

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}
\left(
\eta
P_{\phi^{(n)}}F(\phi^{(n)})
\right)
}
\tag{32}
$$



By construction,

$$
\boxed{
\phi^{(n)}\in\mathcal Z_N
\Longrightarrow
\phi^{(n+1)}\in\mathcal Z_N.
}
\tag{33}
$$

Hence

$$
\boxed{
\mathcal F_N:
\mathcal Z_N
\to
\mathcal Z_N
}
\tag{34}
$$

holds exactly at every finite iteration.

This is the discrete dynamics that preserves zero closure.

The importance of this formula does not lie merely in a constraint-handling device that stabilises numerical computation. Since self-consistency is the starting condition of this axiom system, a state update is not allowed to leave the admissible set. Equations (32)–(34) show that **a discrete self-map that actually updates the phases without destroying zero closure can be constructed.** That is, one need not temporarily discard zero closure in order to introduce dynamics.

$$
\boxed{
\text{admissible state}
\xrightarrow{\;\text{update of }\phi\;}
\text{admissible state}
}
$$

holds at finite $N$. For a system that places self-consistency among its foundational axioms this is an indispensable verification, and it cannot be replaced by implicit state rewriting.

---

## 10. Self-Consistency Condition

A self-consistent fixed point is

$$
\boxed{
\phi_*
=
\mathcal F_N(\phi_*)
\quad
(\bmod 2\pi)
}
\tag{35}
$$



Since $\phi_*$ already lies on $\mathcal Z_N$, it satisfies zero closure.

At the fixed point the tangential update vanishes, so

$$
\boxed{
P_{\phi_*}F(\phi_*)=0.
}
\tag{36}
$$

that is,

$$
F(\phi_*)
\in
\operatorname{span}
\{g_R(\phi_*),g_I(\phi_*)\}.
\tag{37}
$$

Hence the action is stationary under admissible variations that keep zero closure.

Using Lagrange multipliers, the constrained stationarity condition is

$$
\boxed{
\frac{\partial S_N}{\partial\phi_i}
+
\lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\lambda_I
\frac{\partial C_I}{\partial\phi_i}
=
0.
}
\tag{38}
$$

This is the finite-$N$ discrete Euler–Lagrange equation that satisfies self-consistency and zero closure simultaneously.

---

## 11. The Unconstrained Sector and Local Current Conservation

In the sector where the normal reaction force of zero closure does not contribute to the local current sector, i.e.

$$
\lambda_R=\lambda_I=0
$$

, equation (38) becomes

$$
\frac{\partial S_N}{\partial\phi_i}=0
$$



Hence

$$
\boxed{
\sum_{j\sim i}J_{ij}=0.
}
\tag{39}
$$

This is the discrete continuity equation.

In the general constrained sector,

$$
\sum_{j\sim i}J_{ij}
=
\lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\lambda_I
\frac{\partial C_I}{\partial\phi_i}.
\tag{40}
$$

Thus the zero-closure constraint force can be separated as an explicit source.

---

## 12. Readout Resolution $N$ and Lattice Spacing

Let $L$ be the normalised readout length.

At resolution $N$,

$$
\boxed{
h_N=\frac{L}{N}.
}
\tag{41}
$$

Hence

$$
\boxed{
N\to\infty
\quad\Longleftrightarrow\quad
h_N\to0.
}
\tag{42}
$$

$N$ is not the absolute number of existing elements but the resolution with which the same structure is read.

---

## 13. Continuum Limit of the Discrete Current Density

Define the edge current density in direction $\mu$ as

$$
\boxed{
J_N^\mu
\left(
q+\frac{h_N}{2}e_\mu
\right)
=
\frac{A^2}{h_N}
\sin
\left[
\phi(q+h_Ne_\mu)-\phi(q)
\right].
}
\tag{43}
$$



By Taylor expansion,

$$
\phi(q+h_Ne_\mu)
-
\phi(q)
=
h_N\partial_\mu\phi
+
\frac{h_N^2}{2}
\partial_\mu^2\phi
+
O(h_N^3).
\tag{44}
$$

Also,

$$
\sin u
=
u
-
\frac{u^3}{6}
+
O(u^5).
\tag{45}
$$

Hence

$$
J_N^\mu
=
A^2\partial_\mu\phi
+
O(h_N).
\tag{46}
$$

Thus

$$
\boxed{
N\to\infty
\quad\Longrightarrow\quad
J_\mu
=
A^2\partial_\mu\phi.
}
\tag{47}
$$

Choosing a readout metric $g_{\mu\nu}$,

$$
\boxed{
J^\mu
=
A^2g^{\mu\nu}\partial_\nu\phi.
}
\tag{48}
$$

---

## 14. From Discrete Divergence to the Continuum Continuity Equation

Pairing the positive-direction outflow with the negative-direction inflow at vertex $q$,

$$
\sum_\mu
\left[
J_N^\mu
\left(
q+\frac{h_N}{2}e_\mu
\right)
-
J_N^\mu
\left(
q-\frac{h_N}{2}e_\mu
\right)
\right]
=0.
\tag{49}
$$

Dividing by $h_N$,

$$
\boxed{
\sum_\mu
\frac{
J_N^\mu(q+\tfrac12h_Ne_\mu)
-
J_N^\mu(q-\tfrac12h_Ne_\mu)
}{
h_N
}
=
0.
}
\tag{50}
$$

The Taylor expansion of the central difference is

$$
\frac{
J_N^\mu(q+\tfrac12h_Ne_\mu)
-
J_N^\mu(q-\tfrac12h_Ne_\mu)
}{
h_N
}
=
\partial_\mu J^\mu
+
O(h_N^2).
\tag{51}
$$

Hence

$$
\boxed{
N\to\infty
\quad\Longrightarrow\quad
\partial_\mu J^\mu=0.
}
\tag{52}
$$

This is the continuum Noether-type conservation law.

---

## 15. From Discrete Action to Continuum Action

The discrete action is

$$
S_N
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i).
$$

For small phase differences,

$$
\cos u
=
1-\frac{u^2}{2}
+
O(u^4).
\tag{53}
$$

Dropping the constant term,

$$
S_N
\sim
\frac{A^2}{2}
\sum_{\langle ij\rangle}
(\phi_j-\phi_i)^2.
\tag{54}
$$

In direction $\mu$,

$$
\phi(q+h_Ne_\mu)-\phi(q)
=
h_N\partial_\mu\phi
+
O(h_N^2).
$$

Hence, taking the Riemann sum with the cell volume $h_N^d$ and the metric weight,

$$
\boxed{
S_N
\longrightarrow
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq.
}
\tag{55}
$$

The Lagrangian density is

$$
\boxed{
\mathcal L
=
\frac{A^2}{2}
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi.
}
\tag{56}
$$



---

## 16. Connection to the Standard Partial Differential Field Equations

The Euler–Lagrange equation is

$$
\frac{\partial\mathcal L}{\partial\phi}
-
\partial_\mu
\left(
\frac{\partial\mathcal L}
{\partial(\partial_\mu\phi)}
\right)
=
0.
\tag{57}
$$

Since $\mathcal L$ does not depend on $\phi$ itself,

$$
\frac{\partial\mathcal L}{\partial\phi}=0.
$$

Also,

$$
\frac{\partial\mathcal L}
{\partial(\partial_\mu\phi)}
=
A^2g^{\mu\nu}\partial_\nu\phi
=
J^\mu.
\tag{58}
$$

Hence

$$
\boxed{
\partial_\mu J^\mu
=
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi
\right)
=
0.
}
\tag{59}
$$

If $A$ and $g^{\mu\nu}$ are constant,

$$
A^2
g^{\mu\nu}
\partial_\mu\partial_\nu\phi
=
0.
$$

Hence

$$
\boxed{
\Box_g\phi=0.
}
\tag{60}
$$

For a Euclidean readout,

$$
\boxed{
\nabla^2\phi=0.
}
\tag{61}
$$

For a Lorentzian readout,

$$
\boxed{
\Box\phi
=
-\partial_t^2\phi
+
\nabla^2\phi
=
0.
}
\tag{62}
$$

Thus, in the $N\to\infty$ continuum approximation, the discrete relational phase system connects to the standard partial differential equation of a massless scalar / phase field.

Here the Lorentz signature was not placed by privileging time at the foundation. It appears from the sign of the real form chosen at readout.

---

## 17. Constrained Continuum Field Equation

If zero closure is to be kept explicitly in the continuum limit as well, Lagrange multipliers are added to the action.

$$
\boxed{
S_{\mathrm{c}}
=
S_{\mathrm{cont}}
+
\int
\left[
\Lambda_R C_R
+
\Lambda_I C_I
\right]
d^dq.
}
\tag{63}
$$

Varying,

$$
\boxed{
\partial_\mu
\left(
A^2g^{\mu\nu}\partial_\nu\phi_i
\right)
=
\Lambda_R
\frac{\partial C_R}{\partial\phi_i}
+
\Lambda_I
\frac{\partial C_I}{\partial\phi_i}.
}
\tag{64}
$$

that is,

$$
\boxed{
A^2\Box_g\phi_i
=
-2\Lambda_R\sin2\phi_i
+
2\Lambda_I\cos2\phi_i.
}
\tag{65}
$$

In the sector where the constraint reaction vanishes, this returns to

$$
\boxed{
\Box_g\phi_i=0
}
\tag{66}
$$



Hence the general form is a constrained wave / Laplace type equation, and the standard free field equation appears as its zero-constraint-reaction sector.

---


## 18. Derivation of Local Gauge Structure

In this section we localise the global $U(1)$ phase structure of the preceding sections and derive the edge connection, the face curvature, the covariant derivative, and the Maxwell / Yang–Mills type continuum limit.

### 18.1 Anonymity of the Local Phase Origin

Under the global phase transformation

$$
\phi_i\mapsto\phi_i+\alpha
$$

the relational phase difference

$$
\phi_j-\phi_i
$$

is invariant.

If, however, we allow the phase origin of each vertex to be re-chosen independently,

$$
\boxed{
\phi_i\mapsto\phi_i+\alpha_i
}
\tag{67}
$$

then

$$
\phi_j-\phi_i
\mapsto
\phi_j-\phi_i+\alpha_j-\alpha_i
$$

and the bare phase difference is no longer invariant.

Hence, to compare phases between different vertices, a compensating quantity $\theta_{ij}$ on the edge is required.

We define the locally observable relational phase as

$$
\boxed{
\Delta_{ij}^{(g)}
=
\phi_j-\phi_i-\theta_{ij}
}
\tag{68}
$$



For $\Delta_{ij}^{(g)}$ to remain invariant under the local phase transformation (67), we must have

$$
\boxed{
\theta_{ij}
\mapsto
\theta_{ij}
+
\alpha_j-\alpha_i
}
\tag{69}
$$



This is the discrete gauge transformation law of a $U(1)$ connection.

Hence the edge connection is not a variable added from outside; it becomes necessary as

$$
\boxed{
\text{anonymity of the local phase origin}
\Longrightarrow
\text{edge connection comparing different vertices}
}
\tag{70}
$$



### 18.2 link variable

Exponentiating the edge connection, set

$$
\boxed{
U_{ij}
=
e^{-i\theta_{ij}}
}
\tag{71}
$$



Under the local transformation

$$
X_i
\mapsto
e^{i\alpha_i}X_i
$$

, if

$$
U_{ij}
\mapsto
e^{i\alpha_i}
U_{ij}
e^{-i\alpha_j}
$$

then

$$
X_i^*U_{ij}X_j
$$

is invariant.

Indeed,

$$
X_i^*U_{ij}X_j
\mapsto
e^{-i\alpha_i}X_i^*
e^{i\alpha_i}U_{ij}e^{-i\alpha_j}
e^{i\alpha_j}X_j
=
X_i^*U_{ij}X_j.
$$

Hence the gauge-covariant discrete current is

$$
\boxed{
J_{ij}^{(g)}
=
\operatorname{Im}
\left(
X_i^*U_{ij}X_j
\right)
}
\tag{72}
$$



Substituting the equal amplitude

$$
X_i=Ae^{i\phi_i}
$$

gives

$$
\boxed{
J_{ij}^{(g)}
=
A^2
\sin
\left(
\phi_j-\phi_i-\theta_{ij}
\right)
}
\tag{73}
$$



### 18.3 Gauge-Covariant Discrete Action

Likewise from the real part we obtain

$$
\boxed{
S_N^{(g)}
=
-A^2
\sum_{\langle ij\rangle}
\cos
\left(
\phi_j-\phi_i-\theta_{ij}
\right)
}
\tag{74}
$$



By (68), (69) the argument on each edge is gauge invariant, so

$$
\boxed{
S_N^{(g)}
\text{ is exactly invariant under local }U(1)\text{ transformations}
}
\tag{75}
$$



Varying with respect to $\phi_i$,

$$
\frac{\partial S_N^{(g)}}{\partial\phi_i}
=
-\sum_{j\sim i}
J_{ij}^{(g)}.
$$

so the stationarity condition is

$$
\boxed{
\sum_{j\sim i}
J_{ij}^{(g)}
=
0
}
\tag{76}
$$



This is the discrete form of the local gauge-covariant continuity equation.

### 18.4 Curvature Emerges from Simplex Faces

Consider the triangle $i\to j\to k\to i$.

Let the sum of the connection along the closed path be

$$
\boxed{
\Theta_{ijk}
=
\theta_{ij}
+
\theta_{jk}
+
\theta_{ki}
}
\tag{77}
$$



Under the local transformation,

$$
\theta_{ij}
\mapsto
\theta_{ij}+\alpha_j-\alpha_i,
$$

$$
\theta_{jk}
\mapsto
\theta_{jk}+\alpha_k-\alpha_j,
$$

$$
\theta_{ki}
\mapsto
\theta_{ki}+\alpha_i-\alpha_k.
$$

Hence

$$
\Theta_{ijk}
\mapsto
\Theta_{ijk}
+
(\alpha_j-\alpha_i)
+
(\alpha_k-\alpha_j)
+
(\alpha_i-\alpha_k)
=
\Theta_{ijk}.
$$

Therefore

$$
\boxed{
\Theta_{ijk}
\text{ is gauge invariant}
}
\tag{78}
$$



Hence the simplex hierarchy corresponds as

$$
\boxed{
\text{vertex}
\to
\text{phase field}
}
$$

$$
\boxed{
\text{edge}
\to
\text{connection}
}
$$

$$
\boxed{
\text{face}
\to
\text{curvature}
}
\tag{79}
$$



### 18.5 plaquette action

Take the minimal $2\pi$-periodic action for the face curvature as

$$
\boxed{
S_F
=
\beta
\sum_p
\left(
1-\cos\Theta_p
\right)
}
\tag{80}
$$



Since $\Theta_p$ is gauge invariant, so is $S_F$.

For small plaquettes,

$$
1-\cos\Theta_p
=
\frac12\Theta_p^2
+
O(\Theta_p^4).
\tag{81}
$$

Hence in the continuum limit it goes over to a curvature-squared action.

### 18.6 The Covariant Derivative Emerges in the Continuum Limit

Let the readout spacing be $h_N$.

For an edge in direction $\mu$, set

$$
\theta_{ij}
=
gh_NA_\mu(q)
+
O(h_N^2)
\tag{82}
$$



Also,

$$
\phi(q+h_Ne_\mu)-\phi(q)
=
h_N\partial_\mu\phi
+
O(h_N^2).
$$

Hence

$$
\frac{
\phi_j-\phi_i-\theta_{ij}
}{
h_N
}
\longrightarrow
\partial_\mu\phi-gA_\mu.
$$

Thus the gauge current is

$$
\boxed{
J^\mu_{(g)}
=
A^2
\left(
\partial^\mu\phi
-
gA^\mu
\right)
}
\tag{83}
$$



For the complex field

$$
\Psi=Ae^{i\phi}
$$

, setting

$$
D_\mu
=
\partial_\mu
-
igA_\mu
$$

we have

$$
D_\mu\Psi
=
iAe^{i\phi}
\left(
\partial_\mu\phi-gA_\mu
\right).
$$

Hence

$$
\boxed{
J^\mu_{(g)}
\propto
\operatorname{Im}
\left(
\Psi^*D^\mu\Psi
\right)
}
\tag{84}
$$

, in agreement with the standard covariant derivative.

### 18.7 Continuum Limit of the Face Curvature

Consider a small rectangular plaquette in the $\mu\nu$ plane.

The closed-path sum is

$$
\Theta_{\mu\nu}
=
\theta_\mu(q)
+
\theta_\nu(q+h_Ne_\mu)
-
\theta_\mu(q+h_Ne_\nu)
-
\theta_\nu(q).
\tag{85}
$$

Taylor-expanding (82),

$$
\Theta_{\mu\nu}
=
gh_N^2
\left(
\partial_\mu A_\nu
-
\partial_\nu A_\mu
\right)
+
O(h_N^3).
$$

Hence

$$
\boxed{
\frac{\Theta_{\mu\nu}}{gh_N^2}
\longrightarrow
F_{\mu\nu}
=
\partial_\mu A_\nu
-
\partial_\nu A_\mu
}
\tag{86}
$$



### 18.8 The Maxwell Action as a Limit

From (80) and (81),

$$
S_F
\sim
\frac{\beta}{2}
\sum_p
\Theta_p^2.
$$

Using (86),

$$
\Theta_p^2
=
g^2h_N^4
F_{\mu\nu}F^{\mu\nu}
+
O(h_N^5).
$$

Taking the Riemann sum including the cell volume and the coupling normalisation,

$$
\boxed{
S_F
\longrightarrow
-\frac14
\int
F_{\mu\nu}F^{\mu\nu}
\,d^dq
}
\tag{87}
$$



Hence

$$
\boxed{
\text{edge connection}
+
\text{face curvature}
\Longrightarrow
\text{Maxwell gauge dynamics}
}
\tag{88}
$$



### 18.9 Non-Abelian Generalisation

Let the vertex state be

$$
\Psi_i\in\mathbb C^r
$$

and the local internal basis change be

$$
\Psi_i
\mapsto
G_i\Psi_i,
\qquad
G_i\in SU(r)
\tag{89}
$$



To compare states at different vertices an edge transporter

$$
U_{ij}\in SU(r)
$$

is required, and it must transform as

$$
\boxed{
U_{ij}
\mapsto
G_iU_{ij}G_j^{-1}
}
\tag{90}
$$



The closed-path product

$$
\boxed{
W_{ijk}
=
U_{ij}U_{jk}U_{ki}
}
\tag{91}
$$

transforms as

$$
W_{ijk}
\mapsto
G_iW_{ijk}G_i^{-1}
$$

so

$$
\boxed{
\operatorname{Tr}W_{ijk}
\text{ is gauge invariant}
}
\tag{92}
$$



In the continuum limit, setting

$$
U_{ij}
=
\exp
\left(
igh_NA_\mu^aT^a
\right)
$$

the Baker–Campbell–Hausdorff expansion of the plaquette product yields

$$
\boxed{
F_{\mu\nu}
=
\partial_\mu A_\nu
-
\partial_\nu A_\mu
-
ig[A_\mu,A_\nu]
}
\tag{93}
$$



Thus the non-abelian Yang–Mills term

$$
[A_\mu,A_\nu]
$$

arises from the closed-path product of finite transporters.

### 18.10 What Has Been Derived in This Paper So Far

What this section obtained is

$$
\boxed{
\text{local phase anonymity}
\to
\text{edge connection}
\to
\text{face curvature}
\to
\text{covariant derivative}
\to
\text{Maxwell / Yang--Mills type continuum limit}
}
\tag{94}
$$



Hence the connection to the Standard Model is strengthened by one more step relative to the previous version.


## 19. Connection to the Symmetry Derivation of the Preceding Paper and Closure of the Standard-Model Representation

In this section we connect the symmetry structures obtained in the preceding paper to the local gauge connection, curvature, covariant derivative, and Yang–Mills type dynamics derived in this paper.

In the preceding paper, subtracting one complex constraint from the complex six-axis zero closure

$$
\sum_{n=1}^{6}X_n^2=0
$$

gave

$$
\boxed{
\dim_{\mathbb C}=5
}
\tag{95}
$$

independent complex degrees of freedom.

Reading these five degrees of freedom as a self-consistent Hermitian decomposition

$$
\boxed{
V=V_3\oplus V_2,
\qquad
\dim_{\mathbb C}V_3=3,
\qquad
\dim_{\mathbb C}V_2=2
}
\tag{96}
$$

the decomposition-preserving group is

$$
U(3)\times U(2)
$$



Removing the overall phase redundancy further gives

$$
\boxed{
S(U(3)\times U(2))
\cong
\frac{
SU(3)\times SU(2)\times U(1)
}{
\mathbb Z_6
}
}
\tag{97}
$$



The preceding paper treated this as a conditionally exact connection to the Standard-Model global gauge group.

In this paper, the anonymity of the local phase origin makes an edge connection necessary, the simplex face loop yields curvature, and in the continuum limit

$$
D_\mu
=
\partial_\mu-igA_\mu
$$

and

$$
F_{\mu\nu}
=
\partial_\mu A_\nu-\partial_\nu A_\mu-ig[A_\mu,A_\nu]
$$

appear, as derived in the previous section.

Hence, connecting the global stabiliser of the preceding paper with the local transporter of this paper, we obtain the bridge to local gauge dynamics

$$
\boxed{
S(U(3)\times U(2))
\quad\Longrightarrow\quad
\text{local }
\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6}
\text{ connection}
}
\tag{98}
$$



Below we check in turn how far the internal quantum numbers of one Standard-Model generation are fixed by this five-complex-degree-of-freedom structure.

### 19.1 The $U(1)$ Generator Ratio Is Fixed by the Trace-Zero Condition

On $V=V_3\oplus V_2$, the $U(1)$ generator commuting with $SU(3)$ and $SU(2)$ is block diagonal,

$$
Y
=
\begin{pmatrix}
y_3 I_3 & 0\\
0 & y_2 I_2
\end{pmatrix}
\tag{99}
$$



Since the Lie algebra condition of $S(U(3)\times U(2))$ is that the total trace vanish,

$$
\operatorname{Tr}Y
=
3y_3+2y_2
=
0.
\tag{100}
$$

Hence

$$
\boxed{
3y_3+2y_2=0
}
\tag{101}
$$

and the ratio is uniquely fixed as

$$
\boxed{
y_3:y_2=-2:3
}
\tag{102}
$$



The overall normalisation corresponds to the choice of unit of $U(1)$ charge. Following the Standard-Model convention,

$$
\boxed{
y_3=-\frac13,
\qquad
y_2=\frac12
}
\tag{103}
$$

gives

$$
\boxed{
Y
=
\operatorname{diag}
\left(
-\frac13,-\frac13,-\frac13,
\frac12,\frac12
\right)
}
\tag{104}
$$



The important point is that $-1/3$ and $1/2$ are not input independently. The ratio $-2:3$ is fixed by the $3+2$ decomposition and the trace-zero condition alone; the remaining common factor is the unit normalisation of charge.

Hence

$$
\boxed{
\text{hypercharge ratio}
\quad
3\left(-\frac13\right)
+
2\left(\frac12\right)
=0
}
\tag{105}
$$

follows from the structure of $S(U(3)\times U(2))$.

### 19.2 Five Degrees of Freedom and Simplex Two-Body Relations

Axiom A3 of the preceding paper treats $N$ vertices and all two-body relations as a simplex.

Reading the five complex degrees of freedom as a one-vertex register space

$$
V\simeq\mathbb C^5
$$

the oriented two-body relations are represented by the antisymmetric rank-two tensor

$$
\boxed{
\Lambda^2V
}
\tag{106}
$$



Its complex dimension is

$$
\boxed{
\dim_{\mathbb C}\Lambda^2V
=
\binom52
=
10
}
\tag{107}
$$



On the other hand, the dual representation on the one-body side is

$$
V^*
$$

with

$$
\dim_{\mathbb C}V^*=5.
$$

Hence

$$
\boxed{
V^*\oplus\Lambda^2V
}
\tag{108}
$$

has

$$
\boxed{
5+10=15
}
\tag{109}
$$

complex components.

This does not end as a mere coincidence of dimensions. Decomposing into $S(U(3)\times U(2))$ representations in the next section, it coincides with the left-handed Weyl representation of one Standard-Model generation.

### 19.3 Decomposition of $V^*$

By (96), (103),

$$
V
=
(\mathbf3,\mathbf1)_{-1/3}
\oplus
(\mathbf1,\mathbf2)_{1/2}.
\tag{110}
$$

The dual representation is

$$
\boxed{
V^*
=
(\overline{\mathbf3},\mathbf1)_{1/3}
\oplus
(\mathbf1,\mathbf2)_{-1/2}
}
\tag{111}
$$



In the left-handed Weyl notation of the Standard Model, these coincide with the internal quantum numbers of

$$
\boxed{
d^c
:
(\overline{\mathbf3},\mathbf1)_{1/3}
}
\tag{112}
$$

and

$$
\boxed{
L
:
(\mathbf1,\mathbf2)_{-1/2}
}
\tag{113}
$$



### 19.4 Decomposition of $\Lambda^2V$

The exterior square of a direct sum is

$$
\Lambda^2(V_3\oplus V_2)
=
\Lambda^2V_3
\oplus
(V_3\otimes V_2)
\oplus
\Lambda^2V_2.
\tag{114}
$$

We compute this term by term.

#### 19.4.1 $\Lambda^2V_3$

For the fundamental representation of $SU(3)$,

$$
\Lambda^2\mathbf3
\cong
\overline{\mathbf3}.
$$

Since hypercharge is additive,

$$
Y
=
-\frac13-\frac13
=
-\frac23.
$$

Hence

$$
\boxed{
\Lambda^2V_3
=
(\overline{\mathbf3},\mathbf1)_{-2/3}
}
\tag{115}
$$



This coincides with

$$
\boxed{
u^c
:
(\overline{\mathbf3},\mathbf1)_{-2/3}
}
\tag{116}
$$



#### 19.4.2 $V_3\otimes V_2$

The representation is

$$
(\mathbf3,\mathbf2)
$$

and the hypercharge is

$$
-\frac13+\frac12
=
\frac16.
$$

Hence

$$
\boxed{
V_3\otimes V_2
=
(\mathbf3,\mathbf2)_{1/6}
}
\tag{117}
$$



This coincides with

$$
\boxed{
Q
:
(\mathbf3,\mathbf2)_{1/6}
}
\tag{118}
$$



#### 19.4.3 $\Lambda^2V_2$

For the fundamental doublet of $SU(2)$,

$$
\Lambda^2\mathbf2
\cong
\mathbf1.
$$

The hypercharge is

$$
\frac12+\frac12=1.
$$

Hence

$$
\boxed{
\Lambda^2V_2
=
(\mathbf1,\mathbf1)_1
}
\tag{119}
$$



This coincides with

$$
\boxed{
e^c
:
(\mathbf1,\mathbf1)_1
}
\tag{120}
$$



### 19.5 The One-Generation Representation

Collecting (111), (115), (117), (119),

$$
\boxed{
V^*\oplus\Lambda^2V
=
(\overline3,1)_{1/3}
\oplus
(1,2)_{-1/2}
\oplus
(\overline3,1)_{-2/3}
\oplus
(3,2)_{1/6}
\oplus
(1,1)_1
}
\tag{121}
$$



That is,

$$
\boxed{
V^*\oplus\Lambda^2V
=
d^c\oplus L\oplus u^c\oplus Q\oplus e^c
}
\tag{122}
$$



This coincides with the 15 components of one Standard-Model generation, without a right-handed neutrino component, in the standard left-handed Weyl notation in which right-handed particles are represented by left-handed charge-conjugate fields.

Group-theoretically,

$$
\boxed{
V^*\oplus\Lambda^2V
\cong
\overline{\mathbf5}
\oplus
\mathbf{10}
}
\tag{123}
$$

is isomorphic to the well-known decomposition.

In this paper, however, the grand unified group $SU(5)$ is not placed at the starting point.

The order of derivation is

$$
\boxed{
\text{complex zero closure}
\to
5\text{ independent complex degrees of freedom}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
V^*\oplus\Lambda^2V
}
\tag{124}
$$



### 19.6 Direct Verification of Anomaly Cancellation

We compute the perturbative gauge anomalies directly for the representation content (121).

Everything below is counted in left-handed Weyl fields.

#### 19.6.1 $SU(3)^3$

Since $Q=(3,2)$ has two $SU(2)$ components, there are two $SU(3)$ fundamentals.

$u^c,d^c$ contribute one anti-fundamental each.

Taking the signs of the cubic anomaly coefficients of the fundamental and anti-fundamental as $+1,-1$ respectively,

$$
2-1-1=0.
$$

Hence

$$
\boxed{
\mathcal A_{SU(3)^3}=0
}
\tag{125}
$$



#### 19.6.2 $SU(3)^2U(1)_Y$

We use $T(\mathbf3)=T(\overline{\mathbf3})=1/2$.

$Q$ is an $SU(2)$ doublet and hence has multiplicity 2.

$$
\mathcal A_{SU(3)^2U(1)}
=
2\left(\frac16\right)\frac12
+
\left(-\frac23\right)\frac12
+
\left(\frac13\right)\frac12.
$$

Hence

$$
\mathcal A_{SU(3)^2U(1)}
=
\frac16-\frac13+\frac16
=
0.
$$

that is,

$$
\boxed{
\mathcal A_{SU(3)^2U(1)}=0
}
\tag{126}
$$



#### 19.6.3 $SU(2)^2U(1)_Y$

We use $T(\mathbf2)=1/2$.

Since $Q$ has colour multiplicity 3,

$$
\mathcal A_{SU(2)^2U(1)}
=
3\left(\frac16\right)\frac12
+
\left(-\frac12\right)\frac12.
$$

Hence

$$
\frac14-\frac14=0.
$$

Thus

$$
\boxed{
\mathcal A_{SU(2)^2U(1)}=0
}
\tag{127}
$$



#### 19.6.4 $U(1)_Y^3$

Counting all left-handed Weyl components with multiplicity,

$$
\mathcal A_{U(1)^3}
=
6\left(\frac16\right)^3
+
3\left(-\frac23\right)^3
+
3\left(\frac13\right)^3
+
2\left(-\frac12\right)^3
+
1^3.
\tag{128}
$$

The terms are

$$
\frac1{36},
\qquad
-\frac89,
\qquad
\frac19,
\qquad
-\frac14,
\qquad
1
$$

so that

$$
\frac1{36}
-\frac{32}{36}
+\frac4{36}
-\frac9{36}
+\frac{36}{36}
=
0.
$$

Hence

$$
\boxed{
\mathcal A_{U(1)^3}=0
}
\tag{129}
$$



#### 19.6.5 gravitational--$U(1)_Y$

The sum of hypercharges is

$$
6\left(\frac16\right)
+
3\left(-\frac23\right)
+
3\left(\frac13\right)
+
2\left(-\frac12\right)
+
1.
$$

Hence

$$
1-2+1-1+1=0.
$$

Therefore

$$
\boxed{
\mathcal A_{\mathrm{grav}^2U(1)}=0
}
\tag{130}
$$



#### 19.6.6 $SU(2)$ global anomaly

The total number of left-handed $SU(2)$ doublets, including colour multiplicity, is

$$
3\quad(Q)
+
1\quad(L)
=
4.
$$

which is even, so

$$
\boxed{
N_{\mathrm{doublet}}=4
}
\tag{131}
$$

and the $SU(2)$ global anomaly condition is also satisfied.

From the above,

$$
\boxed{
V^*\oplus\Lambda^2V
\text{: its one-generation representation satisfies Standard-Model anomaly cancellation}
}
\tag{132}
$$

is confirmed directly.

### 19.7 Hypercharge and Anomaly Cancellation Are Not Independent Additional Conditions

We summarise the results so far.

The hypercharge ratio follows from

$$
\boxed{
3y_3+2y_2=0
}
$$

, the trace-zero condition of $S(U(3)\times U(2))$.

The one-generation representation is constructed from the five degrees of freedom and the simplex two-body relations as

$$
\boxed{
V^*\oplus\Lambda^2V
}
$$



Anomaly cancellation for that representation holds by the direct computation (125)–(131).

Hence

$$
\boxed{
\text{hypercharge}
+
\text{anomaly cancellation}
}
\tag{133}
$$

need not be input independently to match the Standard Model.

The derivation chain

$$
\boxed{
5\text{ complex degrees of freedom}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
Y
\to
V^*\oplus\Lambda^2V
\to
\text{anomaly cancellation}
}
\tag{134}
$$

closes.

### 19.8 Lorentz Spinors and Chirality

In the preceding paper, for the Lorentz partial readout of $(x,y,z,t)$, the connection to

$$
SO^+(3,1)
$$

and its double cover

$$
\boxed{
Spin^+(3,1)
\cong
SL(2,\mathbb C)
}
\tag{135}
$$

was obtained.

The minimal Weyl spinor representations of $SL(2,\mathbb C)$ are the two

$$
\boxed{
\left(\frac12,0\right)
}
\tag{136}
$$

and

$$
\boxed{
\left(0,\frac12\right)
}
\tag{137}
$$

, which are complex conjugates of each other.

On the internal-representation side there are likewise

$$
\boxed{
V^*\oplus\Lambda^2V
}
\tag{138}
$$

and its complex conjugate

$$
\boxed{
V\oplus\Lambda^2V^*
}
\tag{139}
$$



Hence the complete candidate including the mirror can be written as the conjugate pair of sectors

$$
\boxed{
\left(\frac12,0\right)
\otimes
\left(
V^*\oplus\Lambda^2V
\right)
}
\tag{140}
$$

and

$$
\boxed{
\left(0,\frac12\right)
\otimes
\left(
V\oplus\Lambda^2V^*
\right)
}
\tag{141}
$$



Writing the Standard Model in left-handed Weyl fields only, the observed one generation is collected as the internal quantum numbers on the side of (140),

$$
d^c,\ L,\ u^c,\ Q,\ e^c
$$



Here, however, there is an important logical boundary.

The Fermi-type double cover / odd-harmonic sector of the preceding paper gives the existence of a fermionic sector, but by itself does not uniquely determine which of (140) and (141) is selected as the self-consistent physical sector.

Hence the last remaining problem for the closure of the Standard-Model representation is

$$
\boxed{
\text{self-consistent selection of the chirality sector}
}
\tag{142}
$$



This is not a problem of "generating chirality from nothing".

It is the problem of determining why one of the already existing conjugate sectors,

$$
\left(\frac12,0\right)
\otimes
(V^*\oplus\Lambda^2V)
$$

and

$$
\left(0,\frac12\right)
\otimes
(V\oplus\Lambda^2V^*)
$$

, is physically selected.

### 19.9 The Last Problem Written in the Language of the Present Axiom System

The preceding paper already contains

1. oriented simplices
2. $\partial^2=0$ and $d^2=0$
3. Lorentz spin double cover
4. the Fermi-type single/double cover classification
5. the self-consistent fixed points of A4 and stabiliser selection



Hence the last problem is not an abstract "origin of left–right asymmetry", but the problem of analysing which of the two sectors is selected as the stable fixed-point sector when

$$
\boxed{
\text{simplex orientation}
+
\text{Lorentz Weyl sector}
+
\text{Fermi-type double cover}
+
\text{A4 self-consistency}
}
\tag{143}
$$

acts on

$$
\boxed{
V^*\oplus\Lambda^2V
}
\quad\text{and}\quad
\boxed{
V\oplus\Lambda^2V^*
}
\tag{144}
$$



That is, it suffices to determine

$$
\boxed{
\mathcal F_N
\text{ stabilises which fixed point under}
\text{orientation reversal / complex conjugation}
}
\tag{145}
$$



Once this selection is fixed to one side,

$$
\boxed{
\text{chirality}
\to
\text{hypercharge}
\to
\text{one-generation representation}
\to
\text{anomaly cancellation}
}
$$

— the derivation chain of the Standard-Model internal representation — closes.

### 19.10 The Range Currently in Agreement with the Standard Model

We organise the present state of affairs combining this paper and the preceding one.

| Structure | Derivation route | Present status |
|---|---|---|
| Lorentz-type signature | complex imaginary axis $it$ | derived |
| $Spin^+(3,1)$ | double cover of the Lorentz group | connected |
| 5 complex degrees of freedom | one zero closure on $\mathbb C^6$ | derived |
| $3\oplus2$ readout | Hermitian decomposition of the 5 d.o.f. | derived readout; generalisation of the selection rule is a separate task |
| $S(U(3)\times U(2))$ | decomposition preservation + removal of overall phase | conditionally exact |
| local connection | anonymity of local phase / basis origin | derived in this paper |
| curvature | simplex face loop | derived in this paper |
| Maxwell / Yang–Mills type action | plaquette limit | connected in this paper |
| hypercharge ratio | $3y_3+2y_2=0$ | derived |
| 15 components of one generation | $V^*\oplus\Lambda^2V$ | derived |
| quantum numbers of $Q,u^c,d^c,L,e^c$ | $3\oplus2$ representation decomposition | derived |
| perturbative anomaly cancellation | direct sum of representations | verified |
| $SU(2)$ global anomaly | 4 doublets | verified |
| Bose/Fermi/mixed sectors | odd/even harmonics, double cover | derived and numerically confirmed in our own papers |
| chirality sector selection | orientation + Weyl + A4 | **last unclosed point** |

Hence, for the internal gauge representation of the Standard Model,

$$
\boxed{
\text{what remains last is the self-consistent selection of the chirality sector}
}
\tag{146}
$$



This does not mean, however, that all phenomena of the Standard Model have been derived.

The Higgs radial mode, Yukawa couplings, the number of generations, the concrete mass hierarchy, quantum corrections and renormalisation are dynamical and quantisation problems separate from the closure of the internal gauge representation.

What this section closed is

$$
\boxed{
\text{Standard-Model gauge group}
+
\text{local gauge geometry}
+
\text{one-generation internal representation}
+
\text{hypercharge}
+
\text{anomaly cancellation}
}
\tag{147}
$$

and no more.


## 20. Theoretical Closure of the Chirality Sector

Up to the previous section, for the internal gauge representation we obtained

$$
5\text{ complex degrees of freedom}
\to
3\oplus2
\to
S(U(3)\times U(2))
\to
Y
\to
V^*\oplus\Lambda^2V
$$

. What remained was how one of the two conjugate Weyl sectors is selected self-consistently.

Here we identify the A/B two-channel selection system already constructed in our earlier papers with the conjugate two Weyl sectors, and rewrite its selection rule as a chirality order parameter.

### 20.1 The Conjugate Two Weyl Sectors

Let the minimal Weyl representations on the Lorentz side be

$$
\left(\frac12,0\right),
\qquad
\left(0,\frac12\right)
$$



On the internal-representation side there are

$$
V^*\oplus\Lambda^2V
$$

and its complex conjugate

$$
V\oplus\Lambda^2V^*
$$



We therefore define the conjugate pair of sectors as

$$
\boxed{
\mathcal H_L
=
\left(\frac12,0\right)
\otimes
\left(
V^*\oplus\Lambda^2V
\right)
}
\tag{148}
$$

$$
\boxed{
\mathcal H_R
=
\left(0,\frac12\right)
\otimes
\left(
V\oplus\Lambda^2V^*
\right)
}
\tag{149}
$$



Writing complex conjugation / orientation reversal as $\mathcal P_\chi$,

$$
\boxed{
\mathcal P_\chi:
\mathcal H_L
\leftrightarrow
\mathcal H_R
}
\tag{150}
$$



### 20.2 Identification with the Earlier A/B Two-Channel System

In the earlier papers, from two complex amplitudes $a,b$ we defined

$$
p_A=|a|^2,
\qquad
p_B=|b|^2
$$

and used

$$
\boxed{
S
=
\frac{p_A-p_B}{p_A+p_B}
}
\tag{151}
$$

as the order parameter of two-state selection.

Under $A\leftrightarrow B$,

$$
\boxed{
S\mapsto-S
}
\tag{152}
$$



Here we identify

$$
\boxed{
A\equiv\mathcal H_L,
\qquad
B\equiv\mathcal H_R
}
\tag{153}
$$



Then

$$
p_L:=p_A,
\qquad
p_R:=p_B
$$

and

$$
\boxed{
S_\chi
:=
\frac{p_L-p_R}{p_L+p_R}
}
\tag{154}
$$

is the chirality order parameter.

Under the mirror transformation,

$$
\boxed{
\mathcal P_\chi:
S_\chi\mapsto-S_\chi
}
\tag{155}
$$



Hence the earlier A/B selection problem is mathematically isomorphic to the selection problem of the conjugate two Weyl sectors.

### 20.3 The Internal Correlation of the Earlier Papers Is Mirror-Odd

In the earlier papers, the sum of cross terms was

$$
C
=
\sum_{m<n}x_mx_n
$$

and the quantity reading the phase difference arising between A and B was defined as

$$
\boxed{
J
=
\operatorname{Im}
\left(
B^{*2}C
\right)
}
\tag{156}
$$



Conjugating all complex quantities,

$$
B^{*2}C
\mapsto
\left(B^{*2}C\right)^*.
$$

Hence

$$
\boxed{
J\mapsto-J
}
\tag{157}
$$



Also,

$$
\Delta\phi_{AB}
=
\arg
\left(
1+\frac{2C}{B^2}
\right)
$$

goes, under complex conjugation, to

$$
\boxed{
\Delta\phi_{AB}
\mapsto
-\Delta\phi_{AB}
}
\tag{158}
$$



Hence

$$
J,
\qquad
\Delta\phi_{AB},
\qquad
S_\chi
$$

are all odd under the same mirror transformation.

That is,

$$
\boxed{
\mathcal P_\chi:
(J,\Delta\phi,S_\chi)
\mapsto
(-J,-\Delta\phi,-S_\chi)
}
\tag{159}
$$



### 20.4 The Minimal Normal Form of Chirality Selection

Let the evolution of $S_\chi$ be

$$
\dot S_\chi
=
F(S_\chi,J)
$$



Since the fundamental equations are mirror symmetric,

$$
\boxed{
F(-S_\chi,-J)
=
-F(S_\chi,J)
}
\tag{160}
$$

must hold.

Expanding analytically near the symmetric point

$$
S_\chi=0,
\qquad
J=0
$$

, mirror-even terms such as a constant, $J^2$, or $S_\chi^2$ are not allowed.

At lowest order,

$$
F(S_\chi,J)
=
\lambda J
+
aS_\chi
+
bS_\chi^3
+
O(J^3,J^2S_\chi,JS_\chi^2,S_\chi^5).
\tag{161}
$$

The nonlinear selection term already used in the earlier papers is

$$
gS(1-S^2)
=
gS-gS^3
$$



Hence the minimal normal form is

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi
\left(
1-S_\chi^2
\right)
}
\tag{162}
$$



This is not an equation to which a new arbitrary selection rule has been added.

- $\lambda J$ is the lowest-order coupling by which the mirror-odd internal phase correlation already defined in the earlier papers maps into the selection variable;
- $gS(1-S^2)$ is the nonlinear selection term already used in the earlier papers.



### 20.5 Spontaneous Chirality Selection at $J=0$

First let

$$
J=0
$$



Then

$$
\dot S_\chi
=
gS_\chi(1-S_\chi^2).
\tag{163}
$$

The fixed points are

$$
\boxed{
S_\chi^*=0,\quad +1,\quad -1
}
\tag{164}
$$



Writing the right-hand side as

$$
f(S)=gS(1-S^2)
$$

we have

$$
f'(S)
=
g(1-3S^2).
$$

Hence

$$
f'(0)=g,
$$

$$
f'(\pm1)=-2g.
$$

For $g>0$,

$$
\boxed{
S_\chi=0
\text{ is unstable}
}
\tag{165}
$$

and

$$
\boxed{
S_\chi=\pm1
\text{ are stable}
}
\tag{166}
$$



Hence, with the fundamental equations left-right symmetric, the spontaneous chirality selection

$$
\boxed{
S_\chi=0
\longrightarrow
S_\chi=+1
\quad\text{or}\quad
S_\chi=-1
}
\tag{167}
$$

occurs.

### 20.6 The Internal Phase Correlation $J$ Determines the Selected Sign

When $J$ is small but nonzero,

$$
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
$$



At the symmetric point $S_\chi=0$,

$$
\left.
\dot S_\chi
\right|_{S_\chi=0}
=
\lambda J.
\tag{168}
$$

Hence, if

$$
\lambda J>0
$$

then $S_\chi$ starts to move in the positive direction, and if

$$
\lambda J<0
$$

then in the negative direction.

Thereafter the nonlinear term with $g>0$ amplifies the difference, giving

$$
\boxed{
\lambda J>0
\Longrightarrow
S_\chi\to+1
}
\tag{169}
$$

$$
\boxed{
\lambda J<0
\Longrightarrow
S_\chi\to-1
}
\tag{170}
$$



Under the mirror transformation,

$$
J\to-J,
\qquad
S_\chi\to-S_\chi
$$

so the equation as a whole is invariant.

Hence no chirality is favoured as a foundational axiom.

The selected sign is determined by the sign of the self-consistent internal phase correlation $J$.

### 20.7 Agreement with the Causal Chain of the Earlier Papers

The earlier papers examined the causal chain in which integer harmonics and phase closure re-close the cross-term correlation over long times, so that

$$
J\neq0
$$

arises,

$$
\Delta\phi_{AB}\neq0
$$

appears,

$$
S\neq0
$$

is reached, and nonlinear selection amplifies it into one state.

Rewritten in the language of chirality, this becomes

$$
\boxed{
\text{Fermi-type integer harmonics / phase closure}
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{cross correlation }C
}
$$

$$
\Downarrow
$$

$$
\boxed{
J
=
\operatorname{Im}(B^{*2}C)
\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
\Delta\phi_{LR}\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
S_\chi\neq0
}
$$

$$
\Downarrow
$$

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
}
$$

$$
\Downarrow
$$

$$
\boxed{
S_\chi\to\pm1
}
\tag{171}
$$



Hence the dynamics of chirality selection need not be newly added; it is obtained by identifying the earlier Fermi-type nonlinear selection mechanism with the conjugate two Weyl sectors.

### 20.8 Completion of the Theoretical Connection to the Standard-Model Internal Representation

Up to the previous section, the hypercharge ratio was fixed from

$$
S(U(3)\times U(2))
$$



$$
V^*\oplus\Lambda^2V
$$

gave the 15 Weyl components of one generation, and anomaly cancellation held.

In this section we further connected the earlier A/B nonlinear selection mechanism to the mirror pair

$$
\mathcal H_L
\leftrightarrow
\mathcal H_R
$$



Hence, for the **Standard-Model internal gauge representation and chirality selection** treated in this paper, the theoretical derivation chain closes as

$$
\boxed{
\text{zero closure}
}
$$

$$
\Downarrow
$$

$$
\boxed{
5\text{ complex degrees of freedom}
}
$$

$$
\Downarrow
$$

$$
\boxed{
3\oplus2
}
$$

$$
\Downarrow
$$

$$
\boxed{
S(U(3)\times U(2))
}
$$

$$
\Downarrow
$$

$$
\boxed{
SU(3)\times SU(2)\times U(1)
\text{ local gauge geometry}
}
$$

$$
\Downarrow
$$

$$
\boxed{
Y
}
$$

$$
\Downarrow
$$

$$
\boxed{
V^*\oplus\Lambda^2V
=
d^c\oplus L\oplus u^c\oplus Q\oplus e^c
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{anomaly cancellation}
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{mirror-odd }J
}
$$

$$
\Downarrow
$$

$$
\boxed{
\text{chirality selection }S_\chi\to\pm1
}
\tag{172}
$$



In this sense, **the theoretical connection to the Standard-Model gauge group, local gauge geometry, one-generation internal representation, hypercharge, anomaly cancellation, and chirality selection is complete**.

"Complete" here refers to the gauge-representation/chirality chain that is the subject of this paper.

The Higgs radial mode, Yukawa couplings, the origin of three generations, the concrete mass hierarchy, and quantum corrections / renormalisation are separate problems and are not included in the range of theoretical closure of this paper.

---

## 21. Numerical Verification Specification — Verifying Chirality Selection with Minimal Changes to the Existing Experimental System

After the theoretical derivation has closed, what remains is the numerical verification of whether the existing Fermi-type A/B experimental system reproduces the theoretical equation (162) under the identification with the conjugate two Weyl sectors.

No new fundamental dynamics is added.

The existing update and scattering rules are left unchanged; the work consists mainly of **redefining readout quantities and adding measurements**.

### 21.1 Reinterpretation of the Existing A/B Variables

The existing code's

$$
A,\qquad B
$$

are read as

$$
\boxed{
A\equiv\mathcal H_L,
\qquad
B\equiv\mathcal H_R
}
\tag{173}
$$



The existing

$$
p_A=|a|^2,
\qquad
p_B=|b|^2
$$

are set to

$$
p_L:=p_A,
\qquad
p_R:=p_B
$$



The existing selection order parameter

$$
\boxed{
S_\chi
=
\frac{p_L-p_R}{p_L+p_R}
}
\tag{174}
$$

is used as is.

### 21.2 Observables to Be Added

Leaving the existing equations of motion unchanged, record the following at each iteration:

$$
C
=
\sum_{m<n}x_mx_n
\tag{175}
$$

$$
\boxed{
J_\chi
=
\operatorname{Im}
\left(
B^{*2}C
\right)
}
\tag{176}
$$

$$
\boxed{
\Delta\phi_\chi
=
\arg
\left(
1+\frac{2C}{B^2}
\right)
}
\tag{177}
$$

$$
S_\chi
$$

and

$$
\boxed{
G_\chi
=
gS_\chi
\left(
1-S_\chi^2
\right)
}
\tag{178}
$$



Further, from finite differences, record

$$
\boxed{
\dot S_{\chi,\mathrm{num}}
=
\frac{
S_\chi(n+1)-S_\chi(n)
}{
\Delta s
}
}
\tag{179}
$$



### 21.3 Direct Verification of the Minimal Normal Form

Fit the theoretical equation

$$
\dot S_\chi
=
\lambda J_\chi
+
gS_\chi(1-S_\chi^2)
$$

to the numerical time series.

That is, with

$$
\boxed{
\dot S_{\chi,\mathrm{num}}
=
\lambda_{\mathrm{fit}}J_\chi
+
g_{\mathrm{fit}}
S_\chi(1-S_\chi^2)
+
\epsilon
}
\tag{180}
$$

, estimate

$$
\lambda_{\mathrm{fit}},
\qquad
g_{\mathrm{fit}}
$$



The condition for theoretical support is

$$
\boxed{
\lambda_{\mathrm{fit}}\neq0,
\qquad
g_{\mathrm{fit}}>0
}
\tag{181}
$$



Further confirm that the residual $\epsilon$ is sufficiently small compared with the two main terms.

### 21.4 mirror test

Prepare a mirror run in which the initial state and all internal phases are complex conjugated.

Let the original run be

$$
\mathcal R
$$

and the mirror run be

$$
\mathcal R^*
$$



Theoretically,

$$
\boxed{
J_\chi^{*}(n)
=
-J_\chi(n)
}
\tag{182}
$$

$$
\boxed{
\Delta\phi_\chi^{*}(n)
=
-\Delta\phi_\chi(n)
}
\tag{183}
$$

$$
\boxed{
S_\chi^{*}(n)
=
-S_\chi(n)
}
\tag{184}
$$

must hold.

On the other hand, mirror-even quantities such as

$$
Q=p_L+p_R,
$$

, the total norm

$$
|C|,
$$

, the localisation degree $L$,

$$
N_{\mathrm{eff}}
$$

etc. must coincide.

Hence the pass condition of the mirror test is

$$
\boxed{
\max_n
|S_\chi^*(n)+S_\chi(n)|
<
\varepsilon_{\mathrm{mirror}}
}
\tag{185}
$$

together with the analogous conditions for $J_\chi,\Delta\phi_\chi$.

### 21.5 Symmetric Initial Condition Test

In theory,

$$
S_\chi=0,
\qquad
J_\chi=0
$$

is a symmetric fixed point.

Hence, for mechanically perfectly symmetric initial conditions, neither side may be chosen arbitrarily apart from floating-point error.

Next, use small seeds differing only in sign,

$$
J_\chi(0)=+\epsilon
$$

and

$$
J_\chi(0)=-\epsilon
$$



The theoretical prediction is

$$
\boxed{
+\epsilon
\Longrightarrow
S_\chi\to+1
}
\tag{186}
$$

$$
\boxed{
-\epsilon
\Longrightarrow
S_\chi\to-1
}
\tag{187}
$$



Lower the absolute value of the seed to about

$$
10^{-k},
\qquad
k=2,\ldots,12
$$

and measure the region in which the sign alone determines the final chirality.

### 21.6 Phase-Closure Breaking Test

In the causal chain of the earlier papers,

$$
\text{integer harmonics / phase closure}
\to
J_\chi
$$

is required.

Hence compare the following four conditions with the same amplitude distribution:

1. integer harmonics + phase closure
2. integer harmonics + phase randomisation
3. non-integer frequencies + phase closure
4. non-integer frequencies + phase randomisation

The condition for theoretical support is that only condition 1 shows a long-lasting

$$
|J_\chi|>0
$$

together with chirality selection.

### 21.7 Causal Order Test

Let the cross-correlations be

$$
R_{J\Delta\phi}(\ell)
=
\operatorname{corr}
\left(
J_\chi(n),
\Delta\phi_\chi(n+\ell)
\right)
$$

$$
R_{\Delta\phi S}(\ell)
=
\operatorname{corr}
\left(
\Delta\phi_\chi(n),
S_\chi(n+\ell)
\right)
$$



The theoretically required order is

$$
\boxed{
J_\chi
\to
\Delta\phi_\chi
\to
S_\chi
\to
|S_\chi|\approx1
}
\tag{188}
$$



Confirm that the correlation maxima appear at positive lag.

### 21.8 Separating the Role of $g$

As the earlier material already emphasised,

$$
\boxed{
J_\chi\text{ generation}
\neq
S_\chi\text{ amplification}
}
\tag{189}
$$



Hence compare

$$
g<0,\quad g=0,\quad g>0
$$



- $g=0$: is a signed small $S_\chi$ generated by the internal correlation alone?
- $g>0$: is that seed amplified to $\pm1$?
- $g<0$: is $S_\chi$ restored to 0 while $J_\chi$ itself remains?

If this separation holds, the theoretical interpretation

$$
\boxed{
J_\chi
=
\text{chirality seed},
\qquad
gS_\chi(1-S_\chi^2)
=
\text{selection amplifier}
}
\tag{190}
$$

is verified directly.

### 21.9 Bose-Type Control

The earlier papers used as a control the fact that one-state selection does not hold for Bose-type linear maps.

The same comparison is kept for the chirality readout.

In the Fermi-type sector,

$$
S_\chi\to\pm1
$$

occurs, whereas in the Bose-type control one confirms that

$$
S_\chi\approx0
$$

or periodic mixing is all that remains.

This verifies that chirality selection is not a generic phenomenon of two-channel exchange as such, but is tied to Fermi-type internal closure and nonlinear response.

### 21.10 Axiom-Preservation Audit — A Necessary Condition for Adoption as Dynamics

The most important point of this verification is not only that the selection of $S_\chi$ is observable. **One simultaneously confirms that the phase update itself does not destroy the foundational axioms.** In this system, where self-consistency is a foundational condition, an update rule that fails to preserve the axioms is not adopted as dynamics however plausible its output.

At minimum, for each iteration $n$, record

$$
\epsilon_C^{(n)}
=
\left|\sum_i X_i^{(n)2}\right|
$$

and confirm that the zero-closure error stays within numerical precision. Further, for finite recurrence $U^N=I$, simplex closure, norm, current conservation, and mirror symmetry, compare the audit quantities already defined in the existing code before and after the update.

The decision principle is clear:

$$
\boxed{
\text{observed selection holds}
\;\land\;
\text{axiom preservation holds}
}
$$

only this case counts as theoretical support. A run in which selection occurs but zero closure or any other admissibility condition is destroyed is rejected. Conversely, if $J_\chi\to\Delta\phi_\chi\to S_\chi$ is reproduced while the axioms are preserved, it is a numerical confirmation that **the discrete phase update functions as a dynamics inside the admissible state space**.

For each run, record

$$
\left|
\sum_iX_i^2
\right|
$$



Further, record

$$
p_L+p_R,
$$

, the total norm

$$
\operatorname{div}_dJ,
$$

, and the conservation error of the gauge current derived in this paper.

Even if chirality selection occurs,

$$
\boxed{
\text{zero closure}
}
$$

$$
\boxed{
\text{global conservation laws}
}
$$

must not be violated.

### 21.11 Policy for Changes to the Existing Code

The principle is not to change the existing dynamics.

Changes are limited to the following three kinds.

1. **Reinterpretation of names**
   - `A -> L_sector`
   - `B -> R_sector`
   - `S -> S_chi`

2. **Addition of observables**
   - `C_cross`
   - `J_chi`
   - `delta_phi_chi`
   - `dS_chi`
   - `normal_form_residual`

3. **Addition of mirror runs**
   - generate the conjugate initial condition of the full complex state
   - apply the identical dynamics
   - compare the transformation laws of odd/even observables

Hence no new interaction term is added to the code in order to verify the theory.

### 21.12 Pass Criteria

The necessary conditions for numerical support of the theoretical connection are the following.

1. In the mirror run,
   $$
   J_\chi\to-J_\chi,\quad
   \Delta\phi_\chi\to-\Delta\phi_\chi,\quad
   S_\chi\to-S_\chi
   $$
   holds.

2. Mirror-even conserved quantities coincide.

3. A persistent $J_\chi$ arises only under the integer-harmonic + phase-closure condition.

4. $J_\chi$ precedes or coincides with $\Delta\phi_\chi$.

5. $\Delta\phi_\chi$ precedes the growth of $S_\chi$.

6. The normal-form fit gives
   $$
   \lambda_{\mathrm{fit}}\neq0,
   \qquad
   g_{\mathrm{fit}}>0
   $$


7. Reversing the sign of the seed reverses the sign of the final $S_\chi$.

8. The same stable chirality selection does not hold in the Bose-type control.

9. The zero-closure, norm, and current conservation errors stay within tolerance.

If all of these hold,

$$
\boxed{
\text{earlier Fermi-type A/B selection system}
=
\text{numerical realisation of Weyl chirality selection}
}
\tag{191}
$$

is supported.

### 21.13 Connection to the Subsequent Paper

The theoretical derivation has closed in this paper.

The subject of the subsequent paper is not to add a new theory but

$$
\boxed{
\text{verification of the derivation of chirality selection on the existing numerical experimental system}
}
\tag{192}
$$



The subsequent paper will carry out

- mirror symmetry test,
- the causal order $J_\chi\to\Delta\phi_\chi\to S_\chi$,
- normal form fit,
- phase-closure breaking,
- the Bose/Fermi control,
- the zero-closure preservation audit,

and verify whether the theoretical connection of this paper is reproduced on the existing dynamics.

## 22. From Discrete Dynamics to Continuum Dynamics

### 22.1 The Self-Map Parameter Is Not Physical Time

The parameter $s$ of this section is not a substitute for physical time. It is a **construction/selection parameter** for reaching a self-consistent solution without destroying zero closure. Physical time $t$ was treated in the preceding paper as an observation direction that appears, by Lorentz readout, from a complex axis of the same kind as the other axes. Hence this paper is not a construction that "identifies the relaxational flow along $s$ with physical time evolution".

This distinction does not weaken the claim about dynamics; rather, it is the core of this axiom system. What is sought is not a rule that rewrites states unconditionally along an external time, but **a map that selects a self-consistent solution preserving the axioms over the whole configuration including the time axis**. When that fixed point is expressed in the Lorentz readout, it is read as a partial differential field equation containing physical time.

### 22.2 The Infinitesimal Iteration Limit

We take the infinitesimal iteration limit of the discrete self-map.

Let the iteration interval be $\Delta s$, with

$$
\eta=\Delta s.
$$

The tangential part of (32) is

$$
\frac{
\phi^{(n+1)}-\phi^{(n)}
}{
\Delta s
}
=
P_\phi F
+
O(\Delta s).
$$

Hence

$$
\boxed{
\Delta s\to0
\quad\Longrightarrow\quad
\frac{d\phi}{ds}
=
P_\phi F.
}
\tag{67}
$$

Taking also the spatial readout to $N\to\infty$,

$$
F_i
=
\sum_{j\sim i}
\sin(\phi_j-\phi_i)
$$

goes over to a second difference.

On a regular lattice,

$$
F_N(q)
=
\sum_\mu
\left[
\sin(\phi(q+h_Ne_\mu)-\phi(q))
+
\sin(\phi(q-h_Ne_\mu)-\phi(q))
\right].
\tag{68}
$$

By Taylor expansion,

$$
\boxed{
F_N(q)
=
h_N^2
\sum_\mu
\partial_\mu^2\phi
+
O(h_N^4).
}
\tag{69}
$$

Hence, renormalising the iteration parameter not as

$$
\tau
=
s\,h_N^2
$$

but, to keep a nontrivial limit, as

$$
d\tau
=
h_N^2\,ds
$$

, we obtain

$$
\boxed{
\frac{\partial\phi}{\partial\tau}
=
P_\phi
\left(
\sum_\mu
\partial_\mu^2\phi
\right).
}
\tag{70}
$$

Writing the constraint reaction with Lagrange multipliers,

$$
\boxed{
\frac{\partial\phi_i}{\partial\tau}
=
\Box_g\phi_i
-
\lambda_R\frac{\partial C_R}{\partial\phi_i}
-
\lambda_I\frac{\partial C_I}{\partial\phi_i}.
}
\tag{71}
$$

This is a continuum relaxational dynamics toward the self-consistent solution.

On the other hand, if one chooses the physical Lorentz readout and reads it as the Euler–Lagrange equation of the action, the hyperbolic PDE (65) becomes the physical field equation.

Hence

$$
\boxed{
\text{iteration parameter }s
\neq
\text{physical time }t
}
\tag{72}
$$



By this distinction, the discrete self-map and the standard continuum field equation coexist without privileging time at the foundation. Accordingly, the relaxational dynamics along $s$ must not be interpreted as physical time evolution itself. $s$ is the map parameter that constructs the admissible self-consistent configuration, and physical time $t$ is contained inside the Lorentz readout of the fixed-point configuration.

---

## 23. The Complete Form of the Derivation Chain

The entire derivation of this paper is

$$
\boxed{
\sum_iX_i^2=0
}
$$

$$
\Downarrow
$$

$$
\boxed{
X_i=Ae^{i\phi_i},
\qquad
\sum_i e^{2i\phi_i}=0
}
$$

$$
\Downarrow
$$

$$
\boxed{
J_{ij}
=
A^2\sin(\phi_j-\phi_i)
}
$$

$$
\boxed{
S_N
=
-A^2
\sum_{\langle ij\rangle}
\cos(\phi_j-\phi_i)
}
$$

$$
\Downarrow
$$

$$
\boxed{
F
=
-\nabla_\phi S_N/A^2
}
$$

$$
\Downarrow
$$

projection onto the zero-closure tangent space:

$$
\boxed{
F_T=P_\phi F
}
$$

$$
\Downarrow
$$

finite retraction:

$$
\boxed{
\phi^{(n+1)}
=
R_{\phi^{(n)}}(\eta F_T)
}
$$

$$
\Downarrow
$$

$$
\boxed{
\mathcal F_N:
\mathcal Z_N
\to
\mathcal Z_N
}
$$

$$
\Downarrow
$$

self-consistent fixed point:

$$
\boxed{
P_{\phi_*}F(\phi_*)=0
}
$$

$$
\Downarrow
$$

unconstrained current sector:

$$
\boxed{
\operatorname{div}_dJ=0
}
$$

$$
\Downarrow
\quad
N\to\infty
$$

$$
\boxed{
J^\mu
=
A^2g^{\mu\nu}\partial_\nu\phi
}
$$

$$
\boxed{
\partial_\mu J^\mu=0
}
$$

$$
\boxed{
\Box_g\phi=0
}
$$



---

## 24. Conclusion

The preceding paper derived many symmetry structures from zero closure, finite recurrence, simplex closure, and self-consistency, but the connections to Noether conservation laws, local gauge dynamics, and chirality selection remained.

In this paper we first constructed the discrete current and discrete action from relational phases and, by projection onto the zero-closure tangent space and retraction, constructed the zero-closure-preserving self-map

$$
\mathcal F_N:
\mathcal Z_N\to\mathcal Z_N
$$



This result is the central conclusion of this paper concerning dynamics. **In the discrete system one can construct a dynamics that changes the state by updating the phases $\phi$ without destroying zero closure.** In many continuum theories the state update is formally given first; in this system, since self-consistency is a foundational condition, such an implicit update is not allowed. Only after the update rule itself is shown to satisfy

$$
\boxed{
\mathcal F_N(\mathcal Z_N)\subseteq\mathcal Z_N
}
$$

can it be adopted as a dynamics of this axiom system. In this paper this was closed analytically for zero closure, and for the other previously derived conditions it was made a mandatory audit condition in the numerical verification.

The finite-$N$ discrete continuity equation goes over, in the high-resolution limit

$$
N\to\infty
$$

, to

$$
\partial_\mu J^\mu=0
$$

, and the continuum action connected to

$$
S_{\mathrm{cont}}
=
\frac{A^2}{2}
\int
g^{\mu\nu}
\partial_\mu\phi
\partial_\nu\phi
\,d^dq
$$



Further, from the anonymity of the local phase origin we obtained the edge connection and from the simplex faces the curvature, connecting to

$$
D_\mu
=
\partial_\mu-igA_\mu
$$

and

$$
F_{\mu\nu}
=
\partial_\mu A_\nu-\partial_\nu A_\mu-ig[A_\mu,A_\nu]
$$



Combining the five complex degrees of freedom and the $3\oplus2$ decomposition of the preceding paper,

$$
S(U(3)\times U(2))
$$

appears, and the trace-zero condition

$$
3y_3+2y_2=0
$$

fixes the hypercharge ratio.

Further,

$$
V^*\oplus\Lambda^2V
$$

decomposes into

$$
d^c,\quad L,\quad u^c,\quad Q,\quad e^c
$$

— the 15 left-handed Weyl components of one Standard-Model generation — and for that representation we directly verified the cancellation of the

$$
SU(3)^3,
\quad
SU(3)^2U(1),
\quad
SU(2)^2U(1),
\quad
U(1)^3,
\quad
\mathrm{grav}^2U(1)
$$

anomalies and the absence of the $SU(2)$ global anomaly.

Finally, we identified the conjugate two Weyl sectors

$$
\mathcal H_L
=
\left(\frac12,0\right)
\otimes
(V^*\oplus\Lambda^2V)
$$

$$
\mathcal H_R
=
\left(0,\frac12\right)
\otimes
(V\oplus\Lambda^2V^*)
$$

with the A/B two channels of the earlier papers.

Under the mirror transformation,

$$
J\to-J,
\qquad
\Delta\phi\to-\Delta\phi,
\qquad
S_\chi\to-S_\chi
$$

so the minimal selection equation near the symmetric point is

$$
\boxed{
\dot S_\chi
=
\lambda J
+
gS_\chi(1-S_\chi^2)
}
$$



For $g>0$,

$$
S_\chi=0
$$

is unstable and

$$
S_\chi=\pm1
$$

are stable.

Hence, without putting any left-right preference into the fundamental equations, chirality is spontaneously selected with the sign of the internal phase correlation $J$ as the seed.

From the above, within the range targeted by this paper, the theoretical connection to the **gauge-representation/chirality structure of the Standard Model** has closed as

$$
\boxed{
\text{zero closure}
\to
\text{Noether conservation laws}
\to
\text{local gauge geometry}
\to
S(U(3)\times U(2))
\to
\text{hypercharge}
\to
\text{one-generation representation}
\to
\text{anomaly cancellation}
\to
\text{chirality selection}
}
$$



Hence the main task remaining for this derivation chain is not the addition of new theoretical terms but the **numerical verification** that, on the existing numerical experimental system,

$$
\boxed{
J_\chi
\to
\Delta\phi_\chi
\to
S_\chi
\to
\pm1
}
$$

is reproduced while mirror symmetry, zero closure, and the conservation laws are preserved.

The subsequent paper will, without changing the existing A/B Fermi-type experimental system, add observables and mirror runs, and carry out the direct fit of

$$
\dot S_\chi
=
\lambda J_\chi
+
gS_\chi(1-S_\chi^2)
$$

, phase-closure breaking, the Bose/Fermi control, the causal order, and the conservation-law audit.

Note that "theoretical connection to the Standard Model" in this paper refers to the chain gauge group, local gauge geometry, one-generation internal representation, hypercharge, anomaly cancellation, and chirality selection. The Higgs radial mode, Yukawa couplings, the number of generations, the concrete mass hierarchy, and quantum corrections / renormalisation remain separate theoretical problems.

## References

1. N. Kihara, "Must the symmetries of physics really be given from the start?" — the preceding paper (see [10]).
2. E. Noether, “Invariante Variationsprobleme,” *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 1918, 235–257.
3. J. E. Marsden and M. West, “Discrete Mechanics and Variational Integrators,” *Acta Numerica*, 10 (2001), 357–514. DOI: 10.1017/S096249290100006X.
4. V. A. Dorodnitsyn, “Noether-type theorems for difference equations,” *Applied Numerical Mathematics*, 39 (2001), 307–321. DOI: 10.1016/S0168-9274(00)00041-6.
5. M. Skopenkov, “Discrete Field Theory: Symmetries and Conservation Laws,” *Mathematical Physics, Analysis and Geometry*, 26 (2023), Article 19. DOI: 10.1007/s11040-023-09459-4.
6. S. Navas et al. (Particle Data Group), "Review of Particle Physics," *Phys. Rev. D* 110, 030001 (2024). DOI: 10.1103/PhysRevD.110.030001.
7. M. E. Peskin and D. V. Schroeder, *An Introduction to Quantum Field Theory*, Addison-Wesley (1995), Chapter 20 (gauge theories with spontaneous symmetry breaking; Standard Model representation content and anomaly cancellation).
8. K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* 10, 2445 (1974). DOI: 10.1103/PhysRevD.10.2445 (link variables, plaquette action, continuum Yang–Mills limit).
9. H. Georgi and S. L. Glashow, "Unity of All Elementary-Particle Forces," *Phys. Rev. Lett.* 32, 438 (1974). DOI: 10.1103/PhysRevLett.32.438 (the $\overline{\mathbf 5}\oplus\mathbf{10}$ one-generation content).
10. N. Kihara, "Symmetry Generation from Zero Closure, Finite Order, and Self-Consistent Geometry — The Single External Parameter $N$ and the Remaining Tasks of Generalization and Dynamics," Zenodo, Concept DOI: 10.5281/zenodo.22028072, Version DOI: 10.5281/zenodo.22028073 (2026). (Cited in the text as [1] / "the preceding paper".)
