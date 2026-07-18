# Emergence of Discrete Born-Type Weights in Iterated Two-Channel Exchange Systems

## A Finite-Order Recurrence Law from Wave-Packet Localization Transfer, Metastable Two-State Dynamics, and Observation Selection

**Version:** English complete manuscript v1<br>
**Date:** July 18, 2026<br>
**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Version DOI:** 10.5281/zenodo.21422471<br>
**Concept DOI:** 10.5281/zenodo.21422470<br>
**Position in the series:** Independent report on finite-order recurrence and Born-type weights in the Wave Information Readout series

---

## Abstract

This study was not initiated to reproduce the Born rule or its two-state representation, the $\cos^2$ law. It began with two different numerical experiments. System A examined whether iterative exchange scattering between a weakly localized wave and a localized wave containing higher harmonics could transfer localization and effective harmonic structure between two channels, producing behavior resembling wave-packet contraction or localization. System B examined two A/B states, an intermediate gray metastable state, its retention under weak readout, and its selection into A or B under strong observation.

Although the two models were built to address different physical questions and used different observation functions, later analysis showed that their linear mixing parts, before normalization and observation operations, share the same two-channel exchange-scattering kernel

$$
U_R=
\begin{pmatrix}
r&t\\
t&r
\end{pmatrix}.
$$

System A applies channel-wise normalization after this linear mixing, while the strong D observation in System B applies a separate back-action map; the complete update maps of the two systems are therefore not assumed to be identical. A full-range sweep of the exchange coefficient $R$ in System B produced extremely sharp peaks at particular values. Eigenvalue analysis showed that the peaks arise from exact recurrence in which the antisymmetric eigenvalue returns to unity after finitely many iterations:

$$
U_R^n=I.
$$

We derive this result without depending on the trigonometric parameterization of the scattering amplitudes. Exchange symmetry determines the symmetric and antisymmetric projectors

$$
P_s=\frac12
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix},
\qquad
P_a=\frac12
\begin{pmatrix}
1&-1\\
-1&1
\end{pmatrix}.
$$

After fixing the global phase so that the symmetric eigenvalue is one, an exchange-symmetric unitary operator can be written as

$$
U=P_s+\zeta P_a,
\qquad |\zeta|=1.
$$

When the finite-order condition selects $\zeta=e^{-2\pi im/n}$, the diagonal and off-diagonal amplitudes in the A/B basis are

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2}.
$$

Consequently, for coprime integers $m,n$, the exchange weight satisfying the finite-order condition is

$$
\boxed{
R_{n,m}=\cos^2\left(\frac{\pi m}{n}\right)
}.
$$

This expression was not obtained by entering a Born-type squared law as a search condition or evaluation function. The projector representation above also shows that the trigonometric form need not be imposed as an independent assumption. The iterative exchange systems modeling wave-packet localization transfer, metastable two-state dynamics, weak readout, and strong observation selection were constructed first; numerical peaks were then found; only afterward was their closed-orbit condition analyzed. The result is therefore not a complete derivation of the Born rule. It shows that, in a closed two-channel wave system exhibiting behavior analogous to phenomena associated with quantum measurement, finitely recurrent exchange weights arise endogenously as a discrete Born-type $\cos^2$ series.

The current code represents the unitary scattering amplitudes using $\sin\theta$ and $\cos\theta$. This is a coordinate representation recoverable from the exchange-symmetric projectors $P_s,P_a$ and the eigenvalues $1,\zeta$, rather than an independent assumption required by the central theorem. No angle was prescribed in advance as the one that must recur exactly. Instead, the finite-order condition selects the discrete Born-angle series

$$
\phi_{n,m}=\frac{\pi m}{n}.
$$

The scattering angle used in the implementation is $\theta=\pi/2-\phi_{n,m}$. This distinguishes the exchange-symmetric two-channel projection structure from the phase series selected by closure: the latter discretizes the admissible weights of the former.

This paper does not focus on numerical correspondences with the fine-structure constant, $E_8$, or any specific physical constant. Its central purpose is to report clearly the mathematical structure connecting contraction-like localization redistribution, metastable superposition, observation selection, finite-order recurrence, and Born-type squared weights through a common unperturbed exchange kernel $U$ followed by model-specific normalization, readout, and selection maps.

**Keywords:** Born rule, two-channel exchange, finite-order recurrence, wave-packet localization, metastable state, observation selection, quantum measurement, unitary operator, discrete phase

---

# 1. Background and Objective

## 1.1 The Born Rule and Its Two-State Representation

In standard quantum mechanics, the probability of projecting a state $|\psi\rangle$ onto an observation state $|A\rangle$ is given by the Born rule

$$
P(A)=|\langle A|\psi\rangle|^2
$$

[5,6].

For a two-state system written as

$$
|\psi\rangle
=
\cos\phi\,|A\rangle
+
e^{i\chi}\sin\phi\,|B\rangle,
$$

the observation probabilities of the two orthogonal states are

$$
P(A)=\cos^2\phi,
\qquad
P(B)=\sin^2\phi.
$$

The Born rule is the fundamental rule connecting quantum theory to observation results. The present work, however, did not begin as a derivation of that rule. The numerical systems were first constructed to study wave-packet localization transfer and observation selection in a metastable two-state system; they did not contain the $\cos^2$ law as a target function.

## 1.2 Research Question

The question addressed here is:

> In an iterated exchange system that models contraction-like localization redistribution, metastable mixing of two A/B states, retention under weak readout, and selection under strong observation, can Born-type squared weights arise from exchange symmetry, unitarity, and finite-order closure without being externally entered as a probability axiom?

This paper does not seek a uniqueness proof of the general Born rule or a complete dynamical derivation of single-trial probabilities. Its specific objective is to decompose an already constructed exchange-symmetric two-channel unitary operator into symmetric and antisymmetric subspaces, derive the endogenous discrete Born-type weights

$$
R_{n,m}=\cos^2\left(\frac{\pi m}{n}\right)
$$

from its finite-order condition, and state precisely the meaning and limits of this result.

## 1.3 The Non-Directed Order of Discovery

The chronological order of the discovery is important to its interpretation.

1. System A was constructed to transfer wave-packet localization and harmonic structure between two channels [1,2].
2. System B was constructed to study two A/B states, a gray metastable state, weak readout, and strong observation selection [3].
3. A local sweep was performed for localization transfer in System A, and a full-range sweep of their common exchange-scattering coefficient $R$ was performed in System B [4].
4. Extremely sharp peaks were observed in the full-range System B sweep.
5. Local high-resolution sweeps and multiprecision calculations were performed.
6. Operator eigenvalues were analyzed, identifying the peaks as finite-order recurrence [4].
7. The finite-order condition was recognized to give the exchange weights $\cos^2(\pi m/n)$.
8. Only afterward was this expression reinterpreted physically as isomorphic to the two-state representation of the Born rule.

Thus, this study did not embed a Born-type $\cos^2$ law in the code and rediscover it in the output. The amplitude structure of unitary two-channel scattering was present in the original models, but neither the phase angle nor the exchange weight selected for exact recurrence was entered in advance. Moreover, this paper removes the trigonometric representation itself from the central assumptions and derives

$$
U=P_s+\zeta P_a,
\qquad |\zeta|=1
$$

from exchange symmetry and unitarity. Imposing the finite-order condition on this representation makes the $\cos^2/\sin^2$ weights a necessary consequence of interference between the symmetric and antisymmetric eigenphases, rather than of a chosen coordinate parameterization.

---

# 2. The Two Starting Models

## 2.1 System A: Wave-Packet Localization Transfer

System A iterates the same exchange scattering between a weakly localized wave on side A and a localized wave containing higher harmonics on side B [1,2].

Let the wave-packet vectors at step $j$ be $A_j,B_j$. The linear exchange part is

$$
\begin{pmatrix}
\widetilde A_{j+1}\\
\widetilde B_{j+1}
\end{pmatrix}
=
U_R
\begin{pmatrix}
A_j\\
B_j
\end{pmatrix}.
$$

The implemented one-step update in System A is not only this linear mixing. If $\mathcal N(x):=x/\|x\|$ denotes normalization of each wave-packet channel, the actual update is

$$
\boxed{
F_R(A,B)
=
\left(
\mathcal N(rA+tB),
\mathcal N(tA+rB)
\right)
}.
$$

Thus, $F_R$ is generally nonlinear. The finite-order theorem in this paper applies directly to the pre-normalization linear numerator

$$
(A,B)\longmapsto(rA+tB,\ tA+rB).
$$

The localization transfer observed in System A is treated as a property of the complete map combining this common linear kernel with channel-wise normalization.

Representative observables are

$$
L,
\qquad
N_{\mathrm{eff}},
\qquad
B_{\mathrm{to}A},
$$

where

- $L$ measures wave-packet localization,
- $N_{\mathrm{eff}}$ is the effective harmonic order, and
- $B_{\mathrm{to}A}$ measures how much of the initial harmonic structure on side B has moved to side A.

At intermediate exchange coefficients, localization and harmonic structure transfer periodically between the two channels. If only one channel is observed at a selected time, a broad wave packet can appear to have changed into a localized wave. This behavior was constructed as redistribution of localization by exchange interference, not as nonunitary contraction accompanied by disappearance of the norm of the total system.

## 2.2 System B: White-Cat, Black-Cat, and Gray-Cat Metastable Interface

System B defines, from two complex amplitudes $a,b$,

$$
p_A=|a|^2,
\qquad
p_B=|b|^2,
\qquad
S=p_A-p_B
$$

[3].

The states are read as follows.

- $S\approx+1$: A-dominant, or white-cat, state.
- $S\approx-1$: B-dominant, or black-cat, state.
- $S\approx0$: balanced A/B, or gray, state.

Among gray states, a state that fluctuates periodically with small amplitude, remains under weak C readout, and is selected into A or B by strong D observation is called a gray metastable phase.

The unperturbed A/B exchange in System B follows $U_R$. The implemented pair normalization is

$$
\mathcal N_{AB}(a,b)
=
\frac{(a,b)}{\sqrt{|a|^2+|b|^2}},
$$

which is the identity under exact unitary evolution. The weak C readout and strong D observation are separate readout and back-action maps applied after the unperturbed exchange. In particular, the D implementation explicitly contains nonlinear back-action that updates the allocation difference $S$ as

$$
S_{j+1}
=
S_j+g_D S_{D,j}(1-S_j^2),
$$

pushing the state toward either A or B. Strong observation selection therefore does not arise from attraction by $U_R$ alone, but from the composition of $U_R$ with D back-action.

This model was not constructed as a direct solution to the quantum-measurement problem. Its purpose was to test whether one iterative exchange system could distinguish

- a metastable state intermediate between two states,
- a weak readout that does not substantially destroy the state, and
- a strong observation that changes the state toward one side.

## 2.3 The Scattering Kernel Shared by the Two Models

Later analysis showed that System A and System B, despite their different state spaces, normalizations, and observation functions, use the same exchange-scattering kernel in the linear mixing stage before normalization, readout, and back-action:

$$
\boxed{
U_R=
\begin{pmatrix}
r&t\\
t&r
\end{pmatrix}
}.
$$

The complete update structures can be decomposed conceptually as

$$
\text{System A}:\quad F_R=\mathcal N_A\circ U_R,
$$

$$
\text{System B, unperturbed}:\quad U_R,
$$

$$
\text{System B, C/D}:\quad
F_C=\mathcal C\circ U_R,
\qquad
F_D=\mathcal D\circ U_R.
$$

Here $\mathcal N_A$ is channel-wise normalization, $\mathcal C$ is weak readout and weak back-action, and $\mathcal D$ is strong selection back-action. System A reads localization and internal harmonic structure of wave-packet vectors, whereas System B reads the A/B allocation difference and metastability of two complex amplitudes.

Therefore, candidate positions of particular $R$ values common to both models originate in the intrinsic structure of the exchange operator, while peak visibility, localization transfer, and A/B selection are determined by the downstream maps and observation functions. This separation between the common kernel that fixes root positions and the model-specific maps that visualize or select phenomena is the central implementation structure of this paper.

---

# 3. Exchange-Scattering Operator

## 3.1 Scattering-Coefficient Coordinates Used in the Numerical Implementation

For a reflection weight $R\in[0,1]$, the implementation defines

$$
\theta(R):=\arcsin\sqrt R.
$$

The transmission amplitude $t$ and reflection amplitude $r$ are

$$
t=e^{i\theta}\cos\theta,
\qquad
r=-ie^{i\theta}\sin\theta.
$$

Then

$$
|t|^2=\cos^2\theta=1-R,
$$

$$
|r|^2=\sin^2\theta=R.
$$

Moreover,

$$
|r|^2+|t|^2=1,
\qquad
r^*t+t^*r=0,
$$

so that

$$
U_R^\dagger U_R=I.
$$

Thus, the linear exchange operator at fixed $R$ is unitary. This trigonometric representation is a convenient coordinate system used in the implementation, but it is not an independent axiom needed to obtain the $\cos^2/\sin^2$ series. Sections 4 and 5 reconstruct the same coefficient representation using only the projector decomposition of an exchange-symmetric unitary operator.

## 3.2 State Norm and Localization

The unitarity of the linear kernel $U_R$ preserves the norm of the complete two-channel state:

$$
\left\|
U_R
\begin{pmatrix}
A\\B
\end{pmatrix}
\right\|^2
=
\left\|
\begin{pmatrix}
A\\B
\end{pmatrix}
\right\|^2.
$$

Therefore, even if one wave packet appears to localize under the pre-normalization linear kernel, this does not mean that the total norm disappears or that the state is dissipatively attracted to a point. Localization, harmonic structure, or channel allocation is redistributed within the conserved total state. Since the complete System A map includes channel-wise normalization, this total-norm equation cannot be transferred unchanged to the full update; localization transfer must be evaluated for $F_R=\mathcal N_A\circ U_R$.

This distinction is essential when discussing the analogy with wave-packet contraction. “Contraction-like behavior” in this paper is not claimed to be identical to post-measurement state update in standard quantum mechanics. It denotes a mathematical analogy in which concentration of localization emerges from the composition of a unitary exchange kernel, channel-wise normalization, and partial observation.

---

# 4. Symmetric and Antisymmetric Eigenmodes

## 4.1 Decomposition into Eigenmodes

Decompose the two-channel state as

$$
X_j:=\frac{A_j+B_j}{\sqrt2},
\qquad
Y_j:=\frac{A_j-B_j}{\sqrt2}.
$$

$X_j$ is symmetric under channel exchange, and $Y_j$ is antisymmetric.

The eigenvalues of the exchange operator are

$$
\lambda_s=r+t,
\qquad
\lambda_a=r-t.
$$

Substituting the definitions of $r,t$ gives

$$
\lambda_s=1,
$$

$$
\lambda_a=-e^{2i\theta}.
$$

Therefore,

$$
X_j=X_0,
\qquad
Y_j=\lambda_a^jY_0.
$$

## 4.2 Projector Representation

Let the channel-exchange operator be

$$
X=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix}
=P_s-P_a.
$$

Exchange symmetry is

$$
[U,X]=0.
$$

Because the eigenvalues $+1,-1$ of $X$ are distinct, any operator $U$ satisfying this commutation relation does not mix the symmetric and antisymmetric subspaces.

Define the symmetric and antisymmetric projectors by

$$
P_s=
\frac12
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix},
\qquad
P_a=
\frac12
\begin{pmatrix}
1&-1\\
-1&1
\end{pmatrix}.
$$

Then

$$
P_s+P_a=I,
\qquad
P_sP_a=0,
$$

and

$$
U=\lambda_sP_s+\lambda_aP_a.
$$

If $U$ is unitary,

$$
|\lambda_s|=|\lambda_a|=1.
$$

Removing the global phase $\lambda_s$ by defining $\widetilde U:=\lambda_s^{-1}U$ gives

$$
\boxed{
\widetilde U=P_s+\zeta P_a,
\qquad
\zeta:=\frac{\lambda_a}{\lambda_s},
\qquad
|\zeta|=1
}.
$$

In the numerical implementation $\lambda_s=1$ from the outset, so below we write $\widetilde U$ as $U_R$. Then

$$
U_R
=
\frac12
\begin{pmatrix}
1+\zeta&1-\zeta\\
1-\zeta&1+\zeta
\end{pmatrix}
=
\begin{pmatrix}
r&t\\
t&r
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

Thus, without separately assuming a trigonometric form, the exchange amplitudes are fixed as the sum and difference of the symmetric eigenphase $1$ and the antisymmetric eigenphase $\zeta$.

Furthermore,

$$
\boxed{
U_R^j=P_s+\zeta^jP_a
}.
$$

This expression shows that the apparently complicated time evolution of the model consists of only

- a fixed symmetric component, and
- an antisymmetric component rotating on the unit circle.

## 4.3 Periodicity of Observables

The allocation difference in System B can generally be written as

$$
S_j=C\cos(j\omega+\varphi_0),
$$

where

$$
\omega=\arg\lambda_a
$$

is determined by the operator, while $C$ and $\varphi_0$ are determined by the initial state and observation function.

Changing the initial amplitude, initial phase, internal harmonic distribution, or wave-packet shape therefore does not change the operator-intrinsic phase-rotation number $\omega(R)$. What changes is how visibly that eigenphase appears in the chosen observable.

---

# 5. Finite-Order Recurrence

## 5.1 Exact Recurrence Condition

In the numerical implementation, the symmetric eigenvalue is exactly one. The condition under which this phase-fixed iterated exchange system restores every initial state after $n$ iterations is

$$
U_R^n=I.
$$

Since the symmetric eigenvalue satisfies

$$
\lambda_s^n=1,
$$

the necessary and sufficient condition is

$$
\lambda_a^n=1
$$

for the antisymmetric eigenvalue. For a general exchange-symmetric $U\in U(2)$, $U^n=e^{i\gamma}I$ is projective recurrence up to a global phase of the physical state. The equation $\widetilde U^n=I$ for $\widetilde U=\lambda_s^{-1}U$ is its phase-fixed representation used below.

## 5.2 Exchange-Symmetric Finite-Order Theorem

**Theorem 5.1 (Discrete squared weights of an exchange-symmetric finite-order map).**  
Let a two-channel operator $U\in U(2)$ satisfy the following conditions.

1. It commutes with the channel-exchange operator $X$: $[U,X]=0$.
2. The global phase is fixed by using $\widetilde U:=\lambda_s^{-1}U$.
3. $\widetilde U$ has a nontrivial fundamental order $n\ge3$: $\widetilde U^n=I$.
4. The antisymmetric eigenphase is $\zeta=e^{-2\pi i m/n}$ with $\gcd(m,n)=1$. Since the conjugate pair $m$ and $n-m$ gives the same weights, choose $1\le m<n/2$ as a representative.

Then

$$
\widetilde U
=
P_s+e^{-2\pi i m/n}P_a
=
\begin{pmatrix}
r_{n,m}&t_{n,m}\\
t_{n,m}&r_{n,m}
\end{pmatrix},
$$

$$
r_{n,m}
=
\frac{1+e^{-2\pi i m/n}}{2},
\qquad
t_{n,m}
=
\frac{1-e^{-2\pi i m/n}}{2},
$$

and the channel weights are

$$
\boxed{
|r_{n,m}|^2
=
\cos^2\left(\frac{\pi m}{n}\right),
\qquad
|t_{n,m}|^2
=
\sin^2\left(\frac{\pi m}{n}\right)
}.
$$

*Proof.* The relation $[U,X]=0$ diagonalizes $U$ on the symmetric and antisymmetric eigenspaces of $X$. After fixing the global phase,

$$
\widetilde U=P_s+\zeta P_a.
$$

The condition $\widetilde U^n=I$ gives $\zeta^n=1$. Since the fundamental order is $n$, write $\zeta=e^{-2\pi im/n}$ with $\gcd(m,n)=1$. Expanding the projectors in the channel basis gives

$$
r_{n,m}=\frac{1+\zeta}{2},
\qquad
t_{n,m}=\frac{1-\zeta}{2}.
$$

Therefore,

$$
|r_{n,m}|^2
=
\frac{(1+\zeta)(1+\zeta^*)}{4}
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
\qquad\square
$$

## 5.3 Agreement with the Implementation Coordinates and Boundary Roots

In the implementation coordinates of Section 3,

$$
\zeta=\lambda_a=-e^{2i\theta}.
$$

Identifying this with $e^{-2\pi im/n}$ gives

$$
\theta
=
\frac{\pi}{2}-\frac{\pi m}{n}
\pmod{\pi},
$$

and once again

$$
R=|r|^2=\sin^2\theta
=
\cos^2\left(\frac{\pi m}{n}\right).
$$

The trigonometric form is therefore not an input to the theorem. It is the implementation coordinate expressing the coefficients obtained from the projector decomposition in terms of $R$.

Theorem 5.1 concerns nontrivial interior roots $0<R<1$. At the boundaries,

$$
\zeta=1:\quad U=I,\quad R=1,
$$

$$
\zeta=-1:\quad U=X,\quad R=0,
$$

which are the trivial closures of order one and two, respectively.

## 5.4 A Discrete Series inside a Continuous Parameter

Formally, $R$ ranges over the continuous interval

$$
R\in[0,1].
$$

The points of exact recurrence, however, are classified by integer pairs $(n,m)$:

$$
(n,m)
\longmapsto
R_{n,m}
=
\cos^2\left(\frac{\pi m}{n}\right).
$$

Thus, within an apparently continuous exchange system there is a countable family of closed orbits labeled by rational phase ratios.

Because rational numbers $m/n$ are dense in the interval, finite-order roots also become dense when arbitrarily high orders are allowed. At each root, however,

$$
U_R^n=I
$$

holds exactly: the orbit is an exact discrete closed orbit, not an approximately periodic one.

---

# 6. Discrete Born-Type Weights

## 6.1 Formal Isomorphism with Standard Two-State Projection

In a standard two-state system,

$$
|\psi\rangle
=
\cos\phi\,|A\rangle
+
e^{i\chi}\sin\phi\,|B\rangle,
$$

the Born rule gives

$$
P(A)=\cos^2\phi,
\qquad
P(B)=\sin^2\phi.
$$

Applying the operator studied here to the channel-basis state $|A\rangle=(1,0)^T$ gives

$$
U_{n,m}|A\rangle
=
r_{n,m}|A\rangle+t_{n,m}|B\rangle.
$$

Consequently, the two transition weights obtained by projecting the output onto the same A/B basis are exactly

$$
W_{A\to A}
=
|\langle A|U_{n,m}|A\rangle|^2
=
|r_{n,m}|^2,
$$

$$
W_{A\to B}
=
|\langle B|U_{n,m}|A\rangle|^2
=
|t_{n,m}|^2.
$$

By Theorem 5.1,

$$
W_{A\to A}(n,m)
=
\cos^2\phi_{n,m},
$$

$$
W_{A\to B}(n,m)
=
\sin^2\phi_{n,m},
$$

where

$$
\phi_{n,m}:=\frac{\pi m}{n}.
$$

Thus, the channel-projection weights of a single exchange acting on a basis input are

$$
\boxed{
W_{A\to A}(n,m)=\cos^2\phi_{n,m},
\qquad
W_{A\to B}(n,m)=\sin^2\phi_{n,m}
}
$$

and are isomorphic to the standard two-state Born rule. This is not merely a visual resemblance between formulas; it is the same operation of projecting a complex amplitude onto an observation basis and squaring its absolute value.

For a general superposed input $\alpha|A\rangle+\beta|B\rangle$, however, the A-channel output is $r\alpha+t\beta$, whose weight includes an interference term. Therefore, $R=|r|^2$ is not itself the probability of observing A for an arbitrary input state; it is the transition weight for remaining in the A channel when the input is the A basis state. A separate frequency law is also required to interpret this projection weight as an experimental frequency over repeated trials.

## 6.2 What Was Assumed and What Emerged

To evaluate the result precisely, we separate the structures already present in the initial definitions from those that emerged only after analysis.

### Structures present in the initial definitions

1. Two-channel complex amplitudes.
2. A unitary scattering condition.
3. Exchange symmetry between the A and B channels.
4. Reading squared amplitude magnitudes in the A/B basis as channel weights.
5. A trigonometric parametrization of the amplitudes $r,t$ as numerical implementation coordinates.

Thus, the structure in which the squared magnitude of a complex amplitude becomes a nonnegative channel weight was already part of the present model. Theorem 5.1 nevertheless shows that the trigonometric parametrization need not be an independent derivational assumption: it can be reconstructed from the projector decomposition of an exchange-symmetric unitary operator.

### Structures not present in the initial definitions

1. Which values of $R$ produce the special peaks.
2. Which phase angles are selected by exact recurrence.
3. That the roots are classified by integer pairs $(n,m)$.
4. That the selected exchange weights form the discrete sequence $\cos^2(\pi m/n)$.
5. That the localization-transfer model and the gray-state metastability index connect to finite-order roots of a common linear kernel.

The nontrivial consequence of this work is therefore not the generation of a square law without assumptions, but rather

$$
\boxed{
\text{finite-order closure}
\Rightarrow
\text{discrete phase sequence}
\Rightarrow
\text{discrete Born-type weight sequence}
}.
$$

## 6.3 Closure-Selected Weights Preceding a Probability Postulate

Ordinarily, $\cos^2\phi$ is introduced as a measurement probability.

Here, the same form first appears as

> the reflection/transmission weight required for a two-channel exchange to close exactly after finitely many iterations.

That is, rather than beginning with

$$
\cos^2 \text{ as a probability},
$$

we first obtain

$$
\cos^2 \text{ as a completely recurrent exchange ratio}.
$$

This order suggests that a Born-type weight may be read not merely as an external probability postulate, but as a recurrence condition for a closed phase system.

The word “selection” here does not mean dynamical attraction. During the sweep, $R$ is not a time-evolving variable; each value of $R$ is fixed while the orbital closure error is evaluated. The finite-order roots are therefore not attractors of the form

$$
R_j\longrightarrow R_{n,m},
$$

but parameters selected by the closure criterion

$$
U_R^n=I.
$$

We call this “closure selection” or a “recurrence-compatible weight.” Dynamical selection between A and B under strong observation is a separate mechanism arising from the D backreaction in System B.

## 6.4 The Discrete Sequence and the Continuous Born Rule

What is obtained directly here is the discrete sequence associated with

$$
\phi_{n,m}=\frac{\pi m}{n}.
$$

As $n$ increases, rational phases become dense in the continuum of angles, and hence the set

$$
\cos^2\left(\frac{\pi m}{n}\right)
$$

is dense in $[0,1]$.

Density, however, is not equivalent to having physically derived the continuous Born rule. To assert a continuous limit as a physical law, one must separately establish

- which orders are physically realizable,
- how well high-order roots survive noise,
- how observation and coherence times select roots, and
- whether experimental frequencies converge to the squared weights.

---

# 7. Relation to Wave-Packet-Contraction-Like Behavior

## 7.1 Redistribution of Localization Rather Than Nonunitary Contraction

In System A, when one wave packet changes from a weakly localized state to a strongly localized state, a partial observation can make the packet appear to contract.

The unnormalized exchange kernel $U_R$ is unitary; when the linear update alone is isolated, the total system preserves

$$
\|A_j\|^2+\|B_j\|^2
=
\text{constant}.
$$

The full System A update implemented numerically is, however,

$$
F_R(A,B)
=
\left(
\mathcal N(rA+tB),
\mathcal N(tA+rB)
\right),
$$

which normalizes each channel separately. Therefore, $F_R$ itself cannot be called a unitary operator.

What occurs at the linear exchange stage is consequently not

$$
\text{contraction accompanied by loss of the total state},
$$

but

$$
\text{redistribution of localization and internal structure within a conserved total state}.
$$

The localization transfer observed numerically in System A is this conservative redistribution as read after channel-wise normalization. Whether the finite-order roots are also exact periodic roots of the complete nonlinear map $F_R$ does not follow automatically from the theorem for $U_R$ and is a System-A-specific proposition requiring separate verification.

## 7.2 Contraction-Like Appearance Through an Observation Section

If an observer reads only one channel or a particular localization index rather than the entire two-channel state, one section of the combined orbit of the exchange kernel and normalization may appear as

- a transition from a diffuse to a localized state,
- selection from an intermediate state into A or B, or
- a sudden reduction in wave-packet width.

The model therefore offers the working hypothesis

$$
\boxed{
\text{apparent contraction}
=
\text{unitary exchange kernel}
+
\text{localization redistribution}
+
\text{channel-wise normalization}
+
\text{restricted observation section}
}.
$$

This is not a claim to have reproduced the state-update rule of standard quantum measurement. It is, however, a concrete mathematical model capable of producing contraction-like observations without postulating instantaneous nonunitary contraction as a fundamental operation.

## 7.3 Recurrent Localization Exchange at Finite-Order Roots

When the exchange weight is

$$
R=R_{n,m},
$$

the linear kernel $U_R$ restores every input exactly after $n$ iterations.

Thus, even if localization becomes concentrated in one channel at an intermediate stage of the unnormalized linear evolution, that change is not irreversible absorption but one phase of a closed orbit. For the full System A orbit including channel-wise normalization, the conditions under which this closure survives must be distinguished numerically.

This structure allows the following to be treated as phases of a single periodic orbit:

- a moment that appears localized,
- a moment that becomes delocalized again,
- a moment at which A and B are interchanged, and
- the moment of exact return to the initial state.

---

# 8. Relation to the White-Cat, Black-Cat, and Gray-Cat System

## 8.1 A Metastable Two-State System

In System B, the allocation difference

$$
S_j=|A_j|^2-|B_j|^2
$$

represents the dominance of the two states A and B.

Because the antisymmetric eigenphase rotates, it generally has the form

$$
S_j=C\cos(j\omega+\varphi_0).
$$

The gray state is therefore not merely a static half-and-half mixture of A and B; it can be represented as a dynamically metastable state in which the A/B allocation difference is exchanged periodically with small amplitude.

## 8.2 Weak Readout and Strong-Observation Selection

A one-step update in System B can be written by separating the unperturbed exchange kernel from the observation maps:

$$
\begin{aligned}
\Psi_{j+1}^{(C)}
&=\mathcal C\!\left(U_R\Psi_j^{(C)}\right),\\
\Psi_{j+1}^{(D)}
&=\mathcal D\!\left(U_R\Psi_j^{(D)}\right).
\end{aligned}
$$

Under weak C readout, the observational action is small and reads the allocation difference without substantially changing the iterated orbit. Strong D observation has an explicit nonlinear backreaction that uses the read value $S_{D,j}$ to update

$$
S_{j+1}
=
S_j+g_D S_{D,j}(1-S_j^2),
$$

thereby driving the system toward a state dominated by either A or B. The direction and irreversibility of A/B selection therefore do not arise from the finite-order kernel $U_R$ alone; they depend on the structure of $\mathcal D$ and on the coupling strength $g_D$.

This construction is not a direct implementation of the standard theories of weak and strong measurement. It nevertheless distinguishes, within the same model,

- a readout that tends to preserve the state, and
- a selective operation that changes the state.

Here, $U_R$ determines the locations of the finite-order roots, $\mathcal C$ reads the gray orbit while disturbing it as little as possible, and $\mathcal D$ selects A or B. This division of roles makes it possible to connect recurrence and observation selection without confusing them as effects of the same operator.

## 8.3 The Context in Which the Born-Type Weights Appeared

It is important that the $\cos^2$ sequence was not discovered from an abstract two-state vector alone.

It emerged while identifying the cause of peaks in numerical experiments that evaluated the finite order of the unperturbed exchange kernel in a system containing

1. two states A and B,
2. an intermediate gray metastable state,
3. preservation under weak readout,
4. A/B selection by a separately implemented strong-observation backreaction, and
5. closed orbits with small long-time drift.

The principal significance of this work is therefore that, in a family of mathematical models resembling the physical context in which a Born-type square law is ordinarily used, the same form of weight emerges from the finite-order condition of the unperturbed kernel, and the C/D observation maps can be connected to that same kernel.

---

# 9. Numerical Peaks and Analytic Roots

## 9.1 Peak Metric

For a time series of length $k$ in System B, we evaluated

$$
\overline S_k
=
\frac1k\sum_{j=0}^{k-1}S_j,
$$

$$
A_k
=
\frac{\max_{j<k}S_j-\min_{j<k}S_j}{2},
$$

$$
D_k
=
\left|
\overline S_{k,\mathrm{second}}
-
\overline S_{k,\mathrm{first}}
\right|.
$$

The principal conditions used to confirm the finite-order roots directly were the initial relative phase

$$
\varphi_0=0\ \text{or}\ \pi,
$$

the code-level initial difference $s_0=0.01$, and hence the allocation difference

$$
S_0=|a_0|^2-|b_0|^2=0.02,
$$

with target amplitude $C=0.02$. Under these conditions, the gray error was defined as

$$
\varepsilon_k
=
|\overline S_k|
+
|A_k-C|
+
D_k,
$$

and the depth as

$$
d_k=-\log_{10}\varepsilon_k.
$$

This metric does not evaluate agreement with a Born probability. It evaluates

- neutrality of the mean,
- preservation of the oscillation amplitude, and
- absence of drift between the first and second halves.

## 9.2 Even-Order Roots

Under the principal condition $\varphi_0=0$ or $\pi$, the allocation difference at a root of Theorem 5.1 is

$$
S_j
=
C\cos\left(\frac{2\pi m}{n}j\right),
\qquad C=0.02.
$$

If the primitive order $n$ is even and $\gcd(m,n)=1$, then $m$ is odd, so the lattice points contain exactly

$$
S_0=C,
\qquad
S_{n/2}=C\cos(\pi m)=-C.
$$

The sum of the cyclotomic phases over one period is also zero.

Therefore, when two complete periods $j=0,\ldots,2n-1$, comprising $2n$ samples, are taken, the current gray metric satisfies exactly

$$
\overline S=0,
\qquad
A=C,
\qquad
D=0.
$$

This does not assert that the gray error vanishes for every initial state and every evaluation window whenever the order is even. Operator recurrence $U_R^n=I$ is independent of the initial state, but the amplitude condition and zero error above are exact results under the experimental conditions $\varphi_0=0$ or $\pi$, $S_0=C=0.02$, and a two-complete-period sampling window.

Consequently,

$$
\varepsilon=0,
$$

and in ideal arithmetic

$$
d\rightarrow+\infty.
$$

What diverges is not a physical energy or amplitude but the negative logarithm of an error metric for a perfectly periodic orbit. Finite-precision calculations display a finite depth because of rounding error.

## 9.3 Odd-Order Roots and Observable Dependence

For an odd primitive order $n$,

$$
U_R^n=I
$$

still holds.

If, however, the finite sampled sequence does not contain both positive and negative extrema as lattice points, a small deficit remains in the present amplitude estimate.

Therefore,

$$
\boxed{
\text{exact operator recurrence}
\neq
\text{an infinitely deep peak in a particular observable}
}.
$$

This shows the need to distinguish the candidate sequence of Born-type weights from how visibly that sequence appears in a particular experiment.

## 9.4 The Two Observed Principal Peaks

The two principal peaks tracked by the full-range and local high-precision sweeps coincide with the following roots of Theorem 5.1 [4].

| Root | Primitive order $n$ | $m$ | $R_{n,m}=\lvert r\rvert^2$ | $T_{n,m}=\lvert t\rvert^2$ |
|---|---:|---:|---:|---:|
| $R_{124,23}$ | 124 | 23 | 0.697177927556659 | 0.302822072443341 |
| $R_{122,23}$ | 122 | 23 | 0.688363946817593 | 0.311636053182407 |

Both have even primitive order and odd $m=23$. Hence, under the principal conditions of §9.2, the positive and negative extrema occur on the sampling lattice and the gray error vanishes in ideal arithmetic. The sharp peaks observed in the numerical sweeps are finite-precision readouts of this even-order cyclotomic recurrence.

What has been identified here is not the fine-structure constant itself, but finite-order roots of the exchange operator. The theorem in this paper does not explain why these roots lie near physical constants; that remains the separate problem retained in the preceding report [4].

---

# 10. What This Study Establishes

## 10.1 Derived Mathematical Consequences

The following statements are derived rigorously in this work.

1. Up to an overall phase, an exchange-symmetric two-channel unitary operator is represented uniquely as

$$
U=P_s+\zeta P_a,
\qquad |\zeta|=1.
$$

2. The exchange amplitudes in the channel basis are

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2},
$$

so an independent trigonometric assumption is unnecessary.

3. After the overall phase has been fixed, the exact recurrence condition $U^n=I$ is equivalent to the antisymmetric eigenphase being an $n$th root of unity.

4. At a nontrivial finite-order root $\zeta=e^{-2\pi im/n}$, the A/B transition weight from an A-basis input is

$$
W_{A\to A}(n,m)=\cos^2\left(\frac{\pi m}{n}\right).
$$

5. The complementary weight is

$$
W_{A\to B}(n,m)=\sin^2\left(\frac{\pi m}{n}\right).
$$

6. As squared magnitudes of amplitude projections onto the A/B basis, these weights use the same mathematical operation and have the same $\cos^2/\sin^2$ form as the standard two-state Born rule.

7. This discrete sequence was not obtained by inserting the Born rule as the objective function; it follows from exchange symmetry, unitarity, and finite-order recurrence.

8. The theorem concerns the common linear kernel $U$ and does not claim identity with the complete maps containing channel-wise normalization in System A or the C/D backreactions in System B.

## 10.2 Physical Implications Suggested by the Result

The study suggests the following possibility:

> Born-type squared weights can appear as transition weights obtained by projecting the finitely recurrent eigenphases of an exchange-symmetric closed two-channel wave system onto its channel basis.

It further suggests:

> Localization concentration that appears as wave-packet contraction can be represented as a composition of a unitary exchange kernel, redistribution of localization, channel-wise normalization, and a restricted observation section. Two-state selection under strong observation can be represented as a separate dynamics obtained by connecting a nonlinear D backreaction to the same exchange kernel.

## 10.3 What This Study Does Not Yet Establish

This paper does not prove

1. the complete Born rule of standard quantum mechanics,
2. a mechanism generating a probabilistic outcome in a single trial,
3. convergence of repeated-trial frequencies to $\cos^2\phi$,
4. uniqueness of a general projection measure in arbitrary-dimensional Hilbert space,
5. irreversibility in an actual quantum measurement,
6. complete dynamics of the observer, environment, or decoherence,
7. identity between localization transfer in System A and standard wave-function collapse,
8. exact closure of the complete System A map including channel-wise normalization at the same finite order as the linear kernel, or
9. unique determination of the same channel weights for a general $U(2)$ operator without exchange symmetry.

Accordingly, this paper does not state that it has “derived the Born rule.” Its claim is instead:

> Discrete Born-type weights emerge from finite-order recurrence.

---

# 11. Relation to the Standard Born Rule

## 11.1 Born's Statistical Interpretation

In 1926, Born introduced the statistical interpretation that reads event probabilities from the amplitude of the wave function in quantum scattering [5]. In modern notation, for a normalized state $|\psi\rangle$ and projector $P_A$,

$$
P(A)=\langle\psi|P_A|\psi\rangle.
$$

For a pure state and a one-dimensional projection,

$$
P(A)=|\langle A|\psi\rangle|^2.
$$

## 11.2 Difference from Gleason's Theorem

Gleason's theorem shows that, in Hilbert spaces of dimension three or greater, an additive probability measure on orthogonal projections has the density-operator form

$$
\mu(P)=\operatorname{Tr}(\rho P)
$$

[6]. It is a strong result concerning the mathematical uniqueness of Born-type measures.

The standard Gleason theorem does not directly include two-dimensional Hilbert spaces. Busch, however, obtained a generalized Gleason-type result with a density-operator representation by imposing additivity on effects, a broader class than projective measurements [7]. Thus, established mathematical frameworks also characterize Born-type measures for two-state systems.

The present study does not address that level of generality. It concerns an iterated two-channel exchange system and asks

> Why are particular squared projection values selected as exchange weights capable of finite recurrence?

The roles are therefore different:

- Gleason-type results constrain the form of a consistent probability measure.
- This work identifies the discrete squared weights selected by the closure condition of an iterated exchange dynamics.

## 11.3 Difference from Other Derivations of the Born Rule

Zurek's envariance argument aims to derive Schmidt-component probabilities $p_k\propto|\psi_k|^2$ from entanglement symmetry in a composite system [8]. It directly addresses the question of why squared amplitudes should be read as probabilities.

The question here is different. Without taking probability as the starting concept, we classify the eigenphases satisfying

$$
U^n=I
$$

for a single exchange-symmetric two-channel operator and show that its channel-projection weights become $\cos^2/\sin^2$. This paper is therefore not an alternative to envariance or Gleason-type measure theorems, but presents another mathematical route by which Born-type projection weights arise from finite-order recurrence.

## 11.4 Connection to Exact-Recurrence Research

Anand et al. classified exact, initial-state-independent recurrence in finite-dimensional Floquet systems through the cyclotomic structure of unitary spectra [9]. The condition

$$
\zeta^n=1
$$

in this paper is the minimal two-channel realization of such exact recurrence. Whereas that preceding work gives an arithmetic classification of recurrence times for general Floquet spectra, exchange symmetry here fixes the eigenvectors to symmetric and antisymmetric modes. Returning their cyclotomic eigenphases to the A/B channel basis yields

$$
\left|
\frac{1\pm\zeta}{2}
\right|^2
=
\cos^2\phi,\ \sin^2\phi.
$$

This mapping from cyclotomic recurrence to channel squared weights is the particular focus of the present paper.

## 11.5 Significance of Restricting the Analysis to a Two-State System

Because this is a two-state system, the Born rule is not derived by directly applying Gleason's theorem.

Two-state systems nevertheless form the minimal models of

- interference,
- two-level spin,
- two-path measurement,
- two-mode optics, and
- A/B observation selection.

Showing that recurrent exchange weights in this minimal model appear as

$$
\cos^2\phi,
\qquad
\sin^2\phi
$$

is a result with a clearly delimited domain of validity before any extension to a general theory.

---

# 12. Falsifiable Tests

## 12.1 The Exchange-Symmetric Subgroup and General $U(2)$

Theorem 5.1 rules out the possibility that the $\cos^2$ form arises only from the scattering-coefficient representation of §3. What is required is not a particular trigonometric coordinate but the A/B exchange symmetry

$$
[U,X]=0.
$$

A general two-level unitary operator, including its overall phase, can be written as

$$
U
=
e^{i\alpha}
\left(
\cos\beta\,I
-i\sin\beta\,\boldsymbol n\cdot\boldsymbol\sigma
\right),
\qquad
|\boldsymbol n|=1.
$$

The finite-order condition quantizes the eigenphase $\beta$ to a rational angle, but the transition weight from the A basis to the B basis is

$$
|\langle B|U|A\rangle|^2
=
\sin^2\beta\,(n_x^2+n_y^2),
$$

which retains the eigenvector axis $\boldsymbol n$. Finite order alone therefore cannot uniquely fix the general-$U(2)$ A/B weight to $\sin^2\beta$.

The theorem here applies to the centralizer of the exchange operator $X=\sigma_x$: the case in which $\boldsymbol n$ is fixed to the exchange axis and the symmetric and antisymmetric eigenvectors have equal weights relative to the A/B basis. This is not a defect of the theorem but the identification of the precise symmetry condition that produces the discrete Born-type weights.

A falsifiable next test is to introduce a controlled perturbation with $[U,X]\ne0$ and measure whether the peak weights change according to the axis factor $n_x^2+n_y^2$. If they do not, an additional structure not explained by the projection mechanism in this paper must be present.

## 12.2 Testing the Frequency Law

The present result gives closed-orbit weights but does not directly generate measurement frequencies.

The next stage must prepare the same A-basis input $|A\rangle$ repeatedly, pass it through the finite-order kernel and an explicit measurement/selection map, and measure the A/B output frequencies

$$
f_{A\to A},
\qquad
f_{A\to B},
$$

testing whether

$$
f_{A\to A}\rightarrow W_{A\to A}(n,m)=|r_{n,m}|^2,
\qquad
f_{A\to B}\rightarrow W_{A\to B}(n,m)=|t_{n,m}|^2.
$$

If this holds, it will provide a more direct connection between finite-order exchange weights and the Born rule as an experimental frequency law. For a general superposed input, a frequency law containing the interference term in $r\alpha+t\beta$ must be tested separately.

## 12.3 Observational Perturbation and Irreversibility

A fully unitary system is recurrent and is not irreversible in principle.

To connect the model to actual observational selection, one must introduce

- apparatus degrees of freedom,
- environmental degrees of freedom,
- phase noise,
- channel loss, and
- post-readout feedback,

and determine the conditions under which a recurrent orbit becomes effectively fixed in either A or B.

## 12.4 Finite-Order Selection Under Noise

Let the displacement from a root be

$$
R=R_{n,m}+\delta R.
$$

The resonance width of each finite-order root can be defined by measuring the operator residual after $n$ iterations,

$$
\mathcal E_n(R)
:=
\|U_R^n-I\|.
$$

One may also introduce per-iteration phase noise

$$
\theta_j=\theta+\xi_j
$$

and evaluate the recurrence fidelity of

$$
U_{j+n-1}\cdots U_j.
$$

This will distinguish which low-order roots among the mathematically dense set are physically robust.

## 12.5 Normalization Audit of System A

The System A implementation contains channel-wise normalization, so its complete update map is not the linear operator $U_R\otimes I$ but a state-dependent nonlinear map.

It is therefore necessary to compare three conditions:

1. no normalization,
2. joint normalization of the entire two-channel state, and
3. separate normalization of each channel.

This comparison must separate the location of the finite-order roots, localization transfer, and the visibility of Born-type weights.

Such an audit is required to connect rigorously the finite-order structure derived for the common linear kernel to the nonlinear wave-packet localization behavior of System A. It has already been confirmed that the current implementation uses channel-wise normalization. What remains unconfirmed is under which conditions that normalization preserves, shifts, or eliminates the finite-order roots of the linear kernel.

---

# 13. Discussion

## 13.1 Principal Insight

The principal insight of this study is that a seemingly continuous wave system with a continuous exchange coefficient $R$ contains exact discrete closed orbits classified by integer pairs $(n,m)$, and that their weights are fixed directly by exchange symmetry.

Moreover, the exchange weights of these discrete closed orbits are

$$
R_{n,m}=\cos^2\left(\frac{\pi m}{n}\right).
$$

Exchange symmetry and unitarity give

$$
U=P_s+\zeta P_a,
$$

while finite-order closure gives $\zeta=e^{-2\pi im/n}$. Hence the structure

$$
\text{discrete phase closure}
\longrightarrow
\text{amplitude projection}
\longrightarrow
\text{squared weight}
$$

follows. The $\cos^2$ form is not a remnant of an assumed scattering coordinate; it is the squared interference obtained when the sum and difference of two projected eigenphases are returned to the A/B basis.

## 13.2 A Finite-Order Origin of Born-Type Weights

In the Born rule, squared weights are given as observation probabilities.

Here, the same squared weights appear as the exchange ratios required to return to an identical state after finitely many iterations.

The following proposition is therefore a theorem within this model:

> In an exchange-symmetric closed two-channel system, transition weights obtained by projecting recurrent cyclotomic eigenphases onto the A/B channels form a discrete Born-type squared sequence.

In this derivation, probability is not posited initially as an unstructured random number. The order is instead

1. a closed phase relation,
2. an integer ratio permitting finite recurrence,
3. projection onto the observation basis, and
4. squaring of the projected component.

What remains unresolved is the mechanism by which these transition weights govern physical repeated-measurement frequencies, not the derivation of the discrete squared weights themselves.

## 13.3 A Common Kernel for Contraction, Metastability, and Selection

Systems A and B were constructed to investigate different phenomena.

After analysis, however,

- concentration of localization that appears as wave-packet contraction,
- periodic exchange between the two states A and B,
- a gray metastable state,
- preservation under weak readout,
- selection under strong observation, and
- discrete Born-type weights

were all connected to an unperturbed kernel having the same decomposition into symmetric and antisymmetric eigenmodes. The complete update maps are nevertheless not identical. The localization representation in System A adds channel-wise normalization $\mathcal N_A$; weak readout in System B adds $\mathcal C$; and strong-observation selection adds nonlinear backreaction $\mathcal D$.

This means that phenomena discussed as separate problems in quantum measurement can be compared, in the minimal two-channel model, as combinations of a “common exchange-phase kernel” and distinct downstream maps. The important point is that the root location, localization representation, weak readout, and strong selection can be connected in the same computational model without conflating them as a single effect.

## 13.4 Exact Central Conclusion

The exact central conclusion of this paper is:

> The common linear kernel of model families independently constructed to study wave-packet localization transfer, metastable two-state dynamics, weak readout, and strong-observation selection is an exchange-symmetric two-channel unitary operator. Finite-order closure of this kernel necessarily yields the discrete channel-projection weights $\cos^2(\pi m/n)$ and $\sin^2(\pi m/n)$ without specifying the Born rule as the target of the search.

This is not a derivation of the complete Born rule. It does show, by a reproducible model and analytic formulas, that discrete Born-type squared weights arise from exchange symmetry, unitarity, finite-order closure, and channel projection. The frequency law, nonlinear closure of System A, and correspondence between D observation and standard measurement dynamics are separate next problems, not parts of this theorem.

---

# 14. Conclusion

This study analyzed two-channel unitary operators commuting with the channel-exchange operator $X$. By exchange symmetry, after fixing the overall phase the operator can be written as

$$
U=P_s+\zeta P_a,
\qquad |\zeta|=1.
$$

Returning to the channel basis gives the amplitudes

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2}.
$$

The finite-order condition $U^n=I$ quantizes the antisymmetric eigenphase to

$$
\zeta=e^{-2\pi im/n},
$$

and yields the finitely recurrent channel-projection weight

$$
\boxed{
|r_{n,m}|^2=\cos^2\left(\frac{\pi m}{n}\right)
}
$$

and its complement

$$
\boxed{
|t_{n,m}|^2=\sin^2\left(\frac{\pi m}{n}\right)
}.
$$

As squared magnitudes of the amplitudes obtained by projecting an A-basis input onto the A/B basis, this $\cos^2/\sin^2$ form has the same mathematical structure as the standard two-state Born rule in quantum mechanics. The trigonometric form is not an independent assumption: it follows from interference between the exchange-symmetric projectors $P_s,P_a$ and a cyclotomic eigenphase.

The important point is that the model was not designed to reproduce this form. It began with two numerical models intended to study wave-packet localization transfer, contraction-like localization, the white-cat/black-cat/gray-cat metastable interface, weak readout, and strong-observation selection. Subsequent eigenvalue analysis of the observed sharp peaks identified the finite-order roots and discrete Born-type weights, which this paper then elevated to an exchange-symmetric finite-order theorem independent of the trigonometric parametrization.

The roles of the implementation have also been separated. The common linear kernel $U$ determines the root locations; $\mathcal N_A\circ U$, including channel-wise normalization, produces the localization representation of System A; and $\mathcal C\circ U$ and $\mathcal D\circ U$ produce weak readout and strong selection, respectively, in System B.

The conclusion of this paper is therefore not a complete derivation of the Born rule.

What has been established is the precise mathematical connection

$$
\boxed{
\begin{gathered}
\text{exchange-symmetric unitary kernel}
+\text{cyclotomic finite-order closure}\\
+\text{A/B channel projection}
\Rightarrow\text{discrete Born-type squared weights}
\end{gathered}
}.
$$

The next stage will examine general $U(2)$ control experiments that break exchange symmetry, reproduction of measurement frequencies, inclusion of apparatus and environmental degrees of freedom, selection of roots under noise, and an audit of System A under different normalization conditions.

---

# References

## Self-Citations

1. Noriaki Kihara, “Experimental Specification v1 for Low-Localization and Harmonic Transfer Readout in Fermion-Like Collisions Using an Exchange-Interference Scattering Matrix,” 2026. [Japanese].
2. Noriaki Kihara, “Preliminary Experimental Summary v1 of the Acceleration Basis and Localization Exchange in Fermion-Like Collisions Using an Exchange-Interference Scattering Matrix,” Zenodo Concept DOI: 10.5281/zenodo.21333766, 2026. [Japanese].
3. Noriaki Kihara, “Preliminary Experimental Summary v1 of C Weak Readout and D Strong-Observation Selection at a White-Cat, Black-Cat, and Gray-Cat Metastable Interface,” Zenodo Concept DOI: 10.5281/zenodo.21353208, 2026. [Japanese].
4. Noriaki Kihara, “Discovery of Finite-Order Resonances in Iterated Exchange Scattering: Identification of the Origin of Peaks Near the Fine-Structure Values 137 and 128 and a Reproducible Wave-Packet Mathematical Model,” Version DOI: 10.5281/zenodo.21421367, Concept DOI: 10.5281/zenodo.21421366, 2026.

## External References

5. Max Born, “Zur Quantenmechanik der Stoßvorgänge,” *Zeitschrift für Physik* **37**, 863–867 (1926). DOI: 10.1007/BF01397477.
6. Andrew M. Gleason, “Measures on the Closed Subspaces of a Hilbert Space,” *Journal of Mathematics and Mechanics* **6**, 885–893 (1957). DOI: 10.1512/iumj.1957.6.56050.
7. Paul Busch, “Quantum States and Generalized Observables: A Simple Proof of Gleason's Theorem,” *Physical Review Letters* **91**, 120403 (2003). DOI: 10.1103/PhysRevLett.91.120403.
8. Wojciech H. Zurek, “Probabilities from Entanglement: Born's Rule $p_k=|\psi_k|^2$ from Envariance,” *Physical Review A* **71**, 052105 (2005). DOI: 10.1103/PhysRevA.71.052105.
9. Amit Anand, Dinesh Valluri, Jack Davis, and Shohini Ghose, “Quantum Recurrences and the Arithmetic of Floquet Dynamics,” *Quantum* **10**, 2074 (2026). DOI: 10.22331/q-2026-04-20-2074.

---

# Appendix A. Classification of Claims

| Claim | Classification |
|---|---|
| Phase-fixed representation $U=P_s+\zeta P_a$ of an exchange-symmetric unitary operator | Derived consequence of Theorem 5.1 |
| Finite-order condition $U^n=I\Leftrightarrow\zeta^n=1$ | Derived consequence of Theorem 5.1 |
| $r=(1+\zeta)/2$, $t=(1-\zeta)/2$ | Derived from the projector decomposition |
| $\lvert r_{n,m}\rvert^2=\cos^2(\pi m/n)$ | Derived consequence |
| $\lvert t_{n,m}\rvert^2=\sin^2(\pi m/n)$ | Derived consequence |
| Trigonometric form is recovered as an implementation coordinate rather than an independent assumption | Derived consequence |
| Isomorphism between finite-order weights for a basis input and two-state Born projection | Derived comparison |
| The $\cos^2$ sequence was not entered as the search target | Fact about experimental and development history |
| Systems A and B share the same unnormalized, unperturbed linear kernel | Implementation fact confirmed in source code |
| System A uses channel-wise normalization and its complete map is nonlinear | Implementation fact confirmed in source code |
| The complete nonlinear map of System A recurs exactly with the same order as the linear kernel | Untested |
| D selection in System B has an explicit nonlinear backreaction | Implementation fact confirmed in source code |
| The gray error vanishes at even-order roots | Derived consequence restricted to $\varphi_0=0$ or $\pi$, $S_0=C=0.02$, and $2n$ samples |
| Infinitely deep peak | Mathematical divergence of the error metric $-\log_{10}\varepsilon$; not an energy divergence |
| Finite order alone does not determine the same A/B weights for general $U(2)$ | Derived consequence of §12.1 |
| Contraction-like wave-packet behavior can be represented by the exchange kernel, normalization, and partial observation | Mathematical finding within the present System A model |
| Standard quantum-mechanical wave-function collapse has been explained | Not derived |
| The complete Born rule has been derived | Not derived |
| Repeated-measurement frequencies converge to $\cos^2$ | Untested |
| Discrete Born-type channel weights originate from finite-order phase closure | Derived for the exchange-symmetric two-channel kernel |
| The physical Born frequency law originates from finite-order phase closure | Physical hypothesis for future testing |

---

# Appendix B. Minimal Reproduction Formulas

The central result of this paper can be reproduced without assuming a trigonometric representation, using only the following formulas.

## B.1 Exchange-Symmetric Projectors

$$
X=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix},
$$

$$
P_s=\frac{I+X}{2},
\qquad
P_a=\frac{I-X}{2}.
$$

## B.2 Exchange-Symmetric Unitary Operator

If

$$
[U,X]=0,
\qquad
U^\dagger U=I,
$$

then, after fixing the overall phase,

$$
\boxed{
U=P_s+\zeta P_a,
\qquad |\zeta|=1
}
$$

can be written.

## B.3 Finite-Order Condition

For a nontrivial primitive order $n\ge3$,

$$
U^n=I
\quad\Longleftrightarrow\quad
\zeta^n=1.
$$

Thus, with $\gcd(m,n)=1$,

$$
\zeta=e^{-2\pi im/n}.
$$

## B.4 Channel Amplitudes and Squared Weights

$$
U
=
\frac12
\begin{pmatrix}
1+\zeta&1-\zeta\\
1-\zeta&1+\zeta
\end{pmatrix}
=
\begin{pmatrix}
r&t\\
t&r
\end{pmatrix},
$$

and hence

$$
r=\frac{1+\zeta}{2},
\qquad
t=\frac{1-\zeta}{2}.
$$

Their squared magnitudes are

$$
\boxed{
|r_{n,m}|^2
=\cos^2\left(\frac{\pi m}{n}\right)
}
$$

and

$$
\boxed{
|t_{n,m}|^2
=\sin^2\left(\frac{\pi m}{n}\right)
}.
$$

These are the discrete Born-type weights of an exchange-symmetric two-channel process that recurs exactly after finitely many iterations.

## B.5 Recovery of the Numerical Implementation Coordinates

The implementation parametrization

$$
t=e^{i\theta}\cos\theta,
\qquad
r=-ie^{i\theta}\sin\theta
$$

is a coordinate representation of the same operator obtained by setting

$$
\zeta=-e^{2i\theta}.
$$

Indeed, the sum and difference in B.4 give

$$
\frac{1+\zeta}{2}
=
\frac{1-e^{2i\theta}}{2}
=
-ie^{i\theta}\sin\theta=r,
$$

$$
\frac{1-\zeta}{2}
=
\frac{1+e^{2i\theta}}{2}
=
e^{i\theta}\cos\theta=t.
$$

At a finite-order root,

$$
\theta
=
\frac{\pi}{2}-\frac{\pi m}{n}
\pmod\pi,
$$

which reproduces the weights in B.4 directly.

---

# Appendix C. Correspondence with the Core Implementations

The core updates are shown here so that readers can verify both that the finite-order roots were not inserted into the code as hidden constants and that the linear kernel is separated from the downstream maps.

## C.1 System A: Channel-Wise Normalization After Linear Exchange

The iterated update in the [System A implementation](../20260715/run_system_A_localization_exchange_R_sweep_preliminary_v1.py) is

```python
a_next = src.normalize(r * a + t * b)
b_next = src.normalize(t * a + r * b)
a, b = a_next, b_next
```

This explicitly contains the common linear kernel

$$
(a,b)\mapsto(ra+tb,\ ta+rb)
$$

and separate normalization of each output channel. The complete localization-transfer map is therefore nonlinear.

## C.2 System B: Unperturbed Finite-Order Kernel

The iterated update in the [System B direct finite-order test](../20260715/run_minimal_system_B_gray_direct_check_v5.py) is

```python
a, b = normalize_pair(r * a + t * b, t * a + r * b)
```

As long as $U_R$ is unitary and the input pair is normalized, the normalization factor in `normalize_pair` is one. Thus, the iteration that determines the peak locations in this direct test is exactly the linear kernel $U_R$; neither $R_{124,23}$ nor $R_{122,23}$ is injected into the update equation as a target value.

## C.3 System B: Nonlinear Backreaction of Strong D Observation

The [strong-D-observation implementation](../20260714/run_gray_cat_d_observation_response_preliminary_v1.py) adds the following backreaction after exchange:

```python
s_next = s + d_gain * s_d * (1.0 - s * s)
a_next = math.sqrt(0.5 * (1.0 + s_next)) * phase_a
b_next = math.sqrt(0.5 * (1.0 - s_next)) * phase_b
```

The finite-order roots therefore arise from the cyclotomic spectrum of $U_R$, while strong selection into A or B arises from a separate D backreaction. The code structure itself implements the paper's two-layer distinction

$$
\boxed{
\text{generation of recurrence roots}
\ne
\text{state selection under strong observation}
}.
$$
