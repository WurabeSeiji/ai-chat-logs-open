# An Observation on the Structural Correspondences between Signal/Control Theory and Quantum Mechanics, Quantum Optics, and Open Quantum Systems

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd. / Faculty of Engineering Science, Osaka University (graduate)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Version**: v6
**Date**: 3 June 2026
**License**: CC BY 4.0
**Concept DOI**: [10.5281/zenodo.20521598](https://doi.org/10.5281/zenodo.20521598)
**Version DOI (v1.0)**: [10.5281/zenodo.20521599](https://doi.org/10.5281/zenodo.20521599)

---

## Character of This Note

**This is an observation paper. It is neither a proof paper nor an assertion paper.**

This note makes no new physical prediction. It proves no new mathematical theorem. It proposes no new interpretation.

What it does is limited to **juxtaposing and observing that ten structural correspondences already exist in the literature** between the mathematical structures already established and socially implemented in signal/control theory and the mathematical structures established in quantum mechanics, quantum optics, and open quantum systems.

These correspondences are facts already recorded in the standard textbooks and classical original papers of both fields. This note merely organizes them across fields.

This note does **not** claim:

- to propose a new interpretation of quantum mechanics
- to claim a physical extension of signal/control theory
- to claim the superiority of either field
- to predict physical constants
- to propose a modification of existing physical theory
- to assert new mathematical theorems

Evaluation and interpretation are left to the reader.

---

## Abstract

Theoretical physics (quantum mechanics, quantum optics, open quantum systems) and engineering (signal/control theory) developed independently in the 20th century. The former aimed to describe microscopic physical phenomena; the latter aimed at the practical realization of communication, measurement, and control.

Yet the mathematical structures underlying the two fields have, at several important points, an exact identity or a strong structural correspondence. In particular, Heisenberg's quantum-mechanical uncertainty $\Delta x \cdot \Delta p \geq \hbar/2$ (Heisenberg 1927) and Gabor's time–band uncertainty $\Delta t \cdot \Delta \omega \geq 1/2$ (Gabor 1946) are the same mathematical theorem, and Gabor himself made the equivalence explicit.

Starting from this known correspondence, this note organizes the structural correspondences of the two fields into ten items. For each item, the standard literature of both fields is cited side by side, and the range in which an exact identity holds and the range that remains a structural parallel are made explicit.

As noteworthy facts, we observe: (i) the same mathematical object was discovered independently in the two fields (Wigner 1932 in quantum mechanics, Ville 1948 in signal theory, the equivalent phase-space quasi-probability distribution); (ii) the Kennard–Robertson-type Fourier uncertainty and Gabor's time–band uncertainty are based on the same Fourier-analytic inequality, and Gabor himself made the correspondence explicit; (iii) the mathematical structures on the signal/control side are already socially implemented and operating daily in engineering domains such as communications, GPS, autonomous driving, radar, MRI, and optical communications.

---

## §1 Introduction

### 1.1 Independent development of the two fields

In the first half of the 20th century, quantum mechanics was systematized by Heisenberg, Schrödinger, Dirac, von Neumann, and others [1, 11, 16, 17]. In the same period, signal/control theory was systematized by Nyquist, Shannon, Wiener, Gabor, Kalman, and others [2, 5, 6, 7, 8].

The two fields **developed independently with different aims, motivations, and ranges of application**. The former aimed at the prediction and understanding of microscopic physical phenomena; the latter aimed at solving practical problems in communication, measurement, and control.

### 1.2 Correspondences of mathematical structures

Yet, juxtaposing the mathematical structures of the two fields, one can confirm that **an exact identity or a strong structural correspondence** holds for several important concepts. This note organizes that observation.

The most explicit example is the identity between Heisenberg's quantum-mechanical uncertainty principle [1] and Gabor's signal-theoretic uncertainty principle [2]. Gabor (1946) [2] explicitly pointed out in his own paper that the two have the same mathematical structure.

### 1.3 Method of this note

This note adopts the following method:

1. Select pairs of concepts observed to have a correspondence from the two fields.
2. Cite, for each pair, the standard literature of both fields (classical original papers or canonical textbooks).
3. Juxtapose the correspondences of mathematical structure in table form.
4. **Make explicit the range in which an exact identity holds and the range that remains a structural parallel.**
5. Add no interpretation or evaluation.

This note is limited to organizing facts obvious to experts in both fields; its novelty lies only in the "juxtaposition and classification" of structure. It makes no new mathematical or physical claim.

### 1.4 Structure of this note

§2 lists the ten structural correspondences in turn. For each item it distinguishes "the exactly identical part" from "the part remaining a structural parallel." §3 summarizes the observations. §4 makes explicit what this note does not claim. §5 lists implications (these too are observations, not claims). §6 gives the conclusion.

---

## §2 Ten Structural Correspondences

### 2.0 Classification of the correspondences

We classify the following ten items by the strength of correspondence.

| # | Concept pair | Classification |
|---|---|---|
| 1 | Uncertainty principle | **Exact mathematical identity** (same theorem up to unit scale) |
| 2 | Time–frequency quasi-probability distribution | **Exact mathematical identity** (identical defining formula) |
| 3 | Paraxial equation / Schrödinger equation | **Conditional exact isomorphism** (identical as a PDE form under the paraxial, monochromatic, scalar approximation) |
| 4 | Sampling / phase-space degrees of freedom | **Strong structural correspondence** (formal correspondence of effective-DOF counting in a finite domain; asymptotic identity of effective DOF) |
| 5 | State-space representation / Hilbert-space picture | **Structural parallel** (linear time evolution, but differences in unitarity and input/output structure) |
| 6 | Observability / CSCO | **Structural parallel** (parallel formal problems, but mechanistic difference of dynamics/kinematics) |
| 7 | Kalman filter / quantum filtering | **Structural parallel** (corresponds via POVM, distinct from projective measurement) |
| 8 | Jones vector / qubit | **Exact mathematical identity** (same two-component complex vector, $U(2)$ action) |
| 9 | Dephasing / phase noise | **Structural parallel** (shares only the mathematics of the $T_2$ process, not decoherence as a whole) |
| 10 | SVD / Schmidt decomposition | **Exact mathematical identity** (same linear-algebra tool) |

Exact mathematical identity holds for the 5 items #1, #2, #3, #8, #10; strong structural correspondence for the 1 item #4; and structural parallel for the 4 items #5, #6, #7, #9.

### 2.1 Uncertainty principle (exact mathematical identity)

| Signal theory | Quantum mechanics |
|---|---|
| **Gabor time–band uncertainty** | **Heisenberg–Kennard–Robertson uncertainty** |
| $\Delta t \cdot \Delta \omega \geq \dfrac{1}{2}$ | $\Delta x \cdot \Delta p \geq \dfrac{\hbar}{2}$ |
| [Gabor 1946] [2] | [Heisenberg 1927] [1], [Kennard 1927] [22], [Robertson 1929] [23] |

The two are the same inequality up to units (natural units $\hbar=1$). The mathematical proofs of both rest on the Cauchy–Schwarz inequality and the unitarity of the Fourier transform, and are structurally identical. Since the Fourier conjugate of the position-$x$ representation is the momentum $p/\hbar=k$, $\sigma_x\sigma_p\geq\hbar/2$ is the same Fourier-analytic proposition as the Gabor inequality.

**Historical note**: Heisenberg (1927) [1] presented the uncertainty principle as an intuitive, operational argument. It was Kennard (1927) [22] who formulated it as the modern strict variance inequality $\sigma_x\sigma_p\geq\hbar/2$, and Robertson (1929) [23] who generalized it to an arbitrary pair of self-adjoint operators. The "Heisenberg uncertainty" referred to in this note means the modern standard form based on these three papers.

Gabor (1946) [2] explicitly stated in his own work that his inequality is "an analog of Heisenberg's uncertainty principle in wave mechanics." That is, **the equivalence was already stated, as of 1946, by the founder on the signal-theory side himself.**

### 2.2 Time–frequency (phase-space) quasi-probability distribution (exact mathematical identity)

| Signal theory | Quantum mechanics |
|---|---|
| **Wigner–Ville distribution** | **Wigner function** |
| $W(t,\omega) = \int s(t+\tau/2)\,\overline{s(t-\tau/2)}\, e^{-i\omega\tau}\, d\tau$ | $W(x,p) = \dfrac{1}{2\pi\hbar}\int \psi^*(x+y/2)\,\psi(x-y/2)\, e^{ipy/\hbar}\, dy$ |
| [Ville 1948] [4] | [Wigner 1932] [3] |

Both are quasi-probability distributions following the same mathematical definition. **Wigner (quantum mechanics, 1932 [3]) and Ville (signal theory, 1948 [4]) are different people**, yet the mathematical objects they discovered independently have the same structure. That the standard name in both fields, the "**Wigner–Ville distribution**," carries both names reflects the history of independent discovery and mathematical identity.

Both distributions have the same properties:

- their marginals give the true probability distributions (position/momentum or time/frequency)
- they can take negative values (quasi-probability, not satisfying positivity)
- covariance under the linear Fourier transform
- the same autocorrelation structure

Both are used daily, **as the same mathematical object**, in signal processing (radar, sonar, time–frequency analysis) and in quantum optics (quantum tomography).

### 2.3 Schrödinger equation and paraxial wave equation (conditional exact isomorphism: PDE form under the paraxial approximation)

| Signal theory (optics) | Quantum mechanics |
|---|---|
| **Paraxial wave equation** | **Schrödinger equation** |
| $i\dfrac{\partial E}{\partial z} = -\dfrac{1}{2k}\nabla_\perp^2 E - \dfrac{k\,\delta n}{n_0}E$ | $i\hbar \dfrac{\partial \psi}{\partial t} = -\dfrac{\hbar^2}{2m}\nabla^2 \psi + V\psi$ |
| Optical fiber, laser propagation, holography | Single-particle non-relativistic quantum mechanics |
| [Saleh & Teich 1991] [12] | [Schrödinger 1926] [16] |

Under the **paraxial, monochromatic, scalar approximation**, the two become **mathematically isomorphic partial differential equations** under the correspondence of variables $(z \leftrightarrow t,\ 1/k \leftrightarrow \hbar/m,\ V/\hbar \leftrightarrow -\,k\delta n/n_0)$. Here, adopting the standard carrier convention $E\propto e^{ikz}$, the refractive-index term appears with a negative sign, and a high-index region ($\delta n>0$, waveguiding) corresponds to a **bound potential well ($V<0$)** (relative minus sign). Note that a repulsive potential ($V>0$) is also a legitimate Schrödinger form; this sign does not affect the isomorphism classification and is a note for the consistency of the explicit dictionary.

This isomorphism is used daily in the design of optical-fiber communications and integrated optics, and is made explicit in standard optics textbooks including Saleh & Teich (1991) [12], Chapters 2–3. On the optics side it is sometimes taught as "borrowing the mathematics of quantum mechanics for optics," and sometimes the reverse, "borrowing the mathematics of optics for quantum mechanics"; in either direction the mathematical content is the same.

**Caution**: the correspondence is, after all, an approximate, mathematical isomorphism. On the optics side $z$ is a propagation distance; on the quantum side $t$ is time. The physical meaning (probability interpretation, Hilbert-space inner-product structure, measurement theory, etc.) is not identical; the correspondence stays at the mathematical form of a linear Schrödinger-type PDE.

### 2.4 Sampling theorem and phase-space degrees of freedom (strong structural correspondence)

| Signal theory | Quantum mechanics |
|---|---|
| **Shannon–Nyquist sampling theorem / time-bandwidth product** | **Phase-space degree-of-freedom counting** |
| A real signal of band $B$ and observation time $T$ is approximately described by an effective DOF $N\sim 2BT$ (for complex baseband, $\sim BT$ independent **complex** samples) | The number of quantum states for a phase-space volume $V\cdot p_{\max}^d$ is $\sim V\cdot p_{\max}^d/(2\pi\hbar)^d$ |
| [Nyquist 1928] [5], [Shannon 1949] [6], [Slepian–Pollak 1961–1964] [13] | [Planck 1906] [24], [Sackur 1911; Tetrode 1912] [25], standard statistical-mechanics textbooks |

The two have a **strong structural correspondence** based on the common framework of the concept of a "**minimal phase-space cell**." In signal theory it appears as the time-bandwidth product; in quantum mechanics as the phase-space volume divided by the Planck unit $h$.

**The form of the DOF count corresponds strongly**, but the two are not strictly the same theorem. The Shannon–Nyquist / Slepian–Pollak-type time–band DOF and the semiclassical number of states (the quantum number of states by the Weyl law) each have a different measure, boundary condition, and physical interpretation. Gabor (1946) [2] introduced this DOF count on the signal-theory side under the name "**logon**," a counting essentially parallel to the phase-space-cell concept on the quantum side.

**Refinement**: on the signal-theory side, the strict formulation of the effective DOF in a finite-time, finite-band signal space is established by the prolate spheroidal wave functions of Slepian, Pollak, and Landau [13]. **Since a nonzero signal that is both time-limited and band-limited does not strictly exist, the signal space is formally infinite-dimensional, and $N\sim 2BT$ is the asymptotic expression of "the effective dimension where the concentrated (prolate) eigenvalues are close to 1."** On the quantum side, the count dividing the phase-space volume by $h^d$ originates in Planck's quantization of blackbody radiation (1906) [24] and the Sackur–Tetrode equation of statistical mechanics (Sackur 1911; Tetrode 1912) [25], and is generalized as the Weyl law for the semiclassical number of states.

Stated information-theoretically, "**the independent DOF of a finite-band, finite-observation-time system is effectively finite**," which is the same kind of framework as the finitization of the Hilbert-space dimension in quantum mechanics (by phase-space volume). The correspondence of the two is a strong structural correspondence of **the same form of DOF counting, not the same theorem**.

### 2.5 State-space representation and Hilbert-space picture (structural parallel, with important differences)

| Control theory | Quantum mechanics |
|---|---|
| **Linear state equation** | **Schrödinger equation (linear version)** |
| State eq.: $\dot{\mathbf{x}}(t) = A\,\mathbf{x}(t) + B\,\mathbf{u}(t)$ | State evolution: $i\hbar\,\dfrac{d|\psi\rangle}{dt} = \hat{H}\,|\psi\rangle$ |
| Output eq.: $\mathbf{y}(t) = C\,\mathbf{x}(t)$ | Observed value: $\langle A\rangle = \langle\psi|\hat{A}|\psi\rangle$ |
| [Kalman 1960] [8] | [von Neumann 1932] [11] |

Both are parallel in that **a linear operator ($A$ or $\hat{H}$) describes a time evolution on a linear space**. However, the following important differences exist:

**Difference (i): presence/absence of unitarity**. In quantum mechanics $\hat{H}$ is Hermitian (self-adjoint), so $e^{-i\hat{H}t/\hbar}$ is a unitary operator and the norm (probability) of the state is conserved. By contrast, the $A$ of control theory need not be skew-Hermitian (purely imaginary spectrum) in general, and **allows eigenvalues with a real part (unstable poles, dissipative modes)**. That is, the conservative time evolution of a closed quantum system and a control system handling open/dissipative systems correspond only in the special case where $A$ is skew-Hermitian and $B\mathbf{u}=0$.

**Difference (ii): bilinearity of the observed value**. The quantum expectation $\langle A\rangle=\langle\psi|\hat{A}|\psi\rangle$ is **bilinear (a quadratic form)** in the state $|\psi\rangle$, whereas the control-theory output $\mathbf{y}=C\mathbf{x}$ is **linear** in the state $\mathbf{x}$. The way observables are extracted differs structurally.

**Difference (iii): presence/absence of a driving input**. The driving input $B\mathbf{u}(t)$ of control theory is the essence of a controlled system, but the closed Schrödinger equation has no free term corresponding to it. Extending to open quantum systems (Lindblad equation [15]) adds a driving term, but this correspondence is restricted to the standard closed Schrödinger equation.

Hence the correspondence of the two is a structural parallel at the level of "linear time evolution on a linear space," and even restricted to finite-dimensional quantum systems (qubit, qudit), **essential differences remain in unitarity, the way observables are extracted, and the treatment of the driving input**. These differences are made explicit in the literature of both fields and indicate the limits of the structural correspondence.

### 2.6 Observability and state distinguishability (structural parallel, with mechanistic difference)

| Control theory | Quantum mechanics |
|---|---|
| **Observability** | **Complete Set of Commuting Observables (CSCO)** |
| Rank condition of the observability matrix $\mathcal{O} = \begin{pmatrix}C\\CA\\CA^2\\\vdots\\CA^{n-1}\end{pmatrix}$ | The maximal set of mutually commuting self-adjoint operators |
| Whether the state is recoverable from the output time series | Whether the state is specifiable from the simultaneous eigenvalues of observables |
| [Kalman 1960] [8] | [Dirac 1930] [17], [von Neumann 1932] [11] |

The two are **parallel formalizations** of "the condition for specifying the state from observation." But there is an important mechanistic difference:

**Difference (i): dynamics vs kinematics**. Kalman's observability condition is the condition for uniquely determining the initial state $\mathbf{x}(0)$ from the output time series $\{y(t)\}_{t\geq0}$, using the **dynamics** of stacking $CA^k$. By contrast, the quantum CSCO is the **kinematic** condition that a commuting observable algebra at a single time uniquely specifies a basis of the state space, requiring no time evolution.

**Difference (ii): single trajectory vs many copies**. Kalman observability reconstructs the state from a deterministic single trajectory. In quantum mechanics, by contrast, a single measurement generally cannot determine the state, and state determination requires **quantum state tomography** (repeated measurement on many copies of the same state), because quantum measurement is inherently probabilistic and measurement disturbs the state.

The correspondence is a meta-level parallel of "the criterion for whether the observable algebra separates the state space," and the concrete content of the mechanism differs in the two fields.

**Note on the correspondence**: From the viewpoint of state recoverability, the quantum counterpart closest to control-theory observability is not the CSCO but quantum state tomography or an informationally complete POVM. The correspondence with the CSCO made in this note is restricted to the kinematic aspect that observables separate the state space.

### 2.7 Optimal estimation and filtering (structural parallel, correspondence via POVM)

| Control theory | Quantum mechanics |
|---|---|
| **Kalman filter** | **Quantum filtering** |
| Predict → observe → update cycle | State update by generalized measurement (POVM) and quantum conditional expectation |
| Minimum-mean-square-error estimation for a linear Gaussian state-space model | Bayes-optimal estimation of the quantum state under classical output measurement |
| [Kalman 1960] [8] | [Davies & Lewis 1970] [18], [Belavkin 1992] [19] |

Both are parallel in being operations that "**optimally update the state using observed information**." But the center of gravity of the correspondence lies not in **projective measurement (von Neumann–Lüders) but in POVM and quantum filtering**:

**Important distinction**: projective measurement (von Neumann–Lüders) involves an irreversible disturbing back-action on the state, which has no counterpart in the Kalman filter. The Kalman filter is a Bayes update on a classical probability distribution and shares no mathematical structure with the state collapse of projective measurement.

The correct bridge between the two passes through POVM (Positive Operator-Valued Measure, Davies–Lewis 1970 [18]) and quantum filtering theory (Belavkin 1992 [19]). Quantum filtering formulates, in the framework of Bayes-optimal estimation, how the quantum state is updated under continuous classical output measurement (homodyne/heterodyne detection, etc.) (the Belavkin equation, the stochastic master equation). **For linear Gaussian systems the correspondence with a Kalman-type filter is strong, while general quantum filtering remains a structural parallel of Bayesian state update** (generalizing as "the quantum version of the Kalman filter" is too strong).

That is, the correspondence of the Kalman filter with "simple projective measurement" is inaccurate, and the correspondence of the Kalman filter with "POVM + quantum filtering" is the mathematically correct parallel.

### 2.8 Jones vector and qubit (exact mathematical identity)

| Signal theory (optics/communications) | Quantum mechanics (optics) |
|---|---|
| **Jones vector (polarization state)** | **qubit state / Bloch sphere** |
| $|\psi\rangle = \begin{pmatrix} E_x \\ E_y \end{pmatrix}$, $E_x, E_y \in \mathbb{C}$ | $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$, $\alpha, \beta \in \mathbb{C}$, $|\alpha|^2 + |\beta|^2 = 1$ |
| Two-component complex vector (4 real DOF) | Two-component complex vector (4 real DOF; 2 real DOF after normalization + global phase) |
| [Jones 1941] [9] | [Bloch 1946] [14], standard quantum-mechanics textbooks |

Both represent the state by **the exact same two-component complex vector**. Up to the normalization condition and the global-phase freedom, both are described as points on the **complex projective line $\mathbb{CP}^1\cong S^2$**, which has the same geometry as the Bloch sphere (quantum mechanics) / Poincaré sphere (polarization).

Concrete correspondences:

- Orthogonal polarization-basis transformations (H/V ↔ ±45° ↔ R/L circular) ↔ qubit basis transformations ($|0\rangle,|1\rangle$ vs $|\pm\rangle$ vs $|\pm i\rangle$)
- **Lossless Jones matrices** (phase plates, rotators, birefringent elements, and other reversible polarization elements) ↔ $U(2)$ **unitary transformations**
- **General Jones matrices including absorption/polarizers** ↔ **non-unitary linear transformations** (a formal correspondence with quantum operations including measurement/loss on the quantum side)
- Stokes parameters ↔ Bloch-vector components
- The **orthogonal polarization basis** used in polarization-division multiplexing (PDM) ↔ the **same polarization basis** used in quantum optics to describe single-photon polarization qubits and two-photon polarization states (caution: PDM is the dual-polarization multiplexing of classical optical communication, not quantum entanglement itself; what is shared is only the choice of polarization basis)

This is an exact mathematical identity restricted to pure polarization states, and is used daily in the implementation of quantum optics and quantum information. However, as above, the physical interpretation of the general lossy Jones matrix and of PDM has no direct identity with quantum unitarity or quantum entanglement.

**Caution**: this item treats the correspondence between the **classical-optics Jones vector** and the **quantum qubit**. The **I/Q complex baseband signal** $\tilde{s}(t)=I(t)+iQ(t)$ in communications engineering is a single complex scalar (2 real DOF), of different dimension from the Jones vector (two-component complex, 4 real DOF). What corresponds directly to the I/Q baseband signal is the **single-mode complex amplitude** (the single-mode amplitude of a classical electromagnetic field), not a qubit. Care must be taken not to conflate the two.

### 2.9 Dephasing and phase noise (structural parallel, only part of decoherence)

| Signal theory | Quantum mechanics |
|---|---|
| **Phase noise / jitter** | **Dephasing (pure phase relaxation)** |
| $\phi(t) = \omega_0 t + \delta\phi(t)$ | Dissipation of the off-diagonal elements of the density matrix |
| Lorentzian linewidth (characteristic width $\sim 1/T$, convention-dependent) | Transverse relaxation time $T_2$ |
| Oscillator theory, laser-linewidth theory ([Schawlow–Townes 1958] [26]) | [Bloch 1946] [14], [Lindblad 1976] [15] |

Both share the mathematical structure of **dissipation of phase information**. Concretely, it is described as a process in which the phase distribution of the complex amplitude diffuses over time, and the Lorentzian (the Schawlow–Townes linewidth in laser-linewidth theory [26]) and the $T_2$ relaxation of the Bloch equation have the same mathematical form.

**On the linewidth convention**: the coefficient changes by convention — full width at half maximum (FWHM) or half width at half maximum (HWHM), angular frequency $\Delta\omega$ or ordinary frequency $\Delta f$ (e.g., for the Lorentzian spectrum of the exponential correlation $g^{(1)}(\tau)\propto e^{-|\tau|/T}$, FWHM $=1/(\pi T)$ Hz, HWHM $=1/(2\pi T)$ Hz). This note states only the qualitative isomorphism that the characteristic width is of order $\sim 1/T$; the exact coefficient correspondence depends on convention.

**Important restriction**: "decoherence" is a broad concept in quantum information/optics and includes, in addition to (i) pure phase relaxation ($T_2$, dephasing), (ii) energy relaxation ($T_1$, transitions from excited to ground state) and (iii) loss of purity of the reduced density matrix due to entanglement with the environment. What strictly shares a mathematical structure with classical phase noise is **only (i), pure phase relaxation (dephasing)**; (ii) and (iii) have no classical counterpart (since classical systems have no discrete energy levels or quantum entanglement).

Hence the correspondence of this item is restricted not to "phase noise = decoherence as a whole" but to "**phase noise = the dephasing part of decoherence**." Within this restriction, **the corresponding pure-dephasing process on the quantum side is described by a Lindblad-form master equation [15]** (this means the quantum counterpart is expressed in Lindblad form, not that classical phase noise itself obeys a Lindblad equation).

### 2.10 SVD and Schmidt decomposition (exact mathematical identity)

| Signal theory (communications) | Quantum mechanics |
|---|---|
| **Singular value decomposition (SVD)** | **Schmidt decomposition** |
| $H = U \Sigma V^\dagger$, $H \in \mathbb{C}^{N_r \times N_t}$ | $|\Psi\rangle_{AB} = \sum_i \sqrt{\lambda_i} \, |a_i\rangle_A \otimes |b_i\rangle_B$ |
| Decompose an arbitrary matrix by orthogonal basis changes and singular values | Decompose a bipartite pure state by Schmidt coefficients and orthogonal bases |
| [Eckart–Young 1936] [20] | [Schmidt 1907] [21] |

The two are **the same tool of linear algebra**. Representing a bipartite pure state on the tensor-product space $\mathcal{H}_A\otimes\mathcal{H}_B$ as a coefficient matrix, its singular value decomposition is the Schmidt decomposition itself. The Schmidt coefficients $\sqrt{\lambda_i}$ are the singular values of the coefficient matrix, and the Schmidt bases $\{|a_i\rangle\},\{|b_i\rangle\}$ correspond to the left/right orthogonal bases of the SVD. The historical origin of the concept of bipartite entanglement in the quantum case goes back to Einstein–Podolsky–Rosen (1935) [10].

**Limits of the correspondence (important)**: what corresponds in this item is only **the SVD = Schmidt decomposition as a tool of linear algebra**. No further physical/engineering identity holds:

- **The MIMO channel matrix $H$ and the coefficient matrix of a bipartite pure state are different objects**: the MIMO $H$ is the matrix representing the propagation channel, and the signal vectors $\mathbf{x},\mathbf{y}$ are the input/output vectors of a single linear channel. These are not "vectors on a tensor-product space $\mathcal{H}_A\otimes\mathcal{H}_B$ of two subsystems." By contrast, the quantum entangled state $|\Psi\rangle_{AB}$ belongs to the tensor-product space of subsystems.
- **The functional forms of the capacity formulas differ**: the MIMO capacity is $C=\sum_i\log_2(1+\rho\sigma_i^2)$ (a logarithmic sum over singular values $\sigma_i$ with SNR $\rho$), whereas the entanglement entropy is $S=-\sum_i\lambda_i\log\lambda_i$ (the Shannon entropy of the squared, normalized Schmidt coefficients $\lambda_i$). The two have different functional forms, and the isomorphism of capacity does not hold.

Hence this item is limited to the observation that "**the SVD, a tool of linear algebra, plays the same role (giving the canonical form of a bipartite system) in the two fields**." Since the structural isomorphism between MIMO channel capacity and quantum entanglement capacity does not hold, it is not claimed here.

---

## §3 Summary of the Observations

All ten items above are facts made explicit in the standard textbooks and classical original papers of both fields. This note merely juxtaposes them across fields.

The following noteworthy facts can be observed:

**(A) Independent discovery by independent people**

A salient example of the same mathematical object being discovered independently in the two fields is the phase-space quasi-probability distribution. Wigner (1932 [3]) in quantum mechanics and Ville (1948 [4]) in signal theory introduced equivalent distributions independently of each other. They are different people, and the standard name in both fields, the "**Wigner–Ville distribution**," reflects their independent discovery. That the same mathematical object was reached independently from different physical/engineering motivations suggests the mathematical depth of the correspondence.

**(B) The same Fourier-analytic inequality, and Gabor's explicit statement**

The Kennard–Robertson-type Fourier uncertainty and Gabor's time–band uncertainty are based on the same Fourier-analytic inequality. Importantly, Gabor himself made the mathematical identity of the two explicit in his 1946 paper [2]. That is, the correspondence between the two fields **was already recognized by the founders of the two fields themselves.**

**(C) Social implementation in engineering**

The mathematical structures on the signal/control side are already socially implemented and operating daily in the following engineering domains:

- OFDM, MIMO, and I/Q demodulation in mobile communications (4G LTE, 5G NR, 6G research)
- the Kalman filter in GPS and satellite positioning
- state estimation in autonomous driving and drone control
- the Wigner–Ville distribution in radar and sonar
- numerical solution of the Schrödinger-type paraxial wave equation in optical-fiber communications
- the Bloch equation in MRI (magnetic resonance imaging)
- phase-noise analysis in oscillator and laser design

These implementations are a **confirmation of the operation of a shared mathematical structure in engineering**. The empirical validity on the quantum side is independently verified separately by quantum experiments (quantum optics, atom interferometry, superconducting qubits, etc.). The two stay at the observation that the same mathematical structure is confirmed in different physical domains; the social implementation of one does not directly guarantee the empirical validity of the other (this distinction is important to avoid a category confusion).

---

## §4 What This Note Does Not Claim (explicit scope)

This note does not claim:

1. **It does not propose a new interpretation of quantum mechanics**: the observed correspondences are, after all, a juxtaposition of mathematical structures, and express no view on the physical content or interpretation of quantum mechanics (Copenhagen, many-worlds, hidden variables, etc.).
2. **It does not claim a physical extension of signal/control theory**: signal/control theory is positioned purely as practical engineering, and no claim is made to re-position it as a new theory of physics.
3. **It does not claim the superiority of either field**: both fields developed independently and are both mature, so we are not in a position to judge superiority.
4. **It does not predict physical constants**: this note only juxtaposes structure; predicting constant values is out of scope.
5. **It does not propose a modification of existing physical theory**: it proposes no change to the mathematical structures of standard quantum mechanics, quantum optics, or open quantum systems.
6. **It does not assert new mathematical theorems**: all the structural correspondences are recorded in the existing literature of both fields. The contribution of this note lies only in juxtaposition and classification.
7. **It does not extend a strong identity to a structural parallel**: following the classification of §2.0, it distinguishes the range in which an exact mathematical identity holds from the range that remains a structural parallel, and makes no claim to extend a parallel to an identity.

The evaluation of this note is left to the reader.

---

## §5 Implications (as observations)

We record the following as implications that may be naturally drawn from these observations. These are not claims but direct consequences of the observations.

**(I) Consistency of the uncertainty principle**

The Gabor inequality established in signal/control theory and the Heisenberg uncertainty principle of quantum mechanics are the same mathematical structure (§2.1). **As long as one accepts the content as a Fourier-analytic inequality, one cannot treat the mathematical cores of the Gabor-type and Kennard-type uncertainties as different things** (this does not claim an identity of physical meaning).

**(II) Operation of the shared mathematical structure in engineering**

The mathematical structures on the signal/control side are already socially implemented and operating daily in engineering domains such as communications, GPS, autonomous driving, radar, and MRI (§3 (C)). This means at least that part of the mathematical structure shared by the two fields is empirically confirmed in engineering. The empirical validity on the quantum side is confirmed separately by quantum experiments; the two stay at the observation that the same mathematical structure is independently confirmed in different physical domains.

**(III) Possibility of inter-language translation**

By cross-referencing the literature of the two fields, it is possible to understand the concepts of each field in the language of the other (especially for the exact-identity items of §2.0). For example, Heisenberg uncertainty can be translated as the Gabor inequality, a qubit as a Jones vector, and the Schrödinger equation as the paraxial optics equation.

**(IV) Educational and communicative implications**

By making explicit the mathematical structures shared by theoretical physics and engineering, it is observed that communication between experts in the two fields, and bridging in interdisciplinary education, become easier. However, for items classified as structural parallels in §2.0, one must also understand the differences in mechanism and limits.

These implications too are not new claims but direct consequences of the observations of §2.

---

## §6 Conclusion

This note juxtaposed and observed, from the classical literature of both fields, the **ten structural correspondences** already established between signal/control theory and quantum mechanics, quantum optics, and open quantum systems. Of the ten, 5 items (§2.1, §2.2, §2.3, §2.8, §2.10) have an exact mathematical identity within the explicitly stated conditions/conventions. §2.4 is not the same theorem but a **strong structural correspondence** of effective-DOF counting in a finite domain. The remaining 4 items (§2.5, §2.6, §2.7, §2.9) remain structural parallels (see the classification table of §2.0).

These are not new claims but facts already made explicit in the standard textbooks and original papers of both fields. The contribution of this note lies in organizing, as a single list, the correspondences scattered across the literature of both fields, and in distinguishing exact identity from structural parallel.

The structural correspondences of the two fields have been recognized gradually throughout the 20th century (Gabor 1946 [2] making the Heisenberg equivalence explicit, the independent discovery by Wigner 1932 [3] and Ville 1948 [4], the quantum-filtering theory since Belavkin 1992 [19], etc.). This note merely binds them into a single observation.

We hope this note is useful as a bridging resource for researchers and practitioners standing at the boundary of the two fields.

When readers a century hence reach a more comprehensive understanding of the correspondences of the two fields, this note may be referred to as a transitional organization. That is enough.

---

## References

[1] Heisenberg, W. (1927). "Über den anschaulichen Inhalt der quantentheoretischen Kinematik und Mechanik." *Zeitschrift für Physik*, **43**, 172–198.
[2] Gabor, D. (1946). "Theory of Communication." *Journal of the Institution of Electrical Engineers – Part III*, **93** (26), 429–457.
[3] Wigner, E. P. (1932). "On the Quantum Correction for Thermodynamic Equilibrium." *Physical Review*, **40**, 749–759.
[4] Ville, J. (1948). "Théorie et applications de la notion de signal analytique." *Câbles et Transmission*, **2** (1), 61–74.
[5] Nyquist, H. (1928). "Certain topics in telegraph transmission theory." *Transactions of the AIEE*, **47** (2), 617–644.
[6] Shannon, C. E. (1949). "Communication in the presence of noise." *Proceedings of the IRE*, **37** (1), 10–21.
[7] Wiener, N. (1949). *Extrapolation, Interpolation, and Smoothing of Stationary Time Series*. MIT Press, Cambridge, MA.
[8] Kalman, R. E. (1960). "A New Approach to Linear Filtering and Prediction Problems." *Journal of Basic Engineering*, **82** (1), 35–45.
[9] Jones, R. C. (1941). "A new calculus for the treatment of optical systems." *Journal of the Optical Society of America*, **31** (7), 488–493.
[10] Einstein, A., Podolsky, B., Rosen, N. (1935). "Can Quantum-Mechanical Description of Physical Reality Be Considered Complete?" *Physical Review*, **47**, 777–780.
[11] von Neumann, J. (1932). *Mathematische Grundlagen der Quantenmechanik*. Springer-Verlag, Berlin.
[12] Saleh, B. E. A., Teich, M. C. (1991). *Fundamentals of Photonics*. Wiley-Interscience, New York.
[13] Slepian, D., Pollak, H. O. (1961). "Prolate spheroidal wave functions, Fourier analysis and uncertainty — I." *Bell System Technical Journal*, **40** (1), 43–63; Landau, H. J., Pollak, H. O. (1961). "… — II." *Bell System Technical Journal*, **40** (1), 65–84 (Parts III–V appeared in the same series, 1961–1978).
[14] Bloch, F. (1946). "Nuclear Induction." *Physical Review*, **70** (7–8), 460–474.
[15] Lindblad, G. (1976). "On the generators of quantum dynamical semigroups." *Communications in Mathematical Physics*, **48** (2), 119–130.
[16] Schrödinger, E. (1926). "Quantisierung als Eigenwertproblem (Erste Mitteilung)." *Annalen der Physik*, **384** (4), 361–376.
[17] Dirac, P. A. M. (1930). *The Principles of Quantum Mechanics*. Oxford University Press, Oxford.
[18] Davies, E. B., Lewis, J. T. (1970). "An operational approach to quantum probability." *Communications in Mathematical Physics*, **17** (3), 239–260.
[19] Belavkin, V. P. (1992). "Quantum stochastic calculus and quantum nonlinear filtering." *Journal of Multivariate Analysis*, **42** (2), 171–201.
[20] Eckart, C., Young, G. (1936). "The approximation of one matrix by another of lower rank." *Psychometrika*, **1** (3), 211–218.
[21] Schmidt, E. (1907). "Zur Theorie der linearen und nichtlinearen Integralgleichungen." *Mathematische Annalen*, **63** (4), 433–476.
[22] Kennard, E. H. (1927). "Zur Quantenmechanik einfacher Bewegungstypen." *Zeitschrift für Physik*, **44**, 326–352.
[23] Robertson, H. P. (1929). "The Uncertainty Principle." *Physical Review*, **34**, 163–164.
[24] Planck, M. (1906). *Vorlesungen über die Theorie der Wärmestrahlung*. Barth, Leipzig.
[25] Sackur, O. (1911). *Annalen der Physik*, **341** (15), 958–980; Tetrode, H. (1912). *Annalen der Physik*, **343** (7), 434–442.
[26] Schawlow, A. L., Townes, C. H. (1958). "Infrared and Optical Masers." *Physical Review*, **112** (6), 1940–1949.

---

## Appendix A: Summary Table of Correspondences

For the reader's convenience, the ten items of §2 are re-listed. The **Classification** column distinguishes exact identity (**[★ exact]**), strong structural correspondence (**[◎ strong]**), and structural parallel (**[△ parallel]**).

| # | Signal/control theory | Quantum mechanics/optics/open systems | Common mathematical structure | Classification |
|---|---|---|---|---|
| 1 | Gabor time–band uncertainty | Heisenberg–Kennard–Robertson uncertainty | Lower bound of the variance product of a Fourier conjugate pair | **[★ exact]** |
| 2 | Wigner–Ville time–frequency distribution | Wigner phase-space quasi-probability | Fourier autocorrelation of a quadratic form (with conventions/sign conventions) | **[★ exact]** |
| 3 | Paraxial wave equation | Schrödinger equation | First-order-in-time, second-order-in-space linear PDE (identical as a PDE form under the paraxial/monochromatic/scalar approximation) | **[★ conditional exact]** |
| 4 | Shannon–Nyquist sampling | Phase-space DOF counting | Formal correspondence of effective-DOF counting in a finite domain (asymptotic identity of effective DOF) | **[◎ strong]** |
| 5 | State/output equations | Schrödinger equation/expectation | Linear time evolution on a linear space (differences in unitarity, etc.) | **[△ parallel]** |
| 6 | Kalman observability | CSCO | State separation by the observable algebra (dynamics/kinematics difference) | **[△ parallel]** |
| 7 | Kalman filter update | Quantum filtering (via POVM) | Bayes-optimal state estimation under observation | **[△ parallel]** |
| 8 | Jones vector (pure polarization) | qubit state / Bloch sphere | Two-component complex vector algebra ($\mathbb{CP}^1\cong S^2$); non-unitary Jones matrices excluded | **[★ exact]** |
| 9 | Phase noise / jitter | Dephasing ($T_2$) | Dissipation of the phase distribution (Lindblad form on the quantum side, not decoherence as a whole) | **[△ parallel]** |
| 10 | SVD (linear-algebra tool) | Schmidt decomposition (bipartite pure state) | Canonical form on a tensor-product space (capacity formulas not isomorphic) | **[★ exact]** |

**Exact identity**: 5 items (#1, #2, #3, #8, #10)
**Strong structural correspondence**: 1 item (#4)
**Structural parallel**: 4 items (#5, #6, #7, #9)

---

## Appendix B: Position of This Note

This note is a standalone observation paper and has no aim of reinforcing or extending any particular research program or body of papers. All its references are classical original papers or standard textbooks of the two fields, and it has no direct citation relation to other work by the author.

We hope this note functions, for experts, educators, and interdisciplinary researchers in the two fields, as a reference for the mutual translation of concepts.

---

Author: Noriaki Kihara / WF System Co., Ltd. / ORCID [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020) / CC BY 4.0
