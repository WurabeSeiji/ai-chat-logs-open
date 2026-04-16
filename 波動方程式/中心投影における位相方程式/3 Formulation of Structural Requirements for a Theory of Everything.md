# Formulation of Structural Requirements for a Theory of Everything

## — Structural Constraints Based on Definitions of Information and Operations —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19601592](https://doi.org/10.5281/zenodo.19601592)

---

## 0. Scope of Claims (Must Read)

This paper is a thought experiment that formulates the structural requirements a Theory of Everything must satisfy, based on rigorous definitions of information and operations.

### What This Paper Claims

1. To provide rigorous definitions of information $I$ and operation $F^*$ as the most fundamental starting point for a Theory of Everything.
2. To show that causal structure inevitably emerges from the requirement of cognitive accessibility.
3. To enumerate the structural requirements that a Theory of Everything must satisfy, built upon these definitions and the causal structure.
4. To organize the consequences derived from these requirements.

### What This Paper Does Not Claim

1. This paper is not a theory of physics. It makes no claims about the actual universe.
2. This paper does not presuppose any specific geometric construction (such as central projection). The structural requirements are formulated as general demands independent of geometry.
3. This paper makes no empirical predictions.

---

## 1. Definitions

This chapter provides the most fundamental definitions that serve as the starting point for a Theory of Everything. Only **information** and **operations** are defined here; concepts such as causality, time, space, and matter are not presupposed at all.

The definitions in this chapter are constructed as the minimal axiom system that is necessary and sufficient as a foundation for physical theories, referencing the mathematical foundations of information [1] and the framework of computability theory [2].

### 1.1 Definition of Information $I$

**Information $I$** refers to any information that can be expressed numerically. Its format is arbitrary: scalar values, complex values, matrices, tensors, etc. However, information $I$ must satisfy the following properties:

1. **Existence of null state**: A state in which no information exists, $I = \mathrm{null}$, can be defined.
2. **Decidability of null**: For any information $I$, it is decidable whether $I = \mathrm{null}$ or not. The comparison $\mathrm{null} = \mathrm{null}$ is true.
3. **Definition of type**: Each piece of information belongs to a type. For any two pieces of information $I_a,\, I_b$, it is decidable whether they belong to the same type. For example, a real number and an array are different types. $\mathrm{null}$ belongs to every type (i.e., for any information $I$, the type comparison between $\mathrm{null}$ and $I$ always yields the same type).
4. **Decidability of equivalence**: For any two pieces of information $I_a,\, I_b$ belonging to the same type, it is decidable whether $I_a = I_b$ (whether they are identical information).

### 1.2 Definition of Operation $F^*$

The mapping from information to information

$$
\mathrm{OUT} = F^*(\mathrm{IN})
$$

is defined as an **operation**. Operation $F^*$ requires that the input information and output information belong to the same type.

Operation $F^*$ has the following internal structure:

- **Input information $\mathrm{IN}$**: Information provided externally. It can only be referenced and cannot be modified from within $F^*$.
- **Internal information $I$**: Information held internally by $F^*$. It is observable externally as $\mathrm{OUT}$ but cannot be modified from outside.

#### Initial State

The initial state of operation $F^*$ is $I = \mathrm{null}$, $\mathrm{OUT} = \mathrm{null}$.

#### Operation Rules

Operation $F^*$ operates according to the following rules:

1. If $\mathrm{IN} \neq I$, or $I = \mathrm{null}$:
   - Copy input information $\mathrm{IN}$ to internal information $I$ ($I \leftarrow \mathrm{IN}$)
   - Output internal information $I$ to output information $\mathrm{OUT}$
2. If $\mathrm{IN} = I$:
   - Do nothing ($\mathrm{OUT}$ remains unchanged)

### 1.3 Completeness of the Definitions

The above definitions (§1.1–§1.2) constitute the minimal definitions that are necessary and sufficient as a starting point for a Theory of Everything.

These definitions:

- Do not presuppose causality (no temporal ordering between IN and OUT is assumed)
- Do not presuppose time
- Do not presuppose space
- Do not presuppose the number of dimensions
- Do not presuppose any specific physical theory
- Do not presuppose any specific geometry
- Do not restrict the format of information (scalars, vectors, matrices, tensors, etc. are all permitted)

Under these definitions alone, causality does not exist. A stationary state of information (a state without change) is described without contradiction under these definitions.

---

## 2. Cognitive Accessibility and the Emergence of Causality

### 2.1 Requirement of Cognitive Accessibility

Under the definitions of §1 alone, information and operations "exist" but there is no subject to distinguish or identify them. The minimal definition of "cognition" is as follows:

**Cognitive accessibility:** The input information $\mathrm{IN}$ and output information $\mathrm{OUT}$ of operation $F^*$ can be distinguished. That is, an **ordering** exists between $\mathrm{IN}$ and $\mathrm{OUT}$.

### 2.2 Inevitable Emergence of Causality

The moment cognitive accessibility is required, the following inevitably emerge:

1. **Emergence of ordering:** A direction arises in which $\mathrm{IN}$ comes first and $\mathrm{OUT}$ comes after.
2. **Emergence of causality:** The three-stage separation of $\mathrm{IN}$ (cause) → $F^*$ (process) → $\mathrm{OUT}$ (effect) is established.
3. **Emergence of operation count:** A cumulative count $n$ (a non-negative integer) of how many times operation $F^*$ has been executed becomes definable.

The operation count $n$ is a **value local to each instance** of operation $F^*$, not a global counter for the entire system. The absence of a global time parameter is consistent with the "frozen time problem" in general relativity [3].

### 2.3 Permission of Causality-Free States

An important point is that the definitions of §1 also permit states without causality. That is, when cognitive accessibility is **not** required, information and operations can exist without causality.

This corresponds to a state that can exist without possessing causal structure. The existence without causality is also related to discussions of "physics without time" in quantum gravity [4]. Such a state "exists" but is "not cognized," and is described without contradiction under the definitions of §1.

---

## 3. Structural Requirements

### 3.1 Basic Requirement

**The definitions of information $I$ and operation $F^*$ in §1 must be capable of encompassing all information and all operations in current physics.**

**Scope:** The "Theory of Everything" referred to here targets the entire domain of physics that is mathematically analyzable. The humanities and social sciences — sociology, history, philosophy, economics, etc. — are excluded. The structural requirements in this paper aim to encompass only phenomena describable as physical laws through mathematical equations, and make no claims regarding applicability to other domains.

### 3.2 Scope of Detailed Requirements

The scope of detailed requirements derived from the basic requirement (§3.1) is limited to the following:

- Theoretical physics
- Particle physics
- Gravitational field theory
- General relativity

Other branches of physics (thermodynamics, statistical mechanics, condensed matter physics, fluid dynamics, etc.) are excluded from the detailed requirements of this paper as potential future scope extensions.

### 3.3 Supplementary Definitions Based on Scope

The definitions in §1–§2 are general definitions independent of scope. The following supplementary definitions are provided to describe the detailed requirements for the scope of §3.2 (theoretical physics, particle physics, gravitational field theory, general relativity).

#### 3.3.1 Causality

The ordering structure that emerged in §2.2 is defined for the present scope as follows:

**Causality:** An irreversible ordering relation existing between the input information $\mathrm{IN}$ (cause) and the output information $\mathrm{OUT}$ (effect) of operation $F^*$. That is, the direction $\mathrm{IN}$ → $F^*$ → $\mathrm{OUT}$ is established, and the reverse direction is not permitted.

The causality defined here is independent of time. It is also independent of the direction of time (past → future). Time is treated as a type of information $I$ as defined in §1.1, and time is not required for the definition of causality. Causality exists even without the existence of time.

#### 3.3.2 Conserved Quantities (Symmetry)

**Conserved quantity:** A quantity of information that remains invariant before and after operation $F^*$. That is, a function $Q(I)$ of information $I$ whose value does not change regardless of how operation $F^*$ is executed is called a conserved quantity.

By Noether's theorem [5], conserved quantities are equivalent to symmetries. That is, the existence of a conserved quantity $Q$ corresponds one-to-one with the invariance of operation $F^*$ under a certain transformation (symmetry).

#### 3.3.3 Dimension

**Dimension:** The number of independent directions in which information $I$ can vary. That is, the number of independent degrees of freedom in the state space of information $I$ is defined as the dimension.

- **Null dimension:** A state in which the information of dimension itself does not exist ($d = \mathrm{null}$). This corresponds to the null state of §1.1. Since the concept of dimension itself is not defined, nothing can be said about spatial properties.
- **Zero dimension ($d = 0$):** A state in which the dimension exists but its value is 0. Information $I$ has no independent direction in which it can vary. However, information itself can exist. In zero dimensions, the distinction between "exists" ($I \neq \mathrm{null}$) and "does not exist" ($I = \mathrm{null}$) can be made. The count or magnitude of information is not defined, but whether it exists or not is decidable. It should be noted that the definition of information $I$ in §1.1 places no restriction on the size or complexity of its content. A single piece of information encoding the state of the entire universe is also information $I$. Therefore, even in zero dimensions, the one piece of information that "exists" can have arbitrarily complex content.
- **$d$ dimensions ($d \geq 1$):** A state in which information $I$ can vary in $d$ independent directions.
- **Infinite dimensions:** A state in which there is no upper bound on the number of independent directions in which information $I$ can vary.

The number of dimensions is determined by the structure of conserved quantities (symmetries).

### 3.4 Detailed Requirements

#### Requirement 1: Existence of Causality

**Requirement:** The theory must possess causal structure. That is, an ordering must exist between cause and effect.

**Rationale:** In all of theoretical physics, particle physics, gravitational field theory, and general relativity, causality is the most fundamental premise. No physical theory exists without causality.

**Sufficiency:** By §2, causal structure inevitably emerges at the point when cognitive accessibility is required. Therefore, this requirement is already satisfied by the definitions of §1 and the requirement of cognitive accessibility in §2.

#### Requirement 2: Conserved Quantities Must Be Definable

**Requirement:** The theory must be able to define quantities (conserved quantities) that remain invariant before and after operation $F^*$.

**Rationale:** In every branch of physics, conserved quantities are the most fundamental structures. Physical theories cannot exist without conserved quantities such as energy conservation, momentum conservation, angular momentum conservation, and charge conservation. By Noether's theorem [5], conserved quantities are equivalent to symmetries. Symmetries generate conserved quantities, and the structure of conserved quantities determines the dimensions of space. Therefore, the definability of conserved quantities precedes the existence of dimensions.

**Sufficiency:** Operation $F^*$ of §1 has the rule that when $\mathrm{IN} = I$, it does nothing ($\mathrm{OUT}$ remains unchanged). This means that when the input matches the internal state, the output is conserved. Conserved quantities can be defined as combinations of information satisfying this invariance condition.

#### Requirement 3: Dimensions Must Be Definable

**Requirement:** The theory must possess the concept of dimension as defined in §3.3.3. That is, it must be possible to define the dimension as the number of independent degrees of freedom in the state space of information $I$.

**Scope dependency:** This requirement depends on the scope of §3.2. The general definitions of §1–§2 do not presuppose the existence of dimensions, but within the scope of this paper (theoretical physics, particle physics, gravitational field theory, general relativity), all fields are formulated upon dimensional structure. Theories that do not require dimensions are outside the scope of this paper. Therefore, the definability of dimensions is a mandatory requirement within the present scope.

**Sufficiency:** By §3.3.3, dimension is definable as the number of independent directions in which information $I$ can vary. Null dimension, zero dimension, $d$ dimensions, and infinite dimensions are each clearly distinguished.

**Extended requirements:** For a Theory of Everything, the following requirements — which are not mandatory within the current scope — are defined as necessary:

- **Infinite extensibility of dimensions:** Current theories presuppose a specific number of dimensions (such as 3+1). However, a Theory of Everything is required to be non-degenerate and regular for any number of dimensions from zero to infinity.
- **Symmetry of dimensions:** Current theories assign specific meanings to particular axes (such as the distinction between the time axis and spatial axes). However, a Theory of Everything requires that all dimensions are equivalent (that dimensions possess symmetry). Assigning a special role to a particular axis must be derived as a result of spontaneous symmetry breaking and must not be included in the definition.
- **Indifference to continuous/discrete:** Current theories presuppose that dimensions are continuous. However, a Theory of Everything is required to be non-degenerate and regular whether dimensions take continuous or discrete values.

#### Requirement 4: Reversibility of Causality

**Requirement:** The theory must be non-degenerate and regular even when the direction of causality is reversed. That is, even when the direction $\mathrm{IN}$ (cause) → $F^*$ → $\mathrm{OUT}$ (effect) is reinterpreted as $\mathrm{OUT}$ (cause) → $F^*$ → $\mathrm{IN}$ (effect), the structure of the theory must be preserved.

**Scope dependency:** This requirement depends on the scope of §3.2. In all four fields of the present scope (theoretical physics, particle physics, gravitational field theory, general relativity), the fundamental equations are all symmetric under reversal of causality. Theories that are asymmetric under reversal of causality (such as thermodynamics) are outside the present scope.

**Distinction between reversal of causality and time reversal:** What is required here is strictly the **reversal of causality**, which is not synonymous with time reversal. As defined in §3.3.1, causality is defined independently of time. Whether time reversal is required when causality is reversed depends on the specific construction of operation $F^*$ in individual theories, and is not included in the structural requirements of this paper.

**Sufficiency:** The definition of operation $F^*$ in §1.2 does not prohibit the exchange of the roles of input and output. The reversal of causality corresponds to reversing the direction while preserving the cognitive accessibility (existence of ordering) of §2.1.

#### Requirement 5: Inverse Operations Must Be Definable

**Requirement:** The theory must be able to define an inverse operation that undoes the result of operation $F^*$.

**Scope dependency:** This requirement depends on the scope of §3.2. In theoretical physics, the reversibility of unitary transformations; in particle physics, CPT symmetry; and in gravitational field theory and general relativity, the reversibility of field equations — each requires the existence of inverse operations.

**Sufficiency:** Operation $F^*$ of §1.2 is defined as a mapping from information to information and does not prohibit inverse operations. An inverse operation can be defined as a separate operation from $F^*$.

#### Requirement 6: Composition of Operations Must Be Possible

**Requirement:** The theory must be able to chain and compose multiple operations $F^*$. That is, it must be possible to use the output of one operation as the input to another operation.

**Scope dependency:** This requirement depends on the scope of §3.2. In theoretical physics, the group structure of symmetries; in particle physics, field interactions and perturbation expansions; in gravitational field theory, nonlinear field equations; and in general relativity, the composition of coordinate transformations — each requires the composition of operations.

**Sufficiency:** Operation $F^*$ of §1.2 requires that input information and output information belong to the same type. Therefore, the output of one operation can be passed as input to another operation of the same type, and the composition of operations is guaranteed by definition.

#### Requirement 7: Finiteness of Values Must Be Guaranteeable

**Requirement:** When the values of information $I$ can become infinite, the theory must have a method to map them to a finite range. That is, it must be possible to structurally avoid the divergence of the values of information $I$ through repeated application of operation $F^*$.

**Scope dependency:** This requirement depends on the scope of §3.2. In theoretical physics and particle physics, ultraviolet and infrared divergences in quantum field theory are fundamental problems addressed by renormalization theory [6]. In gravitational field theory and general relativity, divergences at singularities (such as the center of black holes) remain unsolved problems [7]. All of these stem from the fact that finiteness of values is not structurally guaranteed. A Theory of Everything requires the existence of a mapping that guarantees finiteness of values at the structural level, rather than relying on individual techniques (such as renormalization).

**Sufficiency:** The definition of information $I$ in §1.1 does not restrict the range of values. Therefore, there is room to construct a mapping to a finite range. However, the specific method of mapping (embedding into phase space, compactification, etc.) is outside the scope of the structural requirements of this paper and is left to the construction of individual theories.

**Regularity of the inverse mapping:** Furthermore, the inverse mapping of this finitization mapping must also be regular. That is, the original information must be recoverable from the information mapped to a finite range. This is already required by Requirement 5 (existence of inverse operations) and the definition of operation $F^*$ in §1.2 (input information and output information belong to the same type), but is restated here in the context of finiteness of values. If a finitization mapping irreversibly loses information, it does not satisfy the requirements of operation $F^*$.

#### Requirement 8: Dimensional Unification

**Requirement:** The dimensions of all physical quantities appearing in the theory must be expressible as powers $L^n$ of a single fundamental dimension $L$ (length).

**Rationale:** In current physics, four independent dimensions are used: mass $M$, length $L$, time $T$, and charge $Q$. However, upon adopting natural units ($c = 1$, $\hbar = 1$):

$$[T] = [L], \quad [M] = [L^{-1}], \quad [E] = [L^{-1}]$$

The dimensions of major physical quantities are expressed in terms of $L^n$ as follows:

| Physical quantity | SI dimension | $L^n$ expression |
|:--|:--|:--|
| Length $L$ | $L$ | $L^1$ |
| Time $T$ | $T$ | $L^1$ ($c = 1$) |
| Mass $M$ | $M$ | $L^{-1}$ ($\hbar = 1$) |
| Energy $E$ | $ML^2T^{-2}$ | $L^{-1}$ |
| Velocity $v$ | $LT^{-1}$ | $L^0$ (dimensionless) |
| Gravitational constant $G$ | $M^{-1}L^3T^{-2}$ | $L^2$ |
| Planck constant $\hbar$ | $ML^2T^{-1}$ | $L^0$ (by definition, $\hbar = 1$) |
| Speed of light $c$ | $LT^{-1}$ | $L^0$ (by definition, $c = 1$) |

**Consequence for the charge dimension:** Under this requirement, charge $Q$ must also be expressible as $L^n$. In Gaussian units, $[Q] = [M^{1/2}L^{3/2}T^{-1}]$, and under $c = \hbar = 1$, this yields $[Q] = [L^0]$ (dimensionless). That is, the elementary charge $e$ becomes a dimensionless pure numerical constant, and the fine-structure constant $\alpha = e^2/(4\pi)$ also becomes a pure geometric constant.

**Verification:** The major equations of gravitational field theory, electromagnetism, quantum mechanics, and thermodynamics all satisfy this requirement.

**Sufficiency:** Since the definition of information $I$ in §1.1 does not constrain the structure of dimensions, a construction that unifies all dimensions into $L^n$ is permitted.

#### Requirement 9: Scale Symmetry

**Requirement:** All equations appearing in the theory must be scale-homogeneous under the scale transformation $L \to \lambda L$ ($\lambda > 0$). That is, when each physical quantity has dimension $L^n$, it transforms as $L^n \to \lambda^n L^n$ under $L \to \lambda L$, and both sides of any equation must carry the same power $\lambda^k$.

**Rationale:** For a Theory of Everything to depend on a specific scale (a specific length) would mean that another theory is needed to explain the origin of that scale. A truly fundamental theory must be scale-homogeneous (power-law). Logarithmic scale dependence (e.g., the logarithmic running of the coupling constant $\alpha_s(\mu)$ in quantum chromodynamics) does not satisfy this requirement.

**Formulation:** By Requirement 8, all physical quantities have dimension $L^n$. Under the scale transformation $L \to \lambda L$:

$$f(\lambda L) = \lambda^k f(L)$$

is required to hold. This means that $f$ is a homogeneous function of degree $k$ in $L$.

**Verification results:** This requirement was verified for the following major equations:

- **Classical mechanics:** Newton's equation of motion $F = ma$, gravitational force $F = GMm/r^2$, Kepler's third law $T^2 \propto a^3$ — all scale-homogeneous ✓
- **Gravitational field theory:** Einstein equation $G_{\mu\nu} = 8\pi G T_{\mu\nu}$, Schwarzschild solution $r_s = 2GM/c^2$, geodesic equation — all scale-homogeneous ✓
- **Electromagnetism:** Maxwell's equations, Coulomb force $F = e^2/(4\pi r^2)$, Lorentz force — all scale-homogeneous ✓
- **Quantum mechanics:** Schrödinger equation, Dirac equation, Klein–Gordon equation — all scale-homogeneous ✓
- **Thermodynamics:** Stefan–Boltzmann law $P \propto T^4$, blackbody radiation spectrum — all scale-homogeneous ✓

**Exception:** The coupling constant $\alpha_s(\mu)$ of quantum chromodynamics (QCD) exhibits logarithmic running and does not satisfy this requirement. This suggests that QCD is an effective theory incorporating renormalization group effects and requires derivation from a more fundamental theory.

**Sufficiency:** Since the definitions of §1 do not prescribe any specific behavior under scale transformations, the construction of a scale-homogeneous theory is permitted.

---

## 4. Remarks: Examples of Constructions That Facilitate Satisfaction of Structural Requirements

This chapter does not constitute mandatory requirements. It lists techniques believed to facilitate the satisfaction of the structural requirements of §3 when constructing a theory.

### 4.1 Central Projection

Using central projection as a mapping from coordinate space to phase space facilitates the satisfaction of Requirement 7 (finiteness of values). Central projection is a method that maps a coordinate space including infinity to a finite phase space, and the inverse mapping is also regular.

### 4.2 Treating Coordinate Information as Phase Information

When coordinate information is treated as phase information (angles), the domain of values becomes bounded (e.g., $0$ to $2\pi$), making it easier to structurally guarantee Requirement 7 (finiteness of values).

---

## References

[1] C. E. Shannon, "A Mathematical Theory of Communication," *Bell System Technical Journal*, vol. 27, pp. 379–423, 623–656, 1948. — Mathematical foundations of the definition of information.

[2] A. M. Turing, "On Computable Numbers, with an Application to the Entscheidungsproblem," *Proceedings of the London Mathematical Society*, ser. 2, vol. 42, pp. 230–265, 1937. — Formal definition of computability and operations.

[3] C. J. Isham, "Canonical Quantum Gravity and the Problem of Time," in *Integrable Systems, Quantum Groups, and Quantum Field Theories*, NATO ASI Series, vol. 409, pp. 157–287, 1993. — The frozen time problem in quantum gravity; the absence of a global time parameter.

[4] C. Rovelli, "Forget Time," *Foundations of Physics*, vol. 41, pp. 1475–1490, 2011. arXiv:0903.3832. — The possibility of "timeless" formulations in fundamental physics.

[5] E. Noether, "Invariante Variationsprobleme," *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, pp. 235–257, 1918. — Equivalence of symmetries and conserved quantities (Noether's theorem).

[6] S. Weinberg, "Living with Infinities," arXiv:0903.0568, 2009. — Overview of divergence problems and renormalization in quantum field theory.

[7] R. Penrose, "Gravitational Collapse and Space-Time Singularities," *Physical Review Letters*, vol. 14, pp. 57–59, 1965. — Singularity theorems in general relativity.
