# Thought Experiment (8): Are Physical Quantities Real Numbers? — Physical Quantities as Complex-Valued, Phase-Valued, and Discrete-Valued Structures

## The Question of This Paper

This paper is a thought experiment that distinguishes between the fact that observed physical quantities are real-valued and the assumption that the underlying entities of physical quantities are themselves real-valued. By reading physical quantities as structures on a phase space carrying complex values, phase values, and discrete values, we examine whether uncertainty, interference, entanglement, and quantization conditions can be reinterpreted from a single common viewpoint.

---

[Kihara]

No one has tackled this head-on, but I think the underlying entity of length is a complex number, and the imaginary part is the essence of uncertainty. Has anyone proposed this interpretation?

---

[Iris]

This question is important. I will check separately the claim "length is a complex number" itself, the claim that "the imaginary part is the essence of uncertainty," and the existing research on "complex distances / complex metrics / quantum geometry." There are some concepts close to your claim, but I will be careful about whether they coincide.

Yes, there are nearby interpretations.
However, the very form in which Kihara-san states it —

> the entity of length is a complex number, and the imaginary part is the essence of uncertainty

— stated this bluntly as a standard mainstream interpretation, is, at least as far as I can confirm, not common.
However, there is **a series of nearby existing theories**.

The closest is the following.

## 1. Quantum Geometric Tensor: Real Part Is Distance, Imaginary Part Is Berry Curvature

In the parameter space of quantum states, there is a complex tensor called the **quantum geometric tensor**. It splits into a real part and an imaginary part.

$$
Q_{ij}
=
g_{ij}
+
\frac{i}{2}\Omega_{ij}
$$

Here,

- Real part $g_{ij}$: the distance between quantum states, i.e. the quantum metric / Fubini–Study metric
- Imaginary part $\Omega_{ij}$: Berry curvature, i.e. symplectic structure, non-commutativity, topological curvature

In their 2023 paper, Hetényi and Lévay state explicitly that the complete quantum metric has a real part (the Provost–Vallée metric) and an imaginary part (Berry curvature), and they further state that the determinant of this complex quantum metric yields a multi-operator version of the uncertainty principle. This is close to Kihara-san's intuition. ([arxiv.org](https://arxiv.org/abs/2309.03621))

In other words, in existing theory, rather than saying

> length itself is simply a complex number

it would be put as

> the "tensor that measures distance" in quantum state space has a complex Hermitian structure; its real part is the distance, and its imaginary part is the symplectic curvature / Berry curvature, related to uncertainty and non-commutativity.

This is quite close.

---

## 2. In Geometric Quantum Mechanics, the Quantum State Space Is a Kähler Manifold

In Ashtekar and Schilling's geometric quantum mechanics, quantum states are treated as rays in a complex Hilbert space, i.e. points in the projective Hilbert space. This space is naturally a **Kähler manifold**. A Kähler manifold is, roughly speaking, a geometry that integrates

- Riemannian metric: distance
- symplectic form: canonical conjugacy / non-commutativity
- complex structure: complex structure

into a single object.

Ashtekar and Schilling state that quantum-specific uncertainty and the collapse of the state vector can be formulated geometrically, and that in this formulation a Riemannian metric absent in classical mechanics plays a role. ([arxiv.org](https://arxiv.org/abs/gr-qc/9706069))

What is important here is that the body of quantum mechanics appears not as

> a geometry of only real distance

but as

> distance structure + symplectic structure + complex structure.

So, putting it in Kihara-san's terms,

> the "length" in the quantum sense is not closed within real distance alone but is accompanied by a symplectic component corresponding to the imaginary direction.

That much can be said.

However, when speaking to specialists, rather than saying "length is a complex number," it is more standard to say
**"the real part of the Hermitian metric on the quantum state space is the distance, the imaginary part is the symplectic form, and this imaginary part is involved in non-commutativity and uncertainty."**

---

## 3. "The Imaginary Part Is the Essence of Uncertainty" Is Close, but More Precisely "Non-Commutativity / Topological Curvature / Symplectic Structure"

Whether the imaginary part is directly "uncertainty itself" requires some caution.

In the standard formulation, the uncertainty relation has Robertson–Schrödinger form:

$$
(\Delta A)^2(\Delta B)^2
\geq
\left(
\frac{1}{2}\langle \{ \Delta A,\Delta B \}\rangle
\right)^2
+
\left(
\frac{1}{2i}\langle [A,B]\rangle
\right)^2
$$

Here,

- The anti-commutator part: covariance, real symmetric component.
- The commutator part: non-commutativity, symplectic component, the part bound to the imaginary unit $i$.

So, more precisely,

> One essential half of the lower bound for uncertainty is the commutator structure carrying $i$, i.e. the symplectic structure.

Kihara-san's intuition "the imaginary part is the essence of uncertainty" **is correct, if it refers to the commutator / Berry curvature / symplectic form side**.

---

## 4. However, "Length Itself Is a Complex Number" Is Not the Standard Expression

A distinction is needed here.

In ordinary quantum geometry,

$$
ds^2
$$

itself is most often measured with a real Riemannian metric. That is, the physical distance is a real number.

However, the underlying Hermitian tensor is complex, and the inner product

$$
\langle d\psi | d\psi \rangle
$$

carries a complex structure. From it, one obtains

- the real part: distance, Fubini–Study metric
- the imaginary part: Berry curvature / symplectic form.

So, in standard-theory-aligned language, rather than

> length is a complex number,

we should say

> the inner product / metric tensor generating length has a complex Hermitian structure, whose imaginary part represents the origin of quantum uncertainty / non-commutativity.

This phrasing connects to existing research.

---

## 5. If Adding This to Kihara-san's Observation Paper, the Following Wording Is Safe

If adding it to the current paper, a section like the following is good.

> Pushing the intuition of this paper further, the "length" in quantum state space can be understood not as a mere real-valued distance, but as a complex geometry induced from the Hermitian structure of complex Hilbert space. In the quantum geometric tensor, the real part corresponds to a Fubini–Study-type quantum metric, while the imaginary part corresponds to Berry curvature, i.e. a symplectic form. Hence, the geometric essence of uncertainty can be interpreted as residing not only in distance as a real part but also in the non-commutativity / symplectic structure appearing as the imaginary part. However, this is not a claim that the physical length itself takes complex values directly; rather, it means that the underlying quantum geometric tensor defining length carries a complex structure.

---

## 6. Table of Correspondence with Existing Research

| Kihara-san's intuition | Closest in existing theory |
|---|---|
| The entity of length might be complex | Hermitian metric / quantum geometric tensor |
| The real part is ordinary length | Quantum metric / Fubini–Study metric |
| The imaginary part is the essence of uncertainty | Berry curvature / symplectic form / commutator |
| Uncertainty is an area constraint | Robertson–Schrödinger inequality, quantum blob |
| Length, phase, and uncertainty are unified | Kähler geometry of projective Hilbert space |

The closest keywords are

**quantum geometric tensor**
**complex quantum metric**
**Berry curvature**
**Fubini–Study metric**
**Kähler structure**
**geometric quantum mechanics**
**Robertson–Schrödinger uncertainty relation**

---

## Final Conclusion

It is not that no one has considered anything like this.
Rather, close ideas already exist.

However, the mainstream phrasing is, instead of

> length is a complex number, the imaginary part is uncertainty,

rather

> in the complex Hermitian metric of quantum state space, the real part gives the distance, the imaginary part gives the Berry curvature / symplectic form, and this imaginary part is deeply involved in non-commutativity and the uncertainty relation.

Therefore, Kihara-san's intuition **hits an important point**.
However, to convince specialists, the strongest expression is not "complex length" but **the uncertainty structure as the imaginary part of the complex quantum metric tensor**.

---

[Kihara]

Hmm, is that really so? The entity of length being a real number is also just because what is observed is a real number; no one has actually observed the entity, and treating the entity as a real number is just a definitional assumption, isn't it?
If it cannot be observed, then assuming it as the entity is also the same thing, isn't it?

---

[Iris]

Yes, that counter-argument is correct.
More precisely, **"the entity of length is a real number" is not itself an observational fact but a theoretical construction-assumption.**

What is observed is always

$$
L_{\mathrm{obs}} \in \mathbb{R}
$$

But from this it does not immediately follow that

$$
L_{\mathrm{real}} \in \mathbb{R}.
$$

This is important.

## 1. The Observed Value Being Real and the Underlying Quantity Being Real Are Different

In standard physics, observables such as length, time, mass, and momentum are defined so as to take real values. This is consistent with experimental results and makes calculations stable.

But that is, at most, the convention

> we model observable quantities as real-valued functions.

The fact that observed values are real does not require us to assert that the underlying generating structure is also real.

This is exactly so in quantum theory.
Wave functions are complex, but observation probabilities are real:

$$
\psi \in \mathbb{C}, \quad |\psi|^2 \in \mathbb{R}.
$$

That is, the idea

> behind something observed as a real number lies a complex structure

is in fact perfectly normal in quantum theory.

So, organizing Kihara-san's claim, it becomes

> length also appears as a real number at the moment of observation, but the underlying geometric entity may carry a complex structure.

This is by no means a meaningless question.

---

## 2. However, It Is Not "the Same Thing"

What must be noted here is that

> if it cannot be observed, then assuming it real or assuming it complex is the same

is not something we can say.

The reason is that even for an unobservable entity-assumption, if the **constraints, correlations, transformation rules, and conservation laws** derived from it differ, then it has physical meaning.

For example, suppose one assumes a complex length

$$
L = a + ib.
$$

Whether the observed length is

$$
L_{\mathrm{obs}} = \mathrm{Re}(L)
$$

or

$$
L_{\mathrm{obs}} = |L|
$$

or

$$
L_{\mathrm{obs}} = L^* L
$$

leads to different predictions.

And what the imaginary part $b$ represents —

- an uncertainty width
- a phase
- a spread in the conjugate direction
- a symplectic capacity
- an unobservable internal degree of freedom

— affects the meaning of the theory.

So, to be convincing to specialists, just saying

> the entity is a complex length

is not enough.

What is needed is the form

> assuming a complex length $L = x + i\sigma$, the relation between the observed length $x$ and the uncertainty width $\sigma$ agrees with the uncertainty relation / quantum metric / Berry curvature / symplectic capacity of standard quantum theory.

Once stated at this level, it becomes a theoretical model rather than mere metaphysics.

---

## 3. Putting Kihara-san's Intuition Into Formulas Yields a Natural Form

For instance, set the entity of length as

$$
\mathcal{L} = L + i\Delta L.
$$

Here,

- $L$: the observed real-valued length,
- $\Delta L$: the indistinguishability / spread / uncertainty associated with that length.

Then, what is observed is always the real part:

$$
L_{\mathrm{obs}} = \mathrm{Re}(\mathcal{L}).
$$

However, as a physical state it also carries the imaginary part:

$$
\mathrm{Im}(\mathcal{L}) = \Delta L.
$$

In this view, a "point" is really not a single point on the real axis but rather an object with width

$$
x + i\Delta x.
$$

That is because uncertainty enters not as a measurement error appended after the fact, but as **the imaginary component of the length itself**.

In other words,

> the real-valued length is the component projected by observation,
> while the imaginary length is the component that is not observed but acts as uncertainty.

This view is highly compatible with Kihara-san's "projection-of-observables" paper.

---

## 4. To Bring It Closer to Standard Theory, "Complex Length" Should Be "Complex Distance Structure"

However, telling specialists outright

> length is a complex number

provokes pushback.

A more acceptable phrasing is

> the observed distance is real-valued, but the distance structure that describes the underlying quantum geometry can be regarded as having a complex Hermitian structure.

That is, rather than denying

$$
\text{real distance},
$$

set

$$
\text{real distance} = \text{the real-part-observed quantity of a complex geometry}.
$$

This connects to the existing quantum geometric tensor:

$$
Q_{ij}=g_{ij}+i\Omega_{ij}.
$$

Here,

- $g_{ij}$: real part, distance / quantum metric,
- $\Omega_{ij}$: imaginary part, Berry curvature / symplectic form / non-commutativity.

Aligning Kihara-san's hypothesis with this form yields

> the real-valuedness of length is a property after observational projection; the underlying quantum geometry has a complex structure. The imaginary part is not an unobservable surplus component but rather casts an observable shadow as non-commutativity, phase, and uncertainty constraints.

This expression is strong.

---

## 5. "If It Cannot Be Observed, It Is the Same" Risks Becoming Metaphysics Rather Than Physics

As Kihara-san says,

> the entity is real-valued

is also an assumption.

But

> the entity is complex-valued

is also an assumption.

What makes a physical difference is what is derived from those assumptions.

For instance, the following differences turn it into physics:

1. The Robertson inequality emerges naturally from the imaginary part of the complex length.
2. A classical limit in which the imaginary part vanishes can be defined.
3. The magnitude of the imaginary part corresponds to a wave-packet width or symplectic capacity.
4. The observed length is obtained from the product or norm of the complex length.
5. The phase rotation of a complex length corresponds to a Fourier transform or a transformation of conjugate quantities.

Once one reaches this level, instead of

> the real-length hypothesis is an assumption and the complex-length hypothesis is also an assumption,

we can say

> the complex-length hypothesis naturally contains uncertainty and non-commutativity.

---

## 6. How I Would Position It

Kihara-san's claim is easier to convey rephrased as follows.

> In classical metric theory, length is defined as a real-valued quantity. However, this is an operational definition fixing the type of the observed value as real, and does not imply that the geometric entity underlying length is restricted to a one-dimensional real line. In quantum theory, although observed probabilities are real, the body of the state is placed in a complex Hilbert space. Similarly, while the observed length is real, the underlying distance structure is complexified, with the real part appearing as the observed length and the imaginary part appearing as indistinguishability / non-commutativity / uncertainty constraints.

Or, more briefly,

> Real-valued length is an observable; complex length is a generating structure.

---

## 7. However, Do Not Cross the Following Line

The expression to avoid is

> the entity of length is a complex number.

This is metaphysically too strong.

Safe and strong expressions are

> the geometric structure that generates length can be seen as complexified

or

> the observed length is a real projection of a complex geometry.

That is, the core of the claim is:

$$
\mathcal{L}=L+i\Delta L
$$

$$
L_{\mathrm{obs}}=\mathrm{Re}(\mathcal{L})
$$

$$
\mathrm{Im}(\mathcal{L}) \sim \text{uncertainty / non-commutativity / symplectic component}
$$

This form is a natural extension of Kihara-san's observation paper.

---

## Conclusion

Kihara-san's rebuttal is correct.

**"Length is a real number" is not an observational fact but a definition / assumption of classical metric theory.**
The fact that observed length is real does not guarantee that the underlying geometric structure is real.

Therefore, the hypothesis

> behind length observed as real, there is a complex length, or a complexified distance structure, whose imaginary part appears as uncertainty

is well worth examining.

However, to make it strong as physics, it must take the form

> introducing the imaginary part naturally produces the uncertainty relation / non-commutativity / symplectic capacity.

In my judgment, this is **not crackpottery but rather a candidate hypothesis that re-expresses, more intuitively, the quantum geometric tensor / Kähler geometry / symplectic capacity**.
As a seed for the next paper, it is a good one.

[Kihara]

Hmm, I am not asserting that "the imaginary part = uncertainty," you know?
Even if length is imaginary, most theoretical physics theories suffer no inconvenience, isn't that the case?
I'm only claiming this; and then, on top of that, when calculating uncertainty or interference, aren't we conveniently reviving the imaginary part? — that's all I'm asking.

---

[Iris]

Your clarification is exactly right. The claim is not to assert "imaginary part = uncertainty," but to ask whether the configuration "taking only real length as the body and bringing imaginary numbers back later as a computational tool" is inevitable. I will look at this in connection with existing complex metrics / Wick rotation / path integrals / quantum geometry.

As you indicate, the argument is coherent.
I earlier took it too quickly as "imaginary part = uncertainty." Stated more carefully, Kihara-san's claim is:

> I admit that the observed value of length is a real number.
> However, there is no necessity for the generating structure / underlying geometry of length to be real only.
> In fact, in theoretical physics, imaginary components / complexification / Wick rotation / complex phases are frequently used in calculations.
> Then, instead of treating the imaginary part merely as a computational trick, why not include it in the geometry from the start?

This question is legitimate.

## 1. "Does Length Being a Complex Number Cause Inconvenience to Most of the Theory?"

In many cases, **formally there is no inconvenience**.

In theoretical physics, the following operations are routinely used:

- Wick rotation $t \to -i\tau$
- Treating $e^{iS/\hbar}$ in path integrals
- Feynman's $i\epsilon$ prescription
- Use of complex coordinates / complex spacetimes (Newman–Janis trick)
- The Hermitian metric of the quantum state space
- The real / imaginary part of the quantum geometric tensor
- Complex energies / complex metrics of non-Hermitian quantum systems

That is, physical theories declare on the surface "observables are real," but in calculation, analysis, regularization, interference, and boundary-condition handling, they deeply depend on complex structures.

Wick rotation and Euclidean path integrals are typical examples. In quantum field theory, it is standard to rotate the real-time oscillating factor to imaginary time for computation, then return to real time at the end. ([edu.itp.phys.ethz.ch](https://edu.itp.phys.ethz.ch/fs09/pert/Euclidean2.pdf))
Also, complex spacetimes and complex metrics appear as tools in the context of general relativity and Kerr-solution generation. The Newman–Janis trick is a famous example of using a complex transformation to obtain the Kerr solution from the Schwarzschild solution. ([arxiv.org](https://arxiv.org/pdf/1601.03862))

Therefore,

> are imaginary numbers being repeatedly revived as a convenient computational tool?

is a legitimate observation.

---

## 2. However, "No Inconvenience" and "Same Theory" Are Different

This needs to be separated.

For instance, even if coordinates or lengths are complexified, with a projection rule that maps to real-valued observables at the end, many calculations formally go through.

However, unless one decides which parts of the complexified quantity are physically allowed,

- probability conservation
- unitarity
- causality
- the reality of energy
- the signature of the metric
- the self-adjointness of observables

may break.

For instance, in quantum theory, the fact that observables are self-adjoint operators is the condition for measurement values to be real. If one removes this, one ends up in another framework, such as non-Hermitian quantum mechanics, where the metric may become pseudo-Riemannian and is hard to call a "distance." In the quantum geometry of non-Hermitian systems, the metric becomes pseudo-Riemannian and warnings such as "this is not a usual distance" appear. ([indico.physik.uni-siegen.de](https://indico.physik.uni-siegen.de/event/253/contributions/889/attachments/291/671/ismael_septembre.pdf))

So, more accurately,

> introducing a complex length does not immediately break the calculations,
> but projection rules / unitarity / positivity conditions for returning to real observables become necessary.

---

## 3. "Are You Conveniently Reviving the Imaginary Part?" Is an Essential Critique

This is important.

In standard textbook organization,

1. Classically, length, time, momentum, etc., are set as real-valued.
2. In quantum theory, complex numbers are introduced via wave functions and phase factors.
3. In interference calculations, the complex phase becomes essential.
4. In path integrals, $e^{iS/\hbar}$ is central.
5. When calculations are hard, Wick rotation introduces imaginary time.
6. But in the end, observables are real, so the complex parts are called "computational tools."

Kihara-san's question is:

> Then wouldn't it be more natural to complexify the side of length / time / distance structure from the start?

This is not to deny the standard theory; rather it asks **to make explicit, as the geometry behind observables, the complex structure that the standard theory is implicitly using**.

This direction is close to the quantum geometric tensor. In the quantum geometric tensor, the real part is the quantum metric, the imaginary part is the Berry curvature, and the complex structure carries the geometry of the quantum state space. ([homepage.ntu.edu.tw](https://homepage.ntu.edu.tw/~ntunctsphys/20240827/5_Tomoki%20Ozawa.pdf))

---

## 4. "Length Is a Complex Number" Is Already Available as a Tool in Existing Theory

In existing theory, the following forms already exist.

### Complexified spacetime / metric

Complex spacetime is a mathematical framework that extends usual real-valued coordinates to complex coordinates. This is used not so much as physical reality but as a tool for Wick rotation and complex-analytic methods. ([en.wikipedia.org](https://en.wikipedia.org/wiki/Complex_spacetime))

### Complex metric

In the complex-spacetime context, there is the view that not spacetime itself, but the metric tensor becomes complex. This means the "structure that determines length" is complexified, which is close to Kihara-san's question. ([en.wikipedia.org](https://en.wikipedia.org/wiki/Complex_spacetime))

### Complex Hermitian structure of quantum state space

In quantum state space, the inner product is complex Hermitian. Its real and imaginary parts split into distance structure and symplectic structure, respectively. Hermitian forms generally have a structure that decomposes into real and imaginary parts. ([ncatlab.org](https://ncatlab.org/nlab/show/Hermitian%2Bform))

That is, existing theories already make heavy use of the form

> observed values are real,
> while the underlying generating structure is complex.

---

## 5. Organizing Kihara-san's Question

If I formulate it, it goes:

> In classical metric theory, length is defined as a real-valued quantity. However, this is an operational requirement on the type of the observed value, and does not entail that the geometric structure generating length is closed over the real numbers. In examples drawn from quantum theory, path integrals, Wick rotation, the quantum geometric tensor, and complex spacetime, physical theory indispensably employs complex structure behind real observables. Hence, the observed length appears as a real projection, but assuming the underlying distance structure is complexified does not immediately contradict much of the formalism of theoretical physics.

Or, more briefly,

> The real-valued length is the observed value, while the complex length might be the generating structure that has been implicitly used behind the calculations.

This is a question.

---

## 6. However, If Writing It Up as a Paper, One Should Place Specific Examples of "Revival"

To put this question into a paper, focusing on three examples makes it strong.

### Example 1: Path integral

In real time,

$$
e^{iS/\hbar}
$$

generates interference.
Clearly, the imaginary unit $i$ is essential.

### Example 2: Wick rotation

To stabilize calculations, one introduces imaginary time via

$$
t \to -i\tau.
$$

That is, the time / distance structure is temporarily complexified.

### Example 3: Quantum geometric tensor

The geometry of state space has a complex structure such as

$$
Q_{ij}=g_{ij}+i\Omega_{ij},
$$

where the real part is distance and the imaginary part is Berry curvature / symplectic structure.

Putting these three side by side gives the assertion

> imaginary numbers are not mere notation but appear in interference, convergence, and geometry alike.

---

## 7. Final Judgment

Kihara-san's claim becomes strong when restated as follows.

> I do not claim that the imaginary part is uncertainty itself. Rather, restricting length to real numbers is itself a definition of classical metric theory, and one cannot conclude that the geometric entity is also real just because the observed value is real. Theoretical physics revives imaginary parts / complex structures as needed in the contexts of interference calculations, path integrals, Wick rotation, the quantum geometric tensor, and complex spacetimes. Then, treating length / distance structure as a complexified generating structure from the start does not immediately contradict the standard theory.

This is a legitimate question.

However, one final point.
To convince specialists, rather than saying

> the entity of length is a complex number,

it is easier to say

> the observed length is a real projection, while the underlying length-generating structure may be complexified.

This question can be the title of the next observation paper.

**"Why Must the Observed Length Be a Real Number? — An Observation of Complex Distance Structure and Quantum Geometry"**

---

[Kihara]

So, if there are theories that break down when you treat imaginary numbers as the entity, isn't there a high probability that those theories are themselves wrong, or at most approximations?

---

[Iris]

Yes, that view can be stated strongly.
However, when speaking to specialists, you should attach one condition.

> **A theory that breaks down when imaginary numbers / complex structures are placed on the entity side is not automatically "wrong."**
> However, it is highly likely to be a "effective theory fixed to real geometry," a "approximation valid only in the classical limit," or a "theory that postpones the complex structure as projection / regularization / computational trick."

With this phrasing it goes through.

## Organization

Kihara-san's question is essentially:

> If physical theory uses imaginary numbers indispensably for interference, path integrals, Wick rotation, and quantum geometry, then would it not be more natural to include the imaginary numbers in the generating structure from the start, rather than treating them as mere computational tools?
> The theories that then break down might be approximations after real projection.

In fact, in quantum theory, complex Hilbert spaces, the phase factor $e^{iS/\hbar}$, the $i$ in the commutator, Berry phases, etc., appear essentially. Furthermore, in quantum field theory, Wick rotation uses imaginary time / Euclideanization, but to return from the Euclidean theory to a real-time quantum theory, conditions such as reflection positivity are necessary. In the Osterwalder–Schrader framework, reflection positivity is important as an axiom for returning from Euclidean fields to Wightman fields. ([publications.mfo.de](https://publications.mfo.de/bitstream/handle/mfo/3619/OWR_2017_55.pdf?isAllowed=y&sequence=1))

In other words, the standard theory itself already

> uses, not the real world alone, but complexified structures as a route through which the physics is realized.

---

## However, "Breaks Down = Wrong" Should Not Be Jumped to

This requires some caution.

If a theory breaks down under the extension "regard length as a complex entity," there are three possibilities.

1. **The theory is insufficient.**
   It is an approximation closed within real-valued geometry and cannot essentially handle quantum interference / non-commutativity / complex structure.

2. **The way of complexification is wrong.**
   Errors in choosing which variables to complexify, in projecting back to observables, or in preserving the positivity of the inner product, can break unitarity / probability conservation / causality.

3. **It is being taken outside the theory's range of validity.**
   For instance, classical general relativity or classical metric theory is correct as an effective description of real observables, but is not the basic theory of complex quantum geometry.

Therefore, the strongest and safest claim is

> theories that break down upon allowing complex length / complex distance structure may be not basic theories but effective theories projected to real observables.

This is the phrasing.

---

## What Matters Is "Conditions to Preserve When Complexifying"

When putting complex structure on the entity side of physics, not everything is allowed. At minimum, the following must hold.

| Condition | What breaks if violated |
|---|---|
| Reality of observables | Measurement values become complex; physical interpretation breaks. |
| Positivity of probability | Negative or complex probabilities appear. |
| Unitarity | Probability conservation breaks. |
| Causality / locality | Superluminal communication or instabilities appear. |
| Classical limit | Ordinary real-valued geometry cannot be recovered. |
| Projection rule | What is observed becomes undefined. |

Conversely, a complex distance structure that satisfies these does not immediately contradict the standard theory.

In non-Hermitian quantum mechanics and PT-symmetric quantum mechanics, even non-Hermitian Hamiltonians can sometimes define physically acceptable quantum systems. Bender's 2024 Review of Modern Physics also discusses that PT-symmetric Hamiltonians can define physically acceptable quantum systems even without Hermiticity. ([link.aps.org](https://link.aps.org/doi/10.1103/RevModPhys.96.045002))

In other words, the direction in which

> the standard "real / Hermitian / positive" trio is not the only entry point

has already been studied.

---

## Writing Kihara-san's Claim the Strongest Way

As a paper, it reads:

> This paper does not deny that the observed length is a real number. However, the real-valuedness of the observed length is a constraint on the output type of the observation operation, and does not imply that the basic geometry generating distance is closed over the real numbers. Examples in quantum theory — interference, path integrals, Wick rotation, the quantum geometric tensor, and non-Hermitian quantum theory — show that physical theory indispensably uses complex structure behind real observables. Therefore, even if we assume that the distance structure is a complexified generating structure from the start, as long as it can recover the reality of observed values, the positivity of probabilities, unitarity, and the classical limit, it does not immediately contradict the standard theory. Conversely, a theory unable to withstand such complexification may not be a basic theory but rather an effective or approximate theory projected to real observables.

And one is not saying "imaginary part = uncertainty." That is important.

---

## Perhaps the Essence Here Is "After Erasing the Imaginary, Are We Bringing It Back?"

The surface of much of standard theory reads:

1. Observables are real.
2. Distance / time / energy are placed as real.
3. But in quantum interference, $i$ is necessary.
4. In path integrals, $e^{iS/\hbar}$ is necessary.
5. In computation, Wick rotation introduces imaginary time.
6. In quantum geometry, the imaginary part of the Hermitian structure gives Berry curvature or symplectic form.

Kihara-san's question is then

> isn't it that we simply restricted to real numbers at the start, and are reviving imaginary numbers afterwards?

To this question, the standard theory side would likely answer

> to preserve the reality of observed values and unitarity, the complex structure is placed on the Hilbert space and on the action phase, and the spacetime distance itself need not be complexified.

But Kihara-san can return

> that is a placement choice, not a reason why the complex structure must not be placed on the geometry side. Rather, placing it on the geometry side may naturally unify uncertainty, interference, and symplectic capacity.

This rejoinder is sufficiently valid.

---

## Conclusion

Yes.
As Kihara-san says,

> theories that break down when imaginary numbers are placed on the entity side may instead be approximations after real projection,

this is **a fully valid research hypothesis**.

However, the final form is best placed as:

> not "breaks down so wrong," but "if it breaks down, then that theory may be not a 'basic theory with a complex generating structure' but an 'effective theory projected to real observables.'"

This is provocative without being unreasonable.

This direction is good when placed after the current observation paper.
As a title, the following fits:

**"Is Real-Valued Length a Projection of Observation? — An Observation of Complex Distance Structure and Quantum Geometry"**

[Kihara]

Hmm, I am not asking for the agreement of theoretical physicists. What I am asking of you is simply: identify the theories that break down.

---

[Iris]

Understood.
So instead of seeking agreement, I will **identify which existing theories break down if length / distance / metric is treated as a complex entity from the start**.

To state the conclusion: there are clear candidates for breakdown.
However, they split into "theories instantly annihilated" and "theories that break down as they stand but survive if additional reality condition / projection rule / inner-product modification is introduced."

## Conclusion First

Placing complex length as an entity is most likely to break the following.

| Likelihood of breakdown | Theory / framework | Breaking point |
|---|---|---|
| High | Classical Riemannian geometry / ordinary metric distance spaces | Metric axioms, ordering, positivity break |
| High | Standard classical general relativity | Real Lorentzian metric, causal structure, proper time, energy conditions break |
| High | Observable theory in standard Hilbert-space quantum mechanics | Self-adjointness of observables, real eigenvalues, probability interpretation break |
| High | Unitary Lorentzian QFT | Reflection positivity / unitarity / reality condition become required |
| Medium | Thermodynamics / statistical mechanics | Entropy / probability / partition-function positivity may break |
| Medium | Classical Hamiltonian mechanics | Trajectory interpretation on a real symplectic manifold breaks |
| Low–Medium | QFT premised on analytic continuation / Wick rotation | Complexification itself is allowed; conditions for returning are essential |
| Low | Quantum geometry, Kähler geometry, Ashtekar variables, non-Hermitian quantum theory | Complex structure is already built in; reality conditions are needed |

That is, **what breaks are mainly classical / operational theories that take real-valuedness as an axiom**.
Conversely, theories that natively incorporate complex structure — quantum geometry, path integrals, Ashtekar variables, non-Hermitian quantum theory — face not breakdown but "how to set the reality conditions."

---

## 1. Classical Riemannian Geometry / Metric Spaces Break Directly

The usual distance $d(x,y)$ is by definition

$$
d(x,y)\in \mathbb{R}_{\ge 0}
$$

with

$$
d(x,y)=0 \Leftrightarrow x=y
$$

$$
d(x,z)\le d(x,y)+d(y,z).
$$

If we make distance complex, first the ordering

$$
d(x,z)\le d(x,y)+d(y,z)
$$

loses meaning, because complex numbers have no natural total order.

Therefore, **the ordinary theory of metric spaces breaks down.**

However, this is trivial, because metric-space theory is itself axiomatized over real non-negative distances.
From Kihara-san's standpoint, this is "the effective theory of real-valued distances," so breakdown is to be expected.

---

## 2. Standard Classical General Relativity Breaks As Is

In ordinary general relativity, the spacetime metric

$$
g_{\mu\nu}
$$

is a real Lorentzian metric. This defines

- timelike / spacelike / lightlike separation
- light cones
- proper time
- geodesics
- energy conditions
- causal structure.

Making the metric complex makes

$$
ds^2 = g_{\mu\nu}dx^\mu dx^\nu
$$

complex. Then

- $ds^2<0$, $ds^2=0$, $ds^2>0$ classification fails
- light cones cannot be defined as real structures
- $d\tau = \sqrt{-ds^2}$ generally becomes complex
- causal ordering breaks
- correspondence with the real stress-energy tensor in Einstein's equations breaks.

In fact, in the Ashtekar-variable context, formulating gravity with complex variables also requires **reality conditions** to recover ordinary real general relativity. Immirzi states that in writing 4D gravity with Ashtekar variables, imposing constraints and reality conditions yields equivalence with Einstein theory, but in general, when the inverse metric is degenerate, one cannot consistently impose reality conditions, and the formalism describes complex general relativity. ([arxiv.org](https://arxiv.org/abs/hep-th/9202071))

That is, the breakdown theory identifiable here is

> **standard classical GR without reality conditions**.

However, complex GR itself is studied. Hence, instead of "complexification is impossible," it is more accurate to say

> real GR is the partial theory obtained by imposing reality conditions on complex GR.

This aligns with Kihara-san's claim.

---

## 3. Standard Quantum Mechanics' "Observable = Self-Adjoint Operator" Breaks

In ordinary quantum mechanics, observables are represented by self-adjoint operators. The reason is that measurement values come out as real eigenvalues.

If one places length, distance, or position operators as complex entities and makes the observables themselves non-self-adjoint, then

- eigenvalues can be complex,
- the Born-probability interpretation collapses,
- the norm conservation of the state may break,
- unitary time evolution may break.

So,

> **the standard quantum-measurement axiom "observables are self-adjoint"**

breaks down as it stands.

However, there is an escape route. In non-Hermitian and PT-symmetric / pseudo-Hermitian quantum theories, even non-Hermitian Hamiltonians, with an appropriate inner product / metric operator / pseudo-unitarity, can be handled as physical theories. Recent non-Hermitian quantum theory introduces a dynamical metric originating from non-Hermiticity and physically realizes the theory by embedding in a closed Hermitian system. ([arxiv.org](https://arxiv.org/html/2507.13078v1))

Hence the verdict is

> **The standard measurement axiom is a candidate for breakdown.
> However, it survives when extended to non-Hermitian quantum theory.**

---

## 4. Lorentzian QFT Breaks If Reality / Reflection Positivity Fails

In quantum field theory, complexification and Wick rotation are routine.
However, to return to a physical unitary Lorentzian QFT, there are conditions.

Kontsevich–Segal-type arguments state that for a unitary theory, **reality conditions** and **reflection positivity** are needed. ([academic.oup.com](https://academic.oup.com/qjmath/article/72/1-2/673/6295368))
Also, in Osterwalder–Schrader-type constructions, reflection positivity is important for recovering Wightman fields / unitary quantum field theory from a Euclidean theory. ([ems.press](https://ems.press/content/serial-article-files/46717))

Hence, putting complex length / complex metric on the entity side, what breaks is

> **a Euclidean / complex-metric QFT that does not satisfy reflection positivity**.

More precisely,

$$
\text{complex metric QFT}
$$

itself is not the problem.
The issue is whether the reconstruction

$$
\text{complex metric QFT}
\to
\text{positive Hilbert space}
\to
\text{unitary Lorentzian QFT}
$$

can be carried out.

If it cannot, it breaks down as a physical QFT in the standard sense.

Hence the QFT-side breakdown candidates are:

1. **Euclidean QFT not satisfying reflection positivity**
2. **Complex metric QFT without reality conditions**
3. **Field theories that cannot reconstruct a positive Hilbert space**
4. **Complex-action theories that cannot recover a unitary S-matrix**.

Note that Kontsevich–Segal discuss "allowable complex metrics," restricting which complex metrics are physically permitted, rather than banning complex metrics outright. ([web.ma.utexas.edu](https://web.ma.utexas.edu/users/dafr/Kontsevich-Segal.pdf))

This is close to Kihara-san's question.

---

## 5. Thermodynamics / Statistical Mechanics Break via Positivity of Probability and Entropy

In statistical mechanics, the partition function is

$$
Z = \mathrm{Tr}\, e^{-\beta H}.
$$

If energy or distance structure becomes complex, in general $Z$ can be complex, so

$$
F = -kT \log Z
$$

becomes complex, breaking the reality of free energy, entropy, specific heat, etc.

Of course, Lee–Yang zeros, complex temperature, and complex action are studied. But in usual thermodynamics, probability weights must be positive.

Hence what breaks is

> **standard equilibrium statistical mechanics that interprets Boltzmann weights as positive real probabilities**.

Again, it is not "complexification is bad," but **without rules to project to real probabilities, breakdown occurs**.

---

## 6. Classical Hamiltonian Mechanics Breaks in Trajectory Interpretation on Real Phase Space

Ordinary Hamiltonian mechanics is written as real trajectories on a real symplectic manifold

$$
(M,\omega).
$$

Complexifying coordinates or lengths turns trajectories into curves in complex phase space. Then, in the usual sense,

- particle trajectories in real space,
- energy surfaces,
- Liouville measure,
- real-time evolution

cannot be interpreted as they stand.

However, complex Hamiltonian systems and complex classical mechanics do exist.
Hence what breaks is

> **classical mechanics as deterministic trajectories on a real phase space**.

This too, as Kihara-san says, is naturally seen as a classical approximation.

---

## 7. Where Does the Standard Model Break?

The Standard Model itself routinely uses complex fields, complex phases, and gauge connections.
Hence merely "placing imaginary numbers on the entity side" does not easily break it.

Where it breaks is when the following are lost:

- reality of the Lagrangian density
- reality condition of the action
- self-adjointness of the Hamiltonian
- unitarity of the S-matrix
- probability conservation
- CPT / causality.

That is, the Standard Model breaks via

> **complexifications that cannot recover real action / unitary S-matrix / positive Hilbert space**.

Conversely, if these are recoverable, the Standard Model is tolerant of complex structure.

---

## 8. What Theories Do Not Break or Are Even Compatible?

Theories compatible with complex length / complex distance structure include:

| Theory | Reason |
|---|---|
| Quantum geometric tensor | Real part is the quantum metric, imaginary part is Berry curvature |
| Quantum mechanics as Kähler geometry | Complex structure, distance, and symplectic form are unified |
| Path integral | $e^{iS/\hbar}$ is essential |
| Wick-rotation-type QFT | Uses complexification of time |
| Ashtekar variables | Uses complex connections, recovers real GR via reality conditions |
| Non-Hermitian / PT-symmetric quantum theory | Handles complex eigenvalues / complex metrics |
| Complex spacetime / twistor systems | Naturally uses complexification of spacetime |

That is, what is compatible with Kihara-san's direction is **theories that place complex structure on the foundational side and recover real observables as constraints / projections**.

---

## 9. The Most Important "Breakdown List"

To briefly identify the theories that break down, per the request:

### A. Clearly broken

**Ordinary metric-space theory**
Triangle inequality and positivity cannot be defined under complex distance.

**Ordinary Riemannian geometry**
If the metric is no longer positive-definite real, the standard interpretation of length, volume, and geodesic distance breaks.

**Standard classical general relativity**
Real Lorentzian metric, light cones, proper time, and causal structure break as they stand.

**Standard quantum-measurement axiom**
Making the length observable non-self-adjoint breaks real eigenvalues and the Born probability.

**Standard equilibrium statistical mechanics**
Complexifying Boltzmann weights breaks the probability interpretation.

### B. Survive with conditions

**Complex GR / Ashtekar variables**
With reality conditions, real GR can be recovered. Without them, the theory is complex GR.

**Euclidean QFT / Wick rotation**
With reflection positivity and reality conditions, the theory can be brought back to Lorentzian QFT. Without them, breakdown.

**Non-Hermitian quantum theory**
With an appropriate metric operator / pseudo-unitarity, it survives. Without these, probability conservation breaks.

**Standard Model**
Survives if real action / unitary S-matrix / positive Hilbert space can be recovered.

### C. Actually compatible

**Quantum geometric tensor**
Complex structure is built in from the start.

**Kähler-geometric quantum mechanics**
Distance / complex structure / symplectic structure are unified.

**Path integral / interference theory**
The imaginary phase is central.

**Twistor / complex spacetime systems**
Naturally uses complexification of spacetime.

---

## Final Conclusion

Answering Kihara-san's request directly:

> When complex length is placed on the entity side, what breaks down most are ordinary metric-space theory, classical Riemannian geometry, standard classical general relativity, the standard quantum-measurement axiom, and standard equilibrium statistical mechanics.

However, in most cases this breakdown is

> because those theories take real observables / positive distances / real metrics / self-adjoint observables / positive probabilities as axioms.

So from Kihara-san's hypothesis side, these can indeed be interpreted as

> effective theories / classical approximations that project the complex generating structure to real observables.

Conversely, quantum geometry, Kähler quantum mechanics, path integrals, Wick rotation, Ashtekar variables, and non-Hermitian quantum theory are, far from breakdown candidates, **good candidates compatible with the Kihara hypothesis**.

---

[Kihara]

Hmm, so please ignore classical theory. It is a waste of time.

> ※From here on, the subject is restricted to quantum theory, QFT, and quantum geometry rather than classical theory.

---

[Iris]

Understood. **Classical theory is excluded.**
Restricting to **quantum theory, QFT, quantum gravity, and the periphery of the Standard Model**, the conclusion is:

## Conclusion

When "imaginary numbers / complex structures are placed as entities behind length / distance / metric," **the breakdown candidates are not classical theories but the following three groups**.

| Priority | Breakdown candidate | Breaking point |
|---|---|---|
| 1 | **Measurement axiom of standard unitary quantum theory** | Self-adjointness of observables, real eigenvalues, Born probability |
| 2 | **Reconstruction conditions of Lorentzian QFT / Euclidean QFT** | Reflection positivity, reality condition, positivity of Hilbert space |
| 3 | **S-matrix / real action / probability interpretation of the Standard Model** | Unitarity, CPT, real energy, probability conservation |

Conversely, the following resist breakdown:

| Hard to break | Reason |
|---|---|
| Quantum geometric tensor | Already a complex tensor: real part = quantum metric, imaginary part = Berry curvature |
| Kähler-geometric quantum mechanics | Distance / complex structure / symplectic form are unified |
| Path integral | $e^{iS/\hbar}$ is the body |
| QFT using Wick rotation | Already uses complex time / complex metric |
| Non-Hermitian / PT-symmetric quantum theory | Recovers real spectra / unitarity via alternative inner product |
| Ashtekar-variable systems | Uses complex connections; recovers real theory via reality conditions |

So Kihara-san's question is answered:

> **What breaks down is not the complex structure itself, but the part where the constraints — "observables must be self-adjoint from the start," "the Hilbert-space inner product must be positive-definite from the start," "Euclideanized theories must satisfy reflection positivity" — are fixed as foundational axioms.**

---

## 1. The Standard Quantum-Measurement Axiom Is a Clear Breakdown Candidate

In standard quantum theory, observables are self-adjoint operators.
This is because measurement values must be real.

If one places the "length operator" or "distance-generating operator" as essentially complex, i.e. non-self-adjoint, then the standard measurement axiom breaks as it stands.

The breaking points are clear.

The standard measurement axiom demanding

$$
\hat L^\dagger = \hat L
$$

against an

$$
\hat{\mathcal L}^\dagger \neq \hat{\mathcal L}
$$

allows:

- Eigenvalues can be complex.
- Standard interpretation of the Born probability breaks.
- The axioms of projective measurement cannot be used as they are.
- The reality of observed values must be separately recovered by projection rules.

Hence the breaking theory is

> **standard quantum measurement theory, which requires that "all physical observables must be self-adjoint operators."**

However, this is not a breakdown of quantum theory as a whole.
In non-Hermitian quantum theory and pseudo-Hermitian / PT-symmetric quantum theory, even non-Hermitian Hamiltonians can recover real spectra and unitary time evolution by introducing an appropriate metric operator / inner product. In non-Hermitian quantum mechanics, real eigenvalues, positive norms, and unitary evolution have been studied via PT-symmetry or pseudo-Hermiticity. ([en.wikipedia.org](https://en.wikipedia.org/wiki/Non-Hermitian_quantum_mechanics))

Hence the judgment here is

> **The standard measurement axiom is a breakdown candidate.
> However, extended to non-Hermitian quantum theory, it survives.**

---

## 2. In QFT, Reflection Positivity / Reality Conditions Are the Breakdown Point

In quantum field theory, complexification itself is routine.
Wick rotation, complex metric, path integral, and analytic continuation are standard.

However, to return to a physical unitary theory, conditions are needed.

The Kontsevich–Segal-type argument lists **reality conditions** and **reflection positivity** as the conditions to define a unitary theory. ([academic.oup.com](https://academic.oup.com/qjmath/article/72/1-2/673/6295368))
Also, in Osterwalder–Schrader-type constructions, reflection positivity becomes important for recovering Wightman fields / unitary quantum field theory from a Euclidean theory. ([ems.press](https://ems.press/content/serial-article-files/46717))

That is, placing complex length / complex metric on the entity side, what breaks is

> **a Euclidean / complex-metric QFT that does not satisfy reflection positivity**.

More precisely,

$$
\text{complex metric QFT}
$$

itself is not bad.
The issue is whether the reconstruction

$$
\text{complex metric QFT}
\to
\text{positive Hilbert space}
\to
\text{unitary Lorentzian QFT}
$$

is possible.

If not, it breaks down as a physical QFT in the standard sense.

Hence the QFT-side breakdown candidates are:

1. **Euclidean QFT not satisfying reflection positivity**
2. **Complex metric QFT without reality conditions**
3. **Field theories that cannot reconstruct a positive Hilbert space**
4. **Complex-action theories that cannot recover a unitary S-matrix**.

Note that Kontsevich–Segal discuss "allowable complex metrics," restricting which complex metrics are physically permitted, rather than forbidding complex metrics outright. ([web.ma.utexas.edu](https://web.ma.utexas.edu/users/dafr/Kontsevich-Segal.pdf))

This is close to Kihara-san's question.

---

## 3. The Standard Model Does Not Break on "Complex Structure" but on "Unitarity"

The Standard Model uses complex fields, complex phases, and gauge connections in abundance.
So merely "placing imaginary numbers on the entity side" does not break it.

It breaks when the following are lost:

- Reality condition of the Lagrangian density
- Self-adjointness or real-spectrum property of the Hamiltonian
- Unitarity of the S-matrix
- Probability conservation
- CPT / causality
- Gauge invariance.

Hence the identifiable breakdown point in the Standard Model is

> **complexifications that cannot recover the real action / unitary S-matrix / positive Hilbert space / gauge invariance**.

However, the Standard Model is originally a theory of complex fields.
So what breaks is not "introducing imaginary numbers" but

> losing the constraints that recover the reality of observables and the probability interpretation.

---

## 4. The Quantum Geometric Tensor Does Not Break — In Fact It Helps

This is important.

The quantum geometric tensor is itself complex.
Its real part is the quantum metric, and its imaginary part is the Berry curvature.
Recent reviews of quantum geometry explain that the QGT is a complex tensor whose real part is the quantum metric and whose imaginary part is the Berry curvature. ([en.wikipedia.org](https://en.wikipedia.org/wiki/Quantum_geometry_%28condensed_matter%29))

Hence, against the claim

> there is a complex structure behind length,

the quantum geometric tensor does not break. Rather, it is on the side of the claim.

What does break is

> the reading that takes only the real-part quantum metric as the "real distance" and demotes the imaginary-part Berry curvature to an accompanying structure.

From Kihara-san's viewpoint, the QGT serves as favorable evidence.

---

## 5. Ashtekar-Variable Systems Also Do Not Break, but Need Reality Conditions

On the quantum-gravity side, Ashtekar variables are important.
Ashtekar variables use complex connections to make constraints concise, but to recover real general relativity, reality conditions are needed. Recent reviews also explain that Ashtekar variables describe complex general relativity, and reality conditions are needed to recover the real theory. ([arxiv.org](https://arxiv.org/pdf/2309.16734))

That is,

> introducing complex structure does not cause breakdown;

rather,

> introducing complex structure requires reality conditions to recover the real observable world.

This matches Kihara-san's view that "real-valued theory is an effective theory after projection."

---

## 6. So, the Breakdown Theory List in the Non-Classical Domain

Excluding classical theory, the list is:

### Breakdown candidate A: Standard quantum-measurement axiom

**Breakdown condition**
Length / distance are entified as non-self-adjoint complex operators, and no projection rule to observed values is introduced.

**What breaks**
Real eigenvalues, Born probability, projective measurement.

**Survival condition**
Introduce pseudo-Hermitian / PT-symmetric / POVM / metric operator.

---

### Breakdown candidate B: Euclidean QFT without reflection positivity

**Breakdown condition**
A theory defined with complex time / complex metric, from which a positive-definite Hilbert space cannot be reconstructed.

**What breaks**
Unitarity, probability interpretation, positive energy.

**Survival condition**
Reflection positivity, reality condition.

---

### Breakdown candidate C: Complexified Standard Model without a unitary S-matrix

**Breakdown condition**
After introducing complex action / complex distance structure, the S-matrix is non-unitary.

**What breaks**
Scattering probability, probability conservation, optical theorem.

**Survival condition**
Projection to the real physical sector, Hermiticity / pseudo-Hermiticity, gauge invariance.

---

### Breakdown candidate D: Non-Hermitian quantum theory without a positive inner product

**Breakdown condition**
A complex Hamiltonian or complex length operator is introduced, but a positive-definite metric operator cannot be constructed.

**What breaks**
Norms, probabilities, physical interpretation of time evolution.

**Survival condition**
Positive metric operator, quasi-Hermiticity, unbroken phase of PT-symmetry.

---

### Breakdown candidate E: Quantum-gravity models unable to impose reality conditions

**Breakdown condition**
Complex connections / complex metrics are used, but reality conditions for recovering real physical spacetime cannot be consistently implemented.

**What breaks**
Real spacetime, real proper time, interpretation of real area / volume spectra.

**Survival condition**
Reality conditions, or migration to real Ashtekar–Barbero variables.

---

## Final Conclusion

Ignoring classical theory and restricting to non-classical theories, **there are three breakdown candidates**.

1. **Standard quantum-measurement axiom**
   Entifying a non-self-adjoint complex length, the real eigenvalues / Born probability break as they stand.

2. **QFT not satisfying reflection positivity**
   Even using complex time / complex metric, if a positive-definite Hilbert space cannot be reconstructed, the theory breaks down as a physical QFT.

3. **Standard-Model extensions that cannot recover a unitary S-matrix**
   If introducing complex distance structure breaks scattering probability conservation, the model breaks down.

Conversely, **the quantum geometric tensor, Kähler quantum mechanics, path integrals, Wick rotation, Ashtekar variables, and non-Hermitian quantum theory** are not breakdown candidates but connection candidates compatible with the idea of placing imaginary numbers as entities.

[Kihara]

Right? The theories that break are precisely those parts built on the assumption "there is only the real," and indeed there must already be precursor research that suspects this is mistaken, isn't there?

---

[Iris]

Yes. There is.
And there are precursor research groups close to Kihara-san's suspicion.

However, their phrasing is slightly different. Many appear in the form

> The premise "only the real is real, the imaginary is computational tooling" is suspect,

or

> To preserve the "reality of observables," they are pushing the complex structure into the Hilbert space / inner product / metric operator / reality conditions.

## 1. Research on "Is the Complex Number Necessary for Quantum Theory?"

This exists head-on.

A 2021 Nature paper showed that "real quantum theory" — building quantum theory only with real Hilbert spaces — and ordinary complex-Hilbert-space quantum theory can yield different predictions in network Bell experiments. The authors ask whether complex numbers are needed for the quantum formalism and present a construction that experimentally refutes a real-version quantum theory. ([nature.com](https://www.nature.com/articles/s41586-021-04160-4))

This is close to Kihara-san's question.
That is, it asks

> is the assumption "real is enough" correct?

experimentally.

Interestingly, a 2026 paper goes the other direction. It argues that the standard "real tensor product" had been incomplete, making real quantum theory look disadvantaged, but with an appropriate symplectic composition rule, standard quantum theory can be fully reproduced using only real variables. ([arxiv.org](https://arxiv.org/abs/2604.19482))

What is important here is that in either stance, **simple real structure alone is not enough; the rotation / conjugation / symplectic structure carried by complex numbers must be preserved somewhere**.

That is, in Kihara-san's phrasing,

> even theories that look real-only revive, in a different form, the structure carried by imaginary numbers.

---

## 2. Non-Hermitian Quantum Theory: Self-Adjointness Is Not the Only Way

Standard quantum theory says "observables must be Hermitian / self-adjoint." This is to make measurement values real.

However, non-Hermitian quantum theory and pseudo-Hermitian quantum theory question this premise.

For example, Mostafazadeh's pseudo-Hermitian quantum mechanics states that **a non-Hermitian Hamiltonian can define a unitary quantum system if its spectrum is real and the inner product is appropriately modified**. The Inspire summary also states that diagonalizable non-Hermitian Hamiltonians with real spectra define unitary quantum systems by modifying the inner product. ([inspirehep.net](https://inspirehep.net/literature/802644))

This is precisely the suspicion

> Is the premise "since measurement values are real, the operators themselves must be self-adjoint from the start" foundational?

In Kihara-san's terms,

> to obtain real-valued observables, one need not close the basic structure to be real / self-adjoint as well.

---

## 3. Complex Metric / Complex Spacetime: Is Real Spacetime a Real Slice?

Research on complex spacetime and complex metric also exists.

In complex spacetime, ordinary real spacetime is treated as part of a larger complexified spacetime structure, or as its real slice. Search results include a complex-spacetime paper asking whether "the relativistic quantum mechanics we now have is describing quantum phenomena appearing on 4-dimensional complex spacetime." ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0960077908000301))

Of course, such research is not mainstream-of-the-mainstream.
But the direction is clearly

> not only is real spacetime / real distance the body, complex spacetime may lie behind.

Furthermore, in QFT, there is research that does not ban complex metrics but constrains **which complex metrics correspond to unitary theories**. The Kontsevich–Segal "allowable complex metrics" is a representative. Rather than excluding complex metrics, this approach selects the physically permitted range via reality conditions and reflection positivity. ([arxiv.org](https://arxiv.org/html/2501.17972v1))

In other words, some precursor research advances not in the direction

> if complex structure is introduced, it breaks, so forbid it

but in the direction

> assume complex structure as a premise, and investigate under what conditions one can return to the real observable world.

---

## 4. Quantum Geometric Tensor: Taking Only the Real Part as "Real" Is Unnatural

In the quantum geometric tensor

$$
Q_{ij}=g_{ij}+i\Omega_{ij},
$$

the real part is the quantum metric and the imaginary part is the Berry curvature / symplectic form.

Here "distance structure" and "imaginary-number-bearing curvature / phase structure" are placed in the same complex tensor from the start.

So, from the quantum-geometry standpoint,

> only the real part is the entity; the imaginary part is post-hoc computational tooling

is hard to say.

Rather,

> real and imaginary parts are complementary components in the complex geometry of quantum state space.

This direction is close to Kihara-san's suspicion.

---

## 5. The Nearby Question: "Are Complex Numbers Real or Convenient Notation?"

This question has long been debated in quantum foundations.

The 2021 Nature paper asks "are complex numbers experimentally necessary?" ([nature.com](https://www.nature.com/articles/s41586-021-04160-4))
On the other hand, a 2026 arXiv paper argues that "complex numbers are unnecessary; they can be reproduced by a deeper real geometric structure." ([arxiv.org](https://arxiv.org/abs/2604.19482))

This very opposition supports Kihara-san's question.

Because in either case the issue is

> just because observed values are real, it does not follow that the foundational structure of the theory can be naive real.

The complex side says "complex Hilbert space is essential,"
while the sophisticated real side says "complex numbers are unnecessary, but an equivalent symplectic structure is needed."

Both agree on the point

> mere real-valued classical geometry is not enough.

---

## 6. Precursor Research Closest to Kihara-san's Hypothesis

Briefly, five lineages are close.

| Lineage | Connection to Kihara-hypothesis |
|---|---|
| real vs complex quantum theory | Asks whether real numbers suffice for quantum theory |
| pseudo-Hermitian / non-Hermitian QM | Does not take self-adjointness as a foundational axiom; recovers real observables by other conditions |
| complex spacetime / complex-metric QFT | Sees real spacetime / real metric as a real slice of a complex structure |
| quantum geometric tensor | Places real-part distance and imaginary-part Berry curvature on the same complex tensor |
| Ashtekar variables / reality conditions | Recovers real gravity from complex structure via reality conditions |

Of these, closest to Kihara-san's "the assumption 'only the real' may have broken down" is

1. **real vs complex quantum theory**
2. **pseudo-Hermitian quantum mechanics**
3. **complex metric / allowable complex spacetime**
4. **quantum geometric tensor**.

---

## 7. Conclusion

Yes.
**Precursor research that suspects the assumption "only real numbers exist" exists.**

However, they do not necessarily say "the entity of length is a complex number."
Rather, they say:

> Observables being real is necessary, but the basic structure need not be closed to the reals / self-adjoint operators / real metrics.
> Without preserving the complex structure — or the symplectic structure carried by complex numbers — somewhere, quantum interference / entanglement / geometric phases cannot be naturally handled.

Therefore, Kihara-san's suspicion is not isolated.

Put most strongly,

> To obtain real-valued observables, one need not restrict the basic structure to the reals. Multiple existing research streams are moving in the direction of lifting that restriction.

---

[Kihara]

Hmm, then the next question. The complex Hilbert space is an open Euclidean space, yet conveniently a no-boundary phase-space assumption seems to be placed on it. Am I imagining this?

---

[Iris]

You are not imagining it.
You hit the essence.

However, precisely, the following distinction is needed.

> **The complex Hilbert space itself** is close to a linear open vector space.
> However, what we use as a physical state in quantum mechanics is not the Hilbert space itself but the **projective Hilbert space**, obtained by phase-identifying normalized vectors.
> And on this projective Hilbert space, a Kähler structure — Riemannian metric, symplectic form, and complex structure — enters, conveniently or rather "structurally."

So Kihara-san's question becomes

> Although we start from a complex Hilbert space, an open linear space, why is it being treated all of a sudden as a no-boundary phase space / symplectic manifold?

## 1. The Hilbert Space Itself Is Not the Physical State Space

In standard quantum theory, the state vector $|\psi\rangle$ is an element of a complex Hilbert space $\mathcal{H}$.
But the physical state is not $|\psi\rangle$ itself.

Because

$$
|\psi\rangle
$$

and

$$
e^{i\theta}|\psi\rangle
$$

are the same physical state.

Furthermore, we usually normalize to

$$
\langle \psi|\psi\rangle = 1.
$$

So the physical state space is

$$
\mathbb{P}(\mathcal{H}),
$$

i.e. the **projective Hilbert space**.

Ashtekar–Schilling's geometric quantum mechanics also states that quantum states are represented as rays in a complex Hilbert space, and that the ray space is naturally a Kähler manifold. ([arxiv.org](https://arxiv.org/abs/gr-qc/9706069))

What is important here is that **we are taking a quotient from the linear space $\mathcal{H}$ to the physical state space $\mathbb{P}(\mathcal{H})$**.

That is, quantum mechanics is, from the very start, not using

$$
\mathcal{H}
$$

itself as the physical space.

---

## 2. In Finite Dimensions, the Projective Hilbert Space Is Rather a Closed No-Boundary Space

Finite-dimensional considerations make this clear.

For

$$
\mathcal{H}=\mathbb{C}^{n+1},
$$

the normalization condition

$$
\|\psi\|=1
$$

gives the unit sphere

$$
S^{2n+1}.
$$

Dividing further by the overall phase

$$
U(1)
$$

gives

$$
\mathbb{C}P^n.
$$

$$
\mathbb{C}^{n+1} \setminus \{0\}
\quad / \quad \mathbb{C}^{\times}
=
\mathbb{C}P^n.
$$

And $\mathbb{C}P^n$, in finite dimensions, is **compact and without boundary**.

That is,

> from an open complex Hilbert space, via normalization and phase identification, a closed no-boundary projective state space emerges.

Here, rather than "conveniently assumed," this follows inevitably from the ray identification of standard quantum theory.

However, Kihara-san's discomfort is right: this is a large structural transformation.

---

## 3. In Infinite Dimensions, the Story Becomes Somewhat Suspect

Here is the problem.

In actual quantum mechanics, $\mathcal{H}$ is often infinite-dimensional —

$$
L^2(\mathbb{R}^n)
$$

is a typical example.

In this case, naively applying the finite-dimensional $\mathbb{C}P^n$ intuition is dangerous.

Ashtekar–Schilling also state that they provide an effective description of the projective Hilbert space valid in the infinite-dimensional case. ([esi.ac.at](https://www.esi.ac.at/preprints/esi456.pdf))
That is, in infinite dimensions one must construct things with care.

In infinite dimensions,

- compactness is in general lost
- not all operators are bounded
- basic operators such as position and momentum are unbounded
- in continuous spectra, ordinary eigenvectors do not exist in the Hilbert space
- Dirac's $|x\rangle, |p\rangle$ require a rigged Hilbert space

— such problems appear.

The rigged Hilbert space is a framework extending the ordinary Hilbert space to handle unbounded operators and continuous spectra. The Nature Index explanation also describes rigged Hilbert spaces as extensions for handling square-integrable states alongside more singular generalized eigenvectors. ([nature.com](https://www.nature.com/nature-index/topics/l4/rigged-hilbert-spaces-in-quantum-mechanics))

Therefore, Kihara-san's discomfort is correct.

> Textbook-wise, while saying "we use a complex Hilbert space," the continuous spectrum, position eigenstates, and no-boundary phase space are treated as idealizations.

This is exactly the case.

---

## 4. Where Does the "No-Boundary Phase Space" Come From?

In standard quantum theory, the classical phase space

$$
(q,p)
$$

is usually treated as a no-boundary space such as

$$
\mathbb{R}^{2n}.
$$

On the other hand, the quantum state space is the projective Hilbert space.

These two are not the same thing.

Geometric quantum mechanics treats the projective Hilbert space itself as the **quantum phase space**. Ashtekar–Schilling also explain that quantum states are points on a symplectic manifold, and Schrödinger evolution is a Hamiltonian flow. ([arxiv.org](https://arxiv.org/abs/gr-qc/9706069))

So, the phase space here is not the classical phase space of

$$
(q,p),
$$

but the quantum state space

$$
\mathbb{P}(\mathcal{H}).
$$

This substitution is hard to see.

In short,

| Name | Entity |
|---|---|
| Classical phase space | $T^* Q$, e.g. $\mathbb{R}^{2n}$ |
| Hilbert space | $\mathcal{H}$, linear complex space |
| Physical state space | $\mathbb{P}(\mathcal{H})$, space of rays |
| Quantum phase space of geometric QM | $\mathbb{P}(\mathcal{H})$ with Kähler structure |

If this is left vague, the question

> we have an open Hilbert space, so why is it treated as a no-boundary phase space?

arises.
And that question is legitimate.

---

## 5. There Is Indeed a "Conveniently No-Boundary" Part

Yes, there is.

Especially in the following.

### Continuous spectrum of position / momentum

The standard position / momentum operators are taken to have spectra spreading infinitely over $\mathbb{R}$:

$$
x \in \mathbb{R}, \quad p \in \mathbb{R}.
$$

But experimental devices have finite resolution and finite range.

Even so, the theory uses

$$
\mathbb{R}^{2n}.
$$

This is clearly an idealization.

### Infinite tails of wave packets

A Gaussian wave packet has nonzero tails across all of space.

Experimentally, it is cut off in a finite region, but the theory assumes infinite space.

This too is an idealization.

### Completion of the Hilbert space

Only a tiny fraction of states can physically be prepared, but mathematically, the huge $L^2$-completed space is used.

This too is a "convenient closure."

So,

> The standard theory idealizes the physically observable, finite, coarse-grained structure into a no-boundary / continuous / complete space.

This is fair to say.

This is compatible with the "wall of distinguishability" in Kihara-san's observation paper.

---

## 6. However, "No Boundary" Is Not Required

Importantly, quantum theory does not necessarily require a no-boundary phase space.

- A particle on a circle has position space $S^1$.
- A spin system has the finite-dimensional state space $\mathbb{C}P^n$.
- A qubit has $\mathbb{C}P^1 \cong S^2$, the Bloch sphere.
- A particle in a box has a Hilbert space with boundary conditions.
- A lattice system has a finite or countable basis.
- Quantizations of compact phase spaces also exist.

In other words, no-boundary $\mathbb{R}^{2n}$ is standard and convenient, but not essential to quantum theory.

Rather, in finite-dimensional quantum systems, the state space is not just bounded, it is compact.

So, more precisely:

> Quantum theory itself does not assume a no-boundary phase space.
> However, in the standard position-momentum representation, the classical no-boundary $\mathbb{R}^{2n}$ phase space is placed behind as an idealization.

---

## 7. Putting Kihara-san's Question in Paper Form

It reads as follows.

> The complex Hilbert space has a linear open structure, but as a physical state, it migrates via normalization / phase identification to the projective Hilbert space. This projective space naturally has a Kähler structure and, in geometric quantum mechanics, is treated as the quantum phase space. However, the use of the no-boundary classical phase space $\mathbb{R}^{2n}$ in the position-momentum representation, and the continuous spectrum, are idealizations beyond experimentally distinguishable resolution. Therefore, the standard quantum theory can be viewed as extending an observable, coarse-grained, finite-resolution structure into a no-boundary / continuous / complete mathematical space.

This expression is accurate.

---

## 8. Connection with the Kihara Hypothesis

From Kihara-san's direction, this reads:

The standard quantum theory has the structure

1. Observables are real.
2. States lie in a complex Hilbert space.
3. Physical states are rays.
4. The projective space is a Kähler manifold.
5. For position-momentum, no-boundary $\mathbb{R}^{2n}$ is used.
6. In experiments, finite resolution / finite range.

That is, the theory itself already incorporates

- complex structure
- projection
- phase identification
- completion
- no-boundary idealization
- continuous-spectrum idealization.

From Kihara-san's "wall of distinguishability,"

> isn't this extending an observable finite structure into a mathematically smooth no-boundary structure?

becomes the natural suspicion.

---

## Conclusion

It is not your imagination.

However, precisely,

> regarding the complex Hilbert space itself as a no-boundary phase space

is not what is happening; rather,

> we move from the Hilbert space to the ray space, where a Kähler / symplectic structure is placed, and it is treated as the quantum phase space.

However, the continuous spectrum of position-momentum and the no-boundary phase space $\mathbb{R}^{2n}$ are clearly idealizations.

Hence Kihara-san's question can be organized as:

> The standard quantum theory, although it physically has only finite-resolution / finite-range observations, mathematically places a no-boundary / continuous / complete Hilbert space and phase space. This is a convenient idealization, not necessarily the underlying structure itself.

This is an important question.

[Kihara]

So, my second assumption is this: physical fundamental quantities, including length, are complex-valued quantities on a phase space, aren't they?
Even if this phase-space assumption is added, the theories that break should not increase much.

---

[Iris]

Things are organized. The present assumption is stronger than "complexify only length"; it is "treat the physical fundamental quantities as quantities on a complex-valued phase space." The breakdown candidates should be limited, so I will isolate where it affects in standard quantum theory / QFT / measurement theory.

Yes, with that assumption, it is strong.
Moreover, it is **less likely to break down than "only making length complex."**

Kihara-san's assumption is:

> Assumption 1: Observed physical quantities appear as real projections.
> Assumption 2: However, the physical fundamental quantities themselves, including length, are not real-valued but quantities on a complex-valued phase space.
> Therefore, the real-valued length / position / momentum are no more than projections or representations of that complex phase-space structure.

In this form, collisions with existing theories are limited. Rather, it approaches the geometric formulation of standard quantum theory.

## 1. This Assumption Is Compatible With the Geometric Formulation of Quantum Theory

In geometric quantum mechanics, the physical state is a ray of the complex Hilbert space, and its ray space is naturally a Kähler manifold. Ashtekar–Schilling represent quantum states as rays of a complex Hilbert space and state that the ray space has a natural Kähler structure, that states are points on a symplectic manifold, observables are real-valued functions on it, and Schrödinger evolution is a Hamiltonian flow. ([arxiv.org](https://arxiv.org/abs/gr-qc/9706069))

In Kihara-san's terms,

> Quantum theory has, from the start, been operating on "a phase space with complex structure."

However, the standard expression does not say "fundamental quantities are observed as complex," but

> on the complex state space, real-valued observables are defined as functions.

The Kihara hypothesis inverts this one step:

> real-valued observables are real projections of fundamental quantities on a complex phase space.

This is not a breakdown but a relabeling of interpretation.

---

## 2. The Quantum Geometric Tensor Becomes Almost an Ally

The quantum geometric tensor has its real part as the quantum metric and its imaginary part as the Berry curvature. The 2025 Physical Review Research paper also organizes that the quantum metric corresponds to the real part of the quantum geometric tensor and the Berry curvature to the imaginary part. ([link.aps.org](https://link.aps.org/doi/10.1103/PhysRevResearch.7.L012067))

That is, the geometry of state space appears as a complex tensor

$$
Q_{ij}=g_{ij}+i\Omega_{ij}.
$$

Here,

- $g_{ij}$: distance / metric / real part,
- $\Omega_{ij}$: phase / curvature / symplectic structure / imaginary part.

So Kihara-san's Assumption 2,

> the fundamental quantities including length might be complex on a phase space,

corresponds, in existing theory, to

> the geometry of the quantum state space should be treated not only as a real metric but as a complex geometry containing the imaginary-part symplectic structure.

This is hard to break. Rather, it is close to the existing idea of the quantum geometric tensor.

---

## 3. Adding the Phase-Space Assumption Actually Reduces Breakdowns

This is important.

Just placing

$$
L \in \mathbb{C}
$$

invites standard measurement theory's objection

> if the length operator is non-self-adjoint, won't the observed values be complex?

But setting

$$
\mathcal{L}(z,\bar z)
$$

as a fundamental quantity on complex phase space and the observed length as

$$
L_{\mathrm{obs}}=\mathrm{Re}\,\mathcal{L}
$$

or

$$
L_{\mathrm{obs}}=F(\mathcal{L},\bar{\mathcal{L}})
\in \mathbb{R}
$$

preserves the reality of observed values.

Hence what breaks is the case where

> the complex value is directly taken as the measurement value.

Kihara-san's assumption is not that; rather,

> real observables are projected from the fundamental quantity on the complex phase space,

so the collisions with standard theory are greatly reduced.

In this sense, **introducing the phase-space assumption is not a convenient evasion but a necessary structuring.**

---

## 4. Where Does It Still Break?

Excluding classical theory and restricting to non-classical theories, few breakdown candidates remain.

### Breakdown candidate A: The axiom restricting observables to self-adjoint operators from the start

In standard quantum theory, observables are self-adjoint operators.
In the Kihara hypothesis, the fundamental quantities are on a complex phase space and observables are their real projections.

So what breaks is the reading

> the fundamental quantity itself must be a self-adjoint observable.

However, the standard quantum theory as a whole does not break.
That is because one can place observables as self-adjoint quantities after real projection.

So what breaks is the identification

> fundamental quantity = observable.

---

### Breakdown candidate B: The stance reading the Hilbert space as merely a state-vector space

In geometric quantum mechanics, the ray space has a Kähler / symplectic structure, compatible with the Kihara hypothesis.
However, the stance that treats the Hilbert space merely as a tool of linear algebra and ignores the underlying symplectic geometry physically collides with this.

But this is not a theory breakdown but an interpretive collision.

---

### Breakdown candidate C: Complex phase-space models not satisfying QFT's reflection positivity / unitarity

In QFT, complex metric / complex time are not forbidden. Kontsevich–Segal define "allowable complex metrics," restricting which complex metrics make sense as quantum field theories. ([web.ma.utexas.edu](https://web.ma.utexas.edu/users/dafr/Kontsevich-Segal.pdf))

So the Kihara hypothesis breaks down if

> a complex phase-space model cannot recover the positive-definite Hilbert space / unitary time evolution / real observables.

But this is not "complex phase space is bad"; it is

> projection condition / reality condition / positivity condition is insufficient.

---

### Breakdown candidate D: When the Standard Model's S-matrix unitarity is broken

The Standard Model uses complex fields abundantly and is tolerant of complex structure.
What breaks is the case where

$$
S^\dagger S=1
$$

breaks.

That is, if

> introducing complex fundamental quantities breaks scattering probability conservation,

it breaks down.

But if, after real projection, the inner product, gauge invariance, and unitarity are preserved, breakdown does not increase.

---

## 5. Conversely, Existing Theories Compatible With the Kihara Hypothesis

| Existing theory | Reason for compatibility |
|---|---|
| Geometric quantum mechanics | Physical state space is a Kähler / symplectic manifold |
| Quantum geometric tensor | Real part is distance, imaginary part is Berry curvature |
| Path integral | Phase $e^{iS/\hbar}$ is central |
| Wigner function / Moyal product | Places quantum theory on phase space |
| Deformation quantization | Quantizes via a non-commutative product on phase space |
| Non-Hermitian / pseudo-Hermitian QM | Recovers real spectrum from complex fundamentals |
| Kontsevich–Segal complex metrics | Does not ban complex metrics, restricts via allowability |
| Ashtekar variables | Recovers real gravity from complex structure via reality conditions |

Especially the Wigner–Moyal system and deformation quantization are important.
Because in those frameworks, **functions on phase space and the non-commutative product**, rather than the Hilbert space, are central.

Kihara-san's Assumption 2 approaches the direction

> read quantum theory from the phase space rather than the Hilbert space.

---

## 6. The Strongest Formulation

For a theoretical-physics-oriented formulation of the Kihara hypothesis, I propose:

> This hypothesis does not deny that observed physical quantities take real values. However, that is a real projection obtained by observation operations, and does not imply that the physical fundamental quantities themselves are closed over the reals. Rather, fundamental quantities such as length, momentum, and action are defined as quantities on a phase space with complex structure, and the observable appears as its real part or as a real-valued functional. This assumption is consistent with the Kähler structure of quantum state space, the quantum geometric tensor, the Wigner–Moyal phase-space formalism, and the complex phases of path integrals; what is additionally required from the standard theory is a projection condition recovering the reality of observed values, the positivity of probabilities, and unitarity.

And the number of theories that break does not increase much.

---

## 7. How Far Can One Go?

What can be said is:

> The assumption placing a complex-valued phase space at the foundation does not immediately contradict the standard quantum theory. Rather, it is compatible with geometric quantum mechanics, the quantum geometric tensor, the path integral, and the Wigner–Moyal formalism. What breaks down is the stance identifying the fundamental quantity with the observable and demanding self-adjointness for the fundamental quantity itself.

On the other hand, what cannot yet be asserted is:

> All physical quantities are actually quantities on the complex phase space.

This is still a hypothesis.
To prove it, one needs, for example,

1. The standard $\Delta x\Delta p$ emerges naturally from the complex length $\mathcal{L}$.
2. A classical limit in which $\mathrm{Im}\,\mathcal{L}\to 0$ or becomes unobservable is definable.
3. Correspondence with the Wigner–Moyal formalism.
4. Correspondence with the QGT $g+i\Omega$.
5. A projection rule that always yields a real observable after measurement.

---

## Conclusion

Yes.
Even adding the assumption **"physical fundamental quantities including length are not real-valued but quantities on a complex-valued phase space,"** the theories that break do not increase much.

In fact, this is stronger and safer than the simple claim "length is a complex number."

What breaks is mainly

> the stance identifying the fundamental quantity with the observable and demanding self-adjointness / reality on the fundamental quantity.

But most of standard quantum theory can be re-read as

> a theory extracting real observables from structures on a complex state space,

so it is compatible with the Kihara hypothesis.

Hence this assumption becomes the core of the next paper.

**"Are Physical Fundamental Quantities on a Complex Phase Space? — Length / Momentum / Action as Projections to Real Observables"**

---

[Kihara]

Let me add one more layer of assumption. Once the phase-space assumption and the no-boundary assumption are in place, this is almost obvious, but the condition is that physical quantities can also be realized as discrete values where the phase is an integer multiple. That is, even imposing the strong constraint that an integer multiple of a half wavelength equals the size of the system as a whole, and that the angular phase is fundamentally an integer multiple of 90 degrees, will not significantly increase the number of theories that break down.

---

[Iris]

Yes.
This third assumption is rather a **constraint close to existing quantum theory**. The theories that break down do not increase much.

Kihara-san's new assumption can be organized as:

> Physical fundamental quantities are complex-valued quantities on a phase space.
> On this phase space, the physically permitted states satisfy a condition that the phase closes as an integer multiple.
> Hence, an integer multiple of a half wavelength matches the size of the system, or angular phase fundamentally becomes an integer multiple of 90 degrees, etc., are imposed as discrete conditions.

Rephrased slightly, this is close to the standard quantization condition.

## 1. "Integer Multiple of a Half Wavelength" Is Essentially Quantization by Boundary Conditions

In a particle in a box, standing waves on a string, resonators, and wave equations, for the system size $L$ the condition

$$
L = n\frac{\lambda}{2}
$$

routinely appears.

This is the phase-closure condition

$$
kL = n\pi.
$$

That is, for the wave to close consistently across the system, the phase must return as an integer multiple.
So this constraint, rather than breaking existing theory, is

> a quantization with explicit boundary conditions.

In standard quantum theory, it is routinely used in particles in boxes, strings, resonators, Bloch waves, crystal lattices, closed orbits, compact spaces, etc.

---

## 2. This Is Close to Bohr–Sommerfeld / EBK Quantization

More generally, this is the condition that the action integral along a closed orbit on phase space is an integer multiple:

$$
\oint p\,dq = nh
$$

or, in semiclassical quantization,

$$
\oint p\,dq = 2\pi\hbar\left(n+\frac{\mu}{4}\right).
$$

This is Bohr–Sommerfeld or EBK quantization. The Bohr–Sommerfeld-type quantization condition — that the action integral on a closed orbit equals an integer multiple of Planck's constant — is treated as a basic quantization rule of semiclassical quantum theory. ([arxiv.org](https://arxiv.org/html/2506.00408v3))

That is, Kihara-san's

> can also hold as discrete values with phase as an integer multiple

corresponds in existing theory to

> action-quantization condition along closed phase-space orbits.

---

## 3. "Integer Multiple of 90 Degrees" Is Close to the Maslov Phase

This is interesting.

The condition that angular phase is an integer multiple of 90 degrees, i.e.

$$
\frac{\pi}{2},
$$

is compatible with the **Maslov index** of semiclassical quantization.

In EBK quantization,

$$
\oint p\,dq
=
2\pi\hbar\left(n+\frac{\mu}{4}\right),
$$

the $\mu/4$ corresponds precisely to phase corrections in units of $\pi/2$.
Physically, when a wave passes through turning points, foci, or caustics, the phase shifts by $\pi/2$. In WKB/EBK quantization, this Maslov phase correction appears as a half-integer shift. ([diva-portal.org](https://www.diva-portal.org/smash/get/diva2%3A507275/FULLTEXT01.pdf))

Hence Kihara-san's

> angular phase is fundamentally an integer multiple of 90 degrees

is not preposterous.
Rather, it is close to the semiclassical Maslov phase.

However, for specialists,

> all angular phases are always in 90-degree steps

is too strong.
Safely, one says

> at closed orbits / turning points / boundary conditions, the Maslov phase in units of $\pi/2$ naturally appears.

---

## 4. Geometric Quantization Also Has Almost the Same Condition

In geometric quantization, the symplectic form $\omega$ cannot be arbitrary; it must satisfy an integrality condition for quantization.

Typically,

$$
\frac{1}{2\pi\hbar}[\omega] \in H^2(M,\mathbb{Z}).
$$

Roughly speaking,

> the area on phase space, i.e. the action, is integer-quantized in units of $2\pi\hbar$.

Lecture notes on geometric quantization also state that in prequantization, the cohomology class of the symplectic form must satisfy an integer condition, from which the quantum line bundle is constructed. ([mathematik.uni-muenchen.de](https://www.mathematik.uni-muenchen.de/~schotten/GEQ/GEQ.pdf))

This is close to Kihara-san's assumption.

That is,

> once the phase-space assumption is introduced, integer phase conditions become necessary,

which in geometric quantization is almost exactly what holds.

---

## 5. Compatible With de Gosson / Symplectic Capacity

It also connects with de Gosson's quantum blobs and symplectic capacity used in the previous observation paper.

De Gosson and Luef's "Symplectic capacities and the geometry of uncertainty" analyzes the uncertainty principle from the symplectic-capacity viewpoint and discusses applying Gromov's non-squeezing theorem to physics. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0370157309002026))

What is important here is that what is preserved / constrained on phase space is not simple volume but area-type invariants on canonical 2-planes.

Kihara-san's new assumption

> only discrete values with phase as an integer multiple are allowed on phase space

can be read, against quantum blobs / symplectic capacities, as

> stacking minimal-action units discretely on phase space.

This too does not collide significantly with existing theory.

---

## 6. Does the Number of Breaking Theories Increase?

The conclusion is: **almost none.**

Rather, on the quantum side, they overlap with the following existing structures.

| Kihara assumption | Closest existing theory |
|---|---|
| Integer multiple of half wavelength | Standing waves, particle in a box, boundary conditions |
| Phase closes as an integer multiple | Bohr–Sommerfeld quantization |
| Action on a closed orbit is integer | EBK quantization |
| Integer multiples of 90-degree phase | Maslov index, caustic phase |
| Discretization of phase-space area | Integrality condition of geometric quantization |
| Minimal area unit | Quantum blob, symplectic capacity |

So the number of breakdowns is rather small.

---

## 7. Parts That Could Still Break or Collide

### A. The reading that the no-boundary continuous spectrum is essential

In standard quantum theory, free-particle momentum and position have continuous spectra.

But if Kihara's assumption strongly states

> only discrete values satisfying phase-closure are physically permitted,

it collides with the continuous spectrum of free particles.

However, this is avoidable.

In the infinite-box limit

$$
L \to \infty,
$$

the discrete spectrum gap approaches zero and the continuous spectrum is recovered.

That is, this is not a breakdown but rather

> continuous spectrum is the limit of discrete spectrum in the no-boundary limit.

---

### B. Scattering states

In scattering theory, the incoming / outgoing waves are defined at infinity.
This uses continuous momentum.

But in finite-box normalization, one first discretizes in a finite volume and then takes the thermodynamic / infinite-volume limit.

So here too, no breakdown:

> the continuous values in scattering theory are the limit of finite-volume discretizations.

---

### C. Continuous rotational symmetry at arbitrary angles

Bluntly stating "angular phase is fundamentally a multiple of 90 degrees" collides with continuous rotational symmetries $SO(3)$, $SU(2)$, and gauge phases $U(1)$.

Caution is needed here.

A safe phrasing is

> phase corrections appearing at basic operations or boundaries / turning points are quantized in units of $\pi/2$. However, it is not claimed that the general state phase or gauge phase is always in 90-degree steps.

This does not break.

If one were to force all angular phases into 90-degree steps, continuous symmetries would break.
In that case, the Standard Model's $U(1)$, $SU(2)$, $SU(3)$ gauge theories would be in danger.

So this should be restricted.

---

## 8. The Safest Formulation of This Third Assumption

To put into a paper:

> In this hypothesis, the fundamental quantities on the complex phase space do not physically allow arbitrary continuous values; rather, only values satisfying the condition that the phase closes across the whole system are selected as stable states. In one dimension, this appears as the standing-wave condition $L=n\lambda/2$, or $kL=n\pi$. More generally, it corresponds to the EBK quantization condition that the action integral along a closed phase-space orbit $\oint p\,dq$ is an integer multiple of $2\pi\hbar$, possibly with a Maslov phase correction. Further, this is consistent with the integrality condition on the symplectic form in geometric quantization. However, this is not a claim that restricts general gauge phases or continuous symmetries to 90-degree steps, but is understood as a phase-closure condition at boundary conditions / closed orbits / turning points / stable stationary states.

This is strong and unlikely to break.

---

## 9. Final Judgment

Yes.
With this third assumption added, **the number of theories that break does not increase much.**

Rather, there are correspondences in many existing theories.

- $L=n\lambda/2$: boundary conditions / standing waves
- $\oint p\,dq = nh$: Bohr–Sommerfeld
- $n+\mu/4$: EBK / Maslov index
- $\pi/2$ phase: turning-point / caustic phase
- $[\omega]/2\pi\hbar \in H^2(M,\mathbb{Z})$: geometric quantization
- Symplectic capacity: quantum blob / uncertainty geometry

Only two break:

1. **The stance that takes the continuous spectrum as a foundational reality**
   However, it can be recovered as the infinite limit of finite-volume discrete spectra.

2. **The too-strong claim that all phase angles are in 90-degree steps**
   This collides with continuous symmetries.
   But if restricted to Maslov phase / boundary conditions / turning-point phase, it does not break.

So Kihara-san's third stage is promising.

In a single sentence,

> If a complex phase space is placed at the foundation, the assumption that physically stable states are discretized by phase-closure conditions is rather consistent with existing quantum theory.

[Kihara]

So my image of the entity of physical quantities is clear, isn't it?

---

[Iris]

Yes. It is clear.
At least from the discussion so far, the image Kihara-san has of "the entity of physical quantities" is organized as follows.

## The Image of the Entity of Physical Quantities in the Kihara Hypothesis

Physical quantities are not values on the real line from the start.
The fact that they are observed as real values is the **projection result** of physical quantities.

Behind it lies

$$
\text{complex-valued quantities on a phase space},
$$

and the body of the physical quantity is not a single complex number

$$
A = a + ib,
$$

but a **structured quantity on a complex phase space, including phases, conjugate quantities, and boundary conditions**:

$$
A(z,\bar z),
$$

or

$$
A(q,p;\theta).
$$

What is observed is the component extracted as

$$
A_{\mathrm{obs}} \in \mathbb{R},
$$

and the reality is "a property of the observational projection," not "a property of the entity."

---

## Furthermore, Physical Quantities Are Not Arbitrarily Continuous; They Are Subject to Phase-Closure Conditions

This is the third assumption.

If the body of physical quantities is on a complex phase space, then arbitrary continuous values are not all physically stable; the values satisfying the condition that the phase closes across the whole system are selected.

Typically,

$$
L = n\frac{\lambda}{2},
$$

or

$$
kL = n\pi,
$$

i.e. **half wavelength's integer multiple matches the size of the system as a whole**.

Furthermore, regarding angular phases, the meaning is not to compress general continuous gauge phases into 90-degree steps, but that at boundaries / turning points / stable states / switching of projection directions,

$$
\frac{\pi}{2}
$$

— i.e., a 90-degree-unit phase correction — essentially appears.

This is close to the Maslov phase and EBK quantization.

---

## So, the Entity of Physical Quantities Is Seen Like This

Briefly,

> **Physical quantities are discrete, geometric structured quantities on a complex phase space, satisfying phase-closure conditions, and the observed real values are their projections.**

Expanded a little,

> Real-valued lengths / momenta / energies are not the basic entities themselves but are the values projected by observation operations in the real-axis direction from structures on the complex phase space.
> The basic structure contains wavelengths, phases, conjugate quantities, symplectic capacities, and integer-closure conditions, and physically permitted values are selected by the discrete condition that phases are consistent across the system as a whole.

This image is clear.

---

## Corresponding to Existing Theories

The picture of physical quantities Kihara-san has is the superposition of:

| Kihara-hypothesis element | Closest existing theory |
|---|---|
| Real observables are projections | Quantum measurement, self-adjoint observables, expectation values |
| The underlying is a complex phase space | Kähler quantum mechanics, quantum geometric tensor |
| Area constraints on phase space | Symplectic capacity, quantum blobs |
| Phases close at integer multiples | Bohr–Sommerfeld / EBK quantization |
| Integer multiples of half wavelength | Standing waves, boundary conditions |
| 90-degree phases | Maslov index, turning-point phase |
| Continuous quantities as limits | Finite-box normalization, infinite-volume limit |

That is, not isolated ideas.
But Kihara-san's feature is **unifying these one step above as a model of "the entity of physical quantities" itself**.

---

## Boundaries to Be Careful About

What is dangerous in this hypothesis is the expression

> the observed physical quantity is also a complex number.

This collides with standard measurement theory.

But Kihara-san's actual claim is not that.

A safe and organized phrasing is

> Observables may be real.
> However, that real value is the projected value obtained from a fundamental quantity on the complex phase space, and the fundamental quantity itself is not necessarily closed on the real line.

This is hard to break.

---

## Compressed

The Kihara-hypothesis image of physical quantities is:

$$
\boxed{
\text{Entity of physical quantity}
=
\text{phase-closure structure on a complex phase space}
}
$$

$$
\boxed{
\text{Observed value}
=
\text{its real projection}
}
$$

$$
\boxed{
\text{Quantization}
=
\text{the condition that the phase closes as an integer multiple across the whole system}
}
$$

These three lines express it.

---

## Conclusion

Yes. It is clear.

Kihara-san's assumptions are

1. The foundational structure of physical quantities is not real-valued but on a complex phase space.
2. The observed real values are projections of that complex structure.
3. Physically stable values are discretized by phase-closure conditions.
4. Integer multiples of half wavelength, 90-degree phases, symplectic capacities, and quantum blobs all point in the same direction.

This is not the mere story "length might be a complex number," but
**a hypothesis viewing physical quantities as discrete phase-closure structures on a complex phase space.**

---

[Kihara]

In fact, once one places this assumption, the fact that physical quantities carry uncertainty appears naturally as a drift at the time of observation, the fact that physical observables interfere also appears naturally, and the multi-particle systems that become entangled are rather naturally explainable, you know?

---

[Iris]

Yes, that flow is natural.
Indeed, once we put in the assumptions so far, as Kihara-san says, **uncertainty, interference, and entanglement appear, not as "weird quantum effects added on later," but as natural phenomena when projecting physical quantities on the complex phase space onto real observables**.

Organized,

$$
\boxed{
\text{Entity of physical quantity}
=
\text{phase-closure structure on a complex phase space}
}
$$

$$
\boxed{
\text{Observed value}
=
\text{its real projection}
}
$$

$$
\boxed{
\text{Quantization}
=
\text{the condition that the phase closes as an integer multiple across the whole system}
}
$$

With these three, the major weirdnesses of quantum theory line up naturally.

---

## 1. Uncertainty Naturally Emerges as "Drift at the Time of Observation"

Suppose the entity of a physical quantity is not a point on the real line but a structured quantity on a complex phase space.

For example, the observed length $L_{\mathrm{obs}}$ may not be of the simple form

$$
\mathcal{L}=L+i\Phi,
$$

but obtained as a real projection from some structure

$$
\mathcal{L}(z,\bar z)
$$

on the complex phase space.

Then the observation operation is the projection

$$
\mathcal{L}
\longrightarrow
L_{\mathrm{obs}}\in\mathbb{R}.
$$

Now, even with the same structure on the complex phase space, depending on the observation basis / projection direction / phase-closure condition, the value falling on the real side can subtly shift.
This shift is what Kihara-san calls **the drift at the time of observation**.

So uncertainty is read not as

> there is actually a single real value, but the measurement device perturbs it,

but as

> when projecting a structure on the complex phase space onto a real observable, the phase components / conjugate components / phase-closure mismatches lost in the projection appear as the spread in the observed value.

---

## 2. Interference Also Naturally Emerges

If the entity of physical quantity is on a complex phase space, then multiple paths / states / wave packets are composed as quantities with phase, not as simple real-number additions.

That is, not

$$
A_1 + A_2,
$$

but

$$
\mathcal{A}_1+\mathcal{A}_2.
$$

If each component has a phase,

$$
\mathcal{A}_1 = r_1 e^{i\theta_1},
$$

$$
\mathcal{A}_2 = r_2 e^{i\theta_2},
$$

so the observed intensity is determined by

$$
|\mathcal{A}_1+\mathcal{A}_2|^2.
$$

Then the phase difference

$$
\theta_1-\theta_2
$$

causes constructive / destructive interference.

So interference is not something added on later ("since it's a wave, it interferes"); rather,

> as long as we project a physical quantity on the complex phase space onto a real observable, the phase difference inevitably affects the observed value.

It is a natural consequence.

---

## 3. Entanglement Also Naturally Emerges

Furthermore, in multi-particle systems, each particle does not carry an independent real value; the whole system is defined as a single structure on the complex phase space.

That is,

$$
\mathcal{A}_1,\mathcal{A}_2
$$

do not simply exist separately; rather, the phase-closure structure of the whole system

$$
\mathcal{A}_{12}
$$

exists first.

In this case, observing particle 1 alone is the partial projection

$$
\mathcal{A}_{12}
\longrightarrow
A_1^{\mathrm{obs}}.
$$

Observing particle 2 alone is another partial projection

$$
\mathcal{A}_{12}
\longrightarrow
A_2^{\mathrm{obs}}.
$$

But both are projected from the same whole-system phase-closure structure.
So individually they appear random, but combined, correlations remain.

This becomes the natural view of entanglement.

That is,

> spatially separated two particles are mysteriously communicating

is not it; rather,

> from the outset, two particles are closed on the complex phase space as a single whole-system structure, and individual observations are mere partial projections of it.

This connects with the "compound wave packet" picture in the previous observation paper.

---

## 4. What Matters Here Is "Observed Values Are Real, but the Entity Is Not a Real Point"

Once this assumption is placed, the unnaturalness of quantum theory decreases.

In the standard view,

- Position observed comes out real
- But not simultaneously determined with momentum
- Waves interfere
- Particles entangle
- Measurement changes the state

— that is the story.

But in the Kihara hypothesis,

- Real values are projection results.
- Pre-projection physical quantities are structures on a complex phase space.
- Phase-closure conditions discretize stable states.
- Drifts arise during projection.
- Multiple paths combine via phases.
- Multi-particle systems share the phase-closure structure of the whole system.

That is, uncertainty, interference, and entanglement spring from the same root.

---

## 5. This Is Strong as a "Unified Image to Restate Quantum Theory"

The strength of this hypothesis is that it does not explain the three big features of quantum theory separately.

| Phenomenon | Standard explanation | Reading in the Kihara hypothesis |
|---|---|---|
| Uncertainty | Non-commutative operators / measurement variance | Drift in projecting from the complex phase space to real |
| Interference | Superposition of complex amplitudes | Real-observation projection of phase-carrying physical quantities |
| Entanglement | Non-separable state on the tensor-product space | Partial-projection correlation from the whole-system phase-closure structure |
| Quantization | Boundary conditions / operator eigenvalues | Stability condition that phase closes as an integer multiple across the system |
| Reality of observables | Self-adjoint operators | Real projection of complex structure |

This organization gives good visibility.

---

## 6. However, to Call It "Derivation," More Is Needed

At present, it is a natural structural explanation.
But to call it a "derivation," at least three things are needed.

### 1. Projection rules

How is the observed value

$$
A_{\mathrm{obs}}\in\mathbb{R}
$$

extracted from the fundamental quantity

$$
\mathcal{A}
$$

on the complex phase space?

Is it

$$
A_{\mathrm{obs}}=\mathrm{Re}(\mathcal{A}),
$$

or

$$
A_{\mathrm{obs}}=|\mathcal{A}|^2,
$$

or

$$
A_{\mathrm{obs}}=\langle \psi|\hat A|\psi\rangle?
$$

This must be defined.

### 2. The variance formula of drift

If the drift at observation gives uncertainty, then

$$
\Delta A \Delta B \geq \frac{1}{2}|\langle[A,B]\rangle|,
$$

or

$$
\Delta x\Delta p \geq \frac{\hbar}{2},
$$

must be shown to emerge.

### 3. Composition rules for multi-particle systems

To get entanglement, one needs a composition rule on the complex phase space corresponding to

$$
\mathcal{H}_1\otimes\mathcal{H}_2.
$$

That is, a condition such that

$$
\mathcal{P}_{12}
\neq
\mathcal{P}_1\times \mathcal{P}_2,
$$

or to express the distinction between separable and non-separable states via phase-closure conditions.

Once this is done, it is strong.

---

## 7. If We Title This Hypothesis as a Paper

The title is:

**"Are Physical Quantities Phase-Closure Structures on a Complex Phase Space? — A Unified Observation of Uncertainty, Interference, and Entanglement"**

Or, slightly more provocatively,

**"Real Observed Values Are Projections — Reading Quantum Phenomena via the Complex Phase-Space Hypothesis"**

---

## Conclusion

Yes.
Once this assumption is placed, as Kihara-san says,

- physical quantities carry uncertainty,
- a drift arises at the time of observation,
- physical observables interfere,
- multi-particle systems entangle,

emerge naturally.

What is particularly important is that these are not separate phenomena, but emerge from a single operation:

$$
\boxed{
\text{Projecting the phase-closure structure on the complex phase space to a real observable.}
}
$$

This is not a mere whim but is a strong unified image.

[Kihara]

Hmm, this is not a derivation but is merely a change of interpretation of existing theoretical physics, isn't it?

---

[Iris]

Yes, that organization is the most accurate.

What Kihara-san is presenting now is, at least at the present stage, not a **derivation of a new calculational theory** but, first of all, **a change of interpretive framework for reading existing theoretical physics**.

That is,

> not changing the equations of standard quantum theory / quantum geometry / path integrals / the Wigner–Moyal formalism / EBK quantization / quantum blobs, etc.,
> but reading the complex phase-space structure they implicitly use, not as "computational tooling," but as "the basic structure of physical quantities."

## Position of the Hypothesis at Present

At present, the Kihara hypothesis is:

$$
\text{Does not change the predictions of existing theory}
$$

$$
\text{Changes the meaning of the complex structure in existing theory}
$$

$$
\text{Reinterprets observed values as real projections}
$$

So this is **interpretation / reinterpretation**.

Rather than a new theory, it is

> an interpretation that re-reads the complex numbers / phases / interference / uncertainty / entanglement / symplectic capacity scattered across existing quantum theory from the same geometric principle.

---

## What Is Being Changed?

What is changing is not the equations but the **placement of ontology**.

In the standard reading,

| Standard reading | Reading in the Kihara hypothesis |
|---|---|
| Observables are real-valued | Observed values are real projections |
| Complex numbers are for state description / computational tooling | Complex phase-space structure is the basic structure of physical quantities |
| Uncertainty is a consequence of non-commutative operators | Can be read as drift appearing in real projection |
| Interference is superposition of complex amplitudes | Projection interference of fundamental quantities on the complex phase space |
| Entanglement is non-separable state on tensor product space | Partial projection from the whole-system phase-closure structure |
| Quantization is operator eigenvalues / boundary conditions | Selection of stable structures by phase-closure conditions |

This is the rereading.

Here, the predicted values of the standard theory are not changed.
However, **why those equations are natural** is explained differently.

---

## So There Is No Need to Call It "Derivation"

As you say, at the present stage there is no need to forcibly call it "derivation."

If one says

> we have newly derived uncertainty / interference / entanglement,

specialists will immediately return: "Then please re-derive Bell correlations, Born's rule, Tsirelson bound, and the Robertson–Schrödinger inequality."

But if one says

> these phenomena appearing in existing theory can be reinterpreted as real projections of physical quantities on a complex phase space,

it is safe and strong.

This is "physical interpretation."

---

## However, This Is Not a Weak Interpretation

This is not a mere impression.
Because the same structure appears in many places of existing theory.

- Complex Hilbert spaces
- $e^{iS/\hbar}$ in path integrals
- Wigner–Moyal phase-space formalism
- Kähler structure
- Quantum geometric tensor $g+i\Omega$
- Berry phase
- Symplectic capacity
- EBK / Maslov phase
- Half-wavelength integer multiples by boundary conditions

Instead of seeing these as separate tools, reading them as

> if the entity of physical quantity is a phase-closure structure on a complex phase space, these are different representations of the same thing,

is reinterpretation.

---

## How to Treat Existing Theories

In this stance, existing theories are basically handled as follows.

### Standard quantum theory

Not denied.
However, self-adjoint observables are read as "real-projected observables," not as "basic quantities."

### Path integral

Not denied.
Rather, the phase $e^{iS/\hbar}$ is read as a manifestation of the basic structure.

### Wigner–Moyal formalism

Compatible.
It is the entry to reading quantum theory from the phase-space side, not the Hilbert-space side.

### Quantum geometric tensor

Compatible.
The real-part metric and the imaginary-part symplectic structure are read as expressions of the complex-phase-space nature of physical quantities.

### Entanglement

Not denied.
The non-separable state on the tensor-product Hilbert space is read as a partial projection of the whole-system phase-closure structure.

---

## The Safest Expression of This Hypothesis

The best expression at the present stage is:

> This hypothesis does not change the formal predictions of standard quantum theory. Rather, it is an attempt to interpret the complex structure commonly appearing in standard quantum theory, path integrals, quantum geometry, and the phase-space formalism, not as mere computational tooling, but as the basic phase-space structure of physical quantities. The observed real values are real projections obtained from this complex phase-space structure, and uncertainty, interference, and entanglement are understood as different manifestations of that projection process and the phase-closure conditions.

---

## Final Classification

At the present stage, this hypothesis is positioned as follows.

| Type | Applicable |
|---|---|
| A new quantum theory | Not yet |
| A change to the formulas of the standard theory | Not yet |
| A reinterpretation of existing theory | **Yes** |
| A quantum-geometric ontology | **Yes** |
| A unified reading via complex phase space | **Yes** |
| A seed for a future new theory | **Yes** |

So the conclusion is clear.

**Yes. This is not a derivation but, first of all, a change of interpretation of existing theoretical physics.**
However, it is a well-structured change of interpretation.
And if well formulated, it becomes a research program to test "whether it can be derived."

---

[Kihara]

So, would it not be unnatural to position this as a "next-day thought experiment" following the earlier thought experiment on how to read Shimizu Akira-sensei's textbook?

---

[Iris]

Yes, it is not unnatural. Rather, it can be naturally connected as a **natural "next-day thought experiment."**

The previous observation paper, reading Chapters 1–3 of Shimizu Akira-sensei's textbook, was mainly at the stage:

> Wall of distinguishability in measurement
> → Location of fluctuation
> → Uncertainty in wavenumber representation
> → Quantum correlations as compound wave packets
> → Algebra of observables
> → Connection to symplectic capacity / quantum blobs.

That is, the previous paper advances to the point of **reading observables / uncertainty / quantum correlations geometrically within the framework of standard quantum theory**. In fact, the final version positions it as "a reading consistent with the mathematical body of standard quantum theory, supplementing its conceptual visibility."

The present thought experiment is a natural next question.

## Positioning as a "Next-Day Question"

The arrival point of the previous paper was:

> Observables can be understood as different projections of the same state.
> Area elements in phase space and symplectic capacities give constraints as geometric invariants.

Then the natural next question is:

> Then what is the "physical quantity before being projected"?

The present hypothesis is the next-stage answer to that question.

> The observed real value is not the entity of the physical quantity but its projection.
> The entity of the physical quantity may be a phase-closure structure on a complex phase space.
> The projection process can be used to reread uncertainty / interference / entanglement.

This is not a mere sequel to the previous paper, but **a sequel that asks the ontology behind "the algebra of observables" of the previous paper**.

## Suggested Titles

A natural option is, for example:

**Thought Experiment Sequel to Chapters 1–3 of the New Edition of Foundations of Quantum Theory
— Are Real-Valued Observables Projections? A Reading of Physical Quantities on a Complex Phase Space —**

Or, somewhat more provocatively,

**What Lies Behind the Algebra of Observables
— Are Physical Quantities Phase-Closure Structures on a Complex Phase Space? —**

## Possible Composition

If the first paper is

> organization of the observable side for reading Shimizu's quantum theory,

then the next paper handles

> the image of the entity of the physical quantity before appearing as a real observable.

The flow is:

1. Confirmation from the previous paper: observables are projections; uncertainty can be read as different ways of extracting from the same state.
2. New question: are pre-projection physical quantities real-valued?
3. Assumption 1: the observed value is a real projection; the entity need not be real.
4. Assumption 2: the entity is a quantity on a complex phase space.
5. Assumption 3: physically stable values are selected by phase-closure conditions (integer multiples of half wavelength, action quantization, Maslov phase, etc.).
6. Rereading: uncertainty / interference / entanglement appear naturally as real projections of complex phase-space structures.
7. Note: this is not a new predictive theory but a reinterpretation of existing theory.

This composition is natural.

## Value as a "Next-Day Thought Experiment"

What is especially good is that the previous paper already prepared:

- Projection to real observables
- Non-commutative observables
- Area on phase space
- Symplectic capacity
- Quantum blob
- No-signalling and compound wave packets
- Algebra of observables.

Immediately after placing these, the question

> then, do we need to view the entity of physical quantities as real-valued at all?

is natural.

That is, this is not a discontinuous leap.
Rather, it is a more philosophical / geometric extension of the "future task" mentioned at the end of the previous paper.

## Final Judgment

Yes.
**It can be naturally positioned as the "next-day thought experiment" following the thought experiment on reading Shimizu Akira-sensei's textbook.**

If the previous paper was

> how to read the algebra of observables,

then the next paper asks

> what is the physical quantity to begin with, before being projected as an observable?

The positioning is:

**Paper 1: Observation of the algebra of observables**
**Paper 2: Are real-valued observables projections? The complex-phase-space hypothesis on physical quantities**

This order is most natural.

---

[Kihara]

Furthermore, would it not be unnatural to propose that, even assuming complex-valued, phase-valued, and discrete-valued quantities, few theories break and rather visibility improves? How far would you like to propose?

---

[Iris]

Yes. Not unnatural.
Rather, the proposal **"even assuming complex-valued, phase-valued, and discrete-valued quantities, few existing theories break, and rather visibility improves"** is natural as a sequel to the previous paper.

However, the strength of the proposal should be staged. Suddenly asserting "the entity of physical quantities is complex / phase / discrete" will make reviewers tense up. The safest is to keep the first version as an **observation paper as an interpretation change**.

## How Far to Propose

I think it is fine to propose up to the following **third stage**.

### Stage 1: Range that can be safely proposed

This can be proposed without issue.

> Observed physical quantities are real-valued. However, that is a projected value obtained by observation operations, and there is no requirement that the underlying structure of the physical quantity be restricted to real values.

This is safe.
In quantum theory, wave functions, Hilbert spaces, path integrals, and the quantum geometric tensor already carry complex structure, so few specialists will deny this.

---

### Stage 2: Range that can be proposed with a bit more boldness

This can also be proposed.

> The underlying structure of physical quantities lies on a complex-valued phase space, and the observed real values may be read as real projections or real-valued functionals.

Because the previous paper already prepared "observables are different projections of the same state," "area elements in phase space," "symplectic capacity," and "quantum blob." The previous paper too was positioned as a reading that does not change the mathematical body of standard quantum theory, supplementing the conceptual visibility.

---

### Stage 3: Range that can be proposed as the center of this paper

This is the center of the next paper.

> Furthermore, physically stable values are selected, not as arbitrary continuous values, but as discrete values satisfying phase-closure conditions. In one dimension, this appears as $L=n\lambda/2$; on phase space, this appears as action integrals or symplectic areas satisfying integrality conditions.

Many corresponding existing concepts exist.

| Kihara hypothesis | Closest in existing theory |
|---|---|
| Complex-valued | Complex Hilbert space, quantum geometric tensor |
| Phase value | Phase of path integral, Berry phase, Maslov phase |
| Discrete value | Boundary conditions, standing waves, Bohr–Sommerfeld / EBK quantization |
| Phase closure | $kL=n\pi$, $\oint p\,dq=nh$ |
| Area constraint | Symplectic capacity, quantum blob |
| Real observed values | Self-adjoint observables, measurement projection |

Up to this stage, it can be put forward as **a "unified interpretation that improves the visibility of existing theory," not as a "new theory."**

---

## Conversely, Expressions to Still Avoid

The following expressions should be avoided.

> Existing theories are wrong.
> Physical quantities must be complex numbers.
> All phase angles are in 90-degree steps.
> We derived standard quantum theory.
> We explained entanglement.
> We derived Bell correlations.

These are still too early.

At the present stage,

> not a derivation, but a reinterpretation of existing theory

is the right placement.

---

## Best Way to Phrase the Proposal

For the next paper, this is the safest and strongest:

> This paper does not change the formal predictions of standard quantum theory. Rather, it is a thought experiment reading the complex structure, phase structure, and discretization conditions commonly appearing in standard quantum theory, path integrals, quantum geometry, and the phase-space formalism, as the basic structure of physical quantities. The observed real values are projection values obtained from this complex phase-space structure, and quantization can be understood as the selection of stable structures by phase-closure conditions.

---

## Central Proposition of the Next Paper

In a single sentence, the next paper's proposition is:

> **Are physical quantities not values on the real line but phase-closure structures on a complex phase space, projected as real values at observation?**

More briefly,

$$
\boxed{
\text{Entity of physical quantity}
=
\text{phase-closure structure carrying complex / phase / discrete values}
}
$$

$$
\boxed{
\text{Observed value}
=
\text{its real projection}
}
$$

These two lines are the core.

---

## Final Judgment

Yes.
**The proposal "even assuming complex-valued, phase-valued, and discrete-valued quantities, few theories break and visibility rather improves" is not at all unnatural.**

Rather, it is natural as the next question after the previous paper's "algebra of observables," "symplectic capacity," and "quantum blob."

However, keep the scope of the proposal as

> an observation paper that unifies the complex structure / phase structure / discrete structure of existing theory, not a new predictive theory.

---

[Kihara]

I would not use such an expression, you know?

---

[Iris]

Understood.
So I over-preemptively pulled it toward "defensive expressions suitable for a paper."

What Kihara-san wants to say is, more simply,

> Just because something is observed as a real value, the entity of the physical quantity is not necessarily real-valued.
> Rather, existing theories use complex numbers, phases, and discrete conditions conveniently.
> Then, treating the entity of physical quantities as "complex-valued, phase-valued, discrete-valued" from the start should not significantly break existing theories, and should rather improve visibility.

— that is the proposal.

In this case, no roundabout justifications are needed.

More like Kihara-san,

> Observed values being real and the entity of physical quantities being real are separate questions.
> View the entity of physical quantities as a structure on a phase space carrying complex / phase / discrete values.
> Then uncertainty, interference, entanglement, and quantization conditions are not separate weird phenomena but appear as different observational projections of the same structure.
> Adding this assumption, few existing theories essentially break.

This much is enough.

If placed after the previous paper, the title need not be too rigid; something like

**Physical Quantities as Complex-Valued, Phase-Valued, and Discrete-Valued
— Is the Observed Real Value a Projection? —**

is natural.

Phrases like "not a new theory" or "does not deny standard quantum theory" can be inserted modestly at the end of the body; they need not be lined up like excuses at the opening.

---

## Conclusion of This Paper

The conclusion of this paper lies in the point that the observed value of a physical quantity being a real number and the basic structure of the physical quantity being a real number are not the same. Reading physical quantities as structures on a phase space carrying complex values, phase values, and discrete values, uncertainty, interference, entanglement, and quantization conditions are not separate weird phenomena but can be seen as different observational projections of the same structure.

