# Zero Closure Was Four-Dimensional — Central Projection Survives in the Complex World

## How to compute a universe without solving the many-body problem

My series of papers has always started from a single equation.

x² + y² + z² = R²

Treating everything as real, this is the equation of an arbitrary point on a sphere of radius R. Add variables and it becomes a hypersphere in n+1 dimensions.

This form carries an enormous advantage.

**As long as you move on the sphere, you can compute the motion without solving the many-body problem.**

Why? Because there is only one condition to keep. "Do not change R" — that is all. Beyond that you may move on the sphere however you like. There is no need to chase hundreds or thousands of interactions one by one. This is what central projection means, and it is why I started here.

---

## Extending to complex numbers, and reaching a table of 62 particles

From there the series extended xₙ to complex numbers and took

Σ xₙ² = 0

as the basic equation, on the principle of attaching an imaginary symbol to what cannot be observed and moving it to the left-hand side.

A surprising amount came out of this form.

- From a wave of a single frequency, after a run-up period, a geometric expansion occurs and the system moves autonomously into a metastable state
- During that transition, three directions appear in the system
- The number of waves is not determined inside the system; it is set by the resolution given from outside
- Introducing an interaction with a reflection coefficient reproduces fermionic elastic reflection and the contractive reaction of a wave packet
- There is a special solution near reflectance 0.7, and the exact root of finite-order recurrence sits near α⁻¹ ≈ 137

And finally I reached **a periodic table of 62 particle-like waves, built from the waveform alone**.

---

## But the dynamics hit a wall

Here I ran into a wall.

Σ xₙ² = 0 **contains the many-body problem for three or more bodies**. So the dynamics could not be derived.

It had been so easy in the real case. With Σ xₙ² = R², "conserve R" was the whole of it. The moment I went complex, it became unclear whether that advantage was still alive.

**Can the complex extension be solved the same way as central projection?**

That is the question this paper answers.

---

## The answer: yes. And the form turned out to be a light cone

The conclusion first. Yes — and in a clean form.

Writing xₙ = qₙ + i pₙ and separating Σ xₙ² = 0 into real and imaginary parts, algebra alone gives two equations.

Σ qₙ² = Σ pₙ²
Σ qₙ pₙ = 0

The first is the essential one. **The sum of squares of the real parts equals the sum of squares of the imaginary parts.** This is exactly the same type of condition as the real central projection Σ xₙ² = R². Going complex did not destroy the "map onto a surface".

And separating what can be observed from what cannot, this becomes

x² + y² + z² = t² + R² + Q²

Transposing,

r² − t² − R² − Q² = 0

**A light cone.**

That is, the basic representation is the **four-dimensional** (r, t, R, Q), and zero closure gives a light-cone-type null cone on it. That is the meaning of the title.

And just as in the real case, **so long as you satisfy the single constraint of conserving the right-hand side, you may move freely on the surface.** There is no need to solve the many-body problem term by term. The advantage of central projection was not lost in the complex world.

---

## Finding 1: it should have been n-dimensional, but it was a three-dimensional ellipsoid

Along the way, something I had not expected became clear.

Σ xₙ² = 0 is fundamentally an n-dimensional affair. And yet the state can be represented as **an arbitrary point on a three-dimensional ellipsoid**.

The figure below is an actual computation at resolution N = 16 (a system with 120 relations). You can see the 16 vertices genuinely trying to lie on a single ellipsoid.

![Figure 1](複製_ダンプ版_v1/次元の生成構造/対照実験_N掃引1to20_三系_v2/figures_tau/fig_ellipsoid_tau09487_electron_T40000_d0.1_rep-dump40k16_N16_m_v1.png)

The second panel from the top is the clearest. The red, blue and green axes are the three radii of this ellipsoid. **An n-dimensional many-body system reads as a single point on a three-dimensional surface.**

One more, this time from the geometric side.

![Figure 2](figures_v1/fig2_parallelotope_d4_projection.png)

This is the shadow, cast down into three dimensions, of a four-dimensional parallelotope (a skewed hypercube). Sixteen vertices, 120 segments. Classify the segments by "the dimension of the smallest face containing both endpoints" into four classes, attach signs, and the alternating sum lands exactly on zero. The numbers on the right of the figure are the measured values:

32 − 96 + 96 − 32 = 0

The error is 7×10⁻¹⁵ — machine precision itself. **The 16 vertices lie exactly on a four-dimensional ellipsoid.** Error 3×10⁻¹⁶.

---

## Finding 2: the extra dimensions have not compactified

This is what I most want to say this time.

A system at resolution N = 16 has 15 principal axes. Only the top three are what we see as three dimensions. So where did the other twelve go?

In higher-dimensional theories the usual account is that **the extra dimensions curled up small and became invisible**. Compactification.

**But that is not what happens in this system.**

Look at the bottom panel of Figure 1. It shows how the values of the 15 principal axes change with time τ. Above zero is a real axis, below zero an imaginary one.

At τ ≈ 9000, something happens.

- The top three directions (A, B, C) jump upward. The A axis by a factor of 2.80
- The middle five (D through i) barely change in magnitude
- The lower ones **pass through zero, drop below, and become imaginary**

And here is the decisive point. The four axes that fell to imaginary **did not shrink and vanish. They changed sign while increasing in absolute value.**

- p axis: ×5.57
- o axis: ×2.62
- n axis: ×1.67
- m axis: ×1.12

**They are not shrinking. If anything, they are growing.**

So what is happening is not "the extra dimensions curled up and disappeared" but

**only three directions expanded geometrically**

The top-3 occupancy goes from 0.384 to 0.812. The total absolute sum grows by 2.79. At the very least, one can state flatly that this is not compactification of the "extra directions shrink until they cannot be seen" type.

It looks like an inflationary expansion.

---

## Finding 3: the conservation law holds — but a misreading is forbidden

What is interesting is that even in the midst of this violent expansion, **the signed trace is perfectly conserved**.

That is the black dashed line in panel (c) of Figure 1. Horizontal across the whole of τ. In measurement, the minimum and maximum over the last 968 points agree to eight digits.

An important caution here.

**Do not read "the conservation law holds, therefore no expansion is occurring."**

The grey dotted line in the same figure is the sum of the positive eigenvalues alone. It is not horizontal; it jumps by a factor of 1.90. The total in the real directions grows, the total in the imaginary directions appears from zero and grows, and **it is their difference that is conserved**.

Real and imaginary both grow and cancel within the signed sum. That is what is happening in this system.

And one more thing. Imaginary directions never appear on the vacuum side. In the vacuum all 15 axes are real and perfectly isotropic (all 0.0645). **Imaginary directions appear together with the creation of matter** — that is the conclusion of the control experiment.

---

## What this paper does not claim

Having made strong claims, I state the boundaries too. In this paper I do not claim the following.

- I do not claim that principal axes such as A, B, C correspond to axes of real time or space. There is at present no ground whatsoever for such an identification
- I do not claim that the number 15 of principal axes is the same thing as the 11 dimensions of superstring theory. This number is N−1, arising from the resolution; it is not a number chosen from outside
- I do not claim to have derived the law of dynamics itself. What I derived is the structure of the reduction: which single quantity must be conserved for the motion to close

I write these lines every time so that I do not lose track of how far my own claims are supported.

---

## Still, the road went through

What had blocked me for a long time was that, the moment I went complex, it was unclear whether the greatest advantage of central projection — not having to solve the many-body problem — was still alive.

This time it went through.

**Σ xₙ² = 0 is a light cone on the four-dimensional space (r, t, R, Q).**
**And if you satisfy the single constraint of conserving the right-hand side, you can design the motion without solving the many-body problem.**

From here on it is dynamics. At last I am standing at its entrance.

---

The paper itself (Japanese and English, with tex and PDF included) is published on Zenodo.

Version DOI (this version): 10.5281/zenodo.21902806
https://doi.org/10.5281/zenodo.21902806

Concept DOI (always points to the latest version): 10.5281/zenodo.21902805
https://doi.org/10.5281/zenodo.21902805

All programs that reproduce the figures and numbers are in the repository, listed with md5 hashes. The random seeds are fixed, so you get the same numbers.

#Physics #TheoreticalPhysics #Mathematics #Geometry #ComplexNumbers #CentralProjection #ZeroClosure #ExtraDimensions #Compactification #Inflation #LightCone #IndependentResearcher #ResearchNotes #Zenodo #OpenScience
