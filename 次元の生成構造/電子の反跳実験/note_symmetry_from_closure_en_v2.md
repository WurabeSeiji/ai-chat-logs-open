# Do the Symmetries of Physics Really Have to Be Given from the Start? — Reverse-Engineering Spacetime, Curvature, Quantization, and Gauge Symmetry from Just Four Closure Conditions

Noriaki Kihara
August 2026

Original paper (public version v1.0, published August 20, 2026)

- DOI (always redirects to the latest version): https://doi.org/10.5281/zenodo.22028072
- English PDF (direct download from the GitHub repository): https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E9%9B%BB%E5%AD%90%E3%81%AE%E5%8F%8D%E8%B7%B3%E5%AE%9F%E9%A8%93/closure_axioms_symmetry_derivation_en_public_v1.0.pdf

There is something a little strange about physics.

The deeper the theory goes, the more symmetries it accumulates.

Lorentz symmetry.
Gauge symmetry.
SU(3).
SU(2).
U(1).
Cyclic symmetry.
Permutation symmetry.
Symplectic structure.
Conformal structure.

Every one of them is an important structure supporting modern physics.

But I decided, for once, to think in the opposite direction.

**Do these really have to be posited from the start, each as a separate principle?**

The question investigated in this paper is a rather extreme one.

> **If we posit only a handful of closure conditions, how much of the symmetry that modern theoretical physics demands appears afterwards, automatically?**

The result was stronger than I myself had expected.

The negative signature of spacetime, the curvature radius, finite quantization, Born-type squared weights, simplex geometry, fixed-point symmetries — and, conditionally, even the global gauge group corresponding to the Standard Model's

SU(3)×SU(2)×U(1)

— a remarkably long stretch of structure could be lined up as a single derivation chain from the same axiom system.

Moreover, within the derived range, this system has no continuously adjustable free parameters.

The only essential quantity specified from outside is a discrete integer N.

This article introduces the overall picture for a general audience.

## The Starting Point Is a Single Strange Equation

The first axiom is

Σxₙ² ＝ 0 (n ＝ 1 … M)

with the xₙ being complex numbers.

With real numbers only, the sum of squares can be zero only if everything is zero.

With complex numbers, things are different.

In the smallest example,

x² + y² ＝ 0

gives

y ＝ ±ix

So the imaginary direction is not a tool inserted afterwards for computational convenience.

**It appears as a direction required for a nontrivial square closure to hold.**

This is the starting point of the paper.

## Could "Imaginary" Be the Sign of an Invisible Direction?

Consider, for example, the four components

(x, y, z, iR)

The first axiom is

x² + y² + z² + (iR)² ＝ 0

which becomes

x² + y² + z² ＝ R²

A familiar sphere.

What matters is that the R² on the right-hand side was not added afterwards.

The complex direction iR was there from the start; the observer simply does not read it directly.

Then the coordinate value of that invisible direction appears, from the three visible directions, as a **curvature radius**.

Furthermore, with

(x, y, z, it)

we get

x² + y² + z² + (it)² ＝ 0

and therefore

x² + y² + z² − t² ＝ 0

The Lorentz-type negative sign did not come from imposing a (+, +, +, −) metric from outside.

It is simply

(it)² ＝ −t²

In this view, the negative sign of time and the negative sign of the curvature radius have the same origin.

## Time, Curvature, and Internal Quantity Are Not Different Things from the Start

Going further, with

(x, y, z, it, iR, iQ)

we get

x² + y² + z² ＝ t² + R² + Q²

Here I do not posit t, R, and Q from the start as separate entities called time, curvature, and internal quantity.

In the fundamental structure, they are all unobservable complex directions of the same kind.

Whatever is read as a clock by some observation map is called t.

Whatever is read as a curvature scale is called R.

Whatever is read as an internal quantity is called Q.

In other words,

> **The name does not come first; the name is attached according to how the thing is read.**

That is the order.

This is quite important.

## So Is This a Theory Where You Can Read Anything Any Way You Like?

This is the most easily misunderstood point.

The answer is the opposite.

**This model is quite cramped.**

The fundamental components carry no names, but that does not mean they can be interpreted freely.

An admissible structure must simultaneously satisfy

- the zero closure of complex squares
- a phase closure that returns exactly after finitely many steps
- distance consistency as a simplex
- harmonic structure
- self-consistent fixed points

Conceptually,

S(allowed) ＝ S(zero closure) ∩ S(finite recurrence) ∩ S(simplex) ∩ S(harmonics) ∩ S(self-consistency)

In other words, however many candidates there are, only those passing all the conditions at once survive.

The direction of this theory is not anything-goes.

Rather, it is

**extremely few inputs ＋ extremely strong constraints ⟹ extremely rich structure**

I consider this one of the most important results this time.

## The Second Axiom: "Return Exactly After Finitely Many Steps"

The second axiom is

Uᴺ ＝ I

Apply the action N times, and you are back exactly where you started.

The eigenvalues then satisfy λᴺ ＝ 1, so

λₘ ＝ e^(2πim/N)

That is, we get at once:

- finite recurrence
- discretization of phase
- the cyclic group Zₙ
- cyclotomic structure

The Born-type squared weights found in earlier work,

cos²(πm/N), sin²(πm/N)

also appear when this finite recurrence is read through two channels.

Probability was not put in at the start.

First there is a finitely recurring phase structure, and the Born-type weights emerge as its squared readout.

This is the part treated in detail in a previous article.

## The Third Axiom Turns Relations into "Space"

Next, read the same N as a number of vertices.

The total number of pairwise relations is then

M ＝ N(N−1)/2

We require that all these pairwise distances be consistent as a simplex.

Then geometry stands up from distance relations alone, without positing absolute coordinates from the start.

Moreover, a simplex has a boundary operator with

∂² ＝ 0

and, dualizing,

d² ＝ 0

So from a single simplex condition, the entrances to

- distance geometry
- homology
- cohomology
- discrete differential forms

all appear at the same time.

By this point, it becomes rather doubtful whether the idea of positing space from the start is needed at all.

## The Fourth Axiom: "The World You Generate Coincides with Yourself"

The fourth axiom is

X ＝ F(X)

The self-consistency condition.

Feed the generated state back into the generating rule, and it returns to itself.

For such a fixed point X*, collect exactly the transformations that leave it unchanged:

G(X*) ＝ { g ∈ G ｜ g·X* ＝ X* }

This is the stabilizer.

So even without declaring from the start which group the symmetry is, one can read symmetry as

> **the transformations that do not change a structure that managed to exist self-consistently.**

Here the relation to the second axiom Uᴺ ＝ I also becomes visible.

One is the finite self-recurrence of an action.

The other is the self-recurrence of a state.

Both share, at a higher level, the same structure: **self-closure**.

## And Three Dimensions Turned Out Not to Be an Accident

In the earlier numerical experiment make_parent, we confirmed that when a self-consistent coherent wave is evolved, three principal axes come to dominate selectively out of many candidate directions.

But organizing things this time, three-dimensionality also has a geometric derivation route separate from the numerics.

Splitting the first axiom Σxₙ² ＝ 0 into real and imaginary parts gives

aᵀa ＝ bᵀb, aᵀb ＝ 0

that is,

‖a‖ ＝ ‖b‖, a ⊥ b

an equal-norm orthogonal structure.

There is an orthogonal plane, and there is its normal direction.

Layering phase closure and simplex geometry on top leads to an ellipsoidal structure with three distinguishable principal axes.

On the numerical side: spectral concentration onto three principal axes.

On the geometric side: a three-dimensional ellipsoidal structure.

The same "3" came out along separate roads.

What remains is not whether three directions can be produced.

**It is to turn into a general theorem why this three-direction selection is universal across wide ranges of N, initial conditions, and seeds.**

## From Here, Even the Standard Model's Gauge Group Came into View

This was one of the most surprising parts.

Reading the zero closure through six complex registers,

ΣXₙ² ＝ 0 (n ＝ 1 … 6)

That is six complex degrees of freedom with one complex constraint, so

complex dimension ＝ 6 − 1 ＝ 5

This is simple dimension counting.

And when these five complex degrees of freedom enter a sector preserving a self-consistent Hermitian decomposition

V ＝ V₃ ⊕ V₂

the conservation group is

U(3)×U(2)

Removing the redundancy of the global phase gives

S(U(3)×U(2))

and, as known group theory,

S(U(3)×U(2)) ≅ ［SU(3)×SU(2)×U(1)］/Z₆

This is the same form as the faithful global gauge group of the Standard Model.

A careful distinction is needed here.

**I am not claiming to have completely derived the Standard Model from the first axiom alone.**

What comes out rigorously is the five complex degrees of freedom.

If from there the 3⊕2 Hermitian decomposition is selected self-consistently, the same global gauge group as the Standard Model appears mathematically.

The remaining problem is:

> **Why would self-consistent dynamics select that stabilizer?**

But put the other way around: without positing the Standard Model's three gauge groups separately from the start, we have reached the point where they **can appear simultaneously as the conservation group of a single five-complex-degree-of-freedom structure**.

## And There Was More Than One Road to SU(3)

The internal quantity can also be read more finely, as

Q² ＝ Q₁² + Q₂² + Q₃²

Then an internal triplet (Q₁, Q₂, Q₃) appears.

Give it a Hermitian structure, and a second route to an internal U(3), and then SU(3), becomes visible.

Also, reading (R, Q₁, Q₂, Q₃) as four components connects to the known group theory

SU(4) ⊃ SU(3)×U(1)

Meanwhile, on the four-dimensional side there is

Spin(4) ≅ SU(2)×SU(2)

So a route toward the Pati–Salam-type structure

SU(4)×SU(2)×SU(2)

also comes into view as another readout of the same zero closure.

I do not classify this as "Pati–Salam has been derived".

The important point is different.

**Read differently, multiple routes toward the known grand-unified group structures appear from the same zero closure.**

That is what makes this interesting.

## "Multiple Readings" Is Not a Weakness

Normally one would think:

> If you can change the reading, can't you get anything you want?

But in this model it is the other way around.

The fundamental components have no names.

Personally, I call this anonymity.

But anonymity does not mean freedom.

It must close to zero.

It must recur after finitely many steps.

It must be consistent as a simplex.

It must be consistent with the harmonic structure.

It must return to itself self-consistently.

Only the sectors that pass all of that remain.

So even symmetry may appear not as something chosen from outside, but as

**the automorphisms of the structures that survived the strong closure conditions.**

This inversion of order is the real subject of the paper.

## Why Were So Many Symmetries Needed?

Looking at modern physics, nature has an astonishing number of symmetries.

Usually one assumes each has a deep physical meaning.

I think that is right, of course.

But this model suggests another way to see it.

**Perhaps there are many symmetries not because there is much freedom, but because there is extremely little.**

It cannot exist unless it closes.

It must return to where it started.

The distances must not contradict each other.

Generating itself must yield itself.

Try to satisfy such extremely strong conditions all at once, and whatever structure remains necessarily has high symmetry.

In other words,

**high symmetry ＝ the result not of high freedom, but of high constraint**

That may be what is going on.

I believe this possibility became quite clearly visible this time.

## N Is Not the Number of Elements in the Universe — It Is Closer to the Resolution of a Microscope

This is another easily misunderstood point.

The N appearing in this model does not mean the universe is really made of N elements.

A closer analogy is **the resolution of a microscope**.

The same object looks coarse at low magnification. Raise the magnification, and what looked like one structure resolves into finer parts.

It is the same in this model:

Uᴺ ＝ I

The N of this finite recurrence is treated as a finite readout order expressing at what fineness the closure structure is read.

**N ＝ not the number of existing elements, but the resolution of the readout**

Therefore, increasing N reads the same closure structure ever more finely.

And

M ＝ N(N−1)/2

does not mean the world essentially contains M relations.

It is the number of pairwise relations that become readable when N distinctions are made at that resolution.

Raising a microscope's magnification does not add new components to the specimen.

Likewise, refining N in this model does not add new physical degrees of freedom; it **reads the same closure structure at higher resolution**.

Consequently, this model has no problem of "finding the correct N".

Which N to read at is a matter of observation and selection.

What matters is not the value of N itself, but that the same closure principle operates behind every finite resolution.

Note also that this model does not posit

- an R for fitting curvature
- a t for fitting time
- a Q for fitting internal symmetry
- angles for fitting Born-type squared weights
- continuous coefficients for fitting gauge groups

as independent tuning parameters.

R, t, Q — and N — are treated not as tuning knobs pasted onto fundamental existence from outside, but as quantities appearing on the observation/readout side.

## The Derivation Map — What Comes Out, and How Far

The original paper contains a derivation map: a table grading every structure line by line. Since tables cannot be used here, it is reorganized by status.

**Rigorous (holding by the mathematics of the axioms alone)**

- Isometric orthogonal two-planes (‖a‖ ＝ ‖b‖, a ⊥ b) and complex phase rotation —— A1
- Conservation group O(M, ℂ), indefinite signature O(p, q) in real display, null cone —— A1
- Curvature-radius readout x² + y² + z² ＝ R² and curvature K ＝ 1/R² —— A1 including iR
- The negative sign of time (it)² ＝ −t², the (3, 3) signature of the six-component readout, and O(3, 3) —— A1
- Degree-of-freedom count of six complex axes: complex dimension 6 − 1 ＝ 5 —— A1
- The entry to symplectic structure —— real display of the complex structure
- Cyclic group Zₙ, cyclotomic eigenvalues, Galois symmetry —— A2
- Born-type squared weights —— A2 ＋ two-channel projection
- Complete graph Kₙ, permutation group Sₙ, distance geometry, ∂² ＝ 0, cohomology —— A3
- Stabilizers of fixed points —— A4

**Conditionally rigorous (holding once a readout sector is selected)**

- The Lorentz group O(3,1) and Spin(3,1) —— (x, y, z, t) partial readout
- U(3)×U(2) —— preservation of the Hermitian 3⊕2 decomposition of the 5 complex degrees of freedom
- S(U(3)×U(2)) ≅ ［SU(3)×SU(2)×U(1)］/Z₆ —— plus removal of the global phase (corresponding to the Standard Model's gauge Lie algebra and global gauge group)
- Symmetry reduction G → H —— A4

**Constructed and numerically confirmed in the self-papers**

- The 3-dimensional ellipsoidal structure and the selective dominance of three principal axes —— the self-consistent spectrum of make_parent
- The nesting of the coarse readout x² + y² + z² ＝ t² + R² + Q² and the refined readout Q² ＝ Q₁² + Q₂² + Q₃²
- The Bose/Fermi/mixed (Ermion) classification —— odd/even harmonic ratio and Z₂ parity

**Connections to known mathematics (independent cross-check routes)**

- The complex quadric hypersurface Q⁴ ⊂ ℂP⁵ and its real forms SO(3,3), SO(4,2), Spin(6) ≅ SU(4)
- Internal triplet → internal SU(3); (R, Q₁, Q₂, Q₃) → SU(4) ⊃ SU(3)×U(1); Pati–Salam type SU(4)×SU(2)×SU(2)
- Conformal structure —— null cone and scale

**Not yet derived (what truly remains)**

- General Riemannian dynamics —— R(q) ＋ connection
- Chirality, hypercharge, anomaly cancellation
- The dynamical selection rule among readout sectors and stabilizers, and the selection rule for N itself

## What Still Remains

While strong results came out, the remaining problems also became much clearer.

The main ones are:

1. Turning the universality of the three-principal-axis selection into a general theorem
2. Determining which stabilizer, among the multiple readout sectors, is selected dynamically
3. Deriving general Riemann curvature and geodesic dynamics in closed form from R(q)
4. Deriving local gauge connections from simplex cochains
5. Connecting the Bose/Fermi/mixed classification obtained from odd/even harmonics to the general spin-statistics structure
6. Reproducing chirality, hypercharge, and anomaly cancellation within the same derivation chain

So this paper is not a story of "the Standard Model and gravity are now complete".

But what used to look like separate grand problems has been compressed into a small number of **concrete, not-yet-derived spots**.

I consider that a major step forward.

## Reading a Century of Accumulated Symmetries Backwards

Writing this paper, this was the part I enjoyed most.

Modern theoretical physics discovered its necessary symmetries one at a time, guided by experiment.

Lorentz.

U(1).

SU(2).

SU(3).

Gauge.

Spin.

Statistics.

Curvature.

Historically, that is the natural way to proceed.

But looking at the completed structure from behind, a different question arises:

> **Why does nature use only such well-ordered mathematics, again and again?**

What I attempted this time is the inverse problem.

Instead of inputting many symmetries, posit only

**zero closure ＋ finite recurrence ＋ simplex consistency ＋ self-consistency**

Then the structures familiar from modern theoretical physics appear from it, one after another.

If this direction is right, it may not be that nature is choosing the symmetries.

**It may be that almost no structures capable of closing into existence remain in the first place.**

And that is why what remains looks beautiful.

This paper examines that possibility quite seriously.

## About the Paper

The paper separates, far more rigorously than this article:

- what is derived rigorously within the paper
- connections to known mathematics
- what has been constructed and numerically confirmed in earlier self-papers
- conditional derivations
- what is not yet generalized
- the dynamics not yet derived

In particular, for SU(3)×SU(2)×U(1), it does not settle for the single phrase "it came out"; it separates how far the rigorous group theory extends and where the problem of self-consistent physical selection begins.

This general-audience article is written rather boldly; the paper itself cuts everything into falsifiable form.

The shortest summary of this work is this:

**It is not that many structures came out of few axioms.**

**It is that imposing almost excessively strong closure conditions left few admissible structures — and the way they remained came to resemble the symmetries of modern physics.**

If this reading can be generalized further, part of what physics has called "fundamental laws" may be, rather than independent laws,

**geometric survival conditions for existing self-consistently.**

The paper itself is here.

- "Symmetry Generation from Zero Closure, Finite Order, and Self-Consistent Geometry — The Single External Parameter N and the Remaining Tasks of Generalization and Dynamics" (public version v1.0)
- DOI (always redirects to the latest version): https://doi.org/10.5281/zenodo.22028072
- DOI of this version: https://doi.org/10.5281/zenodo.22028073
- English PDF (direct download from the GitHub repository): https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E9%9B%BB%E5%AD%90%E3%81%AE%E5%8F%8D%E8%B7%B3%E5%AE%9F%E9%A8%93/closure_axioms_symmetry_derivation_en_public_v1.0.pdf

## Related Papers

- "Zero Closure Was Four-Dimensional — 'Central Projection' Survives Even in the Complex World" https://doi.org/10.5281/zenodo.21902805
- "The Periodic Table of Waves v2 — Particle Classification by Winding-Number Address and Observation Clock, and the Unification of Mass, Lifetime, and Splitting via the Clock Field ω(x)" https://doi.org/10.5281/zenodo.21830706
- "Two-Layer Separation of Waves and Fields — Unification of Gauge Fields and Gravitational Fields via a Universal Field-Readout Function" https://doi.org/10.5281/zenodo.21832257

#TheoreticalPhysics #QuantumMechanics #GeneralRelativity #StandardModel #GaugeTheory #Symmetry #ComplexNumbers #Geometry #BornRule #LorentzSymmetry #SU3 #SU2 #U1 #IndependentResearch #Preprint
