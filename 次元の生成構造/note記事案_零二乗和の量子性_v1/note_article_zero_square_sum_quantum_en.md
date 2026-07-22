# The Equation Σxₙ² = 0 Was Born Quantized

Throughout this series, I have used the equation Σxₙ² = 0 as the starting point.

Square the components, add them up, and get zero. At first sight it is a strange equation — it looks as if the only answer should be all zeros. Yet placing this equation as a closure condition, one can derive, from phase relations alone, one theoretical-physics-like phenomenon after another: position-like quantities, acceleration-like quantities, three-dimensional directions, the spontaneous splitting of waves. I have published each of these derivations as a paper.

Recently, however, I realized the order was backwards.

This equation is not a tool for deriving quantum-like phenomena. The equation itself had quantum structure built in from the start.

This article introduces, with as few formulas as possible, the content of an expository note that confirmed this by citing known mathematics alone.

## Start with the strangeness of the equation

In the world of real numbers, a square is never negative. So if you try to satisfy Σxₙ² = 0 with real numbers only, the sole answer is that every component is zero.

If you want a nonzero answer, you inevitably need numbers whose squares are negative — imaginary numbers. For example, 1² + i² = 1 − 1 = 0.

Complex numbers, in other words, are not a tool brought in afterwards; they are what the equation demands from the start. To the old question of why quantum mechanics uses complex numbers, this equation answers: in order to be nonzero and still close to zero.

## Only one more requirement: size has no meaning

The axiom system of this series has one more pillar: only ratios are physical. Multiplying an entire state by any factor gives the same state. I call this scale anonymity.

x, 2x, and 100x are all the same state. Only the ratios between components carry meaning.

Combining these two — Σxₙ² = 0, and only ratios are physical — mathematically fixes the shape of the collection of states.

## The state space turned out to be a closed surface

The set of solutions of Σxₙ² = 0 is what mathematicians call an isotropic cone (null cone): a cone radiating out of the origin.

Now impose "only ratios are physical." Each full line through the origin inside the cone collapses into a single state. What remains is the collection of directions only — mathematicians call this projectivization.

The result is a well-known object: a quadric hypersurface inside complex projective space, the projective quadric. The name is intimidating, but only one property matters:

It is a closed surface of finite extent (a compact manifold).

Not an infinite plane, but a surface closed on itself like a sphere. That is what the state space really is.

## On a closed surface, waves come only in discrete steps

Here a well-established fact connecting geometry and quantum theory takes over.

Think of a drumhead. Because the membrane is clamped at the rim — finite and closed — the sounds it makes are restricted to discrete frequencies: the fundamental, the second harmonic, the third. Nothing in between.

Waves on a closed surface behave the same way. For no other reason than that the surface is finite and closed, the wave patterns that can live on it are discrete, and the number of states at each level is finite.

The theory that organizes this in general is called geometric quantization — standard mathematics, established around 1970. It teaches one thing:

A system whose state space is compact has, with no room for choice, discrete spectra and finite-dimensional state spaces.

In ordinary quantum mechanics, discreteness is obtained by introducing Planck's constant from outside. In this system, the fact that the state space is a closed surface supplies the discreteness by itself.

So this system is not merely quantizable. It is born quantized.

## Opening up the inside of the equation

Split each component into real and imaginary parts, x = q + ip, and expand Σxₙ² = 0. The single equation separates into two real conditions.

Real part: Σqₙ² = Σpₙ²

Imaginary part: Σqₙpₙ = 0

The real part, if you read q as position and p as momentum, says that the position-side energy and the momentum-side energy balance — exactly the equipartition of a harmonic oscillator.

The imaginary part is even more interesting. The quantity Σqₙpₙ is known in classical mechanics as the generator of scale transformations — the quantity that produces uniform stretching and shrinking. Its vanishing means the equation kills the size degree of freedom by itself, from the start.

The requirement that only ratios are physical may not be an independent assumption at all: it may already be written into the imaginary part of Σxₙ² = 0. Opening up the equation reveals a structure that reads that way.

## Harmonics were not a metaphor

This series has used the word "harmonics" (overtones) again and again — the way the sound of a single string decomposes into a fundamental and its overtones, a closed state can be decomposed as finely as you like.

The functions available on the cone coincide, with mathematical exactness, with the family called harmonic polynomials. And harmonic polynomials, restricted to a sphere, are the spherical harmonics — the functions familiar from the shapes of atomic orbitals, classified by integer labels.

Spherical harmonics. The word "harmonics" literally means overtones.

The natural function system of this state space consists of overtones in the literal sense. Not wordplay — a theorem about the state space.

## Where do the integers come from?

Integer quantities keep appearing in this series — winding numbers, the number of times a phase goes around. Why are they exactly integers?

Imagine winding a string around a doughnut. However you slide the string around continuously, the number of times it wraps around cannot change. You cannot turn 3 turns into 2.9 turns; you would have to cut the string.

On a closed surface, mathematics guarantees that winding numbers of this kind are exact integers (integer invariants called Chern classes). Not approximately integer — integer with zero error, unchanged under any continuous deformation and any rescaling.

That discrete integers appear is, on this geometry, not a mystery but a property of closed surfaces themselves.

## The smallest system is a qubit

The case of three components (N = 3) turns out to be especially beautiful.

Here the closed surface of states becomes the Riemann sphere. And every point of this sphere can be written using a pair of complex numbers (a, b) — what mathematicians call a spinor — through the explicit formulas

x₁ = a² − b²
x₂ = i(a² + b²)
x₃ = −2ab

which satisfy x₁² + x₂² + x₃² = 0 identically (a classical construction given by Cartan in 1938). The ratio a : b determines the state.

A sphere whose points are ratios of two complex components: this is mathematically identical to the state space of a qubit — of spin 1/2 — familiar from quantum information.

The minimal unit of this series, the relation between two waves A and B, is isomorphic at the state-space level to a qubit. The structures called spin and SU(2) sit in a position to emerge from this geometry without additional assumptions.

## An honest boundary

What the note claims goes only as far as identifying the mathematical identity of the state space.

It contains no new theorems. All the mathematics used is known, and the references are just eight: the founding papers of geometric quantization (Kostant, Souriau, Berezin), standard textbooks of algebraic geometry, Cartan's spinors. Self-citations: zero — the note is written to be readable independently, with no prior context.

Also, the identification concerns the shape of the collection of states (kinematics); it claims no correspondence with the dynamics of the series' numerical experiments, nor with any real physical system. Even for spin, the claim stops at the indication of a mathematical isomorphism — "in a position to emerge."

Still: the moment you place the two requirements — Σxₙ² = 0, and only ratios are physical — the necessity of complex numbers, discreteness, finite dimensionality, harmonics, integers, and the qubit all follow, with no room for choice. This equation was standing, from the beginning, on the quantum geometry that theoretical physicists know well.

## The paper

Expository note, "The Geometric Identity of the Zero-Square-Sum Constraint under Scale Invariance: Isotropic Cone, Projective Quadric, and Intrinsic Quantum Structure"

Concept DOI (always the latest version):

https://doi.org/10.5281/zenodo.21495305

The paper in Japanese and English, with TeX/PDF, is public.

Related article:

"Why Do Directions Stop at Three, While Waves Multiply on Their Own?"

https://note.com/kiharanoriaki/n/nb7f682a299a6

#physics #mathematicalphysics #quantummechanics #quantization #complexnumbers #geometry #projectivegeometry #spinors #qubit #harmonics #independentresearch #preprint #Zenodo
