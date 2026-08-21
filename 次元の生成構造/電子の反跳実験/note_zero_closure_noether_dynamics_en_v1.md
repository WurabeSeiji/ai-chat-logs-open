# I Was Chasing Noether's Theorem, and "Ordinary Physics" Came Out of Discrete Waves

Noriaki Kihara
August 2026

Original paper (public version v1.0, published 21 August 2026)

- DOI (always resolves to the latest version): https://doi.org/10.5281/zenodo.22040735
- DOI of this version: https://doi.org/10.5281/zenodo.22040736
- English PDF (GitHub repository): https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E9%9B%BB%E5%AD%90%E3%81%AE%E5%8F%8D%E8%B7%B3%E5%AE%9F%E9%A8%93/zero_closure_noether_dynamics_en_public_v1.0.pdf
- Japanese version of this article: https://note.com/kiharanoriaki/n/n91202fa73800

The result surprised me a little.

What I wanted to study this time was dynamics.

Up to now I have taken two basic conditions,

Σ Xₙ² = 0

which I call "zero closure", and

Uᴺ = I

which I call "finite recurrence", and asked how much physical structure emerges from the relations among waves alone, without placing a background spacetime or any ready-made symmetry at the start.

But once the research reaches this point, a big problem remains.

How do things move in that world?

In other words: dynamics.

## The trigger was Noether's theorem

So I went back and started studying Noether's theorem.

It is a very famous theorem in physics.

Put simply:

**where there is a symmetry, there is a conservation law.**

If the laws of physics do not change when you shift time, energy is conserved.

If they do not change when you shift space, momentum is conserved.

If they do not change when you rotate, angular momentum is conserved.

This idea sits very deep in the dynamics of modern physics.

My model, however, has an awkward feature.

Time is not placed at the start.

Neither is space.

What exists, first of all, is only the closed relation

Σ Xₙ² = 0

So where on earth would Noether's theorem come from?

## There was no need to change the amplitudes of the waves

What mattered here were the equal-amplitude waves I have been using all along.

If the amplitudes are equal, the essential quantity that distinguishes the waves is the phase.

So I considered

**advancing only the phases a little, without breaking zero closure.**

At finite resolution, the condition that preserves zero closure then becomes the finite-difference equation

Σ Xₐ² (e^(2iΔφₐ) − 1) = 0

Up to this point no continuous derivative has been used.

No time derivative.

No partial differential equation.

It is purely a discrete update of phases. This is the exact finite formula obtained in the present paper.

## What if the resolution is made infinitely fine?

Here a thought occurs.

What happens if the resolution of this discrete system is made finer and finer,

N → ∞ ?

The finite phase differences then approach a continuous phase gradient.

And a local conservation structure appears:

Σ Xₐ² ∂φₐ = 0

In other words, the order was reversed.

Normally one thinks

continuous spacetime → differential equations → conservation laws.

In this model it is

**discrete closure → finite phase differences → higher resolution → continuous conservation structure.**

This made possible a different view: perhaps continuous physical laws are the limit of a discrete structure seen at very fine resolution.

## And this is where it became unexpected

Originally I was only studying the relation to Noether's theorem.

But as I followed this continuum limit, connection routes began to appear toward the standard dynamics that had looked separate until now.

Furthermore, if one constructs a current from the phase gradient, a connection to Maxwell's equations may lie on the same extension.

This is not yet at the stage where I would say "derived".

Indeed, the paper itself states explicitly that a complete identification with the standard Noether theorem

∂μ Jμ = 0

still requires showing the correspondence between the construction of the current and the variational symmetry.

But at least the path has become considerably clearer.

## What about relativity, then?

The remaining big destination is the theory of relativity.

Here, however, there has been an entrance for some time.

If zero closure is split into real and imaginary components and read as

x² + y² + z² − t² = R² + Q²

a Minkowski-type metric already appears on the left-hand side.

That is, it is not necessarily the case that

**relativity must be glued on afterwards, separately from quantum theory and field dynamics.**

Relativistic spacetime may be derivable as a different "way of reading" the same zero closure.

This is the next important task.

## This time I deliberately approached "current standard physics"

In this paper, for the sake of comparison, I deliberately looked at a very special limit.

One is the continuum limit

N → ∞

The other is the region where the curvature radius R is very large and things look almost flat locally.

This is to match the region where present standard physical theory succeeds with extremely high precision.

This point matters.

I am not saying that the standard theory is wrong.

Rather the opposite.

In regions that are nearly flat and where the continuum approximation holds very well, the standard theory is astonishingly accurate.

The possibility that has come into view is this:

**it may be possible to explain why the standard theory is that accurate, from the discrete structure one level below it.**

## But "normalisation" has to be reconsidered

And here there is also a slightly troublesome issue.

In standard quantum theory one normalises states as

⟨ψ|ψ⟩ = 1

This is so ordinary computationally that one hardly notices it.

In my model, however, R is not merely a nuisance scale; it remains as a quantity related to curvature and to the size of the system.

Dividing by R to normalise may therefore hide geometric information that was originally there.

In an earlier study, normalisation by R could be read as "a similarity transformation mapping a space with curvature radius onto the unit sphere".

I think this is a point that needs considerable care when the correspondence of formulas with the standard theory is examined in future.

## The overall picture as it looks now

When I started this research, I did not think it would come this far.

The starting point was very simple.

What happens if one places the closed relation

Σ Xₙ² = 0 ?

From there I studied waves, studied phases, studied symmetries, and then chased Noether's theorem.

And now a single road has begun to appear:

discrete zero closure
↓
discrete dynamics by phase differences
↓
N → ∞
↓
continuous conservation laws and field dynamics
↓
the standard theory

And beside it there is another road, toward relativity:

Σ Xₙ² = 0 ⟶ x² + y² + z² − t² = R² + Q²

If these two roads join at the same place, things will become quite interesting.

## In one line

If I had to put what was learned this time into one line:

**While chasing Noether's theorem, a road from the continuum limit of discrete zero closure toward standard dynamics came into view.**

And what I will examine next is

**whether this road really continues unbroken to Maxwell's equations, and on to the theory of relativity.**

The paper itself is here.

- "Noether Conservation Laws and Relational Phase Dynamics from Discrete Zero Closure — A Zero-Closure-Preserving Discrete Self-Map, the N→∞ Continuum Field Equations, Local Gauge Geometry, and the Standard-Model One-Generation Representation with Chirality Selection" (public version v1.0)
- DOI (always resolves to the latest version): https://doi.org/10.5281/zenodo.22040735
- DOI of this version: https://doi.org/10.5281/zenodo.22040736
- English PDF (GitHub repository): https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E9%9B%BB%E5%AD%90%E3%81%AE%E5%8F%8D%E8%B7%B3%E5%AE%9F%E9%A8%93/zero_closure_noether_dynamics_en_public_v1.0.pdf
- Japanese PDF (GitHub repository): https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E9%9B%BB%E5%AD%90%E3%81%AE%E5%8F%8D%E8%B7%B3%E5%AE%9F%E9%A8%93/zero_closure_noether_dynamics_ja_public_v1.0.pdf

## Related papers

- Preceding paper: "Symmetry Generation from Zero Closure, Finite Order, and Self-Consistent Geometry — The Single External Parameter N and the Remaining Tasks of Generalization and Dynamics" https://doi.org/10.5281/zenodo.22028072
- "Zero Closure Was Four-Dimensional — 'Central Projection' Survives Even in the Complex World" https://doi.org/10.5281/zenodo.21902805
- "The Periodic Table of Waves v2 — Particle Classification by Winding-Number Address and Observation Clock, and the Unification of Mass, Lifetime, and Splitting via the Clock Field ω(x)" https://doi.org/10.5281/zenodo.21830706

#TheoreticalPhysics #NoetherTheorem #ConservationLaws #Dynamics #QuantumMechanics #StandardModel #GaugeTheory #MaxwellEquations #Relativity #DiscreteSystems #ContinuumLimit #ComplexNumbers #Geometry #IndependentResearch #Preprint
