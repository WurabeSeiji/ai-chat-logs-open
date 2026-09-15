# Not Advancing Time from Outside — the Minimal Complete-Relation Closed System Born from an Inflation Experiment

Published: 2026-09-15 (Japanese note: https://note.com/kiharanoriaki/n/n6f3dfc53df36)

I have published the paper
"A Finite Complete-Relation Closed System with Time Evolution as Vertices: Exact Minimal Closure at N=1, K=4 under U^K=I and Complete-Network Analysis of Even and Odd Harmonics".

This work grew directly out of my earlier numerical experiment on inflation-like rapid expansion in self-consistent closed systems of relational waves.

In that earlier study, when a large number of complex relational waves interact, I observed that a highly symmetric initial state — concentrated near the 90-degree series — undergoes a rapid growth of its transverse component and moves to a quasi-stable state with nearly equal amplitudes and widely dispersed phases.

This time I simplified the question further and asked it again.

Does the operation of "advancing time" really have to be given from outside?

That is the starting point of this paper.

## First: what happened in the earlier inflation experiment

In the preceding experiment, a complex relational wave was placed on every pair of N vertices.

The number of relations is

　M = N(N-1)/2

With N=40 there are 780 relational waves.

In the initial state, these waves were strongly constrained near four directions about 90 degrees apart in the complex plane, forming a cross of narrow bands.

![Initial state at N=40, M=780](英語版図_inflation_en_20260915/en/inflation_N40_phase_step0_cross_en_20260915.png)
The initial state at N=40, M=780. The relational waves concentrate in four directions near the 90-degree series.

When this state interacts, the transverse component grows rapidly over a short interval, and the distribution finally moves to a ring with nearly uniform amplitudes.

![Quasi-stable state at step 500](英語版図_inflation_en_20260915/en/inflation_N40_phase_step500_ring_en_20260915.png)
The quasi-stable state at step 500. The waves spread onto a ring of nearly equal radius, while the anisotropy of the initial 90-degree series remains.

I call this rapid transition "inflation-like".

The word inflation here does not presuppose identity with cosmological inflation. It names the behavior actually observed in the numerical experiment: an order parameter changing greatly over a short interval, from a highly symmetric state to a different quasi-stable state.

I also made a Mexican-hat effective display reconstructed from the measured values.

![Mexican-hat effective display](inflation_mexican_hat_SSB_reconstruction_20260915.png)
A Mexican-hat effective display reconstructed from the measured evolution of the order parameter. It is not a figure drawn by assuming a potential in advance.

The next figure shows, in three dimensions, how the many relational waves move from the initial state to the quasi-stable state.

![3D display of the N=40 relational waves](英語版図_inflation_en_20260915/en/inflation_N40_relation_wave_3D_step500_en_20260915.png)
A 3D display of the N=40 relational waves expanding out of the highly symmetric constrained state.

## Then one thing started to bother me

In that experiment, the states were updated in the order

　state at τ → state at τ + Δτ

with the update step

　Δτ = 2π/N

As a numerical experiment, this is the ordinary way to write it.

But my systems take anonymity seriously: no state and no relation is given a name or a special role from the start.

Then a question arises.

Why is τ alone — neither a vertex nor a relation — a special variable that advances the whole system from outside?

Moreover, deciding in advance that "there are N vertices at a given instant" may itself not be a fundamental structure.

So I inverted the idea.

Instead of advancing time from outside,
place the distinct phase states themselves as ordinary vertices — all of them.
Do not distinguish anything like past, present, and future from the start.

If there are K states, each carrying N objects, treat the whole as

　P = KN

vertices, and relate every pair of them completely. The number of relations is

　M_P = KN(KN-1)/2

No assumption "k is time" enters here. The K distinct phase states are simply treated as vertices of the same kind as everything else.

## Not approximating a huge system — solving the minimal exact system completely

Carrying this idea straight to huge systems makes the number of relations explode.

So in this paper I did not take the path of coarsely approximating a large system.

The opposite.

Keep the same closure axiom and the same complete relational structure fully intact, shrink to the minimal size, and solve all states and all relations completely.

In the preceding inflation experiment, the initial state clearly showed a four-direction structure of the 90-degree series. So as the minimal cycle I adopted

　N = 1
　K = 4

Even with N=1, if all four phase states become vertices, the whole is the complete graph K4 with 4 vertices and 6 relations.

The phase generator is

　U = exp(2πi/4) = i

with

　U^4 = 1

The states compared in this work are

　S_k = (U^k, U^(mk))

I compared m=2, the minimal nontrivial even harmonic, with m=3, the minimal nontrivial odd harmonic.

## The same K4 — but the inside differs between even and odd harmonics

The outer graph structure is the same K4 with 4 vertices and 6 edges in both cases.

But the internal cyclic structures differ.

For m=2,

　ord(U^2) = 2

and the harmonic component alone degenerates to a 2-cycle.

For m=3,

　ord(U^3) = 4

and both the base wave and the harmonic keep all four states.

Even so, the full two-component state closes in 4 states for both m=2 and m=3.

This is the first important result.

![K4 complete network for m=2](K4_complete_network_m2.png)
The K4 complete network of the even harmonic m=2. The outer combinatorial structure closes completely with 4 vertices and 6 relations.

![K4 complete network for m=3](K4_complete_network_m3.png)
The K4 complete network of the odd harmonic m=3. The same K4, but the internal phase, distance, and squared-closure structures differ from m=2.

## And the squared sum showed a decisive difference

Sum the squares of the two wave components over all four states:

　Q_loop(m) = Σ[(U^k)^2 + (U^(mk))^2]
　k = 0,1,2,3

Then, at K=4, exactly:

　m even → Q_loop = 4
　m odd → Q_loop = 0

So in the minimal comparison,

m=2 does not close to squared-sum zero.
m=3 closes to squared-sum zero.

This is not a numerical approximation. It follows directly from geometric series over a finite cyclic group.

For general K,

　Q_loop(m;K) = K × 1[K divides 2] + K × 1[K divides 2m]

So the even/odd difference seen at K=4 is understood, more generally, as a divisibility condition.

I consider this one of the most important points of the paper.

## Movies make it clear what is being computed

This system is not about showing "waves flowing continuously".

A finite number of discrete phase states cycle on a fixed complete network.

So I made forward / reverse MP4 movies for m=2 and m=3.

(Movie: [m=2 forward](K4_harmonic_discrete_phase_m2_20260915.mp4))
m=2 forward. The harmonic component degenerates to 2 states, while the full two-wave state closes in 4 states.

(Movie: [m=3 forward](K4_harmonic_discrete_phase_m3_20260915.mp4))
m=3 forward. Both the base wave and the harmonic keep 4 states, and the whole closes to squared-sum zero over the 4 states.

There are also reverse versions confirming that reading the cycle backwards returns to the same finite state set.

This does not mean "physical time-reversal symmetry has been proven". It visualizes that the same fixed finite cycle can be traversed in the reverse reading order as well.

## What this paper deliberately does not claim

I intentionally separated what is mathematically established from what remains a candidate physical interpretation.

Mathematically established are:

- the complete relational structure of K4
- the finite cycle under U^4=1
- the order difference between even and odd harmonics
- the difference in state-space distances
- the squared-sum closure condition
- the divisibility condition for general K

Still left as hypotheses are:

- reading N=1 as "one photon"
- reading k as physical past, present, and future
- reading the state-space distance as real-space or spacetime distance

This separation is deliberate. Even if the physical interpretation changes in the future, the mathematical results about the finite complete-relation system remain as they are.

## What I consider most important this time

In the earlier inflation experiment, states were updated from outside, one after another.

In this paper, that external updating is removed for once, and the distinct states themselves are all placed as vertices.

In other words, the framing of the problem changes from

"evolve the states in time"

to

"can something like order and time be read out of the relations among all states?"

The K4 of this paper is the minimal exactly closed system for that purpose. It is not a model approximating a huge universe by four vertices. It is the system shrunk to the minimal size at which all states and all relations can be analyzed without dropping a single one, while keeping the same complete relations and closure condition.

And already in this minimal system, even and odd harmonics separate sharply in orders, distance geometry, and squared closure.

The next task is whether the cyclic order itself can be reconstructed from the complete relations alone, after removing the order label k.

If that is possible, we may move from descriptions that give "time" from outside to descriptions that read temporal order out of relations.

## Paper and PDF

Title:
A Finite Complete-Relation Closed System with Time Evolution as Vertices: Exact Minimal Closure at N=1, K=4 under U^K=I and Complete-Network Analysis of Even and Odd Harmonics

Author: Noriaki Kihara
ORCID: 0009-0004-6753-4020

Version DOI:
https://doi.org/10.5281/zenodo.22763330

Concept DOI:
https://doi.org/10.5281/zenodo.22763329

Zenodo record:
https://zenodo.org/records/22763330

Preceding work:
The Mechanism of Inflation-like Rapid Expansion in Self-Consistent Closed Systems of Relational Waves
Version DOI: https://doi.org/10.5281/zenodo.22176949
Concept DOI: https://doi.org/10.5281/zenodo.22112008

#TheoreticalPhysics #QuantumMechanics #Cosmology #Time #ComplexNumbers #DiscreteGeometry #GraphTheory #CompleteGraph #Harmonics #SelfConsistency #Inflation #NumericalExperiment #Physics #ResearchNote
