# Where Did Einstein's "Spooky Action at a Distance" Go? — Read Observation as a Handshake of Emission and Absorption, and the EPR Paradox Stops Being a Paradox

I have been thinking about the EPR paradox, the one Einstein called "spooky action at a distance."

Two entangled particles, even if they are as far apart as the two ends of the universe, behave like this: when one side is observed and its state is fixed (say, black or white comes out black), the state of the other side is fixed at that very moment (white, in this case). A destructive readout fixes even the state of the entangled partner far away. Against this interpretation of quantum mechanics, Einstein, Podolsky and Rosen proposed a thought experiment in 1935: an action at a distance faster than light is absurd, so the quantum-mechanical description must be incomplete. That is the EPR paradox.

This time I looked at the problem from the side of observation.

My theory does not place a background spacetime to begin with. So I had to start from "what is far?" and "what is distance?" From there the premises narrowed down to three, and from the three premises alone, the form of the Born rule, Bell correlations, the Tsirelson bound and no-signaling came out in a few lines. And when those premises are applied to EPR, the paradox disappears.

The work is published on Zenodo as two papers: the main paper "Observation as Simultaneous Emission and Absorption" (Concept DOI: https://doi.org/10.5281/zenodo.23857884) and its supplement "The EPR Paradox Is Not a Paradox" (Concept DOI: https://doi.org/10.5281/zenodo.23860605). This article follows the order of the motivation.

## Observation is an interaction

In quantum mechanics one often says that observation collapses the wave packet.

I checked that phrase earlier with two numerical experiments: whether a spread-out wave acquires localization when it interacts with a localized wave, and how a gray-cat state, with white and black nearly half and half, looks when read weakly and how it changes when observed strongly.

What came out was that observation is an interaction. Depending on the strength of the interaction, it splits into a weak readout that only shifts the phase slightly without destroying the partner's state, and a selective interaction that drops the gray cat into either the white cat or the black cat.

Put in order, the single word "observation" mixes three interactions of different character.

The first switches the partner's discrete state. As in resonant absorption, the energy level jumps and a selection happens. The strong observation that drops the gray cat to white or black is this, and so is a photodetector absorbing a photon.

The second does not switch the partner's discrete state; it only shifts the internal phases of both sides slightly. The shifted phase difference is read out as interference. Detection of a phase shift by off-resonant dispersive coupling, non-demolition photon counting in cavity QED, and dispersive readout of superconducting qubits are of this kind. The weak readout that reads the gray cat without breaking it is this.

The third is an interaction so strong, or so numerous, that the mutual phase relation itself is severed. The information scatters into the environment and can no longer be read as interference. Decoherence.

The textbook "observation collapses" identifies the third with the first and hardly treats the second. In this paper the substance of observation is placed in the second: reading the phase difference while keeping the relation. The first and the third are the boundaries at its two ends.

## What is distance?

Since no background spacetime is placed, distance cannot be brought in first. Only relations are readable.

Between two parties A and B there is a shared viewing angle. The distance from A to B is the size of B divided by the viewing angle. The narrower the viewing angle, the farther; at the same viewing angle, the larger the partner, the farther.

Drawn as a figure: take the internal common tangents of two discs as the "lines of sight," and call their crossing the relation point. The angle at the crossing is the viewing angle, the same seen from A and from B. The distance from the relation point to each disc is proportional to the size of that disc.

[Figure 1 here: fig1_viewing_angle_and_size_en.png]

Three things follow at once from this picture.

Toward oneself the seeing line and the seen line are parallel, and the viewing angle is zero. There is no crossing, so no distance is defined. One cannot observe oneself.

With two bodies alone, the viewing angle is readable, but there is no external ruler for size. Only the angle is readable; the scale of the distance is not fixed.

With three bodies a triangle forms: the three viewing angles fix the shape, and the size ratios fix the asymmetry of the coupling. Even so, scaling the whole by any factor leaves the three angles unchanged, so the overall scale is never readable.

[Figure 2 here: fig2_three_bodies_gauge_en.png]

This is the first premise. Only viewing angles and size ratios are readable; the absolute scale is not.

## What is "far"?

The second premise concerns the nature of an interaction.

When A reads the light or the gravitational wave emitted by B, B simultaneously reads its own emission through A. Emission and absorption are the two ends of one relation. They are not two events ordered in time; they hold as a single relation of zero spacetime interval. A handshake.

Then what is the speed of light? In this picture it is not the speed at which the interaction travels. It is a conversion factor: the offset between the records that the local clocks of A and B keep of this relation, divided by the distance. There is no conflict with relativity: the spacetime interval between two points joined by light is zero, and the delay appears only in each observer's coordinate time.

This way of thinking goes back to the Wheeler–Feynman absorber theory: radiation holds as a relation with an absorber, and a lone charge does not act on itself. Cramer's transactional interpretation drew this as a handshake between an offer wave and a confirmation wave. The difference this time is that the content of the handshake is placed first as a readout rule.

The third premise is that the phase of a relation carries a harmonic order that depends on the kind of wave. Turn the detector setting by θ and the phase of the relation turns by nθ/2. Spin-½ has n = 1, photon polarization n = 2, gravitational-wave polarization n = 4.

## What comes out of the three premises

Since only ratios are readable, dividing the partner's size out of each relation leaves only the direction. The overlap of two relations is the inner product cos φ. The one dimension dropped here is the logarithm of the size ratio, the same thing I wrote about in the previous article: what the Born-rule normalization throws away is the rapidity.

The form of the Born rule. A single readout holds only when the overlap seen from the emitting side and the overlap seen from the absorbing side are both in place. The two are complex conjugates, so the readable quantity is the product of an amplitude and its conjugate. The form "probability is the square of the amplitude" comes out of the handshake.

The correlation. With a source and two detectors, the overall rotation is unreadable, so the result depends only on the relative angle, and E = −cos nθ. For spin-½ this is −cos θ, for photons ±cos 2θ, both the known formulas.

The Tsirelson bound. The readout is bilinear in the directions of the two ends and the lengths are normalized to 1. So the correlation is necessarily an inner product of unit vectors, and the CHSH quantity never exceeds 2√2. The reason it stops at 2√2, neither at the local-hidden-variable 2 nor at the PR-box 4, reduces to these two facts.

[Figure 3 here: fig4_tsirelson_bound_en.png]

No-signaling. The absolute direction is unreadable from one side alone, so the one-sided distribution is uniform. The correlation appears only when the records of both sides are compared.

As a consequence, the bound is 2√2 regardless of the kind of wave, and the detector-angle step giving the maximal violation is π/(4n): 45° for spin-½, 22.5° for photons, 11.25° for gravitational waves. All three agree with standard theory; they are not discriminating points but a check that things fit.

For a localized wave the correlation becomes the autocorrelation of the pulse shape. A single tone gives cos θ, the quantum end; a rectangular pulse gives the probability of overlapping arcs, a local hidden-variable model, the classical end. Localization and non-local correlation trade off against each other.

[Figure 4 here: fig3_harmonics_autocorrelation_en.png]

## Back to EPR

Now apply all this to EPR.

First, I will not talk in terms of the direction of a spin, as is commonly done. The direction of a single spin is not readable in the first place. Look at one side alone and, whatever the axis, it is fifty–fifty; there is zero information about direction. What is read is only the sign against the detector's axis, and no reason is needed for "up here" and "up there" to be the same direction. The correlation is the inner product of the two axes; only the relative angle enters, and its reference is carried by the pair itself.

The original EPR paper is built not on spin but on position and momentum. Measure the position of particle 1 and the position of particle 2 is predicted with certainty; measure its momentum and the momentum of particle 2 is predicted, without touching particle 2. Therefore both the position and the momentum of particle 2 are elements of reality, and quantum mechanics, which cannot describe both at once, is incomplete. Those are the three steps of EPR.

In this framework a pair from one source satisfies q_A + q_B = 0 at the source. This is not a condition; it is the identity expressing that the source tied the two ends as one relation, and it holds ever after. From this one line, "the positions sum to zero" and "the momenta sum to zero" come out at the same time. The two relations that EPR found astonishing to be fixed simultaneously had held, before any measurement, as a single closure.

[Figure 5 here: epr_fig1_closure_pair_en.png]

A single end has no value. The overall phase is an unreadable gauge, so the one-sided distribution is uniform. Neither "the position of B" nor "the momentum of B" exists as a value belonging to B alone. What exists is the relation: the positions sum to zero and the momenta sum to zero.

Reading A destructively does nothing at B. An operation on A acts on A alone; B's state does not change and B's distribution does not change. The correlation is the same whether A or B is read first, and however far apart the two are. The order of two spacelike-separated readouts flips with the observer's frame. So the sentence "the readout at A instantaneously fixed the value at B" has no frame-independent meaning to begin with.

[Figure 6 here: epr_fig2_order_independence_en.png]

What experiments can show is a lower bound: if there were a finite-speed influence, its speed would exceed ten thousand times the speed of light. Moreover, assuming a finite-speed influence has been shown to allow signaling. No finite-speed action at a distance exists. That is as far as experiment goes, and this reading is consistent with it.

That a destructive readout exists, drops the chosen variable to one value and severs the relation with the conjugate variable, is placed as a known fact: resonant absorption, photodetection, and for position and momentum, homodyne detection. Deriving its mechanism from the dynamics of the closure is a separate task, and the paper does not depend on it.

The conclusion itself is not new. That a single end has no value and only correlations are real is in Peres and in Mermin; the relational-quantum-mechanics treatment of EPR, the consistent-histories position, the toy model in which EPR-type perfect correlations arise locally from an epistemic restriction alone, the formulations that fix both ends at once, and the transactional interpretation all point the same way. What is particular to this work is that position is placed as a relational quantity without a background spacetime and the EPR pair is written as an algebraic closure identity, that the absence of a one-sided value is derived from the gauge, and that the Bell side, 2√2, comes out of the same readout rule.

## Summary

Of EPR's three steps, the perfect correlation is correct. But it is not a consequence of measurement; it is the identity of a closure that held at the source. Locality is also correct: reading A does nothing at B. What was wrong was the third step, the attribution: taking a relational quantity to be an attribute of one end.

As a result, the paradox Einstein pointed to does not exist. The quantum-mechanical description is correct, and complete, as a description of relations. Neither Einstein's "incomplete" nor the textbook image of "the other side is fixed instantaneously," the action at a distance, turned out to be needed.

Observation is a handshake of emission and absorption; distance is the ratio of the viewing angle to the partner's size; "far" means a narrow viewing angle. From these three alone, the form of the Born rule, Bell correlations, 2√2 and no-signaling come out, and EPR stops being a paradox.

## Paper information

### Observation as Simultaneous Emission and Absorption (main paper)

"Observation as Simultaneous Emission and Absorption: Deriving the Born Rule, Bell Correlations and the Tsirelson Bound from a Viewing-Angle Readout Rule"

Version: v1.4 (2026-11-01)

Concept DOI:
https://doi.org/10.5281/zenodo.23857884

Version DOI:
https://doi.org/10.5281/zenodo.23857885

English PDF (Zenodo):
https://zenodo.org/records/23857885/files/paper_observation_handshake_bell_tsirelson_en_v1.4.pdf

Paper folder (Japanese and English manuscripts, figures and the program that draws them, draft history and correction list, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E8%A6%B3%E6%B8%AC%E3%81%AF%E6%94%BE%E5%B0%84%E3%81%A8%E5%90%B8%E5%8F%8E%E3%81%AE%E5%90%8C%E6%99%82%E6%88%90%E7%AB%8B_%E3%83%99%E3%83%AB%E7%9B%B8%E9%96%A2%E3%81%A8%E3%83%84%E3%82%A3%E3%83%AC%E3%83%AB%E3%82%BD%E3%83%B3%E9%99%90%E7%95%8C%E3%81%AE%E8%AA%AD%E5%87%BA%E3%81%97%E5%B0%8E%E5%87%BA_20261101

### The EPR Paradox Is Not a Paradox (supplement)

"The EPR Paradox Is Not a Paradox: A Supplement to 'Observation as Simultaneous Emission and Absorption' — Reading the Position–Momentum Pair as a Closure Identity"

Version: v1 (2026-11-01)

Concept DOI:
https://doi.org/10.5281/zenodo.23860605

Version DOI:
https://doi.org/10.5281/zenodo.23860606

English PDF (Zenodo):
https://zenodo.org/records/23860606/files/paper_epr_not_a_paradox_supplement_en_v1.pdf

Paper folder (Japanese and English manuscripts, figures and the program that draws them, numerical identity checks, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/EPR%E3%81%AF%E9%80%86%E8%AA%AC%E3%81%A7%E3%81%AF%E3%81%AA%E3%81%84_%E8%A6%B3%E6%B8%AC%E3%81%AF%E6%94%BE%E5%B0%84%E3%81%A8%E5%90%B8%E5%8F%8E%E3%81%AE%E5%90%8C%E6%99%82%E6%88%90%E7%AB%8B_%E8%A3%9C%E9%81%BA_20261101

The Japanese version of each paper is the authoritative text; both are included in the same Zenodo records.

### Japanese version of this article

https://note.com/kiharanoriaki/n/n40a43ad1e1ed

### Earlier articles (in Japanese)

"Does a wave packet shrink by observation, or gather by interaction?" (2026-07-13)
https://note.com/kiharanoriaki/n/nbc6649e30af3

"Reading the white cat, the black cat and the gray cat in a closed wave system" (2026-07-14)
https://note.com/kiharanoriaki/n/n2901bcd95bba

### Previous article

"Did We Understand Complex Numbers Correctly? — Ask Where the Amplitude Lives, and the Born Rule and the Invariance of the Speed of Light Become Two Sections of One Axiom"
https://note.com/kiharanoriaki/n/n64ad85f070d8

#TheoreticalPhysics #QuantumMechanics #EPRParadox #Entanglement #MeasurementProblem #BellInequality #TsirelsonBound #BornRule #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
