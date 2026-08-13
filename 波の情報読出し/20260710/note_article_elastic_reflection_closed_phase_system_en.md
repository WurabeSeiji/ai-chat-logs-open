# Building particle-like collision and observation from complex waves and interference alone

In ordinary physics, we usually begin with space and time.

Then we place particles in that space. The particles move in time. They collide. A measuring device observes what happened.

For example, two particles, A and B, approach each other, undergo a perfectly elastic collision, and bounce back. A measuring device C observes the before and after.

That is a very natural picture.

But in this experiment, I wanted to avoid assuming that picture from the start.

In other words, I did not want to begin with:

- a pre-existing background space,
- particles already placed inside it,
- motion through an external time,
- and an observer looking in from outside.

Instead, I tried to construct A, B, C, and even the things that look like position and time, from complex waves, phases, and interference.

This is the constructive experiment reported here.

The formal paper and data have been published on Zenodo.

- Concept DOI, always pointing to the latest version: https://doi.org/10.5281/zenodo.21291018
- This version: https://doi.org/10.5281/zenodo.21291020

## What was built?

There are three waves in the simulation.

A and B are localized particle-like waves. "Localized" means that they are not spread uniformly everywhere, but have a peak around a certain phase position. They behave like compact lumps of wave.

C is the observer.

But C is not a metal device placed in an already existing space. C is also a wave. It is implemented as a heavy observing wave with a much larger representative amplitude than A and B.

All three, A, B, and C, are described using complex phases.

When a wave is written as a complex number, it has not only magnitude but also phase: a direction or shift in the complex plane. When two waves are combined, aligned phases strengthen each other, and shifted phases weaken each other. This is interference.

In this simulation, observation is a readout through interference.

The wave C is correlated with A. The wave C is also correlated with B. From those correlations, we read:

- where A appears to be,
- where B appears to be,
- which direction A and B appear to move,
- and whether A and B can still be distinguished after the collision.

The important point is that I do not externally declare, "this coordinate is position" or "this variable is time." Position-like and time-like quantities are read as relative phase relations with the observing wave C.

## How can A and B be distinguished?

If A and B are merely variable names in a program, that is too weak.

So I gave A and B different internal identification oscillations.

Intuitively, A has its own internal rhythm, and B has a different internal rhythm.

Mathematically, using an internal identification phase eta, the identification oscillation is written as:

```text
D_m(eta) = exp(i m eta)
```

A and B are assigned different values of `m`.

This means that after the collision, we can ask: is the wave on the left still A, or is it B? The answer is not obtained from the variable name, but by reading the internal oscillation carried by the wave itself.

This is crucial.

If two particle-like objects approach and then separate, their positions alone can be ambiguous. Several interpretations may look similar:

- A and B passed through each other.
- A and B reflected from each other.
- the labels A and B were exchanged.

By reading the internal identification oscillations, the simulation distinguishes reflection from transmission and label exchange.

## How was the elastic collision constructed?

The collision is not treated as an event at a single mathematical point.

Instead, I define a finite-resolution interaction cell. A and B are judged to enter the cell when they are sufficiently close in both spatial phase and temporal phase.

Inside that cell, the direction readout is reversed.

That is:

- A changes from right-moving to left-moving,
- B changes from left-moving to right-moving.

But only the direction is allowed to change.

The internal identification oscillation of A remains A. The internal identification oscillation of B remains B. The representative amplitudes are preserved. The fermion-like internal core is also preserved.

In this precise sense, the collision is constructed as a perfectly elastic reflection.

## What was confirmed?

In the minimal experiment, A and B reached the interaction cell, reversed their direction readouts, and then separated again.

When read out by the observing wave C, the identification oscillations of A and B were preserved.

So the result was not merely that two position curves moved in a plausible way. Rather, it was read as:

"the wave carrying A's internal oscillation reflected, and the wave carrying B's internal oscillation also reflected."

I also tested control maps:

- simple transmission,
- reflection with label exchange,
- transmission with label exchange.

Only the reflection map satisfied both direction reversal and identification-mode preservation.

I also varied the conditions:

- when the observer C was too light,
- when the update step skipped over the finite cell,
- when the temporal phases did not align,
- when the identification oscillations leaked into each other,
- when the sampling resolution in the internal eta phase was too low.

In these cases, the model failed in clear and separate ways.

That is important. This is not a model that succeeds no matter what. It has conditions under which it works, and conditions under which it breaks.

I also tested repeated collisions. Even after 8 collisions, direction reversal, identification modes, representative amplitudes, and the compensated square-closure condition were preserved.

## Why is this interesting?

The interesting point is not simply that a collision was simulated.

The interesting point is earlier than that.

Without assuming an external background space first, the simulation constructs:

- particle-like localized waves,
- an observing wave,
- position-like and time-like readouts,
- internal identity,
- collision-like interaction,
- and measurement-like correlation readout,

using only complex waves, phases, interference, and internal oscillations.

Instead of saying:

"Particles collided in space,"

the construction says:

"From phase relations and interference among complex waves, a structure can be read as two identifiable localized waves undergoing perfectly elastic reflection."

This does not replace standard physics. It is not a derivation of standard scattering theory. It is a constructive experiment asking how much of the language of particles, position, time, collision, and observation can be rebuilt from a smaller set of ingredients.

## What this does not claim

To avoid misunderstanding, let me state the limits clearly.

- This does not derive standard fermion scattering.
- This does not derive the standard quantum measurement process.
- This does not compute a real physical scattering cross section.
- The direction-reversal rule itself is introduced as a construction rule inside the finite interaction cell.

What was shown is more modest, but still meaningful:

the reversal rule is compatible with complex waves, phase relations, interference, internal identification oscillations, and readout by an observing wave C. It is also distinguishable from transmission and label exchange.

## Summary

In this simulation, two fermion-like localized waves A and B, an observing wave C, and even the position-like and time-like structure in which they appear to exist were not given as an external background.

They were constructed from complex wave phases and interference.

Then a perfectly elastic reflection between A and B was numerically executed, and observation by C confirmed direction reversal, identification-mode preservation, and representative-amplitude preservation.

With only a small set of axioms, plus simple complex waves and interference, it became possible to build particle-like existence, collision, and observation inside one closed system.

This is still a small constructive experiment. But I think it is an interesting step toward thinking about physics without placing background space at the beginning.

----

The formal paper, specification, experiment results, and reproducibility data are available on Zenodo.

- Concept DOI: https://doi.org/10.5281/zenodo.21291018
- Version DOI: https://doi.org/10.5281/zenodo.21291020
- Zenn article: https://zenn.dev/noriaki_kihara/articles/elastic-reflection-closed-phase-system
- Japanese note article: https://note.com/kiharanoriaki/n/n15451632027b

<!-- pdf-links -->
The paper PDFs can be downloaded directly from the public repository.

- back.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/.venv/lib/python3.9/site-packages/matplotlib/mpl-data/images/back.pdf
- elastic_reflection_closed_phase_system_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260710/elastic_reflection_closed_phase_system_en.pdf
- elastic_reflection_closed_phase_system_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260710/elastic_reflection_closed_phase_system_ja.pdf

Repository: https://github.com/WurabeSeiji/ai-chat-logs-open

#Physics #QuantumMechanics #ComplexNumbers #Waves #Interference #Observation #Simulation #IndependentResearch #Zenodo #Science
