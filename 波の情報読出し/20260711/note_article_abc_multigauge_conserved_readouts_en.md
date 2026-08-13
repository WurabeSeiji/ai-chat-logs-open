# Where are mass, momentum, and energy being read from?

I have published a new paper.

This time, the theme is mass, momentum, and energy.

However, this paper does not claim that standard physical mass, standard momentum, or standard energy have already been fully derived.

The question is one step before that.

If we do not assume a background coordinate system first, can quantities that look mass-like, momentum-like, and energy-like be read from interference inside a closed phase system?

This paper investigates that question numerically using a minimal ABC closed phase system.

The formal paper, English translation, TeX/PDF files, scripts, and result data are available on Zenodo.

- Concept DOI, always pointing to the latest version: https://doi.org/10.5281/zenodo.21308049

- This version: https://doi.org/10.5281/zenodo.21308050

- Zenn article: https://zenn.dev/noriaki_kihara/articles/abc-multigauge-conserved-readouts

- Japanese note article: https://note.com/kiharanoriaki/n/nd5d3777a6e48

## Do not place background coordinates first

In the Wave Information Readout series, I have treated space and time not as background coordinates given in advance, but as readout quantities from a closed phase system.

In other words, x, y, z, and t are not placed first.

Instead, spatial and temporal quantities are read from phase differences, interference, reference waves, and readout windows.

The question in this paper is the next one.

If space and time are readouts, can mass, momentum, and energy also be treated as readouts?

This is a dangerous question, in a good sense.

Mass, momentum, and energy are among the most basic quantities in physics.

Therefore, this paper does not say that standard physical quantities have already been derived.

The aim is more modest and more precise:

Can conserved readouts that behave mass-like, momentum-like, and energy-like be constructed inside a closed phase system?

## The ABC closed phase system

The model used in this paper consists of local wave A, local wave B, and observer C.

A and B are local waves that collide.

C is not an external observer.

C is a reference wave inside the same closed phase system, used to read interference correlations with A and B.

The important point is this:

A value obtained from a single gauge is not treated as a completed measurement.

One viewpoint alone is not enough.

A readout is accepted only when the same quantity is reconstructed stably across multiple readout gauges, multiple reference waves, and multiple readout windows.

I call this multigauge interference readout.

## What was read?

The experiment uses three readout quantities.

The correlation gradient in the spatial phase direction is read as:

p_read

This is not standard momentum itself.

It is a spatial phase-gradient readout that behaves momentum-like.

The correlation gradient in the temporal phase direction is read as:

E_read

This is not standard energy itself.

It is a temporal phase-gradient readout that behaves energy-like.

Finally, the amplitude-squared residual that remains stable across multiple gauges is read as:

R_read

This is not standard mass itself.

It is a stable residual readout that behaves mass-like.

## R is hard to measure

One interesting point in this paper is the role of R.

The spatial phase gradient p is relatively easy to read because it appears as reversal or relative-gradient change.

The temporal phase gradient E can also be read by changing the temporal readout window.

But R is different.

R is stable precisely because its variation is small.

At the same time, because its variation is small, it is hard to measure.

In this sense, a mass-like quantity is intrinsically difficult to read.

This is an interesting feature of the construction.

The component that changes significantly is read as t.

The component that remains stable is read as R.

From this viewpoint, even the distinction between t and R is not a fixed naming given in advance. It is determined by readout stability.

In this paper, that t/R separation was also checked inside the multigauge readout structure.

## What happened in the symmetric collision?

First, I tested a single collision of A and B with equal amplitudes.

The result was that p_read, E_read, and R_read were reconstructed from multiple gauges with very high precision.

The maximum errors were:

- p_read: 2.5202062658991053e-14

- E_read: 2.2315482794965646e-14

- R_read: 4.440892098500626e-16

Under this condition, p was reversed, while E and R were preserved.

I also tested eight repeated collisions.

The p reversal, E and R preservation, identification oscillation preservation, and compensated closure were maintained.

So this was not a one-shot coincidence.

The readout structure survived repeated collisions.

## Under asymmetric R, simple reversal fails

Next, I tested the case where A and B have different R values.

Under equal-amplitude conditions, a simple reversal:

q becomes -q

looks like reflection.

However, when R_A and R_B are different, this simple reversal breaks the conserved readout.

In particular,

R_A p_A plus R_B p_B

is not preserved.

This is not a failure of the model.

Rather, it is a diagnosis.

If R is read as a mass-like quantity, then the collision map itself must be generalized into an R-weighted form.

This was the important turning point in the experiment.

## The R-weighted conservation map

I then constructed a generalized collision map that preserves:

R_A p_A plus R_B p_B

and

R_A p_A squared plus R_B p_B squared.

This does not mean that standard momentum conservation or standard energy conservation was assumed in advance.

The question was whether the readout quantities R and p, reconstructed inside the closed phase system, preserve R*p and R*p^2.

Across eight asymmetric amplitude cases, the maximum conservation errors were:

- R*p conservation error: 2.3803181647963356e-13

- R*p^2 conservation error: 1.4086509736443986e-12

The same structure was also tested under non-unit and asymmetric phase gradients, same-direction catch-up collisions, repeated collisions, readout noise, and extreme R ratios.

The ratio R_B/R_A was swept from 0.015625 to 64.0.

In the integration summary, all nine experiments were valid.

Also, none of the judgments used only a single gauge.

## What became visible?

The result suggests that mass, momentum, and energy do not necessarily have to be placed first as external substantial quantities.

At least inside this ABC closed phase system:

- a spatial phase gradient becomes a momentum-like readout,

- a temporal phase gradient becomes an energy-like readout,

- a stable amplitude-squared residual becomes a mass-like readout,

- R*p behaves like a momentum-conservation readout,

- and R*p^2 behaves like a squared-quantity conservation readout.

All of these were constructed consistently from multigauge interference.

I think this is an important result.

The previous papers treated space and time as phase readouts.

This paper extends that direction toward mass-like, momentum-like, and energy-like conserved readouts.

## What this does not claim

To avoid misunderstanding, let me state the limits clearly.

This paper does not fully derive standard physical mass, momentum, or energy.

It does not rederive standard mechanics.

It does not quantitatively predict real particle collisions.

What was confirmed is this:

inside the ABC closed phase system, conserved readouts that look mass-like, momentum-like, and energy-like can be constructed from multigauge interference.

The correspondence map to standard physical quantities remains a future task.

However, the readout-side foundation for constructing that correspondence map has become much clearer.

## Summary

In this numerical experiment, p_read, E_read, and R_read were reconstructed from multigauge interference in an ABC closed phase system without assuming background coordinates first.

In symmetric collisions, p was reversed, while E and R were preserved.

Under asymmetric R, the simple q becomes -q reversal failed to preserve R*p.

After that diagnosis, I constructed a generalized collision map preserving R*p and R*p^2.

The same conserved readout structure was maintained under asymmetric amplitudes, asymmetric phase gradients, repeated collisions, noise tests, and extreme R ratios.

Mass, momentum, and energy need not be placed first as substantial quantities.

They may be read as conserved quantities from a closed phase system.

This paper is a first numerical constructive experiment in that direction.

----

The formal paper, English translation, TeX/PDF files, scripts, and result data are available on Zenodo.

- Concept DOI: https://doi.org/10.5281/zenodo.21308049

- Version DOI: https://doi.org/10.5281/zenodo.21308050

- Zenn article: https://zenn.dev/noriaki_kihara/articles/abc-multigauge-conserved-readouts

- Japanese note article: https://note.com/kiharanoriaki/n/nd5d3777a6e48

<!-- pdf-links -->
The paper PDFs can be downloaded directly from the public repository.

- ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_en.pdf
- ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_ja.pdf
- abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_en.pdf
- abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_ja.pdf
- abc_closed_phase_system_multigauge_conserved_readouts_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_system_multigauge_conserved_readouts_en.pdf
- abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf

Repository: https://github.com/WurabeSeiji/ai-chat-logs-open

#Physics #QuantumMechanics #ComplexNumbers #Waves #Interference #ConservedQuantities #Mass #Momentum #Energy #NumericalExperiment #Simulation #IndependentResearch #Zenodo #Science
