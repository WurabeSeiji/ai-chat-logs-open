# Does a Wave Packet Collapse by Observation, or Gather Through Interaction?

I have published a new paper.

The theme this time is wave-packet localization.

In theoretical physics, wave-packet collapse or localization remains one of the difficult questions.

Does a wave packet shrink because it is observed?

Does it become localized because it interacts with another system?

Or is it already somewhat artificial to separate observation and interaction as two different physical processes?

In this experiment, I tested this question numerically in a closed system that does not assume even background spacetime in advance.

The full paper, English translation, TeX/PDF files, execution script, figures, and result data are available on Zenodo.

Concept DOI:
https://doi.org/10.5281/zenodo.21333766

This version:
https://doi.org/10.5281/zenodo.21333768

Zenn article:
https://zenn.dev/noriaki_kihara/articles/exchange-scattering-localization-exchange

## What Was Tested

The question is simple.

If a broadly spread wave interacts with a localized wave, can the broad wave also become localized?

However, the waves are not placed inside an already existing space.

In this series, space and time are not assumed first. Instead, quantities that look like position, time, and acceleration are read from complex phases, interference, and readout waves.

This experiment follows the same rule.

There are two waves, A and B.

One is a broad, low-order wave.

The other is a higher-order, fermion-like localized wave.

The experiment asks whether localization can be transferred from the localized wave to the broad wave through their interaction.

## First, the Readout Floor

Before claiming localization transfer, I first checked the lower bound of the readout.

If the odd harmonics are reduced as far as possible, do the basic readouts still survive?

Can the position phase still be read?

Can the identifying internal mark still be read?

Does the acceleration-like readout remain stable?

The result was yes.

Even when the odd harmonics were reduced to the minimum, these basic readouts remained available.

This means that the basic information that looks like position or acceleration is not carried only by high odd harmonics.

The high harmonics mainly affect localization.

They sharpen the wave and make it look more particle-like.

## Localization Appeared in Intermediate Scattering

Next, only one side was given high odd harmonics.

A was the broad wave.

B was the localized wave.

These two waves were coupled by a fermion-like scattering matrix.

The important point is that this is not just a simple perfect reflection.

I compared three cases:

complete transmission,

complete reflection,

and intermediate scattering, where part is reflected and part is transmitted.

The result was clear.

With complete transmission, localization was not transferred.

With complete reflection, localization was not transferred either.

Localization exchange appeared only in the intermediate scattering case.

The most visible case is the time evolution at an intermediate reflection rate.

![Waveform evolution under recursive scattering at R=0.70](exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1/exchange_scattering_matrix_R070_waveform_evolution_v1.png)

The upper-left panel is the initial state.

The blue wave A is broad.

The orange wave B is sharply localized.

After one interaction, A becomes sharply localized.

After the second, third, fifth, and tenth interactions, localization does not simply move in one direction. It appears to oscillate between A and B.

After enough recursive interactions, the effective harmonic order and localization index of the two channels become close.

This is not a case where the broad wave simply disappears.

It is also not a case where only the localized wave remains.

The natural reading is that localization is exchanged and redistributed inside a closed two-channel scattering system.

## The Readout Apparatus Was Not the Main Cause

I also checked the effect of observation.

In this series, the readout apparatus is not an outside human observer.

It is another wave placed inside the system.

The target wave and the readout wave interfere, and information is read from that interference.

If the readout apparatus were the main cause of localization, stopping the readout should change the result.

In this experiment, it did not.

Within this model, the localization was therefore better understood as a result of the A-B interaction itself, rather than as an effect caused by the readout apparatus.

## The Boson-Like Control Did Not Localize

I also tested a boson-like complete transmission model.

In that case, the waves overlap.

They also interfere.

But after the interaction, localization did not transfer to the other side.

The same was true for a simple complete reflection model.

Localization exchange did not occur merely because two waves overlapped.

It also did not occur merely because the waves reflected.

In this model, localization redistribution appeared when a broad wave and a localized wave interacted through fermion-like intermediate scattering.

## What Was Found

The result can be summarized carefully as follows.

Even in a closed complex-phase system that does not assume background spacetime first, localization exchange can be numerically constructed while keeping the acceleration-like readout intact.

When a broad wave and a fermion-like localized wave interact through intermediate scattering, there are conditions under which localization is transferred to the broad-wave side.

This localization exchange did not appear under complete transmission or complete reflection alone.

Stopping the readout apparatus did not change the result, so within this model the main cause appears to be the interaction map rather than observation itself.

## What This Does Not Claim

This paper does not solve the general wave-packet collapse problem in quantum mechanics.

It does not completely solve the measurement problem.

It does not compute a real electron or fermion scattering cross section.

The claim is more limited.

In a closed complex-phase system, when a broad wave and a localized wave interact through fermion-like intermediate scattering, it is possible to construct a case where localization is transferred to the other channel.

And in this experiment, that effect is better read as the result of interaction itself, not merely as the effect of a readout apparatus.

That is the scope of the result.

## Summary

In this numerical experiment, I tested a phenomenon that looks like wave-packet localization inside a closed system made only from complex waves and interference, without assuming background spacetime first.

The result was that a broad wave can acquire localization when it interacts with a fermion-like localized wave through intermediate scattering.

It did not happen under complete transmission.

It did not happen under complete reflection.

It did not change when the readout apparatus was stopped.

These controls make the result fairly clean.

It is still a preliminary experiment.

But as an entry point for reading wave-packet localization not as a special act of observation, but as an interaction with a localized wave, I think this is a very interesting result.

----

The full paper, English translation, TeX/PDF files, execution script, figures, and result data are available on Zenodo.

Concept DOI:
https://doi.org/10.5281/zenodo.21333766

Version DOI:
https://doi.org/10.5281/zenodo.21333768

Zenn article:
https://zenn.dev/noriaki_kihara/articles/exchange-scattering-localization-exchange

#Physics #QuantumMechanics #WavePacketCollapse #WaveLocalization #ComplexNumbers #WaveInterference #Observation #Fermion #Simulation #IndependentResearch #Zenodo #Science
