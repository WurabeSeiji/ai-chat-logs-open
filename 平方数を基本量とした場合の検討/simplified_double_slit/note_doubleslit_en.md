# The Double-Slit Puzzle Without Observation or "Wave-Packet Collapse"? — A Particle-Like Wave Packet Interferes as a Wave and Lands as a Particle

The double-slit experiment is the most famous puzzle in quantum mechanics.

A single particle seems to pass through both slits at once, interferes with itself, and builds up a fringe pattern on the screen. Yet, looked at one shot at a time, the particle hits a single point. And the probability of that point matches the wave interference intensity (the square of the probability wave).

Here textbooks usually invoke two special rules:

- The wave suddenly shrinks to a point the instant you observe it (wave-packet collapse)
- The act of observation itself selects the outcome

What these two thought-experiment notes want to show is a much plainer picture.

Without invoking observation or wave-packet collapse at all:

- a particle-like "wave packet" (a localized wave),
- passes through the two slits and interferes as a wave,
- and lands on the screen again as a particle-like localized wave — a single sharp peak.

And on top of that:

- the position where it lands rides exactly on the probability-wave distribution (the cos² shape).

All of this falls out of nothing but exact geometric computation.

---

## Paper 1: A fluctuating source just "shifts" the fringe

Start with a plain point source (a single-wavelength wave).

When the source position fluctuates a little, the far-field fringe on the screen does not change shape — it just shifts left or right. The shift amount is nearly proportional to the source position.

So if the source position fluctuates with some distribution (say a centrally peaked cos² shape), then reading the fringe shift each trial and histogramming it reproduces the original cos² shape exactly.

The input (the source-position distribution) is carried over, shape-preserved, into the output (the fringe-shift distribution). That is all. We did not derive a new probability law; it is just a change of variables (a push-forward).

The key is how you read it:

- read the shift each trial and tally → the cos² shape appears
- conversely, accumulate the brightness of many trials onto one plate → just a blurred single fringe of reduced visibility (no shape; this is the classical theory of optical coherence, the van Cittert–Zernike theorem)

So the "shape of the probability" appears only when you read the per-shot position and stack them. It is not in a single photograph. Plain but essential.

---

## Paper 2: Make the source a particle-like "wave packet"

Now the main act. Instead of a single wavelength, make the source a "sharply peaked wave packet (localized wave)" built by neatly superposing odd harmonics — a particle-like, compact bump.

──────────────────────────
[Paste FIGURE 1 here]
Figure 1: A particle-like wave packet = a localized wave built from superposed odd harmonics. A single sharp peak at the center, zero at the ends. This is the source.
Image file: 平方数を基本量とした場合の検討/simplified_double_slit/fig_paper2_localized_wave_N17.png
──────────────────────────

Pass this through the double slit, and —

- it interferes properly as a wave, and
- lands on the screen again as the same sharp localized peak (a particle-like single fringe).

A sum of odd harmonics, yet projected localized without losing its shape. That is the interesting part.

And when the source position fluctuates, the position where that sharp localized peak appears again rides on cos² (the probability-wave distribution).

──────────────────────────
[Paste FIGURE 2 here]
Figure 2: When the source position fluctuates, the sharp localized peaks (green) from each position line up neatly on the probability wave cos² (yellow curve). No observation, no wave-packet collapse is used. This is exactly the picture: a particle-like localized wave interferes as a wave and lands, as a particle-like localized wave, at a position that follows the probability-wave distribution.
Image file: 平方数を基本量とした場合の検討/simplified_double_slit/fig_oddharm_fluct_L10_W5_lam1p0308_N17_dlam1.png
──────────────────────────

To sum up:

- a particle-like localized wave interferes as a wave and is projected as a particle-like localized wave
- the projected position follows the probability-wave distribution
- and throughout, no observation and no wave-packet collapse were used

The essential part of the double-slit "puzzle" is reproduced plainly — from wave interference plus geometry alone — without bringing in collapse or an observer.

---

## But it is not free (the honest part)

I wrote that breezily, but the moment we go to a localized wave, a "condition" and a "fragility" appear that the single wavelength did not have.

- To localize while preserving shape, the source-to-slit distance must be "aligned" to an exact integer (or half-integer) multiple of the fundamental wavelength. In a formula: √(L²+(W/2)²) = (m/2)·λ₀.
- This alignment is very sensitive: a mere 3% change in the wavelength separates the cases that interfere cleanly from those that scatter into a mess.
- The sharper you make the localization (more harmonics), the narrower the allowed wavelength tolerance, shrinking as 1/(2N).
- When the source moves off the center, the two path differences disagree, and the peak slips a little below the probability-wave curve.

And the decisive point: only the single wavelength (the simplest case) holds unconditionally, with none of these conditions or fragilities. In fact, running the localized-wave method back at "one wavelength" agrees with Paper 1 to machine precision.

So the net finding is actually the other way around:

- shape preservation is real, but conditional, and fragile
- the simplest case — "one wavelength" — is the robust special one

There is a clear trade-off between resolution (sharpness of localization) and robustness (tolerance to fluctuation).

---

## What this is, and what it is not

These are thought experiments (exact computations of a model), not measurements.

- The cos² (the shape of the probability wave) appears because the input distribution is carried over directly; we did not derive the squaring rule of probability or the origin of randomness itself.
- So we do not claim to have "derived the Born rule" or "solved the measurement problem."
- Nor do we overturn existing quantum mechanics or wave optics. This is only an organization of how far the double-slit picture can be reproduced in a plain form, without collapse or an observer.

Even so, that one can trace "particle-like localization → wave interference → particle-like localization → a position following the probability-wave distribution" using only wave interference and geometry, without invoking observation or wave-packet collapse, is, I think, one fresh angle for looking at the double-slit puzzle again.

---

## Source papers (Zenodo, CC BY 4.0, JA/EN + reproduction code)

- Paper 1 (positional fluctuation and the push-forward of the fringe shift): https://doi.org/10.5281/zenodo.21035808
- Paper 2 (interference of the localized wave; alignment condition and fragility): https://doi.org/10.5281/zenodo.21035831
- Figures and full computation code (GitHub): https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/平方数を基本量とした場合の検討/simplified_double_slit
- Japanese version of this note: https://note.com/kiharanoriaki/n/n65be6bf06c9b
