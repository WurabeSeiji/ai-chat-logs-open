# Maybe curvature did not disappear, but was woven into the closed wave itself

I have published a new paper.

This time, the theme is a little strange.

Normally, if a wave is placed in a curved space, it should feel the curvature. Its phase should shift. Its interference should change. Even a reflection condition that was perfectly aligned should develop some leakage.

So why can many physical models treat a sufficiently local region as flat, and still agree with observation with very high precision?

Is the effect of curvature simply small enough to ignore?

Or is it possible that a stable wave already contains the curvature effect inside its own internal phase structure?

This paper numerically investigates that question inside my Wave Information Readout series.

The central condition is the all-positive zero closure:

Sum of x_n squared equals zero.

In symbols, I write it as:

Sigma x_n^2 = 0

This is not a conjugate norm.

It is not the sum of |x_n|^2. Instead, each component x_n is squared directly, and all terms are added with positive signs. The total is required to be zero.

The formal paper and data are available on Zenodo.

- Concept DOI, always pointing to the latest version: https://doi.org/10.5281/zenodo.21304039

- This version: https://doi.org/10.5281/zenodo.21304040

- Zenn article: https://zenn.dev/noriaki_kihara/articles/curved-closure-stationary-wave

## What was investigated?

In the previous paper, I constructed perfect elastic reflection without directly commanding an external direction flip.

Instead of writing:

q becomes minus q

as an external rule, the reflection was generated from a fermion-like internal inverse-phase core and interference between direct and exchange paths.

In other words, reflection became a readout produced by internal inverse phase, direct and exchange paths, even-odd channel interference, and a relative phase of pi.

The question in this paper is the next one.

What happens if that exchange-interference reflection is placed inside a curved local cell?

If curvature produces a small relative phase leakage, does the perfect reflection break?

Or does the condition for remaining as a closed stationary wave absorb the leakage into internal phase and recover the reflection?

## Curvature is not removed at the beginning

In this paper, I do not assume that curvature can simply be ignored.

Instead, I explicitly introduce curvature-induced relative phase leakage.

Take a closure pair:

x_m and i x_m.

If curvature induces a relative phase leakage, called delta_K,m, the square sum of the pair generally no longer vanishes.

In words, the square of x_m plus the square of the phase-leaked i x_m leaves a residual proportional to 1 minus exp(i 2 delta_K,m).

So the curvature effect does appear.

The closure condition breaks, and transmission leakage appears in the reflection readout.

This point is important.

The result is not that curvature has no effect.

Curvature does have an effect. It is visible in the transient state.

## The surprise comes after that

However, when the internal phase is re-selected so that the curvature leakage delta_K,m is cancelled by an internal correction beta_K,m, the closure condition is recovered.

That is, when delta_K,m plus beta_K,m equals zero, the closed stationary wave is restored.

At that point, curvature has not disappeared.

Rather, the curvature-induced relative phase leakage has been absorbed into the internal phase configuration of the wave that can stably exist.

In other words, the curvature effect did not simply become invisible from the outside. It may have been woven into the internal condition required for the wave to exist as a closed wave.

That was the interesting part.

It suggests another way to read local flatness.

Instead of saying only, "the region is small, so curvature can be ignored," we can also ask whether the stable readout selects only stationary modes that are already closed including the curvature effect.

## What happened in the numerical experiments?

In the minimal experiment, when curvature-induced relative phase leakage was introduced in the transient state, the closure-pair RMS became 1.2319416790092972e-02, and the transmission leakage became 1.1503183254481797e-01.

So the reflection was no longer perfect.

After re-selection into a closed stationary wave, however, the closure-pair RMS returned to 9.4283259783636047e-19, and the transmission leakage returned to 0.0.

I also swept eight curvature relative phase models and seven correction freedoms.

For maximum curvature relative phase 1.2, the uncorrected case left maximum transmission leakage 1.6202719613622976e-01.

With full correction, however, the closure-pair RMS returned to 7.8949412793793227e-19, and transmission leakage returned to 0.0.

The same structure also appeared when the residual curvature phase was inserted back into the one-sided local exchange-interference scattering map.

Without correction, transmission leakage appears.

After re-selection as a closed stationary wave, the perfect-reflection readout is recovered.

## Why is this important?

The central condition is:

Sigma x_n^2 = 0.

Again, this is not the conjugate norm Sigma |x_n|^2.

Each component is squared directly, and all terms are summed with positive signs.

For the minimal pair, the condition becomes:

A^2 plus (iA)^2 equals zero.

No external negative coefficient is inserted. The sign reversal inside the square sum is generated internally by the complex phase i.

In this paper, that axiom was not merely a conservation condition.

It played at least four roles at once:

- a condition for the existence of nontrivial complex waves,

- a detector of curvature-induced relative phase leakage,

- a re-selection condition for closed stationary waves,

- and a stability condition for perfect-reflection readout.

One condition ended up doing all of that.

That is the largest surprise in this result.

## What this does not claim

To avoid misunderstanding, let me state the limits clearly.

This paper does not quantitatively predict real spacetime curvature effects.

It does not derive general relativity.

It does not replace standard quantum theory.

It also does not claim that curvature is always unobservable.

Actually, the opposite is true in this construction.

When curvature-induced relative phase leakage is introduced, closure residuals and transmission leakage appear in transient states.

What was confirmed is that, once the system is re-selected into a stable closed stationary wave, that leakage is absorbed into internal phase and the external readout recovers perfect reflection.

## Summary

In this simulation, an odd-harmonic complex wave placed in a curved local cell first loses closure when curvature-induced relative phase leakage is introduced.

Perfect reflection also develops leakage.

However, when the system is re-selected into a closed stationary wave satisfying Sigma x_n^2 = 0, the curvature-induced phase leakage is absorbed into internal phase, and the perfect-reflection readout is recovered.

Curvature was not erased.

It may have been woven into the internal phase configuration required for a closed wave to exist.

This is still a numerical constructive experiment internal to the axiom system, before constructing a direct correspondence to standard theory.

Even so, it suggests a way to think about local flatness not only as an approximation, but also as a readout of stable closed stationary waves.

----

The formal paper, English translation, TeX/PDF files, scripts, and result data are available on Zenodo.

- Concept DOI: https://doi.org/10.5281/zenodo.21304039

- Version DOI: https://doi.org/10.5281/zenodo.21304040

- Zenn article: https://zenn.dev/noriaki_kihara/articles/curved-closure-stationary-wave

- Japanese note article: https://note.com/kiharanoriaki/n/n2389460836cf

#Physics #QuantumMechanics #ComplexNumbers #Waves #Interference #Curvature #LocalFlatness #NumericalExperiment #Simulation #IndependentResearch #Zenodo #Science
