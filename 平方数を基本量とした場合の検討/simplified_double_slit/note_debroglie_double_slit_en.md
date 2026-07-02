# Is the Electron a Particle or a Wave? —— "reproducing" the de Broglie wavelength and Planck's constant in a double-slit thought experiment

The electron is a particle and yet also a wave. This is one of the most famous stories in quantum mechanics.

The relation that assigns a wavelength to a particle electron is de Broglie's:

　λ = h / p

where λ is the wavelength, p the momentum, and h Planck's constant. Here that h makes its appearance. This prediction that "the electron has a wavelength" was actually confirmed by electron diffraction off a nickel crystal in the Davisson–Germer experiment of 1927.

This article is a short report on "reproducing" the wave nature of the electron, as a thought experiment, in the simplest classical wave model of a double slit. The formal version with equations is published on Zenodo in Japanese and English.

・Concept DOI (always the latest): https://doi.org/10.5281/zenodo.21109902
・This version: https://doi.org/10.5281/zenodo.21109903

## What I was curious about

The electron's wavelength λ is on a scale much larger than the "classical size" of the electron. So what if we think of the electron not as a point but as a lump of wave with internal structure finer than its own wavelength? Would it still interfere?

So I modeled the electron not as a point but as a spatially localized wave. The recipe is simple: sum only the odd-numbered waves (fundamental, 3rd, 5th, ...).

　S(φ) = cos φ + cos 3φ + cos 5φ + …

This makes a localized lump of wave with a single sharp peak at the center. The more terms you add, the sharper (narrower) the peak becomes. I tested, in a simulated way, whether this "lump packed with structure finer than λ" can interfere in a double slit.

## What I found (1): survival

I computed the interference exactly, at physical scale (5 cm from source to slits, slit separation 5 micrometers).

It turned out that for the lump of wave to make clean interference, an "alignment" condition is needed in which the phases of the harmonics line up. At wavelengths that do not satisfy it, the lump collapses and no clear fringes form. Conversely, at wavelengths that do line up, it "survives" as a train of single sharp peaks.

What matters is the spacing of the surviving fringes. No matter how many harmonics you add, the fringe spacing is set by the fundamental wavelength λ alone, independent of the number of harmonics. What the harmonics change is only the sharpness (the resolution) of the peaks.

In other words, even if you pack in structure finer than λ, the wavelength scale observed in the interference stays at the fundamental wavelength λ. This is the main content of this thought experiment.

## What I found (2): Planck's constant

For each electron, fluctuating its position and wavelength a little each time (position by half a wavelength, wavelength by about one percent), I computed the interference exactly many times and read the wavelength back from the fringes.

The product p × λ of the momentum p and the read-off wavelength λ agreed beautifully with Planck's constant h. So h = pλ does "reproduce" h.

But let me be honest here. This agreement actually happens because we put h in at the start, as λ = h/p. It is a self-consistency (a near-tautological check). Put differently, the h we put in is simply coming back out. So this is not an independent derivation that re-measures h. It is only a simulated reproduction that the model closes consistently.

What is genuinely new lies in "finding (1)": that even when you pack in finer structure, interference survives at the fundamental wavelength.

## What this does not claim

To avoid misunderstanding, let me state this clearly.

・This does not derive the de Broglie relation λ=h/p from first principles (λ=h/p is the input we assume).
・This does not derive the uncertainty relation either.
・It is only a thought experiment / reproduction within a classical wave model, and does not replace quantum mechanics itself.

Even so, being able to check in a visible form the single point that "even with localized structure finer than the electron's wavelength, the observed interference survives at the de Broglie wavelength scale" is, I think, an interesting result.

――――
(For the formal version with equations and figures, see the Zenodo record above and the Zenn article https://zenn.dev/noriaki_kihara/articles/debroglie-localized-double-slit .)

#Physics #QuantumMechanics #Electron #deBroglie #PlanckConstant #DoubleSlit #Waves #IndependentResearch #Zenodo #Science
