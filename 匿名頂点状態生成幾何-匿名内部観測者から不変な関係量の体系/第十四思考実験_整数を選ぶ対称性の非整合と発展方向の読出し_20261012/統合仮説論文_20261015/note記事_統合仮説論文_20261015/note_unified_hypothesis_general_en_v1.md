# Reading Out Mass, Charge, Velocity and Spin from Residents of a Circle, with No Time and No Space Put In - the Integrated Hypothesis Paper

The Standard Model and general relativity put time and space in first.

On that stage they place particles and fields, and give the values of mass and charge as inputs.

This paper reverses the order.

Without putting in any background time or space, it asks what can be read out, one stage at a time, from nothing but the relations and interference of complex waves on a closed circle.

The claim is a single one.

**Without a background spacetime, from the relations and interference of complex waves alone, which are not even travelling waves of odd and even harmonics, time, space, mass, charge, translational velocity and spin can all be read out.**

The paper takes two dialogues between me and an AI (Claude) as its primary sources, and is written in the exact order in which I asked, the AI erred, and I corrected it until the path was fixed.

This article follows the main points in that order. The figures are those of the paper; the interactive figures are HTML files placed on Zenodo.

## The setting: residents of a circle

Consider residents who live on a circle of radius R. No background time or space is assumed.

The residents cannot see the center of the circle, nor its inside, nor its outside. "Radius R" is a word used from outside; it is not a quantity the residents can measure from within.

The only clue left to the residents is this: go once around, and you are back where you started.

Now place a wave on this circle.

The direction in which the amplitude of a transverse wave would be placed is the direction perpendicular to the circle, but that is a direction of the background space, which the residents do not have. So the amplitude is placed in the imaginary direction.

The system has only two elements: the phase θ and the imaginary amplitude A.

Because the circle is closed, the allowed waves are discrete. Call a phase at which the amplitude is zero a node, and let n be the number of nodes in one lap. Since n counts half-wavelengths, it may be odd.

- n = 0 is the pair of states in which the whole circle carries a uniform +iA or −iA
- A wave with odd n flips its sign after one lap and returns to itself after two laps
- A wave with even n returns after one lap

The distinction between waves that return after two laps and waves that return after one has the same structure as the distinction between spinors and vectors, that is, between fermions and bosons.

Superpose the waves from n = 1 to 18 with equal amplitude, and a sharp peak stands at a single point. The appearance of the peak and the fixing of the origin θ₀ are the same event.

The sum of the odd harmonics alone shows its peak once every two laps; the sum of the even harmonics alone, once every lap. The sum of both cancels at the end of the first lap and returns at the end of the second.

[Insert image]

Figure file:
`fig00_circle_residents_waves.png`

Suggested caption:
"Top: waves drawn on the circle with the imaginary amplitude vertical. n = 0 is the pair of states +iA and −iA. n = 1 flips sign after one lap (blue) and returns on the second lap (orange, dashed). n = 2 returns after one lap. Bottom: sums of equal-amplitude harmonics drawn over two laps. Odd harmonics only (n = 1 to 17) peak once every two laps; even only (n = 2 to 18) peak once every lap; both together (n = 1 to 18) cancel at the end of lap 1 and return at the end of lap 2."

Interactive figure (opens in a browser; the maximum order and θ₀ can be moved with sliders):

Circle: imaginary amplitude (vertical) and its slope
https://zenodo.org/records/23386447/files/circle_imag_amplitude_radial_log_gradient.html?download=1

## One condition: the sum of true squares

Only one condition is imposed. The sum of the true squares of the components (squares without conjugation) equals R².

　Σxₙ² = R²

At first I ran a numerical experiment that solved this condition with absolute squares and a real radius. The results looked tidy, but they were wrong.

Read with true squares, R is imaginary; it is a quantity of the same kind as the components, an amplitude. It is not a radius placed outside. There is no outside of the circle.

Moving the term across gives

　Σxₙ² + (iR)² = 0

which is the condition of the light cone.

The residents cannot read curvature. But the length (the geodesic length) of one tick marked at equal intervals of phase differs from place to place.

For the full sum with N = 18, one tick in the interval of the peak is 1.245 times the average. For odd harmonics only, the interval opposite the peak (the antipode) shrinks to 0.835 times.

[Insert image]

Figure file:
`fig2_gauge_interval_lengths.png`

Suggested caption:
"Ratio of geodesic length per tick, for ticks at equal phase intervals (360°/N). Full sums N = 2, 4, 18; odd only N = 17; even only N = 18. Only the interval of the peak is long (1.245 times); for odd only, the antipodal interval is short (0.835 times)."

[Insert image]

Figure file:
`fig3_resident_ruler.png`

Suggested caption:
"The residents' ruler. Top: ticks at equal intervals of phase. Bottom: the same ticks placed by geodesic length. The same single tick has different lengths in different places."

## Stage 1: reading out mass

The waves so far are standing waves; no travelling wave can be made. The readout of mass began from this question.

The key is a phase that does not oscillate. A phase is assigned to each place and stays pinned there without moving.

The vacuum expectation value of the Higgs field in the Standard Model does not oscillate; it sits statically throughout space. The pinned phase has the same form.

When a composite wave (a particle) changes its position within this background, it picks up the difference from the pinned phase, a beat.

- A component that does not couple to the background picks up no phase difference and advances with ratio 1. This corresponds to light
- A component that couples to the background is slowed by the amount needed to realign its phase each time it changes position. This corresponds to mass

Three things are read from one beat at the same time. The strength of the beat that remains even at rest is mass; the direction and size of the beat's slope are the direction and speed of motion, that is, space; the number of beats accumulated while changing position is the time the particle experiences, its proper time.

Here the AI calculated that the phase velocity of the beat exceeds the speed of light. That calculation brought in a uniform ruler from outside. The residents' ruler is stretched where the peak is.

A check that compared only the widening of the spatial tick with ω/k, however, did not hold. The comparison was ill-posed.

The answer was the condition itself.

　The light cone Σxₙ² = 0.

Of the total rotation budget R, the part that rotates in spatial directions is k, and the part that goes into the beat with the pinned background is m.

　R² = k² + m²

The speed is k/R and the inverse Lorentz factor 1/γ is m/R. Time dilation and length contraction both come from this one ratio.

This was checked numerically on the circle with N = 18 components, placed with weight a for the right-handed rotation and b for the left-handed one.

- v² + (1/γ)² = 1 holds with a residual of 10⁻¹⁵. From the form of the coefficients it is an identity
- The speed is set by the logarithm of the weight ratio: v = tanh ln(a/b)
- Composing two speeds reproduces the relativistic velocity-addition law to six decimal places. A simple sum gives values above 1
- The mass 1/γ is read as the fraction of the cross term 2ab between the two counter-rotating terms

The background corresponding to the Higgs field did not have to be put in from outside. Taking true squares, it appears as the cross term of the counter-rotating terms.

[Insert image]

Figure file:
`fig_lightcone_readout_s1.png`

Suggested caption:
"Reading Σxₙ² = R² on the circle as a light cone. Left: v = k/R and 1/γ = m/R read from the field lie on tanh ln(a/b) and sech ln(a/b). Center: every point lies on the circle v² + (1/γ)² = 1, with the photon-like case (b = 0) and rest (a = b) at the two ends. Right: the speed read from the product of weight ratios agrees with the relativistic addition law, while the simple sum exceeds 1."

Interactive figures:

Circle: action of phase rotation and the light-cone split (translation slider, odd and even N separately)
https://zenodo.org/records/23386447/files/circle_phase_rotation_action_lightcone_v2.html?download=1

Circle residents: the special-relativity thought experiment (observer and particle)
https://zenodo.org/records/23386447/files/sr_observer_lightcone_s1.html?download=1

## Stage 2: why three dimensions

The residents of the circle are one-dimensional. Next consider residents of the four-dimensional spherical shell (S³). For residents living on its surface, the world is three-dimensional.

Line up three conditions.

- Uniform translations form a group with no fixed point. The only spheres that satisfy this are S¹ and S³
- Waves propagate without leaving a tail (Huygens' principle). This holds only in odd dimensions
- There is a curvature the residents can read, and a mass gap arises. This holds for every sphere except S¹

Only S³ satisfies all three at once.

On S³ a curvature term is attached to the advance of phase. The term that plays the role of mass comes from the very curvature the residents can read, and R is fixed by the wave itself through Σxₙ² = R².

Slicing S³, the even components raise a peak of the same sign at the point opposite the peak (the antipode), and the odd components raise one of the opposite sign. Only when both are present is the position fixed at a single point.

Interactive figures:

Slices of the four-dimensional shell S³: amplitude (shading, opacity) and slope
https://zenodo.org/records/23386447/files/s3_shell_slice_amplitude_gradient.html?download=1

Four-dimensional shell S³: action of phase rotation and the light-cone split
https://zenodo.org/records/23386447/files/s3_phase_rotation_action_lightcone.html?download=1

## Spin, charge, observation

For residents of S³ the following can also be read.

Spin. The odd/even distinction is directly spin 1/2 versus integer spin. Translation on S³ is multiplication by a quaternion from the left, which is the generator of spin 1/2. The cross term has spin 0.

Being localized at one point and being a spin eigenstate are incompatible. A spin eigenstate is localized at two points, and the opposite spin on the far side can be observed by no one.

Charge. S³ carries a bundle of circles called the Hopf fibres. Their winding number is the charge, and its integer nature gives quantization.

Observation. Particle and observer are both residents, and observation is the cross term between the two waves. Being readable, radiating, and having energy distributed are the same event. Uncertainty and quantization follow as consequences.

Interactive figure:

Four-dimensional shell S³: the special-relativity thought experiment with an observer
https://zenodo.org/records/23386447/files/sr_s3_observer_lightcone.html?download=1

This is the main result of the paper.

## Beyond: the electron, α, the four-dimensional lattice, E₈

What happens when the electron is read in this system?

The electron is chosen as the test case, and the fine-structure constant α is introduced. The numbers from the previous article are used here.

- The root 23/124, at which the exchange of two waves returns after 124 steps. It gives 137.036042915…, a value near α⁻¹
- The number of unit cells of a four-dimensional lattice inside a ball of radius 3, 137 = 1 + 8 + 128

The reason the background is four-dimensional is that four is the smallest dimension in which the diagonal has integer length (√4 = 2), and only there do 4-fold and 6-fold symmetry live together in one lattice. It is the same answer, seen from two sides, as why visible space is three-dimensional.

Imposing the rotation 23/124 and rounding to the lattice at the same time on the four-dimensional lattice, the state is captured by 23/124 as the radius grows, passing through the staircase of approximating fractions 1/6 → 3/16 → 5/27 → 13/70 → 23/124.

[Insert image]

Figure file:
`te14_E2_fig.png`

Suggested caption:
"The rotation 23/124 and rounding to the lattice imposed together on the four-dimensional Hurwitz lattice. As the radius grows, the state is captured by 23/124 through a staircase of approximating fractions. The α readout at the bottom right uses a withdrawn method of measurement and is not used in the paper."

In units G = c = q₀ = 1 (q₀ is one third of the electron charge), Planck's constant is an output, not an input.

　ħ = 9 α⁻¹

Once α is fixed, the electron mass follows automatically, and the only thing left to verify is α.

The locking rate to 23/124 rises to 40 → 70 → 86 percent as the radius grows, while the neighbouring roots stay at 0 to 27 percent. Roots whose denominator is a multiple of 4 lock more easily, and after a disturbance of one lattice unit, 18 of 20 orbits return.

[Insert image]

Figure file:
`te14_E4_fig.png`

Suggested caption:
"Selection, pull-in and recovery. (A) Radius dependence of the locking rate for 23/124 and its neighbouring roots; roots with denominator a multiple of 4 lock more easily. (B) Width of the pull-in to 23/124. (C) Recovery after a disturbance of one lattice unit, and the readout during the return (lower side)."

The paper reads the small deviation of the observed α from the value made from 23/124 as a consequence of observation being an average that includes states in transition. As I wrote in the previous article, being close but not identical is the content of this hypothesis.

The last symmetry that remains is E₈. E₈ can be written as two Hurwitz quaternions joined by the condition that their difference is a multiple of (1 + i). The antipodal point at 180 degrees is included automatically. Since the way of joining has been identified, the Lamb shift and the g-factor become derivation problems from here on.

## What was verified, what is hypothesis, and the means of verification

The paper keeps these three apart.

Verified (numbers obtained in the paper and in the numerical experiments).

- Odd closes in two laps, even in one. Superposed harmonics localize at a single point
- With true squares, R is imaginary. The geodesic length of one tick is 1.245 times in the interval of the peak
- On the circle the light cone holds as an identity, and velocity addition agrees to six decimals
- The only spheres with a group of uniform translations are S¹ and S³. Huygens' principle holds in odd dimensions
- On the four-dimensional lattice the pull-in rate to 23/124 reaches 86 percent
- Joining two Hurwitz integers gives the 240 roots of E₈

Hypotheses (connections placed during the dialogue).

- The pinned phase corresponds to the Higgs field
- Charge is the winding number of the Hopf fibres (the numerical separation of winding number and k has not been done)
- Observation is the cross term between particle and observer
- The reason for locking into that window is that integer-selecting structures lie nearby
- Everything from the beginning of the universe to the background radiation lines up as stages of the same system

Means of verification.

- The answer table for the hydrogen atom made in the Eleventh Thought Experiment; whether unrecorded exchanges give 8170 MHz and 1057.85 MHz (the Lamb shift)
- The difference between the measured α, 137.035999177, and the value from the g-factor of a single electron, 137.035999166
- Whether the low multipoles of the background radiation carry a systematic trace of 4-fold symmetry

## What to check next

- Obtain radiation, absorption, recoil and attraction in one calculation, as the cross term between two localized waves
- Separate the charge winding number from the momentum k numerically
- Enumerate the roots with denominator a multiple of 4 and search for those whose N(R) falls near 137
- Under the E₈ joining condition, rotate 23/124 and 101/124 together and see whether the pairs that can lock are narrowed down

## How the record was kept

The primary sources of this paper are the full texts of two dialogues between me and the AI. Transcripts are placed in the GitHub folder.

The order of the chapters is the order in which I asked, the AI erred, and I corrected. The AI's errors are not erased; they are kept, together with the corrections, as part of the derivation.

The first draft was written with α as its subject, and I rejected it. α is one stage of the development, not the subject of the paper.

This is not a finished theory. What was verified, what is hypothesis and the means of verification are kept apart, and the whole is fixed in a form that can be taken over and tested. It is an alternative interpretive hypothesis of the standard theory.

## About the paper

This article is based on the following paper.

### The integrated hypothesis paper

"Reading Out Time, Space, Mass, Charge, Translational Velocity and Spin from the Relations and Interference of Complex Waves without a Background Spacetime
― An Integrated Hypothesis in the Order of the Thought Experiment, from Residents of a Closed Circle to Residents of a Four-Dimensional Spherical Shell, and on to α, the Four-Dimensional Lattice and E₈ ―"

Version: v1.1 (2026-10-15. In the PDF of v1.0 some links to the interactive HTML figures in the appendix were missing, so a corrected version was published the same day. The text is the same.)

Concept DOI (always to the latest version)
https://doi.org/10.5281/zenodo.23384928

Version DOI (fixed to v1.1)
https://doi.org/10.5281/zenodo.23386447

English PDF (Zenodo)
https://zenodo.org/records/23386447/files/unified_hypothesis_paper_en_v1.1.pdf

Folder of the paper (Japanese and English manuscripts, figures, transcripts of the dialogues; GitHub)
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%8D%81%E5%9B%9B%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E6%95%B4%E6%95%B0%E3%82%92%E9%81%B8%E3%81%B6%E5%AF%BE%E7%A7%B0%E6%80%A7%E3%81%AE%E9%9D%9E%E6%95%B4%E5%90%88%E3%81%A8%E7%99%BA%E5%B1%95%E6%96%B9%E5%90%91%E3%81%AE%E8%AA%AD%E5%87%BA%E3%81%97_20261012/%E7%B5%B1%E5%90%88%E4%BB%AE%E8%AA%AC%E8%AB%96%E6%96%87_20261015

The Japanese full text is included in the same Zenodo record. The Japanese version is the original; the English version is a translation.

### Interactive figures (HTML; download and open in a browser)

Circle: imaginary amplitude (vertical) and slope
https://zenodo.org/records/23386447/files/circle_imag_amplitude_radial_log_gradient.html?download=1

Circle: action of phase rotation and the light-cone split (first version)
https://zenodo.org/records/23386447/files/circle_phase_rotation_action_lightcone.html?download=1

Circle: action of phase rotation and the light-cone split (translation slider, odd and even N separately)
https://zenodo.org/records/23386447/files/circle_phase_rotation_action_lightcone_v2.html?download=1

Circle residents: the special-relativity thought experiment (observer and particle)
https://zenodo.org/records/23386447/files/sr_observer_lightcone_s1.html?download=1

Slices of the four-dimensional shell S³: amplitude and slope
https://zenodo.org/records/23386447/files/s3_shell_slice_amplitude_gradient.html?download=1

Four-dimensional shell S³: action of phase rotation and the light-cone split
https://zenodo.org/records/23386447/files/s3_phase_rotation_action_lightcone.html?download=1

Four-dimensional shell S³: the special-relativity thought experiment with an observer
https://zenodo.org/records/23386447/files/sr_s3_observer_lightcone.html?download=1

### The Japanese version of this article

https://note.com/kiharanoriaki/n/ne3480f0f3931

### The previous article

"Why Don't the Three 137s Agree Exactly? - Search for Pairs That Are Close but Never Lock: the Fourteenth Thought Experiment"
https://note.com/kiharanoriaki/n/n2b247f7b0bc5

#TheoreticalPhysics #ThoughtExperiment #LightCone #HiggsField #Spin #FineStructureConstant #137 #E8 #QuantumMechanics #Relativity #IndependentResearch #OpenScience #Preprint
