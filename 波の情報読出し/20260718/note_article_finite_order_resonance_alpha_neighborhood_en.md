# Why Did Two Sharp Valleys Appear Near α⁻¹ = 137 and 128?

Japanese version:  
https://note.com/kiharanoriaki/n/na16b6a4e5ff2

The fine-structure constant α is approximately 1/137.

More precisely, its inverse α⁻¹ is about 137 in the low-energy regime. At higher energies, the electromagnetic coupling changes, and α⁻¹ moves toward the region around 128.

This article is about a numerical mystery: why did extremely sharp responses appear near these two famous values?

In the previous article, we introduced two numerical experiments. One produced behavior resembling wave-packet collapse. The other was the white-cat, black-cat, and gray-cat experiment.

Neither model was originally designed to reproduce the fine-structure constant.

Nevertheless, the same region of the exchange coefficient appeared in both models, even though they had been built for different purposes. When we later scanned the exchange coefficient across a wide range, extremely sharp responses appeared near regions associated with α⁻¹ = 137 and 128.

Their positions barely moved when we changed the initial conditions, the harmonic components, or even the experimental model.

Why did the responses always appear in the same places?

The central part of that mystery has now been solved.

## The Starting Point Was an Experiment That Looked Like Wave-Packet Collapse

In the first experiment, we repeatedly exchange-scattered a broad wave A with a localized wave B containing additional harmonics.

Part of the localization and harmonic structure initially carried by B moved into A. As a result, A also approached a localized state.

The purpose was to test whether behavior resembling observational wave-packet collapse could be represented as an interaction between waves.

One surprising result was the model's low sensitivity to the detailed harmonic structure.

Odd harmonics worked.

Even harmonics worked.

Their phases and amplitudes could also be changed.

As long as an additional harmonic structure existed, localization transfer appeared, and the central exchange-coefficient region did not move very much.

## The Same Region Appeared in the White-Cat, Black-Cat, and Gray-Cat Experiment

In the second experiment, we gave two states A and B the intuitive names white cat and black cat.

A state dominated by A was called a white cat. A state dominated by B was called a black cat. A metastable mixed state between them was called a gray cat.

This model tested how long the gray-cat state could be retained under weak readout and whether strong observation moved it toward the white-cat or black-cat state.

Its purpose and observables were different from those of the wave-packet localization experiment.

Even so, the same exchange-coefficient region again became important.

When the two experiments were compared afterward, they were found to share the same two-channel exchange-scattering kernel.

## A Full-Range Scan Revealed Two Main Responses

We varied the exchange coefficient R over a broad range and ran the same calculation at every point.

Several candidates appeared, but two responses stood out.

![Full-range scan of the exchange coefficient R](../20260715/system_B_full_R_sweep_full_range_depth_v1.png)

Under the correspondence formula used in this study, one response lies near the low-energy region α⁻¹ ≈ 137.

The other was a candidate region connected to the high-energy value near 128. When the search was later extended to higher orders, another exact root was found extremely close to α⁻¹ ≈ 128.946.

The crucial point is that neither 137 nor 128 was placed in the code as a search target.

The full-range scan specified only the minimum R, maximum R, and step size. The same update rule and evaluation function were applied to every grid point. The numbers 137 and 128, the fine-structure constant, and the exact roots identified later were not included in the evaluation function.

Yet the sharp responses appeared in the same places.

## The Cause Was an Operator That Returns Exactly After a Finite Number of Steps

To determine the cause, we separated the scattering operator into symmetric and antisymmetric components and analyzed its eigenvalues.

At each sharp response, the phase of the antisymmetric component returned exactly to its initial value after a finite number of iterations.

If one exchange-scattering operation is denoted by U, then a finite-recurrence point satisfies

Uⁿ = I.

Here, I is the identity operator, which leaves a state unchanged.

In other words, after n scattering steps, the entire operator has completed an exact cycle and returned to the identity.

This condition shows that the recurrent exchange weights form the discrete sequence

Rₙ,ₘ = cos²(πm/n),

where n and m are integers.

The exchange coefficient R had appeared to be a continuously adjustable parameter. In fact, it contained exact closed orbits indexed by pairs of integers.

## Why Did the Position Stay Fixed When the Initial Conditions Changed?

This result answers the largest question left by the earlier experiments.

The response position was not determined by the initial state.

It was not determined by what was initially placed in waves A and B. It was determined by whether the exchange operator U itself satisfied

Uⁿ = I.

Therefore, as long as the same operator is used, changing the initial waveform, harmonics, phases, or amplitudes does not move the finite-order root.

The wave-packet-collapse-like experiment and the white-cat, black-cat, and gray-cat experiment use different later-stage normalizations and observation rules. Before those later operations, however, they share the same exchange-scattering kernel.

That is why the same exchange-coefficient region appeared in different experiments.

## Is It a Valley or a Peak?

Some descriptions call the response a sharp valley, while others call it a sharp peak.

They refer to the same phenomenon.

Let ε be a residual measuring the deviation from exact recurrence. At an exact root,

ε = 0.

If ε is plotted directly, the root appears as a deep valley.

For visibility, however, we defined the depth as

d = −log₁₀ ε.

As ε approaches zero, d becomes larger. In this representation, the same point appears as a high peak.

In ideal arithmetic, ε is exactly zero, so the depth d theoretically tends to infinity.

This is not a divergence of energy. It means only that the recurrence residual becomes zero and the displayed depth becomes unbounded.

![Two exact roots calculated with 80-digit arithmetic](two_physical_roots_multiprecision_v1/80digit_delta1e-16_comparison_v1.png)

When the calculation precision was increased to 50 and then 80 digits, the response positions did not move. The valleys became narrower, while their peaks in the depth representation became higher.

## What We Found Was Not α as a Scattering Cross Section

The explicit motivation for pursuing this investigation was the hope that it might explain the fine-structure constant.

We considered whether the exchange-scattering strength could be connected to α as an electromagnetic coupling or scattering-cross-section parameter.

However, the sharp responses identified here were not directly α as a scattering cross section.

Their cause was finite-order recurrence of the exchange operator.

These two mechanisms have not been derived as identical.

Therefore, this study does not claim to have derived the fine-structure constant.

That does not make the discovery disappear.

The cause of the sharp, initial-condition-independent responses has been identified analytically. A valley first observed unexpectedly in a numerical experiment has been traced back to an exact closure condition of the operator.

## The Mystery We Solved and the Mystery That Remains

The questions answered by this study are clear.

Why did different experiments produce responses at the same positions?

Why did those positions remain fixed when the initial waveforms and harmonic conditions changed?

Why did the responses become sharper as the numerical precision increased?

The answer is that the common exchange operator closes exactly after a finite number of iterations.

But a larger mystery remains.

When the finite-order roots are translated using the correspondence hypothesis of this study, why do they lie near the physically famous values α⁻¹ ≈ 137 and α⁻¹ ≈ 128.946?

Is this only an arithmetic coincidence?

Or do the fine-structure constant and finite-order recurrence share an underlying geometric, algebraic, or physical structure?

This has not yet been resolved.

We were looking for α.

Instead, we found finite-order recurrence that is independent of the initial state.

The next question is why its exact roots lie so close to α.

That remains an open problem.

## Paper and Data

Paper:

“Discovery of Finite-Order Resonance in Iterated Exchange Scattering: Identifying Sharp Peaks near Fine-Structure-Constant Inverse Values 137 and 128 with a Reproducible Wave-Packet Model”

Concept DOI, always pointing to the latest version:  
https://doi.org/10.5281/zenodo.21421366

Version DOI for this release:  
https://doi.org/10.5281/zenodo.21421367

The source code, precision-sweep data, figures, Markdown manuscripts, TeX sources, and PDFs are available in the public record.

#FineStructureConstant #137 #128 #Alpha #TheoreticalPhysics #MathematicalPhysics #WavePacket #WavePacketCollapse #FiniteOrderRecurrence #QuantumRecurrence #TwoStateSystem #UnitaryOperator #ExchangeScattering #NumericalSimulation #IndependentResearch #Preprint #Zenodo
