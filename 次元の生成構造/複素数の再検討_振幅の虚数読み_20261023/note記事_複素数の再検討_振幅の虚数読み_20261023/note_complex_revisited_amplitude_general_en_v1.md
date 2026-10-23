# Did We Understand Complex Numbers Correctly? - Ask Where the Amplitude Lives, and the Born Rule and the Invariance of the Speed of Light Become Two Sections of One Axiom

In the previous article I introduced a crystallization hypothesis built with quaternions and octonions. This time, before going further, I go back.

Back to complex numbers.

My research reads out spacetime and particles from nothing but the relations of a set of nameless complex numbers. I have used the structured number systems in turn: from the reals to the complex numbers, from the complex numbers to the quaternions, from the quaternions to the octonions. Along the way, one doubt arose.

Have I overlooked a symmetry that is already there at the stage of complex numbers?

If I have, then a property found with quaternions or octonions is taken for something that first appeared at that stage. If the attribution is wrong, the answers to the questions "why quaternions?" and "why octonions?" are wrong too.

The check led to three claims.

**The amplitude cannot live in the radial direction.** It lives in an unreadable direction, orthogonal both to the circle and to the radius. This is the origin of the Q axis in the spacetime equation I have been using.

**If that amplitude is read as imaginary, the cos² of quantum mechanics (the Born rule) and the γ of relativity (the Lorentz factor) come out as two sections of one and the same expression.** What the normalization that makes it look like a probability throws away is the rapidity, a measure of velocity. The Born rule and the invariance of the speed of light are not two separate axioms; they are two sections of one axiom.

**What is added on the way up from complex numbers to quaternions and octonions is less than I thought.** What is added is space, and the distinction between mass and charge. Everything else was already there at the stage of complex numbers.

As mathematics, all of this consists of trivial identities. What is new is on the side of the physical reading. In this article I follow the motivation in order, with figures.

## There is only one rule of judgement

For the check I set one rule.

A property may be said to belong to the quaternions, say, only when it has been shown that it cannot be stated with complex numbers alone. It may be said to belong to the octonions only when it has been shown that it cannot be stated with quaternions alone.

For each claim, I first try to write it one stage down. If it can be written there, it belongs to the lower stage.

It is a simple rule, but when applied mechanically it sent several things that my earlier papers had credited to the quaternions back down to the complex numbers.

## 1. Where does the amplitude live?

Draw a complex number z in the complex plane, and the amplitude is simply the distance from the origin, the radius. I had never found that strange.

But if the state is taken to lie on a circle of radius R, then R is the radius of curvature of that circle, and ought to be a different quantity from the amplitude.

The radius of curvature cannot be seen by moving along the circle. It appears as an amount of rotation only after going once round a closed path that sweeps an area (holonomy). It is the same as the roundness of the Earth: you cannot tell by walking straight, only by the angular mismatch after walking round a triangle.

If the amplitude is placed in the radial direction, the curvature changes every time the amplitude changes. Two states of different mass would sit on spheres of different curvature and could not coexist at the same place.

[Insert Figure 1: fig2_amplitude_radius_vs_Q_axis_en.png]

Figure caption: Left, a circle in the plane. An amplitude placed radially is the same number as the radius of curvature R and cannot be distinguished from it. Right, a cylinder. The amplitude lives in the unreadable direction Q, orthogonal both to the circle and to R.

Only one place remains: the axis orthogonal both to the circle and to R, standing out of the plane of the circle. The circle becomes a cylinder.

This axis cannot be read from the directions that are readable as space. By my convention of writing unreadable directions with imaginary units, (iQ)² = −Q².

I have long been using an equation that combines spacetime and the internal space,

x² + y² + z² + (it)² + (iR)² + (iQ)² = 0.

Of its six terms, the three x, y, z are the tangential directions readable as space, and the three t, R, Q are the unreadable normal directions. The Q axis had been placed there before, but its origin became clear only now. Q is the direction in which the amplitude lives.

## 2. Amplitude and phase exchange places

Read the amplitude inside the exponent. Writing a complex relation u as u = e^w, we have w = φ + iθ, where φ is the logarithm of the amplitude and θ is the phase.

Multiply w by i, and φ and θ exchange places. Amplitude and phase are conjugate: they are exchanged by a 90° rotation in the exponent plane.

[Insert Figure 2: fig1_exponent_plane_duality_en.png]

Figure caption: The exponent plane w = a + iθ. Multiplication by i gives i·w = −θ + ia, exchanging the amplitude exponent and the phase.

This one rotation appears in many places in physics under different names. The conjugacy of particle number and phase. Time and frequency. The Wick rotation that connects the temperature of statistical mechanics with the time of quantum mechanics. The Laplace variable s = σ + iω of electric circuits.

Which axis is called the imaginary part is a choice of where to stop the axis in this rotation; it is not decided on the side of the numbers.

## 3. For some reason, the Born rule and the Lorentz factor come out together

Here is the centre of this article.

Put the exponent of the amplitude inside the same exponent as the phase, with an imaginary unit attached. That is what it means to read the amplitude as imaginary.

Then, for the hyperbolic functions of half of w = φ + iθ, the following identities hold.

|cosh(w/2)|² = (cosh φ + cos θ)/2

|sinh(w/2)|² = (cosh φ − cos θ)/2

Add the two and you get cosh φ; subtract and you get cos θ.

[Insert Figure 3: fig3_complex_angle_sections_gamma_en.png]

Figure caption: The two sections of w = φ + iθ. The section φ = 0 is cos²(θ/2), the form of the Born rule. The sum of the two is cosh φ independently of θ, and equals 1 only at φ = 0.

Look only at the section φ = 0, that is, when the amplitudes agree, and the two are cos²(θ/2) and sin²(θ/2), with sum 1. This is the form of the Born rule of quantum mechanics.

Look only at the section θ = 0, that is, when the phases agree, and the two are cosh²(φ/2) and sinh²(φ/2), with sum cosh φ. This is exactly the Lorentz factor γ of relativity. That a boost is a rotation by an imaginary angle, and that rapidities compose by addition, has been the reading since Minkowski.

One and the same complex exponential gives the probability of quantum mechanics when read on the real section, and the upper limit of velocity of relativity when read on the imaginary section.

Since cosh x = cos(ix), the circular and the hyperbolic functions are the same function read after a rotation by i. It is not strange that one of them appears; it was only natural that both appear.

And to make a general relation, one whose amplitudes do not agree, look like a probability, one has to divide by the sum, cosh φ. The one dimension discarded by that division is the rapidity. The moment you divide, the cause of the thinned-out interference becomes unreadable.

## 4. The denominator in the previous paper was the Born-rule normalization

In the previous paper I wrote the coupling strength between two states as k = (1 − u)/(1 + u). Computing its denominator gives

|1 + u|² = 2e^φ(cosh φ + cos θ) = 4e^φ|cosh(w/2)|².

The denominator that I declared last time as "divide by the composite q₁ + q₂" was the same division as the Born-rule normalization.

Here is why, in the previous paper, I insisted so stubbornly that division be done only once, at the end, with the denominator declared. Divide midway, and you cannot notice that the information φ has vanished.

[Insert Figure 4: fig4_k_selfdual_denominator_en.png]

Figure caption: Left, the coupling strength |k| = |tanh(w/2)|. The line where |k| = 1 is θ = ±π/2, independently of the amplitude mismatch φ. Right, the denominator |1 + u|².

Moreover, the coupling strength is exactly 1 when cos θ = 0, independently of the amplitude mismatch. Last time I identified this relation |k| = 1 with a relation on the light cone. It is also the point where the map of coupling and the map of distance exchange and coincide.

The light cone was the self-dual point of coupling and distance.

## 5. Which section to read, and where to divide, is an axiom

Up to here, everything is a result of algebra, not of physics.

What enters as physics is only the choice: which section to read, and where to divide to normalize. Quantum mechanics chooses φ = 0 and calls it probability; relativity looks at θ = 0 and calls it the upper limit of velocity.

Try to derive this choice from dynamics, and the dynamics is itself already written on some section, so you end up explaining the choice by presupposing it. As long as the observer is part of the system, this self-reference cannot be avoided.

So the choice is placed as an axiom. It cannot be derived.

But the two axioms that existing physics places separately, the Born rule and the invariance of the speed of light, become two sections of one axiom: which section of the same complex exponential to read, and where to divide. Two reduce to one.

## 6. In complex numbers, amplitude and phase cannot be told apart

Now I use the rule of judgement.

If, in complex numbers, the exponent of the amplitude is read as imaginary, we get e^{i(a + θ)}, and a is simply added to θ. An amplitude read as imaginary cannot be told from a shift of the phase. There is only one imaginary unit.

So what remains of the imaginary reading at the stage of complex numbers is only the reversal of the sign of the square (the sign of the Wick rotation). That amplitude and phase cannot be told apart is itself a property of the complex stage.

## 7. Quaternions: they can be told apart, and space is born

In quaternions the amplitude can be given an imaginary unit m different from that of the phase. Since m and n are different, orthogonal directions, amplitude and phase can be told apart.

Moreover, the rotation by m and the rotation by n do not commute.

e^{ma}e^{nθ} − e^{nθ}e^{ma} = 2 sin a sin θ (m × n)

The two imaginary readings, of amplitude and of phase, generate through the order of multiplication a third direction m × n. The three fill the three-dimensional imaginary part, that is, space.

[Insert Figure 5: fig5_quaternion_commutator_third_axis_en.png]

Figure caption: In quaternions the rotations about the amplitude axis m and the phase axis n do not commute, and the commutator generates a third direction. The axis of the composite rotation reverses the sign of its m × n component with the order of multiplication.

This is the root of the coupling strength becoming a three-dimensional vector in the previous paper. The core of that paper, that two standing waves can memorize a three-dimensional position, was a consequence of the non-commutativity of the two imaginary readings, of amplitude and of phase.

Let me add one honest note. The coupling strength of the previous paper, k = −n tan(α/2), is the sign-reversed version of a known map, the modified Rodrigues parameters, long used in attitude dynamics in aerospace. The rediscovery of the map, and its positioning as a principle for reading out position without a background space, are written separately.

One more. Until the previous paper I had written that the light cone comes from the true square of quaternions. Reviewed with the rule of judgement, the sign reversal is already in the complex numbers; what the quaternions add is that the light cone is not a pair of points but a cone with an orientation.

## 8. Octonions: mass and charge are told apart

In quaternions, even if the amplitude is given two units, one in the direction of R and one in the direction of Q, two imaginary units always fall within an associative subalgebra, and which of them is called mass is only a renaming.

In octonions, the associator first stands with the triple of R, Q and phase: (ab)c and a(bc) no longer agree. There, an R-type amplitude and a Q-type amplitude become algebraically different things.

Why can mass and charge be measured as distinct? Are they separate entities, or different readouts of the same thing? This question receives an algebraic answer: because the associator first stands with the third imaginary reading.

The spacetime-side quaternion a = t + xi + yj + zk and the internal-side quaternion b = R + Q₁i + Q₂j + Q₃k repeat the same 1 + 3 split.

Here I adopted the reading that mass is the spacetime curvature R itself and charge is the internal curvature Q itself. That mass curves spacetime is general relativity itself, and charge as internal geometry is a Kaluza–Klein type of reading.

Two readings, not yet checked, follow from this. Charge is the direction and sign of the internal vector, hence discrete and quantized. Mass is the magnitude of a scalar, and can take continuous values. This agrees with charge appearing only in integer multiples of the elementary charge, while mass has a hierarchy but no visible rule of quantization.

The hypothesis written in the previous paper, "mass is the square of the amplitude," is refined by this reading to "mass is the scalar part R of the internal side," and the Q part of the amplitude is read not as mass but as charge. I stated this in the paper explicitly as a correction of the previous one.

## 9. The attribution ledger

All of this is collected in a ledger sorted by the rule of judgement.

[Insert Figure 6: fig6_attribution_ledger_en.png]

Figure caption: The attribution ledger. From left: what is already in the complex numbers, what the quaternions add, what the octonions add.

Already in the complex numbers: the duality of amplitude and phase; imaginary amplitude and the light cone; cos² and γ; the self-dual point; finite order Uⁿ = I and the counting of m/n; the division trap; that amplitude and radius of curvature are different things.

Added by the quaternions: the distinction of amplitude and phase; the third direction, that is, space; the light cone with an orientation; that two standing waves memorize a position; that condensates can be made.

Added by the octonions: the associator; the distinction between mass-type and charge-type amplitude; the split form and |a| = |b|; that the order of binding changes the result.

I also recorded in the paper three cases where I myself nearly misattributed during the check. I was about to write the orthogonality of amplitude and curvature as a property of the four-dimensional shell. I was tempted to read γ as a consequence of the quaternionic signature. And I had credited the light cone to the quaternions.

## 10. Prior work, and what only this paper has

Each of the parts of this paper has prior work.

Reading the phase as a rotation in a plane of space is Hestenes' real Dirac theory. Reading a boost as a rotation by an imaginary angle is Minkowski and Penrose–Rindler. That the choice of imaginary unit selects the kind of mechanics was argued by Kisil. The exchange of temperature and time is Connes–Rovelli. The conjugacy of number and phase is Dirac, Susskind–Glogower and Pegg–Barnett. The standpoint of reading probability as a counting of rational numbers is Palmer and 't Hooft. The map of the coupling strength is the modified Rodrigues parameters. The attempt to reconstruct geometry from relations is the Penrose–Szabados line. Octonions and the Standard Model are Singh and Furey.

In the paper these are cited at the stage where each belongs, and none is presented as a finding of this paper.

What only this paper has are the following four points.

That cosh φ = γ is the one dimension thrown away by the Born-rule normalization, and that the Born rule and the invariance of the speed of light are two sections of one axiom.

That the self-dual point of the duality of amplitude and phase is the light cone.

That raising the imaginary reading from complex numbers to quaternions to octonions produces first the distinction from the phase and then the distinction between mass-type and charge-type, through the associator.

That the amplitude is housed in the normal direction Q of the spherical shell.

## 11. What remains, and the form of falsification

The identities used in the paper, 11 of them, were all checked numerically, and the programs are published with the paper.

The rule of judgement itself has the form of a falsification. If a readout of three directions can be built with complex numbers alone, the claim that space is added by the quaternions is wrong. If the distinction between mass and charge comes out algebraically with quaternions alone, the claim that it is added by the octonions is wrong. If a consistent reading can be built with the amplitude placed radially, Section 1 falls.

The remaining tasks are to actually write the coupling strength on the internal side and check how the mass ratio and the relation of charges connect with the spacetime side, to define the condensate, and the dynamics. Dynamics is not treated in this paper.

## Summary

I started from the naive question of where the amplitude lives.

The amplitude cannot live on the radius; it lives in an unreadable direction orthogonal both to the circle and to the radius. Read that direction as imaginary, and already at the stage of complex numbers the form of the Born rule and the Lorentz factor come out as two sections of the same expression, and the denominator of the previous paper turns out to be the same division as the Born-rule normalization. Which section to read and where to divide is an axiom, and two axioms reduce to one.

In complex numbers, amplitude and phase cannot be told apart. The quaternions tell them apart, and non-commutativity gives birth to space. The octonions split the amplitude into two types and tell mass and charge apart.

The answer to the question whether I had understood complex numbers correctly came out like this.

What lies at the stage of complex numbers is more than I thought; what the quaternions and the octonions add is more limited than I thought; and that limited part is space, and the distinction between mass and charge.

For details, please see the PDF on Zenodo.

## Paper information

### Did We Understand Complex Numbers Correctly?

"Did We Understand Complex Numbers Correctly? Where the Amplitude Lives, the Imaginary Reading, and Attribution across ℂ → ℍ → 𝕆"

Version: v1 (2026-10-23)

Concept DOI:
https://doi.org/10.5281/zenodo.23585436

Version DOI:
https://doi.org/10.5281/zenodo.23585437

English PDF (Zenodo):
https://zenodo.org/records/23585437/files/complex_revisited_amplitude_imaginary_en_v1.pdf

Folder of the paper (Japanese and English manuscripts, figures and the programs that make them, numerical checks of the identities, citation tables, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E8%A4%87%E7%B4%A0%E6%95%B0%E3%81%AE%E5%86%8D%E6%A4%9C%E8%A8%8E_%E6%8C%AF%E5%B9%85%E3%81%AE%E8%99%9A%E6%95%B0%E8%AA%AD%E3%81%BF_20261023

The Japanese original is also included in the same Zenodo record.

### Japanese version of this article

https://note.com/kiharanoriaki/n/n5a455168899b

### Previous article

"Placing Nothing, Spacetime Is Born, the Universe Inflates, and Mass and Charge Appear - A "Crystallization" Hypothesis Built with Quaternions and Octonions"

https://note.com/kiharanoriaki/n/nfe20fe0a7dbe

#TheoreticalPhysics #ComplexNumbers #Quaternions #Octonions #BornRule #LorentzFactor #LightCone #MassAndCharge #Spacetime #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
