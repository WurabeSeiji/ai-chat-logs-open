# Where Does 1/137 Come From? - Rereading Time, Mass, Charge and Symmetry from the Light-Cone Equation: the Twelfth Thought Experiment

In the previous article, the Eleventh Thought Experiment, I made the answer key first: for four two-body systems including the hydrogen atom, I computed what general relativity and electromagnetism should give.

That computation took a few numbers as input.

One of them was the fine-structure constant α.

α ≈ 1/137

Anyone who has dipped into physics has probably seen this number at least once.

Last time I put it in as the measured value, without a second thought.

I had been saying that I start from unnamed states without a background spacetime, and yet I walked right past the most famous constant of all.

So this time I stopped and thought.

What is α, in the first place?

## α is, first of all, a ratio of speeds

In textbooks, α is introduced as the constant that expresses the strength of the electromagnetic force.

But if we go back in history, a simpler face appears.

In the Bohr model of the hydrogen atom, the ratio of the speed v of the electron on the innermost orbit to the speed of light c is exactly α.

　v / c = α ≈ 1/137

In other words, the electron in the hydrogen atom goes around at 1/137 of the speed of light.

1/137 can first be read as a ratio of speeds.

Once a ratio of speeds appears, relativity comes next.

The Lorentz factor, which tells how much the clock of something moving at speed v slows down, is

　γ = 1 / √(1 − v²/c²)

and putting v/c = α into it gives γ = 1 / √(1 − α²).

Seen in this form, it is natural that α appeared when Sommerfeld put relativistic corrections into the hydrogen atom.

## An equation I have liked for a long time

Let me tell you something a little personal here.

When I was in high school, I came to know the Lorentz transformation and the light-cone equation.

The Lorentz transformation was an unsightly formula, tangled with square roots and fractions.

The light-cone equation, by contrast, was surprisingly simple.

　x² + y² + z² − t² = 0

(written in units where the speed of light is 1.)

The places light travels to are the points that satisfy this equation.

If we attach the imaginary number i to time and write it as it,

　x² + y² + z² + (it)² = 0

and everything becomes a sum of squares.

Squares add up to zero. That is all the equation says.

I was drawn to this beauty, and I have remembered it ever since.

The Lorentz transformations are determined afterwards, as the transformations that do not change the form of this equation.

What comes first is the light-cone equation; the Lorentz transformations are merely the transformations that keep it.

This thought experiment is the story of following that order all the way through.

## Add one axis, and a slow particle also rides on the light cone

A particle moving slower than light travels inside the light cone.

It does not ride on the surface of the cone.

But there is a trick here.

If τ is the time ticked by the particle's own clock (its proper time), then

　t² − r² = τ²

holds, where r is the distance the particle has traveled.

Moving τ² to the left gives

　r² + τ² + (it)² = 0

This has exactly the same form as the light-cone equation, only with one more axis.

So a particle slower than light, too, rides on the surface of the light cone in the space where one proper-time axis has been added.

A right triangle appears, with r and τ as its two legs and t as its hypotenuse.

The Lorentz factor is nothing but a ratio of the sides of this triangle, γ = t / τ.

[Insert image]

Figure used:
`fig01_lightcone_central_projection.png`

Caption:
"The light cone and central projection are the same right triangle. Left: the ordinary light cone. The path of a particle slower than light (red) lies inside the cone. Center: adding one proper-time axis τ puts the same path on the surface of the cone. A right triangle forms with distance r (blue) and proper time τ (green) as legs and time t (red) as hypotenuse. Right: the same triangle appears in the picture of central projection: the distance on the tangent plane (blue), the radius of the sphere (green) and the distance from the center (red)."

## The same triangle was also the triangle of central projection

I have been writing papers using a geometry called central projection.

It is the operation of drawing a line from the center of a sphere to a point on a plane tangent to the sphere, and mapping the point to where that line meets the sphere.

The distance ℓ from the center to the point on the plane then satisfies

　ℓ² = |x|² + R_c²

where |x| is the distance on the plane and R_c is the radius of the sphere (the curvature radius).

Placed next to the light-cone triangle above,

- the distance r corresponds to the distance |x| on the plane
- the proper time τ corresponds to the curvature radius R_c
- the time t corresponds to the distance ℓ from the center.

Remove the imaginary i from the light-cone equation, and you get the Pythagorean equation of central projection.

The central projection I have been using all along turned out to have the same form as the light-cone equation.

## Eight axes, with the names taken off

From here we go one step further.

Let us split the content of the proper time τ into several more directions.

　x² + y² + z² + R² + Q₁² + Q₂² + Q₃² + (it)² = 0

Eight axes in all.

x, y, z are the three directions of ordinary space.

R and Q₁, Q₂, Q₃ are directions that cannot be measured directly the way space can.

This equation, too, is the light-cone equation, in eight dimensions. Only t carries the imaginary i.

What matters here is not to give the axes names from the start.

We do not decide in advance that t is time, R is mass and Q is charge.

We think that all there really is are eight directions without names.

The names time, mass and charge may have been attached afterwards to the way an observer read each direction.

In fact, if we read momentum as the gradient along each direction, a relation among the energy E, the momentum p, the mass m and the charge q comes out:

　E² = p² + m² + q²

When the charge is zero, this is E² = p² + m², familiar from relativity.

The idea of reading the momentum along a hidden direction as charge is also found in the five-dimensional theory of Kaluza and Klein, a hundred years ago.

## Directions you cannot read are seen only as a length

An observer cannot necessarily distinguish and read all eight directions.

For example, consider an observer who can read only t.

To this observer, the remaining seven directions cannot be told apart.

What it sees is only the length that gathers the squares of the seven,

　ρ² = x² + y² + z² + R² + Q₁² + Q₂² + Q₃²

The seven degrees of freedom have not disappeared.

Because they cannot be distinguished, they appear only as one length.

[Insert image]

Figure used:
`fig02_observers_symmetry.png`

Caption:
"Three observers read the same eight directions. Left: blue marks the directions that can be read directly, gray the ones that cannot be told apart. From the top: an observer reading the space x, y, z, an observer reading Q₁, Q₂, Q₃, and an observer reading only the time t. Right: the number of readable directions and the number of transformations that leave things unchanged when the indistinguishable directions are rotated (the dimension of the rotation group). The fewer directions can be read, the larger the symmetry."

## Not being able to tell apart: that was symmetry

Rotate the indistinguishable directions all the way around, and the length the observer sees does not change.

From the observer's point of view, before and after the rotation are the same.

In mathematics, this is called symmetry.

If we write the result of observation as P, and applying a transformation G to a state U still gives P(U) = P(GU), the observer cannot tell the two apart.

That is,

**symmetry is the mathematical expression of not being able to tell things apart**

can be read.

Symmetry may not be a rule imposed on nature afterwards.

There are directions the observer cannot name. That very fact is what we see as symmetry.

Conversely, when observation becomes more precise and directions that could not be told apart become distinguishable, the symmetry becomes smaller.

What physics calls symmetry breaking might be read as things that could not be distinguished becoming distinguishable.

## Split them into 3 and 2, and the shape of the Standard Model appears

Finally, consider an observer who, like us, reads the space x, y, z.

The directions it cannot tell apart are the five t, R, Q₁, Q₂, Q₃.

Suppose observation distinguishes these five only as two groups:

- the three Q₁, Q₂, Q₃
- the two t, R

Within the group of three nothing can be told apart, and within the group of two nothing can be told apart either.

The candidate symmetry that then appears is

　SU(3) × SU(2) × U(1)

This is the same shape as the symmetry of the Standard Model of particle physics.

SU(3) of color charge, SU(2) of the weak force, and U(1) connected with electromagnetism.

And this structure of splitting five into 3 and 2 has the same shape as the split in the grand unified theory (SU(5)) proposed by Georgi and Glashow in 1974.

[Insert image]

Figure used:
`fig03_split_3_2.png`

Caption:
"Splitting the five indistinguishable directions into three and two. The blue block (mixing within Q₁, Q₂, Q₃) is U(3), the orange block (mixing within t, R) is U(2). The gray blocks (mixing between the three and the two) can be told apart by observation. Adding the condition that aligns the overall phase gives the same shape as the Standard Model symmetry SU(3) × SU(2) × U(1)."

However, this is not a derivation.

It goes only as far as the shape of the symmetry appearing as a candidate.

## Back to α once more

Let us go back to α, where we started.

In the hydrogen atom, α = v/c.

At first it looked like a mere ratio of speeds.

But in the view developed so far, both v and c may be quantities obtained by reading an unnamed higher-dimensional gradient from different directions.

If so, α may be an angle or a ratio that arises when a higher-dimensional state is mapped onto a certain direction of observation.

Which direction over which direction it is has not been decided yet.

Here I honestly leave it unresolved.

## What I am not saying

Let me make clear what this thought experiment does not say.

- I am not saying that I derived the symmetry of the Standard Model. It goes only as far as the shape appearing as a candidate. The sum of squares (the light-cone equation) and the squared absolute values used in quantum mechanics are different things, and connecting them needs one more assumption.
- I am not saying that I derived the value of α (1/137). Nor have I decided which ratio of directions α is.
- The eight axes, the axis that carries the imaginary i, and the split into 3 and 2 were put in by hand as assumptions.
- This is an organization of ideas, without computation.

## What I found

I found three things.

First. 1/137 can first be read as the ratio of the speed of the electron in the hydrogen atom to the speed of light.

Second. The Lorentz factor of relativity is nothing but a ratio of sides of the light-cone equation with one axis added. A particle slower than light, too, rides on the light cone once one axis is added. And that triangle was the same as the triangle of the central projection I have been using all along.

Third. If we give the axes no names and think only about which directions an observer can distinguish, symmetry appears as the expression of not being able to tell things apart. Splitting into 3 and 2, the same shape as the Standard Model appears as a candidate.

The equation in which squares add up to zero, which I have been calling the "zero closure", was not a special hypothesis.

It was the light-cone equation itself.

Next, with this view, I would like to check by computation whether, keeping the same rule of interaction and changing only the directions that are read, a really different physics comes into view.

## About the paper

This article is based on the following paper.

### The twelfth thought experiment

"The Twelfth Thought Experiment: From the Fine-Structure Constant to Observation Maps and Gauge Symmetry
- Starting from the Fine-Structure Constant, Rereading Spacetime, Mass, Charge, Color Charge and Gauge Symmetry from the Light-Cone Equation and the Choice of Observable Axes -"

Version: v1.0 (2026-10-08)

Concept DOI (always to the latest version)
https://doi.org/10.5281/zenodo.23226258

Version DOI (fixed to v1.0)
https://doi.org/10.5281/zenodo.23226259

English PDF (Zenodo)
https://zenodo.org/records/23226259/files/thought_experiment_12_observation_mapping_gauge_symmetry_en_v1.0.pdf

English PDF (direct download from GitHub)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%8D%81%E4%BA%8C%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E5%BE%AE%E7%B4%B0%E6%A7%8B%E9%80%A0%E5%AE%9A%E6%95%B0%E3%81%8B%E3%82%89%E8%A6%B3%E6%B8%AC%E5%86%99%E5%83%8F%E3%81%A8%E3%82%B2%E3%83%BC%E3%82%B8%E5%AF%BE%E7%A7%B0%E6%80%A7%E3%81%B8_20261008/thought_experiment_12_observation_mapping_gauge_symmetry_en_v1.0.pdf

Folder of the paper (figures and the program that generates them; GitHub)
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%8D%81%E4%BA%8C%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E5%BE%AE%E7%B4%B0%E6%A7%8B%E9%80%A0%E5%AE%9A%E6%95%B0%E3%81%8B%E3%82%89%E8%A6%B3%E6%B8%AC%E5%86%99%E5%83%8F%E3%81%A8%E3%82%B2%E3%83%BC%E3%82%B8%E5%AF%BE%E7%A7%B0%E6%80%A7%E3%81%B8_20261008

The Japanese full text and the three figures (SVG and PNG) are included in the same Zenodo record. The Japanese version is the original; the English version is a translation.

### The Japanese version of this article

https://note.com/kiharanoriaki/n/nb193975c2a11

### The previous article

"General Relativity and Electromagnetism Could Be Computed without a Background Spacetime - Making the Answer Key First, for Four Two-Body Systems Including the Hydrogen Atom: the Eleventh Thought Experiment"
https://note.com/kiharanoriaki/n/n028199130bda

### The design document of the research project

"Design the Search Before Searching for the Answer - rebuilding the search for the foundations of physics as a single project"
https://note.com/kiharanoriaki/n/n4c66305cb6cb

#TheoreticalPhysics #FineStructureConstant #Relativity #LightCone #Symmetry #GaugeTheory #StandardModel #ParticlePhysics #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
