# The Conjugate Complex Number Is Really Just a Sum of Squares —— rereading a basic tool of theoretical physics, simply

In physics and engineering, one expression always shows up when we write down waves and oscillations: the complex number. In particular, the operation of multiplying a complex number by its conjugate appears everywhere, from quantum mechanics to signal processing.

This article introduces just one point, as gently as possible: that this operation on conjugate complex numbers is in fact equal to a very simple sum of squares. This is not a story about a new law of physics. It is a short observational note that merely rearranges already-known mathematics into a clearer view.

Japanese and English versions are published on Zenodo.

・Concept DOI (always the latest): https://doi.org/10.5281/zenodo.21126212
・This version: https://doi.org/10.5281/zenodo.21126213

## A complex number and its conjugate

A complex number bundles two numbers, a real part and an imaginary part, into one.

　Z = x + iy

Here i is the imaginary unit, with the property i² = −1. x is the real part, y the imaginary part.

The conjugate of Z flips only the sign of the imaginary part.

　conjugate = x − iy

Multiplying a complex number by its conjugate gives:

　(x + iy)(x − iy) = x² + y²

The imaginary cross-terms cancel neatly, and what remains is only x² + y², the sum of the squares of two real numbers. The imaginary unit i has vanished.

This is exactly the squared distance of the point (x, y) from the origin, i.e. the Pythagorean theorem itself. Writing the radius as ρ,

　x² + y² = ρ²

If we write a complex number by its size ρ and its direction (phase) θ, then Z = ρ(cosθ + i sinθ); but multiplying by the conjugate makes the direction θ vanish, leaving only the square of the size ρ.

In other words, multiplying by the conjugate is the operation that drops the directional information and keeps only the squared magnitude.

## Extending to many axes, on an equal footing

What is interesting is that this extends directly to many dimensions.

Take many real coordinates and bundle them two at a time into complex numbers.

　Z₁ = x₁ + i x₂ , Z₂ = x₃ + i x₄ , …

Multiplying each by its conjugate and summing gives

　x₁² + x₂² + x₃² + x₄² + … = R²

So however many conjugate complex numbers we bundle, in the end it is just the sum of the squares of all the coordinates = a single squared radius.

The important thing here is that all axes are on an equal footing. No axis is treated specially. An axis that does not appear in the complex plane we are currently looking at (call it an invisible axis) is simply added with the same plus sign, on the same footing. The radius R is always the single "square root of the sum of all the squares," and its meaning never becomes twofold.

This is the main point of the article. The conjugate complex number, a familiar expression in theoretical physics, is in fact just another way of writing an equal-footing sum of squares.

## One more step (but this is not physics)

Finally, one more step. Let me state clearly in advance that this is a game of algebra, not a claim about physics.

Above, we bundled real coordinates. Now let us deliberately take one coordinate to be an imaginary value. For example, set some coordinate to i z. Then its square is

　(i z)² = i² z² = − z²

and a square with a minus sign appears. The minus did not come from us deciding a sign by hand; it came out naturally from the property i² = −1 of imaginary numbers.

Then, setting the sum of all the squares to zero, a form appears in which only one side has moved to the minus side:

　x² + y² + … = z²

The point here is that this form comes from an operation different from the conjugate multiplication (which extracts a magnitude and is always plus). For real coordinates the two operations agree, but only for an imaginary coordinate do they diverge. The minus arises only there.

To repeat: this is an observation about rewriting algebra. I have no intention of saying that the z here is time, or that this formula derives relativity or the light cone. The coordinate symbols are geometric labels only; I make no claim that they are the same as the physical dimensions of space or time. To proceed there would require entirely separate physical assumptions. This article stops just short of that.

## Summary

・Multiplying a complex number by its conjugate drops the direction (phase) and keeps only the squared magnitude. The result is x² + y², a plain Pythagorean sum of squares.
・This extends to many dimensions with all axes kept on an equal footing. Bundling conjugate complex numbers = a single squared radius.
・Taking one coordinate to be imaginary produces a minus square because i² = −1. But this is a consequence of algebra, not a claim of physics (time, relativity, metric).

Rewording the conjugate complex number, so familiar in theoretical physics, into the simplest language of a sum of squares turns out to be surprisingly clarifying —— that was all this small note was about.

――――
(This is the fourth note in the "square-quantity readout" series. For the formal version with equations, see the Zenodo record above and the Zenn article https://zenn.dev/noriaki_kihara/articles/complex-norm-square-readout .)

#ComplexNumbers #Mathematics #Physics #TheoreticalPhysics #QuantumMechanics #Pythagorean #LinearAlgebra #Relativity #Spacetime #WickRotation #IndependentResearch #Zenodo #Science
