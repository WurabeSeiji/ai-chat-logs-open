# Read Two Complex Numbers as One Quaternion, and Spacetime (1+3) Appears - A Map from a Set of Complex Numbers to Spacetime

In the previous article, "Do Symmetries Carry Dimensionless Numbers Rather Than Conserved Quantities?", I introduced a catalogue of symmetries and magic numbers.

While putting that catalogue in order, one big question was left open.

Without placing a background spacetime, how do we read out spacetime (1+3)(t, xyz) from nothing but a set of complex numbers {zₙ}?

And further, how do we map it to (1+3)(1+3)(t, xyz)(R, Q₁Q₂Q₃), which includes a spacetime that cannot be observed?

In the course of the study, I noticed that two complex numbers can be used to make a quantity called a quaternion q.

Looking a little further, this quaternion turned out to have a symmetry that splits it into 1+3 when it is read through multiplication.

As a result, the following map turned out to be possible.

{zₙ} (2N complex numbers) ↦ {qₙ} (N quaternions) ↦ {(t, xyz)ₙ} (N points of spacetime)

The only new thing I placed is a way of reading: read two complex numbers as one quaternion.

No background space, no metric, no derivative and no gauge field are placed.

In this article I introduce what is inside, with figures.

[Insert Figure 1: fig1_mapping_chain_en.png]

Figure caption: The chain of maps. The step that pairs 2N complex numbers into N quaternions is a full projection that loses nothing. The step that reads the quaternions by central projection is a partial projection, in which only part maps.

## 1. A quaternion is just a four-dimensional vector. Only multiplication differs

From two complex numbers z₁, z₂ we can make a quaternion q = z₁ + z₂ j.

In real components, q = q₀ + q₁i + q₂j + q₃k, with i² = j² = k² = ijk = −1.

Looking only at addition and size, a quaternion is just a four-dimensional vector.

The four components are on an equal footing, and there is no special direction.

The difference appears only when we multiply.

When the same unit is multiplied from both sides, one component (the real part) does not move, and the remaining three components turn while staying indistinguishable from each other.

A quaternion splits into 1 and 3 when it is read through multiplication.

What matters here is that which component is called time and which is called space is not a colour attached to the quaternion from the start.

Just as with the real and imaginary parts of a complex number, names are attached afterwards to the properties of the readout.

## 2. The quaternion circle: a three-dimensional sphere made of stacked spheres

Just as a complex number is written by its size and angle, a quaternion can be written

q = |q|(cos θ + n sin θ).

n is just a three-dimensional unit vector.

Multiplied as a quaternion, however, n² = −1, and it works like the imaginary unit i.

The imaginary unit, which was only ±i for complex numbers, fills a whole sphere for quaternions.

The quaternions of size 1 form the three-dimensional sphere S³, which I call the quaternion circle.

Once the angle θ is fixed, the real part is fixed at cos θ, and the remaining three components move on a sphere of radius sin θ.

The quaternion circle is these spheres stacked over the values of the real part.

Think of the circles of latitude on the Earth, raised by one dimension.

[Insert Figure 2: fig2_s3_latitude_slices_en.png]

Figure caption: Cutting the quaternion circle S³ at real part q₀ = cos θ gives a sphere of radius sin θ. θ = 0 is the north pole (a point), and θ = 90° is the equator (the largest sphere, of radius 1). The lower row shows the position of the cut on the meridian.

## 3. The form of a quaternion is already a central projection

This is what surprised me most this time.

A quaternion with positive real part can be written q = q₀(1 + v/q₀).

v is the imaginary part, with three components.

1 + v/q₀ is a point on the hyperplane q₀ = 1 tangent at the north pole, and its coordinate is the ratio u = v/q₀.

Projecting it from the origin onto the sphere gives the point cos θ + n sin θ on the quaternion circle.

In other words,

u = v/q₀ = n tan θ.

The polar form of a quaternion turned out to be the same thing as the central projection (gnomonic projection) from the hyperplane tangent at the north pole.

The name of this research series, "gnomonic mapping", had been inside the multiplication structure of quaternions from the start.

[Insert Figure 3: fig3_gnomonic_projection_en.png]

Figure caption: Central projection onto the tangent plane at the north pole. A circle of latitude at angle θ maps to a circle of radius tan θ on the tangent plane. The cross-section on the right also shows that the point at θ = 120° passes through the opposite side of the origin and maps to the same point as −q.

Several things can be read directly from this correspondence.

Multiplying the whole by any factor maps to the same point. Only ratios can be read, and no ruler is fixed.

The equator, where the real part is 0, goes to infinity.

q and −q map to the same point.

Great circles on the sphere map to straight lines on the tangent plane.

## 4. The light cone comes out of the true square

Square the quaternion as it is, without taking the conjugate.

q² = (q₀² − |v|²) + 2q₀v

The real part q₀² − |v|² places the square of one component and the squares of three components with opposite signs.

This is signature (1, 3), that is, the Minkowski form.

No metric was placed from outside. It comes out of the true square.

For q of size 1, q² = cos 2θ + n sin 2θ.

If θ is less than 45°, q is inside the light cone; at exactly 45°, on the light cone; above 45°, outside.

[Insert Figure 4: fig4_light_cone_true_square_en.png]

Figure caption: The light cone and the true square. At θ = 45°, q lies on the light cone, and its square q² comes to the equator, where the real part is 0. The green circles are where the quaternion circle meets the light cone.

## 5. Treat the four components equally, and the true square does not close

Let N quaternions close by zero closure (the sum is 0, and the sum of true squares is also 0).

For complex numbers, the condition treating the real and imaginary parts equally and the true-square closure were the same thing.

For quaternions they are not.

If the condition treating the four components completely equally (the 13 isotropic conditions) is satisfied, the sum of true squares always deviates from 0.

If the spread on the time side is 1, the spread on the space side is 3, and the difference remains.

The 16-cell {±1, ±i, ±j, ±k} is an example.

It is the four-dimensional version of the complex cross {±1, ±i}; it is isotropic, but the sum of its true squares is −4.

[Insert Figure 5: fig5_16cell_isotropic_vs_true_square_en.png]

Figure caption: The 16-cell. It spreads isotropically, but the spreads on the time side and the space side are 1 : 3, and the sum of true squares is −4, not 0.

What preserves the true-square closure is only three-dimensional rotation by the same unit from both sides, and real scaling of the whole.

The 1+3 split comes from this symmetry.

## 6. The most symmetric form: each one lies on the light cone

Then what is the most symmetric configuration satisfying the true-square closure?

At the highest symmetry, all q are equivalent.

Then the real part of the true square takes the same value for every q.

Since the sum is 0, the real part of each q is 0; that is, each one lies on the light cone.

The answer is

q = ±e^(π n/4) = ±(1 + n)/√2,

with n placed at the vertex directions of a regular polyhedron.

Then the square is q² = n, and the eighth power returns to 1 (q⁸ = 1).

The most symmetric is the configuration of 24 using the 12 directions of the regular icosahedron.

Seen through central projection, the 24 map to 12 points, and all 12 lie on the image of the light cone (the sphere of radius 1), forming the vertices of an icosahedron.

[Insert Figure 6: fig6_most_symmetric_icosahedral_en.png]

Figure caption: The most symmetric form q = ±e^(π n/4) (icosahedral, 24 quaternions). Left: the image under central projection, all on the sphere of radius 1. Right: the sum of true squares when the overall angle θ is shifted; it is 0 only at 45° and 135°.

## 7. The spatial directions stop at three

The number of directions that N quaternions can open as space is min(3, N−2).

For N = 3 there is only one plane, for N = 4 up to two directions, and for N ≥ 5 three directions, where it stops.

It stops because the imaginary part of a quaternion has only three components.

Even though the degrees of freedom of the state grow with N, the dimension of readable space does not.

The smallest number that opens three directions isotropically is N = 5.

N = 5 also appeared as a special place in my earlier numerical experiments on inflation.

## 8. This map is a partial projection

Looking at the chain of maps step by step:

The step from complex numbers to quaternions is a full projection that loses nothing.

The central-projection step is a partial projection in which only the three components of the ratio v/q₀ can be read.

The size, the ± sign, the equator and the global rotation escape the readout.

In the previous catalogue I wrote that, since spacetime is read stably, the map must be at least a partial projection.

This map is a concrete form that satisfies that condition.

## 9. Figures you can rotate

A four-dimensional figure can only be shown as cross-sections in still images.

I have prepared interactive figures (HTML) that can be rotated with sliders and dragging.

Open the links, and they run directly in your browser.

Slices of the quaternion circle, decomposition into real and imaginary parts, central projection, light cone:
https://wurabeseiji.github.io/ai-chat-logs-open/interactive/quaternion_circle_gnomonic.html

The 16-cell configuration and the closure checks:
https://wurabeseiji.github.io/ai-chat-logs-open/interactive/quaternion_16cell.html

The most symmetric form (icosahedral, in the same format as the quaternion-circle figure):
https://wurabeseiji.github.io/ai-chat-logs-open/interactive/quaternion_circle_icosa.html

The most symmetric form (icosahedral, in the same format as the 16-cell figure):
https://wurabeseiji.github.io/ai-chat-logs-open/interactive/quaternion_light_cone_icosa.html

## 10. Next tasks

There are three things I do not yet know.

The first is what is seen when the localization made by superposing equal-amplitude harmonics (all, odd, even) is mapped to (1+3). Localization could be made with complex numbers, but whether it reads as a point localized in three spatial directions has not yet been checked.

The second is whether, just as two complex numbers made a quaternion, two quaternions can make an octonion and map to (1+3)(1+3), including the unobservable side. For octonions the associative law of multiplication fails, and new problems appear.

The third is whether highly symmetric quaternions, each interacting with all others, produce an inflation-like expansion. We may see three-dimensional space opening from a plane and stopping at three directions.

To write that I do not know what I do not know is, for now, the answer.

For details, please see the PDF on Zenodo.

## Paper information

### A Map from a Set of Complex Numbers to Spacetime (1+3)

"A Map from a Set of Complex Numbers to Spacetime (1+3): Readout by Quaternions and Partial Projection"

Version: v1.1 (2026-10-18)

Concept DOI:
https://doi.org/10.5281/zenodo.23451494

Version DOI:
https://doi.org/10.5281/zenodo.23452319

English PDF (Zenodo):
https://zenodo.org/records/23452319/files/quaternion_spacetime_map_en_v1.pdf

Folder of the paper (Japanese and English manuscripts, figures and reproduction scripts, citation tables, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E5%9B%9B%E5%85%83%E6%95%B0%E3%81%AB%E3%82%88%E3%82%8B%E6%99%82%E7%A9%BA%E3%81%B8%E3%81%AE%E5%86%99%E5%83%8F_20261018

The Japanese full text and the four interactive figures (HTML) are included in the same Zenodo record. The Japanese version is the original; the English version is a translation.

### The Japanese version of this article

【日本語版の note の URL を記入】

### The previous article

"Do Symmetries Carry Dimensionless Numbers Rather Than Conserved Quantities? - A Catalogue of Symmetries and Magic Numbers for Reading Out Spacetime from a State Space"

https://note.com/kiharanoriaki/n/ne4dcf69af398

#TheoreticalPhysics #Quaternions #Spacetime #LightCone #CentralProjection #GnomonicProjection #ComplexNumbers #Symmetry #RegularPolyhedra #Octonions #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
