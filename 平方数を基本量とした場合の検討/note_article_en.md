# Reading by Squares: When a Curved Sum of Squares Becomes a Straight Line

## Purpose of This Article

I have released a new observational and organizing paper.

This paper does not claim to prove physical laws. It does not claim to replace established physics.

The starting point is very elementary.

Instead of reading a positive quantity x directly, read the square quantity X = x² as the basic quantity.

Then the sum-of-squares constraint, which looks curved,

　　Σ xᵢ² = E

is read on the square-quantity side as

　　Σ Xᵢ = E.

That is, it becomes a linear simplex constraint.

This readout itself is close to known structures, such as the square-root representation connecting probability simplexes and spherical regions. Therefore, I do not call that correspondence a new discovery.

The point of this paper is different. When this square-quantity readout is placed side by side with the curvature-distortion result proved in my earlier Paper 0, one can separate, in very simple geometric terms,

- which equations should not receive a curvature correction,
- and from which equations curvature-correction candidates first become worth considering.

The paper is available here.

- Concept DOI: https://doi.org/10.5281/zenodo.20785539
- Version DOI: https://doi.org/10.5281/zenodo.20785540
- GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/平方数を基本量とした場合の検討

---

## What Paper 0 Proved

The present discussion rests on the self-referenced Paper 0.

Paper 0 rigorously computed how edge length, angle, area, and volume are distorted when a geodesic unit cell of side length 1 is placed in a positively curved constant-curvature space.

The important conclusions are these.

- The geodesic edge length remains exactly 1 by construction.
- The angle deviates from 90 degrees depending on the curvature radius R.
- The area is distorted depending on the curvature radius R.
- The volume distortion depends on the dimension d, and its leading coefficient is d(d−1)/12.

In short, curvature distortion does not first appear in length itself.

It first appears in the area spanned by two directions.

This is the point that connects Paper 0 with the square-quantity readout in the present paper.

【Paper 0 Figure B: Unfolding of a geodesic band】

Image file: 波長空間と周波数空間の双対幾何/paper0_figB_band_unfold.png

In Paper 0, edges are preserved as geodesic lengths. On the other hand, once a surface is spanned, the effect of curvature appears.

For the angle, the following expression is obtained as a function of curvature radius R:

　　θ(R) = arccos(−tan²(1/(2R))).

Using this angle, the area of a geodesic square of side length 1 is

　　A(R) = R² { 4θ(R) − 2π }.

Since the area of a flat unit square is 1, the area-correction coefficient can be defined as

　　k_s(R) = R² [ 4 arccos(−tan²(1/(2R))) − 2π ].

In the flat limit,

　　k_s(R) → 1,

and the leading expansion is

　　k_s(R) = 1 + 1/(6R²) + O(R⁻⁴).

【Paper 0 Figure C: Angle excess】

Image file: 波長空間と周波数空間の双対幾何/paper0_figC_angle_excess.png

For d-dimensional volume, in the small-curvature limit,

　　V_d(R) = 1 + c_d/R² + O(R⁻⁴),

　　c_d = d(d−1)/12 = C(d,2)/6.

This form is important.

The coefficient d(d−1)/12 is proportional to the number C(d,2) of independent pairs of directions. In other words, the d-dimensional volume distortion can be read as the sum of area distortions over the two-dimensional coordinate planes.

This gives the test rule used in the present paper.

**The first entrance to curvature correction is not length, but area.**

---

## What Is Square-Quantity Readout?

In the present paper, on the positive domain, set

　　Xᵢ = xᵢ².

For example, the positive circular arc in two dimensions,

　　x² + y² = 1,

is read on the square-quantity side as

　　X + Y = 1,

which is a line segment.

【Figure 1: Two-dimensional square map】

Image file: figures/fig01_2d_square_map.png

In three dimensions, the positive spherical octant

　　x² + y² + z² = 1

is read on the square-quantity side as

　　X + Y + Z = 1,

which is a triangular simplex.

【Figure 2: Three-dimensional square map】

Image file: figures/fig02_3d_square_map.png

The important caution is that this is not an isometric flattening.

The square map does not preserve distances.

Taking differentials gives

　　dXᵢ = 2xᵢ dxᵢ,

so the way lengths are measured changes.

Therefore, the intrinsic metric of a sphere does not simply become the Euclidean metric of a plane.

The fact used in this paper is much more limited.

A sum-of-squares constraint becomes a linear equation on the square-quantity side.

That is the only point being used.

---

## Placing Simple Equations on the Square-Quantity Side

Place the following equation on the square-quantity side:

　　X = A T.

Now read back by setting

　　X = x²,

　　T = t².

Then

　　x² = A t².

Taking the positive square root gives

　　x = √A t.

This has the same form as uniform motion,

　　x = vt.

Next, place

　　X = 1/2 B T²

on the square-quantity side. Then

　　x² = 1/2 B t⁴,

and therefore

　　x = √(B/2) t².

This has the same form as uniformly accelerated motion,

　　x = 1/2 a t².

【Figure 3: Equations on the square-quantity side and motion-form readouts】

Image file: figures/fig03_motion_readouts.png

Here again, caution is needed.

This is not a physical derivation of velocity or acceleration.

It only says that when a simple equation placed on the square-quantity side is read on the positive square-root side, the resulting form is isomorphic to a familiar form in classical mechanics.

That is all.

---

## Why No Curvature Correction Appears in One-Dimensional Equations

This is where the result of Paper 0 matters.

In Paper 0, the first entrance to curvature distortion was area.

A relation closed by length alone does not contain an area cell spanned by two directions.

Therefore, even if the square-root side gives

　　x = √A t

or

　　x = √(B/2) t²,

there is no reason to immediately multiply it by k_s(R).

This judgment is not determined by whether acceleration appears.

It is determined by whether the equation is closed by one-dimensional quantities, or whether it contains an area quantity spanned by two directions.

A one-dimensional two-body collision makes this point clearer.

In an ordinary one-dimensional perfectly elastic collision, if masses m₁, m₂ and velocities u₁, u₂ are externally assigned, the usual collision formulas follow from conservation of momentum and conservation of energy.

The energy equation contains u².

However, this u² is the square of the same one-dimensional velocity component. It is not an area spanned by two independent directions.

Therefore, there is still no place to insert k_s(R).

This distinction is important.

The entrance to correction is not mass, acceleration, or the fact that there are two bodies.

The entrance to correction is the appearance of an area quantity.

---

## The Centrifugal-Type Model as the First Natural Entrance

Then where does area appear?

As the smallest natural example, consider a centrifugal-type equation.

Classically, it is written as

　　a_c = v²/r.

Here the tangential speed v and the radial length r appear at the same time.

Geometrically, circular motion involves two independent directions: the tangential direction and the radial direction.

Here, for the first time, one has an entrance for considering a two-dimensional area cell rather than length alone.

Therefore, taking the geodesic-square area coefficient from Paper 0,

　　k_s(R) = R² [ 4 arccos(−tan²(1/(2R))) − 2π ],

one minimal test model places it into the centrifugal-type equation as

　　a_c^(R) = k_s(R) v²/r.

【Figure 5: Area-correction coefficient k_s(R)】

Image file: figures/fig05_area_coefficient_ks.png

Here k_s(R) is not a freely adjustable coefficient.

Once R is fixed, it is uniquely determined by the geometric formula of Paper 0.

But again, this must be separated carefully.

The fact that k_s(R) is uniquely determined is not the same as the claim that it must multiply the centrifugal-type equation in precisely this direction.

Depending on whether one reads the quantity as a measure, density, area, or inverse area, k_s(R), 1/k_s(R), or another function may appear.

Therefore, the claim of the paper stops here.

It does not prove centrifugal force.

It only observes that the centrifugal-type equation contains two directions and, unlike one-dimensional equations, has a natural entrance for considering an area correction.

---

## The Most Important Boundary in This Paper

If misread, this paper could easily be taken to say the following.

The square map makes curvature disappear.

k_s(R) comes from the non-isometry of the square map.

Classical mechanics is derived from square-quantity readout.

None of these are claims of the paper.

The actual claim is more modest and more geometric.

- A sum-of-squares constraint can be read as a linear simplex constraint on the square-quantity side.
- Under that readout, a linear equation and a quadratic equation appear in forms isomorphic to uniform motion and uniformly accelerated motion.
- However, an equation closed by one-dimensional quantities alone has no place for the area correction from Paper 0.
- Only when two directions appear and an area cell is involved does k_s(R) first become a candidate for consideration.

This is the organizing principle of the paper.

The square map itself is close to known structures.

The curvature-distortion formula of Paper 0 is also established separately as geometry.

The interest of the present paper is not to mix the two and derive one from the other.

Rather, the interest is to not confuse them.

Separate equations of length alone from equations involving area.

That separation determines where curvature-correction candidates should be considered.

---

## What This Paper Does Not Claim

Finally, let me make clear what this paper does not claim.

- It does not claim to prove physical laws.
- It does not claim to replace classical mechanics, relativity, or quantum theory.
- It does not identify Xᵢ with probability, energy, mass, time, or any particular physical quantity.
- It does not claim that the square map is an isometry.
- It does not claim that a curved space itself disappears into a flat space.
- It does not claim that the way k_s(R) is inserted into the centrifugal-type model is uniquely derived as a law of mechanics.
- It does not claim that the non-isometry of the square map is the cause of k_s(R).
- It does not identify the metric distortion of the square map with the geodesic-square area excess of Paper 0.

This paper is not written to enlarge a discovery.

It is written to draw lines where misreadings easily enter.

Still, I think this organization has value.

The square-quantity readout simplifies equations more than one might expect, but curvature correction does not enter everywhere.

Clarifying where correction does not enter is just as important as finding where it may enter.

---

## Related Links

- Present paper Concept DOI: https://doi.org/10.5281/zenodo.20785539
- Present paper Version DOI: https://doi.org/10.5281/zenodo.20785540
- Paper 0 Concept DOI: https://doi.org/10.5281/zenodo.20680269
- Paper 0 Version DOI: https://doi.org/10.5281/zenodo.20684135
- Zenn technical article: https://zenn.dev/noriaki_kihara/articles/square-quantity-readout-simplexification
- GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/平方数を基本量とした場合の検討

#mathematics #geometry #differential_geometry #observational_paper #square_map #curvature #area_correction #FisherRao_geometry #Zenodo
