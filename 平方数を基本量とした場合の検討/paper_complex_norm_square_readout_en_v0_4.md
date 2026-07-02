# On the Connection Between the Conjugate Complex Norm and Square-Quantity Readout

## An observational extension from phase representation to a multi-dimensional sum-of-squares representation including an invisible axis (pure algebra / non-physical)

**Author**: Noriaki Kihara
**ORCID**: 0009-0004-6753-4020
**Version**: v0.4 (English draft)
**Date**: July 2026
**DOI**: [10.5281/zenodo.21126213](https://doi.org/10.5281/zenodo.21126213) (this version) / Concept DOI: [10.5281/zenodo.21126212](https://doi.org/10.5281/zenodo.21126212)
**Zenodo**: https://zenodo.org/records/21126213
**Position**: The fourth observational note in the "square-quantity readout / odd-harmonic phase structure / common relative phase" series. It treats not central-projection geometry itself, but the algebra of conjugate complex numbers and its geometric interpretation. It does not connect to physics.

---

## 0. Premises and stance (how to read this note)

This note starts **only from conjugate complex numbers and their geometric interpretation**. It takes the stance of **not connecting to physics** (spacetime, special/general relativity, quantum theory, metric signature). Accordingly, we declare the following explicitly.

1. The coordinate symbols $x,\ y,\ z_1,\ z_2,\ z_3,\dots$ that appear below are all **geometric coordinate labels**. **They are not identified with physical spatial or temporal dimensions.**
2. In particular, this note **does not use** a symbol $t$ for time. Even when a coordinate taking imaginary values is introduced, it is not "time" but merely an **imaginary-valued coordinate label** $z_h$.
3. Even if a sign (a minus) of the sum-of-squares-difference type appears in a formula, it **does not assert a metric signature**. The minus is nothing more than an **algebraic consequence** of $i^2=-1$, arising when a coordinate is placed as imaginary.
4. All coordinate axes are treated **on an equal footing**. Here, "equal footing" means **not assigning any particular physical role or metric signature from the outset**. Splitting, in a later section, into real-valued and imaginary-valued components as an algebraic representation (§5) does not violate this equal footing.
5. This note asserts no new physical law and no new mathematical theorem. It is an **observational note** that re-reads the known conjugate complex norm and sum of squares by placing them within the preceding square-quantity-readout series.

This declaration is the design foundation of this note, made to remove the error of v0.1 (assigning a minus sign to a particular axis by hand, and using the radius symbol in two senses).

---

## Abstract

This note is the fourth in the series following linear simplexification by square-quantity readout, the isolated-peak wave by a half-wavelength odd-harmonic sum, and the contrast law in two-copy common-relative-phase superposition.

In the preceding square-quantity-readout paper, reading $X_i=x_i^2$ in the positive region showed that the sum-of-squares constraint $\sum_i x_i^2=E$ can be read as the linear simplex constraint $\sum_i X_i=E$ in square-quantity space. Here, we connect to this square-quantity readout the most basic norm of conjugate complex numbers,

$$
Z\bar Z=x^2+y^2 .
$$

For the complex polar form $Z=\rho e^{i\theta}=\rho(\cos\theta+i\sin\theta)$, putting $x=\rho\cos\theta,\ y=\rho\sin\theta$ gives $Z=x+iy$ and $Z\bar Z=x^2+y^2=\rho^2$. This is a standard fact, and there is no novelty in the formula itself. What this note observes is that this known norm connects naturally to the square-quantity readout as an **operation that bundles the square quantities of two orthogonal real components into a single real norm**.

Next, we generalize this to an **equal-footing multi-dimensional extension**. Bundling real components two at a time into complex numbers gives

$$
\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2 ,
$$

and a coordinate axis not displayed in the complex plane (what we call an invisible axis) simply joins as an **equal-footing axis with the same sign $+$**, $z_1,z_2,\dots$. There is no special axis here.

Finally, to the front-half conjugate norm $\sum_j Z_j\bar Z_j$ (Hermitian, positive-definite) we add the square $(i z_k)^2$ of an imaginary-valued coordinate, and observe the **zero-sum-of-squares form**

$$
\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0 .
$$

Here the conjugate product $X\bar X$ and the non-conjugate square $X^2$ **agree for real coordinates but diverge only for imaginary coordinates** ($X_h\bar X_h=z_h^2$ versus $X_h^2=-z_h^2$), and this divergence is the sole origin of the minus sign. Expanding gives $\sum_j Z_j\bar Z_j=\sum_k z_k^2$, i.e.

$$
\sum_n x_n^2=\sum_k z_k^2 .
$$

The minus sign appearing here is an algebraic consequence of $i^2=-1$, not a hand-placed metric signature. The conjugate-norm type $Z\bar Z$ and the non-conjugate square $X_h^2$ of an imaginary coordinate are not the same operation, and this note keeps the distinction. All of this is a geometric observation on the **algebraic arrangement** of conjugate complex numbers and sums of squares; it makes no claim about physical spacetime or relativity.

---

## 1. Position within the series

This note is not a paper on central-projection geometry itself. The central-projection series takes the curvature radius, geodesic cells, and angle/area/volume corrections as its subjects, and deliberately kept its distance from conjugate-complex representations. Inserting the content of this note directly into the central-projection series would shift the center of gravity of the argument.

This note belongs, rather, to the following series.

1. Linear simplexification by square-quantity readout, and a dimension-by-dimension organization of curvature-correction candidates.
2. The isolated-peak wave of a constant-amplitude odd-harmonic sum on a half-wavelength phase interval.
3. Waveform invariance and the contrast law in two-copy common-relative-phase superposition.
4. This note, i.e. the connection between the conjugate complex norm and square-quantity readout, and its equal-footing multi-dimensional extension.

Following these, this note treats $Z\bar Z=x^2+y^2$ as the minimal contact point linking phase representation, conjugate complex numbers, and square-quantity readout, and generalizes it to an equal-footing multi-dimensional sum of squares.

---

## 2. Complex polar form and sum of squares

### 2.1 Complex polar form

Write a complex number $Z$ as $Z(\rho,\theta)=\rho(\cos\theta+i\sin\theta)$, where $\rho\ge0$ is the radius in the complex plane and $\theta$ the phase angle. Putting $x=\rho\cos\theta,\ y=\rho\sin\theta$ gives $Z=x+iy$. The conjugate is $\bar Z=x-iy$, so

$$
Z\bar Z=(x+iy)(x-iy)=x^2+y^2 .
$$

From the polar form as well, $Z\bar Z=\rho e^{i\theta}\rho e^{-i\theta}=\rho^2$. Hence

$$
Z\bar Z=x^2+y^2=\rho^2 .
$$

### 2.2 Interpretation as square-quantity readout

In the preceding square-quantity-readout paper, reading $X_i=x_i^2$ in the positive region organized the sum-of-squares constraint $\sum_i x_i^2=E$ as the linear constraint $\sum_i X_i=E$. This note's $Z\bar Z=x^2+y^2$ is the **two-component version** of this square-quantity readout, in which two square quantities $x^2,\ y^2$ are bundled into a single real norm as the conjugate product $Z\bar Z$.

### 2.3 Two roles of $i$ (a distinction made in this note)

This note clearly distinguishes two different roles of the imaginary unit $i$. Neglecting this distinction produces errors in the multi-dimensional extension below.

- **Role (1): bundling.** When $Z=x+iy$ with $x,y$ both real, $Z\bar Z=x^2+y^2>0$. Here $i$ is merely a notational device for gathering two real coordinates into one complex coordinate. The norm is always non-negative.
- **Role (2): imaginary-valued coordinate label.** When the coordinate label itself is specified as an imaginary value, we set $X_h=i\,z_h$ ($z_h$ a real label). Then $X_h^2=(i z_h)^2=-z_h^2$, producing a **negative squared contribution**.

Role (1) gives $x^2+y^2$ (positive-definite); role (2) gives $-z_h^2$ (a negative contribution). The two are distinct and must not be treated as continuous with each other. More precisely, role (1) appears in the **conjugate product** $Z\bar Z$ (Hermitian, positive-definite), while role (2) appears in the **non-conjugate square** $X_h^2$ (a complex bilinear form). This difference in **presence or absence of conjugation** becomes decisive in §5.

### 2.4 The location of novelty

$Z\bar Z=x^2+y^2$ is a standard fact. This note does not present this formula as a new theorem. Its aim is to place the known complex norm within the square-quantity-readout series and make explicit the contact point

$$
\text{square-quantity readout}\quad\leftrightarrow\quad\text{conjugate complex norm},
$$

and to generalize it to an equal-footing multi-dimensional sum of squares.

---

## 3. Reading as a phase representation

In the complex polar form, $\theta$ is an angle in the plane and can also be read as the phase angle of a wave. That is,

$$
\operatorname{Re}Z=\rho\cos\theta,\qquad \operatorname{Im}Z=\rho\sin\theta,
$$

and $Z=\rho e^{i\theta}$ is a representation that expands a state with phase $\theta$ into two orthogonal components in the complex plane.

What we treat at this stage is the map $(\rho,\theta)\longmapsto(x,y)$, from a polar form containing a phase angle to two orthogonal components. No coordinate of time or propagation direction is introduced here (as declared in §0, the time symbol $t$ is not used).

---

## 4. Extension to an equal-footing multi-dimensional sum of squares

### 4.1 Bundling by multiple complex numbers

Let the real coordinates be $x_1,x_2,\dots,x_{2m}$ and bundle them two at a time into complex numbers,

$$
Z_j=x_{2j-1}+i\,x_{2j},\qquad Z_j\bar Z_j=x_{2j-1}^2+x_{2j}^2 .
$$

Their sum is

$$
\sum_{j=1}^{m} Z_j\bar Z_j=\sum_{n=1}^{2m} x_n^2 .
$$

For an odd count, add the unbundled real coordinate $x_{2m+1}$ on its own. In either case,

$$
\sum_n x_n^2=R^2
$$

is a positive-definite sum of squares containing all coordinates on an equal footing.

### 4.2 The invisible axis is not a special axis

When a coordinate does not appear in the complex plane currently displayed (say the $Z_1$ plane), we call it an **invisible axis**. But this means only that it is not carried on that plane in the display; its **sign and role are the same as the other axes**. Adding invisible axes $z_1,z_2,\dots$ merely amounts to

$$
\sum_n x_n^2 + z_1^2 + z_2^2 + \cdots = R^2 ,
$$

adding terms of **the same sign $+$ to the left-hand side so that $R$ grows**. The radius $R$ remains, throughout, the single "square root of the sum of squares of all coordinates," and is never defined in two ways.

This is the point that removes the v0.1 error. Taking a particular axis alone and placing it on the other side, or giving it a minus sign, breaks $\sum_n x_n^2=R^2$ and makes $R$ carry two meanings. As long as the axes are treated equally, that breakdown does not occur.

Note that "invisible" here is limited to the **display-level meaning** of not appearing in the chosen complex-plane display. It does not mean a physically unobservable axis, a hidden physical dimension, or an unknown degree of freedom.

### 4.3 Representation by a tilt angle (spherical coordinates)

For three coordinates, using the radius $R$, a tilt angle $\phi$, and the phase angle $\theta$,

$$
x=R\cos\phi\cos\theta,\qquad y=R\cos\phi\sin\theta,\qquad z_1=R\sin\phi ,
$$

we immediately obtain

$$
x^2+y^2+z_1^2
=R^2\cos^2\phi(\cos^2\theta+\sin^2\theta)+R^2\sin^2\phi
=R^2 .
$$

Here $\phi$ is a tilt angle leaving the complex plane, not an additional in-plane rotation phase. The two must be distinguished: an additional in-plane phase $\alpha$ merely gives $x=\rho\cos(\theta+\alpha),\ y=\rho\sin(\theta+\alpha)$, generates no $z_1$ component, and leaves $x^2+y^2=\rho^2$.

---

## 5. Conjugate-norm-type square quantity and the zero-sum-of-squares by an imaginary-valued coordinate

### 5.1 The difference between the conjugate product and the square —— agreement for real coordinates, divergence for imaginary ones

The previous $\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2$ is a **positive-definite norm containing conjugation** (the Hermitian form $\sum_n X_n\bar X_n=\sum_n|X_n|^2$), which is zero only when all components are zero.

Here we clarify the relation between the conjugate product $X\bar X$ and the non-conjugate square $X^2$. **For real-valued coordinates the two agree**: if $x$ is real, $x\bar x=x^2$. **They diverge only for imaginary-valued coordinates**: for $X_h=i\,z_h$ ($z_h$ real),

$$
X_h\bar X_h=(i z_h)(\overline{i z_h})=(i z_h)(-i z_h)=z_h^2\ (>0),
\qquad
X_h^2=(i z_h)^2=-z_h^2\ (<0).
$$

That is, **a negative squared contribution arises "only" when one takes the non-conjugate square of an imaginary-valued coordinate**. This is the sole origin of the minus signs appearing hereafter, and is not a hand-placed metric signature. Accordingly, this note **does not confuse** the conjugate-norm type $Z\bar Z$ with the non-conjugate square $X_h^2$ of an imaginary coordinate.

### 5.2 The zero-sum-of-squares form (conjugate-norm backbone + imaginary-coordinate square)

To the positive-definite square quantity $\sum_j Z_j\bar Z_j$ obtained from the conjugate norm, we add the square $(i z_k)^2$ of imaginary-valued coordinates, and observe the **zero-sum-of-squares form**

$$
\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0,
\qquad Z_j=x_{2j-1}+i x_{2j},\ \ z_k\in\mathbb{R}.
$$

Here the visible part is the conjugate norm (Hermitian, positive-definite) and the invisible part is the square of imaginary-valued coordinates (non-conjugate); the two are different operations. Expanding the left-hand side, $(i z_k)^2=-z_k^2$ gives

$$
\sum_j Z_j\bar Z_j-\sum_k z_k^2=0
\quad\Longrightarrow\quad
\sum_j Z_j\bar Z_j=\sum_k z_k^2 .
$$

Since further $\sum_j Z_j\bar Z_j=\sum_n x_n^2$,

$$
\sum_n x_n^2=\sum_k z_k^2 .
$$

Note that, since the visible coordinates are real, by §5.1 we have $x_n\bar x_n=x_n^2$. Hence this form coincides numerically with the non-conjugate quadratic form written as $\sum_n X_n^2=0$ for all coordinates at once. However, this note takes as primary the above form, written by **separating** the visible part as a conjugate norm and the invisible part as an imaginary-coordinate square, in order not to confuse the $Z\bar Z$ type (conjugate) with the $X^2$ type (non-conjugate).

### 5.3 A non-trivial solution requires a non-real component

If all coordinates are real, $\sum_j Z_j\bar Z_j=\sum_n x_n^2=0$ has only the all-zero solution (a positive-definite norm). Therefore, for the zero-sum-of-squares to be **non-trivial, at least one non-real component is required**. As the minimal organization, this note treats the non-real component only as a purely imaginary-valued coordinate $i z_k$ (in general the form $a_n+i b_n$ is also possible). The non-real component is not an additional assumption but appears as a necessary condition for non-triviality.

### 5.4 The origin of the minus sign, and non-identification

The minus sign appearing here is, as in §5.1, an algebraic consequence of $i^2=-1$ from the square of an imaginary-valued coordinate, not a metric signature assigned to a particular axis by hand. Moreover, this note does not identify $z_k$ with a time coordinate, a curvature radius, or a projection radius (§0, §6).

### 5.5 Contrast with the v0.1 error

v0.1, to the positive-definite three-dimensional sum $x^2+y^2+z^2=r^2$, added by hand a square $t^2$ of a separate symbol to the right-hand side, writing $x^2+y^2+z^2=r^2+t^2$. But this was not an equal-footing dimensional extension (adding same-sign terms to the left-hand side); it mixed **two independent operations**: (i) introducing a new coordinate, and (ii) placing it with a minus sign. As a result $\sum_n x_n^2=R^2$ broke down and $r$ carried two meanings.

In the present form, the minus sign is **derived** from the imaginary-coordinate square ($i^2$), the visible-part conjugate norm $Z\bar Z$ (positive-definite) and the invisible-part imaginary-coordinate square $X_h^2$ are clearly separated, and $R$ is not used in two ways, so this breakdown does not occur.

---

## 6. Readings this note does not bring about

We restate the §0 declaration in terms of concrete formulas.

- This note does not introduce a time coordinate. The $z_k$ on the right-hand side of $\sum_a x_a^2=\sum_k z_k^2$ is an imaginary-valued coordinate label, not time.
- This note does not assert a metric signature $(+,+,\dots,-)$. The negative squared contribution is algebra arising from $i^2=-1$, not a spacetime metric.
- This note derives none of the light cone, Minkowski space, special relativity, or general relativity. Even if $\sum_a x_a^2=\sum_k z_k^2$ formally coincides with those quadratic forms, this note makes no such identification.
- This note does not identify $z_k$ or $R$ with a curvature radius, projection radius, or internal radius.
- This note does not identify the coordinate labels $x,y,z_1,z_2,\dots$ with physical dimensions.
- The "imaginary-valued coordinate" and "invisible axis" in this note are algebraic and display-level specifications of coordinate labels; they do not mean unobservable physical degrees of freedom or hidden physical dimensions.

These are boundary conditions that prevent the algebraic observation from leaping to physical claims.

---

## 7. Relation to the preceding three notes

### 7.1 Relation to the square-quantity-readout paper

Whereas the square-quantity-readout paper treated the individual square-quantity readout $(x_i)\longmapsto(x_i^2)$, this note treats

$$
(x,y)\longmapsto Z=x+iy\longmapsto Z\bar Z=x^2+y^2 ,
$$

a two-component sum-of-squares readout by the conjugate complex norm, and its equal-footing multi-dimensional sum $\sum_j Z_j\bar Z_j=\sum_n x_n^2$.

### 7.2 Relation to the odd-harmonic isolated-peak paper

Whereas the second note treated "what localized shapes can be made on a phase interval," this note treats "how a complex representation carrying a phase angle connects to square-quantity readout." It targets not a concrete odd-harmonic sum $S_N(\varphi)$ but the phase representation $Z=\rho e^{i\theta}$ itself.

### 7.3 Relation to the relative-phase contrast paper

In the third note, giving two copies a common relative phase produced $\psi_\alpha=2\cos\alpha\,S_N$ and $I_\alpha=4\cos^2\alpha\,I_N$, so that the relative phase appeared as a squared-intensity coefficient without changing the waveform. In this note, the conjugate product $Z\bar Z$ makes the phase vanish, leaving the norm square $\rho^2$. What the two share is that a real quantity appears **when a phase-carrying structure is read as a square quantity**, not the phase itself.

---

## 8. Dangerous readings and safe readings

### 8.1 Dangerous readings (this note does not assert)

- Derived space or spacetime from complex numbers.
- Phase space and real space are identical.
- Derived special/general relativity from a complex norm.
- The negative squared contribution is a spacetime metric signature.
- The coordinate label $z_k$ is time.
- $R$ or $z_k$ is a curvature radius.
- Coordinate labels are physical dimensions.

### 8.2 Safe readings (what this note can say)

1. The complex polar form is a standard map sending a polar form containing a phase angle to two orthogonal components.
2. Its conjugate product is the sum of squares of two orthogonal real components, readable as a form of square-quantity readout.
3. Bundling real components two at a time generalizes to the equal-footing multi-dimensional sum of squares $\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2$.
4. An invisible axis not displayed in the complex plane simply joins as an equal-footing axis with the same sign, and the radius $R$ is unique.
5. The conjugate product $X\bar X$ and the non-conjugate square $X^2$ agree for real coordinates but diverge only for imaginary ones ($X_h\bar X_h=z_h^2$ versus $X_h^2=-z_h^2$). The minus sign arises only from this divergence.
6. Expanding the zero-sum-of-squares form $\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0$—the conjugate-norm backbone plus the imaginary-coordinate square—yields the algebraic relation $\sum_n x_n^2=\sum_k z_k^2$, that the sum of squares of the real coordinates equals the sum of squares of the imaginary coordinate labels. Non-triviality requires at least one imaginary-valued coordinate.
7. All of this is a geometric observation on the algebraic arrangement of conjugate complex numbers and sums of squares, and makes no physical claim.

---

## 9. Conclusion

In this note, as the fourth in the series following square-quantity readout, odd-harmonic phase structure, and relative-phase contrast, we observed the connection between the conjugate complex norm and square-quantity readout, and generalized it to an equal-footing multi-dimensional sum of squares.

For the complex polar form $Z=\rho e^{i\theta}$, putting $x=\rho\cos\theta,\ y=\rho\sin\theta$ gives $Z\bar Z=x^2+y^2=\rho^2$ (a standard fact). The center of this note is to place this known norm within the square-quantity-readout series and, further, to bundle real components two at a time and generalize to the equal-footing multi-dimensional sum of squares

$$
\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2 .
$$

There is no special axis here; an invisible axis also joins as an equal-footing axis with the same sign, and $R$ is unique.

Furthermore, adding the square of an imaginary-valued coordinate to this conjugate-norm backbone, we observed the zero-sum-of-squares form

$$
\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0 .
$$

Since the conjugate product $X\bar X$ and the non-conjugate square $X^2$ agree for real coordinates and diverge only for imaginary ones ($X_h\bar X_h=z_h^2$ vs $X_h^2=-z_h^2$), expanding gives, through $i^2=-1$ for the purely imaginary-valued coordinate, the algebraic relation

$$
\sum_{n} x_n^2=\sum_{k} z_k^2 .
$$

The negative squared contribution appearing here is a consequence of $i^2$, not a hand-placed metric signature, and no double definition of $R$ occurs. The conjugate-norm type $Z\bar Z$ and the non-conjugate square $X_h^2$ of an imaginary coordinate are not the same operation, and this note keeps the distinction.

The above records how the conjugate positive-definite norm $\sum_j Z_j\bar Z_j$ and the non-conjugate square term $(i z_k)^2$ of an imaginary-valued coordinate can be connected within a single zero-sum-of-squares form, together with consistency with square-quantity readout and the equal-footing multi-dimensional sum of squares—an observation limited to elementary algebra and geometric interpretation. This note does not connect to physics (spacetime, relativity, quantum theory, metric signature), and does not identify coordinate labels with physical dimensions.

---

## References

### Self-references

[1] N. Kihara, "Linear Simplexification by Square-Quantity Readout," v0.1-r1, Zenodo DOI: 10.5281/zenodo.20785540, 2026.

[2] N. Kihara, "An Observation on the Isolated Peak Wave of a Constant-Amplitude Odd-Harmonic Sum on a Half-Wavelength Phase Interval and Its Localization," v0.6, Zenodo DOI: 10.5281/zenodo.21073985, 2026.

[3] N. Kihara, "Waveform Invariance and a Contrast Law for the Two-Copy Common-Relative-Phase Superposition of a Half-Wavelength Odd-Harmonic Isolated-Peak Wave," v0.1, Zenodo DOI: 10.5281/zenodo.20923462, 2026.

### Related standard material

[4] T. Takagi, *Kaiseki Gairon* (Introduction to Analysis), 3rd revised ed., Iwanami Shoten, 1961.

[5] Standard textbook material on the polar form of complex numbers, Euler's formula, the complex norm, and standard quadratic forms.

---

## Appendix A: Minimal chain of formulas (revised)

1. In the complex plane, set $Z=x+iy$ (role (1) of $i$: bundling).
2. By the conjugate product, $Z\bar Z=x^2+y^2=\rho^2$.
3. Bundle real components two at a time and add on an equal footing: $\displaystyle\sum_j Z_j\bar Z_j=\sum_n x_n^2=R^2$. An invisible axis merely joins as an equal-footing axis with the same sign.
4. Confirm the difference between conjugate product and square: for real coordinates $x\bar x=x^2$; for imaginary coordinates they diverge, $X_h\bar X_h=z_h^2$ versus $X_h^2=(i z_h)^2=-z_h^2$ (role (2) of $i$). The minus sign arises only here.
5. The zero-sum-of-squares form adding the imaginary-coordinate square to the conjugate-norm backbone: $\displaystyle\sum_j Z_j\bar Z_j+\sum_k (i z_k)^2=0$. Non-triviality requires an imaginary-valued coordinate.
6. Expand and rearrange: $\displaystyle\sum_j Z_j\bar Z_j=\sum_k z_k^2$, i.e. $\displaystyle\sum_{n} x_n^2=\sum_{k} z_k^2$. The minus sign is a consequence of $i^2$. The $Z\bar Z$ type and $X^2$ type are not confused, and no double definition of $R$ occurs.

---

## Appendix B: Candidate connections for a next paper (outside the scope of this note)

This note was limited to observations on algebraic arrangement. Candidates that a next paper might examine (none treated here):

1. Organizing the geometric motivation for introducing imaginary-valued coordinates (a criterion for which axis to place as imaginary).
2. The geometry of the solution set (a complex quadric hypersurface) of the zero-sum-of-squares $\sum_n X_n^2=0$.
3. Projection onto the real subspace and a generalization of the role of the invisible axis.
4. Comparison with the projection radius of the central-projection series (comparison only; no identification).
5. Making explicit the additional assumptions required if one proceeds to a physical interpretation (this note does not proceed there).
