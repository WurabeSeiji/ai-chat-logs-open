# Paper 3: Closed Four-Degree-of-Freedom Structure and Its Correspondence with 4-Dimensional Lattice Counting
## A Geometric Organization from the 5-Component Sum-of-Squares Constraint to the Unit-Cell Counting Region

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v0.2  
**Date**: June 2026  
**DOI (this version)**: 10.5281/zenodo.20589262  
**Concept DOI**: 10.5281/zenodo.20589261  
**Paper 1 (supplemented)**: 10.5281/zenodo.20588037  
**License**: CC BY 4.0

---

## Abstract

In Paper 1 we observed that, in the wavelength space or the frequency space, imposing the 5-component sum-of-squares condition

$$
\sum_{n=1}^{5}x_n^2=R^2
$$

produces a closed geometric structure with four degrees of freedom. In Paper 2 we counted the number of unit 4-cells of side 1 that are fully inscribed in the 4-dimensional ball region of radius $R$,

$$
N_0(R)
=
\#\left\{
k\in\mathbb{Z}^4
\mid
\sum_{i=1}^{4}
\left(|k_i|+\frac12\right)^2
\le R^2
\right\} .
$$

At first sight, the 5-component sum-of-squares constraint of Paper 1 defines a 4-dimensional hypersphere inside a 5-dimensional space, whereas the counting of Paper 2 is an interior-cell counting on a 4-dimensional ball region. This paper is a short supplement that organizes this correspondence without physical interpretation.

The point of this paper is simple. The 5-component sum-of-squares constraint defines a closed object with four degrees of freedom. To actually carry out lattice-cell counting, those four degrees of freedom are treated as a 4-dimensional coordinate region. Hence the counting of Paper 2 does not count the values of $\lambda_n$ or $\nu_n$ themselves as physical quantities; it counts the number of 4-dimensional lattice cells corresponding to the closed four-degree-of-freedom structure of radius $\Lambda$ or $\mathcal{N}$.

Moreover, for a point that already satisfies the sum-of-squares constraint, a spherical projection onto the sphere of the same radius leaves the value unchanged. That is, in the wavelength space

$$
\lambda'=\Lambda\frac{\lambda}{\|\lambda\|}=\lambda
$$

and in the frequency space

$$
\nu'=\mathcal{N}\frac{\nu}{\|\nu\|}=\nu .
$$

Using this fact, we confirm that the projection is not a transformation of values but a geometric description for reading a point that satisfies the sum-of-squares constraint as a point on a closed four-degree-of-freedom structure of constant radius.

---

## Keywords

5-component sum-of-squares constraint, four degrees of freedom, 4-dimensional lattice, unit-cell counting, spherical projection, wavelength space, frequency space, observational model

---

## 1. Purpose

In Paper 1 we imposed the 5-component sum-of-squares condition in the wavelength space or the frequency space.

In the wavelength space,

$$
\sum_{n=1}^{5}\lambda_n^2=\Lambda^2 ,
$$

and in the frequency space,

$$
\sum_{n=1}^{5}\nu_n^2=\mathcal{N}^2 .
$$

These each give a closed constraint with four degrees of freedom in the 5-component space.

In Paper 2, on the other hand, we counted the number of unit 4-cells of side 1 fully inscribed in the 4-dimensional ball region of radius $R$.

The purpose of this paper is to make explicit, as concisely as possible, the geometric organization required to connect these two.

This paper does not treat any connection with physical spacetime, mass, energy, momentum, gravity, or physical constants. What we treat is only the following correspondence.

$$
\text{5-component sum-of-squares constraint}
\quad\longrightarrow\quad
\text{four-DOF structure}
\quad\longrightarrow\quad
\text{4-dimensional lattice-cell counting}
$$

---

## 2. The 5-Component Sum-of-Squares Constraint

In general, for a 5-component vector

$$
x=(x_1,x_2,x_3,x_4,x_5)\in\mathbb{R}^5 ,
$$

imposing

$$
\sum_{n=1}^{5}x_n^2=R^2
$$

defines the 4-dimensional hypersphere of radius $R$,

$$
S^4_R
=
\left\{
x\in\mathbb{R}^5
\mid
\sum_{n=1}^{5}x_n^2=R^2
\right\} .
$$

This lies inside the 5-component space, but because there is one constraint, the number of degrees of freedom is four.

In this paper we do not interpret this object as a physical space. We treat it simply as a closed geometric structure with four degrees of freedom defined by the 5-component sum-of-squares constraint.

---

## 3. The 4-Dimensional Region for Counting

What we counted in Paper 2 are the unit cells on the 4-dimensional lattice

$$
\mathbb{Z}^4 .
$$

The region used there is

$$
B^4_R
=
\left\{
u\in\mathbb{R}^4
\mid
\sum_{i=1}^{4}u_i^2\le R^2
\right\} .
$$

This is a 4-dimensional ball region, the counting region for judging whether a unit cell is fully inscribed.

Taking the center of a unit 4-cell of side 1 to be

$$
k=(k_1,k_2,k_3,k_4)\in\mathbb{Z}^4 ,
$$

the condition that this cell be fully inscribed in $B^4_R$ is

$$
\sum_{i=1}^{4}
\left(|k_i|+\frac12\right)^2
\le R^2 .
$$

Hence the fully-inscribed cell count is

$$
N_0(R)
=
\#\left\{
k\in\mathbb{Z}^4
\mid
\sum_{i=1}^{4}
\left(|k_i|+\frac12\right)^2
\le R^2
\right\} .
$$

---

## 4. Why Use the 4-Dimensional Ball Region

The 5-component sum-of-squares constraint defines $S^4_R$. The lattice counting, on the other hand, is carried out inside $B^4_R$.

This is not a contradiction.

$S^4_R$ is a closed object with four degrees of freedom, and to actually carry out lattice-cell counting one must treat those four degrees of freedom as a 4-dimensional counting region.

We therefore use, as the counting region,

$$
B^4_R
=
\left\{
u\in\mathbb{R}^4
\mid
\sum_{i=1}^{4}u_i^2\le R^2
\right\} .
$$

This operation does not mean that the values of $\lambda_n$ or $\nu_n$ become physical coordinates inside the 4-dimensional ball. What Paper 2 counts is the 4-dimensional lattice indices, or the number of unit cells, introduced for a constraint structure that has four degrees of freedom.

That is, the counting of Paper 2 is read as

$$
\text{4-dimensional unit-cell counting}
$$

for

$$
\text{the closed four-DOF structure of radius }R .
$$

---

## 5. Re-reading in the Wavelength and Frequency Spaces

In the wavelength space the radius is

$$
R=\Lambda ,
$$

so the wavelength-side fully-inscribed cell count can be written as

$$
N_0(\Lambda) .
$$

In the frequency space the radius is

$$
R=\mathcal{N} ,
$$

so the frequency-side fully-inscribed cell count can be written as

$$
N_0(\mathcal{N}) .
$$

That is,

$$
R=\Lambda
\quad\Rightarrow\quad
N_0(R)=N_0(\Lambda)
$$

$$
R=\mathcal{N}
\quad\Rightarrow\quad
N_0(R)=N_0(\mathcal{N}) .
$$

Here too we do not interpret $\Lambda$ or $\mathcal{N}$ as physical quantities. They are the radii of the sum-of-squares constraint in the wavelength space or the frequency space.

---

## 6. Viewed as a Spherical Projection

As a general form of spherical projection, the map that sends an arbitrary nonzero vector

$$
y\in\mathbb{R}^5,\qquad y\ne0
$$

onto the sphere of radius $R$ can be written

$$
\Pi_R(y)=R\frac{y}{\|y\|} .
$$

This map is a spherical projection that moves a point onto the sphere of radius $R$ while preserving the direction seen from the origin.

In this paper we consider the same form of projection for the wavelength space and the frequency space.

In the wavelength space,

$$
\lambda'=\Pi_\Lambda(\lambda)
=
\Lambda\frac{\lambda}{\|\lambda\|} ,
$$

and in the frequency space,

$$
\nu'=\Pi_{\mathcal{N}}(\nu)
=
\mathcal{N}\frac{\nu}{\|\nu\|} .
$$

The important point here is that, in the setup of Paper 1, $\lambda$ and $\nu$ already satisfy the sum-of-squares constraints

$$
\sum_{n=1}^{5}\lambda_n^2=\Lambda^2,
\qquad
\sum_{n=1}^{5}\nu_n^2=\mathcal{N}^2 .
$$

Therefore

$$
\|\lambda\|=\Lambda,
\qquad
\|\nu\|=\mathcal{N} .
$$

Hence, in the wavelength space,

$$
\lambda'
=
\Lambda\frac{\lambda}{\|\lambda\|}
=
\Lambda\frac{\lambda}{\Lambda}
=
\lambda ,
$$

that is,

$$
\boxed{\lambda'=\lambda} .
$$

Similarly, in the frequency space,

$$
\nu'
=
\mathcal{N}\frac{\nu}{\|\nu\|}
=
\mathcal{N}\frac{\nu}{\mathcal{N}}
=
\nu ,
$$

that is,

$$
\boxed{\nu'=\nu} .
$$

Thus, for a point that already satisfies the sum-of-squares constraint, projecting onto the sphere of the same radius does not change the component values.

This point is important in the organization of this paper. The projection here is not an operation that transforms the values of $\lambda$ or $\nu$ into other values. The values after projection are

$$
\lambda\to\lambda'=\lambda,
\qquad
\nu\to\nu'=\nu ,
$$

and the values themselves are invariant.

Therefore the spherical projection in this paper is a geometric description for reading a 5-component vector that satisfies the sum-of-squares constraint as a point on a closed four-degree-of-freedom structure of constant radius. It is not a transformation of values, but a description that makes explicit the arrangement of constraint-satisfying points.

## 7. Connection with Paper 2

The counting of Paper 2 is

$$
N_0(R)
=
\#\left\{
k\in\mathbb{Z}^4
\mid
\sum_{i=1}^{4}
\left(|k_i|+\frac12\right)^2
\le R^2
\right\} .
$$

By the organization of this paper, this $R$ can be read as $\Lambda$ on the wavelength side and as $\mathcal{N}$ on the frequency side.

Hence the counting of Paper 2 corresponds, as

$$
N_0(\Lambda)
$$

or

$$
N_0(\mathcal{N}) ,
$$

to the unit-cell counting of the four-degree-of-freedom structure in the wavelength space or the frequency space.

However, this is not a counting of physical space or physical quantities. It is purely a lattice-cell counting for the dual geometry carrying the sum-of-squares constraint.

---

## 8. Positioning of This Paper

This paper is a supplement that makes explicit the connection between Paper 1 and Paper 2.

Paper 1 treated the reciprocal duality and the sum-of-squares constraint of the wavelength space and the frequency space.

Paper 2 counted the fully-inscribed unit-cell count on the 4-dimensional lattice by a radius sweep.

This paper organized, without physical interpretation, the correspondence

$$
\text{5-component sum-of-squares constraint}
\quad\longrightarrow\quad
\text{four-DOF structure}
\quad\longrightarrow\quad
\text{4-dimensional lattice counting}
$$

that lies between them.

---

## 9. Limitations

This paper does not treat the following.

1. Geodesic cell partitions on $S^4_R$  
2. Lattice partitions based on surface area or geodesic distance  
3. Correspondence with physical spacetime  
4. Correspondence with energy, momentum, mass, gravity, etc.  
5. Correspondence with physical constants  
6. Derivation of the uncertainty width $\delta$  

This paper is merely a geometric organization for treating the four-degree-of-freedom structure that appears under the 5-component sum-of-squares constraint as a 4-dimensional lattice-counting region.

---

## 10. Conclusion

In this paper we organized the correspondence between the 5-component sum-of-squares constraint of Paper 1 and the 4-dimensional lattice counting of Paper 2.

The 5-component sum-of-squares constraint

$$
\sum_{n=1}^{5}x_n^2=R^2
$$

defines a closed geometric structure with four degrees of freedom.

For the actual lattice counting, those four degrees of freedom are treated as the 4-dimensional ball region

$$
B^4_R
=
\left\{
u\in\mathbb{R}^4
\mid
\sum_{i=1}^{4}u_i^2\le R^2
\right\} ,
$$

and we count the fully-inscribed number of unit 4-cells,

$$
N_0(R) .
$$

Since one can read $R=\Lambda$ in the wavelength space and $R=\mathcal{N}$ in the frequency space, the counting of Paper 2 can be used as

$$
N_0(\Lambda)
$$

or

$$
N_0(\mathcal{N}) .
$$

Moreover, the spherical projection of the same radius applied to a point already satisfying the sum-of-squares constraint is the identity map,

$$
\lambda'=\lambda,\qquad
\nu'=\nu .
$$

Therefore the projection in this paper is not an operation that transforms values, but a geometric description for reading a point that satisfies the sum-of-squares constraint as a point on a closed four-degree-of-freedom structure.

---

## References

1. Snyder, J. P. (1987). *Map Projections—A Working Manual*. U.S. Geological Survey Professional Paper 1395. U.S. Government Printing Office.

2. Haines, L. M. (2024). Stereographic Projections for Designs on the Sphere. *arXiv:2401.05931*.

3. Coxeter, H. S. M. (1969). *Introduction to Geometry* (2nd ed.). Wiley.

4. do Carmo, M. P. (1976). *Differential Geometry of Curves and Surfaces*. Prentice-Hall.

---

## Note

The references are listed not as physical grounds for this paper but as general background on spherical projection, spherical geometry, and coordinatization. The claim of this paper is limited to making explicit the geometric correspondence between the 5-component sum-of-squares constraint and the 4-dimensional lattice counting; it is not a physical theory.
