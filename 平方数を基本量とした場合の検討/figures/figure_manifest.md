# Square-Quantity Readout Figures

This folder contains publication-ready PNG and SVG figures for the first
derivation section of the square-quantity readout paper. Labels and captions
are in English so the same files can be reused in an English manuscript.

## Figure 1

- Files: `fig01_2d_square_map.png`, `fig01_2d_square_map.svg`
- Role: first visual proof that the positive arc of `x^2+y^2=1` is read as the
  straight simplex edge `X+Y=1` under `(X,Y)=(x^2,y^2)`.
- Caption: Squaring sends a curved positive arc to a linear simplex edge.

## Figure 2

- Files: `fig02_3d_square_map.png`, `fig02_3d_square_map.svg`
- Role: three-dimensional analogue of Figure 1. The positive octant of
  `x^2+y^2+z^2=1` is read as the standard simplex `X+Y+Z=1`.
- Caption: In three dimensions, the positive spherical octant reads as the
  standard simplex.

## Figure 3

- Files: `fig03_motion_readouts.png`, `fig03_motion_readouts.svg`
- Role: algebraic readout example. Equations placed on the square-quantity side
  recover uniform-motion and constant-acceleration forms after the positive
  square-root readout.
- Caption: Elementary motion forms are recovered by square-root readout.

## Figure 4

- Files: `fig04_quadratic_readings.png`, `fig04_quadratic_readings.svg`
- Role: geometric vocabulary for later sections. The sign pattern of a quadratic
  form selects an elliptic, double-cone, or hyperbolic reading.
- Caption: Different root-side geometries can share a simple quadratic-form
  origin.

## Figure 5

- Files: `fig05_area_coefficient_ks.png`, `fig05_area_coefficient_ks.svg`
- Role: exact area correction used when the derivation reaches area-sensitive
  terms. The plotted coefficient is
  `k_s(R)=R^2 [4 arccos(-tan^2(1/(2R))) - 2 pi]`.
- Caption: Exact area coefficient for a unit geodesic square.

## Rebuild

Run:

```bash
python3 make_square_quantity_figures.py
```

from this folder.
