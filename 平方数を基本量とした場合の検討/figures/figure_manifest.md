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

---

# Odd-Harmonic Localization Figure

Belongs to a separate short paper, `paper_odd_harmonic_localization_ja_v0_1.md`.

## Figure 1 (odd-harmonic paper)

- Files: `fig01_odd_harmonic_localization.png`, `fig01_odd_harmonic_localization.svg`
- Role: normalized squared amplitude of the equal-amplitude odd-harmonic sum
  `S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi) = sin((N+1) phi) / (2 sin phi)`
  for `N = 99, 999, 9999`, showing the central peak narrowing as `1/N`.
- Caption: Localization of the equal-amplitude odd-harmonic sum.
- Rebuild: `python3 make_odd_harmonic_figure.py`
- Provenance: supersedes the original hand-made figure kept at
  `source_originals/chatgpt_original_odd_harmonic_figure.png` (which was
  mislabeled "図2"; the regenerated figure is the canonical "図1").
