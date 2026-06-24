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
- Role: squared amplitude (squared norm) `I_N = |S_N|^2` of the equal-amplitude
  odd-harmonic sum
  `S_N(phi) = sum_{m=0}^{(N-1)/2} cos((2m+1) phi) = sin((N+1) phi) / (2 sin phi)`
  for `N = 99, 999, 9999`, NORMALIZED to peak = 1.0 (divide by the array max,
  the value at phi=0). The amplitude is evaluated, then squared, then divided
  by its max -- no shortcut; the closed form is verified against the direct
  term-by-term sum to < 1e-9 first (verify_closed_form), so squaring the closed
  form and squaring the direct sum agree. Shows the central peak narrowing 1/N.
- X-range / domain: the fundamental domain is the HALF-WAVELENGTH interval
  `phi in [-pi/2, pi/2]` (width pi = 180 deg). Phase axis in PERCENT with
  100% = 180 deg, so +/-50% <-> +/-90 deg <-> +/-pi/2. Every odd harmonic is 1
  at the center and 0 at the ends (cos((2m+1)pi/2)=0), so S_N is max at 0 and
  zero at +/-90 deg: a single isolated pulse, no endpoint peaks even when
  squared. S_N is antiperiodic (S_N(phi+pi) = -S_N(phi)), so I_N=|S_N|^2 has
  period pi and the half-wavelength interval is one period. (The earlier
  100% <-> 2*pi mapping plotted the full 360-deg circle = TWO periods, wrongly
  showing two pulses; 360 deg is the wrong fundamental domain.)
- Display window: zoomed to +/-10% (= +/-18 deg); outside that the squared
  amplitude is essentially zero, so the full +/-50% would be flat baseline.
- Caption: Figure 1. Isolated peak wave: the equal-amplitude odd-harmonic sum.
- Labels: ENGLISH, so the figure can be reused as-is in an English manuscript.
- Rebuild: `python3 make_odd_harmonic_figure.py`
- Provenance: supersedes two hand-made source figures kept under
  `source_originals/` — `chatgpt_original_odd_harmonic_figure.png` and
  `user_japanese_odd_harmonic_figure_zu2.png`. Both were mislabeled "図2", and
  the second one, despite its title saying 振幅二乗 (squared), plots the LINEAR
  amplitude `|S_N|/max` (its first sidelobe ~0.22, matching `|sin u/u|`, not the
  ~0.047 of the squared kernel). The script figure is the canonical "図1" and
  is squared (`|S_N|^2`), consistent with the paper's definition `I_N=|S_N|^2`.

## Figure 2 (odd-harmonic paper) — 1/N scaling, honest overlay

- Files: `fig02_odd_harmonic_scaling.png`, `fig02_odd_harmonic_scaling.svg`
- Role: ALL THREE curves (N=99 red, N=999 blue, N=9999 green) overlaid in EVERY
  panel; only the phase window changes (+/-10%, +/-1%, +/-0.1%). In each window
  exactly one main lobe fits (99@+/-10%, 999@+/-1%, 9999@+/-0.1%); the other two
  appear x10 wider (broad flat top) or x10 narrower (thin spike). This is the
  honest demo of W(N) ~ 1/N -- NOT three separately-tuned plots that trivially
  coincide (that earlier version looked "too convenient"). Shared vertical scale
  (peak = 1.0); data outside each window is clipped.
- Numerics (printed by the script, all ratios exactly x10):
  first-zero half-width 100/(N+1)% = +/-1.00 / 0.100 / 0.0100 %;
  1%-level half-width (first descent through 1%, eq 4.2, ~90.7/(N+1)%)
  = +/-0.908 / 0.0908 / 0.00908 %. (NOT the last-crossing of the sidelobe
  envelope, which is a different, larger measure.)
- Labels: ENGLISH. Rebuild: `python3 make_odd_harmonic_scaling_figure.py`.
- Paper: this is 図2 in `paper_odd_harmonic_localization_ja_v0_1.md` (§4.1).
