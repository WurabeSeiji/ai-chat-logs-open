# Chapter 3: Complex-Plane Readout Figures — Three Grids: Step 0, Final Step, and Zoom into the Condensed Center

(Chapter 3 of the reproduction papers for the N=3..40 stage-1+2+3 sweep. Shares the
Concept DOI 10.5281/zenodo.22317635 with the overview paper; Version DOI
10.5281/zenodo.22317636. Equation numbers continue from Chapters 1–2 (Eqs. 1–25).
Code blocks are quoted verbatim; comments inside them remain in Japanese, as in the
originals.)

## 1. Purpose

This chapter describes the process that generates, from the state npz files saved by
the Chapter-2 sweep (the Δτ=2π/N runs), three grid figures reading out the per-edge
complex waves z_e ∈ ℂ on the complex plane — (1) step 0, (2) the final step
(step 500), and (3) a zoom into the largest angular cluster at the final step. The
design intent of the three figures (their division of roles corroborating the start,
end, and end-interior of the inflation figure from the configuration side, the caution
that "the complex plane of the figures ≠ the parent plane Π", and the bridging
identity) is given in **Chapter 2, §2.6 (Eq. 25)**; this chapter fixes the
implementation, execution, and observations. The data are read-only; the dynamics and
stored data are never touched.

## 2. Theoretical Background

The mathematics of states and measurement follows Chapter 2 (Eqs. 22–25). The
figure-specific conventions of this chapter are defined as Eqs. 26–28.

**(Eq. 26) Duplicate counting (notation of degeneracy)** — The set of plotted complex
values {z_e} is partitioned into equivalence classes by **rounding the real and
imaginary parts to the 12th decimal place**; classes with c > 1 members are annotated
"xc" at their position:

    class(w) = ( round(Re w, 12), round(Im w, 12) )

Rounding at 12 digits is coarser than the effective double precision (~16 digits): it
identifies differences at the level of rounding error while keeping physically
significant separations (relative 10⁻⁴ or larger in this series) distinct. In the zoom
figure the rounding is tightened to the 15th decimal (Eq. 28), so that only
machine-precision coincidences are counted as "xn".

**(Eq. 27) Panel scaling convention** — The axis range of each panel is the square
[−1.15r, +1.15r] (equal aspect) with r = max_e |z_e| the overall amplitude of that
panel, and **the tick values display the actual values** (no relabeling by normalized
values). The line segments from the origin to each point are guides for reading the
argument and modulus of each wave simultaneously.

**(Eq. 28) Extraction of the largest angular cluster (algorithm of the zoom figure)** —
The final-step state, in coordinates normalized by the overall amplitude, is grouped by
rounding at **coarse resolution 1/100**; the group with the most members is the zoom
target:

    key(w) = ( round(Re w / amp, 2), round(Im w / amp, 2) ),   amp = max_e |z_e|
    C* = argmax_{key} |{ w : key(w) = key }|
    center c = mean(C*),  spread = max_{w∈C*} |w − c|,  window = 1.4 × spread

The panel title prints |C*|, spread, and amp. The spread quantifies whether the
cluster is an exact single point (exact degeneracy) or a bundle of finite width (the
zoom row of the table in Chapter 2, §2.6).

## 3. Implementation Method

- Plotting program `plot_complex_plane_N3_N40_stage123_v1.py` (bundled in this
  package; read-only). The drawing style extends the existing plotting programs of
  this series (the grid style of
  `complex_plane_readout_step0_step2000_20260904/plot_complex_plane_step0_step2000.py`
  and the cluster-zoom algorithm of
  `自発的分裂予備実験_v1/N40_state_readout_20260904/plot_complex_plane_N40_v1.py`) to
  the 8×5 grid of N=3..40; the algorithms themselves are identical (the style lineage
  is unchanged).
- The input is solely the Chapter-2 product `results/hm_N{N}_den_{N}_states_500.npz`.
  **The den=N (Δτ=2π/N) run is used as the representative for each N** (step 0 is
  identical across denominators, so the representative choice matters only for the
  final step; final configurations for other denominators can be read out likewise
  from the stored npz).
- The output is 3 PNG files. No data are rewritten or regenerated.

## 4. Detailed Design

### 4.1 Overall Flow

```
[initialization]  fix the input folder (results/)
[loop]  draw_grid(step=0)  : draw the step-0 state on each panel for N=3..40 (Eqs. 26, 27)
        draw_grid(step=500): same for the final-step state
        zoom grid          : largest angular cluster (Eq. 28) per panel for N=3..40
[finalization]  savefig each figure; switch off the 2 unused panels
```

### 4.2 Overall Data Flow

- **Input**: `results/hm_N{N}_den_{N}_states_500.npz` × 38 (Chapter-2 products;
  read-only). Fields used: `Z` (501×M), `denominator`, `steps` (for consistency
  asserts)
- **Parameters** (all constants inside the program):
  - rounding digits: 12 (Eq. 26) / 15 (inside the zoom, Eq. 28) / coarse resolution
    2 digits (the 1/100 of Eq. 28)
  - scale factor 1.15 (Eq. 27), window factor 1.4 (Eq. 28)
  - grid 8×5 (38 used, 2 off), dpi=180
- **Output**:
  - `fig_complex_plane_step0_N3_N40_stage123.png`
  - `fig_complex_plane_final_N3_N40_stage123.png`
  - `fig_complex_plane_final_zoom_N3_N40_stage123.png`

### 4.3 Individual Processes

#### 4.3.1 Initialization (loading and consistency checks)

```python
    18	def load(N, step):
    19	    d = np.load(os.path.join(IN, f'hm_N{N}_den_{N}_states_500.npz'))
    20	    assert int(d['denominator']) == N and int(d['steps']) == 500
    21	    return np.asarray(d['Z'][step], dtype=np.complex128)
```

- The assert on line 20 checks that "the intended file (den=N, 500 steps) is being
  read"; on mismatch an AssertionError stops execution (no drawing occurs).

#### 4.3.2 Loop (1): Grid Rendering (Eqs. 26, 27)

```python
    26	    for k, N in enumerate(range(3, 41)):
    27	        ax = axs[k]
    28	        z = load(N, step)
    29	        M = N * (N - 1) // 2
    30	        assert z.size == M
    31	        for w in z:
    32	            ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    33	        ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    34	        cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    35	        for (a, b), c in cnt.items():
    36	            if c > 1:
    37	                ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
    38	                            fontsize=5, color='black')
    39	        r = float(np.abs(z).max())
    40	        lim = r * 1.15 if r > 0 else 1.0
    41	        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
```

- Lines 31–32: segments from the origin; line 33: points; lines 34–38: the duplicate
  annotations of Eq. 26; lines 39–41: the scaling convention of Eq. 27 (tick values
  remain actual values; `ticklabel_format` on line 46 only switches the notation to
  scientific).

#### 4.3.3 Loop (2): Zoom into the Largest Angular Cluster (Eq. 28)

```python
    67	    z = load(N, 500)
    68	    amp = float(np.abs(z).max())
    69	    coarse = {}
    70	    for w in z:
    71	        key = (round(float(w.real) / amp, 2), round(float(w.imag) / amp, 2))
    72	        coarse.setdefault(key, []).append(w)
    73	    mem = max(coarse.values(), key=len)
    74	    zz = np.array(mem)
    75	    c = zz.mean()
    76	    dev = np.abs(zz - c)
    77	    spread = float(dev.max())
    78	    win = spread * 1.4 if spread > 0 else amp * 1e-12
...
    90	    ax.set_title(f'N={N}: {len(zz)} waves, dev={spread:.2e} (|z|max={amp:.2e})', fontsize=7)
```

- Lines 70–73: the grouping and largest-cluster selection of Eq. 28; lines 75–78:
  center, spread, window; line 90: printing member count, spread, and amp in the panel
  title (this printed output is the source of the observed values in §6).

#### 4.3.4 Finalization and Exceptional Behavior (elements not formulated as equations)

- Line 40, `lim = r*1.15 if r > 0 else 1.0`: fallback for the degenerate case of all
  points at the origin (does not occur in these data).
- Line 78, `win = spread*1.4 if spread > 0 else amp*1e-12`: window fallback when the
  cluster has a single member (spread=0; occurs at N=3, where the largest cluster is a
  single wave).
- Lines 50–51 and 91–92: of the 8×5=40 panels, the unused 39th and 40th are set to
  `axis('off')`.
- The only exception paths are the asserts (lines 20, 30); there are no other
  branches.

## 5. Execution Results

### 5.1 Reproduction Commands

```bash
cd N3_N40_stage123_sweep_20260905
python3 plot_complex_plane_N3_N40_stage123_v1.py
# or as the final step of ./run_all.sh
```

### 5.2 Execution Environment

- Python 3.9.6 (`.venv/bin/python3`), numpy 2.0.2, matplotlib (Agg rendering)
- macOS 26.3.1 (arm64)

### 5.3 Execution Time

**Under about 1 minute** for the three figures in total (dominated by 38×2 npz loads
and the rendering of up to 780 segments × 38 panels).

### 5.4 Verification Gates

| Gate | Pass condition | Measured | Verdict |
|---|---|---|---|
| G1: input consistency | asserts denominator==N, steps==500, z.size==M pass at every load | passed for all 38 N × 3 figures (no AssertionError) | **PASS** |
| G2: completion | `ALL DONE` printed, 3 PNGs generated | confirmed | **PASS** |

(The bit-level lineage of the input npz themselves is guaranteed by Chapter 2 G1 —
`Z[0]` of all 228 runs bit-identical to the static parents.)

### 5.5 Data

| Item | Content |
|---|---|
| Input | `results/hm_N{N}_den_{N}_states_500.npz` × 38 (Chapter-2 products; read-only) |
| Output location | package root |
| step-0 figure | `fig_complex_plane_step0_N3_N40_stage123.png` (632,949 bytes) |
| final-step figure | `fig_complex_plane_final_N3_N40_stage123.png` (4,582,367 bytes) |
| zoom figure | `fig_complex_plane_final_zoom_N3_N40_stage123.png` (494,906 bytes) |
| SHA256 | the bundled `SHA256SUMS.txt` is canonical |

### 5.6 Figures

The three figures above (8×5 grids; each panel equal aspect, actual-value ticks). How
to read them, and their correspondence with the inflation figure, follow the table in
Chapter 2, §2.6.

## 6. Execution Analysis (objective report and observations only; the numbers
originate from the in-figure printing = the output of this program)

1. **Step-0 figure**: every panel for N=3..40 shows the star shape of two antipodal
   pairs (4 bundles). The bundles have a radial spread of amplitudes. The duplicate
   annotations of Eq. 26 appear only for some small N (x4 at N=4; x2–x3 for N=5–9,
   etc.); exact duplicates vanish as N grows.
2. **Final-step figure**: in every panel the star has disappeared, replaced by a
   nearly equimodular ring-like arrangement (sparse spokes for small N). The Eq. 26
   (12-digit) duplicate annotations do not appear except for a few small-N cases (x2
   at N=4, 5).
3. **Zoom figure**: the largest angular cluster has 1–10 members (N=3: 1 wave;
   N=39, 40: 10 waves); the spread is of order 10⁻⁷–10⁻⁴ for N≥6 (e.g., at N=40,
   dev=2.04e-04 with |z|max=3.64e-02, relative ~5.6×10⁻³). Only the two-wave clusters
   of N=4, 5 coincide at machine precision (dev 10⁻¹⁶–10⁻¹⁰).
4. Summary of observations (facts only): for all N≥6 the final-step angular clusters
   are bundles of finite width, not condensation onto exactly identical complex
   values. As per the division of roles in Chapter 2 §2.6, this records the fine
   organization of the final state that H⊥/H does not measure.

---
(End of Chapter 3. The overview paper is separate.)
