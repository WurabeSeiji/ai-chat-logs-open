# Stage Decomposition of the Self-Consistent Inflation Mechanism and Its Full-Range Reproduction for N=3..40 — A Reproduction Specification (Overview)

**Author:** Noriaki Kihara<br>
**ORCID:** 0009-0004-6753-4020<br>
**Date:** September 5, 2026<br>
**Version DOI:** 10.5281/zenodo.22317636<br>
**Concept DOI:** 10.5281/zenodo.22317635 (this overview and Chapters 1–3 share this Concept DOI)<br>
**Positioning:** Reproduction specification (an overview plus three chapter papers). Contains no physical interpretation

---

## 1. Purpose

The purpose of this set of specifications is **to fix, in a reproducible form, the
mechanism by which inflation-like development occurs — the phenomenon in which a seed
too small to be measured (dormant fraction H⊥/H ~ 10⁻³⁰) is exponentially amplified by
the dynamics alone over dozens of orders of magnitude until saturation.**

For this phenomenon, observed at N = 40, 300, 1000 in the canonical paper of
2026-07-22 ("Onset and Threefold Classification of Outcomes of Spontaneous Splitting",
Version DOI 10.5281/zenodo.21486234 / Concept DOI 10.5281/zenodo.21486233), the
components of the dynamics were isolated one factor at a time. The confirmed necessary
condition is the simultaneous presence of **stage 2 (amplitude normalization of the
generator) and stage 3 (removal of the cos symmetric part of the generator = real
orthogonal rotation)** (the audit record of the deletion controls is Chapter 2,
Appendix A). That composition (the stage-1+2+3 dynamics) was then swept from static
parents over the whole range N=3..40, fixing its universality (Chapter 2).

This overview gives only the purpose and the inventory of the uploaded structure.
**The main chapter is Chapter 2**; the equation–program–data correspondences, the
verification gates, and the observations are all recorded in the chapter papers.

## 2. Paper Structure (under this Concept DOI)

| Part | File | Content | Equations |
|---|---|---|---|
| Overview (this document) | `paper_overview/overview_stage123_sweep_en.md` | purpose, structure, upload inventory | — |
| Chapter 1 | `paper_ch1_static_parents/ch1_static_parents_en.md` | generation of the static parent data (seed formula, make_parent, zero-closure seed, Z0; bit-identity gate against the July canonical) | Eqs. 1–9 |
| **Chapter 2 (main)** | `paper_ch2_sweep_dynamics/ch2_sweep_stage123_en.md` | definition of the stage-1+2+3 dynamics (mathematics of the old dynamics; differences of the old/new rotation maps; design hypothesis of the denominators; the dormant fraction as indicator; correspondence with the complex-plane figures §2.6) and the N=3..40 sweep. Appendix A = audit table of deletion controls; Appendix B = program lineage | Eqs. 10–25 |
| Chapter 3 | `paper_ch3_complex_plane/ch3_complex_plane_en.md` | implementation and observations of the three complex-plane readout figures (step 0, final step, zoom into the condensed center) | Eqs. 26–28 |

Japanese originals (`*_ja.md`) are included alongside the English versions.

## 3. Zenodo Upload Inventory

### 3.1 Papers

- The overview and Chapters 1–3 as md (Japanese originals and English versions) and
  tex/PDF (generated at publication).

### 3.2 Main Package `N3_N40_stage123_sweep_20260905/` (~476MB)

| Category | Content |
|---|---|
| Programs (6) | `make_static_parents_N3_N40_v1.py` (parent generation), `run_N3_N40_stage123_v1.py` (sweep body), `check_sweep_inputs_v1.py` (input gate), `analyze_sweep_summary_v1.py` (aggregation), `plot_complex_plane_N3_N40_stage123_v1.py` (plotting), `run_n_scaling_lowrank_v1.py` (bit-identical copy of the July canonical engine) |
| One-shot reproduction | `run_all.sh` (parents → sweep → gate → aggregation → plots) |
| Initial data | `parents/`: 38 static-parent npz + ledger `parents_summary.csv` (~636KB) |
| Result data | `results/`: 228 state npz (all 501 steps × all waves saved), timeseries CSV (114,228 rows), summary CSV (228 rows), aggregation JSON, RUN_METADATA |
| Figures (4) | the H⊥/H denominator-control figure (target figure) and the three complex-plane grids (step0 / final / zoom) |
| Ledgers | `README.md` (progress record, chapter registry), `SHA256SUMS.txt` (canonical SHA256 of all files) |

### 3.3 Companion Folders for Reproducibility (the artifacts behind Appendices A and B)

| Folder | Content | Size |
|---|---|---|
| `ChatGPT_denominator_controls_N40_selfcontrol_20260904/` | the complete artifacts of the N=40 single-factor experiments that fixed the stage composition (62 files): self-control, static-parent substitution (stage-1 baseline), stage 2, stage 3, stage-2 deletion, stage-2-at-init-only, σ clock — each with programs / results / figures / README / SHA256SUMS | ~238MB |
| `自発的分裂予備実験_v1_N40対照実験系_20260904/` | the gated reproduction of the July canonical N=40 run (fcurve bit-identical, all rows), the canonical static parent `parent_static_N40_makeparent_20260904.npz`, complex-plane figures, and the inflation figure | ~1.5MB |

### 3.4 Previously Published Items Referenced by Citation (not re-uploaded)

- The July canonical paper, programs, and figures: "Onset and Threefold Classification
  of Outcomes of Spontaneous Splitting in N-body Relation-Wave Closed Systems" v1.0
  (2026-07-22). **Version DOI 10.5281/zenodo.21486234 / Concept DOI
  10.5281/zenodo.21486233.** The record already bundles the reproduction-program zip
  (`run_spontaneous_splitting_largeN_v1.py`, the engine, etc.) and the figures.

## 4. The Backbone of Verification — the Chain of Gates

The reproducibility of this set is connected to the July canonical by the following
bit-level chain (definitions and measurements in the respective chapters):

1. **July canonical ⇔ reproduction run**: the unmodified rerun in
   `自発的分裂予備実験_v1_N40対照実験系_20260904` matches the July fcurve
   bit-identically over all 3,512 rows (GATE1/GATE2).
2. **Reproduction run ⇔ static parent**: the v, g, Z0 of the canonical static parent
   are bit-identical to the initial values of the July run (Chapter 1, G1). The N=40
   parent generated in Chapter 1 is bit-identical to the canonical static parent (the
   npz SHA256 also coincide).
3. **Static parents ⇔ sweep**: the stored `Z[0]` of all 228 runs are bit-identical to
   the static parents' `Z0` for each N (Chapter 2, G1; checked 228 / MISMATCH 0).

The canonical SHA256 of every file is the `SHA256SUMS.txt` bundled with each package.

## 5. Execution Environment (common)

- Python 3.9.6 (venv), numpy 2.0.2 (BLAS/LAPACK: macOS Accelerate), matplotlib
- macOS 26.3.1 (arm64, Darwin 25.x); RNG: numpy PCG64 (`default_rng`)

---
(End of the overview. For all details, see Chapters 1–3 and the bundled README and
SHA256SUMS files.)
