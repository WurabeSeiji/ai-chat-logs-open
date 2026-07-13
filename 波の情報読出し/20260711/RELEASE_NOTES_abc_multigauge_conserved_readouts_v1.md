# Release Notes: ABC Multigauge Conserved Readouts v1/v2

**Date:** 2026-07-11
**Author:** Noriaki Kihara
**Repository:** `ai-chat-logs-open`
**Latest Zenodo Version DOI:** 10.5281/zenodo.21332875
**Zenodo Concept DOI:** 10.5281/zenodo.21308049
**Latest Zenodo Record:** https://zenodo.org/records/21332875
**Previous v1 Version DOI:** 10.5281/zenodo.21308050
**Zenn Article:** `articles/abc-multigauge-conserved-readouts.md`
**note Japanese Article:** https://note.com/kiharanoriaki/n/nd5d3777a6e48
**note English Article:** https://note.com/kiharanoriaki/n/nd10d4b8d627d
**Facebook Japanese Post:** https://www.facebook.com/kihara.noriaki/posts/pfbid02ujTqp1HVzQLNA8xyzv9vARuzUmtx2c47kjHBk6aKVVTanFZAMm69BhR9GU6e6ryAl
**Facebook English Post:** https://www.facebook.com/kihara.noriaki/posts/pfbid02PAsfWWBXcyKm5cciEusyKZWacz9ey1yjwpyMoPX5YFvi52AVYmLjhmBYHtRGchjzl
**X Japanese Post:** https://x.com/NoriakiKihara/status/2075898071238791620
**X English Post:** https://x.com/NoriakiKihara/status/2075898711432106276

---

## Summary

This release publishes the additional Wave Information Readout paper:

> Construction Experiment of Multigauge Interference Readout Conserved Quantities in an ABC Closed Phase System

The paper tests whether mass-like, momentum-like, and energy-like conserved readouts can be constructed from multigauge interference in an ABC closed phase system without assuming background coordinates.

The paper reconstructs `p_read`, `E_read`, and `R_read` from multiple gauges, diagnoses that simple `q -> -q` reversal fails to preserve `R*p` under asymmetric `R`, and constructs an `R`-weighted generalized collision map conserving `R*p` and `R*p^2`.

---

## v2 Update

V2 recalculates the experiment series after changing the formula for the fermion-like reflection map.

The Concept DOI is maintained. The V2 Version DOI is `10.5281/zenodo.21332875`.

The Japanese/English Markdown, TeX/PDF, numerical outputs, and reproducibility bundle were updated to the V2 outputs.

---

## Primary Outputs

| File | Role |
|---|---|
| `ABC閉鎖位相系における多ゲージ干渉読出し保存量の構成実験 v2.md` | Japanese main paper |
| `abc_closed_phase_system_multigauge_conserved_readouts_en.md` | English main paper |
| `abc_closed_phase_system_multigauge_conserved_readouts_ja.tex` | Japanese TeX source |
| `abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf` | Japanese PDF |
| `abc_closed_phase_system_multigauge_conserved_readouts_en.tex` | English TeX source |
| `abc_closed_phase_system_multigauge_conserved_readouts_en.pdf` | English PDF |
| `abc_multigauge_conserved_readouts_zenodo_deposit_v2.json` | Zenodo metadata without token |
| `abc_multigauge_conserved_readouts_zenodo_deposit_v2.json` | Zenodo draft deposit response without token |
| `abc_multigauge_conserved_readouts_zenodo_uploads_v2.json` | Zenodo upload result summary without token |
| `abc_multigauge_conserved_readouts_zenodo_published_record_v2.json` | Zenodo published record response without token |
| `abc_multigauge_conserved_readouts_publication_bundle_v2.zip` | Reproducibility bundle |
| `articles/abc-multigauge-conserved-readouts.md` | Zenn article source |
| `note_article_abc_multigauge_conserved_readouts_ja.md` | note Japanese article source |
| `note_article_abc_multigauge_conserved_readouts_en.md` | note English article source |

---

## Numerical Experiment Coverage

1. Single ABC collision multigauge readout of `p,E,R`.
2. Repeated symmetric ABC collisions.
3. Readout robustness across gauge centers, widths, phases, and gains.
4. Asymmetric `R` diagnostic showing simple reversal failure for `R*p`.
5. `R`-weighted generalized elastic collision map.
6. Non-unit and asymmetric phase-gradient sweep.
7. Repeated generalized collisions.
8. Readout noise robustness with zero-mean gauge perturbations and common-bias detection.
9. Extreme `R` ratio sweep from `R_B/R_A=0.015625` to `64.0`.

---

## Main Findings

- In the single ABC collision, `p,E,R` were reconstructed from multiple gauges with maximum errors `2.5202062658991053e-14`, `2.2315482794965646e-14`, and `4.440892098500626e-16`.
- In eight repeated symmetric collisions, `p` reversal, `E,R` preservation, identification oscillation preservation, and compensated closure were maintained.
- Under asymmetric `R`, simple `q -> -q` reversal was diagnosed as failing to preserve `R*p`; the maximum weighted momentum collision error was `16.000000000000036`.
- The generalized `R`-weighted collision map preserved `R*p` and `R*p^2` across eight asymmetric amplitude cases with maximum errors `2.3803181647963356e-13` and `1.4086509736443986e-12`.
- The generalized velocity sweep preserved `R*p`, `R*p^2`, and relative-gradient reversal across nine non-unit and asymmetric initial conditions.
- Repeated generalized collisions preserved the same quantities across four cases with up to six AB collisions.
- Zero-mean readout noise was cancelled by multigauge averaging, while common readout bias above the detection floor was detected.
- The extreme `R` sweep preserved the generalized map across `R_B/R_A=0.015625` to `64.0`.
- The integration summary reported all nine experiments as `valid`, with no single-gauge-only judgment used.

---

## Scope

This release does not claim a derivation of standard momentum, standard energy, or standard mass.

The confirmed result is internal to the Wave Information Readout series: in the tested ABC closed phase system, `p_read`, `E_read`, `R_read`, `R*p`, and `R*p^2` can be consistently constructed as conserved readouts from multigauge interference.

The correspondence map from these readout quantities to standard physical quantities is left as a subsequent task.

---

## Build Notes

- Markdown was converted to standalone TeX with Pandoc.
- PDF generation was performed under `/tmp/abc_multigauge_readout_build` to avoid Google Drive filesystem issues during LaTeX compilation.
- XeLaTeX was used for both Japanese and English PDFs.
- Japanese CJK font: `HaranoAjiMincho-Regular.otf`.
- Monospace font: `DejaVuSansMono.ttf`.
- The English TeX source includes `xeCJK` support so Japanese file names in appendices render correctly.
- The TeX build succeeded with non-fatal table-width, font-shape, and overfull-line warnings caused mainly by long identifiers and filenames.

---

## Publication Status

Publication has been completed.

- Zenodo version DOI: https://doi.org/10.5281/zenodo.21332875
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.21308049
- Zenodo record: https://zenodo.org/records/21332875
- Zenn article source: `articles/abc-multigauge-conserved-readouts.md`
- note Japanese article: https://note.com/kiharanoriaki/n/nd5d3777a6e48
- note English article: https://note.com/kiharanoriaki/n/nd10d4b8d627d
- Facebook Japanese post: https://www.facebook.com/kihara.noriaki/posts/pfbid02ujTqp1HVzQLNA8xyzv9vARuzUmtx2c47kjHBk6aKVVTanFZAMm69BhR9GU6e6ryAl
- Facebook English post: https://www.facebook.com/kihara.noriaki/posts/pfbid02PAsfWWBXcyKm5cciEusyKZWacz9ey1yjwpyMoPX5YFvi52AVYmLjhmBYHtRGchjzl
- X Japanese post: https://x.com/NoriakiKihara/status/2075898071238791620
- X English post: https://x.com/NoriakiKihara/status/2075898711432106276

Zenodo uploaded files:

1. `ABC閉鎖位相系における多ゲージ干渉読出し保存量の構成実験 v2.md`
2. `abc_closed_phase_system_multigauge_conserved_readouts_en.md`
3. `abc_closed_phase_system_multigauge_conserved_readouts_ja.tex`
4. `abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf`
5. `abc_closed_phase_system_multigauge_conserved_readouts_en.tex`
6. `abc_closed_phase_system_multigauge_conserved_readouts_en.pdf`
7. `abc_multigauge_conserved_readouts_zenodo_deposit_v2.json`
8. `abc_multigauge_conserved_readouts_publication_bundle_v2.zip`
