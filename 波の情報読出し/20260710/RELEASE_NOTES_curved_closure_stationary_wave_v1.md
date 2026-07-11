# Release Notes: Curved Closed Stationary Wave Curvature Renormalization v1

**Date:** 2026-07-11
**Author:** Noriaki Kihara
**Repository:** `ai-chat-logs-open`
**Zenodo Version DOI:** 10.5281/zenodo.21304040
**Zenodo Concept DOI:** 10.5281/zenodo.21304039
**Zenodo Record:** https://zenodo.org/records/21304040
**Zenn Article:** `articles/curved-closure-stationary-wave.md`
**note Japanese Article:** https://note.com/kiharanoriaki/n/n2389460836cf
**note English Article:** https://note.com/kiharanoriaki/n/nda3623c44423
**Facebook Japanese Post:** https://www.facebook.com/kihara.noriaki/posts/pfbid037wN39hUdgVY7bVCb6BWkFK86pYqDedPgVxaCxtJbWGk479ZqyVoc9XggvRPLmstwl
**Facebook English Post:** https://www.facebook.com/kihara.noriaki/posts/pfbid02QieNB6aGk3TyRkcreL2z14w62WEFdRbenP47gHuzMoPd4aS7VDPzazbuNJhEGBo4l
**X Japanese Post:** https://x.com/NoriakiKihara/status/2075801193281106375
**X English Post:** https://x.com/NoriakiKihara/status/2075802233153998999

---

## Summary

This release prepares the additional Wave Information Readout paper:

> Curvature Renormalization and Perfect-Reflection Stability by Curved Closed Stationary Waves

The paper tests how odd-harmonic complex waves satisfying the all-positive zero closure `Sigma x_n^2 = 0` behave in curved local cells. Curvature is represented as relative phase leakage inside the local cell. The numerical experiments show that leakage is not absent: it appears in transient states as closure residuals and transmission leakage. Once the system is re-selected into a closed stationary wave, the leakage is absorbed into internal phase and the complete-reflection readout is recovered.

---

## Primary Outputs

| File | Role |
|---|---|
| `曲率付き閉鎖定常波による曲率繰り込みと完全反射安定性.md` | Japanese main paper |
| `curved_closure_stationary_wave_curvature_renormalization_en.md` | English main paper |
| `curved_closure_stationary_wave_curvature_renormalization_ja.tex` | Japanese TeX source |
| `curved_closure_stationary_wave_curvature_renormalization_ja.pdf` | Japanese PDF |
| `curved_closure_stationary_wave_curvature_renormalization_en.tex` | English TeX source |
| `curved_closure_stationary_wave_curvature_renormalization_en.pdf` | English PDF |
| `run_curved_closure_stationary_wave_v1.py` | Minimal closed stationary wave experiment |
| `run_curved_closure_stationary_wave_broad_sweep_v1.py` | Broad curvature phase model sweep |
| `run_curved_closure_scattering_integration_v1.py` | One-sided scattering integration experiment |
| `curved_closure_stationary_wave_result_v1/` | Minimal experiment outputs |
| `curved_closure_stationary_wave_broad_sweep_result_v1/` | Broad sweep outputs |
| `curved_closure_scattering_integration_result_v1/` | Integrated scattering outputs |
| `curved_closure_stationary_wave_zenodo_metadata_v1.json` | Zenodo metadata without token |
| `curved_closure_stationary_wave_zenodo_deposit_v1.json` | Zenodo deposit response without token |
| `curved_closure_stationary_wave_zenodo_published_record_v1.json` | Zenodo published record response without token |
| `curved_closure_stationary_wave_publication_bundle_v1.zip` | Reproducibility bundle |
| `articles/curved-closure-stationary-wave.md` | Zenn article source |
| `note_article_curved_closure_stationary_wave_ja.md` | note Japanese article source |
| `note_article_curved_closure_stationary_wave_en.md` | note English article source |

---

## Numerical Experiment Coverage

1. Minimal closed stationary wave verification.
2. Curvature relative phase leakage and closure-pair residual detection.
3. Common-factor curvature action as a closure-preserving control.
4. Internal phase re-selection by `beta_K,m = -delta_K,m`.
5. Broad sweep over eight curvature relative phase models.
6. Comparison of seven correction freedoms: `none`, `constant`, `linear`, `affine`, `quadratic`, `cubic`, and `full`.
7. One-sided local exchange-interference scattering integration.

---

## Main Findings

- In the minimal experiment, the transient state produced closure-pair RMS `1.2319416790092972e-02` and transmission leakage `1.1503183254481797e-01`.
- After internal phase re-selection, the stationary state recovered closure-pair RMS `9.4283259783636047e-19` and transmission leakage `0.0`.
- In the broad sweep at maximum curvature relative phase `1.2`, the uncorrected case left maximum transmission leakage `1.6202719613622976e-01`.
- The `full` correction recovered closure-pair RMS `7.8949412793793227e-19`, residual phase `0.0`, and transmission leakage `0.0`.
- In one-sided scattering integration, the uncorrected case produced maximum dynamic transmission leakage `1.6202719613622971e-01`.
- The `full` correction reduced maximum dynamic transmission leakage to `1.6608667989341789e-19`.
- Dynamic scattering matched the two-channel expectation with maximum error `5.551115123125783e-17`; maximum norm error was `4.440892098500626e-16`.

---

## Scope

This release does not claim a quantitative prediction of real spacetime curvature, a derivation of standard quantum theory, or a derivation of general relativity.

The confirmed result is internal to the paper series: under the basic axiom `Sigma x_n^2 = 0`, curvature-induced relative phase leakage is detected as closure failure in transient states, and complete-reflection readout is recovered only after the system re-selects an internally closed stationary phase configuration.

---

## Build Notes

- Markdown was converted to standalone TeX with Pandoc.
- PDF generation was performed under `/tmp/curved_closure_stationary_wave_build` to avoid Google Drive filesystem issues during LaTeX compilation.
- XeLaTeX was used for both Japanese and English PDFs.
- Japanese CJK font: `HaranoAjiMincho-Regular.otf`.
- Monospace font: `DejaVuSansMono.ttf`.
- The TeX build succeeded with non-fatal table-width and font-shape warnings.

---

## Publication Status

Publication has been completed.

- Zenodo version DOI: https://doi.org/10.5281/zenodo.21304040
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.21304039
- Zenodo record: https://zenodo.org/records/21304040
- Zenn article source: `articles/curved-closure-stationary-wave.md`
- note Japanese article: https://note.com/kiharanoriaki/n/n2389460836cf
- note English article: https://note.com/kiharanoriaki/n/nda3623c44423
- Facebook Japanese post: https://www.facebook.com/kihara.noriaki/posts/pfbid037wN39hUdgVY7bVCb6BWkFK86pYqDedPgVxaCxtJbWGk479ZqyVoc9XggvRPLmstwl
- Facebook English post: https://www.facebook.com/kihara.noriaki/posts/pfbid02QieNB6aGk3TyRkcreL2z14w62WEFdRbenP47gHuzMoPd4aS7VDPzazbuNJhEGBo4l
- X Japanese post: https://x.com/NoriakiKihara/status/2075801193281106375
- X English post: https://x.com/NoriakiKihara/status/2075802233153998999
