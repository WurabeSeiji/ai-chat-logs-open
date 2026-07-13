# Release Notes: Fermionic Interference Reflection Engine v1/v2

**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Repository:** `ai-chat-logs-open`  
**Latest Zenodo Version DOI:** 10.5281/zenodo.21332867  
**Zenodo Concept DOI:** 10.5281/zenodo.21295479  
**Latest Zenodo Record:** https://zenodo.org/records/21332867  
**Previous v1 Version DOI:** 10.5281/zenodo.21295480  
**Zenn Article:** `articles/fermionic-interference-reflection-engine.md`  

---

## Summary

This release prepares the additional constructive experiment paper:

> Interference Construction of a Perfect Reflection Map from a Fermionic Inverse-Phase Core

The paper replaces the direct direction-reversal rule used in the preceding AB/C elastic-reflection simulation with a local exchange-interference map controlled by the internal phase of a fermionic inverse-phase core.

The preceding paper established the full AB/C closed-phase-system architecture. This release supplies the interaction engine: a local map in which direct and exchange paths are decomposed into even and odd channels, the internal core phase `Delta_F` is transferred to the odd channel, and the recombined waveform produces reflection and transmission readouts.

---

## v2 Update

V2 recalculates the experiment with the reflection and transmission amplitudes acting as a two-channel scattering matrix on the incident channels.

The Concept DOI is maintained. The V2 Version DOI is `10.5281/zenodo.21332867`.

The Japanese/English Markdown, TeX/PDF, numerical outputs, and reproducibility bundle were updated to the V2 outputs.

---

## Primary Outputs

| File | Role |
|---|---|
| `フェルミオン的逆相核による完全反射写像の干渉構成 v2.md` | Japanese main paper |
| `fermionic_inverse_phase_core_interference_reflection_en.md` | English main paper |
| `fermionic_interference_reflection_engine_ja.tex` | Japanese TeX source |
| `fermionic_interference_reflection_engine_ja.pdf` | Japanese PDF |
| `fermionic_interference_reflection_engine_en.tex` | English TeX source |
| `fermionic_interference_reflection_engine_en.pdf` | English PDF |
| `run_fermionic_interference_reflection_v2.py` | Executed numerical experiment script |
| `fermionic_interference_reflection_result_v2/` | JSON, CSV, plots, and execution report |
| `fermionic_interference_reflection_engine_zenodo_deposit_v2.json` | Zenodo metadata draft without token |
| `fermionic_interference_reflection_engine_publication_bundle_v2.zip` | Reproducibility bundle |
| `articles/fermionic-interference-reflection-engine.md` | Zenn article draft |

---

## Numerical Experiment Coverage

The publication bundle covers the following checks:

1. One-sided incident packet scattering without a mirror-image initial condition.
2. Local-map phase sweep over `Delta_F`.
3. Exchange-interference node formation at complete overlap.
4. Removal of the exchange path as a control.
5. Auxiliary odd-node half-line readout.
6. Identification oscillation readout for `m_A=1` and `m_B=2`.
7. Reversibility checks for `U(pi)^2` and `U(delta)U(-delta)`.
8. Compensated square-closure checks using `x_n, i x_n` pairs.
9. AB/C cell replacement test, replacing direct `q` reversal with `q_out=q_in*(T-R)`.

---

## Main Findings

- A single incident packet placed only on the left side produced full transmission for `Delta_F=0`, half splitting for `Delta_F=pi/2`, and full reflection for `Delta_F=pi`.
- The initial right-side probability was `2.4303500961591473e-89`, so the reflection is not caused by an embedded mirror component.
- The phase sweep reproduced `R(Delta_F)=sin^2(Delta_F/2)` and `T(Delta_F)=cos^2(Delta_F/2)` with maximum error `5.551115123125783e-16`.
- The maximum norm error was `6.661338147750939e-16`.
- The exchange-interference node at `Delta_F=pi` reached diagonal relative norm `7.4987989133092880e-33`.
- Removing the exchange path destroyed the node, leaving `4.9999999999999994e-01`.
- `U(pi)^2` returned to the initial waveform with relative error `4.8214412843789535e-11`, within the verdict threshold `1e-10`.
- `U(delta)U(-delta)` remained reversible to maximum relative error `1.4547631339456792e-16`.
- The maximum compensated square-closure residual was `1.2143074258005e-17`.
- The AB/C replacement test generated `q_A: 1 -> -1` and `q_B: -1 -> 1` without a direct `q=-q` instruction, while preserving identification modes.

---

## Scope

This release does not claim to derive standard fermion scattering, S-matrix scattering theory, or the standard quantum measurement process.

The confirmed result is that the direct direction-reversal rule of the preceding simulation can be replaced by a conservative, reversible, phase-selective local exchange-interference map. The remaining lower-level constructive question is why the internal core phase is transferred to the exchange-path phase; in this release, that transfer is the local coupling rule being tested.

---

## Build Notes

- Markdown was converted to TeX with Pandoc.
- PDF generation was performed under `/tmp/fermionic_interference_reflection_build` to avoid Google Drive filesystem issues during LaTeX image embedding.
- XeLaTeX was used for both Japanese and English PDFs.
- Japanese CJK font: `HaranoAjiMincho-Regular.otf`.
- Monospace font: `DejaVuSansMono.ttf`.
- The TeX build succeeded with non-fatal table-width and font-shape warnings.

---

## Publication Status

Publication has been completed.

- Zenodo version DOI: https://doi.org/10.5281/zenodo.21332867
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.21295479
- Zenodo record: https://zenodo.org/records/21332867
- GitHub final commit: recorded in the v2 publication commit
- Zenn article source: `articles/fermionic-interference-reflection-engine.md`
