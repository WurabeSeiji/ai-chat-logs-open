# Release Notes: Elastic Reflection in a Closed Phase System v1

**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Repository:** `ai-chat-logs-open`  
**Zenodo Version DOI:** 10.5281/zenodo.21291020  
**Zenodo Concept DOI:** 10.5281/zenodo.21291018  
**Zenodo Record:** https://zenodo.org/records/21291020  
**Zenn Article:** `articles/elastic-reflection-closed-phase-system.md`

---

## Summary

This release publishes the constructive experiment paper:

> Constructive Experiment on Elastic Reflection of Two Fermionic Local Waves in a Closed Phase System Without Assuming Background Space

The work constructs and numerically tests a finite-resolution conservative reflection map for two identifiable fermionic local waves in a closed phase system, without assuming an external background space as the starting point.

The release includes the Japanese and English main papers, supporting axiom/specification/result documents, generated TeX/PDF files, figures, JSON/CSV result files, and simulation scripts.

---

## Primary Outputs

| File | Role |
|---|---|
| `背景空間を仮定しない閉じた位相系におけるフェルミオン的二局所波の完全弾性反射の構成実験.md` | Japanese main paper |
| `constructive_experiment_elastic_reflection_closed_phase_system_en.md` | English main paper |
| `elastic_reflection_closed_phase_system_ja.tex` | Japanese TeX source |
| `elastic_reflection_closed_phase_system_ja.pdf` | Japanese PDF |
| `elastic_reflection_closed_phase_system_en.tex` | English TeX source |
| `elastic_reflection_closed_phase_system_en.pdf` | English PDF |
| `basic_axiom_system_v3_en.md` / `基本公理系.md` | Axiom-system references |
| `elastic_collision_simulation_spec_v1_en.md` / `完全弾性衝突シミュレーション仕様書 v1.md` | Simulation specification |
| `elastic_collision_simulation_experiment_results_v1_en.md` / `完全弾性衝突シミュレーション実験結果 v1.md` | Integrated result report |
| `elastic_reflection_closed_phase_system_publication_bundle_v1.zip` | Reproducibility bundle |

---

## Numerical Experiment Coverage

The publication bundle includes the following experiment groups:

1. Minimal complete elastic collision.
2. Identification-mode robustness.
3. Observer-capacity sweep.
4. Cell-resolution and update-step sweep.
5. Reflection/transmission/label-exchange control maps.
6. Asymmetric-condition tests.
7. Observation-perturbation tests.
8. Repeated-collision tests.
9. Internal identification-phase `η` readout-resolution tests.

---

## Main Findings

- The reflection map preserves identification modes, representative amplitudes, fermionic cores, and compensated square closure.
- Reflection is distinguishable from transmission and label-exchange maps by jointly reading direction reversal and internal identification modes.
- Repeated collisions preserve the same invariants over eight AB collisions.
- Observer insufficiency, cell overshooting, temporal-cell mismatch, identification-mode leakage, observation perturbation, and `η` aliasing appear as distinct failure conditions.

---

## Scope

This release does not claim to derive standard fermion scattering, S-matrix scattering theory, or the standard quantum measurement process.

The direction-reversal rule is introduced as a finite-resolution interaction-cell construction rule. The result shown here is that this rule is compatible with the closed phase-system construction, internal identification modes, internal observation, and compensated square closure, and that it can be distinguished from control maps.

---

## Build Notes

- Markdown was converted to TeX with Pandoc.
- PDF generation was performed under `/tmp` to avoid Google Drive filesystem issues during LaTeX image embedding.
- XeLaTeX was used for both Japanese and English PDFs because the local LuaLaTeX environment failed on `luaotfload` cache initialization in this session.
- Japanese CJK font: `HaranoAjiMincho-Regular.otf`.
- Monospace font: `DejaVuSansMono.ttf`, selected to preserve Greek symbols inside code spans.
