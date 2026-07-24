# Release Notes: 無名等振幅複合波モデル基本公理系 v9

Release date: 2026-07-24

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21315735
- Latest Version DOI (v9.1): https://doi.org/10.5281/zenodo.21522310
- Initial v9 Version DOI: https://doi.org/10.5281/zenodo.21522219
- Latest Zenodo record: https://zenodo.org/records/21522310

## Zenn

- https://zenn.dev/noriaki_kihara/articles/basic-axiom-system-v9

## Publication Status

- The existing Concept DOI was preserved through the Zenodo `newversion` action.
- The Japanese and English Markdown sources contain the new Version DOI.
- The Japanese and English TeX and PDF files were generated from the DOI-fixed Markdown sources.
- The Japanese PDF has 27 pages; the English PDF has 21 pages.
- Six v9 files were uploaded after removing the six inherited v8 files from the draft.
- Every uploaded file was validated against its local size and MD5 checksum before publication.
- The initial v9 publication and the corrected v9.1 publication both returned HTTP 202.
- Both the Version DOI and Concept DOI resolve successfully.

## v9.1 Formatting Correction

The initial v9 archive rendered Markdown hard-line-break markers as visible backslash characters in the title metadata block of the generated TeX and PDF files.

Version v9.1 corrects those markers to actual line breaks. The axiom text, numbering, mathematics, and claims are unchanged. The Concept DOI remains unchanged, and v9.1 is its latest version.

## Main Change

Version 9 adds the following axiom between Axiom 0.6 and Axiom 1:

```text
Axiom 0.7: Anonymous Partial-Projection Existence
```

Let the scale-quotiented candidate space satisfying the adopted mandatory axioms other than Axiom 0.7 be

$$
\mathscr S_{\mathrm{ax}}.
$$

The axiom states that the physically realized state space is the image of at least one nontrivial admissible partial projection:

$$
\boxed{
\exists\mathcal P_*
\in
\mathfrak P(\mathscr S_{\mathrm{ax}})
\quad\text{such that}\quad
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*.
}
$$

Admissibility requires independence from component and formulation names, independence from absolute scale, consistency with foundational closure and image composition, readout consistency, and prohibition against using a post-projection physical or symmetry-group name as a projection condition.

## Numbering Preservation

- Axioms 0, 0.5, 0.6, and 1–17 retain their existing numbers.
- Their texts, together with the working axioms, are unchanged from v8.
- No previously published citation number is shifted by the v9 insertion.

## Scope

Axiom 0.7 does not assert:

```text
uniqueness of the admissible partial projection
an already completed selector functional or dynamics
unique enforcement of a particular symmetry group by the foundational axioms
```

No specific symmetry-group name, eight-component assumption, two-quartet structure, or derived number 240, 248, or 30 is placed in the axiom.

The conditional identification of the $E_8$ lattice under a positive-definite eight-component readout, two-quartet decomposition, local $D_4$ transition lattices, and common-center neutrality projection remains in the independent companion paper:

- Concept DOI: https://doi.org/10.5281/zenodo.21521899
- Version DOI: https://doi.org/10.5281/zenodo.21521900

## Published Files

- `無名等振幅複合波モデル基本公理系v9_純化定義論文.md`
- `anonymous_equal_amplitude_composite_wave_model_basic_axiom_system_v9_pure_definition_en.md`
- `basic_axiom_system_v9_pure_definition_ja.tex`
- `basic_axiom_system_v9_pure_definition_ja.pdf`
- `basic_axiom_system_v9_pure_definition_en.tex`
- `basic_axiom_system_v9_pure_definition_en.pdf`

## Build Validation

- TeX conversion and PDF compilation were performed in `/private/tmp`.
- Both PDFs were compiled twice with LuaLaTeX.
- The English monospaced font was explicitly selected to preserve the Greek lambda glyph in code blocks.
- Final logs contain no LaTeX errors, undefined control sequences, missing characters, or overfull boxes.
- PDF text extraction confirms the Version DOI and Concept DOI in both language editions.
