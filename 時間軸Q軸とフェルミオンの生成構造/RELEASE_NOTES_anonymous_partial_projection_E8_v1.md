# Release Notes: 無名部分投射による物理的対称性選択とE8分枝 v1

Release date: 2026-07-24

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21521899
- Version DOI: https://doi.org/10.5281/zenodo.21521900
- Zenodo record: https://zenodo.org/records/21521900

## Zenn

- https://zenn.dev/noriaki_kihara/articles/anonymous-partial-projection-e8

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Paper source directory: `時間軸Q軸とフェルミオンの生成構造/`
- Zenn article source: `articles/anonymous-partial-projection-e8.md`

## Publication Status

- A new Zenodo preprint record was created and published as v1.0.
- The Japanese and English Markdown, TeX, and PDF files were generated.
- Both PDFs contain the Version DOI and Concept DOI.
- The Japanese PDF has 23 pages; the English PDF has 17 pages.
- Six files were uploaded to Zenodo.
- Every uploaded file was validated against its local size and MD5 checksum before publication.
- The published record reports Version DOI `10.5281/zenodo.21521900` and Concept DOI `10.5281/zenodo.21521899`.

## Main Result

The paper introduces the anonymous partial-projection existence axiom

$$
\mathscr S_{\mathrm{phys}}
=
\operatorname{Im}\mathcal P_*,
\qquad
\mathcal P_*
\in
\mathfrak P(\mathscr S).
$$

It does not place a symmetry-group name in the foundational axioms. Under the explicit branch assumptions of a positive-definite rank-eight readout, a two-quartet decomposition, local $D_4$ transition lattices, and a common-center neutrality projection, only

$$
(0,0),\qquad(v,v),\qquad(s,s),\qquad(c,c)
$$

survive. Their diagonal gluing gives a positive-definite rank-eight even unimodular lattice and therefore

$$
\Lambda\cong\Lambda_{E_8}.
$$

The derived consequences are 240 roots, $\dim\mathfrak e_8=248$, and the existence of a Coxeter element satisfying $C^{30}=I$.

The result is conditional. The foundational axioms do not uniquely force $E_8$; the selector responsible for the eight-component, local-$D_4$, common-center branch remains an open problem.

## Published Paper Files

- `無名部分投射による物理的対称性選択の公理的枠組み_E8分枝の条件付き同定_独立論文_v1.md`
- `anonymous_partial_projection_physical_symmetry_selection_E8_en.md`
- `anonymous_partial_projection_physical_symmetry_selection_E8_ja.tex`
- `anonymous_partial_projection_physical_symmetry_selection_E8_ja.pdf`
- `anonymous_partial_projection_physical_symmetry_selection_E8_en.tex`
- `anonymous_partial_projection_physical_symmetry_selection_E8_en.pdf`

## Zenodo Records

- `anonymous_partial_projection_E8_zenodo_metadata_v1.json`
- `anonymous_partial_projection_E8_zenodo_deposit_v1.json`
- `anonymous_partial_projection_E8_zenodo_before_publish_v1.json`
- `anonymous_partial_projection_E8_zenodo_upload_manifest_v1.tsv`
- `anonymous_partial_projection_E8_zenodo_uploads_v1.json`
- `anonymous_partial_projection_E8_zenodo_publish_response_v1.json`
- `anonymous_partial_projection_E8_zenodo_published_record_v1.json`

## Build Notes

- The TeX conversion and PDF compilation were performed in a local `/private/tmp` directory.
- The Japanese PDF was compiled twice with XeLaTeX and explicit Japanese fonts.
- The English PDF was compiled twice with pdfLaTeX.
- Final logs contain no LaTeX errors, undefined control sequences, overfull boxes, or missing characters.
