# Release Notes: Discrete Born-Type Weights from Finite-Order Recurrence v1.0

Release date: 2026-07-18

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21422470
- Version DOI: https://doi.org/10.5281/zenodo.21422471
- Zenodo record: https://zenodo.org/records/21422471

## Zenn

- https://zenn.dev/noriaki_kihara/articles/finite-order-recurrence-born-type-weights

## note

- Japanese: https://note.com/kiharanoriaki/n/n427655620ae0
- English: https://note.com/kiharanoriaki/n/nb4340afa4a9e

## Facebook

- Japanese: https://www.facebook.com/kihara.noriaki/posts/pfbid02xTby2au4pp2z2urarBygGQvuXCbur2qJRbRGVQJfaeyXFN8BcHFHouEQqB6MyJwGl

## X

- Japanese: https://x.com/NoriakiKihara/status/2078344339659669578
- English: https://x.com/NoriakiKihara/status/2078344911926292864

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Japanese manuscript: `波の情報読出し/20260718/finite_order_recurrence_discrete_born_type_weights_ja.md`
- English manuscript: `波の情報読出し/20260718/finite_order_recurrence_discrete_born_type_weights_en.md`
- Zenn article source: `articles/finite-order-recurrence-born-type-weights.md`
- note article sources: `波の情報読出し/20260718/note_article_quadratic_closure_finite_recurrence_born_probability_wave_{ja,en}.md`
- Facebook post sources: `波の情報読出し/20260718/fb_quadratic_closure_finite_recurrence_born_probability_wave_{ja,en}.md`
- X post sources: `波の情報読出し/20260718/x_quadratic_closure_finite_recurrence_born_probability_wave_{ja,en}.md`

## Main Result

For an exchange-symmetric two-channel unitary operator, fixing the global phase gives

```text
U = P_s + zeta P_a.
```

Exact finite-order recurrence `U^n = I` selects `zeta = exp(-2 pi i m/n)`. Returning to the A/B basis gives

```text
r = (1 + zeta)/2,
t = (1 - zeta)/2,
|r|^2 = cos^2(pi m/n),
|t|^2 = sin^2(pi m/n).
```

The discrete Born-type weights therefore arise from exchange symmetry, cyclotomic finite-order closure, and channel projection. The trigonometric scattering parametrization is recovered as a coordinate representation rather than imposed as an independent assumption.

## Scope and Limits

- The paper derives discrete Born-type transition weights, not the complete Born rule.
- It does not derive single-trial randomness or repeated-measurement frequencies.
- The theorem concerns the common linear exchange kernel.
- System A's channel-wise normalization and System B's weak-readout and strong-observation maps are separated as downstream operations.
- The previously observed peaks `R_(124,23)` and `R_(122,23)` are identified as even-order cyclotomic recurrence roots.

## Publication Files

- Japanese and English Markdown
- Japanese and English TeX
- Japanese and English PDF
- Zenodo metadata and publication audit files
- Zenn explanatory article
- Japanese and English note articles
- Japanese and English Facebook and X post sources

## Build Notes

- TeX and PDF generation was performed under `/tmp`, outside the Google Drive working tree.
- Japanese PDF: XeLaTeX, 28 A4 pages.
- English PDF: pdfLaTeX, 28 A4 pages.
- Markdown table delimiters inside absolute-value formulas were replaced by `\lvert` and `\rvert` to preserve the intended table structure.
- The claim-classification tables use fixed-width TeX columns to prevent overflow.
