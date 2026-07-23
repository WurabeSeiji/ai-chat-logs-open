# Release Notes: Iterated Exchange-Scattering Finite-Order Resonance v1.0

Release date: 2026-07-18

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21421366
- Version DOI: https://doi.org/10.5281/zenodo.21421367
- Zenodo record: https://zenodo.org/records/21421367

## Zenn

- https://zenn.dev/noriaki_kihara/articles/finite-order-resonance-alpha-neighborhood

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Japanese manuscript: `波の情報読出し/20260718/反復交換散乱における有限位数共鳴の発見_微細構造定数137・128近傍ピークの原因特定と再現可能な波束数理モデル_ja.md`
- English manuscript: `波の情報読出し/20260718/finite_order_resonance_iterated_exchange_scattering_en.md`
- Zenn article source: `articles/finite-order-resonance-alpha-neighborhood.md`

## Main Result

The sharp peaks previously observed near exchange-weight readouts corresponding to inverse fine-structure-constant values 137 and 128 were identified as finite-order recurrence roots of the two-channel unitary exchange operator.

The exact roots are

```text
R_(n,m) = cos^2(pi m/n).
```

The principal low-energy peak is `R_(124,23)`, the former high-energy candidate is `R_(122,23)`, and the even root `R_(620,117)` gives `N(R)=128.947864735670559`, which differs from `alpha^(-1)(M_Z^2)=128.946 +/- 0.015` by `0.1243 sigma` under the paper's stated correspondence.

The paper does not claim a derivation of `alpha`. It identifies a different exact recurrence phenomenon and leaves as the central open problem why extremely sharp finite-order roots occur near the two physical `alpha^(-1)` values.

## Reproducibility

- The uniform full-range sweep contains no `alpha`, 137, 128, or analytic root in its evaluation function.
- The source code implementing the complex scattering amplitudes, iterative update, gray error, analytic eigenvalues, and multiprecision recurrence is preserved.
- Root-centered sweeps at `1e-10` and `1e-12`, together with 50- and 80-digit calculations, reproduce the analytic roots.
- Ideal even fundamental-order roots have zero gray error and formally infinite depth; this denotes vanishing residual, not divergent energy.

## Publication Files

- Japanese and English Markdown
- Japanese and English TeX
- Japanese and English PDF
- Five principal figures
- Reproduction bundle containing source code, summary data, root-centered datasets, and figures

## Build Notes

- TeX and PDF generation was performed under `/tmp`, outside the Google Drive working tree.
- Japanese PDF: XeLaTeX, 23 A4 pages.
- English PDF: pdfLaTeX, 23 A4 pages.
- All five figures are embedded in both PDFs.
- LuaLaTeX was not used because the execution environment rejected its font-cache path; XeLaTeX produced the Japanese release PDF without missing glyphs.

