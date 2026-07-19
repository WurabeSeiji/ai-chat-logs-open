# Release Notes: AB二体閉鎖位相系の調和閉鎖による逆二乗則 v1

Release date: 2026-07-19

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21441081
- Version DOI: https://doi.org/10.5281/zenodo.21441082
- Zenodo record: https://zenodo.org/records/21441082

## Zenn

- https://zenn.dev/noriaki_kihara/articles/inverse-square-phase-harmonic-closure

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Paper source directory: 次元の生成構造/加速度逆二乗則機構_追加論文_v1/
- Zenn article source: articles/inverse-square-phase-harmonic-closure.md

## Publication Status

- A new Zenodo preprint record was created and published as v1.0.
- Japanese and English Markdown, TeX, and PDF files were generated.
- Both PDFs contain the Version DOI and Concept DOI.
- The Japanese PDF is 22 pages; the English PDF is 18 pages.
- The two figures are embedded in both PDFs.
- The eight-column numerical summary table is placed on a landscape page in both PDFs.
- Fourteen files were uploaded to Zenodo.
- Every uploaded file was validated against the local size and MD5 checksum before publication.
- The Version DOI and Concept DOI both resolve to the published record.

## Main Result

The published two-body AB experiment already contained the nonzero second-order relation

$$
\Delta^2\chi_s=-\omega_d^2\chi_s.
$$

The present paper connects the future-phase-position acceleration map

$$
\alpha_n=R|\omega_n|^2
$$

to the harmonic closure

$$
|\omega_n|\Delta\theta_n=\Omega.
$$

It follows that

$$
\alpha_n
=\frac{R\Omega^2}{\Delta\theta_n^2}.
$$

Under the specific conditions of the reported experiment, the inverse-square law was established.

Generalization to arbitrary closed systems, harmonic arrangements, nonharmonic updates, and distance mappings remains untested. The paper does not identify the acceleration-like readout with standard gravity.

## Main Paper Files

- AB二体閉鎖位相系における未来位相位置加速度写像と調和閉鎖による逆二乗則_完成論文_v1.md
- inverse_square_law_future_phase_position_acceleration_harmonic_closure_en.md
- inverse_square_phase_closure_ja.tex
- inverse_square_phase_closure_ja.pdf
- inverse_square_phase_closure_en.tex
- inverse_square_phase_closure_en.pdf

## Reproducibility Files

- figures/未来位置中心回転写像と調和位相逆二乗機構_v1.png
- figures/未来位置中心回転写像と調和位相逆二乗機構_v1.svg
- figures/既存AB加速度時系列と二階差分整合_v1.png
- figures/既存AB加速度時系列と二階差分整合_v1.svg
- tables/既存AB加速度発生集計_v1.csv
- tables/既存AB加速度発生集計_v1.md
- make_existing_ab_acceleration_evidence_v1.py
- inverse_square_phase_closure_source_v1.zip

## Zenodo Records

- inverse_square_phase_closure_zenodo_metadata_v1.json
- inverse_square_phase_closure_zenodo_deposit_v1.json
- inverse_square_phase_closure_zenodo_before_publish_v1.json
- inverse_square_phase_closure_zenodo_upload_manifest_v1.tsv
- inverse_square_phase_closure_zenodo_uploads_v1.json
- inverse_square_phase_closure_zenodo_publish_response_v1.json
- inverse_square_phase_closure_zenodo_published_record_v1.json

## Build Notes

- TeX conversion and PDF compilation were performed under a local /private/tmp directory rather than in the Google Drive working tree.
- A writable task-specific LuaTeX font cache was used under the same temporary directory.
- LuaLaTeX was run twice for each language.
- The wide eight-column numerical table was placed in a landscape environment without changing its values.
