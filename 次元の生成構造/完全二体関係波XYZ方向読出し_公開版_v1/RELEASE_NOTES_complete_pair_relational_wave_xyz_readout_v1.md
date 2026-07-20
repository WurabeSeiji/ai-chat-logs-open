# Release Notes: 完全二体関係波から読み出されるXYZ三方向 v1

Release date: 2026-07-20

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21454789
- Version DOI: https://doi.org/10.5281/zenodo.21454790
- Zenodo record: https://zenodo.org/records/21454790

## Zenn

- https://zenn.dev/noriaki_kihara/articles/complete-pair-relational-wave-xyz-readout

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Japanese paper: `次元の生成構造/完全二体関係波の反対称生成子から生じる三方向構造_完成論文_v1.md`
- Publication directory: `次元の生成構造/完全二体関係波XYZ方向読出し_公開版_v1/`
- Zenn article source: `articles/complete-pair-relational-wave-xyz-readout.md`

## Publication Status

- A new Zenodo preprint record was created and published as v1.0.
- Japanese and English Markdown, TeX, and PDF files were generated.
- Both PDFs contain the Version DOI and Concept DOI.
- The Japanese PDF is 18 pages; the English PDF is 18 pages.
- The five figures are embedded in both PDFs.
- Eighteen files were uploaded to Zenodo.
- Every uploaded file was validated against the local size and MD5 checksum before publication.
- The Version DOI and Concept DOI both resolve to the published record.

## Main Result

The number of complete pairwise relational waves for $N$ bodies is

$$
M=\binom{N}{2}.
$$

The numerical experiment used 32 trials and 720 steps for each configuration.

- AB: one relational wave, generator rank 0, stationary state.
- ABC: three active relational waves, generator rank 2, kernel dimension 1, one rotation plane and one invariant normal direction.
- ABCD: six active relational waves, generator rank 6, kernel dimension 0, and three rotation planes.

The paper reads the ABC rank-two plane as the phase-readable XY axes and its invariant normal as Z. In ABCD, three rotation planes provide the three uniquely readable XYZ spatial directions; the other three internal relational directions are not uniquely oriented.

Across all configurations, the maximum quadratic-closure error was $1.92\times10^{-13}$, the maximum absolute-square drift was $2.42\times10^{-13}$, and the maximum label-permutation covariance error was $1.47\times10^{-13}$. No stepwise state normalization, observation damping, or absolute background axis was used.

The results for five or more bodies and the physical-axis identities of residual directions are interpretations beyond the directly tested AB, ABC, and ABCD configurations.

## Main Paper Files

- `完全二体関係波の反対称生成子から生じる三方向構造_完成論文_v1.md`
- `complete_pair_relational_wave_xyz_readout_en.md`
- `complete_pair_relational_wave_xyz_readout_ja.tex`
- `complete_pair_relational_wave_xyz_readout_ja.pdf`
- `complete_pair_relational_wave_xyz_readout_en.tex`
- `complete_pair_relational_wave_xyz_readout_en.pdf`

## Reproducibility Files

- `run_ab_abc_abcd_complete_pair_relation_network_preliminary_v1.py`
- `ab_abc_abcd_complete_pair_relation_network_preliminary_result_v1.json`
- `ab_abc_abcd_complete_pair_relation_network_trial_summary_v1.csv`
- `ab_abc_abcd_complete_pair_relation_network_body_summary_v1.csv`
- `ab_abc_abcd_complete_pair_relation_network_selected_series_v1.csv`
- `ab_abc_abcd_complete_pair_relation_network_preliminary_report_v1.md`
- `figures/complete_pair_relation_wave_count_v1.png`
- `figures/ABC_three_physical_relation_waves_v1.png`
- `figures/generator_plane_normal_structure_v1.png`
- `figures/ABC_relation_wave_conservation_v1.png`
- `figures/ABC_one_plane_one_normal_conservation_v1.png`
- `complete_pair_relational_wave_xyz_readout_source_v1.zip`

## Zenodo Records

- `complete_pair_relational_wave_xyz_readout_zenodo_metadata_v1.json`
- `complete_pair_relational_wave_xyz_readout_zenodo_deposit_v1.json`
- `complete_pair_relational_wave_xyz_readout_zenodo_before_publish_v1.json`
- `complete_pair_relational_wave_xyz_readout_zenodo_upload_manifest_v1.tsv`
- `complete_pair_relational_wave_xyz_readout_zenodo_uploads_v1.json`
- `complete_pair_relational_wave_xyz_readout_zenodo_publish_response_v1.json`
- `complete_pair_relational_wave_xyz_readout_zenodo_published_record_v1.json`

## Build Notes

- Markdown-to-TeX conversion and PDF compilation were performed in `/private/tmp`, not in the Google Drive working tree.
- Pandoc standard Markdown input was used so that standalone equality signs inside display mathematics were not misread as setext headings.
- XeLaTeX was run twice for each language.
- The wide claim-classification table was placed on a landscape page.
