# Release Notes: 交換重み G_R = 1 - R の選別と微細構造定数対応候補 v1

Release date: 2026-07-16

## DOI

- Concept DOI: https://doi.org/10.5281/zenodo.21396760
- Version DOI: https://doi.org/10.5281/zenodo.21396761
- Zenodo record: https://zenodo.org/records/21396761

## Zenn

- https://zenn.dev/noriaki_kihara/articles/exchange-weight-alpha-correspondence

## GitHub

- Repository: https://github.com/WurabeSeiji/ai-chat-logs-open
- Release commit: `3623fef Publish exchange weight alpha correspondence paper`
- Zenn article source: `articles/exchange-weight-alpha-correspondence.md`

## Publication Status

- Concept DOI was maintained for this research line.
- Version DOI was assigned and published on Zenodo.
- Japanese and English Markdown, TeX, and PDF files were generated.
- The Zenodo upload manifest and upload responses were saved in the release directory.
- Zenn article source was added and pushed to GitHub.

## Main Paper

- `交換散乱係数R集中と微細構造定数対応候補の数値実験 v1.md`
- `exchange_weight_alpha_correspondence_numerical_experiment_ja.md`
- `exchange_weight_alpha_correspondence_numerical_experiment_en.md`
- `exchange_weight_alpha_correspondence_numerical_experiment_ja.tex`
- `exchange_weight_alpha_correspondence_numerical_experiment_en.tex`
- `exchange_weight_alpha_correspondence_numerical_experiment_ja.pdf`
- `exchange_weight_alpha_correspondence_numerical_experiment_en.pdf`

## Build Notes

- PDF generation was performed outside the Google Drive working tree under `/tmp/tex_compile_exchange_weight_alpha`.
- LuaLaTeX failed in this environment because `luaotfload` could not create a writable font cache.
- The released PDFs were therefore generated with XeLaTeX.
- Generated PDFs:
  - Japanese PDF: 17 pages
  - English PDF: 18 pages

## Uploaded Data

- `exchange_weight_alpha_correspondence_zenodo_metadata_v1.json`
- `exchange_weight_alpha_correspondence_zenodo_upload_manifest_v1.tsv`
- `exchange_weight_alpha_correspondence_zenodo_uploads_v1.json`
- `exchange_weight_alpha_correspondence_zenodo_extra_uploads_v1.json`
- `exchange_weight_alpha_correspondence_zenodo_publish_response_v1.json`
- `exchange_weight_alpha_correspondence_zenodo_published_record_v1.json`
- `交換散乱係数R集中実験群_全体計画と評価方法_v1.md`
- `系統A_局在性交換R近傍斉一スイープ実験仕様書_v1.md`
- `系統B_灰色猫準安定界面R近傍スイープ実験仕様書_v1.md`
- `System B 全域Rスイープ候補帯一覧.md`
- `System B 低エネルギー標準α近傍R感度一覧.md`
- `System B 高エネルギー標準α近傍R感度一覧.md`
- `System B 中間R近傍R感度一覧.md`
- `System B 極低エネルギー標準α近傍R感度一覧.md`
- `標準理論でのα.md`
- `run_system_A_localization_exchange_R_sweep_preliminary_v1.py`
- `run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py`
- `run_minimal_system_B_gray_direct_check_v5.py`

## Zenodo Upload Summary

- Main upload: 36 files uploaded successfully.
- Extra manifest upload: 3 files uploaded successfully.
- Large intermediate CSV files were not included in the Zenodo upload bundle.

## Uploaded Figures

- `system_A_localization_exchange_R_sweep_result_v1/odd_kernel_N_1_2_3_5_15_63_gap_depth_distribution_overview_v1.png`
- `system_B_gray_cat_metastable_R_sweep_result_v1/system_B_odd_kernel_N_1_2_3_5_15_63_gray_depth_distribution_overview_v1.png`
- `system_B_full_R_sweep_ranked_depth_bar_v1.svg`
- `system_B_full_R_sweep_full_range_depth_v1.svg`
- `system_B_full_R_sweep_band_detail_rank_01_v1.svg`
- `system_B_full_R_sweep_band_detail_rank_02_v1.svg`
- `system_B_low_alpha_R_sensitivity_depth_v1.svg`
- `system_B_high_alpha_R_sensitivity_depth_v1.svg`
- `system_B_mid_R_sensitivity_depth_v1.svg`

## Summary

This release publishes a numerical experiment on the selection of the state-exchange weight `G_R = 1 - R` in closed exchange maps.

The study compares two independent experimental systems: one based on localization exchange in an AB harmonic packet model, and one based on gray-cat metastable selection in an AB allocation model.

The main result is that both systems indicate a stable exchange-weight neighborhood corresponding to the low-energy fine-structure constant within the effective numerical precision, while the broader sweep also exposes a high-energy-neighborhood candidate and weaker auxiliary bands.

The result is not claimed as a first-principles derivation of the fine-structure constant. It is positioned as evidence that a closed state-exchange map can select an intrinsic dimensionless exchange weight.
