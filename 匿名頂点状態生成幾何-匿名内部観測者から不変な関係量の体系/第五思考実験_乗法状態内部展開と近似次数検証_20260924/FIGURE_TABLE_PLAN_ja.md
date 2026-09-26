# Paper 5 figure/table plan

No physics simulation was rerun for these figures. All plots and tables are derived from the existing audit JSON/CSV files.

## Recommended main-text figures

1. `figure00_hidden_state_audit_structure` — central motivation: apparent two-state model versus audited five-state persistent state.
2. `figure00b_interaction_constraints` — explicit computational rules and why Method B was selected.
3. `figure05_methodB_strict_heatmap` — 12-case strict verification of Method B.
4. `figure06_angular_order_convergence` — angular truncation convergence N=12,16,20,24.
5. `figure07_integration_order_convergence` — numerical order 4,6,8 convergence.
6. `figure08_integration_order_9case_summary` — transfer of order improvement across nine cases.

## Recommended supplementary figures

- `figure02_five_state_reconstruction_errors`
- `figure03_full_inline_equivalence`
- `figure04a_methodA_vs_B_internal`
- `figure04b_methodA_vs_B_observables`
- `figure05b_methodB_savedC1`
- `figure01_experiment_flow` (archive/history figure, not essential for main paper)

## Recommended main-text tables

- `table01_state_audit_summary.csv` — what is state, derived quantity, or candidate multiplicative replacement.
- `table02_approximation_orders.csv` — separate PN order, angular truncation order, integration order, and internal microstep depth.
- `table03_key_numerical_results.csv` — compact numerical evidence used in the conclusions.

## Storage policy

Both SVG and PNG are generated locally. PNG files are not intended for Google Drive upload. SVG, plotting programs, source CSV tables, and this index are the archival versions.
