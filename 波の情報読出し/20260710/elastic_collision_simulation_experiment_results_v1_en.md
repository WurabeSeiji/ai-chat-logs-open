# Complete Elastic Collision Simulation Experiment Results v1

**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Status:** Integrated experiment results for the constructive experiment paper  

---

## 0. Conclusion

The complete elastic reflection map for two distinguishable fermionic local waves in a closed phase system was numerically tested through nine experiment families.

The basic experiment confirmed:

- arrival at the finite-resolution collision cell;
- reversal of direction readouts `q_A,q_B`;
- preservation of identification oscillation modes `m_A,m_B`;
- preservation of representative amplitudes;
- no label exchange;
- no identification-mode crosstalk;
- preservation of compensated square closure.

Additional sweeps separated validity conditions and failure modes:

- the observer `C` must be sufficiently heavy;
- the computational step must not skip the finite cell;
- spatial and temporal cells must be satisfied simultaneously;
- identification-mode leakage has a finite tolerance;
- finite sampling of the identification phase `η` can produce aliasing;
- reflection is distinguishable from transmission and label exchange;
- repeated collisions preserve the same structure.

These results do not derive standard scattering theory. They show that the proposed finite-resolution reflection map is compatible with localization, internal identification, internal observation, and compensated square closure in the working axiom system.

---

## 1. Experiment List

| No. | Experiment | Directory |
|---:|---|---|
| 1 | Basic complete elastic collision | `elastic_collision_simulation_result_v1/` |
| 2 | Identification oscillation preservation and robustness | `elastic_collision_label_robustness_result_v1/` |
| 3 | Heaviness condition of observer C | `elastic_collision_observer_sweep_result_v1/` |
| 4 | Cell resolution and time step | `elastic_collision_cell_resolution_sweep_result_v1/` |
| 5 | Reflection/transmission/label-exchange controls | `elastic_collision_control_maps_result_v1/` |
| 6 | Asymmetric conditions | `elastic_collision_asymmetry_sweep_result_v1/` |
| 7 | Observation perturbation | `elastic_collision_observation_perturbation_result_v1/` |
| 8 | Multiple collisions | `elastic_collision_multi_collision_result_v1/` |
| 9 | Readout resolution of `η` identification oscillation | `elastic_collision_eta_resolution_sweep_result_v1/` |

---

## 2. Basic Complete Elastic Collision

### 2.1 Conditions

The minimal symmetric experiment used:

| Quantity | Value |
|---|---:|
| `A_A,A_B` | `1,1` |
| `A_C` | `1000` |
| `N_{h,\chi}^A,N_{h,\chi}^B` | `99,99` |
| `N_{h,\tau}^A,N_{h,\tau}^B` | `99,99` |
| `m_A,m_B` | `1,2` |
| `q_A,q_B` | `+1,-1` |
| `delta_s` | `0.01` |

### 2.2 Verdict

| Item | Result |
|---|---|
| Collision cell reached | `true` |
| Collision step | `19` |
| Final step | `40` |
| Direction reversal | `true` |
| Identification mode preserved A | `true` |
| Identification mode preserved B | `true` |
| Label mode swapped | `false` |
| Label mode crosstalk | `false` |
| Representative amplitude preserved | `true` |
| Closure residual absolute value | `0.0` |
| Elastic collision map valid | `true` |

### 2.3 Figures

![Basic collision trajectory](elastic_collision_simulation_result_v1/elastic_collision_trajectory_v1.png)

![Identification mode readout](elastic_collision_simulation_result_v1/elastic_collision_identification_modes_v1.png)

---

## 3. Identification Oscillation Preservation and Robustness

The purpose was to test whether the internal identification modes remain readable when leakage is added to the identification oscillation.

### 3.1 Result

| Item | Value |
|---|---:|
| Total cases | `36` |
| Valid cases | `20` |
| Invalid cases | `16` |
| First failing leakage rate | `0.35` |

The map remains valid when the dominant identification mode is still separable. When leakage becomes large enough to break identification purity, the identification mode can no longer be treated as preserved.

![Identification mode purity](elastic_collision_label_robustness_result_v1/label_robustness_purity_v1.png)

![Identification mode verdict](elastic_collision_label_robustness_result_v1/label_robustness_detection_v1.png)

---

## 4. Heaviness Condition of Observer C

The purpose was to test how large the observer representative amplitude `A_C` must be for `C` to be treated as quasi-static.

### 4.1 Result

| Item | Value |
|---|---:|
| Total cases | `13` |
| Valid cases | `7` |
| Invalid cases | `6` |
| First valid `A_C` | `20` |

For `A_C < 20`, the heavy-observer condition failed. For `A_C >= 20`, the simulation satisfied the quasi-static observer condition under the tested parameters.

![Observer sweep conditions](elastic_collision_observer_sweep_result_v1/observer_sweep_conditions_v1.png)

![Observer sweep validity](elastic_collision_observer_sweep_result_v1/observer_sweep_validity_v1.png)

---

## 5. Cell Resolution and Time Step

The purpose was to test the relation between the update step and the finite-resolution cell width.

### 5.1 Result

| Item | Value |
|---|---:|
| Total cases | `60` |
| Valid cases | `57` |
| Invalid cases | `3` |
| Off-grid valid cases | `27` |
| Off-grid invalid cases | `3` |

The robust practical condition is:

```math
\Delta s \le \epsilon_\chi^{AB}.
```

When the update step is too large relative to the cell width, the local waves can pass over the cell without satisfying the cell condition at a sampled step.

![Cell resolution sampling condition](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sampling_condition_v1.png)

![Off-grid validity](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_validity_d0_0_203_v1.png)

---

## 6. Reflection, Transmission, and Label-Exchange Control Experiments

The purpose was to check whether the successful result is specific to reflection, or could also be produced by transmission or label exchange.

### 6.1 Result

| Map | Direction reversal | Identification preservation | Verdict |
|---|---|---|---|
| `reflection` | yes | yes | valid |
| `transmission` | no | yes | invalid |
| `label_exchange_reflection` | yes | no | invalid |
| `transmission_with_label_exchange` | no | no | invalid |

Only the reflection map satisfied direction reversal and identification-mode preservation simultaneously.

![Control map trajectories](elastic_collision_control_maps_result_v1/control_maps_trajectories_v1.png)

![Control map verdicts](elastic_collision_control_maps_result_v1/control_maps_verdict_v1.png)

---

## 7. Asymmetric Conditions

The purpose was to test whether amplitude asymmetry, harmonic-order asymmetry, or temporal-rate asymmetry breaks the construction.

### 7.1 Result

| Item | Value |
|---|---:|
| Total cases | `10` |
| Valid cases | `7` |
| Invalid cases | `3` |

Amplitude difference and harmonic-order difference alone did not break the map. Failures were caused by the loss of simultaneous satisfaction of the spatial and temporal cell conditions.

![Asymmetry cell gaps](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_cell_gaps_v1.png)

![Asymmetry verdicts](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_verdict_v1.png)

---

## 8. Observation Perturbation

The purpose was to test finite readout disturbances caused by observation.

### 8.1 Result

| Item | Value |
|---|---:|
| Total cases | `8` |
| Observation-model valid cases | `5` |
| Collision-map valid cases | `7` |
| Invalid cases | `3` |

When the perturbation exceeded the localization width of observer `C`, the observation model failed. If the AB finite-cell condition remained satisfied, however, the collision map itself could still be valid.

![Observation perturbation thresholds](elastic_collision_observation_perturbation_result_v1/observation_perturbation_thresholds_v1.png)

![Observation perturbation verdicts](elastic_collision_observation_perturbation_result_v1/observation_perturbation_verdict_v1.png)

---

## 9. Multiple Collisions

The purpose was to test whether the construction survives repeated application.

### 9.1 Result

| Item | Value |
|---|---:|
| Target AB collisions | `8` |
| Actual AB collisions | `8` |
| Wall reflections | `14` |
| Direction reversal at each AB collision | `true` |
| Identification modes preserved | `true` |
| Representative amplitudes preserved | `true` |
| Closure preserved | `true` |
| Overall validity | `true` |

The repeated-collision test shows that the result is not merely a one-shot demonstration. The map remains repeatable under the tested closed-interval dynamics.

![Multiple collision trajectory](elastic_collision_multi_collision_result_v1/multi_collision_trajectory_v1.png)

![Multiple collision closure residual](elastic_collision_multi_collision_result_v1/multi_collision_closure_v1.png)

---

## 10. Readout Resolution of the `η` Identification Oscillation

The purpose was to test finite sampling of the internal identification phase `η`.

### 10.1 Result

| Item | Value |
|---|---:|
| Total cases | `88` |
| Valid cases | `59` |
| Invalid cases | `29` |
| Aliasing cases | `29` |
| Non-aliasing failures | `0` |
| Minimum `η` sample count valid for all tested mode pairs | `64` |

All invalid cases were aliasing cases. No non-aliasing failures were observed. Therefore, the failures in this test are readout-resolution failures, not failures of the collision map itself.

![Eta resolution validity](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_validity_v1.png)

![Eta resolution purity](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_purity_v1.png)

---

## 11. Validity Conditions

The experiments support the following validity conditions for this simulation system.

| Condition | Meaning |
|---|---|
| Observer heaviness | `A_C` must be large enough for quasi-static readout |
| Finite-cell arrival | `A,B` must satisfy spatial and temporal cell conditions simultaneously |
| Time step | `Δs` must not skip the finite cell |
| Identification separability | `m_A,m_B` must be resolvable at the chosen `η` sampling |
| Identification purity | leakage must not invert the dominant identification mode |
| Direction reversal | `q_A,q_B` must reverse inside the interaction cell |
| Representative amplitude | `A_A,A_B` must be preserved |
| Closure | compensated square closure must remain intact |

---

## 12. Failure Conditions

The following failure conditions were separated.

| Failure condition | Observed form |
|---|---|
| Observer insufficiency | quasi-static condition fails when `A_C` is too small |
| Cell skipping | update step skips over the finite-resolution cell |
| Temporal mismatch | spatial proximity occurs without temporal-cell coincidence |
| Observation perturbation | readout perturbation exceeds the width of `C` or the AB cell |
| Identification leakage | dominant identification mode is no longer cleanly readable |
| `η` aliasing | finite sampling cannot distinguish different identification modes |
| Transmission | direction readout does not reverse |
| Label exchange | identification oscillations are swapped |

---

## 13. Unevaluated Items and Next Steps

The following are not evaluated as derivations from first principles in these experiments:

- unique derivation of the direction-reversal rule from the first axioms;
- connection to standard fermion scattering;
- connection to physical cross sections;
- connection to standard quantum measurement theory;
- full three-body dynamical solution of `A,B,C`;
- direct physical-space interpretation of `χ` and `τ`.

These are outside the scope of the present constructive experiment.

---

## 14. Reference Files

| Experiment | Report | Result JSON | CSV |
|---|---|---|---|
| Basic collision | [report](elastic_collision_simulation_result_v1/elastic_collision_report_v1.md) | [json](elastic_collision_simulation_result_v1/elastic_collision_result_v1.json) | [timeline](elastic_collision_simulation_result_v1/elastic_collision_timeline_v1.csv) |
| Identification robustness | [report](elastic_collision_label_robustness_result_v1/label_robustness_report_v1.md) | [json](elastic_collision_label_robustness_result_v1/label_robustness_result_v1.json) | [csv](elastic_collision_label_robustness_result_v1/label_robustness_cases_v1.csv) |
| Observer sweep | [report](elastic_collision_observer_sweep_result_v1/observer_sweep_report_v1.md) | [json](elastic_collision_observer_sweep_result_v1/observer_sweep_result_v1.json) | [csv](elastic_collision_observer_sweep_result_v1/observer_sweep_cases_v1.csv) |
| Cell resolution | [report](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_report_v1.md) | [json](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_result_v1.json) | [csv](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_cases_v1.csv) |
| Control maps | [report](elastic_collision_control_maps_result_v1/control_maps_report_v1.md) | [json](elastic_collision_control_maps_result_v1/control_maps_result_v1.json) | [csv](elastic_collision_control_maps_result_v1/control_maps_cases_v1.csv) |
| Asymmetry sweep | [report](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_report_v1.md) | [json](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_result_v1.json) | [csv](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_cases_v1.csv) |
| Observation perturbation | [report](elastic_collision_observation_perturbation_result_v1/observation_perturbation_report_v1.md) | [json](elastic_collision_observation_perturbation_result_v1/observation_perturbation_result_v1.json) | [csv](elastic_collision_observation_perturbation_result_v1/observation_perturbation_cases_v1.csv) |
| Multiple collisions | [report](elastic_collision_multi_collision_result_v1/multi_collision_report_v1.md) | [json](elastic_collision_multi_collision_result_v1/multi_collision_result_v1.json) | [timeline](elastic_collision_multi_collision_result_v1/multi_collision_timeline_v1.csv) |
| `η` resolution | [report](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_report_v1.md) | [json](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_result_v1.json) | [csv](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_cases_v1.csv) |
