# Preliminary Summary of Harmonic Readout and c=1 Area Sweep in an AB Two-Body Closed Phase System v4

**Date:** 2026-07-15
**Author:** Noriaki Kihara
**Position:** Summary of the AB two-body preliminary experiment group in the Wave Information Readout series
**Version DOI:** 10.5281/zenodo.21374317
**Concept DOI:** 10.5281/zenodo.21318696

V4 implements the fermion-like recoil map as a direct A/B two-channel scattering matrix acting on the incident channels.

---

## Abstract

This summary consolidates six preliminary experiments performed on an AB two-body closed phase system. V4 keeps the same scope and claims and incorporates a fermion-like recoil map into the same AB acceleration-readout frame as a direct A/B two-channel scattering action:

1. one-angle circumferential phase harmonic readout,
2. parameter sweep of the one-angle circumferential phase harmonic readout,
3. internal `c=1` calibration and `chi-tau` area sweep,
4. parameter sweep of the internal `c=1` calibration,
5. inverse-area compensation diagnosis in the `chi-tau` plane, and
6. extended sweep for native inverse-area readout, and
7. direct two-channel fermion-like recoil harmonic readout protocol.

The preceding ABC multigauge interference readout paper showed that, in a minimal ABC closed phase system that does not place a background space first, the measuring device itself can be defined as a complex phase wave inside the same system, and mass-like, momentum-like, and energy-like conserved readouts can be constructed from interference.

The present experiment group takes one step further. It verifies that an acceleration-like readout can be represented by harmonic readout in the AB two-body relation itself.

However, the AB two-body system alone could not determine whether the acceleration-like readout changes with the position phase difference, or distance, according to an inverse law or an inverse-square law.

The central achievement of this experiment group is not merely the display of harmonic motion. No external standard force, spring equation, gravitational equation, or Coulomb equation is introduced. No individual forces `f_A` and `f_B` are introduced. Using only the label-free two-body relation `f_AB`, an acceleration-like readout was constructed.

At the same time, the experiment clarified a boundary. In the AB two-body system, the measuring device is on the measured object itself. Therefore, although relative-distance change can be read, there is no independent gauge that can decide whether the change follows

```text
L_AB
1 / L_AB
1 / L_AB^2
```

The gauge used to measure the interval is itself affected by the same closed phase system. Even after adding a temporal phase difference in addition to the spatial position phase difference, and even after imposing an internal `c=1`-like calibration to reproduce a time-development-like effect, the conclusion did not change.

Thus, the determination of the distance exponent remains a task for a separate ABC three-body experiment with an independent metrical gauge.

**Keywords:** AB two-body closed phase system, label-free two-arc relative phase, harmonic readout, acceleration-like readout, c=1 internal calibration, chi-tau area, inverse-area diagnosis, native readout

---

## 1. Experiment List

| No. | Experiment | Purpose | Main result |
|---:|---|---|---|
| 1 | One-angle harmonic readout | Test whether `D_AB`, `V_AB`, and `f_AB` can be read label-free | Readable |
| 2 | One-angle parameter sweep | Robustness against initial deviation, period, and readout leakage | Harmonic oscillation detected in all cases |
| 3 | `c=1` area sweep | Separate `s` from `tau_read` and test whether a `chi-tau` plane stands | Area detected under independent `tau` |
| 4 | `c=1` parameter sweep | Test whether `c=1` is sufficient for area formation | Not sufficient |
| 5 | Inverse-area compensation diagnosis | Test whether `1/A_chi_tau` remains natively on the closure-compensation side | Not detected |
| 6 | Native inverse-area extended sweep | Search broad conditions for `alpha≈2` | Not detected natively |
| 7 | V4 fermion-like recoil protocol | Re-test acceleration-like readout by placing an A/B two-channel scattering matrix near the AB center | `q_out_factor` is not used as an operator; harmonic readout and the `chi-tau` surface are re-read from the outgoing channel difference |

---

## 2. One-Angle Circumferential Phase Harmonic Readout

### 2.1 Meaning of the Experiment

The AB pair was read as a closed phase system without placing an observer `C`.

The experiment did not introduce

```text
f_A
f_B
```

as independent forces. It read only

```text
f_AB
```

as relational compensation.

No standard gravitational equation, standard Coulomb equation, or standard spring equation was used.

### 2.2 Main Result

The main numerical results were:

| Quantity | Value |
|---|---:|
| `case_count` | `32` |
| `observer_C_used` | `False` |
| `f_A_or_f_B_used` | `False` |
| `max_Q_closed_abs` | `0.0` |
| `oscillation_detected_all_cases` | `True` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `readout_decay_monotonic_all_cases` | `True` |
| `readout_off_decay_max_abs` | `6.49e-17` |
| `readout_strong_decay_min_abs` | `4.0004e-4` |

The result was:

```text
Protocols F and B differ as internal representations,
but they are completely degenerate in D_AB and V_AB.
```

This means that reflection-type and pass-through-type protocols cannot be distinguished by label-free readout.

It was also confirmed that

```text
no decay occurs when readout is off,
and the envelope decays more strongly as the readout becomes stronger.
```

This is a falsification test showing that a readout wave can affect the envelope of the system when information is read out from the closed phase system.

### 2.3 Figures

| AB two-body geometry | Main observation summary |
|---|---|
| <img src="AB二体問題の図化_fAB_v1.png" width="420"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_observation_summary_v1.png" width="520"> |

| Harmonic state | Readout leak response |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_harmonic_state_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_readout_leak_response_v1.png" width="520"> |

| Protocol degeneracy | Phase-difference scaling |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_protocol_degeneracy_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_observation_figures_v1/ab_two_body_one_angle_phase_difference_scaling_v1.png" width="520"> |

---

## 3. One-Angle Parameter Sweep

### 3.1 Main Result

| Quantity | Value |
|---|---:|
| `sweep_configuration_count` | `210` |
| `case_summary_count` | `420` |
| `period_count` | `5` |
| `deviation_count` | `7` |
| `leak_count` | `6` |
| `max_Q_closed_abs` | `0.0` |
| `label_free_protocol_degenerate_all_cases` | `True` |
| `oscillation_detected_all_cases` | `True` |
| `readout_off_decay_max_abs` | `8.26e-17` |
| `decay_abs_monotonic_by_leak_all_grids` | `True` |
| `max_normalized_f_AB_projection_error` | `2.58e-4` |

This sweep confirmed that the one-angle version is stable under weak readout, while a strong readout wave distorts the `f_AB` projection.

| Leakage and projection error | Leakage and decay |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_leak_f_error_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_leak_decay_v1.png" width="520"> |

| Period sweep | Selected time series |
|---|---|
| <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_period_v1.png" width="520"> | <img src="ab_two_body_one_angle_harmonic_readout_parameter_sweep_preliminary_result_v1/ab_two_body_one_angle_parameter_sweep_selected_series_v1.png" width="520"> |

---

## 4. c=1 Internal Calibration and chi-tau Area Sweep

### 4.1 Meaning of the Experiment

In the one-angle harmonic readout, distance attenuation according to phase difference was not read.

Therefore, the experiment step `s` was not treated as time itself. Instead, an independent temporal phase candidate

```text
tau_read
```

was introduced.

The purpose was to test whether

```text
chi_read and tau_read form an independent plane.
```

### 4.2 Main Result

| Quantity | Value |
|---|---:|
| `case_summary_count` | `288` |
| `power_candidate_count` | `48` |
| `max_Q_closed_abs` | `0.0` |
| `disabled_max_area` | `0.0` |
| `locked_max_area` | `7.11e-15` |
| `independent_min_area` | `0.0024367633602631385` |
| `c1_readout_off_max_epsilon_c_abs` | `2.33e-15` |
| `c1_area_sweep_detected_all_cases` | `True` |
| `tau_is_step_used_any` | `False` |
| `external_c_used_any` | `False` |

The following was established.

```text
tau disabled does not create area.
tau locked does not create area.
tau independent creates chi-tau area.
```

Thus, an independent temporal phase can make the `chi-tau` plane stand.

However, this alone does not generate distance attenuation.

### 4.3 Figures

| `chi-tau` surface | Area sweep |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_area_sweep_v1.png" width="520"> |

| `c=1` calibration error | Power candidates |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_error_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_power_candidate_v1.png" width="520"> |

| Readout leakage |
|---|
| <img src="ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1/ab_two_body_c1_internal_calibration_readout_leak_v1.png" width="620"> |

---

## 5. c=1 Internal Calibration Parameter Sweep

### 5.1 Main Result

| Quantity | Value |
|---|---:|
| `sweep_case_count` | `756` |
| `readout_off_case_count` | `252` |
| `rank2_readout_off_count` | `246` |
| `c1_surface_like_readout_off_count` | `13` |
| `c1_locked_like_readout_off_count` | `8` |
| `min_c_error_readout_off` | `1.11e-15` |
| `max_area_readout_off` | `0.19685536479742288` |

This clarified that the following three conditions must be required simultaneously:

```text
c=1
rank_chi_tau = 2
A_chi_tau != 0
```

Being close to `c=1` alone does not mean that a `chi-tau` surface has been formed.

| c-error heatmap | Area heatmap |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_c_error_heatmap_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_area_heatmap_v1.png" width="520"> |

| Phase-difference response | Readout leakage |
|---|---|
| <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_phase_response_v1.png" width="520"> | <img src="ab_two_body_c1_internal_calibration_parameter_sweep_preliminary_result_v1/ab_two_body_c1_parameter_sweep_readout_leak_v1.png" width="520"> |

---

## 6. Inverse-Area Compensation Diagnosis

### 6.1 Purpose

After the `chi-tau` area was formed, the experiment tested whether

```text
1 / A_chi_tau
```

would naturally remain on the closure-compensation side.

The important point is not to construct `1/A_chi_tau` as post-processing.

That is a constructed control and must be distinguished from native readout.

### 6.2 Main Result

| Quantity | Value |
|---|---:|
| `diagnostic_case_count` | `288` |
| `area_valid_case_count` | `144` |
| `fit_count` | `576` |
| `native_fit_count` | `132` |
| `derived_fit_count` | `84` |
| `native_positive2_count` | `0` |
| `constructed_reciprocal_positive2_count` | `2` |

The result is clear.

```text
If 1 / A_chi_tau is constructed, alpha≈2 appears.
However, alpha≈2 does not appear in native readout.
```

| Constructed inverse-area control | Native candidates |
|---|---|
| <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_constructed_control_v1.png" width="520"> | <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_native_candidates_v1.png" width="520"> |

| Alpha comparison |
|---|
| <img src="ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1/ab_two_body_chi_tau_inverse_area_alpha_comparison_v1.png" width="620"> |

---

## 7. Native Inverse-Area Extended Sweep

### 7.1 Main Result

| Quantity | Value |
|---|---:|
| `sweep_case_count` | `1323` |
| `area_valid_case_count` | `1260` |
| `c1_surface_like_case_count` | `126` |
| `fit_count` | `4158` |
| `native_fit_count` | `1056` |
| `native_positive2_count` | `0` |
| `c1_native_positive2_count` | `0` |
| `constructed_reciprocal_positive2_count` | `198` |

Across a broad parameter range,

```text
native inverse-area scaling was not detected.
```

On the other hand, the constructed control does produce `alpha≈2`.

Thus the boundary found here is:

```text
The inverse-square form can be constructed.
However, it has not yet been read natively.
```

| Reference curves | Alpha scan |
|---|---|
| <img src="ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_reference_curves_v1.png" width="520"> | <img src="ab_two_body_chi_tau_native_inverse_area_extended_sweep_preliminary_result_v1/ab_two_body_native_inverse_area_extended_alpha_scan_v1.png" width="520"> |

---

## 8. V4 Addition: Harmonic Readout with a Fermion-Like Recoil Map

### 8.1 Purpose

In the V1 one-angle harmonic readout, the relative phase development of the AB pair was read continuously, and Protocol F/B degeneracy was confirmed under label-free readout.

V4 keeps this context unchanged and incorporates a fermion-like recoil map near the central interaction region.

However, V4 does not use the compressed implementation that multiplies `chi_read` by `q_out_factor`.

Instead, it explicitly defines the A/B incident channels and applies the scattering matrix directly.

```text
psi_A_in =  (chi_pass + i eta_pass) / 2
psi_B_in = -(chi_pass + i eta_pass) / 2
Psi_out  = S_Delta Psi_in
chi_read = Re(psi_A_out - psi_B_out)
eta_read = Im(psi_A_out - psi_B_out)
```

The test is whether the acceleration-like harmonic readout from the AB two-body relation `f_AB` can still be evaluated in the same form when the recoil map is inserted as A/B two-channel scattering.

This section does not test:

```text
AB asymmetrization
low-localization limit
instantaneous contraction
harmonic transfer
localization transfer
```

Those questions belong to a separate wave-packet focusing experiment line.

The present V4 evaluates the fermion-like recoil protocol within the current AB two-body acceleration-readout assumptions, in addition to the transmission-type central readout.

### 8.2 Recoil Map

Let the right-moving and left-moving incident channels be

```math
\Psi^{\mathrm{in}}
=
\begin{pmatrix}
\psi_+^{\mathrm{in}}\\
\psi_-^{\mathrm{in}}
\end{pmatrix}.
```

Let `Delta_F` be the fermion-like exchange phase. Define the transmission and reflection amplitudes by

```math
t_\Delta
=
e^{i\Delta_F/2}
\cos\frac{\Delta_F}{2}
```

```math
r_\Delta
=
-i e^{i\Delta_F/2}
\sin\frac{\Delta_F}{2}.
```

The two-channel scattering matrix is

```math
S_\Delta
=
\begin{pmatrix}
t_\Delta & r_\Delta\\
r_\Delta & t_\Delta
\end{pmatrix},
```

and the outgoing channels are read as

```math
\Psi^{\mathrm{out}}
=
S_\Delta
\Psi^{\mathrm{in}}.
```

That is,

```math
\psi_+^{\mathrm{out}}
=
t_\Delta \psi_+^{\mathrm{in}}
+
r_\Delta \psi_-^{\mathrm{in}}
```

```math
\psi_-^{\mathrm{out}}
=
r_\Delta \psi_+^{\mathrm{in}}
+
t_\Delta \psi_-^{\mathrm{in}}.
```

For the complete-recoil condition,

```math
\Delta_F=\pi,
```

we obtain the equivalent condition

```text
transmission amplitude t_Delta = 0
reflection amplitude r_Delta = 1
```

In this limit the outgoing channels are read as an exchange of the incident channels.

### 8.3 Pre-Experiment Criteria

The V4 additional protocol checks at least the following items.

| Item | Criterion |
|---|---|
| Two-channel scattering | `S_Delta` acts on the A/B incident channels without using `q_out_factor` as an operator |
| Unitarity | Channel norm error remains within numerical error |
| Closure residual after recoil | `Q_closed` does not break down |
| Harmonic readout | The harmonic structure of `D_AB`, `V_AB`, and `f_AB` is maintained |
| Readout-leak response | The difference between readout-off and strong-readout conditions follows the same discipline as V1 |
| Protocol comparison | Transmission, display-recoil, and fermion-like recoil protocols are all recorded |
| `chi-tau` surface | Independent `tau_read` area readout is not broken |

### 8.4 Main Result

| Quantity | Value |
|---|---:|
| `harmonic_case_count` | `48` |
| `c1_case_count` | `144` |
| `scattering_protocol_count` | `3` |
| `fermionic_delta_f` | `3.141592653589793` |
| `fermionic_reflection_rate` | `1.0` |
| `fermionic_transmission_rate` | `3.749399456654644e-33` |
| `fermionic_q_out_factor` | `-1.0` |
| `q_out_factor_applied_any` | `False` |
| `fermionic_full_two_channel_scattering_used_all_cases` | `True` |
| `fermionic_max_scattering_unitarity_error` | `4.163336342344337e-17` |
| `max_Q_closed_abs` | `0.0` |
| `fermionic_regular_cell_harmonic_consistent_nonstrong_modes` | `True` |
| `fermionic_strong_readout_perturbs_harmonic_projection` | `True` |
| `fermionic_max_f_AB_projection_error_regular_nonstrong` | `2.2756072064286935e-6` |
| `fermionic_max_f_AB_projection_error_regular_strong` | `9.072436576166898e-6` |
| `fermionic_reflection_event_cell_count_total` | `720` |
| `label_free_pass_vs_fermionic_match_all_cases` | `True` |
| `label_free_display_vs_fermionic_match_all_cases` | `True` |
| `readout_off_decay_max_abs` | `3.818134371020083e-18` |
| `readout_strong_decay_min_abs` | `4.0004000533409265e-4` |
| `fermionic_tau_disabled_max_area` | `0.0` |
| `fermionic_tau_locked_max_area` | `8.881784197001252e-16` |
| `fermionic_tau_independent_min_area` | `0.00375268307733001` |
| `fermionic_c1_readout_off_max_epsilon_c_abs` | `2.3314683517128287e-15` |
| `fermionic_c1_area_sweep_detected_all_cases` | `True` |
| `tau_is_step_used_any` | `False` |
| `external_c_used_any` | `False` |
| `f_A_or_f_B_used_any` | `False` |

The result is read as follows.

```text
In the complete-recoil limit, q_out_factor=-1 is obtained as a diagnostic quantity.
However, V4 does not multiply chi_read by q_out_factor as an operator.
```

V4 acts with the two-channel scattering matrix on the A/B incident channels, and re-reads `chi_read` and `eta_read` from the outgoing channel difference.

For all fermion-like recoil cases, `full_two_channel_scattering_used=True`, and the maximum unitarity error was `4.16e-17`.

Under that condition, the label-free readouts `D_AB` and `V_AB` matched the V1 transmission-type readout.

This is consistent with the Protocol F/B degeneracy confirmed in V1.

The display-only foldback protocol `display_reflection` was also recorded as a control.

The display-recoil and fermion-like recoil protocols also matched under the label-free `D_AB` and `V_AB` readouts.

For harmonic readout, projection consistency was maintained under non-`readout_strong` conditions, as in V1.

Strong readout produced projection distortion, but this is the same readout-wave perturbation observed in V1 and is not a breakdown specific to the recoil map.

The `chi-tau` surface also maintained the following structure:

```text
tau disabled: no area
tau locked: no area
tau independent c1: area present
```

Therefore, even when the fermion-like recoil map is incorporated into AB two-body acceleration readout as A/B two-channel scattering, the acceleration-like harmonic readout and the `chi-tau` area readout are not broken within the present assumptions.

### 8.5 Figures

| Protocol comparison | Channel state |
|---|---|
| <img src="ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/ab_two_body_fermionic_reflection_protocol_comparison_v4.png" width="520"> | <img src="ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/ab_two_body_fermionic_reflection_channel_state_v4.png" width="520"> |

| Readout decay | `chi-tau` path |
|---|---|
| <img src="ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/ab_two_body_fermionic_reflection_readout_decay_v4.png" width="520"> | <img src="ab_two_body_fermionic_reflection_harmonic_readout_preliminary_result_v4/ab_two_body_fermionic_reflection_c1_chi_tau_v4.png" width="520"> |

---

## 9. Overall Judgment

### 9.1 Established

| Item | Judgment |
|---|---|
| Read `f_AB` as a single AB relation | Established |
| Read harmonic oscillation without placing `f_A` or `f_B` | Established |
| Confirm label-free degeneracy of Protocol F/B | Established |
| Confirm disappearance of decay when the readout wave is stopped | Established |
| Read `chi-tau` area through independent `tau_read` | Established |
| Confirm `c=1` as a candidate necessary condition | Established |
| Post-processed `1/A_chi_tau` shows `alpha≈2` | Established |
| Apply the fermion-like recoil map as an A/B two-channel scattering matrix to the incident channels | Established |
| Do not use `q_out_factor` as an operator | Established |
| Keep the two-channel scattering unitarity error within numerical error | Established |
| Confirm that label-free `D_AB` and `V_AB` match V1 even after introducing fermion-like recoil | Established |
| Maintain the independent-`tau_read` `chi-tau` surface even after introducing fermion-like recoil | Established |

### 9.2 Not Established

| Item | Judgment |
|---|---|
| Independent metrical determination of distance exponent in AB two-body system alone | Not established |
| Formation of temporal phase surface from `c=1` alone | Not established |
| Native inverse-law readout | Not detected |
| Native inverse-square readout | Not detected |
| Automatic conversion from `chi-tau` area to `1/A` compensation | Not detected |
| Correspondence to standard gravity | Not established |

---

## 10. Interpretation

The main achievements of the AB two-body system are the following two points:

```text
1. Using only the AB relation f_AB, an acceleration-like harmonic readout can be constructed.
2. The AB two-body system alone lacks an independent gauge for reading the distance exponent.
```

The first point is a positive result.

In the preceding ABC multigauge interference readout experiment, mass-like, momentum-like, and energy-like readouts were reconstructed from interference correlations.

The present experiment shows, as a next step, that an acceleration-like readout can be constructed from the harmonic AB two-body relation.

This is not the result of inserting a standard force from outside.

It is a harmonic displacement appearing as a closure readout of the two-body relation `f_AB`.

The second point is a boundary condition.

Inside the AB pair, A and B can be read only through one another.

Therefore, relative-distance change can be read, but there is no independent gauge that fixes whether the change rate should be measured as

```text
proportional
inverse proportional
inverse square
```

This is not a failure of namelessness. Rather, it is a consequence of namelessness.

What cannot be read must not be declared read by post-processing.

In the `chi-tau` area sweep, an independent temporal phase allowed an area readout to stand.

However, this did not automatically provide a distance-exponent readout.

In addition, if `1/A_chi_tau` is constructed in post-processing, an inverse-square form naturally appears. But this is not native readout.

The judgment of this summary is therefore:

```text
Acceleration-like readout was confirmed.
Distance-exponent readout cannot be determined in the AB two-body system alone.
```

The V4 additional protocol did not change this judgment.

In the complete-recoil limit, the diagnostic quantity `q_out_factor=-1` was obtained, but V4 does not use it as an operator.

After the A/B two-channel scattering matrix acted on the incident channels, the label-free readouts `D_AB` and `V_AB` matched the V1 transmission-type readout.

This means that the V1 result, in which recoil-type and transmission-type protocols are degenerate under AB two-body nameless readout, was confirmed again using the two-channel scattering-matrix recoil map.

---

## 11. Connection to the Next Experiment

The next stage following this summary is the ABC three-body system.

However, the third wave `C` is not an external absolute observer.

`C` is itself part of the closed phase system and generates

```text
f_AC
f_BC
f_ABC
```

Therefore, the ABC experiment must test:

```text
whether C can be used as an independent metrical gauge,
whether C contaminates the main AB readout,
whether f_ABC can be used as representative time,
and whether f_AB, f_AC, and f_BC can be separated as circumferential-direction candidates.
```

---

# References

## Self-References

1. Noriaki Kihara, "Basic Axiom System of the Nameless Equal-Amplitude Composite Wave Model v4", Version DOI: `10.5281/zenodo.21316620`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. Noriaki Kihara, "Construction Experiment of Multigauge Interference Readout Conserved Quantities in an ABC Closed Phase System", Version DOI: `10.5281/zenodo.21332875`, Concept DOI: `10.5281/zenodo.21308049`, 2026.
3. Noriaki Kihara, [Definitional Supplement on Label-Free Two-Arc Relative Phase and Harmonic Readout in an AB Two-Body Closed Phase System](AB二体閉鎖位相系におけるラベルなし二弧相対位相と調和読出しに関する定義補足.md), 2026.
4. Noriaki Kihara, [Experiment Specification for One-Angle Circumferential Phase Harmonic Readout in an AB Two-Body Closed Phase System v1](AB二体閉鎖位相系における一角度円周位相調和読出し実験仕様書%20v1.md), 2026.
5. Noriaki Kihara, [Experiment Specification for c=1 Internal Calibration and Spatial-Phase/Temporal-Phase Area Sweep Readout in an AB Two-Body Closed Phase System v1](AB二体閉鎖位相系におけるc=1内部較正と空間位相・時間位相面積スイープ読出し実験仕様書%20v1.md), 2026.
6. Noriaki Kihara, [Definitional Supplement on Internal Closure of Self-Terms and Separation of N-Body External Readout in Closed Complex Phase Waves](閉鎖複素位相波における自己項の内部閉鎖とN体外部読出し分離に関する定義補足.md), 2026.

## External References

The external references are used only as minimal standard background on wave readout, phase, and observational probability, not as derivational premises of this paper.

7. Max Born, "Zur Quantenmechanik der Stossvorgaenge", *Zeitschrift fuer Physik* 37, 863-867, 1926. DOI: `10.1007/BF01397477`.
8. Y. Aharonov and D. Bohm, "Significance of electromagnetic potentials in the quantum theory", *Physical Review* 115, 485-491, 1959. DOI: `10.1103/PhysRev.115.485`.
