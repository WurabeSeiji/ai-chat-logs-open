# Construction Experiment of Multigauge Interference Readout Conserved Quantities in an ABC Closed Phase System

**Subtitle:** Numerical construction of mass-like, momentum-like, and energy-like quantities using `p_read`, `E_read`, `R_read`, and `R*p`, `R*p^2` conserving maps
**Date:** 2026-07-11
**Author:** Noriaki Kihara
**Position:** Additional paper in the Wave Information Readout series
**Version DOI:** 10.5281/zenodo.21308050
**Concept DOI:** 10.5281/zenodo.21308049

---

## Abstract

This paper numerically tests whether conserved readouts corresponding to mass-like, momentum-like, and energy-like quantities can be constructed from multigauge interference in an ABC closed phase system that does not assume background coordinates in advance.

In this series, the variables `chi` and `tau`, which correspond to space and time, have been treated not as external spacetime coordinates but as readout quantities from a closed phase system. This paper extends the same policy to mass-like, momentum-like, and energy-like quantities.

For an ABC system consisting of local waves `A`, `B`, and an internal observer wave `C`, a single gauge value is not accepted as a measurement. Instead, interference correlations are acquired from multiple readout gauges. The correlation gradient in the spatial phase direction is read as

```text
p_read
```

the correlation gradient in the temporal phase direction is read as

```text
E_read
```

and the amplitude-square residual that remains stable across multiple gauges is read as

```text
R_read
```

The quantities `p,E,R` in this paper are not standard physical momentum, energy, and mass themselves. They are multigauge interference readout quantities internal to the closed phase system. The purpose of this paper is not to claim identity with standard quantities, but to test whether readout structures that behave like conserved quantities can be constructed inside the ABC closed phase system.

In the single ABC collision, `p,E,R` were reconstructed from multiple gauges with maximum errors `2.5202062658991053e-14`, `2.2315482794965646e-14`, and `4.440892098500626e-16`, respectively. Across eight repeated collisions, `p` reversal, `E,R` preservation, label-mode preservation, and compensated closure were maintained.

Under asymmetric `R` conditions, the simple `q -> -q` reversal was detected to break `R*p` conservation. A generalized collision map conserving `R_A p_A + R_B p_B` and `R_A p_A^2 + R_B p_B^2` was therefore constructed. Across eight asymmetric amplitude cases, conservation was confirmed with maximum `R*p` error `2.3803181647963356e-13` and maximum `R*p^2` error `1.4086509736443986e-12`.

Additional verification covered nine non-unit and asymmetric initial phase-gradient cases, four repeated-collision cases, readout-noise robustness, and an extreme `R`-ratio sweep from `R_B/R_A=0.015625` to `64.0`. The integration summary reported all nine experiments as `valid`, and no judgment used single-gauge-only measurement.

Within the numerical constructive scope of this paper, the ABC closed phase system therefore supports mass-like, momentum-like, and energy-like conserved readouts constructed consistently through multigauge interference.

**Keywords:** multigauge interference readout, ABC closed phase system, all-positive zero closure, mass-like readout, momentum-like readout, energy-like readout, generalized elastic collision, `R*p` conservation, `R*p^2` conservation

---

## 1. Introduction

### 1.1 Background

This series adopts namelessness, all-positive zero closure, and non-trivial existence as its basic axioms. The central closure condition is

```math
\sum_n x_n^2=0.
```

This condition is not the conjugate norm

```math
\sum_n |x_n|^2.
```

Each component is squared directly, and all terms are summed with positive signs.

The preceding definitional supplement, "Definitional Supplement on Readout Multiplicity of the All-Positive Zero Closure," kept Axiom 1 unchanged and organized how the same closure condition can have multiple readout representations. When a radial representation is read as

```math
a^2+b^2=\rho^2,
```

the first-principle layer does not introduce an externally negative metric expression

```math
a^2+b^2-\rho^2=0.
```

Instead, it reads the relation as an all-positive zero-closure representation:

```math
a^2+b^2+(i\rho)^2=0.
```

Thus quantities that appear as radius, time, mass, energy, or momentum are not placed first as axes with intrinsic names. They are treated as readout representations.

### 1.2 Question

The question of this paper is:

> In an ABC closed phase system without assuming background coordinates in advance, can conserved readouts that look like mass-like, momentum-like, and energy-like quantities be constructed from multigauge interference?

The important point is that a single gauge value is not accepted as a measurement.

Just as readouts corresponding to `chi` and `tau` were constructed from interference correlations, readouts corresponding to `p,E,R` must also be reconstructed from multiple reference waves, multiple readout windows, and multiple gauges.

### 1.3 Minimal Claim

This paper limits its claim to the following scope.

1. Reconstruct `p_read`, `E_read`, and `R_read` from multigauge interference inside an ABC closed phase system.
2. Confirm that symmetric collision reverses `p_read` while preserving `E_read` and `R_read`.
3. Detect that simple reversal breaks `R*p` conservation under asymmetric `R` conditions.
4. Construct a generalized collision map conserving `R*p` and `R*p^2`.
5. Verify whether readout conservation is maintained under repeated collisions, readout-gauge changes, readout noise, and extreme `R` ratios.

---

## 2. What This Paper Does Not Claim

This paper does not claim the following.

| Not claimed | Reason |
|---|---|
| Derivation of standard momentum, standard energy, or standard mass | `p,E,R` are readout quantities internal to the closed phase system |
| Complete re-derivation of standard mechanics | A correspondence map must be constructed separately |
| Quantitative prediction of real particle collisions | This is a numerical constructive experiment on an internal phase model |
| General impossibility of single-gauge measurement | This paper requires multigauge interference readout as its measurement condition |
| Internal interference generation of the `R`-weighted generalized map | This paper constructs the map from conservation conditions and verifies it by multigauge readout |
| Treating `R` as a primitive mass axis | `R` is a readout name for the stable residual remaining across multiple gauges |
| Assumption of background spacetime coordinates | `chi,tau` are readout variables internal to the phase system |
| Measurement validity under arbitrary noise | The noise test is limited to a controlled comparison of zero-mean gauge fluctuations and common bias |
| Generalized collision validity under arbitrary conditions | The verification range is limited to the numerical conditions of this paper |

The claim of this paper is that a readout structure that behaves like conserved quantities can be constructed from multigauge interference inside an ABC closed phase system.

---

## 3. Basic Structure

### 3.1 All-Positive Zero Closure

Axiom 1 of the basic axiom system v2 is

```math
Q(x)=\sum_n x_n^2=0.
```

This is not a condition that first introduces an external negative-sign metric.

For the minimal compensating pair

```math
A,\qquad iA,
```

one obtains

```math
A^2+(iA)^2=0.
```

The negative sign is not an external coefficient. It arises from the square of the internal phase:

```math
i^2=-1.
```

### 3.2 Readout Multiplicity

The same closure condition

```math
\sum_n x_n^2=0
```

can have multiple representations under local readout.

For example, if a local representation reads

```math
a^2+b^2+(i\rho)^2=0,
```

then `a,b,rho` are not intrinsically different kinds of components at the first-principle layer. They are labels assigned by a readout window, a reference wave, or a projection.

The same applies to `p,E,R` in this paper.

`p` is not momentum itself. It is a correlation-gradient readout in the spatial phase direction.

`E` is not energy itself. It is a correlation-gradient readout in the temporal phase direction.

`R` is not mass itself. It is an amplitude-square readout that remains stable across multiple gauges.

---

## 4. ABC Closed Phase System

### 4.1 Local Waves A and B

The local waves `A,B` have spatial phase `chi`, temporal phase `tau`, internal label phase `eta`, representative amplitude, and internal label modes.

The phase centers of `A,B` are denoted by `chi_A,chi_B`, and their temporal phase centers by `tau_A,tau_B`.

The internal label modes are

```text
m_A=1,
m_B=2.
```

### 4.2 Observer Wave C

The observer wave `C` is not an external observer.

It is a reference wave inside the same closed phase system and is used to read interference correlations with local waves `A,B`.

Each gauge changes the readout center, readout width, phase shift, reference-wave gain, and related settings.

### 4.3 Single Gauge Values Are Not Measurements

This paper does not treat a numerical value obtained from a single gauge as a completed measurement.

Measurement requires that the same readout quantity be reconstructed across multiple gauges.

The measurement target is therefore not an individual gauge value, but a quantity that remains stable across a gauge family.

---

## 5. Multigauge Interference Readout

### 5.1 Spatial Phase-Gradient Readout

For a small displacement `h_chi`, the interference correlation in the spatial phase direction is read as

```math
p_{\mathrm{read}}
=
\frac{\arg\left(O(\chi+h_\chi)/O(\chi-h_\chi)\right)}
{2h_\chi}.
```

This is not standard momentum itself. It is a correlation gradient in the spatial phase direction.

### 5.2 Temporal Phase-Gradient Readout

For a small displacement `h_tau`, the interference correlation in the temporal phase direction is read as

```math
E_{\mathrm{read}}
=
\frac{\arg\left(O(\tau+h_\tau)/O(\tau-h_\tau)\right)}
{2h_\tau}.
```

This is not standard energy itself. It is a correlation gradient in the temporal phase direction.

### 5.3 R Readout

The readout amplitude is calibrated for each gauge and read as an amplitude square:

```math
R_{\mathrm{read}}
=
\gamma_g A_{\mathrm{read}}^2.
```

Here `gamma_g` is the readout gain of the gauge.

`R_read` is required to remain stable when multiple gauges are varied.

### 5.4 t/R Separation

`t` and `R` are not named axes assumed in advance.

A component read as strongly varying and continuous is labeled `t`, while a component read as weakly varying and stable is labeled `R`.

This paper uses the working indicator

```math
\frac{\operatorname{Var}(R_{\mathrm{read}})}
{\operatorname{Var}(t_{\mathrm{read}})}
```

to evaluate `t/R` separation.

When this ratio is sufficiently small, `R` can be treated as a stable readout separated from temporal variation.

---

## 6. Basic Readout in Symmetric ABC Collision

### 6.1 Single Collision

The first experiment uses the equal-amplitude condition

```text
A_A=A_B=1
```

and performs a single ABC collision.

The execution script is:

```text
run_abc_multigauge_interference_readout_v1.py
```

The results were:

| Quantity | Value |
|---|---:|
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `2.2315482794965646e-14` |
| `R_max_abs_error` | `4.440892098500626e-16` |
| `p_reflection_error_A` | `3.3306690738754696e-16` |
| `p_reflection_error_B` | `3.3306690738754696e-16` |
| `E_preservation_error_A` | `0.0` |
| `E_preservation_error_B` | `0.0` |
| `R_preservation_error_A` | `0.0` |
| `R_preservation_error_B` | `0.0` |
| `separation_ratio_time` | `3.6838616474030686e-30` |

The verdict was:

```text
multigauge_measurement_valid: true
single_gauge_only_used: false
```

### 6.2 Multiple Collisions

The same readout mechanism was then applied to eight AB collisions.

The execution script is:

```text
run_abc_multigauge_interference_readout_multi_collision_v1.py
```

The results were:

| Quantity | Value |
|---|---:|
| `ab_collision_count` | `8` |
| `wall_reflection_count` | `7` |
| `p_max_abs_error` | `2.5202062658991053e-14` |
| `E_max_abs_error` | `3.341771304121721e-13` |
| `R_max_abs_error` | `5.639932965095795e-14` |
| `max_p_reflection_error` | `4.440892098500626e-16` |
| `max_E_preservation_error` | `0.0` |
| `max_R_preservation_error` | `0.0` |
| `separation_ratio_time` | `2.7083289874897587e-28` |

The verdict was:

```text
multi_collision_multigauge_valid: true
single_gauge_only_used: false
```

### 6.3 Readout-Gauge Robustness

Readout robustness was tested with five gauge families, varying readout centers, widths, phases, and gains across a total of 130 gauges.

The execution script is:

```text
run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

The results were:

| Quantity | Value |
|---|---:|
| `case_count` | `5` |
| `total_gauge_count` | `130` |
| `max_p_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_E_abs_error_all_cases` | `3.0331293032759277e-13` |
| `max_R_abs_error_all_cases` | `1.5765166949677223e-14` |
| `max_R_gauge_std_all_cases` | `5.288392122597181e-15` |
| `max_separation_ratio_time_all_cases` | `1.3338999651354898e-27` |

The verdict was:

```text
robustness_sweep_valid: true
single_gauge_only_used: false
```

---

## 7. Failure Diagnosis of Simple Reversal Under Asymmetric R

### 7.1 Problem Setting

Under equal-amplitude conditions, the simple reversal

```text
p_A -> -p_A
p_B -> -p_B
```

for `p_A=+1`, `p_B=-1` works as a reflection readout.

However, when `R_A` and `R_B` differ, the same simple reversal does not necessarily preserve

```math
R_A p_A+R_B p_B.
```

Therefore, under asymmetric `R` conditions, one must diagnose whether simple reversal breaks the conserved readout.

### 7.2 Results

The execution script is:

```text
run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

The results were:

| Quantity | Value |
|---|---:|
| `case_count` | `8` |
| `asymmetric_case_count` | `7` |
| `individual_multigauge_valid_all_cases` | `true` |
| `weighted_energy_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `equal_case_weighted_momentum_preserved` | `true` |
| `asymmetric_cases_detect_weighted_momentum_failure` | `true` |
| `max_weighted_p_collision_error` | `16.000000000000036` |

These results detect that, under asymmetric `R` conditions, simple reversal breaks `R*p` conservation.

This is not a failure of the model. It is a diagnostic showing that the generalized map in the next section is required.

---

## 8. R-Weighted Generalized Collision Map

### 8.1 Conservation Conditions

Under asymmetric `R` conditions, the conserved readouts are required to be

```math
P_R
=
R_Ap_A+R_Bp_B
```

and

```math
K_R
=
R_Ap_A^2+R_Bp_B^2.
```

Here `P_R` is a momentum-like readout, and `K_R` is a conserved readout of squared phase gradient.

They are not identified with standard physical momentum or energy.

### 8.2 Map

The map that preserves `P_R` and `K_R` while reversing the relative phase gradient is

```math
p_A'
=
\frac{R_A-R_B}{R_A+R_B}p_A
+
\frac{2R_B}{R_A+R_B}p_B
```

and

```math
p_B'
=
\frac{2R_A}{R_A+R_B}p_A
+
\frac{R_B-R_A}{R_A+R_B}p_B.
```

This map does not substitute result values such as `R,T`. It is a local map used to verify whether the conservation relations hold for `R_read` and `p_read` read from multigauge interference.

What is confirmed in this paper is that this conservation-condition-based local map preserves `R*p` and `R*p^2` on the multigauge interference readouts. Whether this map itself is internally generated by an exchange-interference mechanism remains a subsequent question connecting to the preceding paper.

### 8.3 Verification

The execution script is:

```text
run_abc_multigauge_generalized_elastic_collision_readout_v1.py
```

Eight amplitude conditions were verified.

The results were:

| Quantity | Value |
|---|---:|
| `case_count` | `8` |
| `individual_readout_valid_all_cases` | `true` |
| `generalized_P_R_preserved_all_cases` | `true` |
| `generalized_K_R_phase_preserved_all_cases` | `true` |
| `E_tau_R_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `max_P_R_conservation_error` | `2.3803181647963356e-13` |
| `max_K_R_phase_conservation_error` | `1.4086509736443986e-12` |

The verdict was:

```text
generalized_elastic_collision_readout_valid: true
single_gauge_only_used: false
```

---

## 9. Non-Unit and Asymmetric Phase-Gradient Sweep

### 9.1 Purpose

The preceding section mainly used `p_A=+1`, `p_B=-1`.

This section uses non-unit and asymmetric initial phase gradients and includes both head-on and same-direction catch-up collision conditions.

### 9.2 Results

The execution script is:

```text
run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

Nine initial conditions were verified.

| Quantity | Value |
|---|---:|
| `case_count` | `9` |
| `collision_reached_all_cases` | `true` |
| `P_R_preserved_all_cases` | `true` |
| `K_R_phase_preserved_all_cases` | `true` |
| `relative_gradient_flipped_all_cases` | `true` |
| `E_tau_R_preserved_all_cases` | `true` |
| `R_total_preserved_all_cases` | `true` |
| `max_P_R_conservation_error` | `2.8910207561239076e-13` |
| `max_K_R_phase_conservation_error` | `1.5258905250448151e-12` |
| `max_relative_flip_error` | `1.2434497875801753e-14` |

The verdict was:

```text
velocity_sweep_generalized_collision_valid: true
single_gauge_only_used: false
```

---

## 10. Multiple Collisions of the Generalized Map

### 10.1 Purpose

This section tests whether the generalized collision map that holds for a single collision remains valid under repetition.

Wall reflection is an auxiliary condition used only to make the particles encounter each other again. Conservation is judged only immediately before and after each AB collision.

### 10.2 Results

The execution script is:

```text
run_abc_multigauge_generalized_elastic_collision_multi_collision_v1.py
```

Four conditions were executed, each with six AB collisions.

| Quantity | Value |
|---|---:|
| `case_count` | `4` |
| `max_ab_collision_count` | `6` |
| `max_wall_reflection_count` | `8` |
| `P_R_preserved_each_collision_all_cases` | `true` |
| `K_R_preserved_each_collision_all_cases` | `true` |
| `relative_flip_each_collision_all_cases` | `true` |
| `E_tau_R_preserved_each_collision_all_cases` | `true` |
| `R_preserved_each_collision_all_cases` | `true` |
| `max_P_R_error` | `3.552713678800501e-14` |
| `max_K_R_error` | `1.056932319443149e-13` |
| `max_relative_flip_error` | `2.220446049250313e-14` |

The verdict was:

```text
generalized_multi_collision_valid: true
single_gauge_only_used: false
```

---

## 11. Readout-Noise Robustness

### 11.1 Purpose

Measurement requires multigauge interference rather than a single gauge value.

Therefore, if the readout side fluctuates, zero-mean gauge fluctuations should be averaged out, while common bias across all gauges should be detected.

This section adds deterministic pseudo-noise to readout rows after state simulation.

This does not claim measurement validity under arbitrary noise. It is a controlled experiment testing whether zero-mean gauge fluctuations are canceled by multigauge averaging and whether common gauge bias is detected.

### 11.2 Results

The execution script is:

```text
run_abc_multigauge_generalized_elastic_collision_noise_robustness_v1.py
```

Four cases, two noise modes, and six noise levels were tested.

| Quantity | Value |
|---|---:|
| `case_count` | `4` |
| `noise_mode_count` | `2` |
| `noise_level_count` | `6` |
| `max_gauge_count` | `116` |
| `zero_mean_multigauge_valid_all` | `true` |
| `common_bias_detection_floor` | `1e-10` |
| `common_bias_detected_all_above_floor` | `true` |
| `zero_mean_max_p_mean_abs_error` | `1.4477308241112041e-13` |
| `zero_mean_max_R_mean_abs_error` | `3.552713678800501e-14` |
| `zero_mean_max_K_R_error` | `8.322231792590173e-13` |

The verdict was:

```text
noise_robustness_valid: true
single_gauge_only_used: false
```

---

## 12. Extreme R-Ratio Sweep

### 12.1 Purpose

If `R_read` behaves as a mass-like quantity, conserved readouts must be checked under strongly asymmetric `R` ratios.

This section sweeps from

```text
R_B/R_A = 0.015625
```

to

```text
R_B/R_A = 64.0.
```

### 12.2 Results

The execution script is:

```text
run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v1.py
```

The results were:

| Quantity | Value |
|---|---:|
| `case_count` | `12` |
| `max_R_dynamic_range` | `64.0` |
| `min_R_ratio_B_over_A` | `0.015625` |
| `max_R_ratio_B_over_A` | `64.0` |
| `P_R_preserved_all_cases` | `true` |
| `K_R_phase_preserved_all_cases` | `true` |
| `relative_gradient_flipped_all_cases` | `true` |
| `max_P_R_conservation_error` | `6.465938895416912e-13` |
| `max_K_R_phase_conservation_error` | `1.2789769243681803e-12` |
| `max_relative_flip_error` | `2.6867397195928788e-14` |

The verdict was:

```text
extreme_R_sweep_valid: true
single_gauge_only_used: false
```

---

## 13. Integration Summary

Nine experiments were integrated.

The execution script is:

```text
run_abc_multigauge_readout_integration_summary_v1.py
```

The integrated verdict was:

```text
experiment_count: 9
all_experiments_valid: true
single_gauge_only_used_any: false
integration_summary_valid: true
```

The experiment list is:

| Experiment | Purpose | Valid |
|---|---|---|
| `single_collision_multigauge_readout` | Read `p/E/R` by multigauge interference in a single ABC collision | `true` |
| `multi_collision_multigauge_readout` | Maintain `p/E/R` readout across repeated symmetric ABC collisions | `true` |
| `readout_robustness_sweep` | Test stable `p/E/R` reconstruction across readout-gauge configurations | `true` |
| `asymmetric_amplitude_diagnostic` | Detect that simple reversal breaks conservation under asymmetric R | `true` |
| `generalized_elastic_collision_readout` | Read a generalized map conserving `R*p` and `R*p^2` under asymmetric R | `true` |
| `generalized_velocity_sweep` | Verify the generalized map under non-unit and asymmetric phase gradients | `true` |
| `generalized_multi_collision` | Repeatedly apply the generalized map to multiple AB collisions | `true` |
| `generalized_noise_robustness` | Confirm cancellation of zero-mean readout noise and detection of common bias | `true` |
| `generalized_extreme_R_sweep` | Test the generalized map and readout under extreme R ratios | `true` |

---

## 14. Evaluation

The results are classified in the style of the exploratory physicist role.

| Target | Classification | Verdict |
|---|---|---|
| `p_read` is reconstructed from multigauge interference | Numerically constructed consequence | Retained |
| `E_read` is reconstructed from multigauge interference | Numerically constructed consequence | Retained |
| `R_read` is read as a stable residual across multiple gauges | Numerically constructed consequence | Retained |
| Symmetric collision reverses `p` and preserves `E,R` | Numerically constructed consequence | Retained |
| Under asymmetric R, simple reversal breaks `R*p` conservation | Numerically constructed consequence | Retained |
| The `R*p`, `R*p^2` conserving map holds under asymmetric R | Numerically constructed consequence | Retained |
| Zero-mean readout-noise components are averaged out by multigauge averaging | Numerically constructed consequence | Retained |
| Common readout bias is detected | Numerically constructed consequence | Retained |
| Identity with standard mass, standard momentum, and standard energy | Correspondence-map task | Not claimed |

---

## 15. Discussion

### 15.1 Mass, Momentum, and Energy Are Not Placed First

This paper does not place mass, momentum, and energy first.

The starting point is the closed phase system, local waves, observer wave, gauge family, and interference correlation.

From these, the quantities

```text
p_read
E_read
R_read
```

are read.

Thus the order of the paper is not

```text
Assume mass, momentum, and energy.
```

Instead, it is

```text
Test whether conserved quantities behaving as mass-like, momentum-like, and energy-like readouts emerge from interference readout.
```

### 15.2 R Is Difficult to Measure

`p_read` is a spatial phase gradient and is relatively easy to read through sign reversal or relative-gradient reversal.

`E_read` is a temporal phase-direction gradient and can be read by changing temporal windows.

In contrast, `R_read` is a stable residual.

It is stable precisely because its variation is small. For the same reason, it is difficult to measure.

The `t/R` separation in this paper is a working indicator used to test this point numerically.

### 15.3 A Single Gauge Is Not Enough

All experiments in this paper reported

```text
single_gauge_only_used: false
```

This is important.

This result does not prove the general impossibility of single-gauge measurement. It shows that, when the measurement condition is set as multigauge interference reconstruction, that condition alone consistently reads `p,E,R` and `R`-weighted conserved quantities.

`p,E,R` are not quantities directly visible as single local values. They are quantities stably reconstructed through multiple gauges.

Therefore, mass-like, momentum-like, and energy-like readouts must be treated as multigauge interference reconstructions rather than single-gauge measurements.

### 15.4 From Simple Reversal to a Generalized Conservation Map

Under equal-amplitude conditions, simple direction reversal appears conservative.

Under asymmetric `R` conditions, however, simple reversal breaks `R*p` conservation.

This break is not a failure of the model. Rather, if `R` is read mass-like, it diagnoses the need to generalize the conservation map to an `R`-weighted one.

This paper confirmed that the generalized map preserves `R*p` and `R*p^2`.

The generalized map here is a local map constructed from conservation conditions. The result of this paper is therefore the construction and verification of a readout map that behaves conservatively under asymmetric `R` conditions. Its internal interference generation is separated as the next question.

### 15.5 Connection to Standard Theory Is the Next Task

This construction does not directly derive standard physical momentum, energy, or mass.

However, the following structure was obtained:

```text
Spatial phase gradient -> momentum-like readout
Temporal phase gradient -> energy-like readout
Stable amplitude-square residual -> mass-like readout
R*p conservation -> relation that looks like momentum-like conservation
R*p^2 conservation -> relation that looks like squared-quantity conservation
```

Therefore, the readout-side foundation for constructing a correspondence map to standard theory has been obtained.

---

## 16. Conclusion

This paper numerically tested whether conserved readouts corresponding to mass-like, momentum-like, and energy-like quantities can be constructed from multigauge interference in an ABC closed phase system.

In a single ABC collision, `p_read`, `E_read`, and `R_read` were reconstructed from multiple gauges. `p` was reversed, and `E,R` were preserved. The maximum readout errors were `2.5202062658991053e-14` for `p`, `2.2315482794965646e-14` for `E`, and `4.440892098500626e-16` for `R`.

Across eight repeated symmetric collisions, `p` reversal, `E,R` preservation, internal label-mode preservation, and compensated closure were maintained.

Under asymmetric `R` conditions, simple `q -> -q` reversal was detected to break `R*p` conservation. A generalized collision map conserving `R*p` and `R*p^2` was then constructed. Across eight asymmetric amplitude cases, conservation was confirmed with maximum `R*p` error `2.3803181647963356e-13` and maximum `R*p^2` error `1.4086509736443986e-12`.

Additional verification covered non-unit and asymmetric initial phase gradients, same-direction catch-up collisions, multiple collisions, readout noise, and extreme `R`-ratio sweeps. The integration summary reported all nine experiments as `valid`, with no single-gauge-only judgment.

Within the numerical constructive scope of this paper, the ABC closed phase system consistently constructs conserved readouts corresponding to

```text
p_read
E_read
R_read
R*p
R*p^2
```

through multigauge interference.

This is not an identification with standard physical quantities. It is, however, a numerical construction for treating mass-like, momentum-like, and energy-like quantities as conserved readouts from a closed phase system rather than as primitive external entities.

---

# Appendix A. Executed Programs and Outputs

## A.1 Single ABC Multigauge Readout

```text
python3 run_abc_multigauge_interference_readout_v1.py
```

Output:

```text
abc_multigauge_interference_readout_result_v1/
```

Main files:

| Type | File |
|---|---|
| Report | [abc_multigauge_interference_readout_report_v1.md](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_report_v1.md) |
| JSON | [abc_multigauge_interference_readout_result_v1.json](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_result_v1.json) |
| Gauge CSV | [abc_multigauge_interference_readout_gauge_sweep_v1.csv](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_gauge_sweep_v1.csv) |
| p/E/R figure | [abc_multigauge_interference_readout_invariants_v1.png](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_invariants_v1.png) |
| t/R separation figure | [abc_multigauge_interference_readout_tr_separation_v1.png](abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_tr_separation_v1.png) |

## A.2 Symmetric Multiple Collisions

```text
python3 run_abc_multigauge_interference_readout_multi_collision_v1.py
```

Output:

```text
abc_multigauge_interference_readout_multi_collision_result_v1/
```

## A.3 Readout-Gauge Robustness

```text
python3 run_abc_multigauge_interference_readout_robustness_sweep_v1.py
```

Output:

```text
abc_multigauge_interference_readout_robustness_sweep_result_v1/
```

## A.4 Asymmetric Amplitude Diagnostic

```text
python3 run_abc_multigauge_interference_readout_asymmetric_amplitude_sweep_v1.py
```

Output:

```text
abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/
```

## A.5 Generalized Elastic Collision Map

```text
python3 run_abc_multigauge_generalized_elastic_collision_readout_v1.py
```

Output:

```text
abc_multigauge_generalized_elastic_collision_readout_result_v1/
```

## A.6 Asymmetric Velocity Sweep

```text
python3 run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1.py
```

Output:

```text
abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/
```

## A.7 Generalized Multiple Collisions

```text
python3 run_abc_multigauge_generalized_elastic_collision_multi_collision_v1.py
```

Output:

```text
abc_multigauge_generalized_elastic_collision_multi_collision_result_v1/
```

## A.8 Readout-Noise Robustness

```text
python3 run_abc_multigauge_generalized_elastic_collision_noise_robustness_v1.py
```

Output:

```text
abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1/
```

## A.9 Extreme R-Ratio Sweep

```text
python3 run_abc_multigauge_generalized_elastic_collision_extreme_R_sweep_v1.py
```

Output:

```text
abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v1/
```

## A.10 Integration Summary

```text
python3 run_abc_multigauge_readout_integration_summary_v1.py
```

Output:

```text
abc_multigauge_readout_integration_summary_result_v1/
```

---

# Appendix B. Executed Verification Notes

| Type | File |
|---|---|
| Definitional supplement | [全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md](全正符号ゼロ閉鎖の読出し多重性に関する定義補足.md) |
| Specification policy | [現在チャットメモ_多ゲージ干渉読出し仕様方針.md](現在チャットメモ_多ゲージ干渉読出し仕様方針.md) |
| Single readout verification | [ABC完全弾性衝突における多ゲージ干渉読出し数値検証メモ_v1.md](ABC完全弾性衝突における多ゲージ干渉読出し数値検証メモ_v1.md) |
| Symmetric repeated collision | [ABC多ゲージ干渉読出しの複数回衝突検証メモ_v1.md](ABC多ゲージ干渉読出しの複数回衝突検証メモ_v1.md) |
| Readout-gauge robustness | [ABC多ゲージ干渉読出しの読出し器頑健性スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しの読出し器頑健性スイープ検証メモ_v1.md) |
| Asymmetric amplitude diagnostic | [ABC多ゲージ干渉読出しの非対称振幅診断スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しの非対称振幅診断スイープ検証メモ_v1.md) |
| Generalized map | [ABC多ゲージ干渉読出しによる一般化弾性衝突写像検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突写像検証メモ_v1.md) |
| Asymmetric velocity | [ABC多ゲージ干渉読出しによる一般化弾性衝突の非対称速度スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の非対称速度スイープ検証メモ_v1.md) |
| Generalized multiple collision | [ABC多ゲージ干渉読出しによる一般化弾性衝突の複数回衝突検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の複数回衝突検証メモ_v1.md) |
| Noise robustness | [ABC多ゲージ干渉読出しによる一般化弾性衝突の読出しノイズ頑健性検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の読出しノイズ頑健性検証メモ_v1.md) |
| Extreme R ratio | [ABC多ゲージ干渉読出しによる一般化弾性衝突の極端R比スイープ検証メモ_v1.md](ABC多ゲージ干渉読出しによる一般化弾性衝突の極端R比スイープ検証メモ_v1.md) |
| Integration summary | [ABC多ゲージ干渉読出し実験群の統合サマリー_v1.md](ABC多ゲージ干渉読出し実験群の統合サマリー_v1.md) |

---

# References

## Self-Citations

1. Noriaki Kihara, "Basic Axiom System v2 of the Nameless Equal-Amplitude Composite Wave Model," 2026-07-10.
2. Noriaki Kihara, "Definitional Supplement on Readout Multiplicity of the All-Positive Zero Closure," 2026-07-11.
3. Noriaki Kihara, "Construction Experiment of Complete Elastic Reflection of Fermionic Bilocal Waves in a Closed Phase System Without Assuming Background Space," Version DOI: `10.5281/zenodo.21291020`, Concept DOI: `10.5281/zenodo.21291018`, 2026.
4. Noriaki Kihara, "Interference Construction of a Perfect Reflection Map from a Fermionic Inverse-Phase Core," Version DOI: `10.5281/zenodo.21295480`, Concept DOI: `10.5281/zenodo.21295479`, 2026.
5. Noriaki Kihara, "Curvature Renormalization and Perfect-Reflection Stability by Curved Closed Stationary Waves," Version DOI: `10.5281/zenodo.21304040`, Concept DOI: `10.5281/zenodo.21304039`, 2026.

## External References

External references are used only as minimal background for conservation laws, interference phase, and geometric phase. They are not used as derivational grounds for this paper.

6. E. Noether, "Invariante Variationsprobleme," *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 235-257, 1918.
7. S. Pancharatnam, "Generalized theory of interference, and its applications," *Proceedings of the Indian Academy of Sciences A*, 44, 247-262, 1956.
8. M. V. Berry, "Quantal phase factors accompanying adiabatic changes," *Proceedings of the Royal Society of London A*, 392, 45-57, 1984. DOI: `10.1098/rspa.1984.0023`.
