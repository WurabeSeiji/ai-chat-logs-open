# Preliminary Summary of Distance-Exponent Readout by Independent Metric C and Relational Compensation Decomposition in an ABC Closed Phase System

**Date:** 2026-07-12
**Author:** Noriaki Kihara
**Position:** Summary of the ABC distance-exponent preliminary experiment group in the Wave Information Readout series
**Version DOI:** 10.5281/zenodo.21318701
**Concept DOI:** 10.5281/zenodo.21318700

---

## Abstract

This summary consolidates a group of preliminary ABC closed phase system experiments introduced after the AB two-body system failed to determine the distance exponent independently.

The target experiments are:

1. distance-exponent readout of AB relational compensation using an independent metric C,
2. C-position control for distance-exponent readout, and
3. relational compensation decomposition.

The conclusion is:

```text
To address the distance-exponent question left open by the AB two-body problem,
the measuring wave C was placed inside the same closed phase system.

Even in this ABC three-body system, the acceleration-like AB relational compensation can be read.

However, no inverse-law or inverse-square-law dependence on position phase difference was observed.
```

In valid gauge cases, the recovered exponent was mainly

```text
alpha ≈ -1.
```

Since the experiments classify the readout as

```text
I_AB(L) proportional to L^{-alpha},
```

`alpha≈-1` means

```text
I_AB(L) proportional to L.
```

That is, the result is proportional-type.

This is interpreted as a consequence of the fact that the present experiment is essentially a one-dimensional harmonic-oscillation model made of two localized waves imitating localized particles.

Within this model, effects of relative distance are not received, or even if they are received, they cannot be observed.

In particular, although C was introduced as an independent metric gauge, C is itself part of the closed phase system and therefore cannot be an absolute external gauge.

Thus, even after adding C, the proportional type remained within the current model.

This does not generally deny inverse-square behavior.

The narrower claim that was not supported is:

```text
An inverse-law or inverse-square-law relation appears natively
from the present one-dimensional AB harmonic model plus the C metric gauge.
```

The possibility of extending the position phase degrees of freedom to three or four dimensions was considered. However, as long as the same localized two-wave harmonic model is used, any higher-dimensional configuration returns to the same structure when projected onto a geodesic or two-dimensional section. Therefore, this experiment group did not proceed to three- or four-dimensional extensions at this stage.

**Keywords:** ABC closed phase system, independent metric C, distance exponent, relational compensation decomposition, `f_AB`, `f_AC`, `f_BC`, `f_ABC`, harmonic readout, inverse-square non-detection

---

## 1. Experiment List

| No. | Experiment | Purpose | Main result |
|---:|---|---|---|
| 1 | C-gauge valid-window test | Test whether `C` can serve as an independent metric for AB distance exponent | Valid window exists |
| 2 | C-position control | Test whether symmetric/asymmetric/phase-inverted C changes the exponent | Proportional type remained |
| 3 | Relational compensation decomposition | Separate `f_AB`, `f_AC`, `f_BC`, and `f_ABC` | Separation rule supported |

---

## 2. Distance-Exponent Readout Using Independent Metric C

### 2.1 Meaning of the Experiment

In the AB two-body system alone, it was possible to read that

```text
relative distance changed.
```

However, there was no independent gauge for deciding whether the change should be measured as proportional, inverse proportional, or inverse-square.

Therefore, a third wave `C` was placed and tested as a metric gauge for reading the position-change and temporal-phase-change quantities of AB.

The following were not implemented:

```text
1/L_AB
1/L_AB^2
1/A_chi_tau
F = G m_A m_B / L_AB^2
```

Thus, if an inverse-law or inverse-square-law form appears, it must be a result of C-gauge readout rather than a constructed implementation.

### 2.2 Main Result

| Quantity | Value |
|---|---:|
| `config_count` | `160` |
| `gauge_valid_count` | `41` |
| `resolution_valid_count` | `60` |
| `clock_valid_count` | `100` |
| `disturbance_valid_count` | `101` |
| `inverse_like_alpha1_count` | `0` |
| `inverse_square_like_alpha2_count` | `0` |
| `proportional_like_alpha_minus1_count` | `41` |
| `alpha_min` | `-1.0000000000000002` |
| `alpha_max` | `-0.9569628567496087` |

The classification by relational time was:

| Time candidate | Classification |
|---|---|
| `tau_ABC` | `proportional_like_alpha_minus1 = 41` |
| `tau_AB` | `other = 41` |
| `tau_AC` | `proportional_like_alpha_minus1 = 41` |
| `tau_BC` | `proportional_like_alpha_minus1 = 41` |

`tau_AB` deviated from `tau_ABC` by up to about `0.659`.

Therefore, reading the distance exponent requires explicitly specifying

```text
what is read as time.
```

### 2.3 C-Gauge Valid Window

In this experiment, `R_C` was also treated as the metric-cell density of C.

That is,

```text
small R_C
= low frequency of C
= long wavelength of C
= wide cell interval of C.
```

If `R_C` is too small, C cannot read the AB position-change quantity across multiple cells.

If `R_C` becomes large, resolution improves, but the C-derived relations `f_AC`, `f_BC`, and `f_ABC` contaminate the main AB readout.

Therefore, C has a valid window:

```text
not too coarse, not too strong.
```

### 2.4 Figures

| C-gauge eligibility | Alpha candidates |
|---|---|
| <img src="abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_gauge_eligibility_map_v1.png" width="520"> | <img src="abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_alpha_candidates_v1.png" width="520"> |

| Reference curve | Filter diagnosis |
|---|---|
| <img src="abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_reference_curve_v1.png" width="520"> | <img src="abc_c_gauge_ab_distance_exponent_preliminary_result_v1/abc_c_gauge_ab_distance_exponent_validity_filters_v1.png" width="520"> |

---

## 3. C-Position Control

### 3.1 Meaning of the Experiment

The third wave `C` is not an external observer. It is part of the closed phase system.

Thus, placing C simultaneously generates

```text
f_AC
f_BC
f_ABC.
```

The experiment changed the position of C to test whether

```text
the distance exponent changes,
or only the valid gauge window narrows due to contamination from C.
```

### 3.2 C-Position Modes

The compared positions were:

```text
symmetric
symmetric_pi_flip
a_side_small
b_side_small
a_side_large
b_side_large
a_side_large_pi_flip
b_side_large_pi_flip
```

### 3.3 Main Result

| Quantity | Value |
|---|---:|
| `config_count` | `160` |
| `gauge_valid_count` | `70` |
| `gauge_valid_position_pair_count` | `23` |
| `max_pair_abs_alpha_difference` | `0.11776992863447344` |
| Main classification under `tau_ABC` | `proportional_like_alpha_minus1 = 70` |

The number of valid gauges by C-position mode was:

| C position | Valid count |
|---|---:|
| `symmetric` | `12` |
| `symmetric_pi_flip` | `12` |
| `a_side_small` | `11` |
| `b_side_small` | `11` |
| `a_side_large` | `6` |
| `b_side_large` | `6` |
| `a_side_large_pi_flip` | `6` |
| `b_side_large_pi_flip` | `6` |

Large asymmetric C placement reduced the number of valid cases.

It is safer to read this not as the appearance of a new stable distance exponent, but as:

```text
bias in f_AC and f_BC narrowed the validity conditions of the C gauge.
```

### 3.4 Figures

| Valid count by C position | Alpha by C position |
|---|---|
| <img src="abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_valid_counts_v1.png" width="520"> | <img src="abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_alpha_candidates_v1.png" width="520"> |

| A-side/B-side symmetry |
|---|
| <img src="abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1/abc_c_gauge_c_position_pair_symmetry_v1.png" width="620"> |

---

## 4. Relational Compensation Decomposition

### 4.1 Meaning of the Experiment

In the ABC three-body system, the relations separate into at least four:

```text
f_AB
f_AC
f_BC
f_ABC
```

The main target is `f_AB`.

`f_AC` and `f_BC` are the C-derived two-body relations.

`f_ABC` is the common relation of the whole ABC system.

The experiment tested the following separation:

```text
representative time: tau_ABC
circumferential-direction candidates: f_AB, f_AC, f_BC
common-mode or central-direction candidate: f_ABC
```

Especially important is the discipline:

```text
Use f_ABC as representative time.
Do not directly add f_ABC to the AB circumferential direction.
```

### 4.2 Main Result

| Quantity | Value |
|---|---:|
| `config_count` | `240` |
| `decomposition_valid_count` | `180` |
| `AB_dominant_valid_count` | `48` |
| `non_AB_dominant_count` | `132` |
| `inverse_or_inverse_square_in_AB_dominant_count` | `0` |

The `tau_ABC` distance-exponent classification in AB-dominant cases was:

```text
proportional_like_alpha_minus1: 48
```

The native `f_AB` classification was:

```text
proportional_like_alpha_minus1: 180
```

The C-bias and C-asymmetry terms were read mainly as

```text
constant_like_alpha0.
```

This indicates that the C-derived projection bias did not appear as a stable inverse-law or inverse-square term with respect to AB distance `L_AB`.

### 4.3 Direct-Addition Control for f_ABC

In the erroneous direct-injection control where `f_ABC` was directly added to the AB circumferential direction, the result was:

```text
other: 63
proportional_like_alpha_minus1: 117
```

This shows that treating `f_ABC` as a causal term in the AB circumferential direction muddies the classification.

At this stage, therefore, the following separation is appropriate:

```text
Use f_ABC as representative time tau_ABC.
Do not directly add f_ABC to the AB circumferential direction.
```

This does not mean that `f_ABC` is ignored.

`f_ABC` is recorded as the common mode or central-direction compensation of the whole system.

### 4.4 Figures

| Decomposition valid count | Pair alpha |
|---|---|
| <img src="abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_valid_counts_v1.png" width="520"> | <img src="abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_pair_alpha_v1.png" width="520"> |

| Native / pair comparison | C-contamination diagnosis |
|---|---|
| <img src="abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_native_vs_pair_alpha_v1.png" width="520"> | <img src="abc_c_gauge_relation_decomposition_preliminary_result_v1/abc_c_gauge_relation_decomposition_contamination_v1.png" width="520"> |

---

## 5. Overall Judgment

### 5.1 Established

| Item | Judgment |
|---|---|
| C gauge has a valid window | Established |
| Too small `R_C` makes cells too coarse | Confirmed |
| Too strong or biased C causes gauge contamination | Confirmed |
| Discipline of reading `tau_ABC` as representative time | Supported |
| Need to record `tau_AB`, `tau_AC`, and `tau_BC` as auxiliary diagnostics | Confirmed |
| Separation of `f_AB`, `f_AC`, `f_BC`, and `f_ABC` | Supported |
| Discipline of not directly adding `f_ABC` in the circumferential direction | Supported |

### 5.2 Not Established

| Item | Judgment |
|---|---|
| Inverse-law form appears from C gauge alone | Not detected |
| Inverse-square form appears from C gauge alone | Not detected |
| Changing C position alone changes the exponent law | Not detected |
| C-derived terms become inverse-law or inverse-square terms of AB distance | Not detected |
| Present minimal AB harmonic model becomes standard-gravity type | Not established |

---

## 6. Interpretation

The most important consequence of this experiment group is:

```text
Even after returning to an ABC three-body system, the current minimal model is read as proportional-type.
```

This is not a negative failure.

It is a boundary condition: the acceleration-like readout obtained in the AB two-body system was re-read using the independent metric C, but it was not transformed into an inverse-law or inverse-square-law form.

The AB two-body system lacked an independent gauge.

The ABC three-body system introduced the independent metric C, but since C itself is part of the closed phase system,

```text
f_AC
f_BC
f_ABC
```

arise.

Even after separating these terms, AB-dominant cases preserved the proportional type.

Therefore, the current model essentially behaves as

```text
a one-dimensional relative-phase displacement model.
```

Within that scope, inverse-law and inverse-square-law forms do not arise naturally.

The result shows that the distance-exponent problem cannot be solved merely by adding C.

Increasing the strength of C improves metric resolution, but increases contamination by `f_AC`, `f_BC`, and `f_ABC`.

Weakening C reduces contamination, but the cell resolution becomes insufficient for reading changes in the position phase difference.

Thus, C gauge has a valid window, but that window was not a mechanism for converting the present one-dimensional harmonic model into an inverse-square form.

---

## 7. Deferred Problems

### 7.1 The Inverse Square Is Not Generally Denied

The present experiment denies only the narrow claim:

```text
If C gauge and relational decomposition are added to the present minimal AB harmonic model,
an inverse-law or inverse-square form appears natively.
```

This claim was not supported.

However, this does not show that inverse-square forms are absent from closed phase systems in general.

### 7.2 Reason for Not Performing Three- or Four-Dimensional Extensions

To explain the inverse-square form in a standard way, a structure such as area sweep or spherical-shell sweep seems necessary.

Therefore, experiments increasing position-phase degrees of freedom to three or four dimensions were considered.

However, the present experiment system is a model in which two localized waves harmonically oscillate in the relative phase direction.

Even if this structure is made higher-dimensional, projection onto a geodesic or two-dimensional section returns to the same two-body harmonic readout as in this experiment.

Therefore, within this model, three- or four-dimensional extension is unlikely to newly produce inverse-law or inverse-square behavior.

For this reason, the higher-dimensional extension experiment was not performed in this series.

### 7.3 Circular Waves, Spherical Waves, and Limited Connection to the Double-Slit Papers

In the observational principle of this series, however, what can be observed is a localized wave. A merely conceptual extended circular wave cannot simply be placed as a real entity.

If one assumes a circular or spherical-shell wave, it appears easy to construct an area sweep or inverse-square form.

However, such a wave is broadly extended in both spatial and temporal phase directions, and it becomes unclear by which local readout its position phase and temporal phase are measured.

Therefore, at present, the circular-wave or spherical-shell-wave model is suspended not because it is impossible to experiment with, but because it may contain an excessive reality assumption relative to the observational principle.

When referring here to the double-slit thought-experiment series, the scope of citation must be limited.

The first paper in that series showed that, when a position fluctuation `P(y)` is given to a stationary single-wavelength point light source, each far-field interference fringe retains its shape and shifts only by the geometric path difference; repeated trials push the distribution of fringe shifts forward from `P(y)`.

As a concrete example, when `P(y)=cos^2` is placed, the `cos^2` form is preserved under a near-axis linear map, and non-near-axis deviation appears as several-percent geometric nonlinearity.

However, that paper treats an observation mode in which the fringe shift is read in each trial. It is separate from the visibility-reduction mode that integrates many-trial intensities into one image.

The second paper in that series extended this pushed-forward readout to a localized odd-harmonic source.

There, a localized odd-harmonic source can preserve its shape and appear in the double-slit far field only under alignment conditions, and under position fluctuation the shape preservation becomes fragile due to off-axis scattering.

The single-wavelength case `N=1` has no such fragility and was confirmed as a special case that agrees with the first paper to machine precision.

Thus, the implication that can be correctly carried from the double-slit papers to the present paper is limited to the following research direction.

From this viewpoint, rather than placing an extended probability wave or spherical wave directly as real, it is closer to the observational principle of this series to examine whether

```text
many observations of localized waves,
initial-position or phase fluctuation,
and the push-forward of that fluctuation to the observation side
```

can appear as a broadened statistical image or area distribution.

This is not a proposition directly proven by the double-slit papers.

What those papers showed is only that a position-fluctuation distribution is pushed forward into a distribution of fringe shifts, and that extension to localized sources adds alignment conditions and fragility of shape preservation.

Therefore, instead of

```text
assuming an extended wave,
```

the next question should be:

```text
Can many-trial, many-initial-condition readout of localized waves appear statistically as area broadening?
```

That question exceeds the scope of this experiment series and is deferred here.

---

## 8. How This Series Should Be Closed

It is reasonable to close this series with the following results.

```text
In the AB two-body system, harmonic readout and chi-tau area can be read.
However, there is not enough independent gauge structure to determine the distance exponent.

In the ABC three-body system, C gauge and relational decomposition can be introduced.
However, in the present minimal model the proportional type is preserved,
and inverse-law or inverse-square forms do not appear natively.

Therefore, the inverse-square problem must be redesigned not as a simple extension
of this one-dimensional closed phase model, but as a problem of observational statistics,
area readout, localized-wave ensembles, and push-forward maps of initial fluctuations.
```

---

# References

## Self-References

1. Noriaki Kihara, "Basic Axiom System of the Nameless Equal-Amplitude Composite Wave Model v4", Version DOI: `10.5281/zenodo.21316620`, Concept DOI: `10.5281/zenodo.21315735`, 2026.
2. Noriaki Kihara, "Construction Experiment of Multigauge Interference Readout Conserved Quantities in an ABC Closed Phase System", Version DOI: `10.5281/zenodo.21308050`, Concept DOI: `10.5281/zenodo.21308049`, 2026.
3. Noriaki Kihara, "Preliminary Summary of Harmonic Readout and c=1 Area Sweep in an AB Two-Body Closed Phase System", Version DOI: `10.5281/zenodo.21318697`, Concept DOI: `10.5281/zenodo.21318696`, 2026.
4. Noriaki Kihara, [Experiment Specification for Distance-Exponent Readout of AB Relational Compensation Using Independent Metric C in an ABC Closed Phase System v1](ABC閉鎖位相系における独立計量CによるAB関係補償の距離指数読出し実験仕様書%20v1.md), 2026.
5. Noriaki Kihara, [Definitional Supplement on Internal Closure of Self-Terms and Separation of N-Body External Readout in Closed Complex Phase Waves](閉鎖複素位相波における自己項の内部閉鎖とN体外部読出し分離に関する定義補足.md), 2026.
6. Noriaki Kihara, "Double-Slit Interference Thought Experiment with a Source Having Position Fluctuation: Push-Forward of Source Position Distribution to Fringe-Shift Distribution (Shape Preservation)", Version DOI: `10.5281/zenodo.21035809`, Concept DOI: `10.5281/zenodo.21035808`, 2026.
7. Noriaki Kihara, "Double-Slit Interference Thought Experiment with a Localized Odd-Harmonic Source: Shape Preservation Is Conditional and Fragile, and the Single-Wavelength Case N=1 Is a Robust Special Case", Version DOI: `10.5281/zenodo.21035831`, Concept DOI: `10.5281/zenodo.21035830`, 2026.

## External References

The external references are used only as minimal standard background on distance laws, relativistic time, and phase readout, not as derivational premises of this paper.

8. Isaac Newton, *Philosophiae Naturalis Principia Mathematica*, 1687.
9. Albert Einstein, "Die Grundlage der allgemeinen Relativitaetstheorie", *Annalen der Physik* 49, 769-822, 1916. DOI: `10.1002/andp.19163540702`.
10. Y. Aharonov and D. Bohm, "Significance of electromagnetic potentials in the quantum theory", *Physical Review* 115, 485-491, 1959. DOI: `10.1103/PhysRev.115.485`.
