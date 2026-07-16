# Selection of the State-Exchange Weight G_R = 1 - R and Candidate Correspondence with the Fine-Structure Constant: A Numerical Experiment v1

**Date:** 2026-07-15  
**Author:** Noriaki Kihara  
**Position:** Wave Information Readout series, exchange-weight selection experiment summary  
**Version DOI:** 10.5281/zenodo.21396761  
**Concept DOI:** 10.5281/zenodo.21396760  

---

## Abstract

This paper investigates, by numerical experiment, under what conditions a fermion-like exchange-scattering map in a closed phase system selects a sharp value when the degree of state exchange is read through the representative quantity `G_R=1-R`. Here `R` is the reflection coefficient scanned in the numerical experiments; the main readout is not reflection itself, but the amount by which the state moves into the other channel.

The central subject is not the agreement with the fine-structure constant itself. The central subject is whether a closed exchange map constructed from a small set of axioms can spontaneously select a representative weight for the degree of state exchange.

There are two starting points.

First, in System A, the low-localization and harmonic-transfer readout experiment, a low-localization wave and a wave with higher harmonics were placed in a two-channel exchange-scattering map. The experiment observed that the harmonic structure moved into the other channel and that the low-localization wave became localized. This system is based on self-references [3] and [4], with the `R`-sweep specification for this paper defined in self-reference [7]. In that experiment, `R=0.70` appeared empirically as an effective scattering point.

Second, in System B, the white-cat, black-cat, and gray-cat metastable-interface experiment, the condition under which an AB two-component distribution enters a gray metastable phase was explored. This system is based on self-reference [5], with the `R`-sweep specification for this paper defined in self-reference [8]. When exchange scattering controlled by `R` is introduced, the gray-metastable landscape may develop sharp structures along the `R` axis.

In this paper, System A is used to test whether harmonic degrees of freedom are involved in the selection of the exchange weight. System B is then used to sweep `R` over a broad range and to read candidate bands of `G_R=1-R` as a representative quantity for state exchange.

In System A, a single base harmonic did not produce a sharp exchange-weight concentration. Once at least one additional harmonic degree of freedom was present, the evaluation points concentrated near `R=0.697177879...`, or equivalently near `G_R≈0.302822121`. This concentration was robust under odd harmonics, even harmonics, mixed parity packets, phase shifts, wavelength shifts, and extremely small harmonic amplitudes.

In System B, the gray-metastable readout was used to scan from `R=0.686602902...` to `R=0.702465367...` with `Delta R=1.0e-7`. Seven candidate bands were detected. The deepest band appeared at `R=0.697177902556148`. This agrees, within the effective numerical precision of the v5 experiment, with

```text
R_low = 0.697177879231003
```

obtained from the low-energy fine-structure constant `alpha(0)` by

```text
R = 1 - sqrt(4 pi alpha).
```

Converting this deepest band back to `1/alpha` gives `137.036020287643`, corresponding to the low-energy fine-structure constant `1/alpha(0)=137.035999177 +- 0.000000021`. The second candidate band appeared at `R=0.688363902556148`, corresponding to `1/alpha=129.394062925467`, which lies near the effective coupling at the Z-boson mass scale, `1/alpha(M_Z^2)=128.946 +- 0.015`.

This paper does not claim to derive the fine-structure constant. It claims the following numerical fact: in exchange-scattering maps of closed phase systems, when a metastable interface with harmonic degrees of freedom is investigated, candidate bands of exchange weight corresponding to two electromagnetic-coupling regions referenced in the standard theory appear without inserting `alpha` into the scan.

---

## 1. Motivation

### 1.1 System A: Low-Localization and Harmonic-Transfer Readout

System A is based on the low-localization and harmonic-transfer readout experiment described in self-references [3] and [4]. In this paper, that system is reused through the `R`-neighborhood uniform-sweep specification in self-reference [7].

In this system, an exchange-interference scattering matrix is used to repeatedly scatter an AB two-body wave in which one side contains higher harmonics.

The experiment is not a simple reflection test. A is placed as a low-localization wave, B is placed as a localized wave containing higher harmonics, and the map

```text
(A', B') = S_R (A, B)
```

is iterated. The readout asks whether the harmonic distribution of B moves into A and whether the localization of A increases.

The main evaluated quantities are:

```text
L        : localization index
N_eff    : effective harmonic order
B_to_A   : transfer of the initial B-side harmonic distribution into A
```

In this experiment, an effective localization-exchange point appeared near `R=0.70`, or equivalently near `G_R≈0.30`. The first question is therefore:

```text
Do harmonic degrees of freedom participate in the selection of the exchange weight G_R?
```

### 1.2 System B: White-Cat, Black-Cat, and Gray-Cat Metastable Interface

System B is based on the white-cat, black-cat, and gray-cat metastable-interface experiment in self-reference [5]. In this paper, that system is reused through the `R`-neighborhood sweep specification in self-reference [8].

In this system, the white-cat, black-cat, and gray-cat metastable interface is constructed as an AB two-component closed system.

The AB distribution is represented by complex amplitudes

```math
a,\quad b
```

and the probabilities

```math
p_A=|a|^2,\qquad p_B=|b|^2
```

are read through

```math
S=p_A-p_B.
```

The condition

```text
S ≈ 0
```

together with a small oscillation, but not a fully intrinsic gray phase, is classified as the gray metastable phase.

Introducing the exchange-scattering coefficient `R` into this system allows us to test at which exchange weight `G_R=1-R` the gray metastable phase tends to appear. The second question is therefore:

```text
Is R≈0.7, or equivalently G_R≈0.3, merely a tuning value?
Or is it a reaction point intrinsic to the gray metastable interface?
```

### 1.3 Connection Made in This Paper

System A corresponds to self-references [3], [4], and [7], and reads harmonic degrees of freedom and localization exchange.

System B corresponds to self-references [5] and [8], and reads the gray metastable interface and the exchange-weight landscape.

This paper connects the two:

```text
System A [3,4,7]:
  When harmonic degrees of freedom are present,
  candidate exchange-weight bands appear in localization exchange.

System B [5,8]:
  On the gray metastable interface,
  sweep R over a broad range and read candidate exchange-weight bands.
```

This connection tests under what closed-system conditions the representative state-exchange weight `G_R=1-R` is sharply selected.

---

## 2. Exchange-Scattering Coefficient R and Conversion to alpha

In this series, the reflection coefficient of the exchange-scattering matrix is denoted by `R`. However, the subject of this paper is not `R` itself.

The numerical experiments in this paper evaluate neither the reflection amplitude nor the transmission amplitude. They evaluate how much of the state is transferred into the other channel by exchange.

In System A, the transfer of the initial B-side harmonic distribution into A is read as `B_to_A`. In System B, entry into the gray metastable interface under repeated exchange scattering is read through `gray_error` and `gray_depth`.

Thus, the directly evaluated quantity is not the amount remaining on the reflection side, but the degree to which state exchange is established.

Therefore, using `R` as the sweep coordinate, this paper defines

```math
G_R := 1 - R
```

as the representative exchange weight used to compare that exchange-establishment degree along the `R` axis.

On the standard-theory side, in rationalized natural units, if the dimensionless electromagnetic coupling amplitude is denoted by `e`, then

```math
\alpha = \frac{e^2}{4\pi}.
```

Hence

```math
e = \sqrt{4\pi\alpha}.
```

This paper compares the representative exchange quantity `G_R` with the standard-theory dimensionless electromagnetic coupling amplitude `e`:

```math
G_R = e.
```

This gives

```math
1 - R = \sqrt{4\pi\alpha},
```

and therefore

```math
R = 1 - \sqrt{4\pi\alpha}.
```

This is not an arbitrary monotonic map inserted only to fit the number. It is based on the correspondence hypothesis that the degree of state exchange directly evaluated in the experiments is read through the representative quantity `G_R=1-R` on the `R` axis and compared with the standard-theory dimensionless coupling amplitude `sqrt(4 pi alpha)`.

The hierarchy relative to the implementation of the scattering matrix must be kept explicit.

In the implementation, the matrix amplitudes of the reflection and exchange channels are `r` and `t`, and

```text
R = |r|^2
T = |t|^2 = 1 - R
```

are used. Therefore the matrix amplitudes themselves are `|r|=sqrt(R)` and `|t|=sqrt(1-R)`.

This paper does not identify the standard scattering amplitude `|t|` itself with the standard-theory electromagnetic coupling amplitude `e`. If one set `|t|=e`, then

```math
\sqrt{1-R}=\sqrt{4\pi\alpha}
```

and the corresponding expression would be

```math
R = 1 - 4\pi\alpha.
```

The observed landscape in this paper appears sharply on the side of `G_R=1-R`, the representative quantity for the degree of state exchange, rather than on this `|t|` correspondence.

This paper does not claim to derive this correspondence from first principles. The unresolved question is why a closed-phase exchange-scattering map may select, not the matrix amplitude `|t|`, but the representative state-exchange weight `G_R=1-R` near the electromagnetic coupling amplitude of the standard theory. What this paper tests is how sharply this correspondence appears in the numerical experiments.

For the low-energy fine-structure constant, the CODATA 2022 recommended value is used:

```text
1/alpha(0) = 137.035999177 ± 0.000000021
R_low      = 0.697177879231003
```

For the effective coupling at the Z-boson mass scale, the value including hadronic vacuum polarization is used:

```text
1/alpha(M_Z^2) = 128.946 ± 0.015
R_MZ           = 0.687822933884774
```

The v5 numerical experiment used here is estimated to have an effective precision of roughly seven digits in `R`. Thus, for the low-energy side, the numerical precision of the experiment sets the comparison limit, while for the high-energy side the four- to five-digit effective precision of the standard-theory reference value sets the comparison limit.

---

## 3. System A [3,4,7]: Harmonic Degrees of Freedom and Localization Exchange

### 3.1 N-Series

First, the one-sided harmonic condition was tested directly.

The A side was set to the base wave `A=[1]`, and the B side was varied as

```text
B=[1], [2], [3], [5], [15], [63].
```

The same-order control `N=1` had no sharp decision point with respect to the exchange weight. In contrast, for the one-sided harmonic condition `N>=2`, the quantities `R_star_L`, `R_star_N`, `R_star_transfer`, and `R_star_joint` concentrated near `R=0.697177879128`. In terms of exchange weight, this corresponds to `G_R=1-R≈0.302822121`.

![Exchange-weight landscape of the System A N-series](system_A_localization_exchange_R_sweep_result_v1/odd_kernel_N_1_2_3_5_15_63_gap_depth_distribution_overview_v1.png)

**Figure 1.** System A N-series. In the single-base-harmonic control, the landscape along the `R` axis is broad. When an additional harmonic is placed on one side, the evaluation points concentrate near `R_137`, or `G_R≈sqrt(4πalpha(0))`.

Representative values are:

| Condition | `R_star_joint` | Reading |
|---|---:|---|
| `A=[1], B=[1]` | `0.6` | Same-order control. No sharp exchange-weight concentration |
| `A=[1], B=[2]` | `0.697177879128` | Additional harmonic present |
| `A=[1], B=[3]` | `0.697177879128` | Additional harmonic present |
| `A=[1], B=[5]` | `0.697177879128` | Additional harmonic present |
| `A=[1], B=[15]` | `0.697177879128` | Additional harmonic present |
| `A=[1], B=[63]` | `0.697177879128` | Additional harmonic present |

The result indicates that the key condition is not odd harmonics themselves, but the existence of at least one additional harmonic degree of freedom.

### 3.2 Odd, Even, and Mixed Harmonic Packets

Next, the B-side harmonic packet was varied.

```text
even packet:
  B=[1,2,4,6]

alternating packet 1:
  B=[1,2,3,4,5]

alternating packet 2:
  B=[1,3,4,5,6]
```

All cases concentrated at `R_star_joint=0.697177879128`.

| Condition | `R_star_joint` |
|---|---:|
| even packet `B=[1,2,4,6]` | `0.697177879128` |
| mixed packet `B=[1,2,3,4,5]` | `0.697177879128` |
| mixed packet `B=[1,3,4,5,6]` | `0.697177879128` |

Thus, this exchange-weight concentration is not a special property of odd harmonics.

### 3.3 Insensitivity to Amplitude, Phase, and Wavelength

The amplitude, phase, and wavelength of the additional B-side harmonic were then varied.

For the amplitude test, `R_star_joint` was maintained even when the additional harmonic weight was strongly reduced.

| additional harmonic weight | `R_star_joint` | Reading |
|---:|---:|---|
| `0.5` | `0.697177879128` | maintained |
| `0.1` | `0.697177879128` | maintained |
| `0.001` | `0.697177879128` | maintained |
| `0.0001` | `0.697177879128` | maintained |
| `0.000001` | `0.697177879128` | some indicators fluctuate, but joint value maintained |
| `0` | `0.6` | no exchange-weight concentration |

When the phase was shifted by 10% or 30%, `R_star_joint=0.697177879128` was maintained.

When the wavelength was shifted by 3%, 10%, or 30%, the same `R_star_joint` was maintained.

This shows that the exchange-weight concentration is not determined only by harmonic amplitude, phase alignment, or exact wavelength matching.

The reading adopted in this paper is:

```text
The concentration of the exchange weight does not require oddness.
It does not require exact amplitude matching.
It does not require exact phase matching.
It does not require exact wavelength matching.

The presence of at least one additional harmonic degree of freedom
sets up the exchange-weight landscape for localization exchange.
```

---

## 4. System B [5,8]: Exchange-Weight Landscape of the Gray Metastable Interface

System B uses the white-cat, black-cat, and gray-cat metastable interface of self-reference [5] as its base, and sweeps `R` according to the specification in self-reference [8].

The state is classified by

```text
S = p_A - p_B.
```

The quality of the gray metastable phase is evaluated by

```text
gray_error = |S_mean| + |S_amp - S_amp_target| + S_drift + phase_penalty
gray_depth = -log10(gray_error)
```

where

```text
S_amp_target = 0.02.
```

In the initial broad and local sweeps, several candidate points such as `R=0.683`, `R=0.700`, and `R=0.697177879128` appeared. This indicates that the gray metastable interface has a multi-peak exchange-weight landscape along the `R` axis rather than a single peak.

![Initial exchange-weight landscape of System B](system_B_gray_cat_metastable_R_sweep_result_v1/system_B_odd_kernel_N_1_2_3_5_15_63_gray_depth_distribution_overview_v1.png)

**Figure 2.** Initial exchange-weight landscape of System B. The gray metastable interface has several local peaks rather than a single value of `R`.

However, the initial experiment used a coarse `R` range and might over-detect candidate points. Therefore, a minimized v5 implementation was used to scan a broad range with a uniform step.

---

## 5. Full-Range R Sweep with the Minimal v5 Implementation

### 5.1 Experimental Conditions

The v5 implementation removes unnecessary outputs for the investigation and reads only the gray-metastable depth for a specified `R`.

The full-range sweep conditions were:

| Item | Value |
|---|---:|
| `min_R` | `0.68660290255614798` |
| `max_R` | `0.70246536756843059` |
| `delta_R` | `0.0000001` |
| `n_R` | `158625` |
| `steps` | `1024` |
| `min_steps` | `256` |
| `early_stop_patience` | `20` |
| `phi_mode` | `zero` |

This range includes the neighborhood of `R_MZ` on the high-energy side and sufficiently covers `R_low` on the low-energy side.

### 5.2 Candidate Bands

In the full-range sweep, candidate points were grouped into continuous bands. Seven candidate bands were obtained.

| Rank | R_start | R_end | peak_R | depth | normalized depth | converted `1/alpha` |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | `0.697060302556148` | `0.697625902556148` | `0.697177902556148` | `9.08320492077` | `1.000000` | `137.036020287643` |
| 2 | `0.688146602556148` | `0.688539502556148` | `0.688363902556148` | `8.84756342227` | `0.974057` | `129.394062925467` |
| 3 | `0.692118902556148` | `0.692927402556148` | `0.692852802556148` | `5.84993869466` | `0.644039` | `133.203841605880` |
| 4 | `0.701163902556148` | `0.701587402556148` | `0.701489902556148` | `5.83931992499` | `0.642870` | `141.023604736761` |
| 5 | `0.686752802556148` | `0.687177802556148` | `0.687009702556148` | `5.81816843379` | `0.640541` | `128.276799088738` |
| 6 | `0.698514602556148` | `0.699071002556148` | `0.698893702556148` | `5.74438302514` | `0.632418` | `138.602220120334` |
| 7 | `0.692068702556148` | `0.692108002556148` | `0.692075802556148` | `4.61388675601` | `0.507958` | `132.532450378510` |

![Candidate band ranking](system_B_full_R_sweep_ranked_depth_bar_v1.svg)

**Figure 3.** Candidate bands detected in the full-range sweep. Candidate bands 1 and 2 are deep, while the remaining candidate bands are clearly shallower.

### 5.3 Full-Range Figure

In the full-range figure, candidate bands 1 and 2 stand out sharply.

![Full-range R sweep](system_B_full_R_sweep_full_range_depth_v1.svg)

**Figure 4.** Full-range R sweep. Candidate band 1 is near the point corresponding to low-energy `alpha(0)`, while candidate band 2 is near a high-energy effective-coupling region.

### 5.4 Candidate Band 1

The representative point of candidate band 1 is

```text
R_peak = 0.697177902556148.
```

The converted value from the standard-theory low-energy constant is

```text
R_low = 0.697177879231003.
```

The difference is

```text
Delta R ≈ 2.33e-8,
```

which is within the seven-digit effective precision of the v5 experiment.

![Candidate band 1 detail](system_B_full_R_sweep_band_detail_rank_01_v1.svg)

**Figure 5.** Detail of candidate band 1. A sharp depth appears near `R_low`.

### 5.5 Candidate Band 2

The representative point of candidate band 2 is

```text
R_peak = 0.688363902556148.
```

Converting this to `1/alpha` gives

```text
1/alpha = 129.394062925467.
```

The standard-theory reference value at the Z-boson mass scale is

```text
1/alpha(M_Z^2) = 128.946 ± 0.015
R_MZ = 0.687822933884774.
```

Candidate band 2 does not coincide with the standard value within the seven-digit precision of the v5 experiment. However, considering that the effective precision of the high-energy standard-theory coupling is roughly four to five digits, it has the following proximity:

```text
relative difference in 1/alpha ≈ 0.35%
relative difference in R       ≈ 0.08%
```

Therefore, candidate band 2 is not treated as identical to `alpha(M_Z^2)`, but is retained as a high-energy-side connection candidate.

![Candidate band 2 detail](system_B_full_R_sweep_band_detail_rank_02_v1.svg)

**Figure 6.** Detail of candidate band 2. It appears as a deep candidate band distinct from the low-energy side.

---

## 6. Sensitivity Near the Low-Energy and High-Energy Values

On the low-energy side, the candidate band appears stably near the standard value.

![Low-energy alpha neighborhood](system_B_low_alpha_R_sensitivity_depth_v1.svg)

**Figure 7.** R sensitivity near low-energy `alpha(0)`. The center of the candidate band remains near the standard value.

On the high-energy side, the result is better read as the appearance of a nearby candidate band rather than as a sharply fixed point at the standard value itself.

![High-energy alpha neighborhood](system_B_high_alpha_R_sensitivity_depth_v1.svg)

**Figure 8.** R sensitivity near high-energy `alpha(M_Z^2)`. The point is not as sharp as the low-energy fixed point, but a nearby candidate band appears.

In the intermediate region, no clear fixed point corresponding to a standard value was confirmed.

![Intermediate R neighborhood](system_B_mid_R_sensitivity_depth_v1.svg)

**Figure 9.** Intermediate region between the low-energy and high-energy sides. The candidate band is not strongly fixed at the center and is difficult to read as a standard correspondence point.

The comparison is therefore read as follows:

```text
Low-energy side:
  A strong candidate band exists at the R corresponding to alpha(0).

High-energy side:
  A deep candidate band exists near alpha(M_Z^2), but identity is not claimed.

Intermediate region:
  The correspondence candidate is weak.
```

---

## 7. Discussion

### 7.1 Why Are Harmonics Necessary?

In System A, the result was insensitive to the detailed amplitude, phase, and wavelength of the additional harmonic.

However, when the additional harmonic was exactly zero, the concentration near `R_137`, or `G_R≈sqrt(4πalpha(0))`, disappeared.

This suggests that the exchange-weight concentration may be related not to the detailed shape of the harmonic waveform, but to the existence of the additional degree of freedom itself.

In other words, in exchange scattering of closed phase systems, a two-component system with only base waves and a system with at least one additional harmonic degree of freedom may not be observing the same reaction space.

### 7.2 The Exchange-Weight Landscape Is Not Single-Peaked

The full-range sweep of System B did not produce a single candidate band.

Candidate band 1 strongly corresponds to low-energy `alpha(0)`.

Candidate band 2 corresponds to `1/alpha≈129.394`, and lies near the effective coupling at the Z-boson mass scale, `1/alpha(M_Z^2)=128.946±0.015`.

The remaining candidate bands are shallower and are not strongly connected to known standard-theory values at this stage.

This structure shows that `R≈0.7`, or `G_R≈0.3`, is not merely a single tuned value, but part of a multi-peak exchange-reaction landscape in the metastable interface of the closed system.

### 7.3 Candidate Correspondence with alpha, Not a Derivation of alpha

The experiments in this paper do not derive the fine-structure constant of the standard theory.

The essence of this paper is not the agreement with the fine-structure constant itself, but the observation that a closed exchange map built from a small set of axioms spontaneously selects a representative state-exchange weight. The correspondence with the fine-structure constant is treated as an observed result suggesting that the intrinsically appearing exchange weight may be related to the electromagnetic coupling in the real world.

Nevertheless, without inserting `alpha` from the outside, a uniform scan of `R` produced its deepest candidate band at the exchange weight corresponding to low-energy `alpha(0)`.

This suggests that the combination of exchange scattering, harmonic degrees of freedom, and a metastable interface in a closed phase system may not be unrelated to the electromagnetic coupling constant observed in the standard theory.

The conclusion of this paper is therefore limited to:

```text
In the exchange-weight landscape of exchange scattering in a closed phase system,
candidate bands that may correspond to the low-energy alpha
and the high-energy effective alpha of the standard theory appear.
```

---

## 8. Claims Not Made in This Paper

This paper does not claim:

```text
that the fine-structure constant alpha has been derived;
that the standard theory has been replaced;
that candidate band 2 is exactly identical to alpha(M_Z^2);
that every fermion-scattering coefficient concentrates at the same exchange weight;
that the geometric origin of the harmonic degree of freedom has been fully explained.
```

What this paper shows is the numerical fact that alpha-correspondence candidates appear in exchange-scattering experiments of closed phase systems.

---

## 9. Conclusion

This paper used the exchange-scattering coefficient `R` as the sweep coordinate and investigated, through two systems, the selection of the representative state-exchange weight `G_R=1-R`: a localization-exchange model and a gray-cat metastable-interface model.

In System A, a single base harmonic did not produce a sharp exchange-weight concentration. Once at least one additional harmonic degree of freedom was present, the evaluation point of localization exchange concentrated near `R=0.697177879...`, or equivalently near `G_R=1-R≈0.302822121`. This concentration was robust under odd harmonics, even harmonics, mixed parity packets, reduced amplitude, phase shifts, and wavelength shifts.

In System B, the minimal v5 implementation of the gray-cat metastable interface swept

```text
R = 0.686602902... to R = 0.702465367...
```

with

```text
Delta R = 1.0e-7.
```

Seven candidate bands appeared. The deepest candidate band was

```text
R = 0.697177902556148.
```

This agrees, within the seven-digit effective precision of the v5 experiment, with

```text
R_low = 0.697177879231003
```

obtained from the CODATA 2022 low-energy fine-structure constant by

```text
R = 1 - sqrt(4 pi alpha).
```

The second candidate band was

```text
R = 0.688363902556148,
```

corresponding to `1/alpha=129.394062925467`. This is not exactly identical to the standard reference value of `alpha(M_Z^2)`, but is retained as a high-energy effective-coupling candidate.

Thus, in a closed phase-system exchange-scattering map, when a metastable interface with harmonic degrees of freedom is read, candidate bands of the representative state-exchange weight corresponding to the fine-structure constant appear spontaneously.

---

## References

### Self-References

1. Noriaki Kihara, [Basic Axiom System v4](../20260710/基本公理系%20v4.md), 2026.
2. Noriaki Kihara, [Interference Construction of a Complete Reflection Map by a Fermion-Like Opposite-Phase Kernel v2](../20260710/フェルミオン的逆相核による完全反射写像の干渉構成%20v2.md), 2026.
3. Noriaki Kihara, [Experiment Specification for Low-Localization and Harmonic-Transfer Readout in Exchange-Interference Scattering-Matrix Fermion-Like Collisions v1](../20260713/交換干渉散乱行列フェルミオン的衝突における低局在性・倍音移乗読出し実験仕様書%20v1.md), 2026.
4. Noriaki Kihara, [Preliminary Summary of Low-Localization and Harmonic Transfer in Exchange-Interference Scattering-Matrix Fermion-Like Collisions v1](../20260713/交換干渉散乱行列フェルミオン的衝突における低局在性・倍音移乗予備実験総括%20v1.md), 2026.
5. Noriaki Kihara, [Preliminary Summary of C Weak Readout and D Strong-Observation Selection at the White-Cat, Black-Cat, and Gray-Cat Metastable Interface v1](../20260714/白猫黒猫灰色猫準安定界面におけるC弱読出しとD強観測選択予備実験総括%20v1.md), 2026.
6. Noriaki Kihara, [Overall Plan and Evaluation Method for the Exchange-Scattering Coefficient R Concentration Experiment Group v1](交換散乱係数R集中実験群_全体計画と評価方法_v1.md), 2026.
7. Noriaki Kihara, [System A: Experiment Specification for Uniform R-Neighborhood Sweep in Localization Exchange v1](系統A_局在性交換R近傍斉一スイープ実験仕様書_v1.md), 2026.
8. Noriaki Kihara, [System B: Experiment Specification for R-Neighborhood Sweep of the Gray-Cat Metastable Interface v1](系統B_灰色猫準安定界面R近傍スイープ実験仕様書_v1.md), 2026.
9. Noriaki Kihara, [alpha in the Standard Theory](標準理論でのα.md), 2026.

### External References

10. Peter J. Mohr, David B. Newell, Barry N. Taylor, Eite Tiesinga, "CODATA Recommended Values of the Fundamental Physical Constants: 2022", arXiv:2409.03787, 2024. https://arxiv.org/abs/2409.03787
11. National Institute of Standards and Technology, "Fundamental Physical Constants", https://physics.nist.gov/constants
12. Alexander Keshavarzi, Daisuke Nomura, Thomas Teubner, "The muon g-2 and alpha(M_Z^2): a new data-based analysis", arXiv:1802.02995, 2018. https://arxiv.org/abs/1802.02995

---

## Appendix A. Main Outputs Used for Reproduction

| Type | File |
|---|---|
| System A N-series | `system_A_localization_exchange_R_sweep_result_v1/odd_kernel_N_1_2_3_5_15_63_report_v1.md` |
| System A even packet | `system_A_localization_exchange_R_sweep_result_v1/even_packet_A1_B1_2_4_6_report_v1.md` |
| System A mixed packet 1 | `system_A_localization_exchange_R_sweep_result_v1/alternating_packet_1_A1_B1_2_3_4_5_report_v1.md` |
| System A mixed packet 2 | `system_A_localization_exchange_R_sweep_result_v1/alternating_packet_2_A1_B1_3_4_5_6_report_v1.md` |
| System A zero-amplitude control | `system_A_localization_exchange_R_sweep_result_v1/amp000_B12/system_A_amp000_B12_custom_packet_A1_B1-2-w-1-0_Rdefault_C256_report_v1.md` |
| System B full-range candidate bands | `System B 全域Rスイープ候補帯一覧.md` |
| System B low-energy neighborhood | `System B 低エネルギー標準α近傍R感度一覧.md` |
| System B high-energy neighborhood | `System B 高エネルギー標準α近傍R感度一覧.md` |
| System B intermediate region | `System B 中間R近傍R感度一覧.md` |
