# Experimental Specification for the White-Cat, Black-Cat, and Gray-Cat Metastable Interface in AB-C-D Stages v1

**Subtitle:** Closed-system numerical experiment on the AB two-body metastable interface, weak C readout, and strong D observation for white-cat, black-cat, and gray-cat branching  
**Date:** 2026-07-14  
**Author:** Noriaki Kihara  
**Series:** Wave Information Readout / staged experiment for macroscopic cat-like metastable states  
**Target implementation:** Python numerical experiment  
**Version DOI:** 10.5281/zenodo.21353209  
**Concept DOI:** 10.5281/zenodo.21353208

---

## 0. Summary

This experiment does not assume in advance that an equal-allocation state of the white-cat state `A` and the black-cat state `B` is merely an unselected superposition.

First, the AB two-body system alone is used to distinguish three phases.

```text
gray eigen phase:
  A/B remains fixed near 0.5/0.5.

gray metastable phase:
  A/B oscillates weakly near 0.5/0.5 but does not naturally converge to one side.

natural selection phase:
  A or B is selected even without C or D.
```

Then C is added to test whether the metastable state can be read without selection. Finally D is added to test whether observation selects either A or B.

The first checkpoint is not the action of C or D.

```text
The first checkpoint is to find the interface between the gray eigen phase,
gray metastable phase, and natural selection phase in AB alone.
```

---

## 1. Purpose

The experiment checks the following stages in order.

1. Classify the gray eigen phase, gray metastable phase, and natural selection phase in AB alone.
2. Add C and test whether a gray eigen or metastable state can be read without selection.
3. Add D and test the conditions under which observation selects A or B.
4. Confirm whether a gray eigen phase remains gray even under D.
5. In large-amplitude regions where A/B is already separated, confirm that D agrees with the C readout.

D is not assumed to always cause collapse.

If a gray cat has already become an eigen state, it does not need to fall into white or black under D.

If A/B is already clearly separated, the D result is treated as a readout of an already existing A/B branching, not as a selection created by D.

---

## 2. State Variables

At minimum, the following variables are recorded.

| variable | meaning |
|---|---|
| `a` | complex amplitude of white-cat A |
| `b` | complex amplitude of black-cat B |
| `p_A = |a|^2` | white-cat allocation |
| `p_B = |b|^2` | black-cat allocation |
| `Q = p_A + p_B` | total A/B amount |
| `S = p_A - p_B` | A/B selection order variable |
| `S_amp` | oscillation amplitude of `S` |
| `S_drift` | long-time drift of `S` |
| `C_A, C_B` | A/B readout by C |
| `D_A, D_B` | A/B readout by D |
| `L_A, L_B` | localization of A/B |
| `N_eff_A, N_eff_B` | effective harmonic order of A/B |
| `E_total` | total conserved quantity or conservation registry quantity |

The central checks are:

```text
whether Q is conserved;
whether S is fixed, weakly oscillating, or naturally growing;
whether C and D change the phase classification of S.
```

---

## 3. AB Two-Body Metastable Interface Experiment

## 3.1 Purpose

The AB system alone is used to test what kind of state the gray cat is.

```text
C = off
D = off
```

No observation-induced selection is included at this stage.

## 3.2 Initial State

The reference initial values are

```math
a_0 = \frac{1}{\sqrt{2}},
```

```math
b_0 = \frac{e^{i\phi_0}}{\sqrt{2}}.
```

Then

```math
Q_0 = |a_0|^2 + |b_0|^2 = 1,
```

```math
S_0 = |a_0|^2 - |b_0|^2 = 0.
```

When a small initial bias is inserted,

```math
p_{A,0} = \frac12 + s_0,
```

```math
p_{B,0} = \frac12 - s_0.
```

## 3.3 AB Exchange Map

The minimal AB exchange map is

```math
\begin{pmatrix}
a_{k+1}\\
b_{k+1}
\end{pmatrix}
=
U_\epsilon
\begin{pmatrix}
a_k\\
b_k
\end{pmatrix},
```

```math
U_\epsilon
=
\begin{pmatrix}
\cos\epsilon & i\sin\epsilon\\
i\sin\epsilon & \cos\epsilon
\end{pmatrix}.
```

This map alone gives a simple unitary exchange oscillation.

If needed, a weak restoring term or weak nonlinear term is added to test closed-system stability. However, a term that causes one-sided natural convergence without C or D is not used for the cat-selection experiment.

## 3.4 Sweep Parameters

At minimum, the following parameters are swept.

| parameter | meaning |
|---|---|
| `epsilon` | strength of AB exchange oscillation |
| `phi_0` | initial relative phase |
| `s_0` | small initial A/B bias |
| `noise_amp` | numerical fluctuation or small disturbance |
| `K` | number of evolution steps |

When the gray cat is treated as a harmonic metastable state, the amplitude of `S` is kept small.

The target condition is

```text
S_amp < S_gray_limit.
```

The initial guide values are:

```text
S_gray_limit = 0.05
S_amp target = about 0.01 to 0.03
```

## 3.5 Phase Classification

The AB two-body results are classified into three phases.

### Gray Eigen Phase

```text
S_mean ≈ 0
S_amp -> 0
S_drift ≈ 0
Q ≈ 1
```

In this case, the equal A/B allocation is treated as a new gray-cat eigen state.

### Gray Metastable Phase

```text
S_mean ≈ 0
0 < S_amp < S_gray_limit
S_drift ≈ 0
Q ≈ 1
```

This phase is the main target for the D selection experiment.

### Natural Selection Phase

```text
abs(S) grows without C or D
S -> +1 or S -> -1
```

In this phase, observation-induced selection cannot be claimed.

---

## 4. C Readout Experiment

## 4.1 Purpose

C is added to a gray eigen or gray metastable state found in AB.

```text
C = on
D = off
```

C reads the A/B allocation. At this stage, C must not push the system into one side.

## 4.2 C Readout Quantities

C reads:

```math
C_A \approx p_A,
```

```math
C_B \approx p_B.
```

The error is recorded as

```math
E_C = |C_A-p_A| + |C_B-p_B|.
```

## 4.3 C Acceptance Conditions

For the gray metastable phase:

```text
E_C < tol_C
abs(S_afterC - S_beforeC) < epsilon_C
no selection
```

For the gray eigen phase:

```text
C_A ≈ 0.5
C_B ≈ 0.5
S remains near 0
```

For the large-amplitude region, C reads the current A/B bias.

```text
C_A > C_B if S > 0
C_B > C_A if S < 0
```

---

## 5. D Observation Experiment

## 5.1 Purpose

D is then added.

```text
C = optional
D = on
```

D is a strong observation map used to test whether the system falls into either A or B.

The result is classified by state phase.

## 5.2 D in the Gray Eigen Phase

For the gray eigen phase, the test asks whether the gray state remains gray under D.

The expected judgment is:

```text
D_A ≈ 0.5
D_B ≈ 0.5
S remains near 0
```

In this case, not falling into white or black is the result.

## 5.3 D in the Gray Metastable Phase

For the gray metastable phase, the test asks whether D selects one side.

The judgment is:

```text
D_A ≈ 1, D_B ≈ 0
```

or

```text
D_A ≈ 0, D_B ≈ 1
```

If the same fall occurs without D, it is not judged as D-induced selection.

## 5.4 D in the Large-Amplitude Separated Region

When A/B is already clearly separated, D only needs to read the current A/B bias.

The D readout should agree with the C readout:

```text
sign(D_A - D_B) = sign(C_A - C_B)
```

If this agreement is obtained, D is judged to have read an already existing A/B bias, rather than creating the selection.

---

## 6. Four-Stage Main Experiment

## Stage 1: AB Metastable Interface Search

```text
C = off
D = off
```

Output:

```text
gray eigen phase
gray metastable phase
natural selection phase
```

## Stage 2: C Readout Confirmation

```text
C = on
D = off
```

Output:

```text
C_read_error
C_induced_bias
C_selection_triggered
```

## Stage 3: D Observation Confirmation

```text
C = off or recorded
D = on
```

Output:

```text
D_result
D_induced_selection
D_vs_C_agreement
```

## Stage 4: C-D Comparison

The same AB initial condition is compared under C readout and D observation.

Checks:

```text
gray eigen phase:
  C and D both read 0.5/0.5.

gray metastable phase:
  C reads but does not select.
  D is tested for one-sided selection.

large-amplitude separated region:
  C and D read the same A/B bias.
  This is not treated as observation-created selection.
```

---

## 7. Acceptance Conditions

The minimum acceptance conditions are:

1. AB alone can classify the gray eigen phase, gray metastable phase, and natural selection phase.
2. In the gray metastable phase, C can read the A/B allocation without selecting one side.
3. In the gray eigen phase, D does not make the system fall into white or black.
4. In the gray metastable phase, D-induced selection can be judged by comparison with a no-D control.
5. In the large-amplitude separated region, D agrees with the C readout.
6. `Q` or the total conservation registry remains conserved within tolerance.

---

## 8. Failure Conditions

The experiment does not proceed to later C/D stages if:

1. All AB-only conditions naturally fall into selection.
2. `Q` is not conserved.
3. C always causes selection.
4. C cannot read the A/B allocation.
5. D/no-D differences cannot be classified.

---

## 9. Result Fields

## 9.1 Stage 1: AB Metastable Interface Search

Output:

```text
gray_cat_ab_metastable_interface_preliminary_result_v1/
```

Result:

```text
total_cases = 5600
gray_eigen = 1248
gray_metastable = 733
large_oscillation = 470
natural_selection = 1070
unstable_or_drifting = 2079
```

## 9.2 Stage 2: C Readout Window

Output:

```text
gray_cat_c_readout_window_preliminary_result_v1/
```

Result:

```text
total_cases = 1080
C_window_count = 247
C_informative_window_count = 144
C_nonzero_backaction_window_count = 114
```

## 9.3 Stage 3: D Observation Response

Output:

```text
gray_cat_d_observation_response_preliminary_result_v1/
```

Result:

```text
total_cases = 7056
D_induced_selection_count = 2016
gray_kept_eigen = 1062
white_selected = 1244
black_selected = 772
```

## 9.4 Stage 4: D Selection Boundary

Output:

```text
gray_cat_d_selection_boundary_preliminary_result_v1/
```

Result:

```text
boundary_points = 438
selection_possible_boundary_points = 438
no_selection_boundary_points = 0
min_D_gain_overall = 0.0225
max_min_D_gain_overall = 0.065
```

