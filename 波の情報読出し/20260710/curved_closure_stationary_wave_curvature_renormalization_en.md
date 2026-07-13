# Curvature Renormalization and Perfect-Reflection Stability by Curved Closed Stationary Waves v2

**Subtitle:** Numerical constructive experiment on the detection, re-selection, and exchange-interference reflection recovery of curvature-induced relative phase leakage under the all-positive zero closure `Σx_n^2=0`
**Date:** 2026-07-11
**Author:** Noriaki Kihara
**Position:** Additional paper in the Wave Information Readout series
**Version DOI:** 10.5281/zenodo.21332874
**Concept DOI:** 10.5281/zenodo.21304039

V2 recalculates the same conditions after changing the formula for the fermion-like reflection map.

---

## Abstract

This paper numerically constructs how an odd-harmonic complex wave satisfying the all-positive zero closure

```math
Q(x)=\sum_n x_n^2=0
```

is stabilized inside a curved local cell.

This closure condition is not introduced here as an auxiliary assumption for the curvature test. It is Axiom 1 of the existing basic axiom system of this paper series, together with namelessness and non-trivial existence. In that system, it was already placed as the existence condition for non-trivial complex phase waves. The present paper applies Axiom 1 to curved local cells and exchange-interference reflection.

In the preceding paper, a local exchange-interference map was constructed from a fermionic inverse-phase core, and a direction-reversal output corresponding to complete elastic reflection was obtained without using an external `q -> -q` instruction. The present paper asks whether that exchange-interference map breaks under curvature-induced relative phase leakage, and whether it is recovered by re-selection into a closed stationary wave.

The curvature effect is not treated as a quantitative prediction of real spacetime curvature. Instead, curvature-induced relative phase leakage inside a local cell is modeled as

```math
\delta_{K,m}.
```

If this leakage is inserted into a closure pair

```math
x_m,\qquad ix_m,
```

then

```math
x_m^2+\left(e^{i\delta_{K,m}}ix_m\right)^2
=
x_m^2\left(1-e^{i2\delta_{K,m}}\right),
```

and a closure residual appears in general. If an internal phase re-selection `\beta_{K,m}` satisfies

```math
\delta_{K,m}+\beta_{K,m}=0,
```

the closure condition is recovered.

In the minimal experiment, the transient state with curvature-induced relative phase leakage produced closure-pair RMS `1.2319416790092972e-02` and transmission leakage `1.1503183254481797e-01`. After internal phase re-selection, the stationary state recovered closure-pair RMS `9.4283259783636047e-19` and transmission leakage `0.0`.

In the broad verification, eight curvature-relative-phase models were tested against seven correction freedoms: `none`, `constant`, `linear`, `affine`, `quadratic`, `cubic`, and `full`. At maximum curvature relative phase `1.2`, the uncorrected case left maximum transmission leakage `1.6202719613622976e-01`, whereas the `full` correction recovered closure-pair RMS `7.8949412793793227e-19` and transmission leakage `0.0`.

In the integrated one-sided scattering verification, the residual curvature phase obtained from the broad sweep was returned to the local exchange-interference map. The uncorrected case produced maximum dynamic transmission leakage `1.6202719613622971e-01`, whereas the `full` correction reduced it to `1.6608667985580024e-19`. The dynamic scattering matched the two-channel expectation with maximum error `5.551115123125783e-17`, and the maximum norm error was `4.440892098500626e-16`.

Within the numerical constructive scope of this paper, curvature effects are therefore not absent. They appear in transient states as closure residuals and transmission leakage. However, once the system is re-selected into a closed stationary wave satisfying `Σx_n^2=0`, the curvature-induced relative phase leakage is absorbed into the internal phase configuration, and the complete-reflection readout is recovered.

**Keywords:** all-positive zero closure, curved closed stationary wave, curvature renormalization, odd harmonics, exchange interference, complete reflection, local constant-curvature cell, internal phase re-selection

---

## 1. Introduction

### 1.1 Background

This series adopts namelessness, all-positive zero closure, and non-trivial existence as its basic axioms. In particular, all-positive zero closure is Axiom 1 of the basic axiom system v2. It is not an auxiliary assumption added for the present paper. The central closure condition is

```math
\sum_n x_n^2=0.
```

This is not the conjugate norm

```math
\sum_n |x_n|^2.
```

It is a closure condition in which each component is squared as it is and summed with all-positive signs.

In the minimal two-component case,

```math
A^2+(iA)^2=0,
```

so sign reversal is generated internally without using an external negative coefficient. Complex phase is therefore not an ornament; it appears as the phase algebra required for non-trivial closure.

The preceding paper, "Interference Construction of a Perfect Reflection Map from a Fermionic Inverse-Phase Core," constructed a local exchange-interference map from an internal inverse-phase core and two exchange paths. It generated a direction-reversal readout corresponding to complete reflection without using an external `q -> -q` instruction.

The present paper investigates whether this exchange-interference reflection remains stable when placed in a curved local cell and exposed to curvature-induced relative phase leakage.

### 1.2 Question

The question of this paper is:

> When curvature-induced relative phase leakage occurs in a curved local cell, can the all-positive zero closure `Σx_n^2=0` detect the leakage and recover the complete-reflection readout through internal phase re-selection into a closed stationary wave?

The important point is that curvature effects are not ignored at the outset. Instead, this paper first confirms that curvature-induced relative phase leakage produces closure residuals and transmission leakage, and then tests whether those residuals disappear after re-selection into a closed stationary wave.

---

## 2. What This Paper Does Not Claim

This paper does not claim the following.

| Not claimed | Reason |
|---|---|
| Quantitative prediction of real spacetime curvature effects | No standard spacetime metric, field equation, or experimental unit system is used |
| Curvature is always unobservable | Closure residuals and transmission leakage appear in transient states |
| A unique physical formula for the curvature relative phase `δ_K` | It is used here as a local relative-phase model for verification |
| Derivation of standard quantum theory or general relativity | This is a numerical constructive experiment on an internal axiom system |
| Derivation of a continuous interaction Hamiltonian | This paper treats local maps and existence conditions of closed stationary waves |

The claim of this paper is a minimal numerical construction: a closed stationary wave satisfying `Σx_n^2=0` detects curvature-induced relative phase leakage and recovers the exchange-interference reflection readout through internal phase re-selection.

---

## 3. Basic Axiom and the Closure Null Cone

### 3.1 All-Positive Zero Closure

This section restates Axiom 1 of the basic axiom system v2. The present paper uses this axiom as its starting point, but it is not introduced here ad hoc to make curvature renormalization work.

Axiom 1 sets the closure condition as

```math
Q(x)=\sum_{n=1}^N x_n^2=0.
```

The set of closed states is denoted by

```math
\mathcal N
=
\{x\mid Q(x)=0\}.
```

`\mathcal N` is not the set of states with zero conjugate norm. It is the complex closure null cone defined by the all-positive sum of squares.

### 3.2 Closure Pair

The minimal closure pair is

```math
x_m,
\qquad
ix_m.
```

Then

```math
x_m^2+(ix_m)^2=0.
```

For an odd-harmonic complex wave, `m` is the harmonic label and

```math
h_m=2m+1
```

is used.

---

## 4. Curved Local Cell and Relative Phase Leakage

### 4.1 Two Types of Curvature Action

Curvature action is divided into two types.

The first is common-factor action:

```math
x_m\mapsto g_{K,m}x_m,
\qquad
ix_m\mapsto g_{K,m}ix_m.
```

In this case,

```math
(g_{K,m}x_m)^2+(g_{K,m}ix_m)^2=0,
```

and the closure pair is preserved.

The second is relative-phase action:

```math
x_m\mapsto x_m,
\qquad
ix_m\mapsto e^{i\delta_{K,m}}ix_m.
```

In this case,

```math
x_m^2+(e^{i\delta_{K,m}}ix_m)^2
=
x_m^2\left(1-e^{i2\delta_{K,m}}\right),
```

and a closure residual appears in general.

### 4.2 Closed Stationary Wave

A wave that remains stable in a curved local cell is called a closed stationary wave `x_K`.

It satisfies

```math
Q(x_K)=0
```

and

```math
\mathcal U_K x_K=e^{i\alpha_K}x_K.
```

Here `\mathcal U_K` is the local evolution map after passing through the curved local cell, and `e^{i\alpha_K}` is a common phase removed by external readout.

---

## 5. Internal Re-Selection of Curvature Relative Phase

For relative phase leakage `\delta_{K,m}`, introduce an internal phase re-selection `\beta_{K,m}`.

If

```math
\delta_{K,m}+\beta_{K,m}=0
```

holds, then

```math
e^{i(\delta_{K,m}+\beta_{K,m})}=1,
```

and the closure pair is recovered.

Thus curvature-induced relative phase leakage is not erased. It is re-selected into an internal phase configuration that can exist as a closed stationary wave.

---

## 6. Connection to Exchange-Interference Reflection

In the preceding exchange-interference reflection paper, the effective phase was

```math
\Delta_{\mathrm{eff}}=\Delta_F,
```

and the pure inverse-phase core `\Delta_F=\pi` produced complete reflection.

In this paper, including curvature-induced relative phase leakage,

```math
\Delta_{\mathrm{eff},m}
=
\Delta_F+\delta_{K,m}+\beta_{K,m}.
```

The complete-reflection condition is

```math
\Delta_{\mathrm{eff},m}=\pi.
```

Without curvature correction,

```math
\Delta_{\mathrm{eff},m}=\pi+\delta_{K,m},
```

so transmission leakage appears.

If

```math
\beta_{K,m}=-\delta_{K,m},
```

then

```math
\Delta_{\mathrm{eff},m}=\pi
```

is recovered, and the complete-reflection readout returns.

---

## 7. Numerical Experiments

### 7.1 Experiment 1: Minimal Closed Stationary Wave Verification

The executed script is:

```text
run_curved_closure_stationary_wave_v2.py
```

The output directory is:

```text
curved_closure_stationary_wave_result_v2/
```

The curvature relative phase leakage is set as

```math
\delta_{K,m}=\kappa h_m.
```

For maximum `\kappa=0.012` and maximum curvature phase `1.1879999999999999`, the following values were obtained.

| State | Closure-pair RMS | Transmission leakage |
|---|---:|---:|
| flat | `0.0000000000000000e+00` | `0.0000000000000000e+00` |
| conformal | `9.4283259783636047e-19` | `0.0000000000000000e+00` |
| transient | `1.2319416790092972e-02` | `1.1503183254481797e-01` |
| stationary | `9.4283259783636047e-19` | `0.0000000000000000e+00` |

In the relaxation sequence, the final closure-pair RMS was `2.0285753500943141e-09`, and the final transmission leakage was `2.4915322496409796e-15`.

All verdict flags were `true`.

### 7.2 Experiment 2: Broad Sweep of Curvature Phase Models

The executed script is:

```text
run_curved_closure_stationary_wave_broad_sweep_v2.py
```

The output directory is:

```text
curved_closure_stationary_wave_broad_sweep_result_v2/
```

The eight curvature relative phase models were:

```text
linear
quadratic_area
cubic_high
quartic_edge
alternating_linear
sinusoidal_loop
mixed_smooth
rippled_random_like
```

The seven internal correction freedoms were:

```text
none
constant
linear
affine
quadratic
cubic
full
```

The correction-wise aggregate at maximum curvature relative phase `1.2` was:

| Correction | Max closure-pair RMS | Max transmission leakage | Max residual phase |
|---|---:|---:|---:|
| none | `1.3992789770439718e-02` | `1.6202719613622976e-01` | `1.2000000000000000e+00` |
| constant | `1.2323357129304359e-02` | `1.1632365832184562e-01` | `1.1885985027360557e+00` |
| linear | `1.2317721991993871e-02` | `1.1624067781252551e-01` | `1.2122605013795469e+00` |
| affine | `1.2315221532542749e-02` | `1.1619972192010980e-01` | `1.2235371932078578e+00` |
| quadratic | `1.2300337923101595e-02` | `1.1598359742226144e-01` | `1.2740926312051919e+00` |
| cubic | `1.2278927640523116e-02` | `1.1567739316276086e-01` | `1.3311608455010591e+00` |
| full | `7.8949412793793227e-19` | `0.0000000000000000e+00` | `0.0000000000000000e+00` |

Common-factor curvature preserved closure throughout the sweep. Without correction, non-trivial curvature leakage was detected. With `full` correction, closure and complete reflection were recovered for all phase models. With limited corrections, model-dependent residuals remained.

### 7.3 Experiment 3: Integrated Verification with One-Sided Scattering

The executed script is:

```text
run_curved_closure_scattering_integration_v2.py
```

The output directory is:

```text
curved_closure_scattering_integration_result_v2/
```

The residual curvature relative phases from Experiment 2 were returned to the local exchange-interference map for one-sided incident scattering.

The correction-wise aggregate was:

| Correction | Max dynamic transmission leakage | Max closure-pair RMS | Max expectation error |
|---|---:|---:|---:|
| none | `1.6202719613622971e-01` | `1.3992789770439718e-02` | `5.5511151231257827e-17` |
| constant | `1.1632365832184560e-01` | `1.2323357129304359e-02` | `1.3877787807814457e-17` |
| linear | `1.1624067781252549e-01` | `1.2317721991993871e-02` | `1.3877787807814457e-17` |
| affine | `1.1619972192010980e-01` | `1.2315221532542749e-02` | `1.0408340855860843e-17` |
| quadratic | `1.1598359742226147e-01` | `1.2300337923101595e-02` | `2.7755575615628914e-17` |
| cubic | `1.1567739316276089e-01` | `1.2278927640523116e-02` | `2.7755575615628914e-17` |
| full | `1.6608667985580024e-19` | `7.8949412793793227e-19` | `1.6608667985580024e-19` |

The dynamic scattering matched the two-channel expectation with maximum error `5.551115123125783e-17`, and the maximum norm error was `4.440892098500626e-16`.

---

## 8. Classification of Results

The results are classified as follows.

| Target | Classification | Verdict |
|---|---|---|
| `Σx_n^2=0` requires non-trivial complex closure | Derived consequence | retained |
| Curvature relative phase leakage produces closure residuals | Numerically constructed consequence | retained |
| Common-factor curvature action does not break closure | Numerically constructed consequence | retained |
| Internal phase re-selection recovers closure | Numerically constructed consequence | retained |
| Internal phase re-selection recovers complete reflection | Numerically constructed consequence | retained |
| Real-spacetime local flatness is explained by this structure | Connection to existing theory | not claimed |

---

## 9. Discussion

### 9.1 Curvature Has Not Disappeared

The result of this paper is not that curvature has no effect.

In the transient state with curvature-induced relative phase leakage, closure residuals and transmission leakage clearly appear.

Curvature effects are therefore detectable in this model.

### 9.2 Closed Stationary Waves Absorb Curvature Internally

The important case is the state after curvature-induced relative phase leakage has occurred and the system is re-selected as a closed stationary wave.

If the internal phase correction `\beta_{K,m}` satisfies

```math
\delta_{K,m}+\beta_{K,m}=0,
```

then closure residuals and transmission leakage disappear.

This does not mean that curvature is absent. It means that the curvature relative phase is renormalized into the internal phase configuration of the closed stationary wave that can stably exist.

### 9.3 Relation to Local Flatness

In the standard local-flatness approximation, curvature effects are ignored in a sufficiently small region.

The construction in this paper is different.

Curvature effects are first inserted, closure is observed to break, and then the residual disappears from readout only after the state is re-selected as a closed stationary wave.

Thus, within the internal axiom system of this paper, a readout that looks locally flat can be interpreted not as the result of ignoring curvature, but as the result that only curvature-inclusive stationary interference modes that close are stably read out.

This is not an identification with standard theory. It is a working hypothesis for constructing a correspondence map.

### 9.4 Connection to Complete Elastic Reflection

The complete-reflection map constructed in the preceding paper was an exchange-interference map with an internal inverse-phase core:

```math
\Delta_F=\pi.
```

In the present paper, adding curvature-induced relative phase leakage gives

```math
\Delta_{\mathrm{eff},m}
=
\pi+\delta_{K,m}+\beta_{K,m}.
```

Without correction, `\delta_{K,m}` remains and transmission leakage appears. Under `full` correction, `\beta_{K,m}=-\delta_{K,m}`, so

```math
\Delta_{\mathrm{eff},m}=\pi
```

is recovered.

At that point, the one-sided scattering experiment also recovers `R=1,T=0`.

---

## 10. Conclusion

This paper numerically constructed how an odd-harmonic complex wave satisfying the all-positive zero closure `Σx_n^2=0` is stabilized inside a curved local cell.

In the minimal experiment, the transient state with curvature-induced relative phase leakage produced closure-pair RMS `1.2319416790092972e-02` and transmission leakage `1.1503183254481797e-01`. After internal phase re-selection, the stationary state recovered closure-pair RMS `9.4283259783636047e-19` and transmission leakage `0.0`.

In the broad verification, eight curvature relative phase models and seven correction freedoms were compared. Without correction, maximum transmission leakage `1.6202719613622976e-01` remained. Under `full` correction, closure-pair RMS `7.8949412793793227e-19` and transmission leakage `0.0` were recovered.

In the integrated one-sided local exchange-interference scattering verification, the uncorrected case produced maximum dynamic transmission leakage `1.6202719613622971e-01`, while the `full` correction reduced maximum dynamic transmission leakage to `1.6608667985580024e-19`.

Therefore, within the numerical constructive scope of this paper, curvature effects are not ignored. Curvature-induced relative phase leakage appears in transient states as closure residuals and transmission leakage. However, once the system is re-selected into a closed stationary wave satisfying `Σx_n^2=0`, the curvature-induced relative phase leakage is absorbed into the internal phase configuration, and the complete-reflection readout by exchange interference is recovered.

Thus `Σx_n^2=0` functions not merely as a conservation condition, but as an existence condition for non-trivial complex waves, a detection condition for curvature-induced relative phase leakage, a re-selection condition into closed stationary waves, and a stability condition for the complete-reflection readout.

---

# Appendix A. Executed Programs and Outputs

## A.1 Minimal Verification

```text
python3 run_curved_closure_stationary_wave_v2.py
```

Output:

```text
curved_closure_stationary_wave_result_v2/
```

Main files:

| Type | File |
|---|---|
| Report | [curved_closure_stationary_wave_report_v2.md](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_report_v2.md) |
| JSON | [curved_closure_stationary_wave_result_v2.json](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_result_v2.json) |
| sweep CSV | [curved_closure_stationary_wave_sweep_v2.csv](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_sweep_v2.csv) |
| relaxation CSV | [curved_closure_stationary_wave_relaxation_v2.csv](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_relaxation_v2.csv) |
| sweep figure | [curved_closure_stationary_wave_sweep_v2.png](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_sweep_v2.png) |
| relaxation figure | [curved_closure_stationary_wave_relaxation_v2.png](curved_closure_stationary_wave_result_v2/curved_closure_stationary_wave_relaxation_v2.png) |

## A.2 Broad Verification

```text
python3 run_curved_closure_stationary_wave_broad_sweep_v2.py
```

Output:

```text
curved_closure_stationary_wave_broad_sweep_result_v2/
```

Main files:

| Type | File |
|---|---|
| Report | [curved_closure_stationary_wave_broad_report_v2.md](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_report_v2.md) |
| JSON | [curved_closure_stationary_wave_broad_sweep_result_v2.json](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_sweep_result_v2.json) |
| sweep CSV | [curved_closure_stationary_wave_broad_sweep_v2.csv](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_sweep_v2.csv) |
| control CSV | [curved_closure_stationary_wave_broad_conformal_control_v2.csv](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_conformal_control_v2.csv) |
| aggregate figure | [curved_closure_stationary_wave_broad_aggregate_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_aggregate_v2.png) |
| closure heatmap | [curved_closure_stationary_wave_broad_closure_heatmap_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_closure_heatmap_v2.png) |
| leakage heatmap | [curved_closure_stationary_wave_broad_leakage_heatmap_v2.png](curved_closure_stationary_wave_broad_sweep_result_v2/curved_closure_stationary_wave_broad_leakage_heatmap_v2.png) |

## A.3 One-Sided Scattering Integration

```text
python3 run_curved_closure_scattering_integration_v2.py
```

Output:

```text
curved_closure_scattering_integration_result_v2/
```

Main files:

| Type | File |
|---|---|
| Report | [curved_closure_scattering_integration_report_v2.md](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_report_v2.md) |
| JSON | [curved_closure_scattering_integration_result_v2.json](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_result_v2.json) |
| CSV | [curved_closure_scattering_integration_v2.csv](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_v2.csv) |
| Figure | [curved_closure_scattering_integration_v2.png](curved_closure_scattering_integration_result_v2/curved_closure_scattering_integration_v2.png) |

---

# References

## Self-Citations

1. Noriaki Kihara, "Basic Axiom System v2 for the Nameless Equal-Amplitude Composite Wave Model," 2026-07-10.
2. Noriaki Kihara, "Paper 0: Distortion of the Geodesic Unit Cell in Positively-Curved Constant-Curvature Space — Exact Evaluation of Edge, Angle, Area, and Volume," Version DOI: `10.5281/zenodo.21303433`, Concept DOI: `10.5281/zenodo.20680269`, 2026.
3. Noriaki Kihara, "Interference Construction of a Perfect Reflection Map from a Fermionic Inverse-Phase Core," Version DOI: `10.5281/zenodo.21332867`, Concept DOI: `10.5281/zenodo.21295479`, 2026.

## External References

4. H. S. M. Coxeter, *Regular Polytopes*, 3rd ed., Dover, 1973.
5. S. Pancharatnam, "Generalized theory of interference, and its applications," *Proceedings of the Indian Academy of Sciences A*, 44, 247-262, 1956.
6. M. V. Berry, "Quantal phase factors accompanying adiabatic changes," *Proceedings of the Royal Society of London A*, 392, 45-57, 1984. DOI: `10.1098/rspa.1984.0023`.
