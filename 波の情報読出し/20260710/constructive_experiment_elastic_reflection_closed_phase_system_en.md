# Constructive Experiment on Elastic Reflection of Two Fermionic Local Waves in a Closed Phase System Without Assuming Background Space

**Subtitle:** Finite-Resolution Interaction Cells and Internal Observation in the Anonymous Equal-Amplitude Composite-Wave Model  
**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Status:** Constructive experimental paper  
**Version DOI:** 10.5281/zenodo.21291020  
**Concept DOI:** 10.5281/zenodo.21291018  

---

## Abstract

This paper starts from anonymity, all-positive-sign zero closure, nontrivial existence, anonymous equal-amplitude composite waves, and normalized equal-amplitude odd-harmonic localized waves. Within a closed phase system that does not assume background space a priori, we construct and numerically examine a finite-resolution conservative map corresponding to complete elastic reflection of two distinguishable local waves with fermionic cores.

The local waves `A` and `B` carry spatial phase, temporal phase, a readout quantity for direction of propagation, representative amplitude, and an identification oscillation mode on an internal identification phase `η`. The observer `C` is implemented as a heavy local wave with sufficiently large representative-amplitude capacity, serving as a quasi-static readout reference that supplies an effective curvature radius near `A` and `B`. A collision is not defined as a point event, but as entry into a finite-resolution cell in spatial and temporal phase. Inside the interaction cell, a map is applied that reverses the direction readout while preserving the identification oscillation, representative amplitude, fermionic core, and compensated square closure.

We examined the minimal experiment, robustness of identification oscillations, observer capacity, cell resolution, control experiments for reflection/transmission/label exchange, asymmetric conditions, observation perturbations, repeated collisions, and readout resolution of the identification phase `η`. The reflection map was distinguished from transmission and label-exchange maps. Even after eight repeated collisions, direction reversal, identification oscillations, representative amplitudes, and closure residuals were preserved. Observer insufficiency, cell overshooting, temporal-cell mismatch, identification-mode leakage, and `η` readout aliasing were separated as failure conditions.

These results do not claim to derive standard fermion scattering, S-matrix scattering theory, or the standard quantum measurement process. They show, as a numerical constructive experiment, that a repeatable conservative map corresponding to complete elastic reflection of two distinguishable local waves can be constructed from closed phase relations, localized waves, internal identification oscillations, and internal observation alone.

**Keywords:** anonymous equal-amplitude composite wave, all-positive-sign zero closure, finite-resolution cell, fermionic local wave, identification oscillation, complete elastic reflection, internal observation

---

## 1. Introduction

In standard particle scattering, background space, asymptotic states, interaction potentials, and scattering amplitudes occupy the center of the theoretical construction. This paper does not adopt that construction. Instead, it asks whether localized waves, identification oscillations, an observer, and a finite-resolution cell can be constructed inside a closed phase system so that two local waves may be read as having undergone complete elastic reflection.

The question of this paper is:

> Can complete elastic reflection of two distinguishable fermionic local waves, together with their identification before and after observation, be constructed using only the interior of a closed phase system without positing background space a priori?

To address this question, the paper constructs the following sequence:

1. Anonymity and all-positive-sign zero closure are placed at the first-principle layer.
2. Anonymous equal-amplitude composite waves are used as the basic structure of local waves.
3. Localization in spatial and temporal phase is constructed by normalized equal-amplitude odd harmonics, rather than by an external window function.
4. The two local waves `A` and `B` are assigned different identification oscillations on an internal identification phase `η`.
5. The heavy observer `C` is treated not as background space, but as a quasi-static curvature-radius generator and readout reference.
6. Collision is defined not as a point event, but as a reversible map across a finite-resolution cell.
7. The reflection map is distinguished from transmission and label-exchange maps.

External references are used minimally as coordinate axes for comparison with standard contexts, not as derivational foundations of the present construction. Pauli [Pauli1925] and Dirac [Dirac1926] are cited for the standard background of fermions and exclusion; Lippmann and Schwinger [LippmannSchwinger1950] for a representative standard formulation of scattering theory; Rovelli [Rovelli1996] as a nearby relational context for observation; and Shannon [Shannon1949] for finite sampling and aliasing.

---

## 2. Claims Not Made in This Paper

This paper does not claim the following.

| Claim not made | Reason |
|---|---|
| A derivation of standard fermion scattering | No potential scattering, S-matrix, or asymptotic states are used |
| A derivation of the standard quantum measurement process | Observation is defined as a correlation readout inside this axiom system |
| Calculation of physical collision cross sections of real particles | No cross section, experimental units, or standard Hamiltonian is introduced |
| Differential-energy conservation in a standard wave equation | The conserved structure is compensated square closure based on Axiom 1 |
| Direct observation of a collision instant as a point event | Collision is treated as a map before and after a finite-resolution cell |

The claim of this paper is limited to constructive possibility: whether a conservative map inside a closed phase system can be constructed from a small set of axioms.

---

## 3. Basic Axiom System

The experiment is based on the Anonymous Equal-Amplitude Composite-Wave Model, Basic Axiom System v3.

### 3.1 Anonymity

Basic components are not given individual names, privileged axes, privileged signs, or privileged types. Indices are convenient labels and are not intrinsic names of components.

### 3.2 All-Positive-Sign Zero Closure

The closure condition is:

```math
\sum_{n=1}^{N}x_n^2=0
```

This is not the conjugate norm

```math
\sum_n |x_n|^2.
```

The closure uses a square-register form `x_n^2`, not `x_n \bar{x}_n`.

### 3.3 Nontrivial Existence

The closure system is not restricted to the all-zero solution.

```math
\exists n,\quad x_n\neq 0
```

### 3.4 Closure by Compensating Phase

In the experimental implementation, each wave coefficient `x_n` is paired with a compensating phase `i x_n`.

```math
x_n^2+(i x_n)^2=0
```

Thus the squared sum of each component pair is zero.

---

## 4. Construction of Localized Waves

### 4.1 Normalized Equal-Amplitude Odd-Harmonic Localized Wave

Let the highest odd-harmonic order be `N_h=2K-1`. Define the normalized equal-amplitude odd-harmonic localized wave by

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u).
```

It can also be written in closed form as

```math
S_{N_h}(u)=\frac{\sin((N_h+1)u)}{(N_h+1)\sin u}.
```

At the central phase `u=0`,

```math
S_{N_h}(0)=1.
```

Therefore, a waveform with representative amplitude `A`,

```math
\Psi_{N_h}(u)=A S_{N_h}(u),
```

satisfies

```math
\Psi_{N_h}(0)=A.
```

Increasing the highest odd-harmonic order therefore does not cause the representative amplitude at in-phase alignment to diverge.

### 4.2 Spatial and Temporal Localized Waves

Let `χ` be spatial phase and `τ` be temporal phase. A local wave `P` is given by

```math
P(\chi,\tau)
=
A_P
S_{N_{h,\chi}^{P}}(\chi-\chi_P)
S_{N_{h,\tau}^{P}}(\tau-\tau_P)
e^{i\phi_P}.
```

Localization is produced by the odd-harmonic structure itself, not by an external window function.

---

## 5. Identification Oscillation and Fermionic Core

### 5.1 Fermionic Core

A two-component opposite-phase structure is read as a fermionic core.

```math
N=2,\qquad \Delta_{12}=\pi
```

```math
e^{i\theta_1}+e^{i\theta_2}=0
```

This structure is not identified with the standard antisymmetric many-particle wave function. In this paper, "fermionic" is a working classification for a local wave that contains an internal opposite-phase core.

### 5.2 Internal Identification Phase `η`

To avoid distinguishing the local waves `A` and `B` merely by variable names, an internal identification phase `η` is introduced.

The identification oscillation is

```math
D_{m_P}(\eta)=e^{i m_P\eta}.
```

`A` and `B` carry different identification oscillation modes.

```math
m_A\neq m_B
```

In the minimal implementation,

```math
m_A=1,\qquad m_B=2.
```

### 5.3 Local Wave Including Identification Oscillation

The local wave including the identification phase is

```math
P(\chi,\tau,\eta)
=
A_P
S_{N_{h,\chi}^{P}}(\chi-\chi_P)
S_{N_{h,\tau}^{P}}(\tau-\tau_P)
D_{m_P}(\eta)
e^{i\phi_P}.
```

Thus the identifier is not an externally attached name, but an internal oscillation mode embedded in the waveform.

### 5.4 Readout of Identification Oscillation

For a readout mode `m`, define the reference wave

```math
D_{\mathrm{read},m}(\eta)=e^{-im\eta}.
```

The identification readout correlation is

```math
O_{P,m}^{(j)}
=
\frac{1}{(2\pi)^3}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
P(\chi,\tau,\eta)
C_{\mathrm{read},P}^{(j)}(\chi,\tau)
D_{\mathrm{read},m}(\eta)
d\chi d\tau d\eta.
```

The identification mode is read as the value of `m` that maximizes `|O_{P,m}^{(j)}|`.

---

## 6. Observer `C` and Relative-Phase Readout

### 6.1 Heavy Observer

The total system is treated as a closed system consisting of `A`, `B`, and `C`.

```math
\sum_A x_n^2+\sum_B x_n^2+\sum_C x_n^2=0
```

The observer `C` is sufficiently heavy and is read as a curvature-radius generator.

```math
R_C^2=-\sum_C x_n^2
```

Thus, near `A` and `B`,

```math
\sum_A x_n^2+\sum_B x_n^2=R_C^2
```

is used.

### 6.2 Reference Wave for Relative-Phase Readout

The observer body `C^{(j)}` is distinguished from the inverse-phase reference wave `C_{\rm read}^{(j)}` used for relative-phase readout.

```math
C_{\rm read}^{(j)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C^{(j)})
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(j)})
e^{-i\phi_C}.
```

Then

```math
P(\chi,\tau)C_{\rm read}^{(j)}(\chi,\tau)
\propto e^{i(\phi_P-\phi_C)}
```

reads the relative phase with respect to the external reference.

---

## 7. Finite-Resolution Interaction Cell and Reflection Map

### 7.1 Finite-Resolution Cell

A collision is defined not as a point event, but as entry into a finite-resolution cell.

```math
|\chi_A-\chi_B|<\epsilon_\chi^{AB}
```

```math
|\tau_A-\tau_B|<\epsilon_\tau^{AB}
```

The cell width is fixed by the width from the central peak of the odd-harmonic localized wave to the first zero.

```math
\epsilon_\chi^{AB}
=
\frac{\pi}
{\min(N_{h,\chi}^{A},N_{h,\chi}^{B})+1}
```

```math
\epsilon_\tau^{AB}
=
\frac{\pi}
{\min(N_{h,\tau}^{A},N_{h,\tau}^{B})+1}
```

### 7.2 Computational Order Parameter

The numerical implementation introduces a computational order parameter `s`. It is the order of the computational procedure and is not the readout temporal phase `τ` itself.

```math
s\neq \tau
```

Spatial phase is updated by

```math
\chi_A(s+\Delta s)=\chi_A(s)+q_A(s)v_\chi\Delta s
```

```math
\chi_B(s+\Delta s)=\chi_B(s)+q_B(s)v_\chi\Delta s
```

and temporal phase by

```math
\tau_A(s+\Delta s)=\tau_A(s)+\omega_A\Delta s
```

```math
\tau_B(s+\Delta s)=\tau_B(s)+\omega_B\Delta s
```

### 7.3 Complete Elastic Reflection Map

Inside the finite-resolution cell, apply the map

```math
\mathcal C:
(A',B',C')\mapsto(A'',B'',C'').
```

The direction readout is reversed:

```math
q_A''=-q_A'
```

```math
q_B''=-q_B'
```

while identification oscillation, representative amplitude, and fermionic core are preserved:

```math
m_A''=m_A',\qquad m_B''=m_B'
```

```math
A_A''=A_A',\qquad A_B''=A_B'
```

Applying the map twice returns the system to the original state:

```math
\mathcal C^2=\mathrm{id}.
```

---

## 8. Numerical Experiments

### 8.1 Basic Parameters

| Quantity | Value | Meaning |
|---|---:|---|
| `A_A` | `1` | Representative amplitude of local wave A |
| `A_B` | `1` | Representative amplitude of local wave B |
| `A_C` | `1000` | Representative amplitude of observer C |
| `N_{h,chi}^A` | `99` | Highest odd-harmonic order of A in spatial direction |
| `N_{h,chi}^B` | `99` | Highest odd-harmonic order of B in spatial direction |
| `N_{h,chi}^C` | `999` | Highest odd-harmonic order of C in spatial direction |
| `N_{h,tau}^A` | `99` | Highest odd-harmonic order of A in temporal direction |
| `N_{h,tau}^B` | `99` | Highest odd-harmonic order of B in temporal direction |
| `N_{h,tau}^C` | `999` | Highest odd-harmonic order of C in temporal direction |
| `chi_A^{(0)}` | `-0.2` | Initial spatial phase of A |
| `chi_B^{(0)}` | `+0.2` | Initial spatial phase of B |
| `tau_A^{(0)},tau_B^{(0)}` | `-0.2,-0.2` | Initial temporal phases of A and B |
| `q_A^{(0)},q_B^{(0)}` | `+1,-1` | Initial direction readouts |
| `m_A,m_B` | `1,2` | Internal identification oscillation modes |
| `delta_s` | `0.01` | Computational order step |

### 8.2 Experimental Set

| No. | Experiment | Purpose |
|---:|---|---|
| 1 | Basic complete elastic collision | Test whether the map holds under minimal conditions |
| 2 | Robustness of identification oscillation | Test tolerance to identification-mode mixing |
| 3 | Observer C condition | Test the heaviness condition of C |
| 4 | Cell resolution | Test relation between step size and finite cell width |
| 5 | Reflection/transmission controls | Distinguish reflection from transmission and label exchange |
| 6 | Asymmetric conditions | Test amplitude, harmonic, and temporal asymmetries |
| 7 | Observation perturbation | Test effects of observation-induced phase perturbations |
| 8 | Multiple collisions | Test repeated preservation |
| 9 | `η` readout resolution | Test aliasing conditions of identification oscillations |

---

## 9. Results

### 9.1 Basic Complete Elastic Collision

In the minimal experiment, `A` and `B` reached the collision cell, direction readouts reversed, and identification oscillations and representative amplitudes were preserved.

| Item | Result |
|---|---|
| Collision-cell arrival step | `19` |
| Final step | `40` |
| Reversal of `q_A,q_B` | `true,true` |
| Preservation of `m_A,m_B` | `true,true` |
| Label exchange | `false` |
| Crosstalk | `false` |
| Compensated square-closure residual | `0.0` |
| Overall verdict | `true` |

**Figure 1. Trajectories and direction readouts in the basic collision**

![Basic collision trajectory](elastic_collision_simulation_result_v1/elastic_collision_trajectory_v1.png)

**Figure 2. Readout of identification oscillation modes**

![Identification mode readout](elastic_collision_simulation_result_v1/elastic_collision_identification_modes_v1.png)

### 9.2 Robustness of Identification Oscillation

When leakage was introduced into the identification oscillation, the identification mode was preserved for small leakage. As leakage increased, the purity condition failed.

| Item | Value |
|---|---:|
| Total cases | `36` |
| Valid cases | `20` |
| Invalid cases | `16` |
| First failing leakage rate | `0.35` |

**Figure 3. Identification oscillation purity**

![Identification purity](elastic_collision_label_robustness_result_v1/label_robustness_purity_v1.png)

**Figure 4. Identification preservation verdict**

![Identification preservation verdict](elastic_collision_label_robustness_result_v1/label_robustness_detection_v1.png)

### 9.3 Heaviness Condition of Observer C

When the representative-amplitude capacity of observer `C` was too small, the quasi-static condition failed.

| Item | Value |
|---|---:|
| Total cases | `13` |
| Valid cases | `7` |
| Invalid cases | `6` |
| First valid `A_C` | `20` |

**Figure 5. Observer condition quantities**

![Observer conditions](elastic_collision_observer_sweep_result_v1/observer_sweep_conditions_v1.png)

**Figure 6. Validity by observer capacity**

![Observer validity](elastic_collision_observer_sweep_result_v1/observer_sweep_validity_v1.png)

### 9.4 Cell Resolution and Time Step

When the update step was coarser than the finite cell width, the collision cell could be skipped.

| Item | Value |
|---|---:|
| Total cases | `60` |
| Valid cases | `57` |
| Invalid cases | `3` |
| Off-grid valid cases | `27` |
| Off-grid invalid cases | `3` |

**Figure 7. Sampling condition**

![Sampling condition](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sampling_condition_v1.png)

**Figure 8. Validity under off-grid conditions**

![Off-grid validity](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_validity_d0_0_203_v1.png)

### 9.5 Control Experiments: Reflection, Transmission, and Label Exchange

Reflection, transmission, and label-exchange maps were compared.

| Map | Verdict |
|---|---|
| reflection | valid |
| transmission | invalid |
| label_exchange_reflection | invalid |
| transmission_with_label_exchange | invalid |

Only the reflection map simultaneously satisfied direction reversal and preservation of identification oscillations.

**Figure 9. Control-map trajectories**

![Control map trajectories](elastic_collision_control_maps_result_v1/control_maps_trajectories_v1.png)

**Figure 10. Control-map verdicts**

![Control map verdicts](elastic_collision_control_maps_result_v1/control_maps_verdict_v1.png)

### 9.6 Asymmetric Conditions

Amplitude differences and harmonic-order differences alone did not break the map. The main cause of failure was loss of simultaneous satisfaction of the spatial and temporal cells.

| Item | Value |
|---|---:|
| Total cases | `10` |
| Valid cases | `7` |
| Invalid cases | `3` |

**Figure 11. Cell gaps under asymmetric conditions**

![Asymmetric cell gaps](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_cell_gaps_v1.png)

**Figure 12. Asymmetry verdicts**

![Asymmetry verdicts](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_verdict_v1.png)

### 9.7 Observation Perturbation

When observation perturbation exceeded the localization width of `C`, the observation model failed. However, as long as the AB cell condition was satisfied, the collision map itself could still hold.

| Item | Value |
|---|---:|
| Total cases | `8` |
| Observation-model valid cases | `5` |
| Collision-map valid cases | `7` |
| Invalid cases | `3` |

**Figure 13. Observation perturbation thresholds**

![Observation perturbation thresholds](elastic_collision_observation_perturbation_result_v1/observation_perturbation_thresholds_v1.png)

**Figure 14. Observation perturbation verdicts**

![Observation perturbation verdicts](elastic_collision_observation_perturbation_result_v1/observation_perturbation_verdict_v1.png)

### 9.8 Multiple Collisions

Within a closed interval, wall reflections were allowed and AB collisions were repeated.

| Item | Value |
|---|---:|
| Target AB collisions | `8` |
| Actual AB collisions | `8` |
| Wall reflections | `14` |
| All identification modes preserved | `true` |
| `q` reversed at each collision | `true` |
| Representative amplitudes preserved | `true` |
| Closure preserved | `true` |
| Overall verdict | `true` |

**Figure 15. Trajectories in repeated collisions**

![Multiple collision trajectory](elastic_collision_multi_collision_result_v1/multi_collision_trajectory_v1.png)

**Figure 16. Closure residual in repeated collisions**

![Multiple collision closure residual](elastic_collision_multi_collision_result_v1/multi_collision_closure_v1.png)

### 9.9 Readout Resolution of `η` Identification Oscillation

When the internal identification oscillation was read with a finite number of samples, aliasing occurred if the mode difference was a multiple of the sample count.

| Item | Value |
|---|---:|
| Total cases | `88` |
| Valid cases | `59` |
| Invalid cases | `29` |
| Aliasing cases | `29` |
| Non-aliasing failures | `0` |
| Minimum `η` sample count valid for all mode pairs | `64` |

All invalid cases were due to `η` readout aliasing. This is not a failure of the collision map itself, but insufficient readout resolution for the identification oscillation.

**Figure 17. Validity of `η` readout resolution**

![Eta validity](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_validity_v1.png)

**Figure 18. `η` identification purity**

![Eta purity](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_purity_v1.png)

---

## 10. Discussion

### 10.1 Distinguishing Reflection from Transmission

If one observes only position trajectories, reflection and transmission can often be difficult to distinguish. This paper tracks the identification oscillations `m_A,m_B` and the direction readouts `q_A,q_B` simultaneously.

Transmission and label-exchange maps can resemble reflection at the level of position readout, but are rejected by the conditions of direction reversal or identification-mode preservation.

Thus the reflection criterion is not merely an apparent trajectory. It includes the condition:

```text
the local wave carrying its identification oscillation is pushed out in the reversed direction
```

### 10.2 Local Interaction Without Background Space

This paper does not place background space at the start. Position is treated not as an external spatial coordinate, but as a relative-phase readout with respect to the observer.

The heavy observer `C` is not fully solved as a dynamical third body. Instead, it is absorbed as a quasi-static reference that supplies an effective curvature radius near AB. This yields a local linear approximation from inside the closed system.

The variables `χ` and `τ` used in this paper are not a priori physical spacetime coordinates in which local waves are placed. They are readout variables for spatial and temporal phase in a closed phase system, and serve as computational domains for expressing positional and temporal order through relative phase with the observer. Thus, when this paper says that no background space is assumed, it does not deny the use of a parameter domain for describing phase relations. It means that an external spacetime stage existing independently of particles and observers is not taken as the starting point.

### 10.3 Nature of the Conserved Quantity

The conserved quantity in this paper is not the differential energy of a standard wave equation. The basis of conservation is the compensated square register associated with all-positive-sign zero closure.

```math
\sum_n x_n^2=0
```

By pairing each component with `x_n` and `i x_n`,

```math
x_n^2+(i x_n)^2=0
```

holds. In the multiple-collision experiment, this closure residual remained zero after repeated collisions.

Therefore, zero closure residual does not mean that a closure law was accidentally discovered by numerical simulation. It confirms that the defined compensated square-closure structure was preserved under the collision map and repeated updates.

### 10.4 Meaning of Failure Conditions

The model does not succeed under all conditions. The map or its observation fails when the observer is too light, the cell is skipped, the spatial and temporal cells are not simultaneously satisfied, the identification oscillation is mixed too strongly, or aliasing occurs in `η` readout.

This means that the construction is not an arbitrary model that succeeds by definition. It has clear validity conditions and failure conditions.

### 10.5 Status of the Reversal Map

The direction-reversal map itself is given as a constructive rule inside the finite-resolution interaction cell. It is not uniquely derived from the first principles alone. What this experiment shows is that this reversal rule is compatible with identification oscillations, representative amplitudes, fermionic cores, internal observation, and compensated square closure, and that it constructs a repeatable conservative map distinguishable from transmission and label-exchange maps.

---

## 11. Validity and Failure Conditions

### 11.1 Validity Conditions

| Condition | Content |
|---|---|
| Finite-cell arrival | `A,B` simultaneously satisfy the spatial and temporal cells |
| Update step | `delta_s` does not skip over the cell width |
| Observer condition | `C` is sufficiently heavy to be treated as a quasi-static curvature-radius generator |
| Identification oscillation | `m_A,m_B` are separable at the readout resolution |
| Identification purity | Leakage does not invert the dominant mode |
| Reflection map | `q_A,q_B` reverse, while `m_A,m_B` and representative amplitudes are preserved |
| Closure | The compensated square register has zero residual |

### 11.2 Failure Conditions

| Failure condition | Experimental manifestation |
|---|---|
| Observer too light | The quasi-static condition fails for `A_C < 20` |
| Update step too coarse | Collision cell is not reached under off-grid high-resolution conditions |
| Temporal-cell mismatch | Spatial proximity does not imply entry into the interaction cell |
| Large observation perturbation | Observation model fails when C width is exceeded; map fails when AB cell is exceeded |
| Identification leakage | Purity condition fails as leakage increases |
| `η` aliasing | Modes become indistinguishable when mode difference is a multiple of sample count |
| Transmission or label exchange | Rejected by direction reversal or identification-mode preservation |

---

## 12. Conclusion

Starting from a small set of working axioms and definitions--anonymity, all-positive-sign zero closure, nontrivial existence, anonymous equal-amplitude composite waves, and normalized equal-amplitude odd-harmonic localization--this paper constructed a finite-resolution conservative map corresponding to complete elastic reflection of two distinguishable fermionic local waves in a closed phase system without assuming background space a priori.

The two local waves carried different oscillation modes on an internal identification phase. Upon entering a finite-resolution interaction cell, their direction readouts were reversed while their identification oscillations, representative amplitudes, and fermionic cores were preserved. Using spatial, temporal, and identification-phase correlation readout by a heavy observer, the reflection map was distinguished from transmission and label-exchange maps.

Validity and failure conditions were obtained for observer amplitude capacity, update step, temporal-cell alignment, identification purity, and sampling resolution in the identification phase. In repeated collisions, direction reversal, identification oscillation, representative amplitude, and compensated square closure were preserved.

These results do not constitute a derivation of standard fermion scattering or of the quantum measurement process. They are a numerical constructive experiment showing that a repeatable conservative map corresponding to complete elastic reflection of distinguishable local waves can be constructed from closed phase relations, localized waves, internal identification oscillations, and internal observation, without placing background space as an external stage.

---

## Acknowledgments

This paper is based on the Anonymous Equal-Amplitude Composite-Wave Model Basic Axiom System v3, the Complete Elastic Collision Simulation Specification v1, and the Complete Elastic Collision Simulation Experiment Results v1, all prepared on 2026-07-10.

---

## References

### Internal Documents

[KiharaAxioms2026] Noriaki Kihara, [Anonymous Equal-Amplitude Composite-Wave Model Basic Axiom System v3](basic_axiom_system_v3_en.md), 2026-07-10.

[KiharaSpec2026] Noriaki Kihara, [Complete Elastic Collision Simulation Specification v1](elastic_collision_simulation_spec_v1_en.md), 2026-07-10.

[KiharaResults2026] Noriaki Kihara, [Complete Elastic Collision Simulation Experiment Results v1](elastic_collision_simulation_experiment_results_v1_en.md), 2026-07-10.

### External References

[Pauli1925] W. Pauli, "Über den Zusammenhang des Abschlusses der Elektronengruppen im Atom mit der Komplexstruktur der Spektren," *Zeitschrift für Physik*, 31, 765-783, 1925.

[Dirac1926] P. A. M. Dirac, "On the Theory of Quantum Mechanics," *Proceedings of the Royal Society of London. Series A*, 112, 661-677, 1926.

[LippmannSchwinger1950] B. A. Lippmann and J. Schwinger, "Variational Principles for Scattering Processes. I," *Physical Review*, 79, 469-480, 1950.

[Rovelli1996] C. Rovelli, "Relational Quantum Mechanics," *International Journal of Theoretical Physics*, 35, 1637-1678, 1996. arXiv:quant-ph/9609002.

[Shannon1949] C. E. Shannon, "Communication in the Presence of Noise," *Proceedings of the IRE*, 37(1), 10-21, 1949.

[Noether1918] E. Noether, "Invariante Variationsprobleme," *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 235-257, 1918. English translation: "Invariant Variation Problems."

---

## Appendix A. List of Experiment Result Files

| Experiment | Report | Result JSON |
|---|---|---|
| Basic complete elastic collision | [report](elastic_collision_simulation_result_v1/elastic_collision_report_v1.md) | [json](elastic_collision_simulation_result_v1/elastic_collision_result_v1.json) |
| Robustness of identification oscillation | [report](elastic_collision_label_robustness_result_v1/label_robustness_report_v1.md) | [json](elastic_collision_label_robustness_result_v1/label_robustness_result_v1.json) |
| Observer C condition | [report](elastic_collision_observer_sweep_result_v1/observer_sweep_report_v1.md) | [json](elastic_collision_observer_sweep_result_v1/observer_sweep_result_v1.json) |
| Cell resolution | [report](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_report_v1.md) | [json](elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_result_v1.json) |
| Control maps | [report](elastic_collision_control_maps_result_v1/control_maps_report_v1.md) | [json](elastic_collision_control_maps_result_v1/control_maps_result_v1.json) |
| Asymmetric conditions | [report](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_report_v1.md) | [json](elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_result_v1.json) |
| Observation perturbation | [report](elastic_collision_observation_perturbation_result_v1/observation_perturbation_report_v1.md) | [json](elastic_collision_observation_perturbation_result_v1/observation_perturbation_result_v1.json) |
| Multiple collisions | [report](elastic_collision_multi_collision_result_v1/multi_collision_report_v1.md) | [json](elastic_collision_multi_collision_result_v1/multi_collision_result_v1.json) |
| `η` resolution | [report](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_report_v1.md) | [json](elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_result_v1.json) |
