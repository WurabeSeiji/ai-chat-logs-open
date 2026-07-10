# Complete Elastic Collision Simulation Specification v1 for Fermionic Local Waves A,B and a Heavy Observer C

**Date:** 2026-07-10  
**Author:** Noriaki Kihara  
**Status:** Simulation specification based on the Anonymous Equal-Amplitude Composite-Wave Model Basic Axiom System v3  

---

## 0. Confirmed Axiom System Assumed

This specification assumes the following confirmed axiom system.

- Anonymous Equal-Amplitude Composite-Wave Model Basic Axiom System v3
- Axiom 0: Anonymity
- Axiom 1: All-positive-sign zero closure

```math
\sum_n x_n^2=0
```

- Axiom 2: Nontrivial existence
- Axiom 3: Anonymous equal-amplitude composite wave
- Chapter 6: Normalized equal-amplitude odd-harmonic localized wave

This specification does not claim to derive standard scattering theory, standard particle collisions, or standard quantum measurement. Its purpose is to convert the complete elastic collision of local waves with fermionic cores into a numerically simulable finite-resolution map inside the confirmed axiom system v3.

---

## 1. Purpose of the Simulation Experiment

### 1.1 Purpose

The purpose is to test whether two particle-like local waves `A` and `B`, each carrying a fermionic core, satisfy the following properties before and after a finite-resolution cell in a central-projection-like local neighborhood defined by a heavy observer `C`.

1. Preservation of identification oscillation
2. Preservation of representative amplitude
3. Preservation of the fermionic core
4. Reversal of the direction readout
5. Quasi-static behavior of observer `C`
6. Stability of the curvature radius based on all-positive-sign zero closure

Complete elastic collision here is not a point event. It is defined as a reversible map before and after a finite-resolution cell.

```math
(A(+;L_A),B(-;L_B))
\longmapsto
(A(-;L_A),B(+;L_B))
```

### 1.2 Claims Not Made

This specification does not claim:

- derivation of standard fermion scattering;
- derivation of the standard quantum-mechanical measurement process;
- calculation of physical cross sections of real particles;
- differential-energy conservation in a standard wave equation;
- direct observation of a collision instant as a point event.

The target is a finite-resolution map internal to the confirmed axiom system v3.

---

## 2. Geometric Model

### 2.1 Central Projection Model

The total system is treated as a closed system consisting of a heavy observer `C` and particle-like local waves `A` and `B`.

```math
\sum_A x_n^2+\sum_B x_n^2+\sum_C x_n^2=0
```

The observer `C` is sufficiently heavy and is read as a curvature-radius generator.

```math
R_C^2=-\sum_C x_n^2
```

Therefore, in the neighborhood of `A` and `B`,

```math
\sum_A x_n^2+\sum_B x_n^2=R_C^2
```

is used.

When the interaction region of `A` and `B` is sufficiently small compared with `R_C`, a local linear approximation is used.

```math
\ell_{AB}\ll R_C
```

### 2.2 Spatial Phase, Temporal Phase, and Identification Phase

Each local wave has spatial phase `χ`, temporal phase `τ`, and internal identification phase `η`.

```math
A=(A_A,N_{h,\chi}^A,N_{h,\tau}^A,\chi_A,\tau_A,\phi_A,m_A,q_A)
```

```math
B=(A_B,N_{h,\chi}^B,N_{h,\tau}^B,\chi_B,\tau_B,\phi_B,m_B,q_B)
```

```math
C=(A_C,N_{h,\chi}^C,N_{h,\tau}^C,\chi_C,\tau_C,\phi_C)
```

where:

- `A_A,A_B,A_C`: representative amplitudes;
- `N_{h,χ}`: highest odd-harmonic order in the spatial direction;
- `N_{h,τ}`: highest odd-harmonic order in the temporal direction;
- `χ`: spatial/position phase;
- `τ`: temporal phase;
- `φ`: internal reference phase;
- `η`: internal identification phase;
- `m_A,m_B`: identification oscillation modes;
- `q_A,q_B`: direction readout quantities.

### 2.3 Central Projection Diagram

The following diagram shows the concept that a heavy observer `C` generates curvature radius `R_C` and local waves `A,B` are placed in its neighborhood.

![Central projection model](完全弾性衝突_中心投影モデル.png)

---

## 3. Local Waveforms Used

### 3.1 Normalized Equal-Amplitude Odd-Harmonic Localized Wave

Let the highest odd-harmonic order be

```math
N_h=2K-1
```

The normalized equal-amplitude odd-harmonic localized wave is

```math
S_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}\cos((2m+1)u)
```

or, in complex representation,

```math
Z_{N_h}(u)=\frac{1}{K}\sum_{m=0}^{K-1}e^{i(2m+1)u}.
```

A localized wave with representative amplitude `A` is

```math
\Psi_{N_h}(u)=A S_{N_h}(u)
```

or

```math
\Psi_{N_h}(u)=A Z_{N_h}(u).
```

At the central phase `u=0`,

```math
S_{N_h}(0)=1
```

and therefore

```math
\Psi_{N_h}(0)=A.
```

Thus increasing `N_h` preserves the representative amplitude at in-phase alignment.

### 3.2 Spatial, Temporal, and Identification Localized Wave

For `P in {A,B,C}`, use spatial phase `χ` and temporal phase `τ`. For `A,B`, also use identification phase `η`.

The observer body is

```math
C(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C)
e^{i\phi_C}.
```

The local waves `A,B` are

```math
P(\chi,\tau,\eta)
=
A_P
S_{N_{h,\chi}^{P}}(\chi-\chi_P)
S_{N_{h,\tau}^{P}}(\tau-\tau_P)
D_{m_P}(\eta)
e^{i\phi_P},
\qquad P\in\{A,B\}.
```

The identification oscillation is

```math
D_{m_P}(\eta)=e^{i m_P\eta}.
```

The readout reference for identification mode `m` is

```math
D_{\rm read,m}(\eta)=e^{-im\eta}.
```

In the minimal implementation,

```math
m_A=1,\qquad m_B=2.
```

These are not intrinsic names attached from outside; they are different internal oscillation modes on the identification phase.

---

## 4. Simulation Assumptions

### 4.1 Assumption 1: C Is Sufficiently Heavy

The observer `C` has a representative-amplitude capacity much larger than that of local waves `A` and `B`.

```math
A_C\gg A_A,A_B
```

Thus `C` can be treated as a quasi-static curvature-radius generator.

```math
R_C^2=-\sum_C x_n^2
```

### 4.2 Assumption 2: Local Linear Approximation Near AB

The distance between `A` and `B` is read as

```math
\ell_{AB}=R_C|\chi_A-\chi_B|.
```

When

```math
\ell_{AB}\ll R_C
```

holds, the phase arc near AB is approximated linearly.

### 4.3 Assumption 3: Position Phase of C Is Almost Unchanged by Observation

The change in the position phase of `C` due to observation is negligible compared with changes in `A` and `B`.

```math
|\delta\chi_C|\ll |\delta\chi_A|,|\delta\chi_B|
```

Similarly for temporal phase:

```math
|\delta\tau_C|\ll |\delta\tau_A|,|\delta\tau_B|.
```

### 4.4 Assumption 4: Positive Feedback Divergence Is Suppressed

The positive-feedback loop created by mutual readout between observer `C` and local waves `A,B` is suppressed by the sufficiently large representative-amplitude capacity of `C`.

Conceptually,

```math
G_{\rm loop}\sim\frac{A_A+A_B}{A_C}
```

and the assumption is

```math
G_{\rm loop}\ll1.
```

### 4.5 Assumption 5: Observation Is Also Localized in the Temporal Direction

The observer `C` is not a continuous reference wave in the temporal direction. It is localized in temporal phase by odd-harmonic structure.

For readout center `j`,

```math
C^{(j)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(j)})
e^{i\phi_C}.
```

Thus the temporal phase of observation is defined at finite resolution by `τ_C^{(j)}` and `N_{h,τ}^C`.

For relative-phase readout, distinguish the observer body `C^{(j)}` from the inverse-phase reference wave `C_{\rm read}^{(j)}`:

```math
C_{\rm read}^{(j)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(j)})
e^{-i\phi_C}.
```

This is not the introduction of a conjugate norm. It defines a reference wave with inverse phase as a readout basis. Then

```math
P(\chi,\tau)C_{\rm read}^{(j)}(\chi,\tau)
\propto
e^{i(\phi_P-\phi_C)}
```

reads the relative phase to the external reference.

---

## 5. Initial Conditions

### 5.1 Representative Values

Recommended initial values for the minimal implementation are:

| Quantity | Recommended value | Meaning |
|---|---:|---|
| `A_A` | `1` | Representative amplitude of local wave A |
| `A_B` | `1` | Representative amplitude of local wave B |
| `A_C` | `1000` | Representative amplitude of heavy observer C |
| `N_{h,\chi}^A` | `99` | Spatial highest odd-harmonic order of A |
| `N_{h,\chi}^B` | `99` | Spatial highest odd-harmonic order of B |
| `N_{h,\chi}^C` | `999` | Spatial highest odd-harmonic order of C |
| `N_{h,\tau}^A` | `99` | Temporal highest odd-harmonic order of A |
| `N_{h,\tau}^B` | `99` | Temporal highest odd-harmonic order of B |
| `N_{h,\tau}^C` | `999` | Temporal highest odd-harmonic order of C |
| `χ_C` | `0` | Position-phase reference of C |
| `τ_C^{ref}` | `0` | Temporal-phase reference of C |
| `φ_C` | `0` | Internal reference phase of C |
| `Δs` | `0.01` | Computational order step |
| `v_χ` | `1` | Update speed of position phase |
| `ω_A,ω_B` | `1` | Temporal-phase update rates of A and B |
| `s_max` | `10000` | Maximum number of update steps |

### 5.2 Initial Configuration of A and B

Place `A` and `B` near `C` in the local linear approximation:

```math
\chi_A^{(0)}=-d_0
```

```math
\chi_B^{(0)}=+d_0
```

```math
q_A^{(0)}=+1
```

```math
q_B^{(0)}=-1
```

The identification modes are different:

```math
m_A\neq m_B.
```

Here `q_A=+1` means that A moves toward B, and `q_B=-1` means that B moves toward A.

Recommended value:

```math
d_0=0.2.
```

### 5.3 Initial Temporal Phase

The center temporal phase of pre-collision observation is

```math
\tau^{(0)}=-\tau_0.
```

The collision-cell neighborhood is

```math
\tau^{(1)}=0.
```

The center temporal phase of post-collision observation is

```math
\tau^{(2)}=+\tau_0.
```

Recommended value:

```math
\tau_0=0.2.
```

For the minimal symmetric experiment,

```math
\tau_A^{(0)}=\tau_B^{(0)}=-\tau_0
```

and

```math
\omega_A=\omega_B.
```

Thus the zero temporal-phase difference between `A` and `B` is preserved:

```math
\tau_A(s)-\tau_B(s)=0.
```

The observer `C` is localized in time. Therefore, the same `C` is read with different temporal centers for the pre-collision and post-collision observations:

```math
\tau_C^{(0)}=-\tau_0,\qquad
\tau_C^{(2)}=+\tau_0.
```

The inverse-phase reference waves for relative-phase readout are

```math
C_{\rm read}^{(0)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(0)})
e^{-i\phi_C}
```

and

```math
C_{\rm read}^{(2)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(2)})
e^{-i\phi_C}.
```

This does not mean that two observers are introduced. It means that the same heavy observer `C` is read at different temporal phase centers.

### 5.4 Initial Curvature Radius

Conceptually, define

```math
R_C^2=-\sum_C x_n^2.
```

In the minimal numerical implementation, choose units so that

```math
R_C=A_C.
```

With the recommended initial value,

```math
R_C=1000.
```

The initial distance is

```math
\ell_{AB}^{(0)}=R_C|\chi_A^{(0)}-\chi_B^{(0)}|.
```

For the local linear approximation, use the normalized phase-coordinate distance

```math
|\chi_A^{(0)}-\chi_B^{(0)}|=2d_0=0.4.
```

No conversion to physical length is performed in this specification.

### 5.5 Computational Order Parameter and Position-Phase Update

Introduce the computational order parameter `s`. It is the order of computation and is not the readout temporal phase `τ`.

```math
s\neq\tau
```

The minimal position-phase update rule is

```math
\chi_A(s+\Delta s)=\chi_A(s)+q_A(s)v_\chi\Delta s
```

```math
\chi_B(s+\Delta s)=\chi_B(s)+q_B(s)v_\chi\Delta s.
```

If temporal phase is updated,

```math
\tau_A(s+\Delta s)=\tau_A(s)+\omega_A\Delta s
```

```math
\tau_B(s+\Delta s)=\tau_B(s)+\omega_B\Delta s.
```

This update brings `A` and `B` from the initial configuration toward the finite-resolution cell.

---

## 6. Initial Observation

### 6.1 Observation Formula

Observation is defined as spatial-temporal correlation between local waves `A,B` and observer `C`.

The integration domain is the closed phase domain

```math
\chi,\tau\in[-\pi,\pi)
```

and the observation value is normalized by `(2π)^2`.

```math
O_A^{(0)}
=
\frac{1}{(2\pi)^2}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
A(\chi,\tau) C_{\rm read}^{(0)}(\chi,\tau)
d\chi d\tau
```

```math
O_B^{(0)}
=
\frac{1}{(2\pi)^2}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
B(\chi,\tau) C_{\rm read}^{(0)}(\chi,\tau)
d\chi d\tau.
```

This product is not a conjugate norm. It is a correlation readout matched to the square-register side of Axiom 1.

For relative-phase readout,

```math
A(\chi,\tau)C_{\rm read}^{(0)}(\chi,\tau)
\propto
e^{i(\phi_A-\phi_C)}
```

reads the relative phase with respect to the external reference.

When necessary, an orthogonal compensating component `iC_{\rm read}^{(0)}` may be used for I/Q readout. In the minimal implementation, first record `O_A^{(0)}` and `O_B^{(0)}`.

### 6.2 Identification-Oscillation Readout

The local wave identification is not read from a variable name. It is read from the internal identification phase `η`.

Use

```math
D_{\rm read,m}(\eta)=e^{-im\eta}.
```

The identification-mode correlation is

```math
O_{P,m}^{(j)}
=
\frac{1}{(2\pi)^3}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
P(\chi,\tau,\eta)
C_{{\rm read},P}^{(j)}(\chi,\tau)
D_{\rm read,m}(\eta)
d\chi d\tau d\eta.
```

For `P in {A,B}` and readout time `j in {0,2}`.

The detected identification mode is

```math
m_P^{\rm read,j}
=
\arg\max_m |O_{P,m}^{(j)}|.
```

The identification purity is

```math
\Gamma_P^{(j)}
=
\frac{|O_{P,m_P}^{(j)}|}
{\sum_m |O_{P,m}^{(j)}|}.
```

The minimal pass condition is

```math
m_A^{\rm read,0}=m_A,\qquad m_B^{\rm read,0}=m_B
```

and, after collision,

```math
m_A^{\rm read,2}=m_A,\qquad m_B^{\rm read,2}=m_B.
```

Label exchange is detected when

```math
m_A^{\rm read,2}=m_B,\qquad m_B^{\rm read,2}=m_A.
```

### 6.3 Quantities Read at Initial Observation

The initial observation reads

```math
(\chi_A^{\rm read,0},\tau_A^{\rm read,0},m_A^{\rm read,0},q_A^{\rm read,0})
```

and

```math
(\chi_B^{\rm read,0},\tau_B^{\rm read,0},m_B^{\rm read,0},q_B^{\rm read,0}).
```

Expected readouts:

```math
\chi_A^{\rm read,0}<\chi_B^{\rm read,0}
```

```math
q_A^{\rm read,0}=+1,\qquad q_B^{\rm read,0}=-1
```

```math
m_A^{\rm read,0}=m_A,\qquad m_B^{\rm read,0}=m_B.
```

### 6.4 Observation Update

After observation,

```math
(A,B,C^{(0)})\mapsto(A',B',C').
```

By the heavy-observer assumption,

```math
C'\simeq C^{(0)},\qquad
\chi_C'\simeq\chi_C,\qquad
\tau_C'\simeq\tau_C^{(0)}.
```

Finite-resolution phase perturbations are allowed for `A` and `B`.

```math
\chi_A'=\chi_A+\delta\chi_A^{\rm obs},\qquad
\tau_A'=\tau_A+\delta\tau_A^{\rm obs}
```

```math
\chi_B'=\chi_B+\delta\chi_B^{\rm obs},\qquad
\tau_B'=\tau_B+\delta\tau_B^{\rm obs}
```

The perturbation widths are limited by the localization width of observer `C`.

```math
|\delta\chi_A^{\rm obs}|,|\delta\chi_B^{\rm obs}|
\lesssim
\epsilon_\chi^C
```

```math
|\delta\tau_A^{\rm obs}|,|\delta\tau_B^{\rm obs}|
\lesssim
\epsilon_\tau^C
```

where

```math
\epsilon_\chi^C=\frac{\pi}{N_{h,\chi}^C+1},\qquad
\epsilon_\tau^C=\frac{\pi}{N_{h,\tau}^C+1}.
```

---

## 7. Collision Interaction Calculation

### 7.1 Finite-Resolution Cell

Collision is defined as entry into a finite-resolution cell, not as a point event.

Spatial cell:

```math
|\chi_A'-\chi_B'|<\epsilon_\chi^{AB}
```

Temporal cell:

```math
|\tau_A'-\tau_B'|<\epsilon_\tau^{AB}.
```

The odd-harmonic localized wave can be written as

```math
S_{N_h}(u)=\frac{\sin((N_h+1)u)}{(N_h+1)\sin u}.
```

The phase width from the central peak to the first zero is

```math
u_0=\frac{\pi}{N_h+1}.
```

Thus, in the minimal implementation,

```math
\epsilon_\chi^{AB}
=
\frac{\pi}{\min(N_{h,\chi}^A,N_{h,\chi}^B)+1}
```

```math
\epsilon_\tau^{AB}
=
\frac{\pi}{\min(N_{h,\tau}^A,N_{h,\tau}^B)+1}.
```

### 7.2 Collision Test

After the initial observation, repeat the update rule of Section 5.5 along the computational order parameter `s`. At each step, `A` and `B` are judged to have entered the interaction cell when both conditions hold:

```math
|\chi_A'-\chi_B'|<\epsilon_\chi^{AB}
```

```math
|\tau_A'-\tau_B'|<\epsilon_\tau^{AB}.
```

Then the complete elastic collision map is applied. If the condition is not satisfied within `s_max`, terminate as

```text
collision_cell_not_reached
```

### 7.3 Complete Elastic Collision Map

The collision map is

```math
\mathcal C:
(A',B',C')\mapsto(A'',B'',C'').
```

It satisfies the following.

#### 7.3.1 Direction Reversal

```math
q_A''=-q_A'
```

```math
q_B''=-q_B'
```

#### 7.3.2 Preservation of Identification Oscillation

```math
m_A''=m_A'
```

```math
m_B''=m_B'
```

The internal identification oscillation mode is preserved as part of the local wave, rather than copied as an external label variable.

#### 7.3.3 Representative Amplitude Preservation

```math
A_A''=A_A'
```

```math
A_B''=A_B'
```

#### 7.3.4 Fermionic Core Preservation

The internal fermionic cores of `A` and `B` are preserved.

```math
C_F^A{}''=C_F^A{}'
```

```math
C_F^B{}''=C_F^B{}'
```

Here `C_F` denotes a fermionic core with internal phase difference `Δ=π`.

#### 7.3.5 Quasi-Static Preservation of C

```math
C''\simeq C'
```

```math
\chi_C''\simeq\chi_C'
```

```math
\tau_C''\simeq\tau_C'
```

#### 7.3.6 Reversibility

Applying the same map twice returns the system to its original state.

```math
\mathcal C^2={\rm id}
```

### 7.4 Continued Motion After Collision

After applying the collision map, continue the update with the same computational order parameter `s`.

After direction reversal,

```math
\chi_A(s+\Delta s)=\chi_A(s)+q_A''v_\chi\Delta s
```

```math
\chi_B(s+\Delta s)=\chi_B(s)+q_B''v_\chi\Delta s.
```

For temporal phase,

```math
\tau_A(s+\Delta s)=\tau_A(s)+\omega_A\Delta s
```

```math
\tau_B(s+\Delta s)=\tau_B(s)+\omega_B\Delta s.
```

Continue this update until the neighborhood of the post-collision observation center `τ_C^{(2)}` is reached.

Thus complete elastic collision is read not only as reversal of `q`, but also as the fact that the local waves carrying their identification oscillations move away in the reversed directions.

---

## 8. Post-Collision Observation

### 8.1 Observation Formula

Use the same normalized observation map as in the initial observation. For relative-phase readout, use `C_{\rm read}^{(2)}` with the post-collision temporal center.

```math
O_A^{(2)}
=
\frac{1}{(2\pi)^2}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
A''(\chi,\tau) C_{\rm read}^{(2)}(\chi,\tau)
d\chi d\tau
```

```math
O_B^{(2)}
=
\frac{1}{(2\pi)^2}
\int_{-\pi}^{\pi}
\int_{-\pi}^{\pi}
B''(\chi,\tau) C_{\rm read}^{(2)}(\chi,\tau)
d\chi d\tau.
```

Identification modes are read by the same `η`-mode correlation as in Section 6.2.

### 8.2 Quantities Read After Collision

Read

```math
(\chi_A^{\rm read,2},\tau_A^{\rm read,2},m_A^{\rm read,2},q_A^{\rm read,2})
```

and

```math
(\chi_B^{\rm read,2},\tau_B^{\rm read,2},m_B^{\rm read,2},q_B^{\rm read,2}).
```

Expected results:

```math
q_A^{\rm read,2}=-q_A^{\rm read,0},\qquad
q_B^{\rm read,2}=-q_B^{\rm read,0}
```

```math
m_A^{\rm read,2}=m_A^{\rm read,0},\qquad
m_B^{\rm read,2}=m_B^{\rm read,0}
```

```math
\chi_A^{\rm read,2}<\chi_B^{\rm read,2}
```

```math
|\chi_A^{\rm read,2}-\chi_B^{\rm read,2}|>\epsilon_\chi^{AB}.
```

---

## 9. Judgment Criteria

### 9.1 Conditions for Complete Elastic Collision

The complete elastic collision map is judged to hold when the following conditions are satisfied before and after the finite-resolution cell.

#### Condition 1: Preservation of Identification Oscillation Mode

```math
m_A^{\rm read,2}=m_A^{\rm read,0}
```

```math
m_B^{\rm read,2}=m_B^{\rm read,0}
```

and not exchanged:

```math
m_A^{\rm read,2}\neq m_B^{\rm read,0},\qquad
m_B^{\rm read,2}\neq m_A^{\rm read,0}.
```

Identification purity is sufficiently high:

```math
\Gamma_A^{(0)},\Gamma_A^{(2)},\Gamma_B^{(0)},\Gamma_B^{(2)}
>
\Gamma_{\rm min}.
```

#### Condition 2: Direction Reversal

```math
q_A^{\rm read,2}=-q_A^{\rm read,0}
```

```math
q_B^{\rm read,2}=-q_B^{\rm read,0}
```

#### Condition 3: Separation After Collision

```math
\chi_A^{\rm read,2}<\chi_B^{\rm read,2}
```

```math
|\chi_A^{\rm read,2}-\chi_B^{\rm read,2}|>\epsilon_\chi^{AB}
```

#### Condition 4: Representative Amplitude Preservation

```math
A_A^{\rm read,2}\simeq A_A^{\rm read,0}
```

```math
A_B^{\rm read,2}\simeq A_B^{\rm read,0}
```

#### Condition 5: Quasi-Static C

```math
\chi_C^{\rm read,2}\simeq \chi_C^{\rm read,0}
```

```math
A_C^{\rm read,2}\simeq A_C^{\rm read,0}
```

The temporal phase centers of observation are judged against their specified values, not as quasi-static preservation:

```math
\tau_C^{\rm read,0}=-\tau_0
```

```math
\tau_C^{\rm read,2}=+\tau_0.
```

#### Condition 6: Maintenance of Closure

The closure condition

```math
\sum_A x_n^2+\sum_B x_n^2+\sum_C x_n^2=0
```

is maintained within numerical tolerance.

### 9.2 Failure Conditions

The map is judged invalid if any of the following occurs:

1. identification modes `m_A,m_B` are exchanged or lost;
2. `q_A,q_B` do not reverse;
3. `A,B`, carrying their identification modes, do not separate outside the cell after collision;
4. representative amplitude changes significantly;
5. position phase or representative amplitude of `C` changes significantly;
6. observer temporal centers do not match `τ_C^{(0)}=-τ_0` and `τ_C^{(2)}=+τ_0`;
7. positive feedback causes `A_C,A_A,A_B` to diverge;
8. closure `Σ x_n^2=0` is not maintained;
9. observation cell is too wide to distinguish reflection from transmission;
10. observation cell is too narrow and winding ambiguity dominates;
11. the collision cell is not reached within the maximum update step count;
12. identification-mode crosstalk exceeds the permitted threshold.

---

## 10. Implementation Procedure

### 10.1 Steps

1. Set parameters `A_A,A_B,A_C,N_{h,χ},N_{h,τ},χ,τ,φ,m,q,Δs,v_χ,ω,s_max`.
2. Define `S_{N_h}(u)`.
3. Compute finite-resolution widths `ε_χ^{AB}, ε_τ^{AB}, ε_χ^C, ε_τ^C`.
4. Generate `A(χ,τ,η)`, `B(χ,τ,η)`, `C^{(0)}(χ,τ)`, `C^{(2)}(χ,τ)`, `C_read^{(0)}(χ,τ)`, and `C_read^{(2)}(χ,τ)`.
5. Compute initial observations `O_A^{(0)}`, `O_B^{(0)}` using `C_read^{(0)}`.
6. Compute identification-mode readouts `O_{A,m}^{(0)}`, `O_{B,m}^{(0)}`.
7. Apply observation perturbations and obtain `A'`, `B'`, `C'`.
8. Update `A'` and `B'` along `s` and test finite-cell entry.
9. If the cell is not reached within `s_max`, terminate as `collision_cell_not_reached`.
10. If the cell is reached, apply `C`.
11. Continue movement after collision in the reversed directions.
12. Compute post-collision observations `O_A^{(2)}`, `O_B^{(2)}` using `C_read^{(2)}`.
13. Compute post-collision identification readouts `O_{A,m}^{(2)}`, `O_{B,m}^{(2)}`.
14. Judge identification preservation, direction reversal, post-collision separation, representative-amplitude preservation, quasi-static behavior of C, and closure.

### 10.2 Pseudocode

```text
Input parameters:
  A_A, A_B, A_C
  Nh_chi_A, Nh_tau_A
  Nh_chi_B, Nh_tau_B
  Nh_chi_C, Nh_tau_C
  chi_A, chi_B, chi_C
  tau_A, tau_B, tau_0
  phi_A, phi_B, phi_C
  m_A, m_B
  q_A = +1, q_B = -1
  delta_s, v_chi
  omega_A, omega_B
  s_max

Minimal symmetric time phase:
  tau_A = -tau_0
  tau_B = -tau_0
  omega_A = omega_B

Observer readout centers:
  tau_C0 = -tau_0
  tau_C2 = +tau_0

Define S_Nh(u):
  K = (Nh + 1) / 2
  S = (1/K) * sum_{m=0}^{K-1} cos((2m+1)u)

Define identification mode:
  D_mP(eta) = exp(i * m_P * eta)
  D_read_m(eta) = exp(-i * m * eta)

Define cell widths:
  epsilon_chi_AB = pi / (min(Nh_chi_A, Nh_chi_B) + 1)
  epsilon_tau_AB = pi / (min(Nh_tau_A, Nh_tau_B) + 1)
  epsilon_chi_C = pi / (Nh_chi_C + 1)
  epsilon_tau_C = pi / (Nh_tau_C + 1)

Define local waves:
  A(chi,tau,eta) = A_A S_NhA_chi(chi-chi_A)
                   S_NhA_tau(tau-tau_A)
                   exp(i*m_A*eta) exp(i*phi_A)
  B(chi,tau,eta) = A_B S_NhB_chi(chi-chi_B)
                   S_NhB_tau(tau-tau_B)
                   exp(i*m_B*eta) exp(i*phi_B)
  Cread0(chi,tau) = A_C S_NhC_chi(chi-chi_C)
                    S_NhC_tau(tau-tau_C0) exp(-i*phi_C)
  Cread2(chi,tau) = A_C S_NhC_chi(chi-chi_C)
                    S_NhC_tau(tau-tau_C2) exp(-i*phi_C)

Initial observation:
  O_A0 = normalized integral of A * Cread0 over chi,tau
  O_B0 = normalized integral of B * Cread0 over chi,tau
  O_A_m0 = normalized integral of A * Cread0 * D_read_m over chi,tau,eta
  O_B_m0 = normalized integral of B * Cread0 * D_read_m over chi,tau,eta
  read m_A0, m_B0 from argmax over m

Observation update:
  chi_A' = chi_A + delta_chi_A_obs
  chi_B' = chi_B + delta_chi_B_obs
  tau_A' = tau_A + delta_tau_A_obs
  tau_B' = tau_B + delta_tau_B_obs
  q_A' = q_A
  q_B' = q_B
  m_A' = m_A
  m_B' = m_B
  A_A' = A_A
  A_B' = A_B
  C' approximately equals C0

Propagate toward collision cell:
  step = 0
  while abs(chi_A' - chi_B') >= epsilon_chi_AB
        or abs(tau_A' - tau_B') >= epsilon_tau_AB:
      if step >= s_max:
          return collision_cell_not_reached
      chi_A' = chi_A' + q_A' * v_chi * delta_s
      chi_B' = chi_B' + q_B' * v_chi * delta_s
      tau_A' = tau_A' + omega_A * delta_s
      tau_B' = tau_B' + omega_B * delta_s
      step = step + 1

Collision map:
  chi_A'' = chi_A'
  chi_B'' = chi_B'
  tau_A'' = tau_A'
  tau_B'' = tau_B'
  q_A'' = -q_A'
  q_B'' = -q_B'
  m_A'' = m_A'
  m_B'' = m_B'
  A_A'' = A_A'
  A_B'' = A_B'
  C'' approximately equals C'

Post-collision propagation:
  post_step = 0
  while abs(chi_A'' - chi_B'') <= epsilon_chi_AB
        or min(tau_A'', tau_B'') < tau_C2:
      if post_step >= s_max:
          return post_collision_propagation_not_completed
      chi_A'' = chi_A'' + q_A'' * v_chi * delta_s
      chi_B'' = chi_B'' + q_B'' * v_chi * delta_s
      tau_A'' = tau_A'' + omega_A * delta_s
      tau_B'' = tau_B'' + omega_B * delta_s
      post_step = post_step + 1

Post observation:
  O_A2 = normalized integral of A'' * Cread2 over chi,tau
  O_B2 = normalized integral of B'' * Cread2 over chi,tau
  O_A_m2 = normalized integral of A'' * Cread2 * D_read_m over chi,tau,eta
  O_B_m2 = normalized integral of B'' * Cread2 * D_read_m over chi,tau,eta

Judgement:
  identification modes preserved?
  identification modes not swapped?
  direction readouts reversed?
  A and B separated after collision?
  representative amplitudes preserved?
  observer C quasi-static?
  observer readout centers matched?
  closure maintained?
```

---

## 11. Required Outputs

### 11.1 Numerical Table

At minimum, output:

| Item | Initial | After observation | After collision | Post observation |
|---|---:|---:|---:|---:|
| `χ_A` | | | | |
| `χ_B` | | | | |
| `τ_A` | | | | |
| `τ_B` | | | | |
| `q_A` | | | | |
| `q_B` | | | | |
| `m_A` | | | | |
| `m_B` | | | | |
| `A_A` | | | | |
| `A_B` | | | | |
| `χ_C` | | | | |
| `τ_C` | | | | |

### 11.2 Figures

1. Central projection diagram of the initial configuration
2. Spatial localized waveforms of `A,B,C`
3. Temporal localized waveforms of `A,B,C`
4. Evolution of `χ_A,χ_B` before and after collision
5. Reversal of `q_A,q_B`
6. Before/after comparison of observation correlations `O_A,O_B`
7. Identification-mode readout before and after collision

### 11.3 Judgment Log

```text
label_mode_detected_A_initial: true/false
label_mode_detected_B_initial: true/false
label_mode_detected_A_final: true/false
label_mode_detected_B_final: true/false
label_mode_preserved_A: true/false
label_mode_preserved_B: true/false
label_mode_swapped: true/false
label_mode_lost: true/false
label_mode_cross_talk: true/false
q_reversed_A: true/false
q_reversed_B: true/false
amplitude_preserved_A: true/false
amplitude_preserved_B: true/false
observer_C_quasi_static: true/false
observer_time_centers_valid: true/false
separated_after_collision: true/false
collision_cell_reached: true/false
collision_cell_not_reached: true/false
post_collision_propagation_completed: true/false
closure_preserved: true/false
elastic_collision_map_valid: true/false
```

---

## 12. Core of This Specification

### 12.1 The Collision Instant Is Not Directly Observed

Collision is treated as a map before and after a finite-resolution cell, not as a point event.

```math
|\chi_A-\chi_B|<\epsilon_\chi^{AB}
```

```math
|\tau_A-\tau_B|<\epsilon_\tau^{AB}
```

### 12.2 Observer C Is Also Localized in the Temporal Direction

Observer `C` is not a continuous clock. For readout center `j`,

```math
C^{(j)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(j)})
e^{i\phi_C}.
```

For pre-collision and post-collision observations, use

```math
\tau_C^{(0)}=-\tau_0,\qquad \tau_C^{(2)}=+\tau_0.
```

Thus the temporal phase of the observation is determined by temporal odd-harmonic localization of `C` and the readout center `τ_C^{(j)}`.

For relative-phase readout, use the inverse-phase reference wave:

```math
C_{\rm read}^{(j)}(\chi,\tau)
=
A_C
S_{N_{h,\chi}^{C}}(\chi-\chi_C)
S_{N_{h,\tau}^{C}}(\tau-\tau_C^{(j)})
e^{-i\phi_C}.
```

Then

```math
A(\chi,\tau)C_{\rm read}^{(j)}(\chi,\tau)
\propto
e^{i(\phi_A-\phi_C)}.
```

### 12.3 Identification Oscillation Is Embedded in the Local Wave

The identification of `A` and `B` is not a variable name. It is an internal oscillation mode:

```math
D_{m_P}(\eta)=e^{i m_P\eta}.
```

As the local center `χ_P(s)` moves, this identification oscillation remains part of the same waveform.

The collision map changes `χ_P` and `q_P`, but does not change `D_{m_P}`.

Therefore, the test checks that the local wave carrying its identification oscillation is pushed out in the reversed direction.

### 12.4 Higher Odd Harmonics Do Not Increase the Conserved Quantity

Increasing the highest odd-harmonic order `N_h` is not an operation that externally increases the conserved quantity.

```math
S_{N_h}(0)=1
```

Thus

```math
\Psi_{N_h}(0)=A.
```

It refines phase structure with the same representative amplitude inside the closure condition of Axiom 1.

### 12.5 C Is Absorbed as a Heavy Curvature-Radius Generator

When

```math
A_C\gg A_A,A_B,
```

`C` can be treated not as a fully dynamical third body, but as a quasi-static background providing an effective curvature radius near AB.

```math
R_C^2=-\sum_C x_n^2
```

```math
\ell_{AB}\ll R_C
```

Under this condition, the AB interaction can be approximated as a local linear two-body problem.

---

## 13. Next Steps

Based on this specification, implement:

1. generator function for `S_{N_h}(u)`;
2. spatial, temporal, and identification-localized waves `A,B,C^{(0)},C^{(2)},C_read^{(0)},C_read^{(2)}`;
3. SVG output for the central projection diagram;
4. initial observation correlations;
5. identification-mode readout;
6. observation perturbations;
7. finite-resolution cell collision test;
8. complete elastic collision map;
9. continued motion after collision;
10. post-collision observation correlations;
11. judgment log;
12. PNG/SVG figure generation.
