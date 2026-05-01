# D1 v2.1: Particle Data Structure Specification of the Central Projection Framework — Frequency Interpretation and Dynamics Extension (Revised)

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd.
**ORCID**: 0009-0004-6753-4020
**Version**: v2.1
**Date**: May 2026
**DOI (v2.1)**: [10.5281/zenodo.19944488](https://doi.org/10.5281/zenodo.19944488)
**DOI (v2.0)**: [10.5281/zenodo.19941250](https://doi.org/10.5281/zenodo.19941250)
**Concept DOI**: [10.5281/zenodo.19941249](https://doi.org/10.5281/zenodo.19941249) (always points to latest)
**Predecessor**: v1.0 (April 2026, unpublished) / v2.0 (May 2026 published)

---

## §0 Position of This Paper

This paper is a design specification, not a mathematical proof. It explicitly defines, as a typed-integer-arithmetic specification, the particle data structures that have been implicitly assumed in the preceding paper series [W1]–[W11], [M1]–[M3], supplementary lectures [Sup1]–[Sup4], and the thought experiment [TE_Raxis].

### Revisions from v2.0 (v2.1)

**Removal of `δδt` field**: The `δδt : DeltaT` field added to the Particle structure in v2.0 (processing-period change, acceleration equivalent) is **removed** in v2.1. Reasons:

1. **Reversibility violation**: A `step(p)` containing `δt += δδt` does not return to its original state when δt's sign is reversed (violates the information conservation requirement of [W3]).
2. **Inconsistency with conservation laws**: Changing δt for a single particle without external interaction means the proper time changes ex nihilo, implicitly breaking energy and momentum conservation.
3. **GR free-particle principle**: A free particle follows a geodesic, and acceleration effects must arise from external forces (interactions). The basic operation for a single particle should not handle acceleration.

All acceleration effects are represented **within interaction processing (between particles or between particle and Mesh)** as changes in σ or δt that satisfy conservation laws. See the separate paper D2 ("Specification of Dynamics Rules of the Central Projection Framework," in preparation).

The corrected Particle structure has **4 fields** `(σ, V, δV, δt)` (§5). The corresponding `δδt_mesh` is also removed from the Mesh structure (§7).

### Major Revisions from v1.0 (introduced in v2.0, continued in v2.1)

1. **Reinterpret σ.k_xyzt as "vibrational frequencies of each axis"** rather than "wave numbers / acceleration" (§1.4)
2. **Add δV (phase displacement) to the Particle structure** (§5). ~v2.0 also added δδt, but this was removed in v2.1.~
3. **Introduce the Mesh structure to separate longitudinal modes (Higgs, gravitational waves) from Particles** (§7)
4. **Eliminate special status of t-axis, fully preserving 6-axis symmetry** (§3, §6)

The data type system (int256-based) and integer-arithmetic principles (no division, no square root, no trigonometric functions) remain unchanged from v1.0.

### Claims of This Paper

1. The boundedness of the universe and the discreteness at the Planck scale together permit all phase quantities to be represented as 256-bit signed integers (§1).
2. The $[L^{-1}]$ dimension of the xyzt axes is naturally interpreted as "vibrational frequency along each axis in phase space," fully preserving 6-axis symmetry (§3).
3. The complete state of a particle can be described by 4 fields: `SignVector` (identity label) + `OneDWave[κ]` (phase state) + `δV` (phase displacement) + `δt` (processing period) (§5).
4. Longitudinal modes such as the Higgs should be represented not as Particles but as vibrations of the `Mesh`, allowing unified treatment of mass generation, gravitational waves, and cosmological vibrations (§7).
5. All operations close under integer addition and multiplication, requiring no division, square root, or trigonometric functions (§6).
6. The basic operation `step(p)` (which only updates V) satisfies all five operational principles: time-reversal symmetry, conservation-law consistency, null transparency, type integrity, and INPUT/OUTPUT exchangeability (§6).

### Out of Scope

- Dynamical formulation of operation rules (physical meaning of phase evolution, concrete forms of interactions) — deferred to a separate paper (D2).
- Numerical derivation of specific physical quantities.
- Hardware implementation details.
- Mathematical proofs of theorems.

### Referenced Papers

- [W3] Structural Requirements for a Theory of Everything (DOI: 10.5281/zenodo.19601592)
- [W8] Combinatorial Structure of the 6-Dimensional Hyperrectangle (DOI: 10.5281/zenodo.19762134)
- [W9] sine-Gordon Equation — Foundation of Topological Solitons (DOI: 10.5281/zenodo.19650966)
- [W10] Four Modes of Shape-Invariant Waves (DOI: 10.5281/zenodo.19709798)
- [W11] Interactions of Shape-Invariant Waves (DOI: 10.5281/zenodo.19763463)
- [M1] Three-Layer Model of Central Projection: R–R₁–R₀ (DOI: 10.5281/zenodo.19691713)
- [TE_Raxis] Reexamination of the 6-Dimensional Sign Vector xyztRQ (DOI: 10.5281/zenodo.19904714)

---

## §1 int256 Boundedness and Dimensional Interpretation

### §1.1 Basis for Boundedness Requirement (unchanged from v1.0)

[W3] requires boundedness as a structural condition for any theory of everything. In the three-layer structure $R$–$R_1$–$R_0$ of [M1], the maximum curvature radius $R_0$ is finite. Combined with discreteness at the Planck scale, this yields an upper bound on the bits needed to represent all phase quantities.

### §1.2 Required Number of Bits

**Spatial axis**: From the observable universe diameter $D \approx 1.892 \times 10^{27}$ m and Planck length $l_P = 1.616 \times 10^{-35}$ m,

$$n_{\text{space}} = D / l_P \approx 1.171 \times 10^{62}, \quad b_{\text{space}} = 207 \text{ bits}$$

**Time axis**: From cosmic time $T \approx 8.710 \times 10^{17}$ s and Planck time $t_P = 5.391 \times 10^{-44}$ s,

$$n_{\text{time}} = T / t_P \approx 1.616 \times 10^{61}, \quad b_{\text{time}} = 204 \text{ bits}$$

**Conclusion**: 256-bit signed integers leave 14-digit margin in space and 15-digit margin in time.

### §1.3 Natural Units and Dimensions

Under $c = 1$, $G = 1$ (geometrized units; $\hbar$ is not adopted),

$$[L] = [T] = [M], \quad [L^{-1}] = [M^{-1}]$$

### §1.4 Frequency Interpretation (Core Revision in v2.0)

In v1.0, the $[L^{-1}]$ dimension of the xyzt axes was interpreted ambiguously as "wave number" or "acceleration." In v2.0, this is unified as **"vibrational frequency along each axis"** as the natural interpretation in phase space.

| Quantity | Dimension | v1.0 Interpretation | **v2.1 Interpretation** |
|---|---|---|---|
| σ.k_x, k_y, k_z | $[L^{-1}]$ | spatial wave number / acceleration | **vibrational frequency of each spatial axis** |
| σ.k_t | $[L^{-1}]$ | temporal frequency | **vibrational frequency of t-axis** (equal to others) |
| σ.R | $[L]=[M]$ | curvature radius | **scale quantity / corresponding period** (dual) |
| σ.Q | $[1]$ | internal symmetry | internal symmetry (unchanged) |
| δt | $[L]=[M]$ | proper time step | **processing step period** (not "time") |
| ~δδt~ | ~$[L]=[M]$~ | ~(not in v1.0)~ | ~added in v2.0, removed in v2.1 (violates conservation/reversibility)~ |

### §1.5 Advantages of the Frequency Interpretation

**(1) Complete preservation of 6-axis symmetry**: All four position-type axes (xyzt) carry independent vibrational modes; the t-axis is on equal footing with the others.

**(2) Elimination of continuum-mechanical concepts**: "Acceleration" is a continuum concept (second time derivative of position) that does not fit naturally into discrete integer arithmetic. "Vibrational frequency" is the natural discrete quantity in phase space.

**(3) Acceleration handled as interaction**: Acceleration-like effects (gravitational fields, inertial frame influences) are not represented for a single particle but **entirely as interaction processing** (between particles or between particle and Mesh), where σ or δt change in a conservation-consistent manner. This is consistent with the GR free-particle principle (free particles follow geodesics). See the separate paper D2.

**(4) Frequency × period = phase**: $\sigma.k \cdot \delta t$ becomes a dimensionless phase increment ($[L^{-1}] \cdot [L] = [1]$), making the phase evolution rule type-consistent.

**(5) Automatic inclusion of special relativity**: The dispersion relation $\sigma.k_t^2 = \sigma.k_x^2 + \sigma.k_y^2 + \sigma.k_z^2 + m^2$ can be naturally written as integer relations among frequencies.

---

## §2 Definition of Meta Types

### §2.1 Scalar Types (unchanged from v1.0)

```
Type         Dimension  Range                     Use
L256         [L]=[M]    |v| ≤ 2^255 - 1          length / time / mass / period
Linv256      [L⁻¹]      |v| ≤ 2^255 - 1          frequency / inverse scale
Scalar256    [1]        |v| ≤ 2^255 - 1          dimensionless
FullPhase    [1]        0 ≤ v ≤ P                0–720°
CyclePhase   [1]        0 ≤ v ≤ P/2              0–360° (position phase)
DeltaPhase   [1]        -P/2 ≤ v ≤ P/2           ±360° (phase displacement)
SpinPhase    [1]        v ∈ {0, P/4, P/2, P}     {0, 180, 360, 720}°
```

Phase constant $P = \lfloor (2^{255}-1)/720 \rfloor \times 720$.

### §2.2 Axis Component Types (renamed and reinterpreted in v2.0)

| Type | Base | Dimension | Physical Meaning (v2.0) |
|---|---|---|---|
| `SpatialFreq` <: Linv256 | int256 | $[L^{-1}]$ | vibrational frequency of x, y, z axes |
| `TemporalFreq` <: Linv256 | int256 | $[L^{-1}]$ | vibrational frequency of t-axis (equal to others) |
| `CurvatureR` <: L256 | int256 | $[L]$ | scale / period (signed square number) |
| `ColorQ` <: Scalar256 | int256 | $[1]$ | internal symmetry label (finite group) |
| `DeltaT` <: L256 | int256 | $[L]$ | processing step period |

**Renaming from v1.0**: `SpatialK` → `SpatialFreq`, `TemporalK` → `TemporalFreq`. The type system itself is unchanged; only the interpretation is clarified.

### §2.3 Operation Rules (unchanged from v1.0)

**Addition**: Only between same types (CyclePhase + DeltaPhase = CyclePhase, etc.).

**Multiplication**:
- `L256 × L256 → L2_256` $[L^2]$
- `Linv256 × Linv256 → Linv2_256` $[L^{-2}]$
- `Linv256 × L256 → Scalar256` $[1]$ (dimensionless)

**Division**: Not used.

---

## §3 Definition of SignVector (Frequency Interpretation Version)

### §3.1 Structure

```
SignVector = (
    x : SpatialFreq,    -- frequency of x-axis vibration [L⁻¹]
    y : SpatialFreq,    -- frequency of y-axis vibration
    z : SpatialFreq,    -- frequency of z-axis vibration
    t : TemporalFreq,   -- frequency of t-axis vibration [L⁻¹] (equal to others)
    R : CurvatureR,     -- scale quantity / period [L]
    Q : ColorQ,         -- internal symmetry label [1]
)
```

### §3.2 Physical Meaning

Each σ component represents either a "vibrational frequency along the corresponding axis" or a "scale quantity":

- σ.k_axis = 0: the axis does not vibrate (no structure in that direction)
- σ.k_axis ≠ 0: the axis vibrates at frequency |k_axis|
- sign(k_axis): direction of phase rotation (+ for clockwise, − for counterclockwise)
- |k_axis|: magnitude of vibrational frequency (discrete integer value)

### §3.3 Invariants

The values of SignVector change in interaction events, but specific combinations define stable modes:

- $\kappa_s = |\{i \in \{x,y,z\} : k_i \neq 0\}|$ — number of non-zero spatial axes
- $\kappa_t = (1 \text{ if } k_t \neq 0 \text{ else } 0)$ — non-zero indicator for t-axis
- $\kappa = \kappa_s + \kappa_t \in \{0, 1, 2\}$ — [W10] Proposition 2.1: $\kappa \geq 3$ is unstable

### §3.4 Dispersion Relation

The mass-shell condition is directly expressed as an integer relation among frequencies:

$$\sigma.k_t^2 = \sigma.k_x^2 + \sigma.k_y^2 + \sigma.k_z^2 + m^2$$

- $m = 0$ (photon, gluon): $k_t^2 = k_x^2 + k_y^2 + k_z^2$
- $m \neq 0$ (massive particle): the square root of the difference corresponds to mass

This automatically incorporates special relativity into the framework.

### §3.5 Correspondence with [W8] Table A (continued from v1.0)

The sign patterns $\{+, -, 0\}^6$ in [W8] Table A correspond to the signs of SignVector. Integer values (axis scale values) correspond to $v_i$ in [Mass Structure] §1.

Representative particles:

| Particle | Set $(x,y,z,t,R,Q)$ | $\kappa$ | Classification |
|---|---|---|---|
| Photon $\gamma$ | $(0,0,0,0,+,0)$ | 0 | gauge boson (free, spherical wave) |
| Gluon $g$ | $(0,0,0,0,0,Q_i \to Q_j)$ | 0 | gauge boson |
| Higgs $H$ | $(0,0,+,0,0,0)$ | 1 | scalar (longitudinal interpretation in §7) |
| $W^+$ | $(+,0,0,0,+,0)$ | 1 | gauge boson |
| $Z^0$ | $(0,+,0,0,+,0)$ | 1 | gauge boson |
| Graviton | $(0,0,0,\pm,\pm,0)$ | 1 | tensor boson |
| Electron $e^-$ | $(+,0,0,+,0,0)$ | 2 | fermion |

---

## §4 Definition of OneDWave

### §4.1 Structure (unchanged from v1.0)

```
OneDWave = (
    θ₀ : CyclePhase,   -- position phase [1]   (S1 §2.1, W9 §8.2)
    A  : Scalar256,     -- amplitude [1]        (S1 §2.1, W9 §8.2)
    ℓ  : DeltaPhase,    -- spread phase [1]     (W9 §8.3)
)
```

### §4.2 Physical Meaning

Each particle has one OneDWave per non-zero position-type axis (at most $\kappa$). Each OneDWave holds the current phase state of vibration in that axis direction.

- $\theta_0$: current phase of vibration (free parameter, evolves over time)
- $A$: amplitude (free parameter)
- $\ell$: spatial extent of localization (dependent on SignVector, [W9] §8.3)

---

## §5 Definition of Particle (v2.1: 4-field structure)

### §5.1 Structure

```
Particle = (
    σ   : SignVector,    -- identity (vibrational frequencies of each axis)
    V   : OneDWave[κ],   -- current phase state of each non-zero position-type axis
    δV  : OneDWave[κ],   -- phase displacement per processing step (added in v2.0)
    δt  : DeltaT,         -- current processing step period
)
```

**Change from v2.0**: The `δδt : DeltaT` field is removed (5 fields → 4 fields). See §0 for reasons.

### §5.2 Physical Meaning of Each Field

#### `σ : SignVector`

The particle's identity. Specifies frequencies of vibration in each axis. Values change in interaction events (formulated as annihilation + creation, but as a data structure σ is updated).

#### `V : OneDWave[κ]`

Current phase state of each non-zero position-type axis. Array length $\kappa$ equals the number of non-zero position-type axes in SignVector ($\leq 2$). Each V[j] holds the (position phase, amplitude, spread) of vibration in the corresponding axis.

Index correspondence:

```
σ.k_x ≠ 0  ⇒  V contains wave_x
σ.k_y ≠ 0  ⇒  V contains wave_y
σ.k_z ≠ 0  ⇒  V contains wave_z
σ.k_t ≠ 0  ⇒  V contains wave_t
```

#### `δV : OneDWave[κ]` (added in v2.0)

A difference structure isomorphic to V. Holds how V changes in one processing step:

```
δV[j] = (
    Δθ₀ : DeltaPhase,    -- phase displacement
    ΔA  : Scalar256,      -- amplitude displacement
    Δℓ  : DeltaPhase,     -- spread displacement
)
```

For free particles, δV[j].Δθ₀ can be computed from σ:

$$\delta V[j].\Delta\theta_0 = \sigma.k_{a(j)} \cdot \delta t$$

However, in interaction events δV may change independently.

#### `δt : DeltaT`

**Processing step period** (NOT time). A universal computation counter.

- $\delta t > 0$: forward operation (causal forward direction)
- $\delta t < 0$: reverse operation (guarantees reversibility)
- $\delta t = 0$: at rest

What is observed as "time" is the accumulation of V[t-axis].θ₀, not δt itself (see §6.4).

#### ~`δδt : DeltaT` (removed in v2.1)~

The Particle structure included a `δδt` field in v2.0, but this is removed in v2.1. Reasons:

- **Reversibility violation**: Executing `δt += δδt` in `step(p)` does not return to the original state when δt's sign is reversed
- **Conservation-law violation**: Single-particle δt change implies energy/momentum change ex nihilo, inconsistent with conservation laws
- **GR free-particle principle**: Acceleration effects must arise from external forces (interactions); the basic operation for a free particle should not handle them

Acceleration effects are represented in the interaction processing (between particles or between particle and Mesh) of §3 onward, where σ or δt change in a conservation-consistent manner. See the separate paper D2.

### §5.3 Free Particle Dynamics Close with δV

Although the second-order closure of physical laws (Newtonian mechanics, GR, QM) initially suggested that δδt was needed, in v2.1:

- Free particles undergo only uniform geodesic motion (GR principle)
- Acceleration is the result of interaction
- The basic operation for a single particle is **only the update of phase velocity (δV) and current phase (V)**

Therefore, **4 fields (σ, V, δV, δt) form a complete representation for free particles.**

### §5.4 Memory Footprint (v2.1)

```
σ      : 6 × 256 = 1,536 bit
V      : κ × 3 × 256 ≤ 1,536 bit (κ ≤ 2)
δV     : κ × 3 × 256 ≤ 1,536 bit
δt     : 256 bit
Total  ≤ 4,864 bit = 608 byte / particle
```

A 32-byte reduction from v2.0's 640 byte (removal of δδt field). Still extremely small compared to v1.0's 224–416 byte (which lacked δV).

---

## §6 Operation Rules

### §6.1 Basic: 1-Step Update (v2.1)

The basic operation `step(p)` consists of only one operation:

```
function step(p : Particle | null) → Particle | null:
    if p == null:
        return null
    for j in range(p.κ):
        p.V[j].θ₀ += p.σ.k_{a(j)} · p.δt
    return p
```

For each non-zero position-type axis $j$:

$$V[j].\theta_0 \mathrel{+}= \sigma.k_{a(j)} \cdot \delta t$$

That is all. `σ`, `δV`, and `δt` are not updated by `step(p)`.

**Change from v2.0**: v2.0 included `δt += δδt`, but this is removed in v2.1 (see §5.2).

### §6.2 Type Consistency Check

```
σ.k_axis · δt : SpatialFreq × DeltaT 
              = Linv256 × L256 
              = Scalar256 [1]                ✓

V[j].θ₀ + (σ.k · δt) : CyclePhase + Scalar256 
                     = CyclePhase [1] (mod P/2)  ✓
```

### §6.2.1 Compliance with the Five Operational Principles (v2.1)

`step(p)` satisfies all five basic operational principles:

1. **Null transparency**: explicitly handled at the function's start ✓
2. **In-place update**: V[j].θ₀ updated directly, no new Particle generated ✓
3. **Type integrity**: complies with D1 §2.3 rules (see §6.2) ✓
4. **INPUT/OUTPUT exchangeability**: function has symmetric read/write structure ✓
5. **Time-reversal reversibility**: complete reversal by sign-flipping δt:

```
First call (δt = +δt₀):  V → V + σ.k · δt₀
Sign-flip δt (δt → -δt₀)
Second call (δt = -δt₀): V → V + σ.k · δt₀ + σ.k · (-δt₀) = V  ✓
```

This structurally guarantees the information-conservation requirement of [W3].

### §6.3 Polarization and Multi-Axis Phase Difference

For particles with multiple non-zero axes (e.g., a polarized directional photon):

```
σ_polarized_light = (+k_x, +k_y, 0, 0, +R_γ, 0)
V = [wave_x, wave_y]
```

Polarization phase difference:

$$\Delta\varphi = V[0].\theta_0 - V[1].\theta_0 \in \text{CyclePhase}$$

| Polarization State | $\Delta\varphi$ |
|---|---|
| Linear (in-phase) | $0$ |
| Linear (anti-phase) | $P/4 = 180°$ |
| Circular (left) | $+P/8 = +90°$ |
| Circular (right) | $-P/8 = -90°$ |
| Elliptical | any integer value |

All exactly representable as int256 integers (integer multiples of $15°$ correspond to integer multiples of $P/48$ without error).

### §6.4 Treatment of Time: δt is a Processing Step, Time is V[t].θ₀

**Important principle**: δt is not time but a universal processing step. What is observed as time is the accumulation of V[t-axis].θ₀ for each particle.

Electron example (σ.k_t ≠ 0):

```
σ_electron = (+k_x, 0, 0, +k_t, 0, 0)
V = [wave_x, wave_t]
1-step:
  V[0].θ₀ += σ.k_x · δt   -- x-axis phase advances (momentum equivalent)
  V[1].θ₀ += σ.k_t · δt   -- t-axis phase advances ("time experienced by electron")
```

Photon example (σ.k_t = 0):

```
σ_photon = (+k_x, +k_y, 0, 0, +R_γ, 0)
V = [wave_x, wave_y]      -- no t-axis component
-- because σ.k_t = 0, no t-axis OneDWave exists in V
-- consequence: zero proper time of photon is structurally guaranteed
```

### §6.5 Acceleration Effects Are Handled in Interaction Processing (v2.1)

In v2.0 we stated "δδt is determined by the background R(r) or by interactions," but in v2.1:

- **Acceleration effects cannot be handled in the basic single-particle operation** (violates conservation laws and reversibility, see §5.2)
- **All acceleration effects are realized in interaction processing** (where σ or δt are exchanged between particles or with the Mesh)

The concrete interaction rules are defined in the separate paper D2. Gravitational time dilation, polarization rotation, and frequency shifts are all expressed as the cumulative result of interaction events that satisfy conservation laws.

### §6.6 Operations Not Required (same as v1.0)

The following operations are not required in this specification:

1. **Division**: $R/R_1$ is converted to $R \cdot R_1^{-1}$ as multiplication
2. **Square root**: a continuum concept; in discrete systems, integer values directly determine physical quantities
3. **Trigonometric functions**: tools of continuum sine-Gordon; replaceable by discrete modulation tables
4. **Floating-point arithmetic**: all fields are int256, no rounding error structurally arises

Where cos/sin would be needed (e.g., longitudinal mode modulation), discrete modulation tables (lookup tables) of finite resolution are used.

---

## §7 Introduction of the Mesh Structure (added in v2.0)

### §7.1 Separation of Longitudinal Modes

In the interpretation of [W8] Table A, the Higgs ($\sigma_H = (0, 0, +k_z, 0, 0, 0)$) differs from other spin-1 bosons in having $\sigma.R = 0$ and $\text{spin} = 0$. Structurally it is natural to treat the Higgs not as a "transverse particle vibrating position phase" but as a **"longitudinal vibration of the mesh itself."**

Similarly, the graviton (gravitational waves) can be interpreted as longitudinal-like vibration of the spacetime metric — a vibrational mode of the Mesh.

Treating these longitudinal modes as Particles strains the data structure (OneDWave is essentially a type for phase rotation and cannot directly represent axial expansion-contraction). This paper introduces a **Mesh structure to manage longitudinal modes separately from Particles**.

### §7.2 Structure of Mesh (v2.1)

```
Mesh = (
    -- background scale
    R_0           : L256,             -- central length of central projection (background)
    
    -- single-axis longitudinal modes (spin 0, SO(3) symmetry breaking)
    H_x           : OneDWave,         -- x-axis longitudinal
    H_y           : OneDWave,         -- y-axis longitudinal
    H_z           : OneDWave,         -- z-axis longitudinal = Higgs
    H_t           : OneDWave,         -- t-axis longitudinal
    
    -- multi-axis longitudinal modes (spin 2, gravitational wave polarizations)
    H_xy_anti     : OneDWave,         -- xy anti-phase = gravitational wave + polarization
    H_xy_quarter  : OneDWave,         -- xy 90° phase = gravitational wave × polarization
    H_yz_anti     : OneDWave,         -- yz anti-phase
    H_xz_anti     : OneDWave,         -- xz anti-phase
    
    -- displacements of longitudinal modes
    δH_x, δH_y, δH_z, δH_t            : OneDWave,
    δH_xy_anti, δH_xy_quarter         : OneDWave,
    δH_yz_anti, δH_xz_anti            : OneDWave,
    
    -- global processing step
    δt_mesh       : DeltaT,
)
```

**Changes from v2.0**:

- Removed `H_xyz_inphase` (xyz in-phase) and replaced with axis-specific single-axis modes `H_x, H_y, H_z, H_t` (the xyz-isotropic interpretation cannot explain the Higgs's generational mass hierarchy; see [W11] §8.5)
- The Higgs is positioned as `H_z` (z-axis longitudinal mode)
- Removed `δδt_mesh` (same reason as Particle: reversibility and conservation)

### §7.3 Physical Identification of Longitudinal Modes

| Mode | Vibration Pattern | Physical Identification | Spin |
|---|---|---|---|
| H_z | z-axis longitudinal (SO(3) breaking) | **Higgs** | 0 |
| H_x, H_y, H_t | each-axis longitudinal | physical identification undetermined (candidates) | 0 |
| H_xy_anti | x stretches while y compresses (anti-phase) | gravitational wave + polarization | 2 |
| H_xy_quarter | xy 90° phase rotation | gravitational wave × polarization | 2 |
| H_yz_anti | yz anti-phase | gravitational wave (other polarization) | 2 |
| H_xz_anti | xz anti-phase | gravitational wave (other polarization) | 2 |

**Higgs = z-axis longitudinal mode (H_z)**: The z-axis selection in W8 Table A's $\sigma_H = (0, 0, +k_z, 0, 0, 0)$ is preserved as the direction of SO(3) symmetry breaking. Third-generation fermions sharing the z-axis have the strongest coupling and acquire the largest mass (consistent with [W11] §8).

**Two polarizations of gravitational waves (+, ×)**: expressed as multi-axis longitudinal modes (H_xy_anti, H_xy_quarter, etc.). SO(3) symmetry is preserved.

### §7.4 Coupling between Mesh and Particle (v2.1)

The vibrational frequency σ.k of a particle is effectively modulated by the state of the Mesh. Mass generation and frequency shifts are expressed as interactions through this modulation (details in the separate paper D2).

Effective frequency (mass generation example):

$$\sigma.k_{\text{effective}}(particle, mesh) = \sigma.k(particle) - \delta k_H(\text{Mesh.H\_z VEV})$$

Third-generation fermions sharing the z-axis couple most strongly with the H_z VEV and acquire the largest mass (consistent with [W11] §8).

### §7.5 Dynamics (v2.1: basic operation is step only)

Mesh update:

```
For each longitudinal mode H_*:
    H_*.θ₀ += (intrinsic frequency) · δt_mesh
```

Particle update (independent of Mesh, basic operation is `step(p)` only):

```
For each particle p:
    step(p)   -- as in §6.1, V only is updated
```

The Mesh-Particle coupling (mass generation, gravitational waves, etc.) is not realized within single-particle `step(p)`, but in interaction processing (defined in D2) where σ or δt are changed.

### §7.6 Memory Footprint (v2.1)

```
R_0    : 256 bit
H_*    : 8 × 3 × 256 = 6,144 bit  -- 4 single-axis + 4 multi-axis = 8 modes
δH_*   : 8 × 3 × 256 = 6,144 bit
δt_mesh : 256 bit
Total  = 12,800 bit = 1,600 byte / entire universe
```

About 1.6 KB for the entire universe. The total memory for $N$ particles is $608N + 1600$ byte (a 32-byte reduction per particle from v2.0).

---

## §8 Verification of W3 Structural Requirements (v2.1)

| W3 Requirement | Satisfied | Basis |
|---|---|---|
| Boundedness | ✓ | int256 represents all phase quantities (§1) |
| Discreteness | ✓ | all fields are integer types (§2–§5) |
| **Information conservation (time-reversal reversibility)** | **✓ (rigorous in v2.1)** | step(p) is fully reversed by sign-flipping δt (§6.2.1) |
| Avoidance of zero division | ✓ | division not used (§6.6) |
| Closure under integer arithmetic | ✓ | only addition/multiplication; cos via discrete modulation table (§7.4) |
| **6-axis symmetry** | **✓ (v2.0+)** | frequency interpretation makes all axes equivalent (§1.4, §3) |
| **Transverse-longitudinal separation** | **✓ (v2.0+)** | Particle vs Mesh (§7) |
| **Structural conservation laws** | **✓ (rigorous in v2.1)** | step(p) does not change σ or δt; conservation is automatic |
| Non-privileging of time | ✓ | t equal to other axes; time is V[t].θ₀ (§6.4) |
| **All five operational principles satisfied** | **✓ (v2.1)** | null transparency, in-place, type integrity, exchangeability, reversibility (§6.2.1) |
| Finite representability | ✓ | ≤ 608 byte per particle, +1.6 KB for whole universe (§5.4, §7.6) |

---

## §9 Incompatibilities with v1.0

### §9.1 Changes in Interpretation

- σ.k unified from "wave number / acceleration" to "vibrational frequency"
- δt clarified from "proper time step" to "processing step period"
- "Time" is a derived quantity appearing as V[t-axis].θ₀

### §9.2 Structural Extension and Refinement

**Added in v2.0** (continued in v2.1):

- δV added to Particle (v1.0's Particle is included as a subset)
- Mesh structure newly introduced (not present in v1.0)

**Added in v2.0 and removed in v2.1**:

- Particle's `δδt` field (removed for conservation/reversibility violations)
- Mesh's `δδt_mesh` field (same reason)
- Mesh's `H_xyz_inphase` mode (cannot explain observed generational mass hierarchy; replaced with axis-specific modes such as `H_z`)

### §9.3 Compatibility with [W8] Table A

The 62-particle classification of [W8] Table A remains valid in v2.1. The Higgs is positioned as an **excitation of the Mesh's H_z mode** (the H_xyz_inphase proposed in v2.0 has been corrected in v2.1, consistent with [W11] §8.5). The total count (61 SM states + 1 graviton = 62) is unchanged.

---

## §10 Conclusion

Essence of v2.1:

- **xyzt are frequencies** (natural interpretation in phase space)
- **6 axes are completely symmetric** (4 position-type + R + Q)
- **Transverse (Particle) and longitudinal (Mesh) are separated**
- **Basic operation step(p) updates only V** (satisfies all 5 principles; acceleration effects are separated into interaction processing)
- **Time is V[t].θ₀** (a position phase, advancing per particle)
- **All operations close under int256 integer arithmetic**

With this structure, the dynamics rules (D2, separate paper) can describe:

- Gravitational time dilation, polarization rotation (particle-Mesh interaction)
- Higgs mechanism (Mesh.H_z VEV coupling to particle's σ.k_z)
- Gravitational wave emission (excitation of Mesh.H_xy_anti mode)
- Coulomb 2-body problem (inter-particle interaction)

all within the same framework. All acceleration effects are realized within interaction processing in a manner that respects conservation laws and reversibility.

This paper is a design specification; these dynamical applications are deferred to separate papers.

---

## References

- [W3] Kihara, N. "Structural Requirements for a Theory of Everything." DOI: 10.5281/zenodo.19601592 (2026).
- [W8] Kihara, N. "Combinatorial Structure of the 6-Dimensional Hyperrectangle." DOI: 10.5281/zenodo.19762134 (2026).
- [W9] Kihara, N. "sine-Gordon Equation — Foundation of Topological Solitons." DOI: 10.5281/zenodo.19650966 (2026).
- [W10] Kihara, N. "Four Modes of Shape-Invariant Waves." DOI: 10.5281/zenodo.19709798 (2026).
- [W11] Kihara, N. "Interactions of Shape-Invariant Waves." DOI: 10.5281/zenodo.19763463 (2026).
- [M1] Kihara, N. "Three-Layer Model of Central Projection: R–R₁–R₀." DOI: 10.5281/zenodo.19691713 (2026).
- [TE_Raxis] Kihara, N. "Reexamination of the 6-Dimensional Sign Vector xyztRQ." DOI: 10.5281/zenodo.19904714 (2026).

---

**Author**: Noriaki Kihara
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**License**: CC BY 4.0
