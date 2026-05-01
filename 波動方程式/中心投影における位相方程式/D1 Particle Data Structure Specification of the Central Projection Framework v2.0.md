# D1 v2.0: Particle Data Structure Specification of the Central Projection Framework — Frequency Interpretation and Dynamics Extension

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd.
**ORCID**: 0009-0004-6753-4020
**Version**: v2.0
**Date**: May 2026
**DOI**: [10.5281/zenodo.19941250](https://doi.org/10.5281/zenodo.19941250)
**Predecessor**: v1.0 (April 2026, unpublished)

---

## §0 Position of This Paper

This paper is a design specification, not a mathematical proof. It explicitly defines, as a typed-integer-arithmetic specification, the particle data structures that have been implicitly assumed in the preceding paper series [W1]–[W11], [M1]–[M3], supplementary lectures [Sup1]–[Sup4], and the thought experiment [TE_Raxis].

### Major Revisions from v1.0

1. **Reinterpret σ.k_xyzt as "vibrational frequencies of each axis"** rather than "wave numbers / acceleration" (§1.4)
2. **Add δV (phase displacement) and δδt (processing-period change) to the Particle structure** (§5)
3. **Introduce the Mesh structure to separate longitudinal modes (Higgs, gravitational waves) from Particles** (§7)
4. **Eliminate special status of t-axis, fully preserving 6-axis symmetry** (§3, §6)

The data type system (int256-based) and integer-arithmetic principles (no division, no square root, no trigonometric functions) remain unchanged from v1.0.

### Claims of This Paper

1. The boundedness of the universe and the discreteness at the Planck scale together permit all phase quantities to be represented as 256-bit signed integers (§1).
2. The $[L^{-1}]$ dimension of the xyzt axes is naturally interpreted as "vibrational frequency along each axis in phase space," fully preserving 6-axis symmetry (§3).
3. The complete state of a particle can be described by 5 fields: `SignVector` (identity label) + `OneDWave[κ]` (phase state) + `δV` (phase displacement) + `δt, δδt` (processing period and acceleration) (§5).
4. Longitudinal modes such as the Higgs should be represented not as Particles but as vibrations of the `Mesh`, allowing unified treatment of mass generation, gravitational waves, and cosmological vibrations (§7).
5. All operations close under integer addition and multiplication, requiring no division, square root, or trigonometric functions (§6).

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

| Quantity | Dimension | v1.0 Interpretation | **v2.0 Interpretation** |
|---|---|---|---|
| σ.k_x, k_y, k_z | $[L^{-1}]$ | spatial wave number / acceleration | **vibrational frequency of each spatial axis** |
| σ.k_t | $[L^{-1}]$ | temporal frequency | **vibrational frequency of t-axis** (equal to others) |
| σ.R | $[L]=[M]$ | curvature radius | **scale quantity / corresponding period** (dual) |
| σ.Q | $[1]$ | internal symmetry | internal symmetry (unchanged) |
| δt | $[L]=[M]$ | proper time step | **processing step period** (not "time") |
| δδt | $[L]=[M]$ | (not in v1.0) | **change in processing period** (acceleration effect) |

### §1.5 Advantages of the Frequency Interpretation

**(1) Complete preservation of 6-axis symmetry**: All four position-type axes (xyzt) carry independent vibrational modes; the t-axis is on equal footing with the others.

**(2) Elimination of continuum-mechanical concepts**: "Acceleration" is a continuum concept (second time derivative of position) that does not fit naturally into discrete integer arithmetic. "Vibrational frequency" is the natural discrete quantity in phase space.

**(3) Clear separation of acceleration effects**: Acceleration-like effects (gravitational fields, inertial frame influences) are represented by the change δδt of the processing period, structurally separated from the particle's intrinsic frequency σ.k.

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

## §5 Definition of Particle (δV and δδt added in v2.0)

### §5.1 Structure

```
Particle = (
    σ   : SignVector,    -- identity (vibrational frequencies of each axis)
    V   : OneDWave[κ],   -- current phase state of each non-zero position-type axis
    δV  : OneDWave[κ],   -- phase displacement per processing step (added in v2.0)
    δt  : DeltaT,         -- current processing step period
    δδt : DeltaT,         -- change in processing period (added in v2.0, ≈ acceleration)
)
```

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

#### `δδt : DeltaT` (added in v2.0)

**Change in processing period** (acceleration equivalent). Determined by background R(r) or interactions:

```
1-step update: δt += δδt
```

Physical interpretation:

- $\delta\delta t = 0$: free particle (constant processing period)
- $\delta\delta t > 0$: processing period elongates (clock slows down in gravity / acceleration)
- $\delta\delta t < 0$: processing period shortens (rising gravitational potential, etc.)

### §5.3 Hierarchy Closes at δδt

Because the basic laws of physics close at second order in time (Newtonian mechanics, GR, quantum mechanics), no higher-order structure (δδδt, etc.) is needed. **The 5 fields (σ, V, δV, δt, δδt) form a complete representation for dynamics.**

### §5.4 Memory Footprint

```
σ      : 6 × 256 = 1,536 bit
V      : κ × 3 × 256 ≤ 1,536 bit (κ ≤ 2)
δV     : κ × 3 × 256 ≤ 1,536 bit
δt     : 256 bit
δδt    : 256 bit
Total  ≤ 5,120 bit = 640 byte / particle
```

This is increased from v1.0's 224–416 byte, but still extremely small for a complete description of any particle in the universe.

---

## §6 Operation Rules

### §6.1 Basic: 1-Step Update

For each non-zero position-type axis $j$:

$$\delta V[j].\Delta\theta_0 \leftarrow \sigma.k_{a(j)} \cdot \delta t$$

$$V[j].\theta_0 \mathrel{+}= \delta V[j].\Delta\theta_0$$

Then update the processing period itself:

$$\delta t \mathrel{+}= \delta\delta t$$

### §6.2 Type Consistency Check

```
σ.k_axis · δt : SpatialFreq × DeltaT 
              = Linv256 × L256 
              = Scalar256 [1]                ✓

V[j].θ₀ + δV[j].Δθ₀ : CyclePhase + DeltaPhase 
                    = CyclePhase [1] (mod P/2)  ✓

δt + δδt : DeltaT + DeltaT 
         = DeltaT [L]                       ✓
```

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
  δV[0].Δθ₀ = σ.k_x · δt    -- x-axis phase displacement (momentum equivalent)
  δV[1].Δθ₀ = σ.k_t · δt    -- t-axis phase displacement (electron's proper time progress)
  V[0].θ₀ += δV[0].Δθ₀
  V[1].θ₀ += δV[1].Δθ₀       -- "time experienced by electron"
  δt += δδt                    -- update processing period (δδt ≠ 0 in gravity)
```

Photon example (σ.k_t = 0):

```
σ_photon = (+k_x, +k_y, 0, 0, +R_γ, 0)
V = [wave_x, wave_y]      -- no t-axis component
δV = [δwave_x, δwave_y]
-- because σ.k_t = 0, no t-axis OneDWave exists in V
-- consequence: zero proper time of photon is structurally guaranteed
```

### §6.5 Acceleration Effect (Background Coupling)

δδt is determined by the background R(r) or by interactions. The concrete decision rule is defined in a separate paper (D2).

Qualitatively:

- In a gravitational field: $\delta\delta t = f(R_{\text{background}}, \sigma)$ — background curvature modulates processing period
- Result: σ.k · δt changes → frequency shift / time delay / polarization phase rotation

All of these are expressed as phase evolution by integer arithmetic.

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

### §7.2 Structure of Mesh

```
Mesh = (
    -- background scale
    R_0           : L256,             -- central length of central projection (background)
    
    -- longitudinal modes
    H_xyz_inphase : OneDWave,         -- xyz in-phase expansion-contraction = Higgs
    H_xy_anti     : OneDWave,         -- xy anti-phase = gravitational wave + polarization
    H_xy_quarter  : OneDWave,         -- xy 90° phase = gravitational wave × polarization
    H_yz_anti     : OneDWave,         -- yz anti-phase = other gravitational wave polarizations
    H_xz_anti     : OneDWave,         -- xz anti-phase = other gravitational wave polarizations
    
    -- displacements of longitudinal modes
    δH_xyz_inphase : OneDWave,
    δH_xy_anti     : OneDWave,
    δH_xy_quarter  : OneDWave,
    δH_yz_anti     : OneDWave,
    δH_xz_anti     : OneDWave,
    
    -- global processing step
    δt_mesh       : DeltaT,
    δδt_mesh      : DeltaT,
)
```

### §7.3 Physical Identification of Longitudinal Modes

| Mode | Vibration Pattern | Physical Identification | Spin |
|---|---|---|---|
| H_xyz_inphase | xyz axes in-phase expansion-contraction | Higgs (scalar) | 0 |
| H_xy_anti | x stretches while y compresses (anti-phase) | gravitational wave + polarization | 2 |
| H_xy_quarter | xy 90° phase rotation | gravitational wave × polarization | 2 |
| H_yz_anti | yz anti-phase | gravitational wave (other polarization) | 2 |
| H_xz_anti | xz anti-phase | gravitational wave (other polarization) | 2 |

**Higgs = spherically symmetric in-phase longitudinal vibration of xyz axes**: the choice of z-axis in W8 is a representative notation; the essence is a 3-axis symmetric longitudinal vibration.

**Two polarizations of gravitational waves (+, ×)**: correspond to H_xy_anti and H_xy_quarter (and other axis selections give other polarizations).

### §7.4 Coupling between Mesh and Particle

The vibrational frequency σ.k of a particle is effectively modulated by the state of the Mesh. Mass generation and frequency shifts are expressed as this modulation.

Effective frequency:

$$\sigma.k_{\text{effective}}(particle, mesh) = \sigma.k(particle) \cdot \text{mod}(mesh)$$

$$\text{mod}(mesh) = 1 + \alpha_H \cdot \text{disc-cos}(H_{\text{xyz\_inphase}}.\theta_0) + \cdots$$

Here `disc-cos` is a finite-resolution discrete modulation table (e.g., 12 stages at 30° intervals) substituting for the continuous trigonometric function.

### §7.5 Dynamics

Mesh update:

```
For each longitudinal mode H_*:
    δH_*.Δθ₀ ← (intrinsic frequency) · δt_mesh
    H_*.θ₀  += δH_*.Δθ₀
δt_mesh += δδt_mesh
```

Particle update (subordinate to Mesh):

```
For each particle p:
    For each non-zero axis j:
        δV[j].Δθ₀ ← σ.k_effective(p, mesh) · δt
        V[j].θ₀  += δV[j].Δθ₀
    δt += δδt
```

### §7.6 Memory Footprint

```
R_0    : 256 bit
H_*    : 5 × 3 × 256 = 3,840 bit
δH_*   : 5 × 3 × 256 = 3,840 bit
δt_mesh, δδt_mesh : 2 × 256 = 512 bit
Total  = 8,448 bit = 1,056 byte / entire universe
```

About 1 KB for the entire universe. For $N$ particles, total memory is $640N + 1056$ byte.

---

## §8 Verification of W3 Structural Requirements (updated in v2.0)

| W3 Requirement | Satisfied | Basis |
|---|---|---|
| Boundedness | ✓ | int256 represents all phase quantities (§1) |
| Discreteness | ✓ | all fields are integer types (§2–§5) |
| Information conservation | ✓ | δt < 0 enables reverse operation (§5.2) |
| Avoidance of zero division | ✓ | division not used (§6.6) |
| Closure under integer arithmetic | ✓ | only addition/multiplication; cos via discrete modulation table (§7.4) |
| **6-axis symmetry** | **✓ (v2.0)** | frequency interpretation makes all axes equivalent (§1.4, §3) |
| **Transverse-longitudinal separation** | **✓ (v2.0)** | Particle vs Mesh (§7) |
| **Natural representation of acceleration** | **✓ (v2.0)** | δδt is explicit (§5.2) |
| Non-privileging of time | **✓ (v2.0)** | t equal to other axes; time is V[t].θ₀ (§6.4) |
| Finite representability | ✓ | ≤ 640 byte per particle, +1 KB for whole universe (§5.4, §7.6) |

---

## §9 Incompatibilities with v1.0

### §9.1 Changes in Interpretation

- σ.k unified from "wave number / acceleration" to "vibrational frequency"
- δt clarified from "proper time step" to "processing step period"
- "Time" is a derived quantity appearing as V[t-axis].θ₀

### §9.2 Structural Extension

- δV and δδt added to Particle (v1.0's Particle is included as a subset)
- Mesh structure newly introduced (not present in v1.0)

### §9.3 Compatibility with [W8] Table A

The 62-particle classification of [W8] Table A remains valid in v2.0. However, the Higgs, instead of being counted as a Particle, is **counted as an excitation of the Mesh's H_xyz_inphase mode** — an interpretive update. The total count (61 SM states + 1 graviton = 62) is unchanged (only the counting interpretation is updated).

---

## §10 Conclusion

Essence of v2.0:

- **xyzt are frequencies** (natural interpretation in phase space)
- **6 axes are completely symmetric** (4 position-type + R + Q)
- **Transverse (Particle) and longitudinal (Mesh) are separated**
- **Acceleration is δδt** (change in processing period)
- **Time is V[t].θ₀** (a position phase, advancing per particle)
- **All operations close under int256 integer arithmetic**

With this structure, the dynamics rules (D2, separate paper) can describe:

- Gravitational time dilation, polarization rotation (free motion + background coupling)
- Higgs mechanism (coupling to Mesh's H_xyz_inphase)
- Gravitational wave emission (excitation of Mesh's H_xy_anti)
- Coulomb 2-body problem (inter-particle interaction)

all within the same framework.

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
