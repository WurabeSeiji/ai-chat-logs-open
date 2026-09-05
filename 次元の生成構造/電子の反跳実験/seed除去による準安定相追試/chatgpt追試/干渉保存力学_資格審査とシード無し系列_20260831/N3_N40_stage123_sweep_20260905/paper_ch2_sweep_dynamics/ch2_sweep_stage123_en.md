# Chapter 2: The Sweep — Stage 1+2+3 Dynamics and the N=3..40 Runs (H⊥/H Denominator-Control Figure)

(Chapter 2 of the reproduction papers for the N=3..40 stage-1+2+3 sweep. Shares the
Concept DOI 10.5281/zenodo.22317635 with the overview paper; Version DOI
10.5281/zenodo.22317636. Equation numbers continue from Chapter 1 (Eqs. 1–9). Code
blocks are quoted verbatim from the programs; comments inside them remain in Japanese,
as in the originals.)

## 1. Purpose

The central purpose of this series of numerical experiments is **to clarify the
mechanism by which inflation-like development occurs — the phenomenon in which a seed
too small to be measured (dormant fraction H⊥/H ~ 10⁻³⁰) is exponentially amplified
by the dynamics alone over dozens of orders of magnitude until saturation.** Since this
phenomenon was observed in the July canonical runs (N = 40, 300, 1000), the series has
isolated, one factor at a time, which components of the dynamics are necessary
conditions for this development (the confirmed result is the stage-1+2+3 composition;
the audit record is Appendix A).

As the culmination, this chapter takes the 38 static parents (N=3..40) generated and
verified in Chapter 1 as initial data, runs the stage-1+2+3 dynamics (phase-only,
imaginary-part-only generator; real orthogonal rotation; fixed Δτ=2π/den) for
N=3..40 × 6 denominators × 500 steps, records at every step the indicator of
inflation-like development, the **dormant fraction H⊥/H** (definition and its
justification as an indicator: §2.5), the total energy H_total, and the zero closure
|zᵀz|/H, and fixes how this development appears across the whole range of N and the
sweep of the clock (denominator) in the target figure
`fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`.

The descriptive purpose of this paper is complete reproducibility and the fixing of
the correspondence among equations, programs, and data; it contains no physical
interpretation beyond the observed facts. The grounds for the stage-1+2+3 composition
are placed in **Appendix A (audit table of deletion controls)** and **Appendix B
(program lineage)** so as not to interrupt the main text.

## 2. Theoretical Background

The dynamics of this chapter takes the dynamics of the July canonical (old) program as
its starting point and reconstructs it on a new architecture, decomposed into three
components: **stage 1, stage 2, and stage 3**. The theoretical background is therefore
described in the order (2.1) the common state space, (2.2) the mathematics of the old
program's dynamics, (2.3) the mathematics of this program's dynamics (stage
decomposition), (2.4) the mathematical differences between the old and new rotation
maps, and (2.5) measured quantities and conservation laws.

### 2.1 State Space and Adjacency Structure (common to old and new)

**(Eq. 10) Line-graph adjacency matrix** — For edges e=(i,j), f=(k,l) of the complete
graph K_N, A_ef = 1 when e and f share a vertex, 0 otherwise, 0 on the diagonal
(M×M, real symmetric). The edge ordering is the upper-triangular order identical to
Eq. 1. The state is z ∈ ℂ^M (M = N(N−1)/2); the edge phase is θ_e = arg z_e.

### 2.2 Mathematics of the Old Program (July Canonical) Dynamics

The old program is defined by the engine `run_n_scaling_lowrank_v1.py` bundled in this
package (a bit-identical copy of the canonical). One step is the composition of the
following three parts.

**(Eq. 11) Old generator (phase-difference sine, real antisymmetric)** — On
vertex-sharing edge pairs (e,f):

    K_ef(θ) = cos θ_e sin θ_f − sin θ_e cos θ_f = sin(θ_f − θ_e)

The generator is rebuilt at every step from the current phases θ = arg z. The
amplitudes |z_e| do not enter. The vertex decomposition K = C Sᵀ − S Cᵀ = W J Wᵀ
(rank ≤ 2N) is the same as Eq. 3 of Chapter 1.

**(Eq. 12) Power-iteration estimate of σ_max (approximate, history-dependent)** — A
real vector wp is carried across steps:

    wp ← −K(K wp)/‖−K(K wp)‖   (3 iterations)
    σ̂_max = ‖K wp‖

Since −K² is positive semidefinite symmetric with largest eigenvalue σ_max², σ̂_max is
a **power-iteration approximation** of σ_max. Only 3 iterations are used, and wp is
inherited from the previous step (warm start). Strictly speaking, the old dynamics is
therefore not a map of z alone but a **system with a hidden state wp**, and σ̂_max(t)
is an approximate sequence depending on the trajectory history.

Engine implementation (`run_n_scaling_lowrank_v1.py`):

```python
   122	    def sigma_max_power(self, wp, iters=3):
   123	        """warm-start 冪反復による σ_max 推定。wp は前ステップのベクトル（更新して返す）。"""
   124	        for _ in range(iters):
   125	            y = self.kmatvec(wp)
   126	            wp = -self.kmatvec(y)
   127	            nrm = np.linalg.norm(wp)
   128	            if nrm == 0.0:
   129	                return 0.0, wp
   130	            wp = wp / nrm
   131	        sig = np.linalg.norm(self.kmatvec(wp))
   132	        return float(sig), wp
```

**(Eq. 13) Old one-step map (Cayley transform, σ-normalized clock)** — With the
constant γ = tan(π/144):

    K̃ = K / σ̂_max
    z ← (I − γ K̃)⁻¹ (I + γ K̃) z

Engine implementation (2N×2N solve via Woodbury):

```python
   134	    def cayley_step(self, z, sigma):
   135	        """z ← (I-γK̃)^{-1}(I+γK̃) z, K̃ = K/σ。Woodbury で O(N^3)。"""
   136	        gn = GAMMA / sigma
   137	        r = z + gn * self.kmatvec(z)
   138	        A2 = (sigma / GAMMA) * self.J + self.G
   139	        rhs = self.wt(r)
   140	        y = np.linalg.solve(A2, rhs)
   141	        return r - self.w(y)
```

**(Eq. 14) Old eigenplane rotation angle** — The eigenvalues of K are purely imaginary
pairs ±iσ_k (σ_k ≥ 0). The Cayley factor for the eigenvalues ±iσ_k/σ̂_max of K̃ is

    (1 + iγσ_k/σ̂_max) / (1 − iγσ_k/σ̂_max)   (modulus 1)

so the one-step rotation angle of eigenplane k is

    φ_k = 2 arctan( γ · σ_k / σ̂_max )

The rotation angle of the fastest plane (σ_k = σ̂_max) is always fixed at
φ_max = 2 arctan(tan(π/144)) = **π/72**. This is an **adaptive (σ-normalized) clock**
that "ticks time by the system's fastest mode", and the spectral dependence of the
rotation angle undergoes **nonlinear compression** by arctan.

**(Eq. 15) Old conserved quantities** — Since K is real antisymmetric, the Cayley
transform O = (I−γK̃)⁻¹(I+γK̃) is a **real orthogonal matrix** (Oᵀ O = I). A real
orthogonal matrix acts identically on the real and imaginary parts of z = x + iy, so

    ‖Oz‖² = ‖z‖²          (norm conservation)
    (Oz)ᵀ(Oz) = zᵀ OᵀO z = zᵀz   (exact conservation of the zero closure zᵀz)

**(Eq. 15') Relative equilibrium (parent)** — The parent v lies on an eigenplane of
K(arg v) (Chapter 1, Eqs. 4–6); on the flow dZ/dt = KZ we have KZ = −iσZ, i.e.,
Z(t) = e^{−iσt} Z(0). All edge phases advance by the same amount, so phase differences
are invariant, hence K(θ(t)) = K(θ(0)). The parent is a relative equilibrium that
"rotates rigidly without changing phase differences".

### 2.3 Mathematics of This Program's Dynamics (Definitions of Stages 1, 2, 3)

The stage numbering follows the definitions of the N=40 single-factor experiment
series (Appendix A). **Stage 1 provides the baseline architecture**, and **stages 2
and 3 are the two modifications applied to it**.

**(Eq. 16) Stage 1: baseline architecture (explicit matrix, spectral map, fixed Δτ
clock)** — The generator is constructed explicitly as a Hermitian matrix H, and one
step is given by the exact spectral exponential map

    H = V diag(w) V†,   z ← V e^{−i Δτ w} V† z

The clock is a **fixed external clock**

    Δτ = 2π/den,   den ∈ { N−2, N−1, N, N+1, N+2 } ∩ ℕ⁺ ∪ { 124 }

(For small N the series becomes special due to the condition den>0. For N=3 the series
is den=1,2,3,4,5,124; den=1 means Δτ=2π.) In the stage-1-only baseline,
H = A∘(z̄⊗z) (amplitude-weighted) is used. This map is memoryless (one step is a pure
function of z and constants); the old hidden state wp and approximate σ estimation do
not exist.

#### 2.3.1 Design Background (Hypothesis) of the Denominator Series and Its
Implementation Correspondence

The denominator series { N−2, N−1, N, N+1, N+2, 124 } of Δτ = 2π/den is not an
invention of this chapter; it inherits, unchanged, the design of the **denominator
control experiments** of this series (2026-09-03; design instruction document
`ChatGPT_denominator_controls_N3_N40_mixedseed_20260903/CLAUDE_CODE_RUN_INSTRUCTION_N3_N40_20260903.md`).
The background splits into two lines.

**(a) den=124 (fixed-clock control)** — Inherited from the step-size convention
Δτ = 2π/L, **L=124** of the interference-preserving dynamics series (e.g., line 27
`L=124` of `.../program/pass1_parents.py` in the qualification-and-seedless series
folder). It is the control with the clock fixed independently of N (legacy convention).

**(b) den ∈ {N−2..N+2} (hypothesis series of a system-scaled clock)** — The
hypothesis connects to the implementation through the rotation angle of Eq. 20,
ψ_k = (2π/den)·σ_k. Measuring σ_max from the parents' σ spectra saved in Chapter 1
(the `sigma` field of the npz; aggregated in `analysis_sweep_summary_v1.json`,
`sigma_max_by_N`):

| N | σ_max | ψ_max/2π (den=N) | ψ_max/2π (den=124) |
|---|---|---|---|
| 3 | 1.414 | 0.471 | 0.011 |
| 5 | 3.742 | 0.748 | 0.030 |
| 10 | 8.928 | 0.893 | 0.072 |
| 20 | 18.894 | 0.945 | 0.152 |
| 30 | 28.898 | 0.963 | 0.233 |
| 40 | 38.905 | 0.973 | 0.314 |

σ_max(N) grows nearly linearly with N (N − σ_max ≈ 1.1). Therefore:

- In the **den ≈ N series**, the one-step rotation of the fastest eigenplane is placed
  at ψ_max ≈ 2π·(1 − O(1)/N), i.e., in the **nearly-one-full-turn (near-commensurate
  with 2π, stroboscopic) regime**. The five points den = N, N±1, N±2 sweep the detuning
  from this commensurability.
- With **den = 124**, all eigenplanes stay in the small-angle regime (ψ_max/2π ≤ 0.314
  within this sweep), a control close to the small-step approximation of the continuous
  flow dz/dt = Kz (**flow-like regime**).

**Implementation correspondence**: den is determined by the series generation on line
38 (`pairs=[(N+o,…) for o in OFFSETS if N+o>0]+[(124,'124')]`) and enters the dynamics
only through `2.0*math.pi/den` on line 27. den is used nowhere else (as Eq. 20 states,
the effect of den is exhausted by the linear scaling of rotation angles and aliasing).

This chapter does not judge this design hypothesis. The den-dependent observations are
recorded as facts in §6 (denominator dependence of onset; distribution of non-crossing
runs).

**(Eq. 17) Stage 2: amplitude normalization** — The input to the generator
construction is replaced by the unit-modulus ẑ:

    ẑ_e = e^{i θ_e} = e^{i arg z_e}
    Ĥ_ef = A_ef · conj(ẑ_e) ẑ_f = A_ef · e^{i(θ_f − θ_e)}

The generator thereby becomes a function of the **phases only**, as in the old
dynamics (Eq. 11). The amplitudes of the state z itself are retained and continue to
evolve as dynamical variables (the normalization happens only inside the generator
construction, at every step).

**(Eq. 18) Real/imaginary decomposition of Ĥ** — Ĥ is Hermitian and decomposes
uniquely into a real symmetric part and a real antisymmetric part:

    Ĥ = S + iK,
    S_ef = A_ef cos(θ_f − θ_e)   (real symmetric)
    K_ef = A_ef sin(θ_f − θ_e)   (real antisymmetric)

This K is **the same matrix as the old generator (Eq. 11)**.

**(Eq. 19) Stage 3: extraction of the imaginary part (removal of the cos symmetric
part) and real orthogonal rotation** — Instead of Ĥ, the generator

    H₃ = i K = i · Im(Ĥ)

is used (i × real antisymmetric = Hermitian). Substituting into the stage-1 map
(Eq. 16):

    z ← exp(−i Δτ · iK) z = exp(Δτ K) z

exp(ΔτK) is the exponential of a real antisymmetric matrix, i.e., a **real orthogonal
matrix**. It therefore has the same conserved quantities ‖z‖², zᵀz as the old dynamics
(Eq. 15; same proof). Since stage 2 makes K phase-only, the relative-equilibrium
argument of Eq. 15' also holds unchanged: along the parent orbit K is time-invariant
and the parent rotates rigidly.

**(Eq. 20) New eigenplane rotation angle** — The eigenvalues w_k of H₃ = iK are real;
corresponding to ±iσ_k of K, w = ±σ_k. The one-step rotation angle of eigenplane k is

    ψ_k = Δτ · σ_k = (2π/den) · σ_k

The rotation angle is **linear** in the eigenvalue, with no upper bound. The phase
factor e^{−iΔτw} wraps around in w with period 2π/Δτ (**aliasing**): for eigenplanes
with Δτσ_k exceeding 2π, the effective rotation angle is Δτσ_k mod 2π.

### 2.4 Mathematical Differences between the Old and New Rotation Maps (the crux of
this chapter's dynamics)

Both are "real orthogonal rotations generated by the phase-only real antisymmetric K",
but **the spectral map of rotation angles and the clock differ mathematically**.

| Item | Old (Eqs. 12–14) | This chapter = stages 1+2+3 (Eqs. 16–20) |
|---|---|---|
| Eigenplane rotation angle | φ_k = 2 arctan(γ σ_k/σ̂_max) | ψ_k = (2π/den) σ_k |
| Spectral dependence | nonlinear compression by arctan (\|φ\| < π) | linear (with mod-2π aliasing) |
| Clock | σ̂_max-normalized (fastest plane fixed at π/72; system-intrinsic) | fixed Δτ = 2π/den (external; den is the sweep parameter) |
| Obtaining σ | approximate σ̂_max by 3 power iterations (warm start, history-dependent) | exact eigenvalues via eigh (no estimation) |
| Realization of the map | Cayley rational form applied via Woodbury solve (O(N³)/step) | exact exponential applied via spectral decomposition (O(M³)/step) |
| Memory | carries hidden state wp (not a pure function) | memoryless (pure function of z) |
| Conserved quantities | ‖z‖², zᵀz (real orthogonal) | same (real orthogonal) |
| Relative equilibrium (parent) | holds (Eq. 15') | holds (consequence of Eq. 19) |

**(Eq. 21) Correspondence in the small-angle limit** — In the regime Δτσ_k ≪ 1 and
γσ_k/σ̂_max ≪ 1, arctan x ≈ x gives

    φ_k ≈ (2γ/σ̂_max) · σ_k = Δτ_eff · σ_k,   Δτ_eff = 2 tan(π/144)/σ̂_max ≈ (π/72)/σ̂_max

That is, the old dynamics approximately coincides with "Eq. 20 with Δτ_eff
re-normalized by σ̂_max at every step". The differences between the two reduce to
(i) the value of Δτ (externally fixed 2π/den vs. the system-intrinsic (π/72)/σ̂_max),
(ii) the presence of arctan compression and aliasing, and (iii) the exactness of σ.
The measured record of this correspondence is the last row of Appendix A (the
σ-normalized-clock reference experiment: late-segment slope 65.8 steps/decade vs. the
July canonical 64.0 steps/decade). The runs of this chapter use the **fixed Δτ
(stage 1)** throughout; the sweep of den measures the effect of this clock difference.

### 2.5 Measured Quantities and Conservation Laws — Definition of the Dormant Fraction
H⊥/H and Its Justification as the Inflation Indicator

**(Eq. 22) Readout plane** — From the initial state z0 = v + δg (after normalization):

    p = Re z0 / ‖Re z0‖,   q' = Im z0 − (Im z0·p) p,   q = q'/‖q'‖

p, q are **fixed** throughout the run (the Gram–Schmidt orthonormal basis of the real
2-dimensional plane Π = span_ℝ(p, q) spanned by the real and imaginary parts of the
initial state). Since δ = 10⁻¹⁵, Π is essentially **the plane spanned by the parent v**.

**(Eq. 23) Measured quantities** — For the state z at each step:

    H⊥/H = ‖ z − p(p·z) − q(q·z) ‖² / ‖z‖²   (energy fraction outside the readout plane Π)
    H_total = z†z
    closure = |zᵀz| / z†z

**(Eq. 23') Why H⊥/H is the "dormant fraction" and does not pick up the parent's
motion** — The orbit of the pure parent (δ=0) is a relative equilibrium (Eqs. 15',
19): Z(t) = e^{−iσt} Z(0). Expanding,

    e^{−iσt}(x + iy) = (x cos σt + y sin σt) + i (y cos σt − x sin σt),   x=Re Z(0), y=Im Z(0)

i.e., the real and imaginary parts always remain within span_ℝ(x, y) = Π. The
complementary orthogonal projection z − p(p·z) − q(q·z) exactly annihilates any complex
combination inside Π, so **H⊥ is identically 0 along the parent orbit**. Hence H⊥/H
measures only "everything **other than** the parent's rigid rotation" — initially the
out-of-Π component of the seed g (H⊥/H(0) ≈ δ² ~ 10⁻³⁰), and afterwards whatever the
dynamics grows from it. Because the seed is energy invisible at the parent's scale, the
July canonical calls this the **dormant fraction**.

**(Eq. 23'') Identity of definition with the July canonical** — The H⊥/H of this
chapter is the same measurement as f(τ) of the July canonical. Canonical program
`自発的分裂予備実験_v1/run_spontaneous_splitting_largeN_v1.py`:

```python
    57	        Zp = Z - p * (p @ Z) - q * (q @ Z)
    58	        htot = float(np.real(np.conj(Z) @ Z))
    59	        f = float(np.real(np.conj(Zp) @ Zp)) / htot
```

(The construction of p, q is also identical; lines 46–48.) The curves of this chapter
therefore plot the same quantity with the same definition as the July inflation figure
(dormant_growth_large_n_v1.png) and can be compared directly.

**Justification as the inflation indicator** — By the above, the time evolution of
H⊥/H quantifies inflation-like development by three points: "(i) the initial value
lies at the seed scale 10⁻³⁰ (as long as the parent is an equilibrium, there is no
jump at the first step); (ii) from there a straight line on the semilog plot — i.e.,
constant-rate exponential amplification — continues; (iii) it saturates at
O(10⁻²–1)." When (i) fails (a jump to 10⁻⁸–10⁻³ at the first step), that is a mismatch
injection and is distinguished from seed amplification (the deletion controls of
Appendix A are discriminated precisely by this (i)). The onset (Eq. 24) is the
indicator of reaching (iii).

Conservation: by Eq. 19 the map is real orthogonal, so both the numerator |zᵀz| and
denominator z†z of closure are theoretically invariant; the time variation of closure
measures only the accumulation of numerical rounding (measured in §6, G3). Norm
conservation guarantees H⊥/H ∈ [0,1].

**(Eq. 24) onset (bookkeeping indicator)** — The first step number at which
H⊥/H > 0.05 (−1 if none). A tally that does not affect the dynamics.

### 2.6 Correspondence with the Complex-Plane Readout — Design Intent of the Inflation
Figure and the Three Complex-Plane Figures

The complex-plane figures of Chapter 3 (step 0, final step, zoom into the condensed
center) are the readout paired with the inflation figure (H⊥/H curves) of this
chapter. The identity that fixes their relation is given here first.

**Caution (two easily confused "planes")** — The complex plane of the figures is "the
plane on which the M per-edge values z_e ∈ ℂ are plotted", whereas the parent plane Π
= span_ℝ(p,q), the reference of H⊥/H, is "a real 2-dimensional subspace of the state
space ℂ^M". They are distinct, but connected by the following identity.

**(Eq. 25) Per-edge representation of states inside Π** — With x = Re z0, y = Im z0,
v ≈ z0 (δ=10⁻¹⁵), x = (v+v̄)/2 and y = (v−v̄)/(2i) give

    z ∈ span_ℂ{x, y} = Π  ⟺  z_e = a·v_e + b·conj(v_e)   (a, b ∈ ℂ constants independent of the edge)

That is, **motion confined to Π can, in the figure, produce only deformations within
the two-complex-parameter family "uniform rotation/scaling of the parent's star
(a·v_e) plus admixture of its mirror image (b·v̄_e)"**. In particular the pure parent
orbit Z(t)=e^{−iσt}Z(0) (Eq. 15') has a=e^{−iσt}, b=0: **the whole star rotates
uniformly about the origin without changing its shape**. Conversely, the appearance in
the figure of shapes not expressible in this family (e.g., phase dispersal in all
directions = a ring) is the configurational expression of the state having left Π.

By this identity, the three complex-plane figures share the roles of corroborating,
from the configuration side, the start point, end point, and interior of the end point
of the inflation figure (H⊥/H curves):

| Figure | Correspondence to the H⊥/H curve | What it reads |
|---|---|---|
| step-0 figure | start of the curve (H⊥/H ≈ 10⁻³¹) | the contents of Π = the parent's shape (two antipodal pairs, 4-bundle star). That the step-1 figure is indistinguishable from step 0 visualizes f(1) staying at seed scale (§2.5 (i)) |
| final-step figure | saturation region of the curve (H⊥/H ~ 0.05–0.46) | departure to shapes not expressible in the family of Eq. 25 (star → ring) = configuration in which out-of-Π components dominate. onset (Eq. 24) roughly corresponds to when the shape breakdown becomes visible |
| zoom into the condensed center | complement of what the curve does not say | H⊥/H measures only "how much left Π". The zoom resolves "the organization at the destination" — whether angular clusters are condensation onto exactly identical complex values (exact degeneracy) or bundles of finite width |

The implementation of the figures (grid rendering, cluster-extraction algorithm,
line-number correspondence) is described in Chapter 3.

## 3. Implementation Method

- The sweep body `run_N3_N40_stage123_v1.py` is self-contained (no engine import). The
  bundled engine `run_n_scaling_lowrank_v1.py` serves only as the definition body of
  §2.2 (the old mathematics) and for the parent generation of Chapter 1; it plays no
  role in this chapter's dynamics loop.
- This program is derived from the ChatGPT-authored canonical
  `run_and_plot_N3_N40_mixedseed_20260903.py` through **a documented chain of minimal
  diffs only** (Appendix B); the dynamics function `one_step` is identical to the
  stage-3 version that passed the three gates in the N=40 single-factor experiments.
- The initial data are the `Z0` of the Chapter-1 static-parent npz files, loaded
  per N (no regeneration, no modification).
- The input gate `check_sweep_inputs_v1.py` verifies post hoc that the stored `Z[0]`
  of all 228 runs is bit-identical to the static parents' `Z0`.
- The numbers in §6 are transcribed from the output of the aggregation program
  `analyze_sweep_summary_v1.py`, `results/analysis_sweep_summary_v1.json` (no manual
  computation or manual tallying).

**Stage ⇔ equation ⇔ implementation correspondence table** (implementation:
`run_N3_N40_stage123_v1.py`):

| Stage | Content | Equations | Implementation |
|---|---|---|---|
| Stage 1 | explicit matrix + exact spectral map | Eq. 16 | lines 21–22 (`H_of`), lines 27–28 (`eigh` and application of the phase factor) |
| Stage 1 | fixed Δτ=2π/den clock and denominator series | Eq. 16 | `2.0*math.pi/den` on line 27; line 38 (series generation) |
| Stage 2 | amplitude normalization (ẑ=e^{iθ}) | Eq. 17 | `np.exp(1j*np.angle(z))` on line 26 |
| Stage 3 | imaginary-part extraction H₃=i·Im(Ĥ) → real orthogonal rotation | Eqs. 18, 19 | `(1j*np.imag(H))` on line 26 |
| — | readout plane, measured quantities, onset | Eqs. 22, 23, 24 | lines 29–30, 31–32, 45–46 |

(Contrast: the old Eqs. 12–14 exist only in engine lines 122–132 and 134–141 and are
never called in this chapter's loop.)

## 4. Detailed Design

### 4.1 Overall Flow

```
for N in 3..40:
  [initialization]  load static parent Z0 → adjacency A (Eq. 10) → readout plane p,q (Eq. 22) → denominator series (Eq. 16)
  for den in series:
    [loop]  t=0..500: record state and measured quantities (Eq. 23); apply one_step (Eqs. 16–19) if t<500
    [finalization (den)]  save state npz; append onset (Eq. 24) etc. to the summary row
[finalization (all)]  write timeseries/summary CSVs → 8×5 grid figure → RUN_METADATA
```

### 4.2 Overall Data Flow

- **Input**: `Z0` from `parents/parent_static_N{N:05d}_makeparent_20260905.npz`
  (38 files; Chapter-1 products; SHA256 canonical in SHA256SUMS.txt)
- **Parameters**:
  - `STEPS = 500` … steps per run (no early stopping; fixed)
  - `OFFSETS = (-2,-1,0,1,2)` … offsets of the denominator series from N (Eq. 16)
  - denominator 124 … the control denominator common to all N
  - dtype: state complex128 / real float64 (fixed by the assert on line 12)
- **Output**:
  - `results/hm_N{N}_den_{den}_states_500.npz` × 228 (`Z`(501×M), `N`, `denominator`, `steps`)
  - `results/timeseries_64bit_with124_N3_N40.csv` (columns: N, series, denominator,
    step, Hperp_frac, H_total, global_closure; 228×501 = 114,228 rows)
  - `results/summary_64bit_with124_N3_N40.csv` (columns: N, series, denominator,
    onset_gt_0.05, initial, step1, final, max; 228 rows)
  - `results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png` (target figure)
  - `results/RUN_METADATA_N3_N40_stage123.json`

### 4.3 Individual Processes

#### 4.3.1 Initialization (`run_N3_N40_stage123_v1.py`)

Constants and dtype fixing:

```python
    11	STEPS=500; OFFSETS=(-2,-1,0,1,2)
    12	assert np.dtype(np.float64).itemsize==8 and np.dtype(np.complex128).itemsize==16
```

Edges and adjacency (Eqs. 1, 10):

```python
    14	def edges(N):
    15	    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
    16	def adjacency(N):
    17	    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    18	    for e in range(M):
    19	        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    20	    return A
```

Parent loading, readout plane (Eq. 22), denominator series (Eq. 16):

```python
    34	for N in range(3,41):
    35	    # 初期データ: 各 N の静的親ファイルの Z0 を使用
    36	    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)
    37	    A=adjacency(N); p,q=plane(z0)
    38	    pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]
```

```python
    29	def plane(v):
    30	    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q
```

#### 4.3.2 Loop (the one-step dynamics — Eqs. 16–20)

```python
    21	def H_of(z,A):
    22	    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
    23	def one_step(z,A,den):
    24	    # 段3の最小変更（唯一の力学変更点）: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    25	    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    26	    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    27	    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    28	    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
```

- First half of line 26, `np.exp(1j*np.angle(z))`: **stage 2** = Eq. 17 (amplitude
  normalization; the generator input is replaced by ẑ=e^{iθ}).
- Second half of line 26, `(1j*np.imag(H))`: **stage 3** = the imaginary-part
  extraction of Eq. 18 → the generator H₃=i·K of Eq. 19 (`np.imag` is the elementwise
  imaginary part; the imaginary part of the Hermitian Ĥ is the real antisymmetric K).
- Lines 27–28: **stage 1** = Eq. 16 (exact spectral decomposition by `np.linalg.eigh`;
  application of the phase factor e^{−i(2π/den)w} with fixed Δτ=2π/den). The
  eigenplane rotation angle is Eq. 20, ψ_k=(2π/den)σ_k, which differs from the old
  Eq. 14 (arctan compression, σ̂_max normalization) (§2.4).

Measurement and recording (Eq. 23):

```python
    31	def metrics(z,p,q):
    32	    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)
```

```python
    40	        z=z0.copy(); vals=np.empty(STEPS+1,np.float64); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1,np.float64); htot=np.empty(STEPS+1,np.float64)
    41	        for t in range(STEPS+1):
    42	            states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
    43	            if t<STEPS: z=one_step(z,A,den)
```

**Stopping conditions and exceptional behavior (elements not formulated as
equations)**:
- There is no early stopping. Every run consists of a **fixed 501** recordings
  (t=0..500) and 500 applications of the map (lines 41–43). onset (Eq. 24) is used
  only for tallying and never stops the dynamics (the `np.flatnonzero(vals>0.05)` on
  lines 45–46).
- For runs at larger N, numpy RuntimeWarnings (divide by zero / overflow / invalid
  value encountered in matmul) may be printed for the matrix product
  `V@(phase*(V†z))`. These are warnings, not exceptions; execution does not stop. In
  this series, bit-identity of reruns on identical inputs (Chapter 1 G1 and this
  chapter's input gate) confirms they do not affect the results.
- For N=3,4 the `N+o>0` filter of Eq. 16 admits small denominators such as den=1
  (Δτ=2π) into the series. There is no special-casing in the program (a consequence of
  the comprehension on line 38). From the viewpoint of Eq. 20, small den means large
  Δτ — the regime where aliasing (Δτσ_k mod 2π) acts strongly.

#### 4.3.3 Finalization (saving, plotting, metadata)

```python
    44	        np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_{den}_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(den),steps=np.int64(STEPS))
    45	        rows.extend((N,label,den,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1)); ix=np.flatnonzero(vals>0.05)
    46	        summaries.append((N,label,den,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[-1]),float(vals.max())))
```

The CSVs, figure, and metadata (lines 48–66) only write out measured values and render
the figure; they do not affect the dynamics. The figure is an 8×5 grid (lines 55–64),
38 panels used and 2 panels off (line 64).

#### 4.3.4 Input Gate (`check_sweep_inputs_v1.py`)

```python
    18	for N in range(3, 41):
    19	    Z0 = np.load(os.path.join(PARENT_DIR, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0']
    20	    dens = [N + o for o in (-2, -1, 0, 1, 2) if N + o > 0] + [124]
    21	    for den in dens:
    22	        p = os.path.join(RESULT_DIR, f'hm_N{N}_den_{den}_states_500.npz')
    23	        same = bool(np.array_equal(np.load(p)['Z'][0], Z0))
    24	        n_checked += 1
    25	        if not same:
    26	            print(f'MISMATCH: N={N} den={den}')
    27	            ok = False
...
    35	sys.exit(0 if ok else 1)
```

It verifies that the `Z[0]` of every stored npz is `np.array_equal` (bit-identical) to
the corresponding static parent `Z0`; a single mismatch yields exit code 1.

## 5. Execution Results

### 5.1 Reproduction Commands

```bash
cd N3_N40_stage123_sweep_20260905
python3 run_N3_N40_stage123_v1.py        # sweep body
python3 check_sweep_inputs_v1.py         # input gate
python3 analyze_sweep_summary_v1.py      # aggregation for the numbers in §6
# or ./run_all.sh (parents → sweep → gate → aggregation → plots)
```

### 5.2 Execution Environment

- Python 3.9.6 (`.venv/bin/python3`), numpy 2.0.2 (BLAS/LAPACK: macOS Accelerate),
  matplotlib (grid rendering)
- macOS 26.3.1 (arm64)

### 5.3 Execution Time

About **45–50 minutes** for the sweep body (38 N × 6 denominators × 500 steps,
including full state saving; measured 2026-09-05; dominated by the M×M `eigh` × 500 ×
6; N=40 alone about 3 minutes). Input gate and aggregation: under 1 minute each.

### 5.4 Verification Gates

| Gate | Pass condition | Measured | Verdict |
|---|---|---|---|
| G1: input identity | `Z[0]` of all 228 npz bit-identical to the corresponding static parent `Z0` | checked 228 / MISMATCH 0 | **PASS** |
| G2: completion | `ALL DONE` printed, exit code 0 | confirmed | **PASS** |
| G3: conserved quantity (observation) | closure = \|zᵀz\|/H stays small over all runs and steps | step-0 max 5.19e-15; global max 3.72e-13 | **PASS** |

(The G3 numbers are `global_closure_step0_max` / `global_closure_all_max` of
`analysis_sweep_summary_v1.json`.)

### 5.5 Data

| Item | Content |
|---|---|
| Folder | `N3_N40_stage123_sweep_20260905/results/` |
| State npz | `hm_N{N}_den_{den}_states_500.npz` × 228 (tens of KB to ~5.9MB each; actual sizes per SHA256SUMS.txt and the bundle) |
| Timeseries CSV | `timeseries_64bit_with124_N3_N40.csv` (~8.7MB, 114,228 rows) |
| Summary CSV | `summary_64bit_with124_N3_N40.csv` (228 rows) |
| Metadata | `RUN_METADATA_N3_N40_stage123.json` |
| Aggregation | `analysis_sweep_summary_v1.json` (source of the numbers in §6) |
| SHA256 | Canonical values for all files in the bundled `SHA256SUMS.txt`. Sweep program: `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567  run_N3_N40_stage123_v1.py` |

### 5.6 Figures

- Target figure: `results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png`
  (8×5 grid; the **vertical axis of each panel is the dormant fraction H⊥/H** [Eq. 23;
  meaning in §2.5] on a semilog scale, horizontal axis is the step; the six denominator
  curves are overlaid. Straight segments on the semilog correspond to inflation-like
  exponential amplification; horizontal segments to saturation.)
- The complex-plane readout figures (step 0, final, zoom) are treated in Chapter 3.

## 6. Execution Analysis (objective report and observations only; numbers from
analysis_sweep_summary_v1.json)

1. Of the 228 runs, **212 crossed H⊥/H > 0.05 within 500 steps** (onset range 45–481).
2. **In all 228 runs, step 1 remained at seed scale**: H⊥/H(0) ∈ [1.92e-33, 1.60e-31],
   H⊥/H(1) ∈ [3.64e-30, 7.32e-25]. There were **zero** mismatch injections (jumps of
   order 10⁻⁸–10⁻³).
3. Final values of the crossing runs: 0.0501–0.458.
4. Breakdown of the 16 non-crossing runs: 15 runs saturated with final 0.025–0.048,
   **just below the 0.05 threshold** (the 124 series and N+2 series for N≥26, and the
   three series of N=39); only one run was actually slow-growing: N=3, den=124
   (final 1.75e-8).
5. The onset of the 124 series shortens from 318 at N=4 as N increases, entering the
   band 89–160 for N≥15 (crossing N only; the table is `onset_by_series['124']` of
   `analysis_sweep_summary_v1.json`).
6. The zero closure |zᵀz|/H stayed at most 3.72e-13 over all 228×501 recorded points
   (initial max 5.19e-15) — a numerical confirmation of the conservation law of
   Eq. 19.
7. Visual observation of the figure: every panel shows the same family of curves —
   "straight exponential amplification from seed scale → saturation in the
   10⁻³–10⁻¹ range".

## Appendix A: Provenance of the Stage Composition — Audit Table of the Single-Factor
Experiments and Deletion Controls

The dynamics of this chapter (stages 1+2+3) is the composition established by the
single-factor experiment series on N=40 with the static parent (the same Z0 as
Chapter 1). The complete set of experiments (programs, data, figures, README,
SHA256SUMS) is contained in `ChatGPT_denominator_controls_N40_selfcontrol_20260904/`
and is **bundled with the upload of this paper**. The numbers originate from the
summary CSVs and READMEs of the respective results directories.

| Composition | Dynamics (generator / clock) | Level of f(1) | 0.05 crossing (500 steps, 6 dens) | Final values / observation | Output directory |
|---|---|---|---|---|---|
| Stage 1 only (baseline) | H=A∘(z̄⊗z) (amplitude-weighted, cos included) / fixed Δτ | 5.3–6.5e-8 (injection) | 0/6 | stable at ceiling ~1.16e-3 | results_staticparent/ |
| Stages 1+2 (stage 3 deleted) | phase-only Ĥ (cos included) / fixed Δτ | 1.4e-9–8.9e-3 (injection, den-dependent) | 6/6 (τ=4–198) | departure up to 0.94–0.9999 | results_staticparent_phaseonly/ |
| Stages 1+3 (stage 2 deleted) | amplitude-weighted i·Im(H) / fixed Δτ | 1.8–2.3e-8 (injection) | 0/6 | slow drift to ~3.5e-3 | results_staticparent_ampimK/ |
| **Stages 1+2+3 (this chapter)** | phase-only i·K / fixed Δτ | **4.0–9.1e-29 (seed scale)** | 4/6 (τ=276–456) | 0.038–0.10, relaxation curve present | results_staticparent_imK/ |
| Stage 2 at init only (reference) | z0 equimodularized once + amplitude-weighted i·Im(H) / fixed Δτ | 2.9–3.4e-3 (injection) | 6/6 (τ=3–6) | 0.96–0.99, immediate departure | results_staticparent_stage2init/ |
| Stages 2+3 + σ clock (reference) | phase-only i·K / σ-normalized Cayley clock | 5.84e-30 (seed scale) | 0/6 (within 500 steps) | late-segment slope 65.8 steps/decade (July canonical: 64.0) | results_staticparent_sigmaclock/ |

Observation (facts only): f(1) remains at seed scale (~10⁻²⁹) **only** in the
composition where stages 2 and 3 are present **simultaneously**; deleting either one
makes f(1) jump to 10⁻⁸–10⁻³. Moving stage 2 to a one-time initialization produces the
largest jump (10⁻³).

## Appendix B: Program Lineage (chain of minimal diffs and SHA256)

The diffs of each step are recorded in full in the READMEs of the respective packages.

| # | Program | Diff from the previous |
|---|---|---|
| 0 | `ChatGPT_denominator_controls_N3_N40_mixedseed_20260903/run_and_plot_N3_N40_mixedseed_20260903.py` (ChatGPT-authored canonical) | — |
| 1 | `ChatGPT_denominator_controls_N40_selfcontrol_20260904/run_and_plot_N40_only_selfcontrol.py` | 2 lines: OUT destination, `range(40,41)` (N=40 output verified bit-identical to canonical) |
| 2 | `…/run_N40_staticparent_v1.py` | initial data replaced by the Chapter-1 static parent Z0 (input only) |
| 3 | `…/run_N40_staticparent_phaseonly_v1.py` | 1 line in one_step: H built from exp(i·arg z) (stage 2) |
| 4 | `…/run_N40_staticparent_imK_v1.py` | 1 line in one_step: H=1j·Im(Ĥ) (stage 3) |
| 5 | `N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py` (this chapter) | loop range(3,41), per-N parent loading, output destination, figure/metadata names only (dynamics unchanged) |

SHA256 (full digits; transcribed from the output of shasum -a 256):

```
c709c56335d4c67373bff9a3ef6414ea17564d5fbf9c10b7bc9c3724ff091b92  run_and_plot_N3_N40_mixedseed_20260903.py
5a07f354e19985dca0f5de89217e2aa22ac511afaa4b2bb4aa5d93e0c7f9706f  run_and_plot_N40_only_selfcontrol.py
cb5a0ab6db9ae5719eac7aee3b539eb04ca09f5ebfcab6dcce36e8a6e727719e  run_N40_staticparent_v1.py
c1d6d2e60e101f0a99585eefdc22990a087c46246aba112042ac97a7fcbd1a71  run_N40_staticparent_phaseonly_v1.py
a67912b77f7f112731c1eac7612f21b464aed7e2c42431f7b3a7afabb2cd051d  run_N40_staticparent_imK_v1.py
1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567  run_N3_N40_stage123_v1.py
```

---
(End of Chapter 2. Chapter 3, "Complex-Plane Readout Figures", follows.)
