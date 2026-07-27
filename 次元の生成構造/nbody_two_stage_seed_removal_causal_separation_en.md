# Causal Separation of the Time Structure of Three-Direction Formation in Closed N-Body Relational Wave Systems via Two-Stage Seed Removal

**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-07-27
**Version DOI:** [10.5281/zenodo.21614403](https://doi.org/10.5281/zenodo.21614403)
**Concept DOI:** [10.5281/zenodo.21614402](https://doi.org/10.5281/zenodo.21614402)

---

## Claims

Paper 7 discovered that a zero-square-closed relational system undergoes geometric rapid expansion after a long low-change region and, through that process, spontaneously forms a metastable three-directional structure that did not exist initially. This is the central experimental fact of this research. The present paper does not question that three-direction formation. The purpose of Paper 8 is to examine whether another bottom, repeating the same formation process one level below, exists inside the seemingly flat region before the rapid expansion; to verify, at the same time, whether the metastable state after the rapid expansion can produce a further rapid-expansion cycle; and to separate and clarify the time structure leading to three-direction formation and the subsequent long-time evolution.

The verification found that, at least for N=5 within double-precision computation, no separate bottom equipped with independent rapid expansion, arrest, and metastabilization exists inside the initial flat region. Even with the explicit seeds removed, the system starts the same geometric growth from an extremely small state and passes the major amplitude levels at almost the same times. The seemingly flat interval is therefore not a stationary phase but a long latent process continuing into the later rapid expansion.

From this latent stage onward, rank_q=4 in the matrix diagnostics already appears intermittently, persists during the geometric rapid expansion, and becomes a clear and sustained structure in the metastable region after the expansion. After the rapid expansion, new directional occupations grow to finite magnitudes and are retained; that a three-directional structure was actually formed is an undeniable observational fact.

At the same time, detailed direction-lineage tracking revealed that the small direction subspace before the rapid expansion was not simply amplified while fixed to the completed three directions. The direction subspace actually rotates, mixes, and reorganizes during the rapid expansion. This does not weaken three-direction formation; rather, it is a new finding showing that the three directions did not exist as pre-fixed axes but were dynamically selected and formed inside the rapid-expansion process.

This paper further verified the possibility that the metastable state formed after the rapid expansion is an intermediate phase preparing the next rapid expansion. In Paper 7, a second seed was added to this metastable state and the subsequent time evolution was observed; within the observation range in the time direction used since Paper 6, no second rapid expansion was found. In this paper, that second seed itself was removed and the natural long-time evolution of the metastable state was observed. As a result, the metastable state settles into a steady oscillatory state within roughly one period, and even when the observation range along the time axis is extended, no sign connecting to a second geometric rapid expansion was observed.

Therefore, at least within the range observed here, the metastable state after the rapid expansion was not confirmed to be a relay point that automatically produces the next rapid expansion. The metastable state persists for a long time, but no clear growth series heading toward a second rapid expansion has been found inside it. This result indicates that the first rapid expansion may not be a universally repeating cycle, and at the same time suggests that the state after three-direction formation possesses an independent stabilization mechanism.

The phenomenon clarified by this paper is a single sequence of time evolution: a long latent process, geometric rapid expansion, direction reorganization during the expansion, spontaneous formation of a metastable three-directional structure, and long-time metastabilization without a second rapid expansion. What remains unexplained is not the existence of the phenomenon. It is the formation and arrest mechanisms: why the latency time is needed, why the geometric rapid expansion begins, why the reorganization during the expansion selects three directions, and why the formed metastable state settles down instead of proceeding to the next rapid expansion.

An important point about the initial state of this experimental system must be stated explicitly. The experiments initialize the N(N−1)/2 relational waves obtained from N entities by make_parent, and this initialization presupposes the zero-square closure condition, the first axiom of this paper series. The initial state is therefore not an unconditional placement of arbitrary continuous values; it is constructed from the outset as a closed relational system.

As examined in a separate note, this closure condition inherently contains quantization. Closable states do not exist continuously without restriction; the closure condition restricts states to specific phase relations. From this viewpoint, quantization is also demanded of the complex phase, and the phase may be discretized in units of the π period. The implementation in this paper does not explicitly impose such a phase-quantization rule from outside.

Nevertheless, the many relational waves generated by make_parent consequently collapsed, in the initial state, onto two orthogonal directions and formed a stationary two-directional state. The experimenter did not specify two axes as an initial condition. What was given was the relational system and the closure condition; the two-directional structure appeared spontaneously as their consequence.

Therefore, the three-direction formation observed in this paper is not an experiment that added a third axis to two artificially placed axes. It is a phenomenon in which a two-directional state, spontaneously established from zero-square closure, transitions to a metastable three-directional structure through the long latent process, the geometric rapid expansion, and the direction reorganization during the expansion. Whether this initial collapse onto two directions is the direct cause of three-direction formation has not yet been identified in this paper. However, that the initial two directions were not introduced from outside but arose from inside a relational system satisfying the first axiom is an indispensable experimental fact for interpreting the phenomenon.

This paper series has also confirmed, in separate experimental systems, that wave interactions have at least two distinct types. One is a bosonic linear interaction in which waves exchange and mix while remaining superposable; the other is a fermionic interaction that produces nonlinear scattering with reflection and transmission. In the latter interaction, localization is exchanged between an extended wave and a localized wave, and a metastable intermediate state is formed.

However, in the experiments so far, this fermionic wave or nonlinear scattering rule did not arise spontaneously from inside the closed relational system. It is explicitly given from outside by the experimenter as the interaction type. Therefore, at the present stage it cannot be said that the geometric rapid expansion and three-direction formation observed in this paper can be explained by spontaneous emergence of the fermionic interaction.

This paper also searched for precursors of the possibility that, in the process of transition from the initial two-directional state to the metastable three-directional structure, a distinction between linear waves and nonlinearly reflecting-scattering waves emerges from inside — that is, the possibility that the separation into bosonic and fermionic states is related to the onset of the rapid expansion or of direction formation. However, within the time series, direction occupations, rank changes, and metastabilization process observed here, no clear sign was found by which such a distinction could be judged to have formed spontaneously.

Therefore, the possibility that the separation of bosonic and fermionic interactions is the cause of the geometric rapid expansion or of three-direction formation is not excluded, but neither is it supported by this experiment. What was confirmed is that, without installing a fermionic scattering rule, the two-directional state established from zero-square closure transitions to a metastable three-directional structure through spontaneous rapid expansion and direction reorganization. Whether the statistical distinction of waves appears afterwards as a result of this formation process, or requires another degree of freedom not yet observable in this experiment, remains unresolved.

---

## 1. System (definitions and axioms)

**Axiom 1 (nontrivial zero-square-sum closure)** $\sum_e Z_e^2=0,\ Z\neq0$.
**Axiom 2 (finite recurrence)** $U^n=I$. Both axioms follow the basic axiom system [1]. The system definition is identical to Paper 7 [2]; the dynamical code reuses the Paper 7 originals fixed by SHA-256 (§10).

A complex amplitude $Z_e$ is placed on each edge $e$ of the complete graph $K_N$ ($M=N(N-1)/2$). The state is closed on

$$Z\in\mathbb C^M,\qquad \lVert Z\rVert^2=1,\qquad Z^{\mathsf T}Z=0.$$

The generator $K$ is a real antisymmetric matrix depending only on the edge phases $\theta_e=\arg Z_e$, coupling vertex-sharing edges through the sine of phase differences:

$$(KZ)_e=\cos\theta_e\!\sum_{e'\sim e}\!\sin\theta_{e'}Z_{e'}-\sin\theta_e\!\sum_{e'\sim e}\!\cos\theta_{e'}Z_{e'}.$$

The update is the Cayley transform after spectral-norm normalization, $Z\leftarrow(I-\gamma\tilde K)^{-1}(I+\gamma\tilde K)Z$ ($\tilde K=K/\sigma_{\max}$, $\gamma=\tan(\pi/144)$). Since $K$ is antisymmetric this is unitary and conserves $\lVert Z\rVert^2$ exactly.

The parent state $v$ is the self-consistent circularly polarized eigenmode (all edges rotating rigidly at one frequency, $v^{\mathsf T}v=0$). The parent plane is $P_1=\mathrm{span}\{\mathrm{Re}\,v,\mathrm{Im}\,v\}$. The eigenmode residual of the double-precision parent is $2.140\times10^{-13}$ for N=5.

### 1.1 Observables

- **Splitting fraction** $f(t)=\lVert Z-\Pi_{P_1}Z\rVert^2/\lVert Z\rVert^2$: total occupation outside the parent plane.
- **crossing**: first passage of the existing Paper 7 threshold $f>0.05$. A reference coordinate, not an event definition.
- **Metastable start**: crossing+3000 (Paper 7 convention).
- **$q_1\ge q_2\ge q_3\ge q_4$**: singular values from the dominant-plane Gram reduction. **rank_q** is the effective rank under the existing relative threshold. rank_q is the rank of $Q=[B_0\mid B_{\mathrm{dom}}]$ (an $M\times4$ matrix) and is structurally capped at 4.
- **Five-component occupation decomposition**: occupations of directions 1–4 (two parent-plane directions plus two additional directions), the residual rotating subspace, and the kernel. The occupations sum to 1 at every step (projection closure error 0).
- **Direction subspace** $D_{34}(t)$: the orthonormal basis of the two additional directions returned by the existing `s4_new_dirs`. The canonical comparison object is not the columns but the projector

$$P_{34}(t)=D_{34}(t)D_{34}(t)^{\mathsf T},\qquad \mathcal O(P,Q)=\tfrac12\operatorname{tr}(PQ),$$

with $\mathcal O=1$ meaning subspace identity. Individual direction 3/4 columns depend on basis choice (column exchange, sign) and are not used for lineage judgment.

- **First-passage coordinates** $\tau_k=\min\{t\mid f(t)\ge 10^{-k}\}$: first passage of each decimal level. The interval amplification rate is defined as $\ln 10/\Delta\tau$. These are coordinates of the amplification process, not event thresholds.

### 1.2 The two explicit seeds

The Paper 7 orbit contains two artificial inputs.

1. **Initial seed**: a zero-closure pair in the parent kernel added at $\delta=10^{-15}$ ($f(0)\approx\delta^2=10^{-30}$).
2. **Metastable seed**: a single transverse perturbation $\varepsilon\eta_\perp$ ($\varepsilon=10^{-8}$) injected once at $t_1=\mathrm{crossing}+3000$ (the perturbation-response measurement of Paper 7 Figure 4).

This paper removes the two independently and separates which observed facts causally depend on the seeds.

### 1.3 Structure of the initial state: two-directionality fixed by zero-square closure [Result]

For $Z=X+iY$,

$$Z^{\mathsf T}Z=\lVert X\rVert^2-\lVert Y\rVert^2+2i\,X\cdot Y,$$

so the zero-square closure $Z^{\mathsf T}Z=0$ is equivalent to

$$\lVert X\rVert=\lVert Y\rVert,\qquad X\cdot Y=0.$$

Every nonzero state satisfying the closure condition necessarily fixes a pair of equal-norm orthogonal directions in the real $M$-dimensional relational space. The two-directional structure is not a choice of initial condition but a consequence of the first axiom.

Measurements of the `make_parent` output are given in Table 1.

*Table 1: Structure of the parent state (originals verified by SHA-256; same PRNG as Stage A2a)*

| Quantity | N=5 (M=10) | N=40 (M=780) |
|:--|:--|:--|
| Nonzero complex components / M | 10 / 10 | 780 / 780 |
| min $\lvert v_e\rvert$ | 0.2887 | 0.0324 |
| rank$[\mathrm{Re}\,v,\ \mathrm{Im}\,v]$ | **2** | **2** |
| $\lVert\mathrm{Re}\,v\rVert=\lVert\mathrm{Im}\,v\rVert$ | $1/\sqrt2$ (error $10^{-15}$) | $1/\sqrt2$ (error $10^{-15}$) |
| $\mathrm{Re}\,v\cdot\mathrm{Im}\,v$ | $-4.4\times10^{-16}$ | $-3.6\times10^{-16}$ |
| $\lvert v^{\mathsf T}v\rvert$ | $1.1\times10^{-15}$ | $2.3\times10^{-15}$ |
| Invariant-plane residual $\lVert KX-\mu Y\rVert$ | $2.0\times10^{-13}$ | $1.1\times10^{-13}$ |
| Decomposition dimensions of the initial $K$ (parent / rotating complement / kernel) | 2 / 4 / 4 | 2 / 78 / 700 |

Table 1 establishes three facts. First, all $M$ relational waves are nonzero; this is not a state keeping only two particular components. Second, the real rank of the whole state is 2: all relational waves are collapsed onto a single two-dimensional rotation plane (an invariant plane of the generator $K$, residual $10^{-13}$). Third, `make_parent` imposes only the closure condition and self-consistency (being an eigenmode of the $K$ built from its own phases); neither a phase-quantization rule nor the coordinates of two axes are given from outside.

The separate examination referred to in the Claims — that the closure condition inherently contains quantization — is the expository note [3]. From the two requirements of zero-square-sum constraint and scale anonymity, the state space becomes a projective quadric in complex projective space (a compact Kähler manifold), and discrete spectra, integer Chern classes, and the qubit isomorphism at $N=3$ follow from citations of known mathematics alone. The implementation here does not impose that quantization rule from outside; the two-direction collapse in Table 1 is a measurement of the double-precision dynamical system.

Whether the initial collapse onto two directions is the direct cause of three-direction formation is not identified (§9).

---

## 2. Experimental design

### 2.1 Conditions and experimental systems

| Condition | Initial seed | Metastable seed | Initial state |
|:--|:--:|:--:|:--|
| A (fully seedless) | OFF | OFF | $Z_0=v$ (no random numbers consumed) |
| B (initial only) | ON | OFF | $Z_0=(v+\delta g)/\lVert\cdot\rVert$ |
| D (Paper 7 equivalent) | ON | ON | same as B; single injection of $\varepsilon\eta_\perp$ at $t_1$ |

- **Preliminary experiment** (N=5, 40, 300, t≤55000): full comparison of conditions A/B/D.
- **Stage A2a** (N=5, every-step recording, t≤5000): high-density observation of the fully seedless orbit, starting from the numerical floor $f(0)=3.275\times10^{-33}$ of $Z_0=v$.
- **Stage A2c / A2d** (N=5): direction-lineage tracking of $P_{34}(t)$ for the seeded / fully seedless orbits.
- **Long-horizon control** (N=5, 40, 300, fully seedless, t≤110000): search for a second rapid expansion.
- **Transition anatomy** (N=5, Paper 7 reproduced orbit): temporal separation of rank_q, occupations, and f.

### 2.2 Hypotheses tested

The hypotheses remaining at the close of Paper 7 were fixed in advance.

- **H-lower-bottom**: a lower metastable phase (a second bottom) with an isomorphic growth-and-arrest process is hidden inside the flat region before the rapid expansion.
- **H-seed-causation**: the explicit seeds determine the occurrence, the timing, or the first generated direction of the splitting.
- **H-scaled-copy**: the small direction subspace before the rapid expansion is a scaled copy of the late three directions.
- **H-repeating-cycle**: the metastable state after the rapid expansion is an intermediate phase preparing the next rapid expansion.

The results reject or fail to support all four, and in doing so establish the time structure stated in the Claims.

---

## 3. Result 1: the full sequence reproduces with the seeds removed [Result]

Condition A (fully seedless) reproduced crossing, rank_q increase, and metastabilization at every N (Table 2).

*Table 2: Conditions A/B/D at all N (t≤55000)*

| N | Cond. | crossing | metastable start | max_f | final_f | final rank_Q | mean_meta(f) | std_meta(f) |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| 5 | A | 1166 | 4166 | 0.9655 | 0.8754 | 4 | 0.8744 | $1.29\times10^{-2}$ |
| 5 | B | 1167 | 4167 | 0.9600 | 0.8053 | 4 | 0.8054 | $2.87\times10^{-2}$ |
| 5 | D | 1167 | 4167 | 0.9600 | 0.8026 | 4 | 0.8034 | $2.88\times10^{-2}$ |
| 40 | A | 2011 | 5011 | 0.1938 | 0.1938 | 4 | 0.1934 | $7.38\times10^{-4}$ |
| 40 | B | 2011 | 5011 | 0.2035 | 0.2035 | 4 | 0.2030 | $1.29\times10^{-3}$ |
| 40 | D | 2011 | 5011 | 0.2035 | 0.2035 | 4 | 0.2030 | $1.29\times10^{-3}$ |
| 300 | A | 4849 | 7849 | 0.0898 | 0.0898 | 4 | 0.0882 | $3.53\times10^{-3}$ |
| 300 | B | 4844 | 7844 | 0.0861 | 0.0861 | 4 | 0.0848 | $2.94\times10^{-3}$ |
| 300 | D | 4844 | 7844 | 0.0861 | 0.0861 | 4 | 0.0848 | $2.94\times10^{-3}$ |

![f comparison N=5](第8論文_二段階seed除去による準安定相の因果分離/figures/fig01_f_compare_N00005.png)

*Figure 1: f(t) for conditions A/B/D at N=5. The crossing difference between seedless (A) and seeded (B/D) is 1 step. N=40 and 300 are isomorphic.*

The crossing differences are 1 step at N=5, 0 steps at N=40, and 5 steps at N=300. The explicit initial seed contributes neither to the occurrence nor to the timing of the splitting (first part of H-seed-causation rejected).

---

## 4. Result 2: no separate bottom inside the latent region [Result]

The fully seedless N=5 orbit (Stage A2a, every-step recording) was inspected over its entire structure from the initial numerical floor $f(0)=3.275\times10^{-33}$ to the main growth.

### 4.1 Sign inspection of every-step differences

For all 1129 differences from step 0 to just before the first passage of $f\ge10^{-2}$ (step 1129):

- positive differences: **1129 (all)**
- negative differences: **0**
- zero differences: **0**

The shape "growth → arrest → lower metastable shelf → regrowth" exists nowhere in the every-step data from the initial floor to the main growth. The flatness in linear display is an appearance caused by the smallness of the observable; in logarithmic display it is one continuous rise (Figure 2).

![seedless log10 f](第8論文_二段階seed除去による準安定相の因果分離/paper8_stage_A2a_seedless_N5/figures/figure02_seedless_log10_f.png)

*Figure 2: $\log_{10}f$ of the fully seedless N=5 orbit. A single continuous rise from the numerical floor.*

### 4.2 First-passage series and amplification rate

The passage steps from $10^{-30}$ to $10^{-23}$ rise as 1, 2, 5, 14, 36, 68, 107, 150, and over all 21 intervals of $10^{-23}\to10^{-2}$ the passage time converges to **46 or 47 steps per decade**. The median interval exponential rate is $4.95\times10^{-2}$ (natural log per step).

Linear regression of $\ln f$ over the 980 points from $f\ge10^{-23}$ (step 150) to $f\ge10^{-2}$ (step 1129) gives

$$\lambda=0.04937,\qquad R^2=0.9999993.$$

The per-step ratio of f is $e^{\lambda}=1.0506$; the amplitude ratio is $e^{\lambda/2}=1.0250$. The main part of the latent region is described by a single geometric law (Figure 3).

![decade rate comparison](第8論文_二段階seed除去による準安定相の因果分離/paper8_stage_A2a_seedless_N5/figures/figure05_decade_growth_rate_comparison.png)

*Figure 3: Mean exponential rate per decade. Seeded and seedless coincide at every level.*

### 4.3 Judgment

**H-lower-bottom is rejected.** The scope of the rejection is N=5, double precision, every-step recording (t≤5000). Structure below the double-precision floor, and unperturbed departure in exact arithmetic, are outside the scope of this paper (§9).

---

## 5. Result 3: the seeds select neither timing nor initial direction [Result]

### 5.1 Timing

Comparing first passages of the seeded and seedless orbits at all levels, the first passage of $f\ge10^{-12}$ is absolute step 662 for both, and the per-decade passage-time difference is at most 1 step over all intervals (Figure 4).

![absolute step comparison](第8論文_二段階seed除去による準安定相の因果分離/paper8_stage_A2a_seedless_N5/figures/figure03_seeded_vs_seedless_absolute_step.png)

*Figure 4: f(t) of seeded and seedless orbits overlaid on absolute steps. They coincide without any translation.*

The initial values are inverted. The seeded orbit starts from $f(0)=1.066\times10^{-30}$ and the seedless orbit from $f(0)=3.275\times10^{-33}$, and the two merge onto the same geometric series within a few hundred steps. The explicit seed $\delta=10^{-15}$ is smaller than the eigenmode residual $2.140\times10^{-13}$ of the double-precision parent. Even in the seeded orbit, what governed the effective minimal transverse component was not the explicit seed.

### 5.2 Initial direction

The same-step direct comparison of $P_{34}(t)$ between the seedless and seeded orbits is given in Table 3 (Stage A2d).

*Table 3: Same-step overlap, seedless vs. seeded*

| Interval | Samples | Overlap median | Overlap minimum |
|:--|--:|--:|--:|
| Pre-crossing, resolvable | 40 | **0.999999987** | 0.999976 |
| crossing–1799 | 128 | 0.999997 | 0.994306 |
| 1800–2500 | 141 | 0.999852 | 0.990155 |
| after 2500 | 500 | 0.265980 | 0.203573 |

From before the rapid expansion until just after the establishment of the three-direction closure, the direction subspaces of the two orbits are identical. **The explicit initial seed does not select the first generated direction subspace** (second part of H-seed-causation rejected).

The divergence after step 2500 (overlap first below 0.95 at step 2690, below 0.5 at step 3185) is a long-time orbit difference separated from the initial selection of the generated direction. After this divergence, rank_q=4 and the occupation levels are retained in both orbits (§6.3).

---

## 6. Result 4: rank_q=4 onset, occupation growth, and the f expansion occur at different times [Result]

### 6.1 rank_q=4 appears at the numerical floor

- Seedless orbit: the first rank_q=4 is at saved step **120**, where $f=2.089\times10^{-24}$; the bracketing saved occupation rows (steps 107 / 125) have direction 3/4 occupations at the $10^{-30}$–$10^{-29}$ level.
- Seeded orbit: the first saved rank_q=4 is at step **265**, $f=3.127\times10^{-21}$, with bracketing direction 3/4 occupations at $4\times10^{-28}$–$8\times10^{-28}$. The number of rank_q=4 saved rows before crossing is 97.

rank_q=4 appears while splitting fraction and direction occupations are at the numerical floor. Therefore

$$\mathrm{rank_q}=4\quad\not\Rightarrow\quad\text{establishment of a direction structure with finite occupation},$$

and the onset of rank_q=4 cannot serve as the time of three-direction establishment (Figure 5).

![q ratios and rank_q](第8論文_二段階seed除去による準安定相の因果分離/paper7_N5_transition_anatomy/figures/figure06_q_ratios_and_rank_q_0_3000.png)

*Figure 5: $q_3/q_1$, $q_4/q_1$ and rank_q at N=5 (steps 0–3000). The early rank_q response and the occupation growth occur at different times.*

### 6.2 Quantities that grow with delays after crossing

After crossing=1167, direction 3/4 occupations, the kernel, and the q ratios do not change at the same speed (Table 4).

*Table 4: Observables after crossing (N=5, actual saved values)*

| step | f | d3 occ. | d4 occ. | $q_3/q_1$ | $q_4/q_1$ |
|--:|--:|--:|--:|--:|--:|
| 1200 | 0.133 | 0.0027 | 0.0035 | 0.279 | 0.016 |
| 1400 | 0.585 | 0.023 | 0.195 | 0.608 | 0.248 |
| 1800 | 0.840 | 0.213 | 0.210 | 0.733 | 0.507 |
| 2500 | 0.801 | 0.319 | 0.045 | 0.752 | 0.549 |

The f crossing, the finite growth of direction occupations, and the growth of the q ratios are separated over hundreds of steps; a single-event description "the three directions were established at one time" is not supported by the data. That direction 3/4 occupations grow to finite magnitudes (the $10^{-1}$ level) after the rapid expansion and are retained is an observational fact, as in Table 4 and Figure 5.

### 6.3 What is metastable is the occupation structure

Conditions A and B diverge as orbits in the metastable region (§5.2), while rank_Q=4, the ceiling-level f, and the finiteness of the direction occupations are preserved in both (Table 2). The three-direction closure is not three fixed axes; it is a dynamical closure that preserves the number of occupied dimensions and the occupation distribution while its internal directions can mix over long times.

---

## 7. Result 5: the direction subspace reorganizes during the rapid expansion [Result]

### 7.1 The pre-expansion direction is not a scaled copy of the late direction

In the fully seedless orbit (Stage A2d), the overlap between the $P_{34}$ that became resolvable before the rapid expansion (resolvability criterion $\min(q_3,q_4)/q_1\ge10^{-6}$; 40 samples at steps 985–1165) and the late representative subspace (top two eigenvectors of the mean projector over steps 1800–2500) is

$$\mathcal O_{\mathrm{early,late}}\ \text{median}=0.1425\qquad(\text{range}\ [0.1210,\ 0.1445]).$$

The seeded orbit (Stage A2c) gives the same level (median 0.1495). A quantity that would approach 1 for identical subspaces stays at 0.14. **The small direction subspace before the rapid expansion is not a scaled copy of the late three directions** (H-scaled-copy rejected; Figure 6).

![seedless direction lineage](第8論文_二段階seed除去による準安定相の因果分離/paper8_stage_A2d_seedless_direction_lineage_N5/figures/figure01_seedless_early_vs_late_lineage.png)

*Figure 6: Early-vs-late direction lineage of the fully seedless N=5 orbit. The overlap between the early $P_{34}$ and the late representative stays at 0.14.*

### 7.2 The reorganization is a real rotation during the expansion

In the fixed observation window steps 900–1400 of the seeded orbit, the sum of maximal principal angles between consecutive-step $P_{34}$ (rotation path length) is 3.871 rad, and the maximum of the consecutive-step maximal principal angle is 0.372 rad. Column exchanges are 0; this is not a relabeling of basis columns but rotation and mixing of the subspace itself. Stages A2c/A2d classify this lineage as **ROTATING_OR_MIXED_LINEAGE**.

### 7.3 Reclassification of Paper 7 Figure 4

The overlap between the transverse perturbation directions $T_\perp$ applied in the metastable region in Paper 7 Figure 4 and the late representative direction is, for the three perturbation seeds,

$$\mathcal O(D_{34}^{\mathrm{late}},T_\perp)=0.0330,\ 0.0301,\ 0.0182.$$

The artificially excitable response channels coincide neither with the naturally occupied three directions nor with the precursor directions before the rapid expansion. What Paper 7 Figure 4 measured was not "sprouts of additional directions arising naturally" but **response channels outside the three-direction closure that are excitable but not naturally occupied**.

---

## 8. Result 6: the metastable state does not proceed to a second rapid expansion [Result]

### 8.1 Removal of the second seed: the metastable oscillation is not a product of the seed

Under condition B (no metastable seed), the temporal variation in the metastable region continued to the final step (std_meta: $2.87\times10^{-2}$ at N=5, $1.29\times10^{-3}$ at N=40, $2.94\times10^{-3}$ at N=300). The difference from condition D (single injection of $\varepsilon=10^{-8}$ at $t_1$) is exact agreement within recording precision at N=40 and 300, and about $10^{-3}$ at N=5 (Figure 7).

![metastable B vs D](第8論文_二段階seed除去による準安定相の因果分離/figures/fig04_metastable_B_vs_D_N00005.png)

*Figure 7: Conditions B and D in the metastable region at N=5. The transverse perturbation does not change the metastable oscillation.*

The metastable oscillation is not a product of the transverse-perturbation seed; it is a consequence of the closed dynamics. After the injection in condition D, no second crossing or rank increase appears within the observation range t≤55000.

### 8.2 Long-horizon control: no sign of a second rapid expansion up to t=110000

The fully seedless orbits were extended to t=110000 at N=5, 40, 300.

- N=40, 300: the raw f, $q_3$, $q_4$ over t=55000–110000 form **no extremum at all** and approach their terminal values monotonically. The peak-to-peak over t=90000–110000 is $f=1.0\times10^{-11}$ at N=40 and $f=2.1\times10^{-11}$ at N=300 (recording-precision level).
- N=5: the metastable oscillation (f range [0.807, 0.966] over steps 4166–55000) decays by t=55000, and over the extension 55000–110000 f is constant at 0.875392 (peak-to-peak $1.0\times10^{-11}$, recording-precision floor). In the extension, steps exceeding the metastable maximum: 0; drops below f=0.05: 0; second crossings: 0.
- The late relaxation cannot be absorbed by a single exponential; a double exponential removes the residual valley (AIC improvements at N=300: f=2848.6, $q_3$=1568.8, $q_4$=2172.0; Figure 8).
- The early metastable region contains raw-value reversals of a damped oscillation (N=40: $q_3$ maximum at t=1650, minimum at t=1700). N=300 is recorded at 100-step intervals, leaving shorter oscillations possibly unresolved.

![long-horizon residual](第8論文_二段階seed除去による準安定相の因果分離/paper7_seedless_natural_figures3_4_v1/outputs/long_horizon_110000/figures/figure_long_horizon_one_exp_residual_x20.png)

*Figure 8: Single-exponential residual of the long-horizon control. The remaining structure is explained by double-exponential relaxation; there is no reversal toward re-amplification.*

Within the observation range **H-repeating-cycle is not supported**. The metastable state is a long-lived phase that, after the initial damped oscillation, relaxes multi-exponentially toward a steady level; it contains no growth series heading toward a second geometric expansion.

### 8.3 No precursor of a bosonic/fermionic distinction appears in the observation range

This series has confirmed in separate experimental systems that wave interactions have two types. In contrast to the linear interaction in which waves exchange and mix while remaining superposable, the nonlinear exchange-interference scattering with reflection and transmission exchanges localization between an extended wave and a localized wave [4] and forms metastable intermediate states (the white-cat/black-cat/gray-cat interface) [5]. In those experiments the scattering rule is given from outside by the experimenter as the interaction type.

The experimental system of this paper has no direct observable that discriminates this type. The update is a single unitary step on the whole state by the state-dependent generator; it has no test direction constructing a per-wave-pair reflection/transmission scattering process. What follows is therefore not a designed test but a passive precursor search over the existing observables (time series, direction occupations, rank, metastabilization process).

The operational content of the search is whether an additional discrete direction beyond three, or a new discrete channel, becomes naturally occupied. If the separation of linear and nonlinearly scattering waves arises internally along with the rapid expansion or direction formation, its first trace appears as natural occupation of additional directions. The result is Table 5.

*Table 5: Residual rotating occupation in the metastable region of the fully seedless natural orbits (occupation outside directions 1–4 and the kernel)*

| N | Dimension of residual rotating subspace | Metastable max | Terminal value |
|--:|--:|:--|:--|
| 5 | 2 | $1.5\times10^{-3}$ | $4.4\times10^{-16}$ |
| 40 | 76 | $4.7\times10^{-6}$ | 0 |
| 300 | 596 | $9.7\times10^{-8}$ | $2.1\times10^{-17}$ |

Since rank_q is structurally capped at 4 (§1.1) and cannot detect additional directions, the judgment rests on this direct measurement. The residual rotating occupation grows at no N and decays to numerical zero (Figure 9). The response channels outside the closure are excitable but not naturally occupied (§7.3). The long-horizon control (t≤110000) likewise shows no growth series of a new channel.

![five-component occupation](第8論文_二段階seed除去による準安定相の因果分離/paper7_seedless_natural_figures3_4_v1/outputs/figures/figure3_compare_N5_N40_N300.png)

*Figure 9: Five-component occupation decomposition of the fully seedless natural orbits (N=5, 40, 300). No direction outside the three directions and the kernel becomes naturally occupied.*

No sign by which a spontaneous separation into bosonic and fermionic states could be judged exists in these observables. That this separation causes the rapid expansion or three-direction formation is not excluded, and not supported. What is confirmed is that the transition from the two-directional state to the metastable three-directional structure occurs without installing a fermionic scattering rule.

---

## 9. Boundary conditions (what this paper does not prove)

1. **Double-precision floor**: the minimal input of the seedless orbit is the eigenmode residual of the parent state ($2.14\times10^{-13}$) and rounding error. Unperturbed departure in exact arithmetic — the physical origin of the first minimal difference — is not separated here. This is not a condition weakening the central claims but the next object to be separated.
2. **N and time range**: the every-step inspection for the lower-bottom rejection (§4) is N=5, t≤5000. The long-horizon control (§8.2) is N=5, 40, 300, t≤110000. Repetition with still longer periods is not logically excluded, but no raw-data reversal indicating such a candidate exists in the present range.
3. **Recording intervals**: q is saved every 5 steps and occupations every 25–100 steps; simultaneity claims are limited to the actually saved records. Sample-to-sample rotation quantities in the lineage analysis are not one-step maxima.
4. **Incomplete control**: the numerical-resolution sweep (control with an artificial quantization operator) has incomplete aggregation. The completed part is a control in which cutting off the minimal amplification source at finite resolution pins the orbit to the parent state; it is not evidence of a second bottom.
5. **Correspondence with physical XYZ**: the correspondence between the formed three-direction closure and physical spatial directions is a separate connection problem handled with projection Gram matrices, principal angles, and interference quantities. What this series establishes is that three directions became readable from the $M$ phase relations, not that they are the xyz of a background space. The physical canonical readout is the interference of direction planes (inner products of normals).
6. **Causality between the initial two directions and three-direction formation**: that the initial state collapses onto two directions as a consequence of the first axiom (§1.3) is a measured fact, but whether this collapse directly causes three-direction formation is not identified. A systematic control from generic closure-only states (`zero_closure_generic` initialization without the eigenmode condition) is the next task.
7. **Inherent phase quantization**: the viewpoint of the expository note [3] that the closure condition inherently contains quantization is not imposed as an external rule in this implementation and is not verified by the present data.
8. **Untestability of the statistical distinction**: the precursor search of §8.3 is a passive search constrained by the absence of a direct observable discriminating the interaction type; it is not an experiment testing spontaneous formation of the statistical distinction. Testing it requires designing an observation system that constructs scattering processes inside the closed relational system (generation of localized wave packets and readout of collisions) — a connection problem between the exchange-interference scattering systems [4,5] and the present system.

---

## 10. Reproducibility

- The dynamical code reuses the Paper 7 originals (`run_n_scaling_lowrank_v1.py` and others) fixed by SHA-256 and imported read-only; the original folders are unchanged. All hashes are recorded in the stage reports and `config/source_file_hashes.json`.
- The Paper 7 N=5 orbit was reproduced bitwise in both the f and q series (Stage A0). The seedless orbit is bitwise identical across two independent executions (Stage A2a exec 1/2). The N=300 five-component decomposition was cross-checked against the existing condition A with all rows of the 27 common columns identical.
- Numerical health: normalization error $\le1.4\times10^{-14}$; zero-square closure $\le1.8\times10^{-10}$ (immediately after the condition-D injection; otherwise $\le1.3\times10^{-13}$); projection closure error 0; occupation-sum error 0.
- Canonical sources of the numbers: Table 1, §4.1, the §4.2 regression, §6.1 seedless, §8.2 N=5, and Table 5 are in `paper8_draft_supporting_analysis_v1/reports/`. Table 2 and §8.1 are in `reports/preliminary_seed_ablation_report.md`. §4.2 and §5.1 are in `paper8_stage_A2a_seedless_N5/reports/`. Table 3 and §7.1 are in `paper8_stage_A2d_seedless_direction_lineage_N5/reports/`. §6 is in `paper7_N5_transition_anatomy/reports/`. §7.2 and §7.3 are in `paper8_stage_A2c_direction_lineage_N5/reports/`. §8.2 is in `paper7_seedless_natural_figures3_4_v1/validation_report.md`.
- Environment: Python 3.9.6, NumPy 2.0.2, macOS arm64. Random numbers via `numpy.random.default_rng` with fixed seeds.

---

## 11. Conclusion

The phenomenon discovered in Paper 7 — the long low-change region, the geometric rapid expansion, and the spontaneous formation of the metastable three-directional structure — reproduces at all of N=5, 40, 300 with both explicit seeds removed. The seeds cause none of the occurrence, the timing, the first generated direction, or the metastable oscillation.

The initial two-directional state at the starting point is not two axes placed by the experimenter. Zero-square closure forces equal-norm orthogonal two-directionality on every nonzero state, and all $M$ relational waves of `make_parent`, while nonzero, collapse spontaneously onto a single rotation plane of real rank 2. The observed sequence of this series begins with a two-direction closure spontaneously established from the axioms and ends with the dynamical formation of a three-direction closure.

The seemingly flat region before the rapid expansion is a single continuous rise with zero negative differences under every-step inspection — a long latent process with no lower bottom inside. The early onset of rank_q=4 is a response at the numerical floor and cannot be identified with direction establishment. The direction subspace rotates, mixes, and reorganizes during the rapid expansion; the three directions are dynamically selected and formed inside the expansion process, not laid down in advance. The formed metastable state relaxes multi-exponentially after the initial damped oscillation and produces no second rapid expansion up to the observation range t=110000.

The open problems are narrowed to three. (i) Why does the geometric growth begin after the long latency — including the physical origin of the first minimal difference below the double-precision floor. (ii) Why does the reorganization during the expansion select three directions. (iii) Why does the formed metastable state settle down instead of proceeding to the next expansion. Direct verifications are: for (i), reproduction of the latent region in exact or high-precision arithmetic; for (ii), comparison of the most unstable eigenplane of the transverse Jacobian in co-rotating coordinates (the 144-step Floquet matrix) with the measured precursor $P_{34}$; for (iii), identification of the relaxation spectrum on the metastable shell.

The possibility that the separation of bosonic and fermionic interactions underlies the phenomenon remains neither excluded nor supported (§8.3). Three-direction formation without an installed fermionic scattering rule is confirmed. Whether the statistical distinction of waves appears afterwards as a result of this formation process, or requires a degree of freedom not observable in this system, is a task for the next experimental series, together with the design of an observation system constructing scattering processes inside the closed relational system.

### 11.1 Relation to prior work

- **Latency followed by geometric growth** shares its time structure with modulational instability, in which a uniform carrier destabilizes exponentially into transverse sidebands from a noise floor [8]. In the present system, however, what grows is not spatial sidebands but the direction subspace itself, and the directions reorganize during amplification (§7).
- The **transverse spectrum of a state-dependent generator** is a toolset shared with the locked-state spectrum analysis of coupled phase oscillators [7]. The present system, however, carries the norm-structural constraint of zero-square closure, and what is read out is direction structure, not synchronization.
- **Long-time metastabilization without a second expansion** (§8) can be read in the vocabulary of the FPUT metastable state [6] and prethermalization in isolated systems [9] — a long-lived nonthermal intermediate state. No recurrence corresponding to the FPUT type appears within t≤110000.
- The prior work closest in claim shape — three directions alone expanding out of a symmetric candidate space — is the dynamical emergence of (3+1)-dimensional spacetime in the IIB (IKKT) matrix model [10], where three directions extend via spontaneous breaking of SO(9,1)/SO(10) symmetry. Mechanism and level both differ: against a path integral over matrix degrees of freedom, this paper is a deterministic closed dynamics of N-body relational waves, and the three directions are not dimensions of a background space but direction structure read out from the interference of relational waves. Whether the shared selection of "3" is another expression of one mechanism is left as an open question.

No prior work corresponds to the central results of this paper — the causal removal of explicit seeds, the rejection of a lower bottom by every-step continuity of the latent region, and the direct measurement of direction-subspace reorganization during the rapid expansion.

---

## References

[1] N. Kihara, Basic Axiom System of the Anonymous Equal-Amplitude Composite Wave Model (purified definition paper), Zenodo. Concept DOI: [10.5281/zenodo.21315735](https://doi.org/10.5281/zenodo.21315735).

[2] N. Kihara, Arrest of Spontaneous Splitting and Creation of Additional Axes in Closed N-Body Relational Wave Systems (Paper 7), Zenodo. Concept DOI: [10.5281/zenodo.21543070](https://doi.org/10.5281/zenodo.21543070).

[3] N. Kihara, The Geometric Identity of the Zero-Square-Sum Constraint and Scale Anonymity — Projectivized Isotropic Cone, Projective Quadric, and Inherent Quantumness (expository note), Zenodo. Concept DOI: [10.5281/zenodo.21495305](https://doi.org/10.5281/zenodo.21495305).

[4] N. Kihara, How Is Localization Exchanged by the Exchange-Interference Scattering Matrix?, Zenodo. Concept DOI: [10.5281/zenodo.21333766](https://doi.org/10.5281/zenodo.21333766).

[5] N. Kihara, Reading the White-Cat/Black-Cat/Gray-Cat Metastable States in a Closed System, Zenodo. Concept DOI: [10.5281/zenodo.21353208](https://doi.org/10.5281/zenodo.21353208).

[6] G. Benettin et al., The Metastable State of Fermi–Pasta–Ulam–Tsingou Models, Entropy **25**(2), 300 (2023).

[7] R. E. Mirollo, S. H. Strogatz, The spectrum of the locked state for the Kuramoto model of coupled oscillators, Physica D **205**, 249 (2005).

[8] I. Daumont, T. Dauxois, M. Peyrard, Modulational instability: first step towards energy localization in nonlinear lattices, Nonlinearity **10**, 617 (1997).

[9] T. Mori, T. N. Ikeda, E. Kaminishi, M. Ueda, Thermalization and prethermalization in isolated quantum systems, J. Phys. B **51**, 112001 (2018).

[10] S.-W. Kim, J. Nishimura, A. Tsuchiya, Expanding (3+1)-dimensional universe from a Lorentzian matrix model for superstring theory in (9+1)-dimensions, Phys. Rev. Lett. **108**, 011601 (2012).
