# Geometric Rapid Expansion Is Specific to Unstable Self-Consistent Closures: Causal Discrimination of Onset Modes via Generic Zero-Closure Initial States

**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-08-05
**Version DOI:** [10.5281/zenodo.21798855](https://doi.org/10.5281/zenodo.21798855)
**Concept DOI:** [10.5281/zenodo.21798854](https://doi.org/10.5281/zenodo.21798854)

---

## Claims

Paper 8 [3] causally removed the two explicit seeds from the time evolution discovered in Paper 7 [2] — long latency, geometric rapid expansion, and spontaneous formation of a metastable three-directional structure — and established that the phenomenon does not depend on the seeds. At the same time, Paper 8 left the causality of the initial state itself unverified, as a stated boundary condition. The experiments had consistently started from the self-consistent circularly polarized eigenmode, i.e., a maximally coherent state in which all relational waves ride on a single eigen-oscillation, and Paper 8 explicitly stated that "a systematic control from generic states satisfying only the closure condition is the next task" (Paper 8 §9.6). The present paper is the execution report of that control.

First, we built a generator that constructs generic zero-closure initial states without contrivance. The raw material is Gaussian white noise, and closure is not imposed afterwards by projection or correction; it holds identically, term by term, through the construction of two orthogonal equal-norm components. The removal of the DC component ($k=0$) and, for even $N$, the Nyquist component ($k=N/2$) is not an additional assumption. The zero square sum of a single mode $k$ is proportional to $\sum_n e^{4\pi ikn/N}$, so a mode with $2k\equiv 0 \pmod N$ can never close on its own. The first axiom already denies these components the qualification to exist; the removal is the enforcement of a theorem (§2).

Second, with the classification criteria fixed before execution, we measured a dichotomy of onset modes (§3–§4). The series starting from self-consistent closures are of the latency-burst type: after a latency of 288–695 steps, $\log f$ grows linearly over about 18 decades (N=5: 0.0494/step, coefficient of determination $R^2=0.9999998$; N=40: 0.0350/step, $R^2=0.998$). The generic zero-closure white harmonic sectors — all 42 of them across N=5 and N=40 — are of the immediate type: zero floor residence, crossing at 11–16 steps, no amplification epoch whatsoever, and an immediate transition to the metastable three-directional structure (rank 4). There is not a single exception.

Third — this is the central claim of the paper. Geometric amplification is not a universal property of arbitrary closed states. It is a phenomenon specific to the instability of the self-consistent closure — the single coherent condensate — and no amplification arises from a sea that already carries harmonic structure. The metastable three-directional structure is a terminal point common to both onset modes; the only fossil that distinguishes the onset mode is the trace of amplification. Inside this model, therefore, a time evolution bearing the trace of 18-decade exponential amplification and an onset from white chaos are incompatible. The onset mode is neither an assumption nor an interpretation; it is a proposition dynamically decidable by the presence or absence of amplification. This proposition is sharpened by the stability measurements of the sixth item into "specific to unstable self-consistent closures" — the title adopts that final form.

Fourth, whether amplification occurs is determined not by the provenance of the state but by its structure alone (§6). The parent state obtained by solving the self-consistency iteration with white noise as its seed (N=40) shows the same latency-burst type as the deterministically constructed control, with matching growth rates (0.0347 vs. 0.0350, under 1% difference). The parent of the same generator at N=5 converged to a different eigen-solution and did not amplify out to 12000 steps. What is established at this stage is therefore the necessity side: if amplification occurs, the initial state is a self-consistent relative equilibrium. This asymmetry is resolved by the stability measurements of the sixth item, which establish that, within the tested range, **an unstable relative equilibrium is the necessary and sufficient condition**. What remains is a theoretical classification predicting, from the structure of a solution, which solutions are stable and which unstable (§8).

Fifth, the dichotomy admits a natural mechanistic reading (§5). The current update rule carries no explicit generation vertex built from amplitude products — a fact mechanically verifiable from the structure of the code. The self-consistent closure is a relative equilibrium of the dynamics; geometric amplification is its unstable departure, and the immediate departure of generic states is the ordinary motion of non-equilibrium states. We state this reading explicitly as the linear no-generation hypothesis and do not claim a proof of a general no-generation theorem. That the white sea behaves as an already-unfolded state from the start, and that only the condensate has room to amplify, are two faces of this reading.

Sixth, we converted this mechanistic reading into measurements through three verification experiments (§5.1–5.3). (i) The one-step residual is $\le2.5\times10^{-15}$ for the 4 parent states and $\ge1.8\times10^{-2}$ for the 82 non-equilibrium states (42 sectors + 40 generic states) — a 13-order-of-magnitude separation with nothing in between. "Self-consistent closure = relative equilibrium" is a direct observable. (ii) The maximal growth rate of the tangent map at the relative equilibria matches the prediction $\lambda_{\max}=\text{rate}_f/2$ exactly at N=5 (0.02468 vs. 0.02468) and within 4% at N=40, and the inert white-origin parent (N=5) turned out to be a stable equilibrium with not a single unstable eigenvalue. Hence, across all 86 tested states, the equivalence **latency exponential amplification ⟺ unstable relative equilibrium** holds without exception. (iii) In the sweep of the perturbation amplitude $\varepsilon$, $t_{\mathrm{cross}}$ is linear in $-\ln\varepsilon$ (N=5: measured/predicted slope ratio 0.996) and saturates at the natural-trajectory value for $\varepsilon\lesssim3\times10^{-13}$. This saturation dynamically measures the location of the intrinsic seed of the seedless trajectory; its value $\varepsilon_{\mathrm{eff}}\approx3.4\times10^{-13}$ (N=5) is of the same order as the parent's eigenmode residual $2.14\times10^{-13}$ — a partial answer to Paper 8's open problem (i), the physical origin of the first minute difference. Furthermore, all 40 generic zero-closure states (direct random samples from the manifold and sector mixtures) are of the immediate type, establishing that the dichotomy is not a peculiarity of single sectors.

Seventh, we record the full process of refutation and correction (§7). The first version of the generator satisfied closure through the cancellation of cross pairs of harmonics $k$ and $N-k$, leaving the self-paired sectors unclosed — a defect. The first closure measurement on the experimental side also used a wrong quadratic form. Both were corrected through control experiments (bit reproduction of the original, correction, re-measurement). The central dichotomy was observed in both the defective first version and the corrected version, surviving the bookkeeping correction. The result is not an artifact of a particular implementation.

Finally, the experiments in this paper were not designed on the basis of prior work. The starting point was an internal question about the naturalness of the initial spectrum; a control was designed independently, a result contradicting the working hypothesis was obtained, and only then was the prior literature re-surveyed — revealing that the standard inflationary initial-condition problem is exactly the point of controversy that this result touches (§10). The process by which a self-consistent closure spontaneously forms out of a white sea (the capture problem) is not observed within the scope of this paper, and the generation of structure from the condensate lies outside the linear dynamics, being treated in the series on inelastic vertices derived from closure conservation [5]. The unresolved part of the onset-mode problem is thereby confined to a single location: the capture from chaos to condensate.

---

## 1. System (Inherited Definitions)

**Axiom 1 (nontrivial zero-square-sum closure)** $\sum_e Z_e^2=0,\ Z\neq0$. **Axiom 2 (finite recurrence)** $U^n=I$. Both axioms follow the foundational axiom system [1].

The dynamical system is identical to Papers 7 [2] and 8 [3]. A complex amplitude $Z_e$ is placed on each edge $e$ of the complete graph $K_N$ ($M=N(N-1)/2$); the state closes on $Z\in\mathbb C^M$, $\lVert Z\rVert=1$, $Z^{\mathsf T}Z=0$. The generator $K$ is a real antisymmetric matrix depending only on edge phases, and the update is a Cayley transform after spectral-norm normalization ($\gamma=\tan(\pi/144)$), preserving the norm exactly. The dynamics code is the Paper 7 original, imported read-only with SHA-256 fixation and never modified (§11).

The observable is the departure fraction from each series' **own initial plane**. From the initial state $v_0$ we form $p=\mathrm{Re}\,v_0/\lVert\mathrm{Re}\,v_0\rVert$ and $q$ ($\mathrm{Im}\,v_0$ orthogonalized against $p$), and define

$$f(t)=\frac{\lVert Z-\Pi_{\mathrm{span}\{p,q\}}Z\rVert^2}{\lVert Z\rVert^2}$$

By construction $f(0)=0$ for every series, so onset modes are compared in the same coordinate. Crossing is the first step at which $f>0.05$, the pre-existing threshold of Paper 7.

## 2. Construction of Generic Zero-Closure Initial States

### 2.1 The generator (white, zero-closure, post-hoc harmonic readout)

The only theoretical input is $N$. Each relational wave $w_m$ is an $N$-sample complex waveform built from the two orthonormal columns $(q_1,q_2)$ obtained by QR decomposition of $N\times2$ Gaussian white noise:

$$w_m=\frac{q_1+iq_2}{\sqrt2}$$

The zero square sum of each row holds identically:

$$\sum_n w_m[n]^2=\lVert q_1\rVert^2-\lVert q_2\rVert^2+2i\,q_1\!\cdot\!q_2=0$$

Closure is not manufactured by corrections; it emerges structurally from two orthogonal equal-norm components. Harmonics are not placed by the generator; they are read out after saving via an N-point DFT.

### 2.2 DC/Nyquist removal is a theorem, not an assumption

The zero square sum of a single mode $k$ is proportional to $\sum_n e^{4\pi ikn/N}$: it vanishes if $2k\not\equiv0\pmod N$ and does not vanish if $2k\equiv0\pmod N$. That is, DC ($k=0$) and, for even $N$, Nyquist ($k=N/2$) can never close alone for any $N$, and under the first axiom they lack the qualification to exist. The corrected generator (v3) projects the DC direction (the all-ones vector) and, for even $N$, the Nyquist direction (the alternating $\pm1$ vector) out of each noise column before QR. Since the projection preserves Gaussianity, this is a rigorous implementation of "white noise on the admissible modes." The strengthened version (v4) further makes the closure debts between $|k|$ pairs (pair-product debt) identically zero, term by term, through direct spectral construction on admissible signed modes. The dynamical experiments of this paper use v3; v4 was used to audit the classification of initial states.

### 2.3 Injected series

Three kinds of initial state are injected into the same dynamics with the same measurement.

1. **control**: the self-consistent circularly polarized eigenmode returned by `make_parent` of Papers 7–8 (deterministic construction).
2. **white-origin parent**: the eigenmode obtained by solving the self-consistency iteration with white noise as its initial value (the v3 generator's `parent_vector`). Its structure — self-consistent closure — is the same as the control; only its provenance (random-noise origin) differs.
3. **white harmonic sectors**: states $Z_k$ obtained by extracting only the harmonic-$k$ content of the white sea as $M$ relational amplitudes (all $k$ with $2k\not\equiv0\pmod N$; 4 for N=5, 38 for N=40, 42 in total). The field closure of each sector is exactly zero (measured: 0.0 for all 42). These are the "generic states satisfying only the closure condition, with no eigenmode condition imposed."

## 3. Classification Criteria (Fixed Before Execution)

The following were fixed before execution, with post-hoc changes prohibited.

- **Latency** $t_{\mathrm{launch}}=\max\{t: f(t)<10^{-20}\}$ (the last time on the floor; 0 if $f\ge10^{-20}$ from the start).
- **Geometric criterion**: least-squares fit of $\ln f$ vs. $t$ in the window $f\in[10^{-20},10^{-2}]$; the window must span at least 6 decades with $R^2>0.99$.
- **Immediate criterion**: $f(2)>10^{-6}$ (finite departure within 2 steps).
- **Classes**: latency-burst type = $t_{\mathrm{launch}}>100$ and geometric / immediate type = $f(2)>10^{-6}$ / anything satisfying neither is recorded as "other" together with its raw data.

The measurement records $f(t)$ at every step ($T=6000$). The only random numbers are the generator seeds (N=5: 2, N=40: 1) and the auxiliary injection seeds (parent 91000, sector $92000+k$), all recorded.

## 4. Result: the Dichotomy of Onset Modes

*Table 1: Classification of all series (canonical sources: `paper8_em9r_profile_N00005_v1.json` / `N00040_v1.json`)*

| Series | Class | Latency | Crossing | Decades | Rate /step | $R^2$ | $f(2)$ |
|---|---|--:|--:|--:|--:|--:|--:|
| control N=5 | latency burst | 288 | 1166 | 17.97 | 0.0494 | 0.9999998 | $2.5\times10^{-29}$ |
| control N=40 | latency burst | 468 | 2011 | 17.99 | 0.0350 | 0.9981 | $2.0\times10^{-30}$ |
| white-origin parent N=40 | latency burst | 695 | 2167 | 17.98 | 0.0347 | 0.9997 | $5.2\times10^{-32}$ |
| white-origin parent N=5 | inert (other) | — (stays on floor) | none ($t\le12000$) | 0 | — | — | $3.0\times10^{-32}$ |
| white sectors N=5 (4) | immediate 4/4 | 0 | 11–16 | — | — | — | $0.95$–$1.6\times10^{-3}$ |
| white sectors N=40 (38) | immediate 38/38 | 0 | 13–14 | — | — | — | $1.25$–$1.46\times10^{-3}$ |

![N=5 rise profiles](figs/fig_em9r_profile_N00005.png)

*Figure 1: Rise profiles at N=5 (semilogy). Gray = control, black = white-origin parent, blue = the 4 white sectors. Only the control shows latency followed by 18 decades of linear growth (geometric amplification). The sectors depart immediately with zero latency.*

![N=40 rise profiles](figs/fig_em9r_profile_N00040.png)

*Figure 2: Rise profiles at N=40. The white-origin parent (black) shows a latency burst of the same form as the control (gray), with growth rates matching to under 1%. All 38 sectors (blue) are immediate.*

Three points stand out.

1. **The presence or absence of an amplification epoch dichotomizes.** The three latency-burst series ride a single linear law of $\log f$ over about 18 decades (best $R^2=0.9999998$). The 42 immediate series contain no interval that could be called exponential growth: zero floor residence, finite departure to $f\sim10^{-3}$ within 2 steps, crossing at 11–16 steps. Not a single intermediate type appeared.
2. **The terminal point is shared.** All 42 white sectors transition after crossing to the metastable three-directional structure (rank 4; late-time $f$: 0.73–0.83 for N=5, 0.995–0.999 for N=40; canonical source: `paper8_em9r_v3_result_*.json`). The three-directional structure of Papers 7–8 is not specific to coherent onsets. What is specific is **the road there** — the trace of amplification.
3. **The burst growth rate does not depend on preparation.** The N=40 growth rates of the control (deterministic construction) and the white-origin parent (random origin) agree at 0.0350 vs. 0.0347. The amplification rate is a property of the dynamics and the equilibrium structure, not of the state's provenance.

## 5. Mechanism: Amplification as Instability of a Relative Equilibrium

**Absence of an explicit generation vertex (a mechanically verifiable fact).** One step of the current engine is a Cayley-type unitary transform by the generator $K(\theta)$ frozen at the current phase configuration $\theta=\arg Z$, and carries no explicit generation vertex built from amplitude products of $Z$. This is verifiable from the structure of the code.

**Linear no-generation hypothesis (interpretation).** Each step is linear with $\theta$ frozen, but the full map $F(Z)=U(\arg Z)\,Z$ is state-dependent and hence nonlinear, so that iteration creates no new content does not follow from the above fact alone. This paper interprets the observed geometric amplification not as generation of new content but as unstable departure from a relative equilibrium. This is a working hypothesis explaining the dichotomy of this paper; no proof of a general no-generation theorem is claimed. A proof of an invariant subspace, or a direct derivation of the update equations of the harmonic coefficients, is the task of the next stage.

Under this hypothesis the dichotomy reads as follows.

- The self-consistent closure is a **relative equilibrium** of the dynamics (fixed up to a uniform phase rotation). On the equilibrium, $f$ does not strictly grow; departure occurs only through the linear instability of the equilibrium growing exponentially from the numerical floor (eigenmode residual $\sim10^{-13}$ plus rounding). This is the identity of the latency and of the ensuing 18 decades of linear $\log f$.
- Generic zero-closure states are not equilibria. One step of the map immediately imparts a finite-angle rotation, so there is neither floor residence nor an exponential interval, and the state leaves its initial plane in $O(10)$ steps. **The amplification epoch is absent because the "small departure from equilibrium" that would have to grow does not exist in the first place.**

Thus the measurement "geometric amplification is specific to the coherent condensate" is not an accidental contrast; it is naturally explained by the difference between relative equilibria and non-equilibrium initial states. From the same viewpoint, the series of negative results up to Paper 8 (no internal harmonic generation, no precursor of statistical distinction) is consistent with the absence of an explicit generation vertex — although deriving them rigorously as corollaries requires proving the hypothesis above.

This reading generates three verifiable predictions. (1) The one-step residual $r(Z_0)=\min_\phi\lVert F(Z_0)-e^{i\phi}Z_0\rVert$ is numerically zero for parent states, finite for non-equilibrium states, with nothing in between. (2) The logarithm of the largest tangent-map eigenvalue $\lambda_{\max}$ at a relative equilibrium equals half the measured burst rate, $\text{rate}_f/2$, since $f$ is an amplitude-squared quantity; the inert parent has no unstable eigenvalue. (3) Injecting a perturbation amplitude $\varepsilon$, the growth time from $f(0)\approx\varepsilon^2$ to $f_{\mathrm{cross}}$ obeys the logarithmic law $t_{\mathrm{cross}}\approx\text{const}+(-2\ln\varepsilon)/\text{rate}_f$. The following three experiments (ONS-1/2/3, criteria and predictions fixed before execution at the head of each script) tested these.

### 5.1 Verification 1: One-Step Residual and Tangent-Map Spectrum (ONS-2)

*Table 2: One-step residuals and tangent spectra (canonical source: `onset_equilibrium_residual_N*.json`)*

| State | Residual $r$ | $\lambda_{\max}$ measured | $\lambda_{\max}$ predicted | Unstable eigenvalues |
|---|--:|--:|--:|--:|
| control parent N=5 | $2.5\times10^{-15}$ | **0.02468** | 0.02468 | 4 |
| control parent N=40 | $7.1\times10^{-16}$ | 0.01824 | 0.01750 | 383 |
| white-origin parent N=40 | $1.0\times10^{-16}$ | 0.01815 | 0.01734 | 389 |
| white-origin parent N=5 (inert) | $1.6\times10^{-16}$ | **0.00000** | $\le10^{-3}$ | **0** |
| white sectors N=5 (4) | $2.1$–$2.6\times10^{-2}$ | — | — | — |
| white sectors N=40 (38) | $1.8$–$2.1\times10^{-2}$ | — | — | — |

Predictions (1) and (2) hold. The residual separates the 4 equilibria from the non-equilibrium states by 13 orders of magnitude, with not one intermediate. $\lambda_{\max}$ matches the prediction exactly at N=5 and runs 4–5% high at N=40 — consistent with the N=40 window fit ($R^2=0.998$) including a nonlinear-correction regime. The ratio of the two N=40 equilibria's $\lambda_{\max}$, 1.005, matches the measured burst-rate ratio 1.009. The decisive entry is the last equilibrium row: **the inert white-origin parent is a stable relative equilibrium (zero unstable eigenvalues).** This resolves the asymmetry of §6 (coherence necessary but not sufficient), and across all 86 tested states (4 equilibria + 42 sectors + 40 generic states),

$$\boxed{\text{latency exponential amplification}\iff\text{unstable relative equilibrium}}$$

holds without exception. The equivalence can be written as a trichotomy: non-equilibrium states depart immediately by a finite angle, stable relative equilibria never depart, and only unstable relative equilibria depart exponentially after a latency. What determines the mode of time evolution is neither the provenance of the state nor its apparent coherence, but the dynamical orbit class to which the initial state belongs.

### 5.2 Verification 2: Logarithmic Law of the Perturbation Amplitude and Measurement of the Intrinsic Floor (ONS-1)

A random perturbation $\varepsilon\eta_\perp$ orthogonal to the parent plane ($\eta$ over 3 seeds, $\varepsilon=10^{-4}$–$10^{-14}$) was injected onto the control parent. $t_{\mathrm{cross}}$ is linear in $-\ln\varepsilon$ (canonical source: `onset_eps_sweep_N*.json`); over $\varepsilon\ge10^{-12}$ the slope is 40.35 at N=5 (prediction $2/\text{rate}_f=40.52$, ratio 0.996, $R^2=0.9875$) and 65.45 at N=40 (prediction 57.16, ratio 1.145, $R^2=0.9985$). The N=40 excess is consistent with the initial growth of a random perturbation being compositionally slower on an equilibrium carrying 383 unstable modes (growth rates spread over 0.0160–0.0182), and the measured per-run window rates, which rise monotonically from 0.027 toward 0.035 as $\varepsilon$ decreases, corroborate this.

At $\varepsilon=10^{-14}$, $t_{\mathrm{cross}}$ saturates for all seeds at the value of the seedless natural trajectory (N=5: 1166, N=40: 2011) — evidence that the injected perturbation fell below the system's intrinsic seed. Extrapolating the line, the intrinsic floor is measured as $\varepsilon_{\mathrm{eff}}=3.4\times10^{-13}$ (N=5) and $2.5\times10^{-13}$ (N=40), of the same order as the parent's eigenmode residual $2.14\times10^{-13}$ (Paper 8 §9.1). That is, **the "first minute difference" that set the latency time of the seedless trajectory is the eigenmode residual of the parent construction itself** — a partial answer, by dynamical measurement, to Paper 8's open problem (i). Both the pre-registered pooled regression including the saturation point ($\varepsilon=10^{-14}$) and the post-hoc analysis excluding the saturation regime are recorded side by side in the canonical source (`onset_eps_sweep_analysis_v1.json`).

### 5.3 Verification 3: Direct Injection of Generic Zero-Closure States (ONS-3)

The white sectors of §4 are 42 generator-derived states, and since the dynamics is nonlinear, the behavior of other points on the closure manifold does not follow from their results. Two ensembles were therefore injected directly (10 states each per N; canonical source: `onset_general_closed_N*.json`).

- **A: direct random samples from the manifold** — $Z=(u_1+iu_2)/\sqrt2$ from the two QR-orthonormal columns of $M\times2$ Gaussian white noise. Closure holds identically (measured $\le5.6\times10^{-16}$).
- **B: random mixtures of admissible sectors** — the 42 sectors superposed with complex Gaussian coefficients, then projected onto the closure manifold by minimal change (SVD singular-value equalization of $[X\,Y]$; projection distances 0.003–0.16 recorded; post-projection closure $\le3.9\times10^{-16}$).

Result: **all 40 states are of the immediate type** (crossing 11–21), with no latency-burst and no "other." As predicted. The dichotomy is not a peculiarity of the single-harmonic-sector construction; it is a consequence of the state structure — equilibrium or not.

## 6. Provenance Independence

The white-origin parent results show that the boundary of the dichotomy is not "random origin vs. deterministic construction." Starting from the same white noise, if the self-consistency iteration converges to an eigenmode (N=40), that state shows the same latency burst at the same rate as the deterministic control. At N=5, meanwhile, the iteration converged to an eigen-solution of a different family and did not amplify out to 12000 steps. That is:

- The first boundary determining the onset mode is whether the state is a **self-consistent relative equilibrium**, and what finally determines the presence of amplification is **the stability of that relative equilibrium** (§5.1). How the state was prepared plays no role at either boundary.
- Among self-consistent closures there are amplifying solutions (unstable equilibria) and inert solutions (stable equilibria). Coherence alone is a **necessary but not sufficient condition** for amplification.

This asymmetry does not weaken the claim of this paper, whose implication runs in a fixed direction: if amplification was observed, the initial state was a self-consistent closure. That there is not one counterexample among the 82 non-equilibrium states (42 sectors + 40 generic states, §5.3) is the measured basis of this implication. Moreover, the tangent-spectrum measurement of §5.1 resolved the asymmetry itself: the inert white-origin parent (N=5) is a stable relative equilibrium (zero unstable eigenvalues), and what finally separates amplification from its absence is whether the state is an equilibrium and, if so, its stability.

## 7. Record of Refutation and Correction

This result was not obtained in a single experiment. The full process is recorded.

1. **Defect of the generator's first version (v2).** Writing the row closure in DFT form, $\sum_n w^2=N\sum_k c_kc_{N-k}$: v2 satisfied closure through the **cancellation of cross pairs** of harmonics $k$ and $N-k$, leaving the self-paired sectors (DC, Nyquist) unclosed. Based on the theorem of §2.2, v3 projects them out at the source of the noise. v3 consumes random numbers in the same order as v2; the parent vector and raw noise are bitwise identical.
2. **Measurement error on the experimental side.** The first E-M9r run used a wrong quadratic form for the sectors' closure measurement and misreported "white sectors are unclosed." Under the correct field closure, the 42 white sectors are legitimate, exactly zero-closed states (conservation verified: drift $\sim10^{-16}$ over 500 steps).
3. **Correction via control experiments.** The generator was copied into a managed area with SHA-256 identity, and the erroneous results were bit-reproduced with identical seeds before being corrected. We do not claim a correction without reproducing the error.
4. **Survival of the dichotomy.** The central dichotomy (latency burst vs. immediate) was observed in both v2 (defective) and v3 (corrected). The bookkeeping error changed the handling of the self-paired sectors but did not change the dichotomy itself.

## 8. Boundary Conditions (What This Paper Does Not Prove)

1. **The theoretical classification of stability is incomplete.** That amplification is decided by the stability of the equilibrium was measured in §5.1 (inert parent = stable equilibrium, λ_max=0). But no theoretical classification is given that predicts, from the structure of a solution (N, family, spectrum), which self-consistent families are stable and which unstable. Why the N=5 white-origin solution converged to a stable one and the N=40 solution to an unstable one remains open.
2. **Seed count.** One generator seed per $N$ (N=5: 2, N=40: 1). On the non-equilibrium side, however, there was no exception among the 42 independent sectors within the same seed plus the 40 generic zero-closure states from independent random streams (§5.3, rng 94000+/95000+). A sweep over generator seeds has not been performed.
3. **Time range.** The profile measurement uses $T=6000$; the v3 main run $T=12000$. The "absence of amplification" of the white sectors is definitive as a finite departure within the first 2 steps, but a second amplification beyond the observation window is not logically excluded.
4. **Range of $N$.** Only N=5 and 40. The N=300 of Paper 8's main series was not run in this control.
5. **The capture problem.** The spontaneous formation of a self-consistent closure out of the white sea is not observed in the linear dynamics of this paper. This is the province of the inelastic generation channels from Paper 9 [5] onward.
6. **Limitation of terminology.** "Inflation" and "amplification" in this paper are phenomenon names of this series (the time structure latency → geometric → metastable of $f$), not cosmological inflation (exponential expansion of the metric) itself. The correspondence of §10 is a structural analogy, not a dynamical derivation in field theory or general relativity.
7. **The linear no-generation hypothesis is unproven.** Of the mechanistic reading of §5, the absence of an explicit generation vertex is a verified fact, and relative equilibrium, instability, the logarithmic law, and the intrinsic floor became measurements in §5.1–5.3. But that iteration of the state-dependent unitary cannot create new harmonic content (preservation of the invariant subspace corresponding to the initial spectral support) remains unproven. The central results — the dichotomy of §4 and the equivalence of §5.1 — are measurements and do not depend on this unproven part.

## 9. Conclusion

We executed the systematic control left by Paper 8 — starting from generic states satisfying only the closure condition. The result is a complete dichotomy of onset modes. Only series starting from self-consistent closures (single coherent condensates) show latency and 18 decades of geometric amplification. The 82 non-equilibrium zero-closure states (42 white sectors + 40 generic states) have, without exception, no amplification epoch and move immediately, by finite-angle rotation, to the metastable three-directional structure. The presence of amplification does not depend on the provenance of the state; it is decided by the structure of the state alone.

The mechanism closed with measurements. The one-step residual separates equilibria from non-equilibria by 13 orders of magnitude; the maximal tangent-map growth rate equals half the measured burst rate (exactly at N=5, within 4% at N=40); the inert parent is a stable equilibrium (zero unstable eigenvalues); the logarithmic law of the perturbation amplitude matches the prediction down to its slope; and its saturation measures the intrinsic floor of the seedless trajectory at the location of the eigenmode residual. Across all 86 tested states, the equivalence — latency exponential amplification ⟺ unstable relative equilibrium — holds without exception.

Together with Papers 7 and 8, the causal ledger closes in the following form. The explicit seeds are not the cause of the phenomenon (Paper 8). The cause is that the initial state belongs to an unstable self-consistent relative equilibrium (this paper). The necessary and sufficient condition for latency exponential amplification is, within the tested range, an unstable self-consistent relative equilibrium. The metastable three-directional structure is a terminal point independent of the onset mode, and only the trace of geometric amplification remains as the fossil showing that the system began from an unstable relative equilibrium — a coherent condensate.

What remains unresolved: (i) the theoretical classification of the correspondence between self-consistent solution families and stability (why the N=5 white-origin solution is stable and the N=40 solution unstable); (ii) the capture from white chaos to condensate; (iii) the inelastic generation of structure (matter) from the condensate — (iii) is already underway in the series deriving the vertex from closure conservation [5].

## 10. Relation to Prior Work (via Re-Survey After the Design)

The control experiments of this paper were not designed as replication or refutation of prior work. The starting point was a question internal to this series — "would not a flat spectrum with all harmonics mixed be a more natural initial state than concentration on a single wavelength at the minimal resolution?" — and the working hypothesis at design time was "the white sea should show the same geometric evolution." The result is the rejection of that working hypothesis (§4). Having obtained a result contrary to expectation, we re-surveyed the literature on initial conditions in inflationary cosmology, and found the following.

**The standard theory does not predict the initial state; it assumes it.** The starting point of slow-roll inflation [6,7] is a nearly homogeneous, potential-dominated inflaton condensate with Bunch–Davies vacuum fluctuations on top. The origin of this initial state is undecided — the "initial conditions problem of inflation" [8,9] — and the observed nearly flat spectrum ($n_s=0.9649$ [14]) is an **output** of inflation, not an initial condition. That is, the standard theory also posits, without deriving it, the onset mode of a coherent condensate. The measurements of this paper show, dynamically and inside this model, that only the state corresponding to that assumption (the self-consistent closure) amplifies.

**The claim that "inflation ignites even from inhomogeneous initial conditions" is a model-dependent numerical result of a few groups.** Numerical-relativity studies [10,11] reported that in large-field (plateau-type) models inflation begins even from gradient-energy-dominated inhomogeneous initial conditions. But the initial data of these computations restrict the field values to the potential plateau before adding inhomogeneity — in the language of this paper, they correspond not to a "white sea" but to a "condensate plus perturbations." The coherence assumption is not removed; it is relocated into the construction of the initial data. If this reading is correct, those results do not contradict the dichotomy of this paper (no amplification from generic states).

**A decision criterion for the critics' intuition.** Against the criticism that inflation demands initial conditions more special than the problem it is supposed to solve [12,13], and the resulting standoff (called a "schism" by the participants themselves [13]), the results of this paper provide a dynamical criterion of the following form. First, as a result inside the model: trajectories with the exponential-amplification time structure did not arise directly from white non-equilibrium states, but only from unstable self-consistent relative equilibria (§4–§5, no exception among 86 states). Translated, as an analogy, into cosmological language: an onset state bearing the trace of exponential amplification cannot be white chaos and must be an unstable coherent condensate — with the proviso that the range of this translation extends only to the structural analogy of this model, and any transfer into field theory is subject to the limitation of §8.6.

The order of emphasis is this. This paper did not draw its hypotheses from these references. An independently designed control returned an unexpected dichotomy, and the post-hoc re-survey revealed that this dichotomy amounts to a decision experiment, within this model, on a point of controversy that has remained open for forty years — is coherence necessary for the onset? Although the model contains no field theory, it arrived at a question of the same structure as the initial-conditions problem and returned a definite answer there (coherence is necessary). This convergence is itself a verifiable outcome of the method of this series: fix the axioms minimally and let the dynamics answer.

## 11. Reproducibility

- The dynamics code is the Paper 7 original (`run_preliminary_seed_ablation_v1.py`), imported read-only (SHA-256: `75a10a5b951302be…`). Generator v3 SHA-256: `d3217579ab54ce29…`. Both hashes are recorded in the `imports` of every result JSON.
- Canonical sources for the generator, classifier, and audits: `make_parent_white_managed_v1/` (v2 control copy, v3, v4, unit tests, classifier control). v2 and v3 are bitwise identical in parent vector and raw noise at identical seeds. The control experiments include SHA-256 identity with the Codex original and bit reproduction of the erroneous results.
- Canonical sources for the experiments: `paper8_em9r_white_harmonics_inflation_v1/` — all JSONs and figures of the first run (v1), the bookkeeping correction (fix_v2), the v3-generator run (v3_result), and the profile measurement (profile; canonical source of Table 1 and Figures 1–2).
- Canonical sources for the mechanism verification: `paper8_onset_mechanism_v1/` — ONS-1 (ε sweep: `onset_eps_sweep_N*.json`, post-hoc analysis `onset_eps_sweep_analysis_v1.json`), ONS-2 (residual and tangent spectrum: `onset_equilibrium_residual_N*.json`, canonical source of Table 2), ONS-3 (generic-state injection: `onset_general_closed_N*.json`) — all scripts, JSONs, and figures. Criteria and predictions are fixed before execution at the head of each script. Random seeds: η=93000+j, direct samples=94000+j, mixtures=95000+j, injection auxiliary=96000+series.
- Commit series: `52c83f13` (initial-state classification) → `aeaf8e57` (first E-M9r) → `d9346e37` (bookkeeping correction) → `e5232885` (v3 re-run) → `94aa76f9` (profile measurement) → `b1cded81` (v4 generator, catalog closure audit) → `e6681744` (discussion note) → `2f443467` (mechanism verification ONS-1/2/3).
- The classification criteria (§3) are stated as pre-registered at the head of the measurement programs and saved in the `criteria` of the result JSONs. All random seeds are recorded (generator: N=5→2, N=40→1; injection auxiliary: parent 91000, sector $92000+k$).
- Environment: Python 3.9.6, NumPy 2.0.2, macOS arm64.

## References

[1] N. Kihara, Foundational Axiom System of the Anonymous Equal-Amplitude Composite Wave Model (Purified Definition Paper), Zenodo. Concept DOI: [10.5281/zenodo.21315735](https://doi.org/10.5281/zenodo.21315735).

[2] N. Kihara, Arrest of Spontaneous Splitting and Creation of Additional Axes in Closed N-Body Relational Wave Systems (Paper 7), Zenodo. Concept DOI: [10.5281/zenodo.21543070](https://doi.org/10.5281/zenodo.21543070).

[3] N. Kihara, Causal Separation of the Time Structure of Three-Direction Formation in Closed N-Body Relational Wave Systems via Two-Stage Seed Removal (Paper 8), Zenodo. Concept DOI: [10.5281/zenodo.21614402](https://doi.org/10.5281/zenodo.21614402).

[4] N. Kihara, The Geometric Identity of the Zero-Square-Sum Constraint and Scale Invariance (Explanatory Note), Zenodo. Concept DOI: [10.5281/zenodo.21495305](https://doi.org/10.5281/zenodo.21495305).

[5] N. Kihara, The Generation Structure of Fermions (Paper 9), Zenodo. Concept DOI: [10.5281/zenodo.21766706](https://doi.org/10.5281/zenodo.21766706).

[6] A. H. Guth, Inflationary universe: A possible solution to the horizon and flatness problems, Phys. Rev. D **23**, 347 (1981).

[7] A. D. Linde, Chaotic inflation, Phys. Lett. B **129**, 177 (1983).

[8] R. Brandenberger, Initial conditions for inflation — A short review, Int. J. Mod. Phys. D **26**, 1740002 (2017). [arXiv:1601.01918](https://arxiv.org/abs/1601.01918).

[9] A. Linde, On the problem of initial conditions for inflation, Found. Phys. **48**, 1246 (2018). [arXiv:1710.04278](https://arxiv.org/abs/1710.04278).

[10] W. E. East, M. Kleban, A. Linde, L. Senatore, Beginning inflation in an inhomogeneous universe, JCAP **09** (2016) 010. [arXiv:1511.05143](https://arxiv.org/abs/1511.05143).

[11] K. Clough, E. A. Lim, B. S. DiNunno, W. Fischler, R. Flauger, S. Paban, Robustness of inflation to inhomogeneous initial conditions, JCAP **09** (2017) 025. [arXiv:1608.04408](https://arxiv.org/abs/1608.04408).

[12] A. Ijjas, P. J. Steinhardt, A. Loeb, Inflationary paradigm in trouble after Planck2013, Phys. Lett. B **723**, 261 (2013). [arXiv:1304.2785](https://arxiv.org/abs/1304.2785).

[13] A. Ijjas, P. J. Steinhardt, A. Loeb, Inflationary schism, Phys. Lett. B **736**, 142 (2014). [arXiv:1402.6980](https://arxiv.org/abs/1402.6980).

[14] Planck Collaboration, Planck 2018 results. X. Constraints on inflation, Astron. Astrophys. **641**, A10 (2020). [arXiv:1807.06211](https://arxiv.org/abs/1807.06211).
