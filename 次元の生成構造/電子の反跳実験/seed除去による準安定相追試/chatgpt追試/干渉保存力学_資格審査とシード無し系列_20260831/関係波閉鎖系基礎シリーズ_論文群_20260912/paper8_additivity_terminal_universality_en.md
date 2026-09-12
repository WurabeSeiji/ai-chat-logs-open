# The Additivity Law of the Terminal Equal Amplitude $a_\infty=\sqrt{H/M}$ —— Separating Terminal Universality from Ignition Dynamics ($N=9$ alone vs $N=7+N=6$, 3 initial-value families)

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 8 of 10)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22729116
**Concept DOI:** 10.5281/zenodo.22729115
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22729116

---

## Abstract

The post-lock terminal state is equal-amplitude (Paper 5). We verify the **additivity law** $a_\infty=\sqrt{H/M}$, in which its radius is determined solely by the extensive variables $(H,M)$, using the mixed sum of triangular numbers $M(9)=M(7)+M(6)$ ($36=21+15$). The mixed system of $N=7$ and $N=6$, normalized to the density $H_i/M_i=1/36$, agrees in terminal ring radius with the standalone $N=9$ system (both giving $a_\infty=\sqrt{1/36}=1/6=0.1666667$). On the other hand, onset, path, and plateau are imprinted by the composition ($N=9$ alone onset $191$ vs mixed union departs at the earliest stage). Furthermore, changing the initial-value family (staged $2\pi/N$ increments / 90-degree analytic floor / 90-degree symmetric v2) changes **the presence or absence of ignition itself** ($N=9$: staged onset $191$ / 90-degree analytic floor onset $140$ / 90-degree v2 is the non-igniting exact fixed point $\sigma=N-1$), yet the terminal value $1/6$ is common across all three families. That is, **terminal universality ($a_\infty=\sqrt{H/M}$) is independent of the symmetry family of the initial values**, whereas **ignition dynamics depends strongly on it**. The composition (which combination of $N$) is imprinted only during the inflation period, and at the terminal state it is indistinguishable from a single system. $(H,M)$ are extensive quantities, $\rho=H/M$ is an intensive quantity, and $a_\infty=\sqrt\rho$ is an intensive terminal state quantity; the formula itself is an exact consequence of equal amplitude plus $H$ conservation, and what is nontrivial is terminal universality. **The transient remembers the composition (the $J_0,\prod_tJ_t$ of Paper 6), while the terminal state forgets the composition in its amplitude via $|z_e|^2\to H/M$, whereas the phase configuration and the orientation of the face are degenerate and selected by the seed (Paper 5)** —— the amplitude becomes a state quantity, the orientation retains the SSB history. Factorization itself is dynamically impossible because $K_N$ is maximally connected (crossing edges necessarily remain), and the terminal agreement is a consequence of equipartition unrelated to decomposability. Classification: numerical verification (additivity law) + mathematical fact (impossibility of decomposition).

---

## 1. Problem

If the terminal state is equal-amplitude, then $a_\infty=\sqrt{H/M}$ **follows exactly** from $H=\sum_{e=1}^M|z_e|^2=Ma_\infty^2$ —— the formula itself is not nontrivial. The nontrivial part that this paper confirms by numerical experiment is **"terminal universality," namely that even different compositions and different initial-value families reach the same equal-amplitude law at the terminal state**. Generalization: if the subsystems have identical density $\rho_i=H_i/M_i=\rho$, then since $H_{\rm tot}=\sum_iH_i=\rho\sum_iM_i=\rho M_{\rm tot}$,

$$\frac{H_{\rm tot}}{M_{\rm tot}}=\rho,\qquad a_{\infty,\rm tot}=\sqrt{H_{\rm tot}/M_{\rm tot}}=\sqrt\rho=a_{\infty,i}.$$

That is, $(H,M)$ are **extensive quantities**, $\rho=H/M$ is an **intensive quantity**, and $a_\infty=\sqrt\rho$ is an **intensive terminal state quantity**; the additivity of $(H,M)$ gives the composition-invariance of $\rho$ ($N=9=7+6$ is one instance). We verify this with a heterogeneous mixture and separate terminal universality from ignition dynamics (the imprint of composition).

## 2. Dynamics and Design (faithful copy + density normalization)

The map is the canonical one_step (SHA `1abf2353…`). The wrapper is a faithful copy of the canonical implementation with only minimal changes to the initial-value source and output names (diff audited, physics formulas unchanged). Using the triangular-number decomposition $M(9)=36=21+15=M(7)+M(6)$, the theoretical-floor parents of $N=7,N=6$ (each with $H=1$) are normalized to the density-matching condition $H_i/M_i=1/36$ ($c_7=\sqrt{21/36}=0.7638$, $c_6=\sqrt{15/36}=0.6455$, total $H=1$). den=N・500 step. The mixed system runs each subsystem independently with its own adjacency and its own den.

## 3. Establishment of the Additivity Law (Experiment 21)

- **The terminal ring radius agrees**: $N=9$ alone $0.1666667$ (min=max), the $N=7$ part of the mixture $0.1666667$, the $N=6$ part $0.16649$–$0.16677$, all equal to the theoretical value $\sqrt{1/36}=0.1666667$.
- **The ignition dynamics is entirely different by composition** (only the terminal state is universal): onset$(>0.05)$ is $191$ for $N=9$ alone vs the mixed union departs at the earliest stage (the union is meaningless, so it is evaluated per subsystem). Per-subsystem onset: $N=7$ part $145$, $N=6$ part $153$ (agreeing with the $153$ of the $N=6\times3$ series = scale-invariant).
- $H$ conservation and closure maintenance are sound ($H$ relative difference $\sim10^{-14}$, maximum closure $\sim10^{-14}$, finite at all steps).

Whether by homogeneous $m$-fold or heterogeneous mixture, the terminal state is indistinguishable from a single system via $a_\infty=\sqrt{H/M}$. The structure (which combination of $N$) is imprinted only on onset, path, and plateau.

## 4. Change of Ignition Dynamics by Initial-Value Family (Experiments 23 & 25) —— comparison of 3 families

The same experiment was re-run changing only the initial-value family:

| $N=9$ | staged ($2\pi/N$ increments) | 90-degree analytic floor (integer K) | 90-degree symmetric v2 (0/90 equal amplitude) |
|---|---|---|---|
| ignition | ignition (onset $191$) | ignition (onset $140$) | **non-ignition** (the $\sigma=N-1$ exact fixed point of $N\equiv1\,\mathrm{mod}\,4$) |
| mixed-part onset ($N=7,N=6$) | $145,\ 153$ | $113,\ 121$ | $5,\ 5$ |
| terminal $a_\infty$ | $1/6$ | $1/6$ | $1/6$ |

**The terminal value $1/6$ is common across all three families**, and only the presence/absence of ignition and the onset change by family. For 90-degree v2, the $N=9$ case is an exact fixed point where the $0^\circ$ family becomes a $(N-1)/2$-regular graph ($\sigma=N-1$), which sits at rest on the equal-amplitude ring and does not ignite (consistent with the qualification screening of Paper 2: a neutral fixed point does not ignite). The 90-degree analytic floor (Paper 3) and staged ignite on an unstable floor.

## 5. Factorization Is Impossible; the Terminal Agreement Is a Consequence of Equipartition (Experiment 22)

The complete enumeration of the triangular-number decomposition $M(N)=\sum M(k)$ (equal partition, 27 two-part cases, 46 three-part cases; the "prime" $N$ having no decomposition with parts $M(4)$ or larger are $3,4,5,6,8$) is the basis for the selection, but **dynamical factorization is impossible in principle**: since $K_N$ is a complete graph, crossing edges necessarily remain and it is maximally connected, and $K_{10}\to3K_6$ violates the degree condition. The inflation mismatch of Experiment 21 (onset $191$ vs $145/153$) is the experimental evidence of dynamical impossibility of decomposition. Therefore, **the terminal agreement is a consequence of equipartition ($a_\infty=\sqrt{H/M}$) unrelated to decomposability**. The hierarchization of a condensate = a single vertex (mass = the $H$ privately held by the closure one level below) is presented as a working hypothesis, and its derivation is left as an open problem for future study.

## 6. Memory and Erasure —— the Transient Remembers the Composition and the Terminal Amplitude Forgets It (unification of Papers 5, 6, 8)

The carrier of information separates between the transient and the terminal state. The **transient** (inflation period) remembers the composition and initial-value family: the initial structure is imprinted in the floor Jacobian $J_0$ and the tangent cocycle $\prod_t J(z_t)$ (Paper 6), determining onset, path, and plateau (§3–4: $N=9$ alone onset $191$ vs mixed $145/153$, and even the presence/absence of ignition changes with staged/analytic floor/v2). The **terminal state** has $|z_e|^2\to H/M$ after nonlinear saturation, and the composition information disappears from this amplitude ($a_\infty=\sqrt{H/M}$ is common across 3 families × mixture). That is,

$$\text{initial structure}\ \to\ \text{memorized in}\ J_0,\ \textstyle\prod_t J_t\ \to\ \text{inflation}\ \to\ \text{saturation}\ \to\ \text{universalization of amplitude information}.$$

Connection with Paper 5 (amplitude forgets, orientation remains): the degenerate vacuum family of Paper 5 possesses both a **common amplitude scale** $|z_e|=\sqrt{H/M}$ and a **seed-dependent phase configuration and orientation of the face (geometric moduli)**. The radial coordinate $r=\sqrt{H_\perp/H}$ of Paper 5 ($0.19$–$0.38$ per vacuum) is the order parameter of the tilt of the face, and it is a **distinct quantity** from the per-edge amplitude $a_\infty=\sqrt{H/M}$ of this paper (they are not identified). Combining the two, the terminal appearance is

$$\boxed{\text{the amplitude scale }|z_e|=\sqrt{H/M}\text{ is universal (becomes a state quantity) / the phase configuration and the orientation of the face are degenerate and selected by the seed (retain the SSB history)}}$$

—— **at the terminal state the amplitude forgets the composition, and the orientation retains the history of SSB**. This is the unification of Paper 5 (degenerate vacuum), Paper 6 (transient inflation), and Paper 8 (terminal state quantity).

## 7. Claim Classification & Open Items

- Classification: numerical verification (additivity law $a_\infty=\sqrt{H/M}$, 3 families × mixture) + mathematical fact (impossibility of factorization). Verdict: retained.
- Note: the terminal equal amplitude has the phase normalization of this engine (stage 2, which kills the amplitude) built in, so it is not presented as an independent corroboration of the norm form (suspicion of circular derivation; independent verification with an amplitude-preserving control dynamics is open).
- Verified (re-run C): the 3-way decomposition of the edge space $1\oplus(N-1)\oplus N(N-3)/2$ (line graph / Johnson scheme) exists in reality and the terminal $\sigma_{\max}\to N-1$ is confirmed for all $N$. However, the terminal state does not ride on the $(N-1)$-dimensional vertex sector but spreads into the residual sector (frac_rest $\to0.93$) —— **"$\sigma=N-1$ (the eigenvalue) ↔ vertex sector (the dimension)" is a coincidence of dimension labels and does not correspond** (rejected).
- Open: quantification of the immediate departure of the mixed system (union onset).

## Related work (structural agreement, not the source of derivation)

- Barycentric decomposition of mean-field / all-pairs-coupled systems ($K_k$ spectrum $\{0,k,\ldots,k\}$): the reduction N-body = conserved barycenter + independent two-body problem.
- Extensive / intensive variables (thermodynamics): $a_\infty=\sqrt{H/M}$ is a state quantity of the extensive variables $(H,M)$.

## References

1. Noriaki Kihara, "The Mechanism of Inflationary Rapid Expansion in Self-Consistent Relational-Wave Closed Systems," Concept DOI 10.5281/zenodo.22112008.
2. Noriaki Kihara, "Onset-Mode Discrimination (Amplification and Unstable Relative Equilibria)," Concept DOI 10.5281/zenodo.21798854.

## Reproduction

Wrappers, parents, and figures are packaged together in each set: `ロールバック対照テスト_正本実装N3_N6_20260908/` (Experiments 20 & 21: `run_single_simplex_N9_v1.py`, `make_normalized_parents_N7N6_v1.py`, `run_mixed_simplex_N7N6_v1.py`, `plot_Hperp_mixed_parts_N7N6_v1.py`, 8 figures), `加法法則_90度系_N9_vs_N7N6_20260909/` (Experiment 23), `加法法則_90度解析床_N9_vs_N7N6_20260909/` (Experiment 25), `二体三体問題_関係波の検討_20260908/` (Experiment 22: `enumerate_triangular_decompositions_v1.py`, `三角数分解一覧_N3_N40.txt`, `分析_二体三体問題_関係波と凝縮体頂点階層_20260908.md`). The test of the correspondence between the vertex sector and $\sigma=N-1$ (re-run C) is in `位相ロック全N確認_相対平衡_20260912/頂点セクター_sigmaN1_20260912/` (`analyze_vertex_sector_sigmaN1_20260912.py`, `results/vertex_sector_summary.csv`, README / REPORT / SHA / run_all; read-only). Each includes README, SHA256SUMS, run_all.sh, and execution logs. The wrapper is a faithful copy of the canonical `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`) with only minimal changes to the initial values and output names.
