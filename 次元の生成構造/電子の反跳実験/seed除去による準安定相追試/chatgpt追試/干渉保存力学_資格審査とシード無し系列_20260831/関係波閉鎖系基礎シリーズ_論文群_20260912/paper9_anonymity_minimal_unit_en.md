# Anonymity and the Minimal Unit of Information — The minimal zero-closure unit is 2 waves; how many waves is the minimal unit of information that survives after quotienting by $S_N$? (Hypothesis / Exploratory Paper)

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 9 of 10 / Exploratory Paper)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22729128
**Concept DOI:** 10.5281/zenodo.22729127
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22729128

---

## Abstract

This paper is a **hypothesis / exploratory paper** (no new runs; it organizes algebraic facts + existing empirical data + literature cross-checks, with claim classification made explicit). We examine, from the standpoint of anonymity, "what can carry information" in a relational-wave closed system. (1) The edge index $m$ is a gauge of the vertex permutation $S_N$, and "the type of wave number $m$" is not a physical question (mathematical fact). (2) The final state is exactly equimodular (relative spread $N=12{:}1.17\times10^{-13}$, $N=40{:}3.54\times10^{-14}$), and under 12-digit rounding there is only one kind of $|z|$ — indistinguishable even by value, and the indistinguishability of identical particles appears as a consequence of the dynamics (numerical verification). (3) The nontrivial solution of closure $\sum z^2=0$ is minimally 2 waves ($z_2=\pm i z_1$); a single wave admits only $z=0$ and cannot carry information (mathematical fact = **minimal closure unit is 2 waves**). However, $\pm$ flips under pair exchange ($S_2$), so on its own it is not an $S_N$-invariant bit. The primitive block size is floor-dependent, being $p=\operatorname{spf}(L)$ of Paper 4 ($L=N/\gcd(N,4)$; only $8\mid N$ reaches the algebraic lower bound 2). Hence **the carrier of information is not the individual wave but the minimal zero-closure block** (Lam–Leung); but "2 waves = minimal closure unit" does not mean "2 waves = anonymous information unit." Analyzing $\mathcal C_k\simeq Q_{k-2}/S_k$ (projective quadric hypersurface $\dim_{\mathbb C}=k-2$), **closure alone does not generate discrete anonymous species** ($k=2$ collapses to 1 orbit under $S_2$; $k\ge3$ is a continuous moduli). The discrete invariant is not $\pm$ but the block degree $p=\operatorname{spf}(L)$ (Paper 4; the closure axiom + the discrete symmetry of the floor quantizes it), and the self-consistent dynamics that quantizes the continuous closure moduli into discrete types is the next question (§5). (4) A structure that records anything beyond $a,b$ on a wave (the band/hair register of prior work) is a **scaffold** not derivable from the closure axiom or from anonymity, and the types read there (boson/fermion) are not cited as a physical claim (downgrade of the verdict). Classification: we distinguish mathematical fact / numerical verification / hypothesis.

---

## 1. Problem

The waves of this system are anonymous (they carry no external index). We examine, using algebra and existing empirical data, what information can reside in an anonymous object and whether distinctions such as particle species can be made at all. This paper does not assert conclusions; it classifies and fixes what is an established fact and what is a hypothesis.

## 2. The edge index is an $S_N$ gauge (mathematical fact)

In the generator $K_{ef}=A_{ef}\sin(\varphi_f-\varphi_e)$, $A$ is the edge-graph adjacency of $K_N$, and the dynamics is equivariant under the vertex permutation $S_N$. The numbering of edges is a product of the choice of vertex labels; permuting them moves to a different wave. Therefore "the type of wave number $m$" is not a physical question, and **what is readable are the $S_N$-invariants** (totals, multisets, the count "how many waves of type $X$"). What is not readable is "which wave" (identification / tracking).

This is the **same principle** as the relational readout of Paper 7 (the absolute normal cannot be recovered from the state alone; what is readable is the relative phase between faces) — physical readout resides not in absolute labels (which wave / which orientation) but in relational quantities invariant under transformation. Adding Paper 8: (i) **kinematic anonymity** = the wave index is an $S_N$ gauge from the outset, (ii) **dynamical anonymization** = at termination even the $|z_e|$ become equal so identification by value also vanishes, (iii) **what remains** = the relational invariants of phase relations, relative arrangement of faces, and the type/multiplicity of closure blocks. This is the relationalism of the whole series.

## 3. The final state is indistinguishable even by value (numerical verification)

The final state is exactly equimodular: relative spread $1.17\times10^{-13}$ for $N=12$ and $3.54\times10^{-14}$ for $N=40$. Under 12-digit rounding there is only one kind of $|z|$ ($66\to1$, $780\to1$). At step0 they are all distinct ($66/66$, $780/780$), so the anonymization is **driven by the dynamics**. That is, the indistinguishability of identical particles appears not as an initial condition but as a consequence of the dynamics (a different expression of the same fact as the equal amplitude of Paper 5 and $a_\infty=\sqrt{H/M}$ of Paper 8, and consistent with the anonymity H-theorem $S\to1$ of Experiment 7).

## 4. The minimal unit of information is 2 waves (mathematical fact + hypothesis)

Nontrivial solutions of closure $\sum z^2=0$:

| Wave count | Solution | Degrees of freedom |
|---|---|---|
| 1 | $z^2=0\Rightarrow z=0$ only | None = cannot carry information |
| 2 | $z_2=\pm i z_1$ (minimal nontrivial closure unit) | all of $z_1$ (scale, phase). The sign $\pm$ flips under pair exchange (see below) |
| 3 | cube-root type $z_k=z\,\omega^k$, etc. | wider solution variety |

**The minimal closure unit is 2 waves (mathematical fact)**: a nonzero closure block requires at least 2 waves; a single wave admits only $z=0$ and cannot carry information.

**However, $\pm$ on its own is not an $S_N$-invariant 2-valued quantity (mathematical fact)**: pair-exchanging ($S_2\subset S_N$) the pair $z_2=+iz_1$ gives $z'_1=z_2,\ z'_2=z_1$ so that $z'_2=-iz'_1$ — **the very $S_2$ action that renders the two waves of the pair anonymous flips the sign**. Hence $+i\leftrightarrow-i$ lie on the same $S_N$ orbit, and to physicalize a two-valued quantity such as right-/left-handed one needs additional relational structure that can compare the signs (propagation direction, orientation, a third relation; the helicity discussion of Paper 4). That is, **the existence of a closure unit $\neq$ the existence of an anonymous information bit**.

**The primitive block size is floor-dependent (Paper 4)**: by the general rule $L=N/\gcd(N,4)$, $p=\operatorname{spf}(L)$ of Paper 4, the primitive block takes $p=2,3,5,7,\ldots$ ($N=10\Rightarrow L=5\Rightarrow p=5$; odd prime $N\Rightarrow p=N$; Lam–Leung). Therefore

$$k_{\min}=2\ (\text{absolute algebraic lower bound})\qquad\text{vs}\qquad k_{\rm primitive}=\operatorname{spf}(L)\ (\text{the minimal unit realized on that floor}),$$

the "minimal **possible** unit" and the "minimal unit actually allowed" differ, and only $8\mid N$ reaches the lower bound 2.

**Consequence (hypothesis) and the next question**: the carrier of information is not the individual wave but the minimal closure block (a single wave cannot close and cannot carry meaning). However, by the above, "2 waves = minimal closure unit" does not mean "2 waves = anonymous information unit." What should really be asked is the minimal $k$ such that discretely distinct types still remain after quotienting the closure solutions by symmetry:

$$\mathcal C_k=\Big\{(z_1,\ldots,z_k):\ \textstyle\sum_{j=1}^k z_j^2=0\Big\}\Big/(\text{global phase, scale, }S_k)$$

what is the minimal $k$ such that discretely distinct types remain (= **minimal anonymous information unit**)? $k_{\min}=2$ is the minimal closure unit. The minimal anonymous information unit is analyzed in §5 ($\mathcal C_k\simeq Q_{k-2}/S_k$). Candidate particle species (block size, multiplicity, way of assembly) can be defined as $S_N$-invariant multisets and do not break anonymity, but $\pm$ on its own is not an invariant bit.

## 5. The closure moduli are continuous — discrete anonymous species require additional self-consistent structure ($\mathcal C_k\simeq Q_{k-2}/S_k$)

The $\mathcal C_k$ of §4 can be analyzed without further computation. Quotienting phase and scale together by $z\sim cz$ ($c\in\mathbb C^\times$), $\sum_{j=1}^k z_j^2=0$ becomes a **projective quadric hypersurface** on $\mathbb{CP}^{k-1}$,

$$Q_{k-2}=\Big\{[z_1:\cdots:z_k]\in\mathbb{CP}^{k-1}:\ \textstyle\sum_j z_j^2=0\Big\},\qquad \dim_{\mathbb C}Q_{k-2}=k-2,$$

and $\mathcal C_k\simeq Q_{k-2}/S_k$.

- **$k=2$**: the projective solutions of $z_1^2+z_2^2=0$ are the 2 points $[1:i],[1:-i]$ ($\dim_{\mathbb C}Q_0=0$), but $S_2$ exchanges the two, so $Q_0/S_2$ is **1 orbit** — **for 2 waves there is closure but no anonymous 2-valued type** (an alternative proof of the $\pm$ flip of §4).
- **$k\ge3$**: $\dim_{\mathbb C}Q_{k-2}=k-2\ge1$. Since $S_k$ is a finite group, taking the quotient does not remove the continuous dimension near a generic point, and $\dim_{\mathbb C}(Q_{k-2}/S_k)=k-2$ — the solution space is still a **continuous moduli**.

Therefore **the closure condition $\sum z^2=0$ alone does not, after quotienting by $S_k$, produce a finite number of discretely distinct types (particle species)** ($k=2$ collapses to 1 species under $S_2$; $k\ge3$ retains continuous degrees of freedom). Hence the question shifts from "how many waves can carry information" to

$$\boxed{\text{which self-consistent dynamics quantizes the continuous closure moduli into discrete types}}$$

**Connection with Paper 4 (the carrier of quantization)**: in Paper 4, the **discrete block size** $p=\operatorname{spf}(N/\gcd(N,4))$ emerged from the cyclic structure of the high-symmetry floor. Although the closure equation alone is a continuous moduli, **imposing the discrete symmetry of the high-symmetry floor quantizes the block size**:

$$\text{closure axiom}\ +\ \text{discrete symmetry of the floor}\ \Longrightarrow\ \text{quantization of block size }p=2,3,5,7,\ldots.$$

At the present stage the invariant that is reliably discrete is not $\pm$ but the **block degree $p$ itself**, and whether a discrete invariant corresponding to a further particle species arises within it is the next problem. Rather than adding, from outside the closure, the register of Paper A (128 declared constants), this paper points in the direction that discrete types must be generated **from within** by closure + floor + dynamical selection rules (orientation, winding, topology, etc.).

## 6. Recording anything beyond $a,b$ on a wave is an anonymity violation (downgrade of the verdict)

- The waves of this system are a single complex number ($z\in\mathbb C^M$, only $a+ib$) = the minimal description.
- Prior work (Paper A) assigns each wave $16\times8=128$ complex numbers (band $k$ / hair $\eta$ registers). $N_n=16,N_\eta=8$ are **declared constants**, derived neither from $N$ nor from the closure axiom $\sum z^2=0$. A wave becomes distinguishable by internal structure = anonymity is broken from outside.
- Therefore the types read within that register (the boson/fermion distinction, "the clock is ticked by fermions," "charge = closure defect") are **properties of the scaffold, not physics derived from the closure axiom**. This series does not cite these as physical claims (downgrade of the verdict). To treat the types as physics, one must derive the band-register-equivalent degrees of freedom from the closure axiom and anonymity (currently not achieved).

## 7. Claim classification / open

- **Mathematical fact**: the $S_N$ gauge (§2), the minimal nontrivial closure unit of $\sum z^2=0$ = 2 waves $z_2=\pm i z_1$ (§4; but $\pm$ flips under pair exchange $S_2$ and on its own is not an $S_N$-invariant bit), the primitive block size $=\operatorname{spf}(L)$ ($L=N/\gcd(N,4)$; only $8\mid N$ reaches the algebraic lower bound 2, §4 · Paper 4), the Lam–Leung decomposition.
- **Numerical verification**: the final state is exactly equimodular (§3).
- **Hypothesis**: the carrier of information = the minimal closure block (not the individual wave), particle species = $S_N$-invariant multisets (block size, multiplicity, way of assembly) (§4).
- Verified (Reproduction B): the test of the block decomposition of the final state was carried out for all $N=3$–$40$ (existing 10000-step terminal states). The terminal state is exactly equimodular (deviation $\sim10^{-13}$) and conserves $\Sigma z^2\approx0$, but **exact $\pm i$ pairs are zero** (the make_parent floor as well; the minimal $\varepsilon$ only decreases asymptotically to $\sim2\times10^{-6}$ at large $N$). That is, minimal-block ($\pm i$ pair) decomposition is **specific to the cyclic-Fourier / integer-K high-symmetry floor ($N=8k$, Paper 4)** and does not stand at the terminus of the data series. Hence "terminal block decomposition = catalog of species" does not hold as a property of the terminal state, and block structure (the type corresponding to species) resides on the high-symmetry-floor side (Paper 4) — consistent with the thrust of this paper (the carrier of information is the block; the type compatible with anonymity resides in the closure structure of the floor).
- **Analytic answer (§5 · mathematical fact)**: closure $\sum z^2=0$ alone does not generate discrete anonymous species — with $\mathcal C_k\simeq Q_{k-2}/S_k$ (projective quadric hypersurface $\dim_{\mathbb C}Q_{k-2}=k-2$), $k=2$ collapses to 1 orbit under $S_2$, and $k\ge3$ is a continuous moduli ($\dim_{\mathbb C}=k-2\ge1$, not removed by the finite group $S_k$). What is discrete is not $\pm$ but the **block degree $p=\operatorname{spf}(N/\gcd(N,4))$** (Paper 4), where the closure axiom + the discrete symmetry of the floor quantizes the block size.
- Open / next research question: what is the self-consistent dynamics that quantizes the continuous closure moduli into discrete types (dynamical selection rules, orientation, winding, topology, etc. should generate them from within)? Does a discrete invariant corresponding to a further particle species arise within the block degree $p$? Reproduction B showed that exact $\pm i$ pairs ($p=2$ blocks) are specific to the high-symmetry floor $8\mid N$ and vanish at the generic terminus (the type corresponding to species resides not in the individual wave but on the block side, compatible with anonymity).

## Related work (structural agreement, not the source of derivation)

- Vanishing sums (Lam–Leung): prime-size decomposition of zero-closure blocks.
- Indistinguishability of identical particles (quantum statistics): in this system it appears as a consequence of the dynamics.

## References

1. Noriaki Kihara, "The mechanism of inflationary rapid expansion in self-consistent relational-wave closed systems," Concept DOI 10.5281/zenodo.22112008.
2. Noriaki Kihara, "Conditions under which a seed gives rise to particles and the lower bound of resolution (Paper A)," Concept DOI 10.5281/zenodo.21874481. (Source of the band/hair register referenced in §5.)
3. T. Y. Lam and K. H. Leung, "On vanishing sums of roots of unity," J. Algebra 224 (2000) 91.

## Reproduction

This paper is a hypothesis / exploratory paper with no new runs. The numerical verification (final-state equimodularity) relies on the existing data `N3_N40_long10000_20260905/` (10000-step sweep, `check_final_states_10000_v1.py`, etc.). The $N$-dependence test of the terminal block decomposition (Reproduction B) is in `位相ロック全N確認_相対平衡_20260912/終端状態ブロック分解_N依存_20260912/` (`analyze_terminal_block_decomposition_20260912.py`, `results/terminal_block_summary.csv`, README / REPORT / SHA / run_all; read-only). The closure factorization is shared with the package of Paper 4. The algebra (the minimal solution of $\sum z^2=0$, Lam–Leung) is self-contained in the main text. The primary record of the study is in the discovery ledger `発見台帳_干渉保存力学_20260905起点.md` (Studies 25, 26, 27).
