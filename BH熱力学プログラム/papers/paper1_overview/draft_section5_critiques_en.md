# 論文1 §5 限界と批判の先回り 試案

---

## §5. Anticipated Limitations and Critiques

This section addresses the principal critiques that we expect the program to encounter. We state each critique in the form a careful reader might raise, and offer the program's response.

### §5.1 Why a four-plus-one-dimensional formulation?

**Critique.** The observed universe has three spatial dimensions and one temporal dimension. Black holes that have been theoretically studied in detail and whose thermodynamic properties are anchored in the literature are 3+1-dimensional. Why does this program take the natural comparison object to be a 4+1-dimensional Tangherlini black hole?

**Response.** The 4+1-dimensional formulation is a consequence of the geometric construction adopted in Paper 2: gnomonic projection from a four-sphere $S^4 \subset \mathbb{R}^5$ onto a four-dimensional subjective chart. The induced metric on the four-dimensional chart, equipped with the curvature radius $R$ as an additional independent geometric parameter, gives rise naturally to a five-dimensional structure $(x, y, z, t, R)$. The black hole metric that arises in this construction is, on dimensional grounds, the Tangherlini analogue rather than the Schwarzschild form.

We do not claim that 4+1 is the correct dimensionality of physical spacetime. We claim that, *given the geometric construction of Paper 2*, the natural comparison is with the 4+1 case, and that the agreement with Tangherlini thermodynamics serves as a non-trivial check on the construction itself. The relation to observed 3+1-dimensional physics — whether through a Randall–Sundrum-style brane localization, an ADD-style large extra dimension, or some other mechanism — is left for future work and is not within the scope of this program.

This positioning is shared with a substantial fraction of the higher-dimensional gravity literature: Tangherlini, Myers–Perry, Randall–Sundrum, and ADD all study generalizations whose physical realization in 3+1 dimensions is a separate question from their internal consistency.

### §5.2 Why not assume the holographic principle?

**Critique.** The holographic principle and AdS/CFT have produced a great deal of progress in black hole thermodynamics, including detailed entropy counts and explicit dynamical descriptions. Why does this program not make use of them?

**Response.** The program is constructed as an independent path of inquiry. Its starting assumptions — central projection from a higher-dimensional spacetime and discreteness at the Planck scale — are logically independent of the holographic principle. Whether the resulting picture is compatible with, supports, contradicts, or is orthogonal to holography is a question that can be examined only after the independent path has been worked out. To assume holography at the outset would obscure that question.

We emphasize that nothing in the program relies on the failure of the holographic principle. If holography is correct, the agreement with Tangherlini thermodynamics derived in Papers 4 and 5 must, by some mechanism, be consistent with the holographic counting; the precise correspondence is an open problem that the program does not attempt to solve.

### §5.3 Why not derive the results fully and quantitatively?

**Critique.** Paper 1 only states hypotheses; the quantitative derivations are deferred to companion papers. The companion papers themselves leave the evaporation-rate exponent ($n = 3/2$ in the simplest treatment, $n = -1$ in the standard Hawking computation) unreconciled. Is the program offering anything beyond a sketch?

**Response.** Yes. The program offers (i) a precise mathematical framework — the central projection construction and the integer-lattice packing problem — that admits rigorous analysis in its own right; (ii) the explicit asymptotic constant $c = 8/(3\pi)$ for the volume deficit, derived analytically and confirmed numerically; (iii) the matching of the $R^3$ scaling and the leading numerical coefficient with the 4+1-dimensional Bekenstein–Hawking entropy; and (iv) the recovery of the standard Tangherlini scaling laws $R \propto M^{1/2}$, $T_H \propto 1/\sqrt{M}$, and lifetime $\propto M^2$.

The unreconciled evaporation-rate exponent is a genuine open problem. The program is explicit that the exponent obtained from the simplest dynamical assumption ($n = 3/2$) does not match the standard semi-classical result ($n = -1$). We mark this as an open problem rather than concealing it, and we conjecture that a more careful treatment of the back-reaction or of the discrete drift mechanism is required to resolve it. This honesty about discrepancies is part of the methodology, not a defect.

### §5.4 Why is the program testable?

**Critique.** The hypotheses are stated as conditional propositions and the program does not predict new observations. In what sense is it scientific?

**Response.** The program is testable in three distinct senses:

1. *Mathematical consistency.* Papers 2 and 3 are entirely mathematical: any error in the central projection metric, the asymptotic constant $c$, or the Jacobi-style closed-form representations of $N(k)$ would falsify the foundational claim. The numerical computation is reproducible from the definition; the analytic derivation is checkable line by line.

2. *Quantitative correspondence.* Papers 4 and 5 produce specific numerical relations — the constant $c = 8/(3\pi)$, the ratio $\Delta/S_{BH} \to 32/(3\pi) \approx 3.40$, the Tangherlini scalings $R \propto M^{1/2}$ — that are either matched or not. They are matched, with the noted exception of the evaporation-rate exponent. Failure of these correspondences would falsify the program's main claim of contact with established thermodynamics.

3. *Open-problem identification.* Paper 6 formalizes the residual discrepancies. The program identifies the evaporation-rate exponent as the principal locus of disagreement with the standard semi-classical computation. This is itself a useful research output: it specifies what would have to be modified, in either the program or the standard computation, to achieve full consistency.

The program does not predict new observable signatures. We do not claim it as an empirical contribution to physics; we claim it as a mathematical and conceptual exploration of what follows from a particular set of geometric assumptions. Whether the exploration turns out to have empirical content is a question for future work that lies outside the scope of the present program.

### §5.5 What is *not* claimed

To prevent misreading, we list explicitly what the program does **not** claim:

- It does not claim that physical spacetime is four-plus-one-dimensional.
- It does not claim that the central projection construction is the unique or correct geometric description of gravity.
- It does not claim to have derived the Bekenstein–Hawking entropy from first principles.
- It does not claim that the evaporation-rate exponent is $3/2$; it identifies the discrepancy with the standard $n = -1$ as an open problem.
- It does not claim to compete with, displace, or refute the holographic, loop, causal-set, or string-theoretic programs.

What it claims is the conjunction of two conditional propositions (Hypotheses 1 and 2), the mathematical structures that follow from them, and the partial quantitative correspondence with established Tangherlini results.

---

**Style notes**:
- Q&A 形式で予想批判を先回り
- 各応答は控えめかつ誠実
- 最後に「主張しないこと」リストで明確化
- §5.5 は §0 風の主張範囲明示として機能
