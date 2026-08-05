# Generation of Fermionic Structure Is Induced, Autocatalytic, and Pair-Correlated: Assumptions and Consequences of a Universal Inelastic Map Acting on Waveforms Alone

**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-08-05　**Version:** v1
**Version DOI:** [10.5281/zenodo.21808092](https://doi.org/10.5281/zenodo.21808092)
**Concept DOI:** [10.5281/zenodo.21808091](https://doi.org/10.5281/zenodo.21808091)

---

## Abstract

For a two-wave relational wave system we design, as an explicit assumption, a single inelastic interaction operation that **takes only the waveforms as input and contains no discriminating branches (no IF statements)**. Within the assumed class (two-wave, pointwise, cubic, common-phase invariant, no self-scattering, strength reusing the existing readout $R=\sin^2\theta$), conservation of the zero-square-sum closure (the axiom) uniquely fixes the relative coefficient, and the vertex is restricted to the single form $\delta a = i R\,(|b|^2 a - b^2 \bar a)$. Measuring the consequences of this assumed map numerically, we find: (1) purely even-harmonic pumps generate no odd content identically (a parity theorem — generation is necessarily an induced process); (2) over a four-decade sweep of the seed fraction $f_0$ the initial growth rate obeys the autocatalytic ignition law $\mathrm{rate}=C f_0^2$ with $C=11.45$ constant to 0.4%; (3) after ignition there is no runaway — the system relaxes oscillatorily to a statistical equilibrium $f^*=0.4690$, close to the fermionic-mask phase-space fraction 0.494; (4) although no antiparticle degrees of freedom are input anywhere, products appear only in the partner bins dictated by the sum rule (the 485 unpredicted empty bins remain at machine zero, ratio $10^{-27}$), as correlated pairs that settle the hair (winding-number) ledger (coherence 0.83), rising simultaneously from exact zero. Furthermore, the uniquified vertex rewrites exactly as $\delta a = -2R\,\mathrm{Im}(\bar b a)\,b$, $\delta b = +2R\,\mathrm{Im}(\bar b a)\,a$ — the driving scalar is the single imaginary part of the relative phase of the two waves, and the vertex is a pointwise two-channel rotation that is **exactly solvable in closed form** (closure and power are conserved identically pointwise; measured drift $4\times10^{-14}$ per 3000 collisions). This paper is not a paper that proves the interaction form correct. We report with the three layers — assumption, conditional derivation, and measurement — strictly separated, and we measure the value of the assumption not by correctness but by consistency and productivity (the same laws are reproduced when embedded in the N-body model). All of the above shows that from a single waveform map with no particle species, no antiparticles, and no selection rules as input, a discrete parity class, autocatalytic generation, statistical saturation, and ledger-closing conjugate pairs emerge as output-side classifications.

---

## 1. The Character of This Paper — Separating Assumption, Derivation, and Measurement

This paper neither derives the generative interaction from the axioms nor proves that it is correct. What this paper does is three things, stated separately from the outset:

- **(a) Assumption**: we posit the form of the inelastic interaction as an explicit assumption (Section 4).
- **(b) Conditional derivation**: within the posited class, we derive what the axioms force (Sections 5, 6).
- **(c) Measurement**: we measure the consequences of the assumed interaction in numerical experiments (Sections 8–10).

What the reader can verify are the derivations (b) and the measurements (c); (a) is not an object of verification but the starting point of this paper. The moment this distinction is blurred, the paper loses its meaning.

## 2. Claims

**First — the most important design principle of this paper — we define the interaction as a single universal operation whose only inputs are the two waveforms.** Following the anonymity principle, the interior of the map contains no special-case branches (IF statements) whatsoever: no branch that discriminates particle species, no case analysis of "if boson," "if fermion," "if seed," and no threshold switches. The map receives only the waveforms $(a,b)$ and returns only waveforms. Even the collision strength is not an externally supplied parameter but is read out from the waveforms themselves (the reflection rate $R=\sin^2\theta$, with $\theta$ the existing waveform readout). Consequently, every "boson," "fermion," "seed," and "partner" appearing in this paper is a **readout-side classification** (a spectral mask), not a **dynamics-side label**. The map does not know what it is generating. This design is what makes all subsequent results nontrivial: the parity selection rule, the autocatalytic ignition, and the pair structure all appear as consequences of a single indiscriminate expression applied identically to arbitrary inputs.

**Second, we enumerate all assumptions.** On top of the existing universal interaction (the elastic part: a two-channel rotation by the waveform readout $\theta$, exactly conserving the combined AB power), we posit the inelastic extension in the following form: (i) two-wave; (ii) pointwise (products at each point of the register lattice); (iii) cubic (lowest nonlinear order); (iv) common-phase invariant; (v) no self-scattering; (vi) no new coupling constant — the strength reuses $R=\sin^2\theta$. **All six items are design choices, not derivations** — although every item is chosen to be consistent with the first design principle (waveforms as sole input, no internal branches), and (vi) is its direct consequence. Each item has a motivation (anonymity, operational algebra, lowest order), but motivation is not derivation.

**Third, within the posited class, Axiom 1 (conservation of the zero-square-sum closure) uniquifies the form.** For the most general form under the assumptions, $\delta a = i(g_1|b|^2 a + g_2 b^2\bar a)$, closure conservation forces $g_2=-g_1$. Hence within the assumed class the vertex is restricted to the single form $\delta a = igR(|b|^2 a - b^2\bar a)$. **This is conditional uniqueness, not uniqueness of the class itself.** It has moreover turned out that the uniquified vertex is exactly solvable in closed form as a pointwise two-channel rotation (Section 5.2) — closure conservation becomes a pointwise identity and no numerical integration is needed. Note also that what closure conservation demands of the strengths is pair symmetry only; the arithmetic-mean rule for mixed collisions is an **additional experimental fact** not derived from closure (the product form also conserves closure — we do not blur this distinction).

**Fourth, the parity theorem.** Under the assumed vertex, no odd-harmonic content is ever generated from purely even-harmonic pumps. Hence **in this model generation is necessarily an induced process**; spontaneous seedless generation does not occur. We prove this as a theorem and verify it by machine. A single expression with no IF statements nonetheless possesses a strict selection rule — the selection rule arises not from branching but from the exponent arithmetic of the expression.

**Fifth, measurement of the ignition law.** Over a four-decade sweep of the seed fraction, $\mathrm{rate} = C f^2$ with $C = 11.45$, constant to 0.4%. The process is autocatalytic, and the ignition time scales as $1/f_0^2$. The coefficient has been calibrated: the initially measured $C=10.4$ contained an integration bias of about 10% from the midpoint method. A RK4 re-run improved the zero-closure drift from $2.9\times10^{-3}$ to $2.8\times10^{-9}$, and the law and its constancy were invariant across both integrators (the history is recorded in Section 11).

**Sixth, measurement of the post-ignition fate.** Total conversion (runaway) does not occur. The vertex is reversible, and the system oscillates back and forth to a statistical equilibrium $f^* = 0.4690$. This value nearly coincides with the phase-space fraction of the fermionic mask, 0.494, and reads as equipartition thermalization.

**Seventh, measurement of the pair structure (three census criteria).** Although antiparticle-like degrees of freedom are input neither in the initial conditions nor in the map: (1) generation appears only in the partner bins dictated by the sum rule, and the 485 unpredicted empty bins remain at machine zero ($10^{-27}$ of the partner-band mean); (2) the partner shows coherence 0.83 with the seed at the position $q^*=+4$ predicted by the hair sum rule (the ledger is machine-derived from the output of a single vertex application — nothing hand-set); (3) the partner band rises from exact zero simultaneously with the sidebands. Conclusion: **inside this model, antiparticles are outputs, not inputs.**

**Eighth, we record the full process of refutation and correction (Section 11).**

**Finally, what this paper does not claim, and the status of the assumption (Sections 13, 14).** This paper does not claim that nature's generative interaction has this form. What raises the value of the assumption is not a proof of correctness but a demonstration of consistency and productivity. The assumption is falsifiable.

## 3. Notation and Foundation — the Elastic Universal Interaction

The foundation is the two-wave relational wave system of Paper 9 [1]. The state is a pair of complex waveforms $(a,b)$ (a $\chi$ lattice of 512 bins × an $\eta$ lattice of 16 points). The existing **elastic** universal interaction applies, at each collision,

$$\theta = \operatorname{atan2}(\sqrt{P_f}, \sqrt{P_b}), \qquad R = \sin^2\theta,$$

$$\begin{pmatrix} a' \\ b' \end{pmatrix} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}$$

where $P_f, P_b$ are the fermionic/bosonic relational quantities of the combined AB spectrum (the fermionic mask = the bin set with even $|k|\ge 4$). Hereafter we call $f = P_f/\lVert(a,b)\rVert^2$, the ratio of fermionic relational power to total power, the fermionic fraction. This rotation exactly conserves the per-bin power of the AB combination — it is **elastic** and has no channel that changes harmonic content. Hence generating fermionic structure from a bosonic sea requires an added inelastic channel. This is the question from which this paper starts.

The readout $\theta$ (and hence $R$) is determined by the waveforms alone. The construction and verification of this readout belong to Paper 9 [1]; in this paper the original code is referenced read-only as an unmodified snapshot (SHA-256 recorded; Section 15).

## 4. Assumptions — Design of the Inelastic Extension

### 4.1 Design principle (Claim One)

The inelastic extension must be a **universal operation** with the same standing as the elastic part: inputs are the waveforms $(a,b)$ only, the output is waveforms only, and there is no internal branch discriminating particle species or situations. The only conditional in the implementation is the numerical skip "do nothing if the strength is zero" (`if r > 0`), which is not a decision rule — the strength $r=R$ is itself a readout of the waveforms. This was confirmed by source audit (test T5).

### 4.2 Assumption list

| # | Assumption | Motivation (not derivation) |
|---|------|--------------------|
| (i) | Two-wave (no vertex of three or more waves) | A collision is a two-wave event |
| (ii) | Pointwise (products at each point of the $\chi\times\eta$ lattice) | Four-wave-mixing sum rules of type $k_1+k_2-k_3=k_4$ hold automatically as lattice arithmetic; translation invariant |
| (iii) | Cubic (lowest nonlinear order) | Vanishes in the dilute limit, reducing to the elastic operation (correspondence principle) |
| (iv) | Common-phase invariant | The convention that the overall phase is meaningless |
| (v) | No self-scattering (no $|a|^2a$-type term) | Operational algebra: "a single wave has no self-scattering" |
| (vi) | No new coupling constant; strength $=R=\sin^2\theta$ | Anonymity: introduce no new free parameter. The complement of $R$ (the transmission rate $1-R$) has been identified with the elementary charge [4] — the strength is a readout quantity that already has physical standing |

These six items are design choices. The motivations make each item **plausible**, but do not **force** it.

### 4.3 Implemented form

One collision = elastic part (unchanged) ∘ inelastic part:

$$a' = a + i\,R\,(|b|^2 a - b^2 \bar a), \qquad b' = b + i\,R\,(|a|^2 b - a^2 \bar b)$$

(simultaneous update, pointwise; the relative sign is fixed by the theorem of the next section. The numerical implementation uses the closed-form exact solution of Section 5.2 — the RK4 substeps are the historical implementation of the calibration process, recorded in Section 11.)

## 5. Theorem 1: Closure Conservation Uniquifies the Coefficient

The axiom is conservation of the zero-square-sum closure $\mathcal{C} = \sum a^2 + \sum b^2 = 0$. For the most general vertex under assumptions (i)–(v),

$$\delta a = i(g_1 |b|^2 a + g_2 b^2 \bar a), \qquad \delta b = i(g_1 |a|^2 b + g_2 a^2 \bar b),$$

the first variation is

$$\delta \mathcal{C} = 2\sum (a\,\delta a + b\,\delta b) = 2i\,(g_1+g_2)\sum\left(|b|^2 a^2 + |a|^2 b^2\right).$$

The sum on the right does not vanish in general, so $\delta\mathcal{C}=0 \iff g_2 = -g_1$. Hence **within the assumed class** the vertex is restricted to the single form of Section 4.3 (the strength absorbs $g_1$ into $R$ by redefinition).

We make three points explicit. First, this is **conditional uniqueness**: with the class (i)–(v) posited, the axiom selects the coefficient; the class itself has not been derived. Second, what closure conservation demands of the strengths is **pair symmetry only**. As the mixing rule when the two waves have different strengths, both the arithmetic mean and the product conserve closure exactly (numerically verified; T5 of [9]). That the actual dynamics is consistent with the arithmetic mean and inconsistent with the product is an **additional experimental fact** not derived from closure. Third, as a historical record, the first version v1 was implemented in the general form ($g_2=+g_1$) and closure drift was observed — the current form is the correction to $g_2=-g_1$ demanded by closure conservation (Section 11).

### 5.1 Corollary 1: the vertex is an exchange flow driven by the relative phase

By an identity transformation (verified to machine precision), the uniquified vertex can be written

$$\delta a = -2R\,\mathrm{Im}(\bar b a)\,b, \qquad \delta b = +2R\,\mathrm{Im}(\bar b a)\,a.$$

The substance of what looked like a cubic nonlinearity is **an exchange flow that reads the single quantity $s=\mathrm{Im}(\bar b a)$ — the imaginary part of the relative phase of the two waves — and transports amplitude toward the other wave**. Under the common phase $(a,b)\mapsto(e^{i\phi}a, e^{i\phi}b)$, $s$ is invariant, so the common-phase invariance of assumption (iv) can be read off explicitly from this form — the scalar driving the interaction is reduced to one antisymmetric relative-phase quantity, and the common phase has no way to enter the expression from the start.

### 5.2 Corollary 2: the vertex is a pointwise two-channel rotation, exactly solvable in closed form

The flow of Corollary 1 is a flow of the two-component vector $(a,b)$ along the tangent direction $(-b, a)$ — that is, a **pointwise two-channel rotation** (angular velocity $2s$). Moreover, that $s$ itself is exactly invariant under this flow can be shown by direct computation: for the flow $\dot a = -2s\,b$, $\dot b = +2s\,a$ ($s$ real),

$$\frac{d}{d\tau}(\bar b a) = \dot{\bar b}\,a + \bar b\,\dot a = 2s\,|a|^2 - 2s\,|b|^2$$

and the right-hand side is real. Hence $\dfrac{ds}{d\tau} = \mathrm{Im}\,\dfrac{d}{d\tau}(\bar b a) = 0$. Since the angular velocity is conserved along the flow, the inelastic part has a **closed-form exact solution** at each lattice point:

$$a' = \cos\varphi\,a - \sin\varphi\,b, \qquad b' = \sin\varphi\,a + \cos\varphi\,b, \qquad \varphi(x) = 2R\,\mathrm{Im}(\bar b a)(x).$$

There are three consequences. First, the closure $a^2+b^2$ and the power $|a|^2+|b|^2$ are **conserved identically pointwise** (invariants of a real rotation) — the geometric meaning of the uniqueness theorem (Section 5 main text) becomes transparent: "the closure-conserving direction is the tangent direction of rotation." Second, numerical integration becomes unnecessary: re-running all experiments of this paper in closed form reproduced everything — 3000-collision closure drift $4.1\times10^{-14}$ and norm drift $1.4\times10^{-14}$ (machine precision), with $f^*=0.4690$, the maximum of $f$ and its position, and $C$ (within 0.5%) all recovered (criteria E1–E3 fixed in advance; all passed). Third, it becomes visible that the elastic and inelastic parts share the **same rotation-generating direction**: one collision is "a global rotation by a finite angle read from the spectrum" ∘ "a pointwise rotation by an angle determined by the local relative phase," and the two parts can be read in a unified way as rotations of the same two-channel relational space. Generation (the appearance of new harmonics) arises from the rotation angle $\varphi(x)$ differing from lattice point to lattice point — the spectrum changes while the pointwise combined power never changes at all. This unified reading (the possibility of deriving both from one generator family) is listed as an unverified problem in Section 13.1.

## 6. Theorem 2: the Parity Theorem — Generation Is Necessarily Induced

**Theorem**: under the assumed vertex, if $a$ and $b$ both contain only even harmonics, the odd-harmonic content remains identically zero.

**Sketch of proof**: the vertex is a pointwise cubic product, so in spectral space it is a convolution composing bins as $k_1 + k_2 - k_3$. Even + even − even = even. Hence the even parity class is closed and no generation term for odd content exists. The selection rule arises not from branching but from the arithmetic of exponent sums. ∎

**Machine verification** (T3): colliding purely even pumps (no seed) at strong coupling $g=0.5$ for 200 collisions, the maximum of the fermionic relational power stays below $10^{-25}$ (machine zero).

**Consequence**: in this model, spontaneous seedless generation does not occur. Generation is necessarily an **induced** (seeded) process — the coexistence of a tiny odd seed and strong even pumps is the condition. This is the strict, discrete-map version of the fact that spontaneous emission does not appear in semiclassical radiation theory (Section 12). We state explicitly that the spontaneous version (seeded by zero-point fluctuations) is outside the scope of this map.

**Positioning**: the even/odd class of harmonic number is the quotient $\mathbb{Z}/2\mathbb{Z}$, and the sum rule $k_1+k_2-k_3$ conserves it. Parity thus appears not as an externally attached quantum number but as **the minimal discrete invariant conserved by the exponent arithmetic of the interaction**. We emphasize: group theory was not an input to the design. The accurate order of events is that when the outputs of the single IF-free expression were classified after the fact, a conserved class of $\mathbb{Z}_2$ type was already there.

## 7. Correspondence Principle and the Ledger of Conserved Quantities

- **T1 control**: at zero strength, bitwise identical to the original pipeline (the inelastic part vanishes completely).
- **T3′ dilute reduction**: the trajectory difference with/without the inelastic part vanishes as $s^2$ in the scale $s$ (measured ratios 3.95, 4.00 per halving — correspondence principle: reduction to the elastic operation in the dilute limit).
- **Conserved quantities**: by the closed form of Section 5.2, the combined AB norm and the zero closure are conserved **identically pointwise**. Measured (3000 collisions, strong occupation $S=8$, closed-form implementation): zero-closure drift $4.1\times10^{-14}$, norm drift $1.4\times10^{-14}$ (machine precision). The drifts observed with the RK2/RK4 implementations of the calibration process ($2.9\times10^{-3}$ / $2.8\times10^{-9}$) are retroactively confirmed to have been pure integrator error (Section 11). The combined $\chi$ momentum is conserved by the sum rule via the translation invariance of the pointwise operation.

## 8. Experiment 1: the Autocatalytic Ignition Law rate = C·f²

**Procedure** (preregistered): inject an odd seed (amplitude ratios $10^{-3}, 10^{-2}, 10^{-1}$ — four decades in the initial seed fraction $f_0$) into two even pumps (scale $S=8$), measure $\mathrm{rate}_0$ from the initial slope of $\ln P_f$ over 200 collisions ($j\le 20$), and test the constancy of $C = \mathrm{rate}_0/f_0^2$.

**Results** (Figure 1):

| Seed amplitude ratio | $f_0$ | $C=\mathrm{rate}_0/f_0^2$ |
|---|---|---|
| $10^{-3}$ | $4.77\times10^{-7}$ | 11.453 |
| $10^{-2}$ | $4.77\times10^{-5}$ | 11.462 |
| $10^{-1}$ | $4.74\times10^{-3}$ | 11.410 |

Over four decades of $f_0$, $C = 11.45$ with spread 0.4%. The growth rate is proportional to the **square** of the seed fraction — an autocatalytic law — and the ignition time scales as $1/f_0^2$: a sea with a small seed fraction does not ignite for astronomical times.

![Figure 1 ignition law](fig1_ignition_law_v1.png)

**Figure 1**: (a) The growth trajectories $(\ln f - \ln f_0)/f_0^2$ collapse onto a single common curve over four decades of $f_0$ (the extension of the $f^2$ law to the whole trajectory). The dashed line is the initial slope $C=11.45$. The common gain curve grows superlinearly with collisions (this superlinearity is common to all seed fractions and is not explained by the growth of the seed itself — identifying the mechanism is recorded as outside the scope of this paper). The largest seed (green) departing upward from the common curve in the second half is the direct appearance of autocatalytic acceleration by the current $f^2$. (b) Constancy of $C$. Gray: before calibration (RK2 integration, systematic bias of about 10%; Section 11).

## 9. Experiment 2: the Post-Ignition Fate — No Runaway, but Statistical Equilibrium

**Procedure** (preregistered): $S=8$, seed amplitude ratio $0.1$, 3000 collisions. Record the zero-closure drift, the norm drift, and the full series of $f$. Criteria: (R1) closure drift below 1/100 of the pre-calibration measurement $2.9\times10^{-3}$; (R2) $f^*$ (mean of the last 500) within $0.469\pm0.02$ with the oscillatory-type verdict retained.

**Results** (Figure 2): total conversion (runaway to $f\to1$) does not occur. $f$ rises steeply at $j\approx450$, passes a maximum $0.607$ ($j=1841$), oscillates, and settles at the statistical equilibrium $f^* = 0.4690$. This value nearly coincides with the phase-space fraction of the fermionic mask, $0.494$ (the ratio of mask bins to all bins), and reads as equipartition thermalization under the reversible vertex (the attribution of the ~5% difference between the two values is in Section 13.1). The closure drift is $4.1\times10^{-14}$ in the closed-form implementation, norm drift $1.4\times10^{-14}$ — machine precision (the calibration RK4 run gave $2.8\times10^{-9}$ / $8.6\times10^{-11}$; both R1 and R2 hold).

![Figure 2 fate](fig2_fate_v1.png)

**Figure 2**: the post-ignition fate. Blue is the RK4 calibration run; gray (RK2, pre-calibration) overlaps the blue almost completely and is invisible — the six-orders-of-magnitude difference in integrator precision does not change the trajectory or the equilibrium value (the physics is integrator-independent). Red dashed line = statistical equilibrium $f^*=0.4690$; green dotted line = mask phase-space fraction $0.494$ (candidate equipartition value).

## 10. Experiment 3: the Pair-Structure Census — Antiparticles Are Outputs, Not Inputs

**Design**: so that the pre-thermalization imprint separates bin by bin, the pumps are single-winding, narrow-band (harmonics $30,32,34$ → raw bins $\{29,31,33,35\}$) × two channels; the seed is single-winding harmonic $21$ (raw bin $22$, amplitude ratio $0.2$); the measurement window is limited to 40 collisions. **Preregistered predictions** (fixed and printed before execution):

- Partner band (the pair-generation output) = (pump bins ⊕ pump bins) − seed bin = $\{36,38,\dots,48\}$ (7 bins)
- XPM sidebands (phase modulation, not pairs) = seed bin ± pump differences = $\{16,18,20,24,26,28\}$ (6 bins)
- **All empty bins other than the above plus the pump-sum band and the seed neighborhood (485 bins) are predicted not to generate** (a machine-zero prediction)
- Hair ($\eta$ winding) ledger: the partner's hair is **machine-derived from the increment of a single vertex application to the initial state** (hand-setting forbidden). The derived value agrees with the sum rule of the pair-generation term $b^2\bar a$, $m^* = 2m_B - m_s = 2(+2)-(+1) = +3$, and the predicted winding offset of the pair correlation is $q^* = m_s + m^* = +4$

**Criteria and results** (all re-measured with RK4; agreement with the pre-calibration measurements to 4 digits):

- **P1 sum-rule exclusivity**: the median growth of the 485 unpredicted bins, $8.2\times10^{-32}$, is $1.4\times10^{-27}$ of the partner-band mean growth $5.8\times10^{-5}$ — machine zero. **Generation appears only in the predicted bins.**
- **P2 the hair ledger**: the $q$-scan coherence (anomalous pair correlation of seed × partner) is $0.832$ at $q^*=+4$ and $7.5\times10^{-8}$ at the uniform average ($q=0$) — not the absence of correlation, but the hair selection rule extinguishing the $q=0$ component. **The partner is born correlated with the seed in the form that settles the hair ledger.**
- **P3 simultaneity of the pair**: the partner band rises from an initial $1.5\times10^{-28}$ (exact zero), and its 40-collision growth is $0.45$ of the sideband growth — simultaneous output of the same vertex.

![Figure 3 census](fig3_census_v1.png)

**Figure 3**: (a) per-bin power growth (40 collisions). Only the predicted bins grow (red = partner band, blue = XPM sidebands, green = pump-sum band), and the 485 unpredicted empty bins (gray, the flat part on the right) stay pinned at machine zero. The gray skirt on the left is the initial structure of pumps and seed (not empty bins). (b) The anomalous pair correlation for each hair winding offset $q$. Only the prediction $q^*=+4$ stands at $0.83$; the others are extinguished below $10^{-8}$ by the selection rule.

**Conclusion**: antiparticle-like degrees of freedom are input neither in the initial conditions nor in the map. Nevertheless, the products appeared at the positions dictated by the sum rule, as correlated pairs carrying the hair ledger, simultaneously from exact zero. **Inside this model, antiparticles are outputs, not inputs.** More precisely, the non-circular definition this model provides is: **a partner (antiparticle candidate) is the conjugate relational state that the generating vertex co-produces in order to close all conservation ledgers (the bin sum rule $k^*=2k_p-k_s$ and the hair sum rule $m^*=2m_B-m_s$)**. It is not that a wave viewed alone is a particle or an antiparticle; the distinction is determined as a role within the generation event — the appearance, at the level of generation phenomena, of this series' design principle that relations precede attributes.

## 11. Record of Refutations and Corrections

Following the series convention (record the full process of hypothesis → refutation → correction), we list every place where this study actually went wrong.

| # | Error | Discovery | Correction |
|---|------|------|------|
| 1 | v1 vertex implemented in the general form $g_2=+g_1$ | Measured zero-closure drift | Derived $g_2=-g_1$ from the first variation of closure (Section 5). The error turned into the discovery of the uniqueness theorem |
| 2 | census v1: measured the pair imprint after thermalization (after 800 collisions) | Three criteria FAIL | Moved the measurement window to 40 collisions before thermalization |
| 3 | census v2: broadband pumps dispersed the partner over 444 bins; correlation measured with the $\eta$-uniform average | Correlation exactly zero | Changed to narrow-band single-winding design + hair-resolved correlation. The "zero" was not decorrelation but the hair selection rule |
| 4 | Hair ledger preregistered by hand calculation (assuming $m_B=-1$) | P2 FAIL (measured $q^*=+4\neq$ prediction) | The measured carrier hair of pump B is $+2$. Replaced the ledger with machine derivation from a single vertex application — what was wrong was the hand-set ledger on the judging side, not the sum rule of the map |
| 5 | Integration bias of RK2 (midpoint) substeps | Closure drift $2.9\times10^{-3}$; RK4 control in the many-body implementation | RK4 re-run: drift $2.8\times10^{-9}$ (six orders better). $f^*$, fate, census invariant to 4 digits; only $C$ calibrated $10.4\to11.45$ (a ~10% systematic error). Subsequently the discovery of the closed-form solution (Section 5.2) made the integrator itself unnecessary (drift $10^{-14}$), finally confirming this family of drifts as pure integrator error |
| 6 | Temporarily claimed "closure conservation selects the arithmetic mean of mixed strengths" | Numerical verification: the product form also conserves closure exactly | Closure demands pair symmetry only. The selection of the mean is demoted to an additional experimental fact of the dynamics (Section 5) |

## 12. Prior Work — Correspondences Found After Independent Study

The vertex of this paper was constructed from an independent axiomatic study (closure conservation) and was not derived by relying on the following prior fields. Noticing that the measured results (partner sum rule, induced-only generation, statistical equilibrium) are isomorphic to known phenomena, we re-surveyed the literature and found the following correspondences. The correspondences are supporting evidence that the assumed class is not physically vacuous; they are not proofs of correctness.

1. **The algebraic form of the vertex**: the coupled nonlinear Schrödinger equations of birefringent optical fibers [5] carry both the XPM term $|A_y|^2A_x$ and the coherent-coupling term $A_y^2 A_x^*$ — the same algebraic class as our $g_1|b|^2a + g_2 b^2\bar a$. In optics, however, the coefficient ratio is fixed by the $\chi^{(3)}$ tensor symmetry of the medium, whereas here closure conservation (the axiom) forces $g_2=-g_1$. Same class, different selection principle.
2. **The partner sum rule**: in seeded four-wave mixing [6] the idler is automatically generated at $\omega_i = 2\omega_p - \omega_s$ with phase-conjugate correlation — strictly isomorphic to our $k^* = 2k_p - k_s$. In optics, phase matching is imposed by the dispersion relation, whereas our sum rule comes from the exponent arithmetic of the pointwise cubic map. We find no counterpart of the hair-inclusive ledger and the machine-zero exclusivity ($10^{-27}$).
3. **Induced processes only**: in semiclassical radiation theory only stimulated emission and absorption appear; spontaneous emission does not [7]. The seed of spontaneous four-wave mixing is quantum vacuum noise. Our parity theorem is the discrete-map version of this fact — not an approximate statement but a strict selection rule with machine verification.
4. **Fermion generation from a bosonic background**: in cosmological preheating, the oscillating inflaton parametrically excites fermions, saturating under Pauli blocking [8]. The phenomenology is parallel but the ontology differs: preheating posits a Dirac field as input. We input no fermionic degrees of freedom — "fermion" is a readout mask. The saturation mechanisms also differ (Pauli blocking vs. equipartition equilibrium of a reversible vertex).
5. **Oscillation and thermalization**: FPU-Tsingou recurrence [10] and Rayleigh-Jeans thermalization/condensation of classical light waves [11] form the lineage in which nonlinear classical waves reach statistical equilibrium via recurrences. Our $f^*\approx$ mask phase-space fraction is consistent with the equipartition reading of this lineage. However, thermalization under an exactly conserved quadratic invariant (the zero closure, drift $10^{-9}$…$10^{-14}$) has no counterpart. The Manley-Rowe relations [12] are a distant analogue of "conservation laws constraining conversion paths," but those follow from Hamiltonian structure, whereas here conservation acts as the axiom that selects the vertex.
6. **No counterpart found for rate = C·f²**: in seeded parametric amplification the gain rate is set by the pump, and the seed enters only as a prefactor. An autocatalytic law in which the growth **rate** itself is proportional to the square of the seed fraction could not be found in the surveyed literature (second-order autocatalysis in chemical kinetics is a formal similarity only). We present it as a novelty candidate and welcome the presentation of prior examples (refutation).

## 13. Outlook — Embedding into the N-Body Model (Verified Connections)

The two-body map of this paper has been embedded into the N-body relational wave system (complete graph, multi-slice register) by a family of reproducible experiments in the public repository. We record only the main points (details in the scripts and commit series [9]).

1. **"When does the interaction happen" is resolved**: the interaction acts at all times; the apparent timing of generation is produced by two factors (condensation gate × seed fraction squared). The measured gate ratio of generation rate under condensation is 2.74.
2. **The uniquification of this paper is the M=2 special case**: on general graphs, closure conservation ⟺ $g_2=-g_1$ and pair symmetry of the mediated strengths (the mediated-vertex uniqueness theorem).
3. **Transplantation of the ignition law**: the many-body system reproduces the power 2.001 with $C$ constant to 0.4%.
4. **The micro-mechanism of the fate**: the generation direction is a phase condition (a 50/50 lottery under random phases), and net matter arises by rectification through the $C f^2$ autocatalysis (80% of the ensemble grows net) — the decomposition of this paper's "oscillatory equilibrium."
5. **Position and motion**: locality lives on the register (waveform-position) axis, and products are pinned at their birth position. Motion is provided by the harmonic-closure dispersion $\omega_k = k\,\omega_1$ [3] (shape-invariant translation demonstrated in 1D/3D).
6. **Transplantation of the census**: the three criteria — sum-rule exclusivity, hair ledger, simultaneity — hold in the many-body system up to the hair-complete version.

### 13.1 Open Problems (working guide — not claimed in this paper)

The following are analysis problems opened by the results of this paper; none of them is claimed here.

1. **Perturbative derivation of the $f^2$ law**: derive analytically, from the vertex and the initial spectrum, why the initial growth rate is the square of the seed fraction. Includes identifying the mechanism of the common superlinear gain of Figure 1(a) (common to all seed fractions — not explained by the current $f^2$ alone). The closed growth law may involve, besides $f$, an alignment quantity of relative phases (e.g. the distribution of $\mathrm{Im}(\bar b a)$).
2. **State-count derivation of the equilibrium value**: can the ~5% difference between $f^*=0.4690$ and the naive bin ratio $0.494$ be explained by the **measure of the reachable phase space** under closure, hair, and momentum conservation? If so, $f^*$ is promoted to a statistical-mechanical formula.
3. **Extraction of the generator algebra**: numerically extract the generators of the elastic rotation, the inelastic pointwise rotation (Section 5.2), and relative-phase operations, and measure whether the commutators close. If the elastic and inelastic parts come from one generator family, assumption (vi) approaches a derivation. Group theory is not input; the operators extracted from the map are classified after the fact.
4. **Reduction of the assumptions**: can (ii) pointwise, (iii) cubic, (iv) common-phase invariant, (v) no self-scattering be re-derived from a single higher principle (local, anonymous, relational, lowest-order)?

## 14. What This Paper Does Not Claim; Falsifiability

This paper does not claim that nature's generative interaction has the form of Section 4. Assumptions (i)–(vi) can be replaced by other choices, and all results of this paper are consequences of the conditional "if this class is assumed." What raises the value of the assumption is not a proof of correctness but the demonstration of **consistency** (no contradiction with the axioms, the existing elastic operation, and the conserved quantities) and **productivity** (the same vertex reproducing the same laws in the N-body model; Section 13).

Falsification conditions: if any of the predictions of this class — sum-rule exclusivity (machine zero of unpredicted bins), the hair ledger (the machine-derived value of $q^*$), the coincidence of the statistical equilibrium $f^*$ with the mask phase-space fraction — is broken in a more constrained experimental system, the class as a whole is rejected. If a prior example of the autocatalytic law of Section 12-6 is presented, the novelty claim is withdrawn.

## 15. Reproducibility

We follow the series convention (self-contained reproduction within the series; code of past papers is imported read-only after snapshot verification).

- **Original**: the Paper 9 execution environment `ab_invariant_theta_toy_v1` (unmodified; snapshot SHA-256 90c7b272…, bit-identity confirmed).
- **Code of this paper** (repository folder `universal_inelastic_map_managed_v1`, under the series folder for dimensional generation structure):
  - `universal_inelastic_map_v1.py` … v1 implementation + T1 control / T2 conservation / T3 parity / T4 seeded gain / T5 IF audit
  - `universal_inelastic_map_v3.py` … final map (strength $=R$) + dilute reduction and ignition sweep
  - `run_ignition_fate_rk4_v2.py` … RK4 fate re-run + $C$ sweep (criteria R1–R3 fixed in advance; calibration record)
  - `run_ignition_fate_exact_v3.py` … verification of the closed-form exact solution (Section 5.2) and re-run of fate and sweep (criteria E1–E3 fixed in advance; all passed)
  - `run_pair_structure_census_v3.py` … final census (with preregistered predictions printed)
  - `make_paper_figures_v1.py` … generation of Figures 1–3 + census RK4 re-measurement
- **Main parameters**: $\chi$ 512 bins × $\eta$ 16, pump harmonics $2\!-\!62$ (even), $S=8$, RK4 substeps $h_{\max}=0.01$ (historical), closed-form implementation for final numbers.
- **Result JSONs**: `ignition_fate_rk4_result_v2.json`, `ignition_fate_exact_result_v3.json`, `pair_structure_census_result_v3.json`, `pair_structure_census_rk4_result_v1.json`, etc. (all numbers in the text agree with these files).

## References

**Within the series (dependence indicated)**

[1] N. Kihara, *The Generative Structure of Fermions* (2026), Zenodo. Concept DOI: 10.5281/zenodo.21766706 (original of the elastic universal interaction, the mask, and the execution environment; Section 3)

[2] N. Kihara, *Geometric Rapid Expansion Is Specific to Unstable Self-Consistent Closures* (2026), Zenodo. Concept DOI: 10.5281/zenodo.21798854 (motivation and the correspondence "amplification ⟺ unstable equilibrium"; Sections 1, 6)

[3] N. Kihara, *Future-Phase Position-Acceleration Map in an AB Two-Body Closed Phase System and the Inverse-Square Law via Harmonic Closure* (2026), Zenodo. DOI: 10.5281/zenodo.21441081 (harmonic-closure dispersion $\omega_k=k\,\omega_1$; Section 13-5 only)

[4] N. Kihara, *Two-Grammar Decomposition of Interaction* (Anonymous Two-Channel Closed Wave System Trilogy, Part I) (2026), Zenodo. Concept DOI: 10.5281/zenodo.21763995 (identification of the transmission rate with the elementary charge; motivation of assumption (vi) only)

**External**

[5] C. R. Menyuk, "Nonlinear pulse propagation in birefringent optical fibers," IEEE J. Quantum Electron. **23**, 174 (1987)

[6] G. P. Agrawal, *Nonlinear Fiber Optics*, 6th ed. (Academic Press, 2019)

[7] P. W. Milonni, "Semiclassical and quantum-electrodynamical approaches in nonrelativistic radiation theory," Phys. Rep. **25**, 1 (1976)

[8] P. B. Greene, L. Kofman, "Preheating of fermions," Phys. Lett. B **448**, 6 (1999); "Theory of fermionic preheating," Phys. Rev. D **62**, 123516 (2000); L. Kofman, A. Linde, A. A. Starobinsky, "Towards the theory of reheating after inflation," Phys. Rev. D **56**, 3258 (1997)

[9] This repository, folder `universal_interaction_manybody_connection_v1` (design policy, the mediated-vertex uniqueness theorem, stages 1–3, GENESIS v1–v3, the 3D demo — all scripts, result JSONs, and the commit series)

[10] E. Fermi, J. Pasta, S. Ulam (with M. Tsingou), "Studies of nonlinear problems I," Los Alamos Report LA-1940 (1955)

[11] K. Baudin *et al.*, "Classical Rayleigh-Jeans condensation of light waves," Phys. Rev. Lett. **125**, 244101 (2020)

[12] J. M. Manley, H. E. Rowe, "Some general properties of nonlinear elements," Proc. IRE **44**, 904 (1956)
