# Paper 1 §2 Review of Existing BH Thermodynamics — Draft (English)

---

## §2. Review of Existing Black Hole Thermodynamics

This section reviews the established results in black hole thermodynamics that the verification program will eventually need to either recover or contact. The review is intentionally compact and neutral. It is organized along the principal lines of historical development, and selected to highlight precisely those results that bear on the hypotheses introduced in §3.

### §2.1 The Schwarzschild metric and the classical description

The simplest non-trivial vacuum solution of the Einstein field equations is the Schwarzschild solution, derived in 1916. In standard $(t, r, \theta, \phi)$ coordinates,
$$ds^2 = -f(r)\, dt^2 + f(r)^{-1}\, dr^2 + r^2 \, d\Omega_2^2, \qquad f(r) = 1 - \frac{2GM}{r}.$$
The function $f(r)$ vanishes at $r_s = 2GM$, defining the Schwarzschild horizon. This horizon plays no special role in the classical mechanics of test particles — geodesics traverse it smoothly in the appropriate coordinate frames — but it acquires central importance in the thermodynamic and quantum-mechanical descriptions reviewed below.

### §2.2 Hawking radiation: thermal evaporation

In 1975, Hawking applied the formalism of quantum field theory in curved spacetime to the Schwarzschild background and derived a striking result: a Schwarzschild black hole emits thermal radiation at temperature
$$T_H = \frac{\hbar \kappa}{2\pi k_B c} = \frac{\hbar c^3}{8\pi G M k_B},$$
where $\kappa$ is the surface gravity at the horizon. The radiation has a Planckian spectrum and carries energy away from the black hole, leading to a finite evaporation lifetime $\tau \sim G^2 M^3 / (\hbar c^4)$.

The temperature is inversely proportional to mass: smaller black holes are hotter and evaporate faster. The endpoint of evaporation — the fate of a Planck-mass black hole — remains an open problem.

### §2.3 Bekenstein entropy

In 1973, Bekenstein argued on thermodynamic grounds that a black hole must possess an entropy proportional to the area of its horizon. Combining his proposal with Hawking's temperature gives the Bekenstein–Hawking entropy
$$S_{BH} = \frac{k_B c^3 A}{4 G \hbar} = \frac{A}{4 \ell_P^2},$$
where $A = 4\pi r_s^2$ is the area of the Schwarzschild horizon and $\ell_P = (G\hbar/c^3)^{1/2}$ is the Planck length. The entropy scales as the square of the mass.

The microscopic origin of $S_{BH}$ — what counts as a black hole microstate — is the principal open problem in the field. The number of microstates required is enormous: a solar-mass black hole has $S_{BH}/k_B \sim 10^{77}$.

### §2.4 The four laws of black hole thermodynamics

Bardeen, Carter, and Hawking (1973) established a complete thermodynamic analogy:

- **Zeroth law.** Surface gravity $\kappa$ is constant on a stationary horizon (analogue of temperature uniformity in equilibrium).
- **First law.** $dM = (\kappa / 8\pi G)\, dA + \Omega\, dJ + \Phi\, dQ$ (analogue of $dU = T\, dS - p\, dV$).
- **Second law.** Horizon area is non-decreasing in classical processes; in quantum processes, the generalized entropy $S_{\mathrm{gen}} = S_{BH} + S_{\mathrm{matter}}$ is non-decreasing.
- **Third law.** A horizon with $\kappa = 0$ (extremal) cannot be reached by any finite physical process.

These laws are presently understood as a robust thermodynamic structure, valid in classical and semi-classical settings.

### §2.5 The holographic principle

't Hooft (1993) and Susskind (1995) proposed the holographic principle: the number of independent quantum-mechanical degrees of freedom inside a region of space is bounded by the area of the boundary in Planck units, not by the volume. The principle is motivated by the Bekenstein–Hawking entropy: if a black hole's entropy is governed by its horizon area, then the maximum entropy that can fit inside any region is also area-bounded.

The principle suggests that a fundamental description of gravitational physics in $D$ spacetime dimensions has degrees of freedom living on a $(D-1)$-dimensional surface. It is the conceptual root of the AdS/CFT correspondence reviewed below, but it is logically more general.

### §2.6 The AdS/CFT correspondence

Maldacena (1997) proposed that Type IIB string theory on $\mathrm{AdS}_5 \times S^5$ is equivalent to $\mathcal{N} = 4$ super Yang–Mills theory on the boundary of $\mathrm{AdS}_5$. The correspondence is a concrete realization of the holographic principle: bulk gravitational physics is encoded in a boundary gauge theory.

In the context of black hole thermodynamics, AdS/CFT gives explicit microstate counts (in the BPS sector) and dual descriptions of black hole horizons in terms of thermal CFT states. The correspondence is the most successful current implementation of the holographic principle.

### §2.7 Higher-dimensional black holes

Tangherlini (1963) generalized the Schwarzschild solution to arbitrary spacetime dimension $D$. In $D = d + 1$ dimensions, the Schwarzschild–Tangherlini metric is
$$ds^2 = -f_T(r)\, dt^2 + f_T(r)^{-1}\, dr^2 + r^2\, d\Omega_{d-1}^2, \qquad f_T(r) = 1 - \frac{\mu}{r^{d-2}},$$
with $\mu$ proportional to mass. The horizon is at $r_h = \mu^{1/(d-2)}$, the surface area of the horizon is the volume of $S^{d-1}(r_h)$, and the Bekenstein–Hawking entropy formula generalizes to $S_{BH} = A_{d-1}/(4\ell_P^{d-1})$.

For our purposes the relevant case is $D = 5$ ($d = 4$), where
$$f_T(r) = 1 - \frac{8GM}{3\pi r^2}, \qquad r_h^2 = \frac{8GM}{3\pi}, \qquad T_H = \frac{1}{2\pi r_h}, \qquad S_{BH} = \frac{2\pi^2 r_h^3}{4 G_5}.$$
This is the natural comparison object for the program developed here.

Myers and Perry (1986) generalized the solution to include rotation; the Myers–Perry solutions are the higher-dimensional analogues of Kerr.

### §2.8 Quantum gravity programs

Several programs treat spacetime as fundamentally discrete at the Planck scale.

*Loop quantum gravity* (Rovelli, Smolin, and others, 1990s) constructs a Hilbert space for canonical quantum gravity on which area and volume operators have discrete spectra with a minimum unit at the Planck scale. The Bekenstein–Hawking entropy is recovered in this framework via the counting of horizon spin-network states.

*Causal dynamical triangulations* (Ambjørn, Jurkiewicz, Loll, 1998 onwards) implements a non-perturbative path integral over discrete causal structures, recovering an emergent four-dimensional spacetime in continuum simulations.

*Causal set theory* (Bombelli, Lee, Meyer, Sorkin, 1987) replaces the smooth spacetime manifold with a discrete partial order representing causal relations between Planck-scale events.

These programs share Premise C of §3.1 (Planck-scale discreteness) but differ widely in their specific implementations. The verification program of this paper does not commit to any particular implementation; it requires only the abstract premise that discreteness is meaningful at the relevant scale.

### §2.9 Synthesis

The established body of black hole thermodynamics has the following structural features:

- An entropy proportional to a geometric quantity (Premise A).
- A semi-classical evaporation mechanism whose origin is local to the horizon (Premise B).
- A range of programs invoking Planck-scale discreteness as a fundamental ingredient (Premise C).

These three premises are common ground; the differences among the programs concern the specific implementation of each. The verification program of this paper extracts these common premises and pursues a particular geometric and combinatorial implementation — central projection from a higher-dimensional structure, and integer-lattice packing of unit cubes — that has not, to our knowledge, been pursued in the literature.

---

**Style notes**:
- 中立的レビュー、各小節は事実記述に徹する
- 自己引用ゼロ
- §2.7 で 4+1 次元 Tangherlini を強調（プログラムの自然な比較対象）
- §2.9 で前提 A, B, C への接続を明示

**未確定事項（後で検討）**:
- 各小節の参考文献の正確な書誌情報
- §2.5 と §2.6 の長さ調整（ホログラフィーは仮定しないので簡潔に）
- §2.8 での量子重力プログラムの選定範囲
