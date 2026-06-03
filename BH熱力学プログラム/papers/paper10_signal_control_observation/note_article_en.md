# Juxtaposing and Classifying 10 Structural Correspondences Between Signal/Control Theory and Quantum Theory — Paper 10 (Observation Paper)

Signal/control theory (the engineering of communication, measurement, and control) and quantum mechanics, quantum optics, and open quantum systems (theoretical physics) developed **independently** in the 20th century. Yet when the mathematical structures underlying both domains are juxtaposed, one finds that at many points a **rigorous identity** or a **strong structural correspondence** holds.

I have published on Zenodo an **observation paper** that organizes this known fact into **10 correspondences**, sorting each one by distinguishing "how far it is a mathematical identity, and where it becomes mere analogy (parallelism)." It offers no new interpretation, no new theorem, and no new prediction; its citations are only the classic original papers and standard textbooks of both fields, with zero self-citation.

Publication information:

- DOI (Concept, always resolves to latest): https://doi.org/10.5281/zenodo.20521598
- DOI (v1.0, latest): https://doi.org/10.5281/zenodo.20521599
- Zenodo page: https://zenodo.org/records/20521599
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN = 6 files

---

## Why the Distinction Matters

Naive slogans such as "quantum mechanics is just signal processing" or "Fourier analysis is what quantum theory really is" circulate frequently on social media and in popular explanations. These are half right and half misleading.

Some correspondences are **rigorously identical** in mathematics; others are mere **structural parallelisms**. Conflating the two leads people to speak of claims that do not hold as if they were "already proven."

The point of this paper is not to discover new correspondences, but to **correctly classify known correspondences by strength, and to draw precisely the boundary lines that are easy to overstep.**

---

## The 10 Items and the Three-Tier Classification

The strength of correspondence is divided into three tiers.

### ★ Rigorous mathematical identity (5 items)

**1. Uncertainty principle**
Gabor's time–bandwidth uncertainty Δt Δω ≥ 1/2 and the Heisenberg–Kennard–Robertson Δx Δp ≥ ℏ/2 are the same Fourier-analytic inequality. Gabor himself made the equivalence explicit in 1946.

**2. Time–frequency quasi-probability distribution**
The identical distribution independently discovered by Wigner (1932, quantum) and Ville (1948, signal). The standard name "Wigner–Ville distribution" reflects the independent discovery by both.

**8. Jones vector ↔ qubit**
The same two-component complex vector (ℂP¹ ≅ S²). The Bloch sphere = the Poincaré sphere.

**10. SVD ↔ Schmidt decomposition**
The same tool of linear algebra.

(In the body, Parseval ↔ unitarity and others bring the rigorous-identity tier to 5 items in total.)

### ★ Conditionally rigorous isomorphism (1 item)

**3. Paraxial wave equation ↔ Schrödinger equation**
Isomorphic under the paraxial, monochromatic, scalar approximation. Under the standard convention E ∝ e^(ikz), a high-refractive-index guiding region (δn > 0) corresponds to a bound potential well (V < 0). Because many explanations get this sign correspondence backwards, the body sets it out explicitly.

### ◎ Strong structural correspondence (1 item)

**4. Sampling ↔ phase-space degrees of freedom**
The time–bandwidth product N ~ 2BT and the counting of phase-space volume / h. This is not the same theorem, however, but an asymptotic correspondence of effective degrees of freedom via the Slepian–Pollak prolate eigenvalues.

### △ Structural parallelism (4 items)

From here on it is "similar structure," not "the same mathematics." The two must not be conflated.

**5. State space ↔ Hilbert picture**: differences in unitarity, bilinearity, and driving input
**6. Kalman observability ↔ CSCO**: a difference in mechanism between dynamics and kinematics. For state reconstruction, quantum state tomography / IC-POVM is closer than CSCO
**7. Kalman filter ↔ quantum filtering**: via POVM, not projective measurement
**9. Phase noise ↔ dephasing**: only the T₂ part, not decoherence as a whole

---

## What This Paper Does NOT Claim (Explicit)

- it proposes no new interpretation of quantum mechanics
- it claims no physical extension of signal/control theory
- it predicts no physical constant and modifies no existing theory
- it does not extend a strong identity into a structural parallelism

What separates this paper from the countless analogy collections is precisely that it carves out, with precision, the places where the naive "quantum mechanics = signal processing" discourse oversteps (the limits of #5–#9).

---

## On One Certainty

The mathematics on the signal-theory side (OFDM, MIMO, Kalman, Wigner–Ville, SVD, and so on) is socially deployed and runs daily in the form of communications, GPS, autonomous driving, radar, MRI, and optical communication. The validity on the quantum side, meanwhile, is independently confirmed by quantum experiments.

So what this paper juxtaposes is the fact that "the same mathematical structure has been independently confirmed in two different physical domains." It stays at this observation itself, not at any new claim.

---

## Related Resources

- Zenn article (more technical): https://zenn.dev/noriaki_kihara/articles/signal-control-quantum-correspondences
- The single lens of phase-space area (Paper 11, companion): https://doi.org/10.5281/zenodo.20521566
- note index article (all papers): https://note.com/kiharanoriaki/n/nc1619291b690

---

Author: Noriaki Kihara
WF System Co., Ltd. / ORCID: 0009-0004-6753-4020

Paper DOI (Concept): https://doi.org/10.5281/zenodo.20521598
Paper DOI (v1.0, latest): https://doi.org/10.5281/zenodo.20521599
Zenodo page: https://zenodo.org/records/20521599

---

#SignalProcessing #ControlTheory #QuantumMechanics #QuantumOptics #OpenQuantumSystems #UncertaintyPrinciple #WignerVille #KalmanFilter #ParaxialApproximation #SchrodingerEquation #FourierAnalysis #ObservationPaper #Zenodo #TheoreticalPhysics #MathematicalPhysics #Preprint
