# Thought Experiments through Chapters 1–3 of Shimizu's Foundations of Quantum Theory — Published on Zenodo

I have published an observational paper on Zenodo that organizes, in academic form, a reading note developed while studying Akira Shimizu's well-regarded Japanese textbook *Foundations of Quantum Theory* (New Edition). The work treats measurement, the uncertainty principle, quantum entanglement, and the algebra of observables in five stepwise thought experiments.

This is NOT a formal research paper. It is a reading note arranged in academic form. It does not modify the mathematical predictions of standard quantum theory; it offers a geometric interpretation of each concept as one reading.

Publication information:

- Title: Quantum Theory as the Algebra of Observables: Thought Experiments on the First Three Chapters of Shimizu's Foundations of Quantum Theory
- DOI: https://doi.org/10.5281/zenodo.20391523
- Zenodo record: https://zenodo.org/records/20391523
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN = 6 files

---

## The source textbook

The starting point of this thought experiment is the following textbook (currently available in Japanese only):

- Title: New Edition: Foundations of Quantum Theory — For an Easy Understanding of Its Essence (新版 量子論の基礎: その本質のやさしい理解のために)
- Author: Akira Shimizu (清水 明), Professor Emeritus, The University of Tokyo
- Publisher: Saiensu-sha (New Physics Library, Supplementary Volume 2)
- First edition: 2003 (Original 2003, New edition 2003)
- ISBN: 978-4781910628

Purchase links (Japan):

- Amazon Japan: https://www.amazon.co.jp/dp/4781910629
- Rakuten Books: https://books.rakuten.co.jp/search?sitem=9784781910628
- Saiensu-sha (publisher direct): https://www.saiensu.co.jp/

Note for international readers: As of this writing, no official English translation is available. International readers proficient in Japanese can purchase through the links above. Amazon Japan ships internationally.

### Why this textbook was chosen

Shimizu's textbook adopts a distinctive structure as a quantum theory textbook.

Chapter 1 develops, in abstract form, complex Hilbert spaces, eigenvalues, eigenvectors, and the algebra of observables. Chapter 2 onward introduces probability amplitudes, the Born rule, the Schrödinger equation, and the uncertainty relation in sequence.

The construction differs from the standard order of physics education (classical mechanics → Schrödinger equation → operatorization). The abstract algebraic discussion of the first chapter does not directly connect to the wave equations in later chapters at first glance, which can puzzle readers.

Through this thought experiment, however, a reading emerged in which this construction presents the essence of quantum theory at the outset — namely, the perspective of understanding quantum theory as the algebra of observables. The present work is a record of that reading.

---

## Five thought experiments

### Thought Experiment I: The identification wall in measurement

Consider a physical quantity A taking real values, measured by a measure based on the concept of distance (returning only discrete values with minimum spacing Δ). We measure the distance L between two points A1 and A2.

- For L ≥ Δ, regardless of the absolute value of Δ, infinite measurements give σ_N = Δ/√(2N) → 0 and convergence to the true value
- For L < Δ, the region is not buried in measurement error; the concept of "two points" itself is undefined

Conclusion: When a measure based on distance is used to measure distance, measurement at scales below the measure's spacing Δ is undefined prior to any invocation of the uncertainty principle. Δ functions not as a resolution but as the threshold of individuation.

### Thought Experiment II: Indistinguishability of the locus of fluctuation

Considering the reverse case (infinite-precision measure + fluctuation on the physical quantity side), as a Bernoulli process the measure-side Δ and the quantity-side Δ produce the same probability distribution observationally.

Conclusion: Whether the observed fluctuation originates from the measuring instrument side or the quantity side cannot be determined in principle from observed data alone. Even with two-sided fluctuation (δ on the quantity, σ on the measure), separating δ and σ requires assumptions external to the observation.

### Thought Experiment III: Wavelength representation of momentum and the uncertainty relation

Using the de Broglie relation p = ℏk, the uncertainty relation Δx Δp ≥ ℏ/2 can be rewritten in the wavenumber representation:

　　Δx · Δk ≥ 1/2

Conclusion: Planck's constant ℏ disappears from both sides of the inequality. What remains is a mathematical inequality of Fourier transforms concerning the product of position and wavenumber. Planck's constant appears as a unit-conversion factor.

### Thought Experiment IV: Quantum correlations as a composite wave packet

A two-particle entangled state can be described as "a state in which one composite wave packet is localized in two regions in space." The observed correlation becomes a geometric consequence of the conservation laws of the entire packet (conservation of an area element in phase space).

String analogy: When a long taut string undergoing baseline vibration has one end suddenly fixed, a soliton-like deformation appears at the opposite end. This is a result of the entire string obeying conservation laws as one system; no information has propagated from one end to the other.

Conclusion: This description does not require information transfer between spatially separated systems. What is conserved is not individual wavelengths but the product of the spreads of conjugate quantities (an area element in phase space).

### Thought Experiment V: The algebra of observables

Pairs of physical quantities satisfying an uncertainty relation ΔA ΔB ≥ |⟨[A,B]⟩|/2 — position and momentum, the three spin axes, polarization, angular momentum, energy and time — all share a common structure.

- Eigenvector: wave packet (one area element in phase space)
- Operator: projection extracting a real value from the area element
- Eigenvalue: result of the projection
- Non-commutativity: two projections viewing the same area element from different directions

Conclusion: Physical quantities satisfying an uncertainty relation are described as different projections of the same wave packet (eigenvector). The product of non-commuting projections is bounded below by the area element of the wave packet (Robertson inequality).

---

## Correspondence with Shimizu's textbook

The framework of complex Hilbert spaces, eigenvalues, eigenvectors, and the algebra of observables introduced in Chapter 1 of Shimizu's textbook can be read consistently with the conclusions of this thought experiment (the geometry of area elements and projections in phase space).

The Born rule, Schrödinger equation, and uncertainty relation introduced in Chapter 2 and beyond are positioned as individual topics developed on top of this algebraic structure.

For readers who have felt puzzled by Shimizu's construction, this thought experiment may serve as a clue to one possible reading.

---

## Related links

The GitHub repository ai-chat-logs-open contains both verbatim records of the thought experiments (including trial-and-error dialogue) and the polished version, in both Japanese and English.

- GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/新版量子論の基礎
- Zenn article (integrated introduction of three papers, in Japanese): https://zenn.dev/noriaki_kihara/articles/quantum-theory-algebra-of-observables
- Japanese note article: https://note.com/kiharanoriaki/n/n2410d4863565

Companion papers (Central Projection Framework):

- Paper 7 (α with 8.7 ppb precision): https://doi.org/10.5281/zenodo.19876200
- Paper 8 (structural correspondence with Wilson lattice gauge theory): https://doi.org/10.5281/zenodo.19881119
- Six-dimensional sign vector xyztRQ thought experiment: https://doi.org/10.5281/zenodo.19904714

---

## Important reservations

The present paper does NOT claim:

- Modification of the mathematical predictions of standard quantum theory
- A complete reformulation of quantum entanglement, spin, or statistics
- Prediction of new physical phenomena

The present paper only records:

- Five observations obtained from reading Chapters 1–3 of Shimizu's textbook
- The structural connections among the identification wall, indistinguishability of the locus of fluctuation, wavenumber representation, composite wave packet, and the algebra of observables
- A verbatim record of the thought process via AI dialogue (including cycles of erroneous answer → correction)

---

## Recommendation: please have the textbook on hand

The arguments in this paper are written so that they can be understood independently without reading Shimizu's textbook. However, for readers who wish to study the foundations of quantum theory seriously, I strongly recommend obtaining the original textbook.

In particular, the construction of Chapter 1 — understanding the essence of quantum theory as the algebra of observables — has a depth not found in other textbooks. The present work is a reading note of that construction.

Purchase links:

- Amazon Japan: https://www.amazon.co.jp/dp/4781910629
- Rakuten Books: https://books.rakuten.co.jp/search?sitem=9784781910628
- Saiensu-sha: https://www.saiensu.co.jp/

Comments and questions are welcome via the note comment section or Zenodo's comment function.

Author: Noriaki Kihara, WF System Co., Ltd., ORCID: 0009-0004-6753-4020

---

#QuantumTheory #QuantumMechanics #ShimizuTextbook #FoundationsOfQuantumTheory #UncertaintyPrinciple #AlgebraOfObservables #ThoughtExperiment #QuantumEntanglement #WavePacket #PhaseSpace #FourierTransform #RobertsonInequality #HilbertSpace #deBroglieRelation #TheoreticalPhysics #MathematicalPhysics #Physics #ReadingNote #Textbook #Zenodo #Preprint #WorkingPaper #IndependentResearcher #AIDialogue
