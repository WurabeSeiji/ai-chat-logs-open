# Thought Experiments through Chapters 1–5 of Shimizu's Foundations of Quantum Theory — Published on Zenodo (Sequel to the Previous Paper)

I have published an observational paper on Zenodo that organizes, in academic form, a reading note developed while studying Akira Shimizu's well-regarded Japanese textbook *Foundations of Quantum Theory* (New Edition). The work treats measurement, the uncertainty principle, quantum entanglement, the algebra of observables, the ontology of physical quantities, and the finite-width structure of particles in seven stepwise thought experiments.

This paper continues from [the previous paper (Thought Experiments through Chapters 1–3)](https://note.com/kiharanoriaki/n/nbfc40cb3cfa3), carrying over its five thought experiments and adding two new ones corresponding to Chapters 4–5: the representation of physical quantities over a complex phase space, and the rereading of particles as rectangular phase-energy windows with a central phase and a finite width.

This is NOT a formal research paper. It is a reading note arranged in academic form. It does not modify the mathematical predictions of standard quantum theory; it offers a geometric and ontological interpretation of each concept as one reading. The Born rule is not denied; it is localized at the final step of a three-tier hierarchy: rectangular body → Fourier partial sum → inner-product projection onto the detection basis.

Publication information:

- Title: Thought Experiments through Chapter 5 of Foundations of Quantum Theory: From the Wall of Identification in Measurement to the Rectangular Phase-Energy Window
- v1 DOI: https://doi.org/10.5281/zenodo.20398527
- Concept DOI (always redirects to latest version): https://doi.org/10.5281/zenodo.20398526
- Zenodo record: https://zenodo.org/records/20398527
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN = 6 files
- AI peer review: Gemini × 1, Grok × 1, ChatGPT × 2 integrated

---

## The source textbook

The starting point of this thought experiment is the following textbook (currently available in Japanese only):

- Title: New Edition: Foundations of Quantum Theory — For an Easy Understanding of Its Essence (新版 量子論の基礎: その本質のやさしい理解のために)
- Author: Akira Shimizu (清水 明), Professor Emeritus, The University of Tokyo
- Publisher: Saiensu-sha (New Physics Library, Supplementary Volume 2)
- First edition: 2003
- ISBN: 978-4781910628

Purchase links (Japan):

- Amazon Japan: https://www.amazon.co.jp/dp/4781910629
- Rakuten Books: https://books.rakuten.co.jp/search?sitem=9784781910628
- Saiensu-sha (publisher direct): https://www.saiensu.co.jp/

Note for international readers: As of this writing, no official English translation is available. International readers proficient in Japanese can purchase through the links above. Amazon Japan ships internationally.

### What Chapters 4–5 cover

In Shimizu's textbook, Chapter 4 introduces a particle moving in one-dimensional space, the wave function, and the probability density. Chapter 5 treats the box potential, the infinite well, and the tunneling effect.

In the previous paper, reading Chapters 1–3 led to the picture of "the algebra of observables" — a geometric viewpoint reading observables as different projections of the same wave packet. The present paper deepens this viewpoint and asks: what is the substance of a physical quantity prior to projection, and what structure does the particle itself possess?

---

## Two newly added thought experiments

### Thought Experiment VI: Are physical quantities real numbers?

The starting question: does the fact that observed values are real numbers imply that the underlying structure of physical quantities is closed over the real-number field? Just as wave functions are complex while observed probabilities are real, can we allow the underlying generating structure to be complex even though what is observed is real?

Three assumptions:

Assumption 1: Observed physical quantities are real-valued, but this is only a constraint on the output type of the observation operation.

Assumption 2: Fundamental physical quantities are defined as quantities over a complex-valued phase space, and real-valued observed values appear as real projections of that complex-phase structure.

Assumption 3: Physically stable values are selected not as continuous arbitrary values but as discrete values satisfying a phase-closure condition across the whole system. One-dimensionally this is the standing-wave condition; in phase space, it corresponds to the Bohr–Sommerfeld / EBK quantization condition for the action integral; more generally, the integrality condition of geometric quantization.

Complex representation of the position phase: representing position x as a phase θ_x = kx, the central position x_0 corresponds to the central phase, and the spatial spread Δx corresponds to the phase width Δθ_x = k Δx.

Conclusion: The substance of a physical quantity can be read not as a single value on the real line but as a phase-closure structure on a complex phase space. Under this assumption, uncertainty, interference, entanglement, and quantization are not separate phenomena but the same structure seen through different observation projections.

### Thought Experiment VII: Are particles rectangular phase-energy windows?

We place the body of the particle as a rectangular phase window: a rectangular structure centered at the central phase θ_x, taking the value 1 in the region of full width Δθ_x and 0 elsewhere.

The observation image as a Fourier partial sum: expanding the rectangular phase window R_x(θ) as a Fourier series of period 2π and defining the low-pass operator L_Λ in the observation bandwidth as "truncation at the N-th harmonic," the observation image appears as the low-order Fourier partial sum S_N(θ). As N grows, the partial sum converges to the rectangle, but Gibbs-type ringing remains at discontinuities.

Ontological inversion:

　　"the particle has a position" → "position-phase energy appears as a particle"

The standard order is "particle → position → wave function → |ψ(x)|^2," whereas in this paper's rereading it is "finite-width position-phase-energy window → particle → observed wave form."

The definition of the particle is a definitional hypothesis: why R_x(θ) is stable, by what dynamics Δθ_x is determined, why the energy density E_0 rides on it — none of this is derived within the scope of this paper. These are the questions answered, from specific Lagrangians, by the Skyrme model, the MIT bag model, the Friedberg–Lee–Sirlin non-topological soliton, the Q-ball, etc. The rectangular window of this paper is merely an alternative starting point for the same questions.

Interaction between particles is an overlap indicator: placing two particles as rectangular phase-energy windows, the overlap region R_1(θ) · R_2(θ) can be read as a candidate interaction kernel, providing a natural entry to a unified representation of interference, tunneling, scattering, binding, and repulsion.

Correspondence with the box potential: in Shimizu's Chapter 5, the infinite well produces standing waves via boundary conditions, while the finite box barrier gives rise to reflection, transmission, and tunneling via boundary connection. The common structure across both — "finite-width region and boundary conditions" — is what this paper's rectangular phase-energy window picture focuses on.

---

## Related work

The central rereading of this paper — "particle = finite-width phase window," "observation = bandwidth limitation," "phase-space area = quantum invariant" — partially resonates with the following lineages.

de Broglie double-solution theory, Madelung hydrodynamics, Bohm pilot wave: the genealogical sources of the intuition "particle = field structure in a finite region." Colin–Durt–Willox (2017) gives the modern review of the double-solution program.

Skyrme model, MIT bag model, Q-ball: lineages that derive the finite-width structure of particles from specific Lagrangians. Skyrme (1961), Chodos et al. (1974), Coleman (1985), etc.

Gabor transform, Slepian–Pollak prolate spheroidal wave functions (PSWF), Hardy's theorem: the rigorous mathematical counterparts of the rectangular window and bandwidth limitation. The Fourier partial sum of this paper is positioned as an intuitive, elementary approximation to these rigorous frameworks.

Coherent states: the lineage of minimum-uncertainty wave packets with central phase and finite width (Gaussian version). Schrödinger (1926), Glauber (1963).

De Gosson quantum blobs and symplectic capacity: the framework for reading the phase-space area element as a quantum invariant. De Gosson (2013), de Gosson–Luef (2009), Gromov's (1985) non-squeezing theorem.

This paper presents no new Lagrangians or dynamical equations. It is an attempt, leaving the mathematical framework of Chapters 1–5 of Shimizu's textbook unchanged, to replace the ontological reading of "particle," "position," and "wave function" on top of that framework.

---

## Correspondence with Shimizu's textbook

The framework of complex Hilbert spaces and the algebra of observables introduced in Chapter 1 of Shimizu's textbook was addressed in the previous paper. Thought Experiment VI of this paper rereads the structure of complex Hilbert space corresponding to Chapter 4 from an ontological perspective as "physical quantities as quantities over a complex phase space."

Thought Experiment VII reconstructs the framework of the one-dimensional particle, the box potential, and the tunneling effect of Shimizu's Chapter 5 as the finite-width structure of the particle.

Through Chapters 1–5, the construction of Shimizu's textbook — algebra of observables → measurement → uncertainty → one-dimensional particle → box potential — can be read consistently with the seven thought experiments of this paper.

---

## Related links

The GitHub repository ai-chat-logs-open contains both verbatim records of the thought experiments (including trial-and-error dialogue) and the polished version, in both Japanese and English.

- GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/新版量子論の基礎
- Zenn article (integrated introduction of three papers, in Japanese): https://zenn.dev/noriaki_kihara/articles/quantum-theory-through-chapter5
- Zenn article (Thought Experiment 8 standalone): https://zenn.dev/noriaki_kihara/articles/are-physical-quantities-real-numbers
- Zenn article (Thought Experiment 9 standalone): https://zenn.dev/noriaki_kihara/articles/particles-and-box-potential
- Previous note article (Chapters 1–3): https://note.com/kiharanoriaki/n/nbfc40cb3cfa3
- Japanese note article: https://note.com/kiharanoriaki/n/n8ffc8e2c9123

Companion papers (Central Projection Framework):

- Paper 7 (α with 8.7 ppb precision): https://doi.org/10.5281/zenodo.19876200
- Paper 8 (structural correspondence with Wilson lattice gauge theory): https://doi.org/10.5281/zenodo.19881119
- Six-dimensional sign vector xyztRQ thought experiment: https://doi.org/10.5281/zenodo.19904714

---

## Important reservations

The present paper does NOT claim:

- Modification of the mathematical predictions of standard quantum theory
- Prediction of new physical phenomena
- Complete extension to regions beyond Chapters 1–5 — Lorentz covariance, locality of fields, positivity, unitarity, scattering amplitudes, spin statistics, gauge symmetry, particle creation and annihilation

The present paper only records:

- Seven observations obtained from reading Chapters 1–5 of Shimizu's textbook
- The structural connections among the wall of identification, indistinguishability of fluctuation locus, wavenumber representation, composite wave packet, the algebra of observables, the complex phase space, and the rectangular phase-energy window
- A verbatim record of the thought process via AI dialogue (including cycles of erroneous answer → correction)
- Clarification of the position of this paper among existing lineages (de Broglie, Skyrme, PSWF, coherent states, de Gosson quantum blobs, etc.)

---

## Recommendation: please have the textbook on hand

The arguments in this paper are written so that they can be understood independently without reading Shimizu's textbook. However, for readers who wish to study the foundations of quantum theory seriously, I strongly recommend obtaining the original textbook.

In particular, the construction of Chapter 1 — understanding the essence of quantum theory as the algebra of observables — and the construction of Chapter 5 dealing with the one-dimensional particle and the box potential have a depth not found in other textbooks. The present work is a reading note of these constructions.

Purchase links:

- Amazon Japan: https://www.amazon.co.jp/dp/4781910629
- Rakuten Books: https://books.rakuten.co.jp/search?sitem=9784781910628
- Saiensu-sha: https://www.saiensu.co.jp/

Comments and questions are welcome via the note comment section or Zenodo's comment function.

Author: Noriaki Kihara, WF System Co., Ltd., ORCID: 0009-0004-6753-4020

---

#QuantumTheory #QuantumMechanics #ShimizuTextbook #FoundationsOfQuantumTheory #BoxPotential #TunnelingEffect #Particle #RectangularPhaseEnergyWindow #FourierPartialSum #OntologicalInversion #ComplexPhaseSpace #PSWF #CoherentState #QuantumBlob #SymplecticCapacity #ThoughtExperiment #PhaseSpace #TheoreticalPhysics #MathematicalPhysics #Physics #ReadingNote #Zenodo #IndependentResearcher #AIDialogue
