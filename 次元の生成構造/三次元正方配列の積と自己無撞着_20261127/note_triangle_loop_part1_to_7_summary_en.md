# I Got This Far With Almost Zero Free Parameters — The Seven-Part Triangle-Loop Product Series (Parts I–VII), Which Began From the Frustration "I Still Can't Read the Table of Spacetime," and the Open Problems for the Next Phase

Let me start with my honest impression.

**Nobody is more surprised than I am, and I am the one building this.**

This morning I published a unified paper (Part I) and a connection hypothesis (Part II). They report that the skeleton of quantum structure, gauge symmetry, gravity, and fermions comes out of a single one-line computational rule: take a three-dimensional cube of real numbers (a rank-3 array T), multiply it by itself around a triangular loop, and normalize. I call this the triangle-loop product Tri(T).

Right after that, however, a strong unease stayed with me.

**"The equations and the symmetries came out, but I still cannot properly read the table of spacetime."**

Which part of this three-dimensional array is time t? Which part is three-dimensional space (x, y, z)? Which part is internal space (electric charge and color charge)? And why does our world look like a smooth 3+1-dimensional spacetime to an observer? When I wrote the first paper, I had not fully understood these points.

I wrote Part III to resolve that frustration, and then kept going, working through the remaining problems one at a time. Before I knew it, **seven papers, Part I through Part VII**, had come together.

When I lined up all seven and looked back, I noticed something I had not expected at all.

---

## I Thought I Would Need Free Parameters Much Sooner

When you build a theory in physics, you normally need several constants put in by hand (free parameters) at an early stage.

- **Planck's constant ħ**, to make the theory quantum
- **The fine-structure constant α ≈ 1/137.036** (the elementary charge), which sets the strength of electromagnetism
- The coupling constants of the weak and strong forces, the particle masses, the number of generations, the direction in which parity is broken, and so on

The Standard Model alone takes about 20 free parameters as manual input. So as I went through this work, I expected that somewhere around Part III or Part IV I would be unable to go further without bringing in some physical constant or assumption from outside.

Yet **even at Part VII, I have not used Planck's constant ħ a single time.**
I have not used the fine-structure constant α either. The gauge group, parity violation, the number three for the generations, the dimension of spacetime: none of these was assumed from outside.

All I use is a unit convention equivalent to c = 1 (or G = 1) in geometric units. That is not a physical free parameter. It is only a conversion of the units on the ruler.

Running **one dimensionless, parameter-free line**, the triangle-loop product Tri(T) on a real three-dimensional array followed by normalization, in a deterministic program with no random numbers was enough for the main skeleton of modern physics to rise in this much detail.

---

## The View That Opened Up Because I Did Not Assume Continuous Spacetime (A Contrast With the Yang–Mills Problem)

In November 2026 there was widely discussed news: an OpenAI system generated a Lean formal proof of more than 1,000 pages for the existence of pure Yang–Mills theory, one of the million-dollar Millennium Prize Problems, and a committee led by Nobel laureate David Gross examined it.

Why did the existence proof and the mass gap of Yang–Mills theory stay unsolved for 50 years? The main reason is that conventional quantum field theory imposes the **continuum assumption** as a rule: spacetime is an infinitely fine continuum ℝ⁴. A theory that runs consistently on a finite lattice meets an enormous mathematical wall when you take the lattice spacing a → 0. The energy diverges to infinity, or the interaction vanishes altogether (triviality).

In this series, by contrast, **the assumption of a continuous spacetime ℝ⁴ is absent from the very beginning.**

- On the side of what exists, there is only the finite, discrete Tri(T) structure itself.
- It looks like a smooth continuous spacetime **because the observer's resolution is finite** (the resolution gauge and logarithmic transform formulated in Part III).

From a position that does not assume a continuum, there are no divergences or infinite limits to struggle with. A one-line discrete computation gave, quite directly, **the non-Abelian gauge symmetry SU(3)_C × SU(2)_L × U(1)_Y, massless long-range propagation of the photon versus non-Abelian confinement of gluons, W/Z mass generation by Higgs condensation, and the bare unification ratio** (1 : 1 : 1, with the bare Weinberg angle sin²θ_W = 3/8).

(A note for clarity: the OpenAI proof was carried out in the arena of pure mathematics, the Wightman axioms on the continuum ℝ⁴, and its premises differ from those of this discrete approach. I do not claim to have solved the Millennium Prize Problem. Still, as a way of asking again whether we are confusing a continuum approximation used for computation with physical observation, the contrast is a very symbolic one.)

---

## What Happened in Parts I–VII, Still With No Parameters

All papers are published on Zenodo (CC BY 4.0) together with deterministic Python code that uses no random numbers, execution logs, and PDFs in English and Japanese.

1. **Part I: Spontaneous Emergence of Quantum, Gauge, Gravitational, and Fermionic Structures** (the unified paper)
   - Iterating Tri(T) automatically selects the octonions O₈ and the time crystal M₈, and quantum interference, SU(3) × U(1), local gravitational distortion, and an event horizon emerge.
   - Concept DOI: https://doi.org/10.5281/zenodo.24622200
   - English PDF direct download: https://zenodo.org/records/24622201/files/triangle_loop_emergent_quantum_gauge_gravity_fermion_unified_full_paper_en.pdf?download=1
2. **Part II: Unified Connection Hypothesis for Gravity, Electromagnetism, QCD, and Dirac Matter**
   - A theoretical map that connects each sector of Part I to the four forces and the matter fields in a single line.
   - Concept DOI: https://doi.org/10.5281/zenodo.24622565
   - English PDF direct download: https://zenodo.org/records/24622566/files/triangle_loop_unified_gravity_em_qcd_dirac_hypothesis_part2_en.pdf?download=1
3. **Part III: Spacetime and Internal-Space Readout Connection Hypothesis** (deciphering the table of spacetime)
   - Resolves the original frustration about the table of spacetime. Formulates the 1+3+4 split of the outer eight dimensions, the past light cone obtained by a logarithmic transform, and the resolution gauge.
   - Concept DOI: https://doi.org/10.5281/zenodo.24625962
   - English PDF direct download: https://zenodo.org/records/24625963/files/triangle_loop_spacetime_readout_connection_hypothesis_part3_en.pdf?download=1
4. **Part IV: Electroweak Symmetry Breaking and Maximal Parity Violation**
   - From G₂ ⊃ SO(4), **an exact 100% : 0% chiral splitting** (maximal parity violation) appears automatically. Demonstrates Higgs mass generation by the Tri(T) dynamics and, with zero tuning parameters, **the bare Weinberg angle sin²θ_W = 3/8** (m_W² / m_Z² = 5/8).
   - Concept DOI: https://doi.org/10.5281/zenodo.24626592
   - English PDF direct download: https://zenodo.org/records/24626593/files/triangle_loop_electroweak_symmetry_breaking_parity_part4_en.pdf?download=1
5. **Part V: Multi-Shell Large-Lattice Gravity, Causal Propagation, and the Spin-2 Unlock of g₀₀(d)**
   - The discovery of an exact closure makes it possible to compute N = 32,768 (four shells), which would take about 280 TB as a dense array, in 9.5 seconds with essentially zero approximation. In addition to a weak-coupling 1/5 decay (1/r type after the logarithmic transform) and a light-cone front of finite speed, it finds a selection rule: **a scalar mass distorts only the spatial components gᵢⱼ, and a spin-2 degree of freedom σₓ is indispensable for unlocking the weak-coupling time delay g₀₀(d) by 12 orders of magnitude.**
   - Concept DOI: https://doi.org/10.5281/zenodo.24627519
   - English PDF direct download: https://zenodo.org/records/24627520/files/triangle_loop_multishell_gravity_causal_part5_en.pdf?download=1
6. **Part VI: Two-Particle Dynamical Scattering, Gauge-Mediated Interactions, and Dynamical Violation of the Spatially Separated CHSH Bell Inequality**
   - Demonstrates collapse under the summed gravity of two particles, local U(1) gauge invariance versus charge splitting mediated by gauge links, and remote spin-2 scattering with causal delay (an exact square law). It further observes that **without ħ, two spatially separated sites spontaneously reach and self-repair the Tsirelson bound S_CHSH = 2√2, and that this is cut off discontinuously, 2√2 → 0, at the moment the event horizon collapses.**
   - Concept DOI: https://doi.org/10.5281/zenodo.24628297
   - English PDF direct download: https://zenodo.org/records/24628298/files/triangle_loop_two_particle_scattering_chsh_part6_en.pdf?download=1
7. **Part VII: The Exceptional Jordan Algebra J₃(O) Crystal, Three-Generation Triality, and Geometric Unification of the Bare Coupling Constants**
   - Within the 27-dimensional J₃(O), it finds that **the off-diagonal three-generation cyclic product J_off (24 dimensions) is an exact crystal with zero error, Tri(J_off) = −(1/64) J_off.** From F₄ ⊃ SO(8) ⊃ G₂ ⊃ SU(3)_C it proves three-generation triality (8_v ⊕ 8_s ⊕ 8_c) and the parameter-free bare gauge unification **α₃ : α₂ : α₁ = 1 : 1 : 1**.
   - Concept DOI: https://doi.org/10.5281/zenodo.24628776
   - English PDF direct download: https://zenodo.org/records/24628777/files/triangle_loop_jordan_three_generations_part7_en.pdf?download=1

---

## Why I Stop Here for Now — I Have Reached the Boundary of Free Parameters

Through Parts I–VII, **the skeleton that is uniquely fixed by the parameter-free law Tri(T) alone** has been almost fully excavated.

Put the other way, the problems that remain from here run straight into **the region that normally cannot be solved without bringing in some free parameter or initial condition**: the elementary charge (α ≈ 1/137.036), the mass couplings (the precise mass ratios of the electron, muon, tau, and the quarks), and so on.

I do not want to cloud the seven papers so far by forcing numbers to fit or by piling on hypotheses. So I close the first phase cleanly here and set out the remaining problems.

---

## Open Problems for the Next Phase

**Problem 1: The absolute values and RG running of the elementary charge α ≈ 1/137.036 and the strong coupling α_s**

- Reached so far: The bare geometric unification ratio (α₃ : α₂ : α₁ = 1 : 1 : 1, bare Weinberg angle sin²θ_W = 3/8) is fixed with no parameters.
- The wall: A first-principles derivation of the absolute values and the scale running that lead to the measured low-energy values (α ≈ 1/137.036, sin²θ_W ≈ 0.2312).

**Problem 2: Precise mass ratios of the three fermion generations and the Yukawa couplings**

- Reached so far: A two-stage breaking on the J_off crystal has been confirmed to amplify an exponential hierarchy of two to three orders of magnitude.
- The wall: The specific ratios and the Koide ratio Q = 2/3 depend on the initial perturbation ratio, and no mechanism that locks them automatically to universal values has been found.

**Problem 3: CKM / PMNS generation mixing angles and CP violation**

- Reached so far: The three-generation triality of J₃(O) (8_v ⊕ 8_s ⊕ 8_c) is established.
- The wall: Deriving the inter-generation mixing and the CP phase from the coupling with the external spacetime and electroweak sectors has not been started.

**Problem 4: Mapping the multi-shell K ≥ 3 continuum limit to the Einstein equations, and the two-particle scattering cross section**

- Reached so far: The K = 4 weak-coupling 1/5 decay, the strong-coupling linear gradient, the 12-order unlock of g₀₀(d) by spin 2 at K = 2, and two-particle scattering at N = 64 have been demonstrated.
- The wall: A direct mapping to the Einstein equations G_μν = 8πG T_μν at K ≥ 3, and the distance dependence of the scattering cross section.

**Problem 5: Measuring multi-shell color confinement in pure SU(3)_C (the distance-dependent potential)**

- Reached so far: The non-Abelian confinement structure in a single shell has been confirmed.
- The wall: A direct measurement of how the energy difference between the color singlet and the color octet behaves with distance d for multiple shells, K ≥ 2, is the next step.

### A Hypothesis for the Next Breakthrough: E₈ Symmetry × Entrainment in the Kuramoto Model

Problems 1 and 2 ask whether the absolute values of the elementary charge α and the mass ratios can be produced with no free parameters. Pure linear symmetry (the representation theory of Lie groups) fixes only ratios, and it does not fix absolute values.

As the hypothesis for the next phase, I am considering the combination of the following two.

1. **E₈ symmetry** (248 dimensions): The F₄ (52 dimensions) of J₃(O) reached in Part VII lies mathematically on the way along F₄ ⊂ E₆ ⊂ E₇ ⊂ E₈. I bring in E₈ as the top-level symmetry that holds spacetime, gauge structure, and the three generations in a single vessel.
2. **Entrainment and synchronization theory of the Kuramoto model**: Entrainment in nonlinear coupled oscillators absorbs differences in initial conditions and makes the whole system converge automatically to a particular order parameter or frequency-locking ratio (a dimensionless absolute value). In fact, the rigid synchronization lock of the background clock at 10⁻¹³ rad seen in Part V, and the self-repairing attractor in Part VI that pulled a product state back to 2√2 in three steps, are already signs that a Kuramoto-type entrainment mechanism is at work inside Tri(T).

If Kuramoto-type entrainment works on the vessel of E₈, and α and the mass ratios fall out uniquely as the convergence values of a nonlinear attractor, not as free parameters, that will be the next big challenge.

For now, I place here all seven papers of Parts I–VII, reached with almost zero free parameters, together with the reproduction code, and mark the end of the first phase.

---

- **All reproduction code and paper sources** (GitHub): https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E7%A9%8D%E3%81%A8%E8%87%AA%E5%B7%B1%E7%84%A1%E6%92%9E%E7%9D%80_20261127
- **List of all papers** (PAPERS): https://wurabeseiji.github.io/ai-chat-logs-open/

#Physics #Mathematics #TheoreticalPhysics #QuantumMechanics #Gravity #ParticlePhysics #StandardModel #GaugeTheory #Octonions #YangMills #UnifiedTheory #Zenodo
