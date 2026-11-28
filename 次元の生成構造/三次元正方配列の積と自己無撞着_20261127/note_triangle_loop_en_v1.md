# Can a Universe Emerge Without Spacetime? — Running "Crystal Melting, Spatial Expansion, and Particle Birth" in Pure Python the Day After the Viral AI Quantum Gravity Paper

[Header Image: note_figs/note_header_triangle_loop.png]

On November 25, 2026, V. Chandrasekaran released a preprint on background-independent quantum gravity (arXiv:2611.34924v1) prominently featuring AI-assisted theoretical physics, drawing widespread attention not only from physicists but also from the general science community.

The core question that has captured so many people's imagination is remarkably simple:

> **Can a universe—spacetime, particles, and gravity—emerge purely from relational laws, without preparing space or time as a pre-existing container?**

I have been exploring this exact "no background spacetime" starting point together with AI. Following yesterday's preprint, I have published a rapid communication on Zenodo reporting interim results extending my earlier complex-number numerical inflation experiments to **octonions** and **3D cubic arrays**, accompanied by **all deterministic source code, raw datasets, and plotting scripts** so that anyone can verify every single claim independently.

- **Paper (Zenodo Concept DOI — always points to the latest version):**  
  https://doi.org/10.5281/zenodo.24614040
- **English PDF (direct download):**  
  https://zenodo.org/records/24614041/files/triangle_loop_product_crystal_selection_melting_rapid_en.pdf
- **Japanese PDF (direct download):**  
  https://zenodo.org/records/24614041/files/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E4%B8%89%E8%A7%92%E5%BD%A2%E3%83%AB%E3%83%BC%E3%83%97%E7%A9%8D%E3%81%AB%E3%82%88%E3%82%8B%E7%B5%90%E6%99%B6%E9%81%B8%E6%8A%9E%E3%81%A8%E8%9E%8D%E8%A7%A3_%E9%80%9F%E5%A0%B1.pdf
- **Full Code & Datasets (GitHub):**  
  https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E7%A9%8D%E3%81%A8%E8%87%AA%E5%B7%B1%E7%84%A1%E6%92%9E%E7%9D%80_20261127
- **Japanese edition of this note article:**  
  https://note.com/kiharanoriaki/n/n6a1c65eb049b

In this article, I walk through the key findings in plain language while keeping a strict, honest boundary between **what was verified by exact computation** and **what remains a physical interpretation (hypothesis)**.

---

## 1. Just One Rule: No Space, No Clock, No Randomness

Conventional physics simulations begin by introducing a 3D spatial grid, a time step t, random noise, and external potentials or operators.

Here, we remove **every single one of them**.

All we prepare is an N × N × N cubic array T of numbers, with three completely anonymous, unlabeled axes. Then we repeatedly apply one self-updating operation—the **triangle-loop product**—to itself:

　Tri(T)[i, j, k] ＝ Σ_{p, q, r} T[i, p, q] × T[q, j, r] × T[r, p, k]

Three factors pass three internal indices (p, q, r) around a closed triangular loop like a relay baton. There is no division, no randomness, and no external operator—only pure multiplication and addition.

What happens when you iterate this single rule?

---

## 2. Finding 1: The Building Blocks—Quaternions (4D) and Octonions (8D)—Select Themselves

In mathematics, the Cayley–Dickson construction doubles dimensions step by step: real numbers (1D) → complex numbers (2D) → quaternions (4D) → octonions (8D) → sedenions (16D).

We placed each algebra's multiplication table (pure integers 0, +1, −1) directly into the cubic array T and iterated Tri(T) using **exact integer arithmetic** with zero rounding error. The result was striking:

- **Complex numbers (2D):** Vanish completely to **0 on step 1**.
- **Quaternions (4D):** All 16 nonzero cells survive as an **exact period-2 crystal** (returning to themselves every 2 steps).
- **Octonions (8D):** All 64 nonzero cells survive as an **exact period-2 crystal**.
- **Sedenions (16D):** Even with zero defects, the structure destabilizes and **self-destructs to 4 cells by step 36**.

We never put "4D" or "8D" in by hand. **The triangle-loop product itself spontaneously selected only quaternions (4D) and octonions (8D) as exact crystals.**

---

## 3. Finding 2: A Single-Cell Scratch Melts the Octonion Crystal in Two Stages and Leaves Persistent "Particles"

Next, we took the perfect octonion crystal (8 × 8 × 8 ＝ 512 cells) and scratched a single cell by a tiny amount (0.0001).

When you do the same test with a standard star-contraction ternary product (the Bhattacharya–Mesner product, which has no closed loop), 100% of defects at N ≥ 8 trigger an avalanche that wipes out the entire array to zero.

Under the triangle-loop product, across **all 1,024 single-cell defects** (512 cells × real or imaginary scratch), **complete extinction occurred 0 times (100% survival)**—and the remnants were extraordinarily structured:

- **Real scratches (all 448 originally zero cells):**  
  After melting, every single case precipitates into an **8-cell complex rotation-matrix pair**. It obeys the exact analytical equation z' ＝ (z₀² ＋ z₁²) z, conserving its internal ratio z₁ / z₀ while rotating on the circle S¹ for **more than 20,000 steps without ever decaying**. Perturbations applied to the external 504 cells are blocked to exact zero in one step or neutrally deflected.
- **Imaginary scratches (448 originally zero cells):**  
  They split strictly **336 : 112 (exact 3 : 1)**. The 336 cases undergo a **two-stage cascade phase transition: octonions (8D) → quaternions (4D, 32 cells) → complex numbers (2D, 8 cells)**, saturating completely at step 100 and remaining frozen through step 20,000. Without any external constraint, they spontaneously satisfy the complex sum-of-squares condition ΣU² ＝ 0 and orthogonal Pauli-matrix projections (I ± σ_y).

---

## 4. Finding 3: Melting a 64×64×64 Crystal Produces "Spatial Expansion" and "Particle Condensation" Simultaneously

We then constructed a globally coupled large crystal of **64 × 64 × 64 ＝ 262,144 cells** (O₈ ⊗ O₈), introduced a single real scratch, and tracked it for 500 steps.

[Figure 1: results/figures_en/fig_n64_512melt_snapshots_en.png]

Figure 1. Evolution of the 64×64 projected intensity map. Ripples spread across the uniform crystal (step 30), an 8×8 grid of dark melting channels carves it into blocks (steps 40–60), and by steps 80–500 it recrystallizes into a stable 512-cell composite survivor arranged as an 8×8 lattice of 2×2 spots.

The 512-cell composite remnant stayed stable all the way through step 500, and the energy across the 64 elements separated into a crisp **three-level discrete energy hierarchy**.

This led to the most exciting question of all:

> **In a system with no coordinates, what is "distance"?**

Our operation is built out of **triangle loops**. Intuitively, a distant object subtends a **razor-thin, elongated triangle (tiny view angle)**, whereas tightly bound nearby components form a **nearly equilateral triangle**.

So we exhaustively measured the side-length ratio (**view-angle ratio ρ ＝ shortest side / longest side**, where 1 is equilateral and smaller means elongated / distant) across all 262,144 triangle loops at every step.

[Figure 2: results/figures_en/fig_n64_viewangle_vs_condensation_en.png]

Figure 2. Particle condensation (top left: E_512), exponent contrast between particle and melted background (top right), and triangle view-angle ratio ρ (bottom left). Synchronized precisely with energy condensing into the 512 cells during steps 65–75, background-side triangles (blue dashed) plunge below 10⁻³⁰⁰, while particle-side triangles (red solid) rebound in a V-shape back to equilateral geometry (ρ ~ 1).

[Figure 3: results/figures_en/fig_n64_viewangle_spectrum_en.png]

Figure 3. Evolution of the triangle view-angle spectrum at the particle representative (left). A single initial group splits sharply into 192 background-connected triangles marching steadily leftward into causal separation, and 64 particle-internal triangles condensing at the equilateral peak on the right.

Without coordinates or an inflaton field, two phenomena happened in exact lockstep:

1. **Background-connected triangles (blue):** The view-angle ratio ρ drops linearly on a log scale and plunges past the double-precision floor 10⁻³⁰⁰ by step 75—**complete causal separation via exponential stretching (inflation-like rarefaction)**.
2. **Particle-internal triangles (red):** After briefly stretching, they execute a **V-shaped rebound back to equilateral geometry (ρ ~ 1)**—**spontaneous particle condensation**.

---

## 5. Comparison with Chandrasekaran (arXiv:2611.34924v1) and Honest Limitations

V. Chandrasekaran's November 25 preprint formulates background-independent quantum gravity starting from 4D vacuum General Relativity on the null cone. A full-text comparison shows deep structural resonance: no fixed background metric, internal corner reference frames with triple gluing consistency, finite-dimensional Hilbert spaces, and subspace separation via null states.

At the same time, this rapid communication contributes distinct strengths:

- **Immediate third-party verifiability:** Whereas Chandrasekaran's preprint contains no code or numerical experiments, we release all 18 deterministic Python scripts and complete datasets that anyone can reproduce in seconds to tens of seconds.
- **Spontaneous selection of dimensions and algebra:** Without assuming 4D spacetime or a gravitational action, exact integer arithmetic shows that only 4D (quaternions) and 8D (octonions) survive as exact crystals.
- **Particle remnants, cascade melting, and bimodal view-angle splitting:** Phenomena corresponding to matter/particles and inflation—outside the scope of Chandrasekaran's vacuum-gravity paper—emerge spontaneously without an inflaton field.

For full scientific honesty, let us also state the current limitations clearly:

- Interpreting these verified numerical structures as "spatial expansion" and "particles" is a **physical hypothesis**, not yet a proven physical derivation.
- Constructing multiple spatially separated independent particles that collide and scatter remains an open problem.
- Quantitative recovery of semiclassical General Relativity or Standard Model mass ratios and coupling constants lies beyond the scope of this rapid communication.

---

## 6. Behind the Scenes: Switching AI Models Unlocked a Half-Day Breakthrough

A brief note on how this work came together:

I initially began this investigation with another AI (Claude Code), which got stuck in dead ends involving 2D matrix multiplication and a misinterpretation of "exponentiation." Switching to **Google Gemini 4 Neon** and asking it to redesign the algebraic operation from scratch unlocked everything—the triangle-loop product, exact integer proofs, analytical theorems, exhaustive 1,024-case classification, 20,000-step tracking, N ＝ 64 large-crystal view-angle measurements, and bilingual publication—in **roughly half a day**.

For full transparency, the record of the earlier AI's detours and failed attempts is also openly preserved in the repository.

---

## 7. Verify It Yourself in Seconds

With standard Python (NumPy and Matplotlib) and zero randomness, you will obtain the exact same numbers:

```bash
git clone https://github.com/WurabeSeiji/ai-chat-logs-open.git
cd ai-chat-logs-open/次元の生成構造/三次元正方配列の積と自己無撞着_20261127

# 1. Exact integer selection of C, H, O (seconds)
python3 tri_loop_product_v0.py

# 2. N=64 large-crystal melting and view-angle spectrum
python3 tri_loop_n64_triangle_viewangle_run_v0.py
python3 plot_tri_loop_n64_viewangle_v0.py
```

If you are curious whether such behavior really arises from a single anonymous multiplication rule, please run the scripts and check it for yourself. Feedback and independent verifications are warmly welcome!

---

#QuantumGravity #TheoreticalPhysics #Physics #Mathematics #Cosmology #AI #Gemini #Python #Octonions #Inflation #BackgroundIndependence #ReproducibleScience #Chandrasekaran
