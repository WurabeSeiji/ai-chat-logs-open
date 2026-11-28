# Can a Universe Emerge Without Spacetime? — Testing Crystal Melting, Spatial Expansion, and Particle Condensation in Pure Python the Day After the Viral AI Quantum Gravity Paper

[Header Image: note_figs/note_header_triangle_loop.png]

On November 25, 2026, V. Chandrasekaran released a preprint on background-independent quantum gravity (arXiv:2611.34924v1) featuring AI-assisted theoretical physics, sparking intense interest among both researchers and the broader public.

The question capturing everyone's imagination is remarkably simple:

> **Can spacetime, particles, and gravity emerge purely from relational rules—without assuming space or time as a pre-existing container?**

I have been exploring precisely this background-free question alongside AI. In response to yesterday's release, I have published a rapid communication on Zenodo reporting interim results extending my earlier complex-number numerical inflation experiments to **octonions** and **3D cubic arrays**—complete with **every deterministic script, raw dataset, and plotting program** so anyone can verify the findings independently.

- **Paper (Zenodo Concept DOI — always points to latest version):**  
  https://doi.org/10.5281/zenodo.24614040
- **English PDF (direct download):**  
  https://zenodo.org/records/24614041/files/triangle_loop_product_crystal_selection_melting_rapid_en.pdf
- **Japanese PDF (direct download):**  
  https://zenodo.org/records/24614041/files/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E4%B8%89%E8%A7%92%E5%BD%A2%E3%83%AB%E3%83%BC%E3%83%97%E7%A9%8D%E3%81%AB%E3%82%88%E3%82%8B%E7%B5%90%E6%99%B6%E9%81%B8%E6%8A%9E%E3%81%A8%E8%9E%8D%E8%A7%A3_%E9%80%9F%E5%A0%B1.pdf
- **Full Code & Datasets (GitHub):**  
  https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E4%B8%89%E6%AC%A1%E5%85%83%E6%AD%A3%E6%96%B9%E9%85%8D%E5%88%97%E3%81%AE%E7%A9%8D%E3%81%A8%E8%87%AA%E5%B7%B1%E7%84%A1%E6%92%9E%E7%9D%80_20261127

In this post, I summarize what happened—strictly separating **verified computational facts** from **physical interpretations (hypotheses)**.

---

## 1. Just One Rule: No Space, No Clock, No Randomness

Conventional physics simulations start by setting up a 3D spatial grid, a time variable t, random noise, and external potentials. Here we strip all of them away.

We start with a single N × N × N cubic array T with unlabeled, completely symmetric axes, and iterate one self-updating operation—the **triangle-loop product**:

　Tri(T)[i, j, k] = Σ_{p, q, r} T[i, p, q] × T[q, j, r] × T[r, p, k]

Three factors pass three internal indices (p, q, r) around a closed triangular loop. No division, no randomness, no external operators—only multiplication and addition.

---

## 2. Finding 1: Only Quaternions (4D) and Octonions (8D) Survive as Exact Crystals

Feeding the exact integer multiplication tables (entries 0, +1, −1) of the Cayley–Dickson algebras into Tri(T) with zero floating-point error yields a stark selection rule:

- **Complex numbers (2D):** Vanish completely to 0 on step 1.
- **Quaternions (4D):** All 16 entries survive as an **exact period-2 crystal** (U² = I).
- **Octonions (8D):** All 64 entries survive as an **exact period-2 crystal** (U² = I).
- **Sedenions (16D):** Self-destruct to 4 cells by step 36 even without any defect.

We did not choose 4D or 8D by hand—the triangle-loop product selected them spontaneously.

---

## 3. Finding 2: A Single-Cell Defect Triggers an 8D → 4D → 2D Cascade and Persistent "Particles"

Adding a tiny perturbation (0.0001) to a single cell of the octonion crystal (8³ = 512 cells × real/imaginary = **1,024 exhaustive cases**) produced **zero complete extinctions (100% survival)**—in sharp contrast to star-contraction (BM) products, which suffer 100% zero extinction at N ≥ 8.

- **Real defects (all 448 zero cells):** Precipitate into an **8-cell complex rotation-matrix pair** obeying the exact closed equation z' = (z₀² + z₁²) z. It conserves the internal ratio z₁/z₀, rotates on S¹ for **over 20,000 steps without decaying**, and shields itself against perturbations from the other 504 cells.
- **Imaginary defects (448 zero cells):** Split strictly **336 : 112 (exact 3 : 1)** into a 32-cell quaternionic intermediate phase and an 8-cell complex fixed point. Tracking to 20,000 steps revealed a **two-stage cascade transition (8D → 4D → 2D)** that saturates completely by step 100 and remains invariant through step 20,000—spontaneously satisfying ΣU² = 0 and Pauli projection structure (I ± σ_y) without any external constraint.

---

## 4. Finding 3: In a 64×64×64 Crystal, Background Expansion and Particle Condensation Emerge Together

Next, we built a globally coupled large crystal of **64 × 64 × 64 = 262,144 cells** (O₈ ⊗ O₈) and introduced a single real defect.

[Figure 1: results/figures_en/fig_n64_512melt_snapshots_en.png]

Figure 1. Evolution of the 64×64 projected intensity map. An 8×8 grid of dark melting channels carves the crystal (steps 40–60), leaving a stable 512-cell composite survivor arranged as an 8×8 lattice of 2×2 spots through step 500.

To measure "distance" inside a coordinate-free system, we measured the **triangle view-angle ratio ρ = L₃ / L₁** (shortest side divided by longest side; 1 = equilateral, smaller = elongated / distant) across all 262,144 triangle loops.

[Figure 2: results/figures_en/fig_n64_viewangle_vs_condensation_en.png]

Figure 2. Particle condensation (top left: E_512), exponent contrast (top right), and triangle view-angle ratio ρ (bottom left). Synchronized with energy condensation into the 512 cells during steps 65–75, background triangles plunge below 10⁻³⁰⁰ while particle-internal triangles rebound in a V-shape back to equilateral geometry (ρ ~ 1).

[Figure 3: results/figures_en/fig_n64_viewangle_spectrum_en.png]

Figure 3. Bimodal splitting of the triangle view-angle spectrum at the particle representative: 192 background-connected triangles march steadily leftward to causal separation, while 64 particle-internal triangles condense at the equilateral peak.

Without coordinates or an inflaton field, **exponential background rarefaction (inflation-like separation)** and **V-shaped particle condensation** occurred simultaneously.

---

## 5. Comparison with Chandrasekaran (arXiv:2611.34924v1) and Honest Limitations

Chandrasekaran's preprint starts from 4D vacuum General Relativity on the null cone without a fixed background metric. While our work shares structural themes (no fixed background, internal reference frames, triple consistency, subspace separation), this rapid communication adds **fully reproducible deterministic code**, **spontaneous selection of 4D/8D algebras**, and **emergent particle-like remnants and bimodal view-angle splitting**.

To be completely transparent about current limits:
- Interpreting these numerical structures as "spatial expansion" and "particles" remains a **physical hypothesis**.
- Multi-particle scattering of spatially separated independent particles is an open problem.
- Quantitative recovery of semiclassical General Relativity or Standard Model constants lies beyond this rapid communication.

---

## 6. Behind the Scenes: Switching AI Models Unlocked a Half-Day Breakthrough

This project initially started with another AI (Claude Code), which got stuck in dead ends involving 2D matrix multiplication and a misinterpretation of exponentiation. Handing the algebraic redesign to **Google Gemini 4 Neon** unlocked the triangle-loop product, exact proofs, exhaustive 1,024-case classification, 20,000-step tracking, N = 64 view-angle measurements, and bilingual publication in roughly half a day. For transparency, the earlier failed logs remain public in the repository as well.

---

## 7. Try It Yourself in Seconds

```bash
git clone https://github.com/WurabeSeiji/ai-chat-logs-open.git
cd ai-chat-logs-open/次元の生成構造/三次元正方配列の積と自己無撞着_20261127

python3 tri_loop_product_v0.py
python3 tri_loop_n64_triangle_viewangle_run_v0.py
python3 plot_tri_loop_n64_viewangle_v0.py
```

Everything is deterministic with zero randomness—feel free to run it and verify every number yourself!

---

#QuantumGravity #TheoreticalPhysics #Physics #Mathematics #Cosmology #AI #Gemini #Python #Octonions #Inflation #BackgroundIndependence #ReproducibleScience
