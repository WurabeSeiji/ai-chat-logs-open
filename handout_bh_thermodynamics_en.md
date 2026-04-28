# BH Thermodynamics Programme: Overview of the Six-Paper Research Programme

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. / Osaka University, Faculty of Engineering Science (B.Eng.)  
**ORCID:** [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)

---

## Three Claims

**Claim 1: Asymptotic constant $c = 8/(3\pi)$ for unit-cube packing in a 4D ball**  
For a four-dimensional ball of radius $R = 2k+1$, the count $N(k)$ of integer-centred unit cubes contained in it gives a volume deficit $\Delta(R) = V_4(R) - N(k)$ whose leading constant is derived in closed form by inclusion–exclusion: $c = 8/(3\pi) \approx 0.84883$. (Paper 3)

**Claim 2: Exact agreement with 4-dimensional Tangherlini black hole thermodynamics**  
Black hole entropy is proportional to the surface area of the event horizon, which for a four-dimensional ball of diameter $R$ scales as $R^3$. Self-consistent determination of the central-projection metric yields $R \propto M^{1/2}$, Hawking temperature $T_H \propto M^{-1/2}$, lifetime $\tau \propto M^{3/2}$, and evaporation rate $dM/dt \propto -M^{-1/2}$. The ratio with the Bekenstein–Hawking entropy is $\Delta(R)/S_{BH} \to 32/(3\pi) \approx 3.40$. (Papers 4 and 5)

**Claim 3: Independent path without holographic assumptions**  
The programme does **not** assume the 't Hooft–Susskind holographic principle, the Maldacena AdS/CFT correspondence, or any specific quantum-gravity programme (LQG/CDT/string theory). It derives the above results from minimal geometric and combinatorial assumptions only — central projection and Planck-scale discreteness. Zero self-citation. (Papers 1–6)

---

## Programme Structure

**Paper 1: Overview and Hypotheses**  
From premises A (entropy is geometric), B (evaporation is horizon-local), and C (spacetime may be discrete at the Planck scale), two conditional propositions are extracted: Hypothesis 1 (Projection) and Hypothesis 2 (Discrete Drift).

**Paper 2: Mathematical Foundations of Central Projection**  
Gnomonic projection from $S^4(R) \subset \mathbb{R}^5$ to a four-dimensional tangent hyperplane. The induced metric $g_{\mu\nu} = (R^2/\ell^2)(\delta_{\mu\nu} - x_\mu x_\nu/\ell^2)$ satisfies $G_{\mu\nu} + (3/R^2)g_{\mu\nu} = 0$, becomes de Sitter under Lorentzian continuation, and Schwarzschild–de Sitter with a mass parameter.

**Paper 3: Mathematical Foundations of Discreteness**  
Packing condition $\sum (|c_i|+1/2)^2 \le R^2$. Exact computation of $N(k)$ for $k = 0..60$ (with $N(60) = 1{,}028{,}515{,}513$). Asymptotic expansion $\Delta(R) = (16\pi/3)R^3 - 6\pi R^2 + O(R)$. Connection to the Lagrange–Jacobi four-square representation count $r_4(N) = 8\sigma(N)$.

**Papers 4 and 5: Physical Connection and Dynamics**  
Constant ratio $32/(3\pi)$ between $\Delta(R)$ and the 4-dimensional Bekenstein–Hawking entropy. From a Friedmann-like relation $E_{\mathrm{total}} \propto R^2$, all Tangherlini scaling laws are recovered.

**Paper 6: Synthesis and Open Problem**  
The principal open problem — how the 4-dimensional results reduce to observed 3+1-dimensional black hole physics — is precisely formulated. Three candidate mechanisms: KK compactification, Randall–Sundrum/ADD brane localization, dS-CFT-style holographic projection.

---

## DOIs of the Six Papers

| Paper | DOI |
|-------|-----|
| 1 (Overview)  | [10.5281/zenodo.19837587](https://doi.org/10.5281/zenodo.19837587) |
| 2 (Projection) | [10.5281/zenodo.19837589](https://doi.org/10.5281/zenodo.19837589) |
| 3 (Packing) | [10.5281/zenodo.19837591](https://doi.org/10.5281/zenodo.19837591) |
| 4 (Entropy) | [10.5281/zenodo.19837593](https://doi.org/10.5281/zenodo.19837593) |
| 5 (Dynamics) | [10.5281/zenodo.19837595](https://doi.org/10.5281/zenodo.19837595) |
| 6 (Synthesis) | [10.5281/zenodo.19837597](https://doi.org/10.5281/zenodo.19837597) |

## Correspondence with 4-dimensional Tangherlini

| Quantity | Central-projection | 4-dimensional Tangherlini |
|----------|--------------------|-----------------|
| Mass–radius | $R \propto M^{1/2}$ | $r_h \propto M^{1/2}$ |
| Hawking temperature | $T_H \propto M^{-1/2}$ | identical |
| Lifetime | $\tau \propto M^{3/2}$ | identical |
| Entropy scaling | $\Delta(R) \propto R^3$ | $S_{BH} \propto R^3$ |
| Ratio | $\Delta/S_{BH} \to 32/(3\pi)$ | — |

---

## What This Research Does *Not* Claim

- That physical spacetime is 4-dimensional
- That higher-dimensional spacetime physically exists
- A first-principles derivation of the Bekenstein–Hawking entropy
- A replacement of holographic or quantum-gravity programmes

The "discrepancy" in the evaporation-rate exponent originally flagged in the programme outline is, in fact, a confusion between 3+1 and 4 conventions; in 4 dimensions the agreement is exact. The genuine open problem is the 4 → 3+1 reduction mechanism.

---

**Zenn article**: https://zenn.dev/noriaki_kihara/articles/bh-thermodynamics-projection  
**Note article (Japanese)**: https://note.com/kiharanoriaki/n/n4524006ef175  
**Note article (English)**: https://note.com/kiharanoriaki/n/n916eac778c1c  
**Source code (GitHub)**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/BH熱力学プログラム  
All papers are available on Zenodo as PDF, LaTeX, and Markdown (Japanese/English bilingual, CC BY 4.0).
