# Central & Spherical Projection: Geometry, the Multi-Axis Model, and a Glance at a Physical Stage (Observational Note)

**Noriaki Kihara** / WF System Co., Ltd. / Faculty of Engineering Science, Osaka University (graduate)
ORCID 0009-0004-6753-4020 / Contact: kihara.noriaki@gmail.com / CC BY 4.0 / 2026-06

> **Observational note.** No new geometric theorem is claimed. A classical map (radial projection) is organized, and the structures appearing on it are observed. No physical derivation is made.

## 1. Spherical projection $\sigma_R$ (radial projection)

$$\sigma_R:\mathbb{R}^{n+1}\setminus\{0\}\to S^n(R),\qquad \sigma_R(x)=\frac{R}{\|x\|}\,x$$

- Identical to the **radial projection** (the deformation retract of Hatcher, *Algebraic Topology*, Ch.0). A $C^\infty$ retraction, idempotent $\sigma_R\circ\sigma_R=\sigma_R$.
- Quotient $(\mathbb{R}^{n+1}\setminus\{0\})/\mathbb{R}_{>0}\cong S^n(R)$. $\ker D\sigma_R=$ radial $\mathrm{span}\{x\}$; $\operatorname{im}D\sigma_R=x^{\perp}=T_{\sigma_R(x)}S^n(R)$ (a submersion of rank $n$). Angle-preserving; direction invariant under the radius scale.

## 2. Central (gnomonic) projection $\Phi_R=\sigma_R|_{\Pi_R}$

Restriction to the tangent hyperplane $\Pi_R=\{x_{n+1}=R\}\cong\mathbb{R}^n$:

$$\Phi_R=\sigma_R\big|_{\Pi_R}:\ \Pi_R\ \xrightarrow{\ \sim\ }\ S^n_+(R)\quad(\text{open upper hemisphere; diffeomorphism})$$

- The pullback metric $g_{\mu\nu}=\dfrac{R^2}{\ell^2}\!\left(\delta_{\mu\nu}-\dfrac{x_\mu x_\nu}{\ell^2}\right)$ ($\ell=\sqrt{R^2+|x|^2}$) satisfies $G_{\mu\nu}+\Lambda g_{\mu\nu}=0$ ($\Lambda=\tfrac{(n-1)(n-2)}{2R^2}$), coinciding intrinsically with the **Beltrami coordinates of de Sitter space**. As $R\to\infty$ it degenerates to the flat (Minkowski) case.
- **The crux**: the contrast between $\sigma_R$ (non-injective, whole sphere) and $\Phi_R$ (injective, upper hemisphere only).

## 3. The multi-axis model

Choosing any of the $n+1$ axes of the background $\mathbb{R}^{n+1}$ as the projection center yields the same metric/curvature structure (**axis equivalence**); $n+1$ subjective spaces coexist, with transitions $T_{A\to B}=\Phi_B^{-1}\circ\Phi_A$. The roles of the "projection-center axis" (inaccessible from inside) and a "subjective coordinate axis" interchange with the observer.

## 4. A glance at a physical stage (observation only)

Taking the central projection of a vacuum (an exterior-free intrinsic system) as the homogeneous 4-sphere $S^4(R_{\mathrm U})$ (sectional curvature $1/R_{\mathrm U}^2$, the Euclidean version of $dS_4$), one observes that "phase-bearing local modes (particle-like states)" $P_a=(\boldsymbol{\nu}_a,R_a)$ with spread $W_a=2R_a$ can be placed on it. No Lorentzian metric, causal structure, or physical derivation is claimed. [Paper 14]

## References (Concept DOI)

- Definition of the spherical projection $\sigma_R$ (foundational map): [10.5281/zenodo.20462569](https://doi.org/10.5281/zenodo.20462569)
- A geometric formulation of 4-dimensional space by central projection: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)
- Composition of central projections (multi-axis; closed form of the composite curvature radius): [10.5281/zenodo.20060728](https://doi.org/10.5281/zenodo.20060728)
- [Paper 14] Central projection of a vacuum universe and particle-like states with a spread phase: [10.5281/zenodo.20543044](https://doi.org/10.5281/zenodo.20543044)

GitHub: github.com/WurabeSeiji/ai-chat-logs-open / Full texts on Zenodo (JA/EN, CC BY 4.0)
