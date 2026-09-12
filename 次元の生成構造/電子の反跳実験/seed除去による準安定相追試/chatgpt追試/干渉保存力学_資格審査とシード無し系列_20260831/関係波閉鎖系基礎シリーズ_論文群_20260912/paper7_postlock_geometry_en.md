# Post-Lock Geometry — The Initial Plane and the Terminal Plane Span 4 Dimensions (Two Nonzero Principal Angles); Plane Switching Is Not an Aggregation Artifact of a Long-Window SVD

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 7 of 10)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22729108
**Concept DOI:** 10.5281/zenodo.22729107
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22729108

---

## Abstract

We fix geometrically what inflation (Paper 6) and SSB lock (Paper 5) do in state space. The map is a real orthogonal rotation $z(t)=Q_t z_0$, and the state remains on a 2-dimensional plane (rank 2) at every instant (the singular values of $[\Re z,\Im z]$ take only two values at every step, and $z^{\mathsf T}z=0$ conservation gives an orthogonal isometric frame). On the floor the floor plane $P_0=\mathrm{span}\{\Re z_0,\Im z_0\}$ is $K$-invariant and the state stays in $P_0$. When SSB leaves the floor, the lock plane $P_{\rm term}$ tilts away from $P_0$, and the principal angles between the two planes are, independently of $N$, **both nonzero** ($N=6{:}27.1^\circ,30.8^\circ$／$N=40{:}18.1^\circ,20.4^\circ$), so the dimension spanned by the two planes is $4$. $P_0$ and $P_{\rm term}$ span 4 dimensions with two nonzero principal angles, and from the same floor the terminal orientation is seed-dependent and not uniquely determined. We do not claim a specific $SO(4)$ rotation structure for the cumulative map (in the real data $M_t$ is confirmed not to preserve a 4-dimensional subspace). The tilt of the lock plane away from the floor is, over all $N=6$–$40$, $14.2^\circ$–$29.0^\circ$ (median $\sim19.2^\circ$). During the transition a second direction rises (transient $s_3/s_1\approx0.10$–$0.20$ = effective rank 4), and after lock the system re-locks to a single plane (tail $s_3/s_1\approx10^{-5}$–$10^{-13}$ = rank 2). By a short-window SVD test we separate this "plane switching" from an appearance of long-window aggregation, showing it is a **real plane tilt**. Classification: derived geometry + measurement (artifact excluded).

---

## 1. Problem

Reading the endpoint merely as a "circle" fails to capture the essential observation of this map — **switching from one plane to another plane**. We quantify this for all $N\ge6$ and separate the appearance (aggregation of a long-window SVD) from the real rotation.

## 2. Instantaneous Structure — Rank 2 (a Single Plane)

The map is a real orthogonal rotation $z(t)=Q_t z_0$ ($Q_t$ real orthogonal, acting identically on $\Re z,\Im z$), and the state stays on a 2-dimensional plane at every instant (the singular values of $[\Re z,\Im z]\in\mathbb R^{2\times M}$ take only two values at every step, and $z^{\mathsf T}z=0$ conservation gives an orthogonal isometric frame). **Formulation as Grassmann motion**: for $z=a+ib$, from $z^{\mathsf T}z=0$ and $\|z\|=1$ we have $\|a\|^2=\|b\|^2=\tfrac12$ and $a^{\mathsf T}b=0$, so $u=\sqrt2\,a,\ v=\sqrt2\,b$ are always an orthonormal real 2-frame $(u,v)\in V_2(\mathbb R^M)$. The global phase $z\mapsto e^{i\theta}z$ merely rotates $(u,v)$ within the same plane, so the object observed — modulo global phase — is the **2-dimensional subspace $P(t)=\mathrm{span}\{u,v\}$ itself** (motion on the Grassmann manifold $\mathrm{Gr}(2,M)$). The projection $\Pi_t=u_tu_t^{\mathsf T}+v_tv_t^{\mathsf T}$ moves exactly, by each step's real orthogonal $Q_t=e^{\Delta\tau K_t}$, as

$$\Pi_{t+1}=Q_t\,\Pi_t\,Q_t^{\mathsf T}$$

On the floor $Q_0P_0=P_0$, so $\Pi_{t+1}=\Pi_t$ (relative equilibrium, fixed); in inflation $\Pi_t\neq\Pi_0$ and the plane moves; at saturation it re-locks to $P_{\rm term}$. Sections 3–5 measure this plane motion $\Pi_0\to\Pi_t\to\Pi_{\rm term}$. On the floor $P_0=\mathrm{span}\{\Re z_0,\Im z_0\}$ is $K$-invariant ($Ka=\mu b,\ Kb=-\mu a$) and the state stays in $P_0$ (principal angles zero). The state-space SVD ratio $s_3/s_1$ of the third/first singular values over the last 200 steps after lock is $10^{-5}$–$10^{-13}$ = a single plane (rank 2/$S^1$).

## 3. The Dimension Spanned by the Initial Plane and the Terminal Plane — 4D Span

When SSB leaves the floor, the lock plane $P_{\rm term}$ tilts away from $P_0$. The principal angles between the initial plane $P_0$ and the terminal plane $P_{\rm term}$ are, independently of $N$, **both nonzero** ($N=6{:}27.1^\circ,30.8^\circ$／$N=40{:}18.1^\circ,20.4^\circ$), so the dimension spanned by the two planes is $4$. The tilt of the lock plane away from the floor ($\arcsin\sqrt{H_\perp/H_{\rm end}}$) is, over all $N=6$–$40$, as follows (excerpted from `plane_switch_allN.csv`):

| $N$ | tilt [°] | transient $s_3/s_1$ | tail $s_3/s_1$ |
|---|---|---|---|
| 6 | 29.0 | 0.172 | $3.8\times10^{-6}$ |
| 8 | 20.3 | 0.135 | $2.8\times10^{-4}$ |
| 16 | 24.5 | 0.172 | $1.1\times10^{-5}$ |
| 22 | 15.7 | 0.130 | $1.1\times10^{-13}$ |
| 40 | 19.3 | 0.132 | $4.5\times10^{-7}$ |
| 64 | 21.7 | 0.021 | $2.5\times10^{-4}$ |

Over all $N$ the tilt is $14.2^\circ$–$29.0^\circ$ (median $\sim19.2^\circ$), neither $0^\circ$ nor $90^\circ$ = isotropy is broken.

**The tilt, the two principal angles, and $H_\perp/H$ are not independent measured quantities but the same Grassmann geometry.** For the initial and terminal frames $U_0=(u_0,v_0),\ U_1=(u_1,v_1)$ and $P_0^\perp=I-U_0U_0^{\mathsf T}$, from the definition of principal angles (the singular values of $U_0^{\mathsf T}U_1$ equal $\cos\theta_i$) we have $\|P_0^\perp U_1\|_F^2=\sin^2\theta_1+\sin^2\theta_2$. Since $z=\tfrac1{\sqrt2}(u+iv)$,

$$\frac{H_\perp}{H}=\frac12\left(\sin^2\theta_1+\sin^2\theta_2\right)=\frac14\big\|\Pi_{\rm term}-\Pi_0\big\|_F^2,\qquad \alpha=\arcsin\sqrt{\frac{\sin^2\theta_1+\sin^2\theta_2}{2}}$$

(from $\Pi=UU^{\mathsf T}$ we get $\|\Pi_1-\Pi_0\|_F^2=2(\sin^2\theta_1+\sin^2\theta_2)$). This is not an additional hypothesis but an **exact geometric identity** from the definitions, and $N=6$'s $(27.1^\circ,30.8^\circ)\to\alpha=29.0^\circ$ and $N=40$'s $(18.1^\circ,20.4^\circ)\to\alpha=19.3^\circ$ agree with the table. That is, the tilt $\alpha$ away from the floor is the **Grassmann distance** $\tfrac14\|\Pi_{\rm term}-\Pi_0\|_F^2$ from the initial plane itself.

![Figure 1: Plane switching (all N≥6). Left = tilt of the lock plane away from the floor vs N; right = transient s3/s1 (rank 4/S³) → tail s3/s1 (rank 2/S¹).](../位相ロック全N確認_相対平衡_20260912/面の乗り換え全N_20260912/fig_plane_switch_allN.png)

**Figure 1** $N$-dependence of plane switching.

## 4. The Initial Plane and the Terminal Plane Span 4 Dimensions

$P_0$ and $P_{\rm term}$ meet at two nonzero principal angles ($N=6{:}27.1^\circ,30.8^\circ$／$N=40{:}18.1^\circ,20.4^\circ$) and $\dim(P_0+P_{\rm term})=4$ — the terminal plane tilts away from the initial plane in two independent directions (these "4 directions" are meant only in the sense of the dimension of the span). **From the same floor this terminal orientation is seed-dependent and not uniquely determined** (Paper 5). We do not claim a specific $SO(4)$ rotation structure for the cumulative map $M_t=\prod_n\exp(\Delta\tau K(\varphi_n))$ — in the real data $M_t$ is confirmed not to preserve this 4-dimensional subspace (the invariance $M_tS=S$ is broken), and $M_t$ is an orthogonal rotation of the full $\mathbb R^M$ that tilts $P_0$ into the 4-dimensional span.

During the transition a second rotation plane rises (transient $s_3/s_1\approx0.10$–$0.20$ = rank 4), but after lock the system re-locks to a single plane (tail $s_3/s_1\approx0$ = rank 2). At $N=64$ the transient $s_3/s_1=0.021$ is small, so the second-plane transient weakens for large $N$. This transitional rank 4 is temporary; after lock it returns to rank 2, and the cumulative map also does not preserve a 4-dimensional subspace (§4) — **no persistent $S^3$ structure is found** (inspected, short-window SVD + discriminating computation).

In the time-varying, non-normal tangent cocycle confirmed in Paper 6, the maximally amplified direction rotates with time while it is amplified. The additional direction of the transient observed in this paper can be read as the geometric manifestation of that tangent dynamics tilting the instantaneous rank-2 plane $P(t)$ with time. Therefore the effective rank 4 of the transient is **not a 4-dimensionalization of the instantaneous state but the local span swept in a short time by the moving 2-plane $P(t)$**, and after saturation it re-locks at rank 2 to a single terminal plane. That is, Paper 6 (in tangent space the amplified direction rotates) and Paper 7 (in state space the 2-plane $P(t)$ tilts and moves to another 2-plane) are the tangent-side and state-side representations of the same phenomenon. Using $H_\perp/H=\tfrac14\|\Pi_t-\Pi_0\|_F^2$ (the §3 identity), Paper 6's initial ramp $H_\perp/H\propto e^{2\gamma t}$ is geometrically $\|\Pi_t-\Pi_0\|_F\propto e^{\gamma t}$ — **inflation is an exponential departure from the initial plane on the Grassmann space $\mathrm{Gr}(2,M)$**.

## 5. Exclusion of the Artifact

A long-window SVD can, if a single plane rotates with time, be aggregated into an apparently higher rank. With the short-window SVD test (`artifact_test_shortwindow_svd`) we separated the fact that, instantaneously staying at rank 2, the initial plane $P_0$ tilts in two directions ($P_0$ and $P_{\rm term}$ span 4 dimensions) — the "plane switching" is not an appearance of long-window aggregation but a **real plane tilt**.

![Figure 2: Short-window SVD test. While preserving instantaneous rank 2, the initial plane tilts in two directions (plane switching is real).](../位相ロック全N確認_相対平衡_20260912/面の乗り換え全N_20260912/fig_artifact_shortwindow_svd.png)

**Figure 2** Artifact exclusion (short-window SVD).

## 6. A Caveat on Readout

When there is initially only a single complex plane, its normal is a degenerate orthogonal complement and is unobservable (conceptual). Only once inflation raises a specific direction in the normal space does the relative phase (tilt) between the initial plane and the new plane become readable. **What can be read out is always a relational quantity (the relative phase between planes)**; the absolute normal direction cannot be recovered from the state alone ($Q_t$ hides the orientation, and there is no conserved quantity pointing to it).

## 7. Claim Classification and Open Items

- Classification: derived geometry (real orthogonal rotation rank 2) + measurement (principal angles, tilt, 4D span, artifact exclusion) + discriminating computation (confirmed that the cumulative map does not preserve a 4-dimensional subspace and is not $SO(4)$). Verdict: the 4D span, the two nonzero principal angles, and the double rotation of the relative configuration are retained; "the dynamics being an $SO(4)$ double rotation on 4 dimensions" is rejected.
- Reference (a different regime): the crystalline terminal states at $N=3,4$ are special solutions in which the plane refreezes at the magic angle $\arccos\sqrt{2/3}=35.264^\circ$ (the general $N\ge6$ lock is $14$–$29^\circ$). The crystal/glass bifurcation is a separate topic of this series.
- Open: the physical identification of the 4D span ($xyzt$ or $xyztRQ$, etc.) and the type of gauge (scale) are outside the scope of this paper (Paper 0 §Open).

## Related Work (Structural Agreement, Not the Source of Derivation)

- The double rotation of $SO(4)$ (two invariant planes, two independent rotation angles): this is a description of the **relative configuration** of $P_0\to P_{\rm term}$, not a property of the cumulative map itself (in the §4 discriminating computation $M_tS\neq S$).
- Relative equilibrium (Marsden & Ratiu): the rigid-body-rotation relative equilibrium of the floor and the lock.

## References

1. Noriaki Kihara, "The mechanism of inflationary rapid expansion in self-consistent relational-wave closed systems," Concept DOI 10.5281/zenodo.22112008.
2. J. E. Marsden and T. S. Ratiu, *Introduction to Mechanics and Symmetry*, Springer (1999).

## Reproduction

The analysis scripts, figures, and CSV are in `位相ロック全N確認_相対平衡_20260912/面の乗り換え全N_20260912/` (`analyze_plane_switch_allN_20260912.py`, `artifact_test_shortwindow_svd_20260912.py`, `plane_switch_allN.csv`, `fig_plane_switch_allN.png`, `fig_artifact_shortwindow_svd.png`). For the locked data, small $N$ and $N=23$–$29$ are 500 steps, and the unsaturated $N$ (22, 30–40) and $N=64$ are 2000 steps. The principal angles of the initial and terminal planes are the step0/step2000 complex planes in `N22_2000step検証走行_/` and `N40_2000step検証走行_/`. The early rank-3 audit is in `independent_rank3_geometry_audit_20260831/`. The dynamics map is identical to the canonical `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`), readout only.
