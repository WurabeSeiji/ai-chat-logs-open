# Analytic Exact Solution of the 90-Degree High-Symmetry Floor — The $\sigma_{\max}$ Eigenmode of the Integer Antisymmetric Generator Fixed by Z4 Phase Labels, and Its Closed-Form Spectrum

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 3 of 10)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22728996
**Concept DOI:** 10.5281/zenodo.22728995
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22728996

---

## Abstract

The 90-degree high-symmetry floor, established in Paper 2 as satisfying the ignition qualification, is here constructed **not** by iterative approximation (make_parent, `iters=1200`, `tol=1e-12`) but **analytically and exactly as an integer eigenvalue problem**, and its closed form is derived for **arbitrary $N$**. Because the floor phases, up to a global phase, are quantized to $\{0,90,180,270\}^\circ$, assigning each edge a Z4 label $k_e\in\{0,1,2,3\}$ makes the generator $K_{ef}=A_{ef}\sin(90^\circ(k_f-k_e))$ an **integer matrix (entries $0,\pm1$)**. For the deterministic seed $k_e=(i+j)\bmod4$, what determines the generator is not mod 4 but only the **parity** $p_e=(i+j)\bmod2$, and the real symmetric matrix $W_{ef}=A_{ef}\sin^2\psi_{ef}$ of Paper 1 becomes a bipartite graph $W_{ef}=A_{ef}\mathbf 1_{p_e\neq p_f}$ ($L(K_N)$ split by the parity of the sum of vertex indices). Therefore the floor can be constructed **in one shot** by $k_e=(i+j)\bmod4\to W\to$ Perron vector $r>0\to z=Dr$ ($D=\operatorname{diag}(e^{i90^\circ k_e})$, Paper 1 $D^{-1}KD=iW$), and **discrete-label iteration is essentially unnecessary**. The floor is the eigenvector of the **negative eigenvalue of largest absolute value** $\mu=-\sigma_{\max}$ of $iK$, and its phase advances as $d\varphi/d\tau=-\mu=\sigma_{\max}>0$, consistent with the phase equation. From the reduced Perron problem of the bipartite structure ($2\times2$ for even $N$, $3\times3$ for odd $N$), $\sigma^2$ (even $N(N-2)$ / odd $N^2-2N-1$), per-edge squared amplitude (rational within each parity class), and the full spectrum ($BB^{\mathsf T}$) all come out in closed form for arbitrary $N$ (**verified symbolically in $N$-general with sympy, and numerically matched over $N=3$–$40$**). Moreover, for general $N$ the number of nonzero spectral values is $O(N)$, and for even $N\ge6$ / odd $N\ge7$ it is exactly $2(N-1)$ $=\operatorname{rank}K=2(N-1)$ with kernel dimension $(N-1)(N-4)/2$ ($N=3,4,5$ are degenerate exceptions): of the $M=N(N-1)/2=O(N^2)$-dimensional relational-wave space, only $O(N)$ dimensions are active. Over all $N=3,\ldots,40$ (38 floors), the one_step relative-equilibrium residual is $\sim10^{-15}$–$1.2\times10^{-14}$ (better than make_parent's $\sim10^{-13}$), the phase lies exactly on the axes, and $\|z\|^2=1$. **Zero-closure $\sum z^2=0$ is not an extra condition imposed on the floor but a general theorem that follows automatically from the mode being a bipartite eigenmode ($\|x\|^2=\|y\|^2$)** (holding not only for the Perron mode but for every nonzero eigenmode), and the connectivity of the bipartite $W$ is also proven for arbitrary $N\ge3$ (recovering the applicability condition of the Perron argument in Paper 1). The floor can be written down completely algebraically as "$\sigma=\sqrt{\text{integer}}$, per-edge amplitude $=\sqrt{\text{rational}}$ (fixed by the parity class), phase $\in\{0,90,180,270\}$." Classification: analytic exact solution (parity bipartite structure + reduced Perron problem) + symbolic (arbitrary $N$) / numerical ($N=3$–$40$) verification.

---

## 1. The Problem — Making the Approximate Floor Exact

The make_parent floor used in Paper 2 is a solution that iteratively approximates the self-consistent eigenmode, with residual $\sim10^{-13}$. Merely copying an approximation stopped by an iteration cap, tolerance cut, and damping is meaningless. The 90-degree floor has the structural property that its phases are quantized to a discrete set (Z4), and using this it can be solved **algebraically and exactly as an integer eigenvalue problem**. This paper constructs that exact solution for all $N$ (including odd $N$) and fixes the closed form.

## 2. Principle — From Z4 Labels to the Integer Antisymmetric Generator

Quantize the floor phases, up to a global phase, into the four quadrants and assign each edge $k_e\in\{0,1,2,3\}$. Since the phase difference is $90^\circ(k_f-k_e)$,

$$K_{ef}=A_{ef}\sin\!\big(90^\circ(k_f-k_e)\big),\qquad (k_f-k_e)\bmod4:\ 0\to0,\ 1\to+1,\ 2\to0,\ 3\to-1,$$

that is, $K$ is an **integer antisymmetric matrix** with entries $0,\pm1$. The floor $z$ is its $\sigma_{\max}$ eigenmode, i.e. the eigenvector of the **negative eigenvalue of largest absolute value** $\mu=-\sigma_{\max}$ of $iK$ (of the form $z=\sigma_{\max}u-iKu$, where $u$ is the real eigenvector of $K^2u=-\sigma_{\max}^2u$). $iK$ is Hermitian with eigenvalues in conjugate pairs $\pm\sigma_j$, and $\sigma_{\max}$ is the largest singular value (spectral radius). The floor takes the branch $\mu=-\sigma_{\max}$ where the phase advances, and $d\varphi/d\tau=-\mu=\sigma_{\max}>0$ is consistent with the phase equation $r\dot\varphi=\sum r_f\sin^2\Delta\varphi\ge0$ (confirmed with real data and code: for N=6 the make_parent floor has $\mu=-4.868$, the integer-K floor $\mu=-4.899$, both the smallest eigenvalue). What the square roots of the eigenvalues of $K^{\mathsf T}K$ fix are the **singular values $\sigma=\sqrt{\text{integer}}$**, whereas each per-edge amplitude $r_e$ is a component of the eigenvector belonging to it (= the Perron vector of $W$), and in this construction it takes the closed form $r_e=\sqrt{\text{rational}}$ (§5).

Note that the complex self-consistent problem $Kz=i\omega z$ is exactly equivalent to the eigenvalue problem $Wr=\omega r$ of the real symmetric non-negative integer matrix $W_{ef}=A_{ef}\sin^2\psi_{ef}$ ($D^{-1}KD=iW$, $D=\operatorname{diag}(e^{i\varphi_e})$), and the existence and uniqueness of a positive-amplitude solution by Perron–Frobenius when $W$ is irreducible is **given in general by Paper 1**. This paper, with that $W$ = the concrete form of the integer $K$, solves $\sigma=\rho(W)$ and each per-edge amplitude in closed form (§5) (the Perron eigenvector of $W$ and the $\sigma_{\max}$ mode of the integer $K$ are the same floor; confirmed $r>0$ and agreement of the principal eigenvector for all $N$).

**Bipartite structure from parity.** What determines $W$ is not mod 4 but only **mod 2** (parity) $p_e=(i+j)\bmod2$. Since $\sin^2[90^\circ(k_f-k_e)]$ is $1$ when $k_f-k_e$ is odd ($p_e\neq p_f$) and $0$ when even ($p_e=p_f$),

$$W_{ef}=A_{ef}\,\mathbf 1_{p_e\neq p_f}.$$

That is, $W$ is the **adjacency matrix of the bipartite graph** obtained by splitting the vertices of $L(K_N)$ (= the edges $\{i,j\}$ of $K_N$) by the parity $p_e$ of the sum of vertex indices. This $W$ is **connected (irreducible)** for arbitrary $N\ge3$ (proof in §8: for the even-vertex set $U$ ($|U|=\lceil N/2\rceil\ge2$) and odd-vertex set $V$, cross-parity edges are bridged by internal edges, and same-parity edges are connected via cross-parity edges). Hence the Perron vector $r>0$ is unique (Paper 1), and $z=Dr$ ($D=\operatorname{diag}(e^{i90^\circ k_e})$) exactly satisfies $Kz=i\sigma_{\max}z$. Therefore the floor can be constructed in one shot from the seed labels by $k\to W\to r\to z=Dr$, and the discrete-label iteration (§3) is unnecessary for constructing the floor (confirmed for all $N=3$–$40$ that the $W$ of the iterative floor agrees with the seed-parity $W$ to machine precision, $\|\text{diff}\|\le3\times10^{-27}$, `parity二部構造_解析的証明_20260912/`). From this bipartite structure the closed forms of §5 and §6 come out for arbitrary $N$.

**Square-closure $\sum_e z_e^2=0$ is a consequence of the bipartite structure (not an extra condition).** For the bipartite $W=\begin{pmatrix}0&B\\B^{\mathsf T}&0\end{pmatrix}$, a nonzero eigenmode $r=(x,y)$ satisfies $By=\sigma x,\ B^{\mathsf T}x=\sigma y$, so $\sigma\|x\|^2=x^{\mathsf T}By=(B^{\mathsf T}x)^{\mathsf T}y=\sigma\|y\|^2$, and since $\sigma\neq0$, $\|x\|^2=\|y\|^2$ (each $1/2$ under normalization). On the other hand, from $z_e=e^{i90^\circ k_e}r_e$ we have $z_e^2=(-1)^{k_e}r_e^2=(-1)^{p_e}r_e^2$ (real), so

$$\sum_e z_e^2=\sum_{p_e=0}r_e^2-\sum_{p_e=1}r_e^2=\|x\|^2-\|y\|^2=0.$$

That is, zero-closure $\sum z^2=0$ is not an extra condition imposed on the floor but **arises automatically from being a bipartite eigenmode** (holding not only for the Perron floor but for any nonzero eigenmode of the bipartite $W$; the real part and the cross-term imaginary part vanish separately as a consequence of $z_e^2$ being real). Confirmed $|\sum z^2|\le10^{-15}$ and $|\|x\|^2-\|y\|^2|\le10^{-15}$ for every nonzero eigenmode over all $N=3$–$40$.

## 3. Construction of the Floor — One Shot from the Seed Labels via Perron (Iteration for Confirmation Only)

By §2 the floor can be constructed **directly without iteration** from the seed labels $k_e=(i+j)\bmod4$ (edge $e=\{i,j\}$, $0\le i<j\le N-1$):

1. Build the bipartite matrix $W_{ef}=A_{ef}\mathbf 1_{p_e\neq p_f}$ using the parity $p_e=(i+j)\bmod2$.
2. Take the Perron vector $r>0$ of $W$ and $\sigma_{\max}=\rho(W)$ ($W$ connected, Paper 1).
3. With $D=\operatorname{diag}(e^{i90^\circ k_e})$, the floor is $z=Dr$. By Paper 1 $D^{-1}KD=iW$, so $Kz=i\sigma_{\max}z$ (relative equilibrium).

We confirmed for all $N=3$–$40$ that this direct construction $z=Dr$ satisfies a one_step relative-equilibrium residual $\le1.2\times10^{-14}$. The conventional generator `make_90deg_floor_analytic_v1.py` obtains the floor from the seed via discrete-label iteration (re-quantize the eigenvector of the eigenvalue of largest absolute value $\mu=-\sigma_{\max}$ of $iK$ into the four quadrants → best-match to the previous labels over the four global-phase choices → stop at a fixed point / gauge-equivalent cycle, adopting the minimum-residual floor), but the **$W$ of the iterative floor agrees with the seed-parity $W$ to machine precision** ($\|\text{diff}\|\le3\times10^{-27}$), and the floor is identical to the direct seed construction up to a **sign gauge** ($z_e\to-z_e$, $k_e\to k_e+2$, $K\to SKS$ with $Kz=i\sigma_{\max}z$ and $W$ invariant) (after gauge alignment $\|\text{diff}\|\le8\times10^{-15}$). Therefore the iteration is unnecessary for constructing the floor and amounts to confirming that the seed is a fixed point (up to the sign gauge). adjacency/one_step are extracted from the canonical source `run_N3_N40_stage123_v1.py` (SHA-checked) (no physics changed). For each floor's verification (`build_analytic_floor`), the following are recorded: $\|z\|^2$, $|\sum z|$, $\Re(a^{\mathsf T}a-b^{\mathsf T}b)$, $a^{\mathsf T}b$, $|\sum z^2|$, the one_step relative-equilibrium residual, the on-axis test, the number of quadrants, the number of amplitude classes, and $\sigma$.

## 4. Results (All $N=3,\ldots,40$ · 38 Floors)

For all floors the phase is exactly on-axis ($0/90/180/270$), closure $\sum z^2=0$ (both the real part and the cross-term imaginary part machine-zero; a consequence of the bipartite structure of §2), $\|z\|^2=1$, and the relative-equilibrium residual is $\sim10^{-15}$–$1.2\times10^{-14}$.

- **$N=3$**: 3 quadrants (T-shape), $|\sum z|=0.707$, $\sigma=\sqrt2=1.414$.
- **$4\mid N$ (4,8,…,40)**: centroid $|\sum z|\sim10^{-15}$ (exactly zero), 2-valued amplitude.
- **$N\equiv2\ (\mathrm{mod}\,4)$ (6,10,…,38)**: centroid $0.236\to0.037$ (decreasing as $N$ increases), 2-valued amplitude.
- **Odd $N$ (5,7,…,39)**: centroid $0.154\to0.019$ (decreasing), 3-valued amplitude. **Odd $N$ are all floors too** (on-axis, zero closure).

An exactly-zero centroid occurs only for $4\mid N$ (where complete $\pm$ symmetry of the quadrants is possible); otherwise it is a residual centroid determined by parity imbalance and the Perron solution (a determined quantity, not a minimization; closed form in §8; an exact zero cannot be attained). Precision comparison (same $N$):

| $N$ | make_parent residual | integer-K analytic residual |
|---|---|---|
| 6 | $2.3\times10^{-13}$ | $5.4\times10^{-16}$ |
| 7 | $3.5\times10^{-14}$ | $1.7\times10^{-15}$ |
| 8 | $8.2\times10^{-14}$ | $5.7\times10^{-16}$ |
| 9 | $8.2\times10^{-14}$ | $1.0\times10^{-15}$ |

All 38 floors are collected in a complex-plane grid figure (Figure 1).

![Figure 1: 90-degree theoretical floors (integer-K σ_max analytic solution, all N=3..40). Phase 0/90/180/270 · exact relative equilibrium (residual ~1e-15) · centroid exactly zero for 4|N.](../90度理論床_整数K解析解_全N_20260909/fig_90deg_floor_complex_planes_analytic.png)

**Figure 1** 90-degree theoretical floors (integer-K $\sigma_{\max}$ analytic solution, $N=3$–$40$).

## 5. Closed Form (From the Parity-Reduced Perron Problem, Arbitrary $N$)

From the bipartite structure of §2, the Perron vector $r$ takes a constant value within each parity class. Counting the number of adjacencies for each class yields a reduced Perron problem, and $\sigma$ and each per-edge amplitude are solved in closed form for **arbitrary $N$** (sympy symbolic verification, `parity二部構造_解析的証明_20260912/results/実行ログ`). Since the integer $K$ is antisymmetric, $K^2$ is integer symmetric (eigenvalues $-\sigma^2$), so $\sigma=\sqrt{\text{integer}}$ and each per-edge squared amplitude is rational.

**Even $N=2m$ (2 classes).** Cross-parity edges ($p=1$, $i+j$ odd, count $N^2/4$) and same-parity internal edges ($p=0$, $i+j$ even, count $N(N-2)/4$). Since each cross-parity edge has $N-2$ same-parity edges adjacent and each same-parity edge has $N$ cross-parity edges adjacent, the class representative values $x$ (cross-parity) and $y$ (same-parity) satisfy

$$\sigma x=(N-2)y,\qquad \sigma y=Nx\ \Longrightarrow\ \sigma^2=N(N-2),$$

and normalizing by $\|z\|^2=1$ gives $x^2=2/N^2$, $y^2=2/(N(N-2))$.

**Odd $N$ (3 classes).** $a=\tfrac{N+1}2$ even vertices, $b=\tfrac{N-1}2$ odd vertices. Cross-parity ($ab=\tfrac{N^2-1}4$), even-side internal ($\binom a2=\tfrac{N^2-1}8$), odd-side internal ($\binom b2=\tfrac{(N-1)(N-3)}8$). The representative values $x$ (cross-parity), $y$ (even-side internal), $z$ (odd-side internal) satisfy

$$\sigma x=\tfrac{N-1}2y+\tfrac{N-3}2z,\quad \sigma y=(N-1)x,\quad \sigma z=(N+1)x\ \Longrightarrow\ \sigma^2=N^2-2N-1,$$

and normalizing gives $x^2=\tfrac{2}{(N-1)(N+1)}$, $y^2=\tfrac{2(N-1)}{(N+1)(N^2-2N-1)}$, $z^2=\tfrac{2(N+1)}{(N-1)(N^2-2N-1)}$.

**$\sigma^2$**: even $N(N-2)$ / odd $N^2-2N-1$. Sequence $N=3$–$20$: $2,8,14,24,34,48,62,80,98,120,142,168,194,224,254,288,322,360$.

**Squared amplitude (normalized $\|z\|^2=1$, value × count) — the class assignment is fixed by parity.**

Even $N$ (2 classes, each class total power $1/2$):

| Squared amplitude | Class (edge) | Count |
|---|---|---|
| $2/N^2$ | cross-parity ($i+j$ odd) | $(N/2)^2$ |
| $2/(N(N-2))$ | same-parity internal ($i+j$ even) | $N(N-2)/4$ |

Odd $N$ (3 classes; $N=3$ has 0 odd-side internal edges = T-shape):

| Squared amplitude | Class (edge) | Count |
|---|---|---|
| $2/((N-1)(N+1))$ | cross-parity ($i+j$ odd) | $(N^2-1)/4$ |
| $2(N-1)/((N+1)(N^2-2N-1))$ | even-side internal (both $i,j$ even) | $(N^2-1)/8$ |
| $2(N+1)/((N-1)(N^2-2N-1))$ | odd-side internal (both $i,j$ odd) | $(N-1)(N-3)/8$ |

In each case the total of the counts is $M=N(N-1)/2$. The closed form is verified symbolically with sympy (arbitrary $N$) and matches the numerical per-edge squared amplitudes over $N=3$–$40$ (error $\le1.2\times10^{-16}$). The centroid is exactly zero for even $N$ because the 2 classes, each with $1/2$ power, cancel completely as $\pm$. For odd $N$ the $\pm$ cancellation is incomplete across the 3 classes, leaving a small residual centroid (closed form in §8).

## 6. Eigen-Oscillation Spectrum and Active Dimension ($O(N)$)

The nonzero eigenvalues of the bipartite structure $W=\begin{pmatrix}0&B\\B^{\mathsf T}&0\end{pmatrix}$ are $\pm\sqrt{\operatorname{eig}(BB^{\mathsf T})}$. Indexing cross-parity edges by an $a\times b$ lattice of (even vertex)×(odd vertex),

$$BB^{\mathsf T}=(a-2)I_a\otimes J_b+(b-2)J_a\otimes I_b+2J_a\otimes J_b,$$

with eigenvalues $4ab-2(a+b)$ ($\times1$), $b(a-2)$ ($\times(a-1)$), $a(b-2)$ ($\times(b-1)$), and the rest $0$. Substituting $a,b$ (even $N$: $a=b=N/2$; odd $N$: $a=\tfrac{N+1}2,b=\tfrac{N-1}2$):

- **Even $N\ge6$**: $\lambda_{\max}^2=N(N-2)$ (multiplicity 1), $\lambda_1^2=N(N-4)/4$ (multiplicity $N-2$).
- **Odd $N\ge7$**: $\lambda_{\max}^2=N^2-2N-1$ (multiplicity 1), $\lambda_1^2=(N-5)(N+1)/4$ (multiplicity $(N-3)/2$), $\lambda_2^2=(N-1)(N-3)/4$ (multiplicity $(N-1)/2$).

Numerically matched for all $N$ (even $N\ge6$ / odd $N\ge7$). $N=3,4,5$ are degenerate and some eigenvalue families vanish ($N=4$: $\lambda_1^2=0$; $N=5$: $(N-5)(N+1)/4=0$).

**Active dimension $O(N)$.** The number of nonzero eigenvalues is $a+b-1=N-1$ (with $\pm$ pairs, $2(N-1)$). Hence

$$\operatorname{rank}K=\operatorname{rank}W=2(N-1),\qquad \dim\ker=M-2(N-1)=\frac{(N-1)(N-4)}2\quad(\text{even }N\ge6,\ \text{odd }N\ge7).$$

Of the $M=N(N-1)/2=O(N^2)$-dimensional relational-wave space, only $O(N)$ dimensions are active under the 90-degree floor's generator; the rest fall into $\ker K$. Small $N$ are exceptions due to degeneracy: $N=3$ (rank 2, kernel 1), $N=4$ (rank 2, kernel 4), $N=5$ (rank 6, kernel 4).

A simple integer-harmonic system (integer multiples of a fundamental) occurs only for $N=3,4$. For $N\ge5$ the ratios are in general square-root ratios. Comparing the raw eigenvalues $\lambda$ across $N$, there are 65 pairs of machine-precision integer ratios of $2$–$8$, and with the phase step included, $q=\lambda/N$ detects 11 pairs (`固有振動数解析_90度系列_N3_N40_20260910/`). This is consistent with the floor being a single eigenmode that co-winds a single carrier (Paper 4: even after factorization, the frequency ratio on the floor is 1:1).

## 7. Relation to the Non-90-Degree Theoretical Floor

| | Phase | Generator | Solution method | Symmetry |
|---|---|---|---|---|
| Cyclic-Fourier floor (Paper 2's (C)) | $2\pi(i+j)/N$ (N-valued) | continuous | distance-class-reduced eigenvector | cyclic symmetry (amplitude constant within a distance class) |
| 90-degree integer-K floor (this paper) | $0/90/180/270$ (4-valued) | **integer matrix** | $\mu=-\sigma_{\max}$ eigenvector | breaks distance-class symmetry (amplitude per parity class) |

The 90-degree integer-K floor is not distance-class symmetric (amplitude per parity class: 2-valued for even $N$, 3-valued for odd $N$), so it cannot be built by distance-class reduction; the correct approach is to solve it as the Perron mode of the bipartite $W$ (= the $\sigma_{\max}$ eigenmode of the integer $K$). The two floors are distinct members, but both satisfy the ignition qualification of Paper 2 (an unstable relative equilibrium).

## 8. Claim Classification · Open Items

- Classification: analytic exact solution (parity bipartite structure + reduced Perron problem) + symbolic (arbitrary $N$, sympy) / numerical ($N=3$–$40$) verification. Verdict: retained.
- **Resolved (given closed form / general proof in this version)**: (1) The edge assignment of the amplitude classes is fixed to be exactly the parity classes themselves (cross-parity = $i+j$ odd / even-side internal = both $i,j$ even / odd-side internal = both $i,j$ odd) (§5). (2) The residual centroid for odd $N$ is also in closed form ($D=N^2-2N-1$):

$$|\textstyle\sum z|=\begin{cases}0&N\equiv0\ (\mathrm{mod}\,4)\\[1mm] \sqrt2/N&N\equiv2\ (\mathrm{mod}\,4)\\[1mm] \sqrt{(N-1)/(2(N+1)D)}&N\equiv1\ (\mathrm{mod}\,4)\\[1mm] \sqrt{(N+1)/(2(N-1)D)}&N\equiv3\ (\mathrm{mod}\,4)\end{cases}$$

Matches the measured values for all $N=3$–$40$ (error $\le10^{-9}$; e.g. $N=5$: 0.154303, $N=6$: 0.235702, $N=7$: 0.140028). (3) Square-closure $\sum z^2=0$ is a consequence of the bipartite eigenmode, a general theorem holding not only for the Perron mode but for every nonzero eigenmode (§2). (4) **The connectivity of the bipartite $W$ is proven in general for arbitrary $N\ge3$ (recovering the open item of Paper 1)**: even vertices $U$ ($|U|=\lceil N/2\rceil\ge2$), odd vertices $V$. Cross-parity edges $(u,v),(u',v')$ are connected by $(u,v)$–$(u,u')$–$(u',v')$ ($(u,u')$ is even-side internal) if $u\neq u'$, and by $(u,v)$–$(v,v')$–$(u,v')$ if $u=u'$ ($v\neq v'$, $|V|\ge2$). A same-parity edge $(u,u')$ connects to its component via $(u,u')$–$(u,v)$ (any $v\in V$). Hence for $N\ge3$ it is connected = irreducible, and Perron–Frobenius is applicable for arbitrary $N\ge3$ (numerically confirmed for $N=3$–$40$ that the number of connected components $=1$). (5) Enumeration of the sign-gauge family: independently for each edge $k_e\to k_e+2t_e$ ($S=\operatorname{diag}((-1)^{t_e})$, $K\to SKS$ with $W$ invariant) gives $2^M$ possibilities, and identifying the all-edge simultaneous flip $S=-I$ (global phase $\pi$) yields a relative sign gauge of $2^{M-1}$.
- Open: within the scope of this paper no structural open items remain. Classification of floors by equivalences other than the sign gauge, such as vertex permutations (organizing the geometric moduli), is outside the scope of this paper.

## Related Work (Structural Agreement · Not the Derivation Source)

- Spectra of line graphs (Cvetković et al., algebraic graph theory): the eigenvalues of $A=L(K_N)$ and the integer antisymmetric operator on it.
- Relative equilibria (Marsden & Ratiu): the floor = a rigid-rotation relative equilibrium.

## References

1. Noriaki Kihara, "The Mechanism of Inflationary Rapid Expansion in Self-Consistent Relational-Wave Closed Systems," Concept DOI 10.5281/zenodo.22112008.
2. D. Cvetković, P. Rowlinson, S. Simić, *An Introduction to the Theory of Graph Spectra*, Cambridge (2010).

## Reproduction

The generator, verification, and data are in `90度理論床_整数K解析解_全N_20260909/` (`make_90deg_floor_analytic_v1.py` = generator, `closed_form_90deg_floor_verify_v1.py` = closed-form sympy exact verification, `parents_90deg_floor_analytic/parent_90deg_floor_N{00003..00040}_analytic.npz` = 38 floors, `summary_90deg_floor_closed_form.csv`, `manifest_90deg_floor_analytic.json`, `fig_90deg_floor_complex_planes_analytic.png`, `分析_90度理論床_整数K解析解_全N_20260909.md`, `run_all.sh`, `SHA256SUMS.txt`). The verification of the arbitrary-$N$ analytic reduction via the parity bipartite structure (§2 · §3 · §5 · §6 · §8) is in `同/parity二部構造_解析的証明_20260912/` (`verify_parity_bipartite_analytic_20260912.py` = read-only verification importing the canonical adjacency/one_step with SHA check + sympy symbolic, `results/parity_bipartite_summary.csv`, `results/実行ログ_20260912.log`, `分析_parity二部構造_解析的一本化_20260912.md`, `README.md`, `run_all.sh`, `SHA256SUMS.txt`). The eigen-oscillation spectrum is in `同/固有振動数解析_90度系列_N3_N40_20260910/`. The equivalence of the complex self-consistent problem ≡ the real symmetric $W$ ($D^{-1}KD=iW$) and Perron–Frobenius are in Paper 1. adjacency/one_step are extracted from the canonical source `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`) with no physics changed.
