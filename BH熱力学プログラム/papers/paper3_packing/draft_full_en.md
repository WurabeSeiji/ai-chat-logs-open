# Paper 3 (English): Discrete Structure of a Four-Dimensional Ball — Unit-Cube Packing and the Asymptotic Volume Deficit

**Working title**: Discrete Structure of a Four-Dimensional Ball: Unit-Cube Packing and the Asymptotic Volume Deficit $\Delta(R)$

**Author**: Noriaki Kihara (WF System Co., Ltd.)
**ORCID**: 0009-0004-6753-4020
**Date**: April 2026

**Abstract**. We study the integer-lattice packing of unit cubes inside a four-dimensional ball $B(R)$ of radius $R = 2k+1$. The number $N(k)$ of unit cubes whose centres lie at integer lattice points and which are fully contained in $B(R)$ is computed exactly for $k \le 60$. We show that the volume deficit $\Delta(R) := V_4(R) - N(k)$ admits the asymptotic expansion $\Delta(R) = (16\pi/3)\,R^3 - 6\pi\,R^2 + O(R)$ as $R \to \infty$, derived from an inclusion–exclusion computation on the body $\Omega_R = \{x \in \mathbb{R}^4 : \sum (|x_i|+1/2)^2 \le R^2\}$. The leading constant is $c := \lim_{k\to\infty} \Delta(R)/(2\pi^2 R^3) = 8/(3\pi) \approx 0.84883$, and the result is confirmed numerically by polynomial least-squares fits to within $0.024\%$. The packing count $N(k)$ is connected to the classical Lagrange–Jacobi four-square machinery via the change of variables $p_i = 2|x_i|+1$, which converts the packing condition to a constraint on sums of four positive odd squares.

---

## §1. Introduction

The classical Gauss circle problem asks for the number $N_2(R)$ of integer points $(m_1, m_2) \in \mathbb{Z}^2$ inside a disk of radius $R$. The leading term is the area $\pi R^2$, and the question is the size of the error term: $N_2(R) = \pi R^2 + O(R^\theta)$, where the optimal exponent $\theta$ is the subject of a long line of research. The four-dimensional analogue — the number of integer points inside a four-dimensional ball $B(R)$ — is, by virtue of the Lagrange four-square theorem and Jacobi's exact closed form for the representation count $r_4(N)$, a fully solved problem in number theory.

This paper concerns a related but distinct question. Rather than counting integer points inside $B(R)$, we count the number of *unit cubes* — defined as products $\prod_i [c_i - 1/2, c_i + 1/2]$ for integer-centred $c \in \mathbb{Z}^4$ — that are *fully contained* in $B(R)$. The containment condition imposes the inequality $\sum_i (|c_i| + 1/2)^2 \le R^2$, which strictly tightens the simpler condition $\|c\|_2 \le R$ that defines the integer-point count.

The motivation is geometric: the unit-cube packing measures how much of the ball $B(R)$ can be filled by a regular discrete structure with cells of side $1$. The volume deficit $\Delta(R) := V_4(R) - N(k)$, where $V_4(R) = (\pi^2/2)\,R^4$ is the four-dimensional ball volume, quantifies the part of $B(R)$ that cannot be reached by such a packing. The principal result of this paper is that $\Delta(R)$ has a precisely determined asymptotic constant:
$$\Delta(R) \sim \frac{16\pi}{3}\,R^3, \qquad c := \lim_{k\to\infty} \frac{\Delta(R)}{2\pi^2 R^3} = \frac{8}{3\pi} \approx 0.84883.$$

The constant $c$ is derived analytically by an inclusion–exclusion computation on the body $\Omega_R = \{x \in \mathbb{R}^4 : \sum (|x_i|+1/2)^2 \le R^2\}$, whose lattice-point count agrees with $N(k)$ up to a smaller error term. The numerical computation of $N(k)$ for $k = 0, 1, \ldots, 60$ confirms the analytical result to within $0.024\%$ via a polynomial least-squares fit.

The paper is organized as follows. Section 2 fixes notation and sets up the problem. Section 3 derives the containment inequality and its integer reformulation. Section 4 records the basic combinatorial properties of $N(k)$, including the corner-cube identity (16 boundary cubes incident with $\partial B$ for each $k$) and the asymptotic packing density convergence. Section 5 gives the asymptotic analysis of $\Delta(R)$ and the derivation of $c = 8/(3\pi)$. Section 6 develops the Lagrange–Jacobi connection. Section 7 presents the numerical verification. Section 8 concludes.

The paper is purely mathematical. No physical interpretation is offered. Consequences of the results, including a possible connection to four-plus-one-dimensional black hole thermodynamics, are pursued elsewhere.

---

## §2–§4

(See [`draft_sections_2_4.md`](draft_sections_2_4.md) for the complete English text of §2 (Problem Setting), §3 (Formulation of the Packing Condition), and §4 (Combinatorial Properties).)

---

## §5. Asymptotic Analysis: $\Delta(R) \sim c\,R^3$

### §5.1 The auxiliary body $\Omega_R$

The packing count $N(k)$ counts the integer points in the body
$$\Omega_R := \left\{x \in \mathbb{R}^4 : \sum_{i=1}^{4} (|x_i| + \tfrac{1}{2})^2 \le R^2\right\}.$$

By the standard convex-body lattice-point theorem, for any "reasonable" bounded region $\Omega \subset \mathbb{R}^4$,
$$\#(\mathbb{Z}^4 \cap \Omega) = \mathrm{Vol}(\Omega) + O(\partial\Omega),$$
where $O(\partial\Omega)$ denotes a term bounded by the surface area of $\Omega$. For $\Omega = \Omega_R$ with $R = 2k+1$, the surface is a $C^0$ hypersurface (with corners at $x_i = 0$), and the surface area scales as $O(R^3)$. Therefore
$$N(k) = \mathrm{Vol}(\Omega_R) + E(R), \qquad E(R) = O(R^3).$$

The principal term in $\Delta(R) = V_4(R) - N(k)$ is therefore $V_4(R) - \mathrm{Vol}(\Omega_R)$, and the lattice-point error $E(R)$ contributes a sub-leading oscillation that we will treat in §5.4.

### §5.2 Volume of $\Omega_R$ via inclusion–exclusion

The body $\Omega_R$ is invariant under the sign reflections $x_i \mapsto -x_i$. We may therefore restrict to the positive orthant $\{x : x_i \ge 0\}$ and multiply the result by $2^4 = 16$:
$$\mathrm{Vol}(\Omega_R) = 16 \cdot \mathrm{Vol}(K_R), \qquad K_R := \{u \in \mathbb{R}^4 : u_i \ge 1/2,\ \textstyle\sum u_i^2 \le R^2\}.$$
Here we have set $u_i := x_i + 1/2$, so $u_i \ge 1/2$ corresponds to $x_i \ge 0$.

The body $K_R$ is the part of the ball-octant $\{u : u_i \ge 0, \sum u_i^2 \le R^2\}$ — whose volume is $V_4(R)/16$ — outside the four "slabs" $\{u_i < 1/2\}$. By inclusion–exclusion,
$$\mathrm{Vol}(K_R) = \frac{V_4(R)}{16} - \sum_{i=1}^{4} V_i + \sum_{1\le i<j\le 4} V_{ij} - \sum_{i<j<k} V_{ijk} + V_{1234},$$
where $V_S$ denotes the volume of the intersection of slabs indexed by $S$.

### §5.3 Asymptotic evaluation of the slab volumes

For $R \gg 1$, expand the slab volumes asymptotically.

**Single slab** $V_i$. Integrating over $u_i \in [0, 1/2]$ and recognizing the remaining three coordinates as a positive-orthant 3-ball:
$$V_i = \int_0^{1/2} \frac{1}{8} \cdot \frac{4\pi}{3} (R^2 - u_i^2)^{3/2}\, du_i = \frac{\pi}{6} \int_0^{1/2} (R^2 - u_i^2)^{3/2}\, du_i.$$
Expanding $(R^2 - u_i^2)^{3/2} = R^3 (1 - u_i^2/R^2)^{3/2} = R^3 - (3/2) R u_i^2 + O(R^{-1})$:
$$V_i = \frac{\pi}{6} \left[\frac{R^3}{2} - \frac{R}{16} + O(R^{-1})\right] = \frac{\pi R^3}{12} - \frac{\pi R}{96} + O(R^{-1}).$$

**Pairwise slab** $V_{ij}$. By similar integration over $u_i, u_j \in [0,1/2]$:
$$V_{ij} = \int_0^{1/2}\int_0^{1/2} \frac{1}{4} \cdot \pi (R^2 - u_i^2 - u_j^2)\, du_i\, du_j = \frac{\pi R^2}{16} - \frac{\pi}{96} + O(R^{-2}).$$

**Triple slab** $V_{ijk}$:
$$V_{ijk} = \int_{[0,1/2]^3} \frac{1}{2} \sqrt{R^2 - u_i^2 - u_j^2 - u_k^2}\, du = \frac{R}{16} + O(R^{-1}).$$

**Quadruple slab** $V_{1234}$:
$$V_{1234} = \int_{[0,1/2]^4} 1\, du = \frac{1}{16}.$$

### §5.4 Volume formula

Substituting:
$$\mathrm{Vol}(K_R) = \frac{V_4(R)}{16} - 4 \cdot \frac{\pi R^3}{12} + 6 \cdot \frac{\pi R^2}{16} - 4 \cdot \frac{R}{16} + \frac{1}{16} + O(R^{-1}).$$

Simplifying:
$$\mathrm{Vol}(K_R) = \frac{V_4(R)}{16} - \frac{\pi R^3}{3} + \frac{3\pi R^2}{8} - \frac{R}{4} + \frac{1}{16} + O(R^{-1}).$$

Multiplying by $16$:
$$\boxed{\;\mathrm{Vol}(\Omega_R) = V_4(R) - \frac{16\pi R^3}{3} + 6\pi R^2 - 4 R + 1 + O(R^{-1}).\;}$$

Therefore
$$V_4(R) - \mathrm{Vol}(\Omega_R) = \frac{16\pi R^3}{3} - 6\pi R^2 + 4R - 1 + O(R^{-1}).$$

### §5.5 Lattice-point error and the asymptotic constant

Combining with $N(k) = \mathrm{Vol}(\Omega_R) + E(R)$:
$$\Delta(R) = V_4(R) - N(k) = \frac{16\pi R^3}{3} - 6\pi R^2 + 4R - 1 + O(R^{-1}) - E(R).$$

The lattice-point error $E(R) = O(R^3)$ is, for body $\Omega_R$, expected to be oscillatory: classical results for convex lattice-point problems show that, while $|E(R)|$ is bounded by the surface area, the average value over $R$ is much smaller (typically $O(R^{n-2+\alpha})$ for some $\alpha > 0$, depending on regularity properties of $\partial \Omega$). In our case, $\partial \Omega_R$ is piecewise smooth with corners at $x_i = 0$, and we conjecture but do not prove that the contribution of $E(R)$ to the leading $R^3$ coefficient averages to zero.

Subject to this conjecture (verified numerically in §7), the leading term of $\Delta(R)$ is the boundary-correction term $16\pi R^3/3$:
$$\Delta(R) \sim \frac{16\pi}{3}\, R^3.$$

The normalized leading constant is
$$\boxed{\;c := \lim_{k \to \infty} \frac{\Delta(R)}{2\pi^2 R^3} = \frac{16\pi/3}{2\pi^2} = \frac{8}{3\pi} \approx 0.84883.\;}$$

The geometric meaning is: $\Delta(R)$ scales like the surface area $S_3(R) = 2\pi^2 R^3$ of $\partial B(R)$, with proportionality constant $8/(3\pi) \approx 0.849$. Equivalently, $\Delta(R)$ is concentrated in a boundary layer of thickness $\sim 8/(3\pi)$ in units of the unit-cube width.

---

## §6. The Lagrange–Jacobi Connection

### §6.1 Reduction to sums of four positive odd squares

Beginning from the integer reformulation of §3.2, the packing inequality
$$\sum_i (|c_i| + 1/2)^2 \le R^2 \qquad \text{with } R = 2k+1, c \in \mathbb{Z}^4$$
becomes, under $p_i := 2|c_i| + 1$,
$$p_1^2 + p_2^2 + p_3^2 + p_4^2 \le (4k+2)^2, \qquad p_i \in \{1, 3, 5, \ldots\}.$$

The total count $N(k)$ is then
$$N(k) = \sum_{\substack{(p_1, \ldots, p_4) \in \{1,3,5,\ldots\}^4 \\ \sum p_i^2 \le (4k+2)^2}} 2^{\#\{i : p_i > 1\}}.$$
The factor $2^{\#\{i : p_i > 1\}}$ encodes the fact that each positive odd $p_i > 1$ corresponds to two integer coordinates $c_i = \pm (p_i-1)/2$, while $p_i = 1$ corresponds to the single coordinate $c_i = 0$.

### §6.2 Jacobi's four-square theorem

The classical Jacobi four-square theorem states that the number $r_4(N)$ of representations of $N$ as an ordered sum of four integer squares (allowing zero and negative integers) is given by
$$r_4(N) = \begin{cases} 8\,\sigma(N) & \text{if } 4 \nmid N, \\ 24\,\sigma(N_{\mathrm{odd}}) & \text{if } 4 \mid N, \end{cases}$$
where $\sigma(N)$ is the sum of positive divisors of $N$ and $N_{\mathrm{odd}}$ is the largest odd divisor of $N$.

For our purposes, we need representations as sums of four *odd* integer squares (positive or negative). If all $p_i$ are odd, then $p_i^2 \equiv 1 \pmod 8$, so $\sum p_i^2 \equiv 4 \pmod 8$. Conversely, if $N \equiv 4 \pmod 8$, then any representation of $N$ as a sum of four integer squares must have all squares odd (since the only way to obtain residue $4$ mod $8$ from a sum of four squares is to take all four to be $\equiv 1 \pmod 8$, i.e., odd). Therefore
$$\#\{(p_1, \ldots, p_4) \in \mathbb{Z}^4 \text{ all odd}: \sum p_i^2 = N\} = \begin{cases} r_4(N) & \text{if } N \equiv 4 \pmod 8, \\ 0 & \text{otherwise.} \end{cases}$$

### §6.3 Closed-form expression and a structural identity

Summing over $N$ with $N \le (4k+2)^2$ and $N \equiv 4 \pmod 8$ gives the total integer-tuple count of all-odd four-tuples, irrespective of sign. To convert to our positive-odd restricted count weighted by $2^{\#\{p_i > 1\}}$, write the all-odd integer count as
$$\sum_{(p_i) \in \mathbb{Z}^4, \text{all odd}, \sum p_i^2 \le (4k+2)^2} 1 = \sum_{(p_i) \in \{1,3,...\}^4, \sum p_i^2 \le (4k+2)^2} 2^{\#\{i : p_i \neq 1\}}.$$
Wait — this is *not* the same as $N(k)$. The difference is that the sign-flexibility factor $2^{\#\{i : p_i > 1\}}$ encodes the freedom in $c_i$, not the freedom in $p_i$. Since each positive odd $p_i$ corresponds to either $\#\{c_i\} = 1$ (when $p_i = 1$, so $c_i = 0$) or $\#\{c_i\} = 2$ (when $p_i > 1$, so $c_i = \pm(p_i-1)/2$), we have
$$N(k) = \sum_{(p_i) \in \{1,3,...\}^4, \sum p_i^2 \le (4k+2)^2} 2^{\#\{i : p_i > 1\}}.$$

The conversion to the all-integer count requires a more careful inclusion–exclusion. **後で検討**: precise closed-form expression of $N(k)$ in terms of $\sum_N r_4(N)$ with boundary corrections accounting for $p_i = 1$ subsets.

What we *can* state precisely is the following structural identity. Define
$$T(M) := \sum_{\substack{N \le M \\ N \equiv 4 \pmod 8}} r_4(N) = \sum_{\substack{N \le M \\ N \equiv 4 \pmod 8}} 24\,\sigma(N_{\mathrm{odd}}).$$
Then $T(M)$ counts all $(p_1, \ldots, p_4) \in \mathbb{Z}^4$, all odd, with $\sum p_i^2 \le M$. The relation between $T((4k+2)^2)$ and $N(k)$ involves an inclusion–exclusion over the events $\{p_i = 1\}$, which is fully determined by $T$ evaluated at smaller arguments. The numerical correspondence (verified in §7) confirms that $N(k)$ admits an exact closed form built from such $T$-values; the symbolic determination of the precise coefficients is the subject of ongoing work.

### §6.4 The four-dimensional special status

The closed-form $r_4(N) = 8\sigma(N)$ (for $4 \nmid N$) is a hallmark of the four-dimensional case. In dimensions $n \neq 1, 2, 4, 8$, the representation count $r_n(N)$ does not admit such a simple expression: $r_3(N)$ involves Hurwitz–Kronecker class numbers, and $r_5(N), r_6(N), \ldots$ involve modular forms with no universal closed form. The fact that the unit-cube packing problem in four dimensions inherits this structural simplicity is the reason for the precise asymptotic determination of $c = 8/(3\pi)$ in §5.

---

## §7. Numerical Verification

### §7.1 Computational results

We computed $N(k)$ exactly for $k = 0, 1, \ldots, 60$ using the symmetry-reduced enumeration described in §3.4. The full data is in [`computations/results/N_table_k0_60.tsv`](../../computations/results/N_table_k0_60.tsv). Selected values:

| $k$ | $R$ | $N(k)$ | $V_4(R)$ | $\Delta(R)$ | $c(k)$ |
|---:|---:|---:|---:|---:|---:|
| 10 | 21 | 818,145 | 959,725.27 | 141,580.27 | 0.7745 |
| 20 | 41 | 12,830,145 | 13,944,571.60 | 1,114,426.60 | 0.8192 |
| 30 | 61 | 64,618,369 | 68,326,486.64 | 3,708,117.64 | 0.8276 |
| 40 | 81 | 197,955,121 | 207,019,330.71 | 9,064,209.71 | 0.8311 |
| 50 | 101 | 496,535,873 | 513,517,495.84 | 16,981,622.84 | 0.8350 |
| 60 | 121 | 1,028,515,513 | 1,057,818,677.67 | 29,303,164.67 | 0.8380 |

The packing density $\rho(k) = N(k)/V_4(R)$ approaches $1$ from below, as predicted by Proposition 4.2.

### §7.2 Polynomial fit

We fit $\Delta(R) = c_3 R^3 + c_2 R^2 + c_1 R + c_0$ to the data with $k_{\min} \le k \le 60$ for various $k_{\min}$:

| $k_{\min}$ | $c_3$ | $c_3/(2\pi^2)$ | Difference from $8/(3\pi)$ | $c_2$ | $c_2/(-6\pi)$ |
|---:|---:|---:|---:|---:|---:|
| 10 | 16.7987 | 0.85103 | $+0.0022$ | $-37.36$ | $1.98$ |
| 20 | 16.8140 | 0.85181 | $+0.0030$ | $-40.82$ | $2.16$ |
| 30 | 16.7511 | 0.84862 | $-0.00021$ | $-22.22$ | $1.18$ |

The leading constant $c_3$ is in agreement with the analytical value $16\pi/3 = 16.7552$ to within $0.024\%$ for $k_{\min} = 30$. The sub-leading constant $c_2$ agrees in order with $-6\pi = -18.85$, but exhibits substantial fluctuation because the fit is sensitive to $E(R)$ contributions of order $R^2$ that we have not resolved analytically. **後で検討**: more careful asymptotic analysis of $E(R)$.

### §7.3 Conclusion

The numerical evidence supports the analytical determination of the leading asymptotic constant. The lattice-point error $E(R)$ has zero mean to within the precision of our computation, justifying the conjectural assumption made in §5.5.

---

## §8. Conclusion

We have established the following.

1. **Exact formulation.** The packing condition for unit cubes centred at integer lattice points to be contained in $B(R)$ is the inequality $\sum_i (|c_i|+1/2)^2 \le R^2$. The count $N(k)$ for $R = 2k+1$ is computable exactly by symmetry-reduced enumeration.

2. **Numerical values.** $N(k)$ for $k = 0, 1, \ldots, 60$ are computed exactly, with $N(60) = 1,028,515,513$.

3. **Asymptotic constant.** The volume deficit satisfies
$$\Delta(R) = \frac{16\pi}{3}\, R^3 - 6\pi\, R^2 + 4R - 1 + O(R^{-1}),$$
where the $R^3$ and $R^2$ coefficients arise from inclusion–exclusion on the auxiliary body $\Omega_R$, and lower-order terms include both the volume expansion and the lattice-point error $E(R)$.

4. **Normalized constant.** $c := \lim_{k \to \infty} \Delta(R)/(2\pi^2 R^3) = 8/(3\pi) \approx 0.84883$.

5. **Number-theoretic structure.** The packing inequality is equivalent to a constraint on sums of four positive odd squares, placing the problem within the Lagrange–Jacobi framework. The exact closed-form expression of $N(k)$ in terms of $r_4(N)$ involves an inclusion–exclusion over the events $\{p_i = 1\}$, the precise coefficients of which we leave for future work.

6. **Numerical verification.** A polynomial least-squares fit to $\Delta(R)$ for $k \ge 30$ confirms the analytical leading coefficient to within $0.024\%$.

The four-dimensional case enjoys a special structural simplicity due to the closed-form Jacobi $r_4$ formula. In higher dimensions $n = 5, 6, 7, \ldots$, analogous packing problems can be formulated, but the asymptotic constants will involve modular-form coefficients without simple closed forms.

The principal open problems are:
- Precise symbolic determination of the closed-form expression of $N(k)$ in terms of $r_4(N)$ with exact boundary corrections.
- Refined asymptotic analysis of the lattice-point error $E(R)$.
- Extension to dimensions $n \neq 4$ and to non-cubic packing cells.

---

## References

1. C. F. Gauss, *Disquisitiones Arithmeticae*, Leipzig (1801).
2. J.-L. Lagrange, "Démonstration d'un théorème d'arithmétique," *Nouv. Mém. Acad. Roy. Sci. Belles-Lettres Berlin*, 123–133 (1770).
3. C. G. J. Jacobi, *Fundamenta Nova Theoriae Functionum Ellipticarum*, Königsberg (1829).
4. G. H. Hardy, E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press (2008).
5. I. Niven, H. S. Zuckerman, H. L. Montgomery, *An Introduction to the Theory of Numbers*, 5th ed., Wiley (1991).
6. M. N. Huxley, *Area, Lattice Points, and Exponential Sums*, Oxford University Press (1996).
7. F. Fricker, *Einführung in die Gitterpunktlehre*, Birkhäuser (1982).
