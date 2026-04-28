# 論文3 §2–§4 試案

**仮タイトル**: Discrete Structure of a Four-Dimensional Ball: Unit-Cube Packing and the Asymptotic Volume Deficit $\Delta(R)$

---

## §2. Problem Setting

### §2.1 The four-dimensional ball

Let $R > 0$ be a real number. The four-dimensional ball of radius $R$ centred at the origin is the closed set
$$B(R) := \left\{x = (x_1, x_2, x_3, x_4) \in \mathbb{R}^4 : \|x\|_2 \le R\right\},$$
where $\|x\|_2 = (x_1^2 + x_2^2 + x_3^2 + x_4^2)^{1/2}$.

Its volume and the volume of its boundary $\partial B(R) \simeq S^3(R)$ are
$$V_4(R) = \frac{\pi^2}{2} R^4, \qquad S_3(R) = 2\pi^2 R^3.$$

We will be concerned with the asymptotic behaviour of certain integer counts associated with $B(R)$, and these volumes will appear repeatedly as natural normalizations.

### §2.2 The unit cube packing problem

A *unit cube* in $\mathbb{R}^4$ centred at a point $c \in \mathbb{R}^4$ is the closed set
$$Q(c) := \left\{x \in \mathbb{R}^4 : |x_i - c_i| \le \tfrac{1}{2}, \;\; i = 1, 2, 3, 4\right\}.$$

The cube has side length $1$ and four-dimensional volume $1$. Its $2^4 = 16$ vertices are the points $c + (\pm 1/2, \pm 1/2, \pm 1/2, \pm 1/2)$.

We restrict attention to cubes whose centres lie at integer lattice points: $c \in \mathbb{Z}^4$. Two questions then arise:

**(a)** For which $c \in \mathbb{Z}^4$ is the unit cube $Q(c)$ fully contained in the ball $B(R)$?

**(b)** How many such cubes are there, as a function of $R$?

We denote the answer to (b) by $N_4(R)$:
$$N_4(R) := \#\left\{c \in \mathbb{Z}^4 : Q(c) \subseteq B(R)\right\}.$$

This is the number of integer-centred unit cubes that fit inside $B(R)$ without protrusion.

### §2.3 Restriction to odd diameter

For combinatorial convenience we restrict $R$ to take the values
$$R_k := 2k + 1, \qquad k \in \mathbb{Z}_{\ge 0}.$$

This choice has two motivations:

1. *Symmetry.* The condition $Q(c) \subseteq B(R)$ is symmetric under $c_i \mapsto -c_i$ for any axis. With $R = 2k+1$, the centre $c = 0$ is included for all $k \ge 0$ (the unit cube around the origin fits in $B(1)$), and the lattice points come in symmetric multiplets without ambiguity.

2. *Integer arithmetic.* Multiplying the packing inequality by $4$ converts all quantities to integers, allowing exact computation.

Throughout this paper, $R$ denotes $R_k = 2k + 1$ unless otherwise stated. We write $N(k) := N_4(R_k)$ for the count.

The restriction is not essential: one may extend to all integer $R$ or to continuous $R$ with minor modifications. The asymptotic results we derive hold without modification in either extension; the only change is that the data points become more closely spaced.

---

## §3. Formulation of the Packing Condition

### §3.1 The containment inequality

A cube $Q(c)$ centred at $c \in \mathbb{R}^4$ is contained in $B(R)$ if and only if all $16$ of its vertices lie in $B(R)$. Among the $16$ vertices, the one farthest from the origin is
$$P_{\max}(c) = \left(|c_1| + \tfrac{1}{2}, |c_2| + \tfrac{1}{2}, |c_3| + \tfrac{1}{2}, |c_4| + \tfrac{1}{2}\right),$$
since the absolute value of each coordinate is maximized by adding $1/2$ in the same sign as $c_i$ (or, if $c_i = 0$, by choosing either sign).

The containment condition $Q(c) \subseteq B(R)$ is therefore equivalent to $\|P_{\max}(c)\|_2 \le R$:

**Lemma 3.1.** *Let $c \in \mathbb{Z}^4$ and $R > 0$. Then $Q(c) \subseteq B(R)$ if and only if*
$$\boxed{\;\sum_{i=1}^{4} \left(|c_i| + \tfrac{1}{2}\right)^2 \le R^2.\;}$$

*Proof.* The vertices of $Q(c)$ are $c + (\epsilon_1/2, \epsilon_2/2, \epsilon_3/2, \epsilon_4/2)$ for $\epsilon_i \in \{-1, +1\}$. The squared norm of such a vertex is
$$\sum_i (c_i + \epsilon_i/2)^2 = \sum_i (c_i^2 + \epsilon_i c_i + 1/4),$$
which is maximized over the choice of $\epsilon$ by taking $\epsilon_i = \mathrm{sgn}(c_i)$ when $c_i \neq 0$ (and either sign when $c_i = 0$). The maximum value is $\sum_i (|c_i| + 1/2)^2$. The condition that all vertices lie in $B(R)$ is therefore that this maximum does not exceed $R^2$. $\square$

### §3.2 Integer reformulation

Setting $p_i := 2|c_i| + 1$, the inequality of Lemma 3.1 becomes
$$\sum_{i=1}^{4} \left(\frac{p_i}{2}\right)^2 \le R^2, \qquad \text{i.e.,} \qquad \sum_{i=1}^{4} p_i^2 \le (2R)^2.$$

For $R = 2k+1$, this is
$$\sum_{i=1}^{4} p_i^2 \le (4k+2)^2,$$
where each $p_i \in \{1, 3, 5, \ldots\}$ is a positive odd integer.

This integer reformulation has two consequences. First, it places the packing problem within the ambit of classical sum-of-squares number theory, which we exploit in Section 6. Second, it allows exact computation in arbitrary-precision integer arithmetic, with no floating-point error.

### §3.3 The count $N(k)$

Combining Lemma 3.1 with the conversion between $c \in \mathbb{Z}^4$ and $(p_1, p_2, p_3, p_4)$ with $p_i \in \{1, 3, 5, \ldots\}$, we have

$$N(k) = \sum_{\substack{c \in \mathbb{Z}^4 \\ \sum (|c_i| + 1/2)^2 \le R_k^2}} 1.$$

Equivalently, in terms of the odd positive integers $p_i$:

$$N(k) = \sum_{\substack{(p_1, p_2, p_3, p_4) \in \{1,3,5,\ldots\}^4 \\ \sum p_i^2 \le (4k+2)^2}} 2^{\#\{i : p_i > 1\}}.$$

The factor $2^{\#\{i : p_i > 1\}}$ accounts for the sign choice in $c_i = \pm (p_i - 1)/2$ when $p_i > 1$ (i.e., $c_i \neq 0$); when $p_i = 1$, we have $c_i = 0$ and there is no sign choice.

### §3.4 Computational considerations

The direct enumeration of $\mathbb{Z}^4$ in the ball $B(R)$ has cost $O(R^4)$. For $R = 2k+1$ and $k$ in the range we consider ($k \le 60$, so $R \le 121$), this is computationally tractable but inefficient.

A more efficient enumeration restricts to ordered non-negative tuples $y_1 \ge y_2 \ge y_3 \ge y_4 \ge 0$, treats sign and permutation multiplicities analytically, and reduces the cost by a factor of $|B_4| = 384$ (the order of the hyperoctahedral symmetry group of the cube). Specifically:

$$N(k) = \sum_{\substack{y_1 \ge y_2 \ge y_3 \ge y_4 \ge 0 \\ \sum (y_i + 1/2)^2 \le R_k^2}} \mu_{\mathrm{sign}}(y) \cdot \mu_{\mathrm{perm}}(y),$$

where
$$\mu_{\mathrm{sign}}(y) = 2^{\#\{i : y_i > 0\}}, \qquad \mu_{\mathrm{perm}}(y) = \frac{4!}{\prod_v |\{i : y_i = v\}|!}.$$

This enumeration is implemented in the accompanying code and produces exact values of $N(k)$ for all $k \le 60$ in under one minute on commodity hardware.

---

## §4. The Sequence $N(k)$: Combinatorial Properties

### §4.1 Initial values

Direct computation gives the initial values:

| $k$ | $R = 2k+1$ | $N(k)$ | $a_k := N(k) - N(k-1)$ |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 |
| 1 | 3 | 137 | 136 |
| 2 | 5 | 1,545 | 1,408 |
| 3 | 7 | 7,281 | 5,736 |
| 4 | 9 | 22,409 | 15,128 |
| 5 | 11 | 53,161 | 30,752 |
| 6 | 13 | 108,081 | 54,920 |
| 7 | 15 | 199,953 | 91,872 |
| 8 | 17 | 337,417 | 137,464 |
| 9 | 19 | 537,409 | 199,992 |
| 10 | 21 | 818,145 | 280,736 |

### §4.2 Boundary cubes

**Proposition 4.1.** *For every $k \ge 0$, the $16$ unit cubes centred at $c = (\pm k, \pm k, \pm k, \pm k)$ have their farthest vertex on the sphere $\partial B(R_k) = \partial B(2k+1)$. These are the "corner cubes" of the configuration.*

*Proof.* For $c = (\pm k)^4$, $|c_i| = k$ for all $i$, so $|c_i| + 1/2 = k + 1/2$. Then
$$\sum_i (|c_i| + 1/2)^2 = 4 \cdot (k + 1/2)^2 = (2k + 1)^2 = R_k^2.$$
By Lemma 3.1, the inequality is satisfied with equality, meaning the farthest vertex lies exactly on $\partial B(R_k)$. There are $2^4 = 16$ such corner cubes (one in each orthant). $\square$

**Remark.** The exact incidence of the corner cubes' vertices with $\partial B(R_k)$ is special to four dimensions and to the choice $R = 2k+1$. In dimension $n$, the analogous condition reads
$$n \cdot (k + 1/2)^2 = R^2,$$
which has an integer solution $R = (2k+1) \sqrt{n}/2$ only when $\sqrt{n}$ is rational, hence only for $n$ a perfect square. The smallest non-trivial dimension with this property is $n = 4$, in which case $\sqrt{n} = 2$ and $R = 2k+1$.

### §4.3 Asymptotic packing density

Define the *packing density*
$$\rho(k) := \frac{N(k)}{V_4(R_k)} = \frac{N(k)}{(\pi^2/2) (2k+1)^4}.$$

This is the fraction of the four-dimensional ball $B(R_k)$ that is occupied by the union of the $N(k)$ unit cubes.

**Proposition 4.2.** *$\lim_{k \to \infty} \rho(k) = 1$.*

*Proof.* The condition $\sum (|c_i| + 1/2)^2 \le R^2$ defines a region $\Omega_R \subset \mathbb{R}^4$. By the standard Gauss-style lattice point theorem for convex bodies, the number of integer points in $\Omega_R$ satisfies
$$N(k) = \mathrm{Vol}(\Omega_R) + O(R^3) \quad \text{as } R \to \infty.$$
Direct computation (deferred to §5) gives $\mathrm{Vol}(\Omega_R) = V_4(R) - O(R^3)$. Therefore $N(k) / V_4(R) \to 1$. $\square$

The numerical values of $\rho(k)$ confirm the convergence:

| $k$ | $\rho(k)$ |
|---:|---:|
| 0 | 0.2026 |
| 1 | 0.3427 |
| 2 | 0.5009 |
| 5 | 0.7358 |
| 10 | 0.8525 |
| 20 | 0.9201 |
| 30 | 0.9457 |
| 40 | 0.9562 |
| 50 | 0.9669 |
| 60 | 0.9723 |

The convergence is not uniformly rapid: at $k = 60$, the density is still only $\approx 0.972$. The slow convergence is governed by the boundary correction whose precise form is the subject of §5.

### §4.4 Volume deficit and its scaling

Define the *volume deficit*
$$\Delta(R) := V_4(R) - N(k).$$

This quantity measures the "missing volume": the difference between the continuous four-dimensional ball volume and the integer count of unit cubes that fit inside.

The principal result of this paper, proved in §5–§6, is:

$$\Delta(R) = \frac{16\pi}{3} R^3 - 6\pi R^2 + O(R) \qquad \text{as } R \to \infty.$$

In normalized form, with $c(k) := \Delta(R) / (2\pi^2 R^3)$, the leading constant is
$$\boxed{\;c := \lim_{k \to \infty} c(k) = \frac{16\pi/3}{2\pi^2} = \frac{8}{3\pi} \approx 0.84883.\;}$$

The numerical values of $c(k)$ approach this limit from below:

| $k$ | $c(k)$ |
|---:|---:|
| 10 | 0.7745 |
| 20 | 0.8192 |
| 30 | 0.8276 |
| 40 | 0.8311 |
| 50 | 0.8350 |
| 60 | 0.8380 |

The convergence is slow, but a polynomial fit confirms the leading constant to within $0.024\%$ for $k \ge 30$ (see §5.4).

### §4.5 Generating function (deferred)

A generating-function representation of $N(k)$ in terms of the Jacobi theta function will be developed in §6. We do not need the explicit form for the asymptotic analysis of §5. The connection serves a separate purpose: it places $N(k)$ within the classical framework of sum-of-squares representations and provides an exact closed form (modulo a sign-correction factor whose precise determination is deferred to §6.3).

---

**(§5–§7 は別途執筆)**

---

**Style notes**:
- 自己引用ゼロ
- 用語は本論文内で完結的に定義
- 命題と証明を学術論文スタイルで記述
- 「The smallest non-trivial dimension with this property is $n = 4$」は数学的記述として論文8 Rev.2 と同様の内容を含むが、引用なしで再導出
- §5（漸近解析）, §6（Jacobi接続）, §7（結論）は別ファイルで執筆予定

**未確定事項（後で検討）**:
- §5の漸近解析の正式な記述（math_report.md の内容を学術論文形式に整形）
- §6の Jacobi 接続の正確な係数
- §1（序論）と §7（結論）
