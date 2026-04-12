# Realization of Perfectly Elastic Collision in the Delay Circuit Model

**Author:** Noriaki Kihara
**Affiliation:** WF System Co., Ltd. (Osaka University, School of Engineering Science, graduated)
**Date:** April 2026
**Category:** Research Note
**DOI:** [10.5281/zenodo.19534353](https://doi.org/10.5281/zenodo.19534353)

---

## Abstract

We extend class $B$ (transverse wave) defined in the previous paper [3] to define class $E$, which holds position information and displacement. After constructing a non-interacting model in which two instances of class $E$ propagate in opposite directions, we define class $F$, which adds a collision radius $r$ per instance and the perfectly elastic collision formula to the operation function $*$. By treating $w_{\max}$ as an analogue of mass, a perfectly elastic collision satisfying conservation of momentum and energy is realized. The definitions of the delay circuit $F_D$ and the NOT operation circuit are not changed.

---

## 1. Introduction

In the previous paper [1], oscillation and propagation patterns were constructed using the delay circuit $F_D$ and the NOT operation circuit. In the previous paper [2], the framework of superclasses, meta-information, and the operation function $*$ was introduced. In the previous paper [3], class $B$ was used to realize a sine wave (transverse wave).

In this paper, we define class $E$ as an extension of class $B$, and treat the case in which two instances propagate in opposite directions on the same open rope model. We first construct a model without interaction (passing through), and then construct a model that realizes a perfectly elastic collision.

---

## 2. Preparation: Summary of Definitions from Previous Papers

This paper presupposes the delay circuit $F_D$ and the NOT operation circuit from the previous paper [1], the superclass, meta-information, and operation function $*$ from the previous paper [2], and class $B$ (transverse wave) from the previous paper [3]. For details of these definitions, refer to the respective previous papers.

---

## 3. Definition of Class $E$

### 3.1 Information and Meta-information of Class $E$

We extend class $B$ ([3] §3.4) and define class $E$ with position information as follows.

**Information of class $E$:**

- $i$: scalar value (input value; identical to class $B$)
- $w$: real value (amplitude; for output; identical to class $B$)
- $x_0$: initial position (integer value)
- $x$: current position (real value)
- $dx$: displacement (integer value; position change per step)

**Meta-information:**

- $n$: counter (number of delays)
- $m$: number of delays corresponding to one wavelength
- $w_{\max}$: maximum amplitude (positive real value)
- color: attribute for identification

### 3.2 Operation Function $*$

$$
x = n \cdot dx + x_0
$$

$$
w = w_{\max} \cdot \sin\!\left(\frac{n}{m} \cdot 360°\right)
$$

Each time $F_D$ is passed, $n$ is incremented and the position $x$ and amplitude $w$ are updated.

---

## 4. Non-interacting Model (Passing Through)

### 4.1 Definition of Two Instances

We define two instances of class $E$.

**$E_1$ (red):**

| Information / Meta-information | Value |
| :----------------------------- | :---- |
| $x_0$                          | 0     |
| $x$                            | 0     |
| $dx$                           | $+1$  |
| $m$                            | 12    |
| $w_{\max}$                     | 2     |
| color                          | red   |

**$E_2$ (blue):**

| Information / Meta-information | Value |
| :----------------------------- | :---- |
| $x_0$                          | 24    |
| $x$                            | 24    |
| $dx$                           | $-1$  |
| $m$                            | 12    |
| $w_{\max}$                     | 1     |
| color                          | blue  |

### 4.2 Sequence of Passing Through

$E_1$ and $E_2$ execute the operation $*$ independently and do not affect each other. The position at each step is shown below.

|   $n$    | $x$ of $E_1$ | $x$ of $E_2$ |
| :------: | :----------: | :----------: |
|    0     |      0       |      24      |
|    1     |      1       |      23      |
|    2     |      2       |      22      |
| $\vdots$ |   $\vdots$   |   $\vdots$   |
|    12    |      12      |      12      |
| $\vdots$ |   $\vdots$   |   $\vdots$   |
|    24    |      24      |      0       |

At $n = 12$, $E_1$ and $E_2$ reach the same position $x = 12$, but since there is no interaction they pass through each other, and at $n = 24$ each reaches the initial position of the other.

### 4.3 Visualization of Passing Through

The spacetime diagram below shows the trajectories of $E_1$ (red) and $E_2$ (blue). The horizontal axis is position $x$ and the vertical axis is the number of delays $n$ (progressing from top to bottom). The size of each point corresponds to the amplitude $|w|$; dark color represents $w > 0$ and light color represents $w < 0$.

![Figure 9: Passing through in the non-interacting model](fig9_pass_through.svg)

**Figure 9: Non-interacting model — passing through of $E_1$ and $E_2$ ($m = 12$)**

At $n = 12$, both pass through the same position $x = 12$, but since there is no interaction they continue straight ahead without any effect on each other.

---

## 5. Model of Perfectly Elastic Collision

### 5.1 Definition of Class $F$

We extend class $E$ (§3) to define class $F$, which has collision functionality. The information and meta-information of class $F$ are identical to those of class $E$, with the following meta-information added:

- $r$: collision radius (a real number with $r > 0$; each instance may have a different value)

#### Collision Detection

For two instances $F_1$ (collision radius $r_1$) and $F_2$ (collision radius $r_2$), a collision is detected when the following condition is satisfied:

$$
(x_1 - x_2)^2 < \max(r_1^2,\; r_2^2)
$$

#### Update of $dx$ at Collision

When a collision occurs, $dx$ is updated using the one-dimensional perfectly elastic collision formula, treating $w_{\max}$ as an analogue of mass:

$$
dx_1^{\text{new}} = \frac{(w_{\max,1} - w_{\max,2}) \cdot dx_1 + 2\, w_{\max,2} \cdot dx_2}{w_{\max,1} + w_{\max,2}}
$$

$$
dx_2^{\text{new}} = \frac{(w_{\max,2} - w_{\max,1}) \cdot dx_2 + 2\, w_{\max,1} \cdot dx_1}{w_{\max,1} + w_{\max,2}}
$$

This update maintains the following two conserved quantities:

- **Analogue of momentum**: $Q = w_{\max,1} \cdot dx_1 + w_{\max,2} \cdot dx_2$
- **Analogue of energy**: $K = w_{\max,1} \cdot dx_1^2 + w_{\max,2} \cdot dx_2^2$

### 5.2 Sequence of Collision

Using the same initial conditions as §4.1, and setting $r_1 = r_2 = 1$ for simplicity.

**Before collision ($n < 12$):**

The trajectories are identical to those of $E_1$ and $E_2$ (see §4.2).

**Collision at $n = 12$:**

From $(x_1 - x_2)^2 < \max(r_1^2, r_2^2)$, a collision is detected. $dx$ is updated:

$$
dx_1^{\text{new}} = \frac{(2-1)(1) + 2(1)(-1)}{2+1} = \frac{-1}{3} = -\frac{1}{3}
$$

$$
dx_2^{\text{new}} = \frac{(1-2)(-1) + 2(2)(1)}{2+1} = \frac{5}{3}
$$

**Verification:**

$$
Q^{\text{new}} = 2 \cdot (-1/3) + 1 \cdot (5/3) = -2/3 + 5/3 = 1 = Q \quad \checkmark
$$

$$
K^{\text{new}} = 2 \cdot (1/9) + 1 \cdot (25/9) = 27/9 = 3 = K \quad \checkmark
$$

**After collision ($n > 12$):**

| $n$ | $x$ of $F_1$ | $x$ of $F_2$ |
| :---: | :---: | :---: |
| 12 | 12.00 | 12.00 |
| 13 | 11.67 | 13.67 |
| 14 | 11.33 | 15.33 |
| 15 | 11.00 | 17.00 |
| 18 | 10.00 | 22.00 |
| 21 | 9.00 | 27.00 |
| 24 | 8.00 | 32.00 |

The heavier $F_1$ with $w_{\max} = 2$ slowly returns in the reverse direction ($dx = -1/3$), while the lighter $F_2$ with $w_{\max} = 1$ is knocked away at high velocity ($dx = 5/3$).

### 5.3 Visualization of Collision

The spacetime diagram below shows the trajectories of $F_1$ (red) and $F_2$ (blue). At the collision point at $n = 12$ (dashed circle), both trajectories bend. The heavier $F_1$ slowly returns to the left after the collision, while the lighter $F_2$ is knocked away to the right at high velocity.

![Figure 10: Perfectly elastic collision model](fig10_elastic_collision.svg)

**Figure 10: Perfectly elastic collision model — collision at $n = 12$, $dx$ is updated**

Comparing with Figure 9 (non-interacting model), the trajectories before the collision ($n < 12$) are identical, but after the collision ($n > 12$) the slopes of the trajectories change to reflect the mass ratio (the ratio of $w_{\max}$).

---

## 6. Conclusion

We extended class $E$ (§3) to define instances with position information and displacement. For the case in which two instances propagate in opposite directions on the same space, we first constructed a non-interacting model (§4, passing through), and then realized a perfectly elastic collision by adding a collision radius $r$ and the perfectly elastic collision formula to the operation function $*$ via class $F$ (§5).

In either case, the definitions of the delay circuit $F_D$ and the NOT operation circuit are not changed. The presence or absence of interaction is expressed solely as a difference in the definition of the operation function $*$.

---

## References

[1] Noriaki Kihara, "Information-Theoretic Organization of Information Transmission," Research Note, April 2026.  
[2] Noriaki Kihara, "Organization of Symmetries Inherent in the Delay Circuit Model," Research Note, April 2026.  
[3] Noriaki Kihara, "Realization of Simple Harmonic Oscillation and Sinusoidal Waves in the Delay Circuit Model," Research Note, April 2026.
