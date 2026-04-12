# Realization of a Wave Packet Collapse Model in the Delay Circuit Model

**Author:** Noriaki Kihara
**Affiliation:** WF System Co., Ltd. (Osaka University, School of Engineering Science, graduated)
**Date:** April 2026
**Category:** Research Note
**DOI:** [10.5281/zenodo.19534355](https://doi.org/10.5281/zenodo.19534355)

---

## Abstract

This paper does not prove the physical mechanism of wave packet collapse. It extends the perfectly elastic collision framework defined in the previous paper [4] and realizes behavior analogous to wave packet collapse on the delay circuit model. By treating the wave packet radius $r$ and the maximum wave height $w_{\max}$ as information, and by defining an operation function $*$ that causes $r$ to contract and $w_{\max}$ to increase upon collision, we construct a process in which a spread-out wave becomes localized through interaction.

---

## 1. Introduction

In the previous paper [4], perfectly elastic collisions were realized by class $F$. In this paper, we extend class $F$ from [4] and define class $G$, in which the wave packet radius $r$ and the maximum wave height $w_{\max}$ are promoted to information.

What this paper demonstrates is the fact that, by adding to $*$ an operational rule under which $r$ contracts and $w_{\max}$ increases upon collision, behavior analogous to wave packet collapse can be constructed. This paper neither theoretically proves actual wave packet collapse nor contains any quantum-mechanical claims.

---

## 2. Definition of Class $G$

### 2.1 Information and Meta-information

In class $F$ of the previous paper [4], we promote $r$ (collision radius) and $w_{\max}$ (maximum wave height), which were meta-information, to information, and define class $G$.

**Information in class $G$:**

- $i$: scalar value (input value)
- $w$: real value (amplitude)
- $x_0$: initial position
- $x$: current position
- $dx$: displacement
- $r$: wave packet radius (real number with $r > 0$)
- $w_{\max}$: maximum wave height (positive real number)

**Meta-information:**

- $n$: counter (number of delays)
- $m$: wavelength (defined via $r$ by setting $m = r$)
- $F$: type (electron, photon, etc.)
- color: attribute for identification

### 2.2 Mass Analogue

The mass analogue is defined as $M = w_{\max} \cdot r$. $M$ is conserved before and after a collision.

### 2.3 Operation Function $*$

The updates to position and amplitude are identical to the previous paper [4]:

$$
x = n \cdot dx + x_0
$$

$$
w = w_{\max} \cdot \sin\!\left(\frac{n}{m} \cdot 360°\right)
$$

### 2.4 Collision Operations

#### Interaction Rules

- Photon and photon: no interaction (pass through each other)
- Photon and electron: interact

#### Wave Packet Collapse

When a collision occurs, the photon's $r$ is updated as follows:

$$
r^{\text{new}} = \min(r_1, r_2)
$$

Since $M = w_{\max} \cdot r$ is conserved, $w_{\max}$ increases as $r$ contracts:

$$
w_{\max}^{\text{new}} = \frac{M}{r^{\text{new}}} = \frac{w_{\max} \cdot r}{r^{\text{new}}}
$$

The wavelength is also updated: $m^{\text{new}} = r^{\text{new}}$

#### Update of $dx$

The update of $dx$ follows the perfectly elastic collision formula in [4] §5.1, using the mass analogue $M = w_{\max} \cdot r$:

$$
dx_1^{\text{new}} = \frac{(M_1 - M_2) \cdot dx_1 + 2\, M_2 \cdot dx_2}{M_1 + M_2}
$$

$$
dx_2^{\text{new}} = \frac{(M_2 - M_1) \cdot dx_2 + 2\, M_1 \cdot dx_1}{M_1 + M_2}
$$

---

## 3. Sequence of Wave Packet Collapse

### 3.1 Definition of Instances

Three instances are defined.

**$G_1$ (electron, red):**

| Information / Meta-information | Value |
|:---|:---|
| $x_0$ | 48 |
| $dx$ | 0 (at rest) |
| $r$ | 6 |
| $w_{\max}$ | 30 |
| $M = w_{\max} \cdot r$ | 180 |
| $F$ | electron |

**$G_2$ (photon, blue):**

| Information / Meta-information | Value |
|:---|:---|
| $x_0$ | 72 |
| $dx$ | $-8$ |
| $r$ | 48 |
| $w_{\max}$ | 1 |
| $M = w_{\max} \cdot r$ | 48 |
| $F$ | photon |

**$G_3$ (photon, green):**

| Information / Meta-information | Value |
|:---|:---|
| $x_0$ | 24 |
| $dx$ | $+8$ |
| $r$ | 48 |
| $w_{\max}$ | 1 |
| $M = w_{\max} \cdot r$ | 48 |
| $F$ | photon |

$G_2$ and $G_3$ are photons with a wide wave packet radius ($r = 48$), and $G_1$ is an electron with a narrow wave packet radius ($r = 6$). $G_1$ does not exist in the initial state; it is assumed to appear at $n = 3$ through external intervention ([2] §3.5).

### 3.2 Before Collision ($n = 0, 1, 2$)

$G_2$ and $G_3$ propagate independently. $G_1$ does not yet exist.

| $n$ | $x$ of $G_2$ | $x$ of $G_3$ |
| :---: | :---: | :---: |
| 0 | 72 | 24 |
| 1 | 64 | 32 |
| 2 | 56 | 40 |

Since $G_2$ and $G_3$ are both photons, they pass through each other without interacting.

### 3.3 $n = 3$: External Intervention and Collision

At $n = 3$, $G_1$ (electron) appears at $x = 48$ through external intervention.

The position of $G_2$ is $x = 48$, and the position of $G_3$ is $x = 48$, so all three reach the same position $x = 48$. A collision is determined from $(x_1 - x_2)^2 < \max(r_1^2, r_2^2)$.

**Wave packet collapse:**

$$
r^{\text{new}} = \min(48, 6) = 6
$$

$$
w_{\max}^{\text{new}} = \frac{48}{6} = 8
$$

$$
m^{\text{new}} = 6
$$

**Update of $dx$ (by symmetry):**

Since $G_2$ and $G_3$ collide with $G_1$ from symmetric positions and with symmetric velocities, the impulses on $G_1$ cancel each other, and $G_1$ remains at rest ($dx_1 = 0$).

$$
dx_2^{\text{new}} = \frac{(48 - 180)(-8) + 2 \cdot 180 \cdot 0}{48 + 180} = \frac{1056}{228} \approx +4.632
$$

$$
dx_3^{\text{new}} = \frac{(48 - 180)(+8) + 2 \cdot 180 \cdot 0}{48 + 180} = \frac{-1056}{228} \approx -4.632
$$

### 3.4 After Collision ($n > 3$)

| $n$ | $x$ of $G_1$ | $x$ of $G_2$ | $x$ of $G_3$ |
| :---: | :---: | :---: | :---: |
| 3 | 48.00 | 48.00 | 48.00 |
| 4 | 48.00 | 52.63 | 43.37 |
| 5 | 48.00 | 57.26 | 38.74 |
| 6 | 48.00 | 61.89 | 34.11 |
| 9 | 48.00 | 75.79 | 20.21 |
| 12 | 48.00 | 89.68 | 6.32 |

$G_1$ (electron) remains at rest, while $G_2$ and $G_3$ (photons) move away in opposite directions as collapsed wave packets ($r = 6$, $w_{\max} = 8$, $m = 6$).

### 3.5 Visualization of Wave Packet Collapse

Below, the shape of the wave packets at each delay count $n$ is shown. The horizontal axis is position $x$ and the vertical axis is amplitude $w$.

$n = 0, 1, 2$: $G_2$ (blue) and $G_3$ (green) are spread throughout space as wide wave packets ($r = 48$, $w_{\max} = 1$).

$n = 3$: $G_1$ (electron, red) appears through external intervention and a collision occurs. The wave packets of $G_2$ and $G_3$ collapse to $r = 6$, $w_{\max} = 8$.

$n = 4, 5, \ldots, 12$: The collapsed $G_2$ and $G_3$ move away in opposite directions as narrow, tall wave packets.

![Figure 11: Wave Packet Collapse](fig11_wave_packet_collapse.svg)

**Figure 11: Wave Packet Collapse — Shape change of wave packets from $n = 0$ to $n = 12$**

---

## 4. Conclusion

We extended class $F$ (perfectly elastic collisions) from the previous paper [4] and defined class $G$, in which the wave packet radius $r$ and the maximum wave height $w_{\max}$ are promoted to information. Under conservation of the mass analogue $M = w_{\max} \cdot r$, we constructed an operational rule in which $r$ contracts and $w_{\max}$ increases upon collision.

We realized a process in which photons with wide wave packets ($r = 48$, $w_{\max} = 1$) collide with an electron having a narrow wave packet ($r = 6$), causing the wave packets to collapse ($r = 6$, $w_{\max} = 8$) and the wavelength to shorten ($m = 6$). This is behavior analogous to wave packet collapse, but it contains no quantum-mechanical claims.

---

## References

[1] Noriaki Kihara, "Information-Theoretic Organization of Information Transmission," Research Note, April 2026.  
[2] Noriaki Kihara, "Organization of Symmetries Inherent in the Delay Circuit Model," Research Note, April 2026.  
[3] Noriaki Kihara, "Realization of Simple Harmonic Oscillation and Sinusoidal Waves in the Delay Circuit Model," Research Note, April 2026.  
[4] Noriaki Kihara, "Realization of Perfectly Elastic Collisions in the Delay Circuit Model," Research Note, April 2026.
