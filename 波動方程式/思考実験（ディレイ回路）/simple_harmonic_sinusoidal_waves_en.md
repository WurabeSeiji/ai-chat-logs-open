# Realization of Simple Harmonic Oscillation and Sinusoidal Waves in the Delay Circuit Model

**Author:** Noriaki Kihara
**Affiliation:** WF System Co., Ltd. (Osaka University, School of Engineering Science, graduated)
**Date:** April 2026
**Category:** Research Note
**DOI:** [10.5281/zenodo.19534349](https://doi.org/10.5281/zenodo.19534349)

---

## Abstract

The delay circuit model constructed in the previous paper [1] could only express square wave oscillation, because it used only the class $A = \{I, \overline{I}\}$ and the inversion (NOT) operation. In this paper, using the framework of superclasses, meta-information, and the operation function $*$ introduced in the previous paper [2], we realize sinusoidal waves (transverse waves), longitudinal waves, and rotational waves by concretely defining the class information and the $*$ operation. In all cases, we show that different wave patterns can be obtained solely by the choice of class definition and operation function $*$, without modifying the definitions of the delay circuit $F_D$ or the NOT operation circuit.

---

## 1. Introduction

In the previous paper [1], we constructed discrete square wave oscillation using only the class $A = \{I, \overline{I}\}$ and the inversion (NOT) operation. In the previous paper [2], we introduced superclasses, meta-information (counter $n$, operation function $*$), and an extended definition of $F_D$.

In this paper, using the framework of [2], we construct models that realize sinusoidal waves (longitudinal waves, transverse waves, and rotational waves) by concretely defining the class $A$ information, meta-information, and operation function $*$.

---

## 2. Preparation: Summary of Definitions from Previous Papers

This paper presupposes the delay circuit $F_D$, the NOT operation circuit, and the properties of information $I$ (null state, equivalence judgment, inversion) defined in the previous paper [1]. We also use the following concepts introduced in the previous paper [2]:

- **Superclass**: A tuple $(A, B, C_1, \ldots, C_m)$ of information belonging to different classes ([2] §3.2)
- **Meta-information**: Counter $n$ (number of delays) and operation function $*$ ([2] §3.2)
- **Extended definition of $F_D$**: When passing through $F_D$, the counter $n$ is incremented, and $*$ is executed if it is not $\mathrm{null}$ ([2] §3.3)

For the details of these definitions, refer to the respective previous papers.

---

## 3. Definition of the Simple Harmonic Oscillation Class $A$

### 3.1 Simple Harmonic Oscillation Circuit in the Previous Paper

In the previous paper [1] §3, a feedback structure was constructed by connecting the delay circuit $F_D$ and the NOT operation circuit: the output of $F_D$ is connected to the input of NOT, and the output of NOT is connected to the input of $F_D$.

![Figure 4: Simple harmonic oscillation circuit](fig4_oscillation_circuit.svg)

**Figure 1: Simple harmonic oscillation circuit (reprinted from [1] Figure 4)**

The output transitions of the simple harmonic oscillation circuit at each step are shown below.

| Number of delays $n$ | OUT of $F_D$ | OUT of NOT |
| :----------: | :-------------: | :-------------: |
| 0 | $\mathrm{null}$ | $\mathrm{null}$ |
| 1 | $I$ | $\overline{I}$ |
| 2 | $\overline{I}$ | $I$ |
| 3 | $I$ | $\overline{I}$ |
| 4 | $\overline{I}$ | $I$ |
| $\vdots$ | $\vdots$ | $\vdots$ |

For $n \geq 1$, the output of $F_D$ alternates between $I$ and $\overline{I}$. This is a discrete oscillation.

Taking information $I$ as a scalar value and defining inversion as $\overline{I} = -I$, the output of $F_D$ in the simple harmonic oscillation circuit is expressed as follows:

$$
x(n) = I \cdot (-1)^{n+1} \quad (n = 1, 2, 3, \ldots)
$$

That is, the output alternates between $I$ and $-I$ with respect to the number of delays $n$. This is a square wave.

### 3.2 Definition of the Square Wave Class $A$

In the framework of the previous paper [2], we redefine the square wave oscillation of §3.1.

- **Class $A$ information**: Scalar value $i$
- **Meta-information**: $\mathrm{null}$ (neither the counter $n$ nor the operation function $*$ is required)

Inversion is given as $\overline{i} = -i$ by the NOT operation circuit. Since the meta-information is $\mathrm{null}$, no additional operation is performed when passing through $F_D$, and the same result as the simple harmonic oscillation circuit of §3.1 is obtained.

$$
x(n) = i \cdot (-1)^{n+1} \quad (n = 1, 2, 3, \ldots)
$$

### 3.3 Propagation and Oscillation via the Square Wave Class $A$

We apply the square wave class $A$ defined in §3.2 (scalar value $i$, meta-information $\mathrm{null}$) to the open rope model of the previous paper [1] §4 and the closed rope model of [1] §5.

Since the meta-information is $\mathrm{null}$, no additional operation is performed when passing through $F_D$. Only the basic operation rule of $F_D$ (copying the input information internally and outputting it after a delay of one step) operates.

Therefore, in the open rope model, information $i$ passes through each $F_D$ one step at a time unchanged, giving the same propagation sequence as [1] §4 (arrival position $k = n$). In the closed rope model, including inversion by the NOT operation circuit, the same oscillation sequence as [1] §5 is obtained (1 oscillation $= 2n$ operations).

That is, the square wave class $A$ exactly reproduces the results of all configurations in the previous paper [1] (simple harmonic oscillation circuit, open rope, closed rope).

### 3.4 Definition of the Transverse Wave Class $B$

We define class $B$ to realize a transverse wave as follows.

**Class $B$ information:**

- $i$: Scalar value (input value; same as class $A$)
- $w$: Real value (amplitude; for output)

**Meta-information:**

- $n$: Counter (number of delays; as defined in [2] §3.3)
- $m$: Number of delays corresponding to the wavelength (a positive integer $1, 2, 3, \ldots$)

**Operation function $*$:**

$$
w = i \cdot \sin\!\left(\frac{n}{m} \cdot 360°\right)
$$

Each time a pass through $F_D$ occurs, $n$ is incremented and the operation function $*$ is executed to update $w$. By observing the output $w$, a sinusoidal oscillation is obtained as the number of delays $n$ advances.

### 3.5 Oscillation Sequence of Class $B$

For $i = 1$ and $m = 12$, the transitions of the output $w$ of $F_D$ in the simple harmonic oscillation circuit are shown below.

| Number of delays $n$ | $n/m \cdot 360°$ | $w = \sin(n/m \cdot 360°)$ |
| :----------: | :---------------: | :------------------------: |
| 0 | 0° | 0.000 |
| 1 | 30° | 0.500 |
| 2 | 60° | 0.866 |
| 3 | 90° | 1.000 |
| 4 | 120° | 0.866 |
| 5 | 150° | 0.500 |
| 6 | 180° | 0.000 |
| 7 | 210° | −0.500 |
| 8 | 240° | −0.866 |
| 9 | 270° | −1.000 |
| 10 | 300° | −0.866 |
| 11 | 330° | −0.500 |
| 12 | 360° | 0.000 |

One period is completed at $n = 12$, and the same pattern repeats thereafter.

![Figure 7: Oscillation sequence of class B](fig7_sine_wave_sequence.svg)

**Figure 7: Oscillation sequence of class $B$ (meta-information $m = 12$, $i = 1$)**

The bar chart shows the discrete output $w$ at each number of delays $n$, and the thin line shows the continuous sinusoidal wave.

## 4. Realization of Longitudinal Waves

### 4.1 Definition of the Longitudinal Wave Class $C$

We define class $C$ to realize a longitudinal wave. The structure of the information, meta-information, and operation function $*$ of class $C$ is exactly the same as that of class $B$ (§3.4).

- **Class $C$ information**: Scalar value $i$, real value $w$
- **Meta-information**: Counter $n$, number of delays $m$ corresponding to the wavelength
- **Operation function $*$**: $w = i \cdot \sin(n / m \cdot 360°)$

The only difference from class $B$ is in the interpretation of the output $w$. In class $B$, $w$ is treated as a displacement perpendicular to the direction of propagation (transverse wave), whereas in class $C$, $w$ is treated as a **displacement in the direction of propagation**. That is, $w$ is interpreted as the deviation of each $F_D$ from its reference position along the horizontal axis.

### 4.2 Visualization of Longitudinal Waves

When class $C$ is applied to the open rope model, each $F_D$ is displaced by $w$ from its reference position in the direction of propagation. As a result, regions where the bars are dense (compression) and regions where they are sparse (rarefaction) appear alternately, yielding the compression-rarefaction pattern of a longitudinal wave.

![Figure 8: Compression-rarefaction pattern of a longitudinal wave](fig8_longitudinal_wave.svg)

**Figure 8: Compression-rarefaction pattern of a longitudinal wave by class $C$ (meta-information $m = 12$, $i = 1$)**

Dashed lines indicate the reference positions of each $F_D$; solid lines indicate the positions after displacement by $w$.

---

## 5. Realization of Rotational Waves

### 5.1 Definition of the Rotational Wave Class $D$

We define class $D$ to realize a rotational wave. The structure of the information and meta-information of class $D$ is exactly the same as that of class $B$ (§3.4).

- **Class $D$ information**: Scalar value $i$, real value $w$
- **Meta-information**: Counter $n$, number of delays $m$ corresponding to the wavelength

The only differences from class $B$ are the operation function $*$ and the interpretation of the output $w$.

**Operation function $*$:**

$$
w = \frac{n}{m} \cdot 360° \mod 360°
$$

The output $w$ is interpreted as an **angle**. Each time a pass through $F_D$ occurs, $n$ is incremented and $w$ is updated. $w$ advances from $0°$ to $360°$ in constant increments, and returns to $0°$ upon reaching $360°$.

### 5.2 Oscillation Sequence of Class $D$

For $i = 1$ and $m = 12$, the transitions of the output $w$ (angle) of $F_D$ in the simple harmonic oscillation circuit are shown below.

| Number of delays $n$ | $w = (n/m \cdot 360°) \mod 360°$ |
| :----------: | :------------------------------: |
| 0 | 0° |
| 1 | 30° |
| 2 | 60° |
| 3 | 90° |
| 4 | 120° |
| 5 | 150° |
| 6 | 180° |
| 7 | 210° |
| 8 | 240° |
| 9 | 270° |
| 10 | 300° |
| 11 | 330° |
| 12 | 0° |
| 13 | 30° |
| 14 | 60° |
| $\vdots$ | $\vdots$ |

One full rotation ($360°$) is completed at $n = 12$, and the same pattern repeats thereafter. If the output $w$ is interpreted as an angle, this represents motion rotating at a constant angular velocity.

---

## 6. Conclusion

Using the framework of superclasses, meta-information, and operation function $*$ introduced in the previous paper [2], we defined classes $B$, $C$, and $D$ to realize transverse waves (sinusoidal waves), longitudinal waves, and rotational waves. The correspondence of each class is summarized below.

| Class | Information | Meta-information | Operation $*$ | Interpretation of output $w$ | Wave realized |
|:---:|:---|:---|:---|:---|:---|
| $A$ | $i$ | $\mathrm{null}$ | $\mathrm{null}$ | — | Square wave |
| $B$ | $i$, $w$ | $n$, $m$ | $w = i \cdot \sin(n/m \cdot 360°)$ | Displacement perpendicular to propagation direction | Transverse wave (sinusoidal wave) |
| $C$ | $i$, $w$ | $n$, $m$ | $w = i \cdot \sin(n/m \cdot 360°)$ | Displacement in propagation direction | Longitudinal wave |
| $D$ | $i$, $w$ | $n$, $m$ | $w = (n/m \cdot 360°) \mod 360°$ | Angle | Rotational wave |

In all cases, the definitions of the delay circuit $F_D$ and the NOT operation circuit ([1] §2.2–§2.4) have not been modified in any way. We have shown that four types of wave patterns, including square waves, can be constructed solely by the choice of class information and operation function $*$.



---

## References

[1] Noriaki Kihara, "Information-Theoretic Organization of Information Transmission," Research Note, April 2026.  
[2] Noriaki Kihara, "Organization of Symmetries Inherent in the Delay Circuit Model," Research Note, April 2026.
