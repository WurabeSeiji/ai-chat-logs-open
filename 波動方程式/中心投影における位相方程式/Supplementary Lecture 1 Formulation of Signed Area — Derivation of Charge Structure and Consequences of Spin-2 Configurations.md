# Formulation of Signed Area — Derivation of Charge Structure and Consequences of Spin-2 Configurations

## — Three Consequences Derived from the Sign Product $P$ of [4] §7 —

**Author**: Noriaki Kihara (WF System Co., Ltd.)

**ORCID**: 0009-0006-7814-2793

**Date**: April 2026

**DOI**: https://doi.org/10.5281/zenodo.19731600

---

## §0 Purpose of This Paper

In the remark of [4] §7, it was stated that "the orientation structure of the signed area $M(\sigma)$ for square faces containing the $R$-axis, and the rigorous definition of conserved quantities based thereon, lie beyond the scope of this paper and will be addressed in a separate work."

This paper responds to that declaration and derives the following three results.

1. The rigorous definition of the signed area $M$ and its relation to the sign product $P$ of [4] Proposition 7.1 (§1)
2. The definition of the fermion area $M_F$ and its properties concerning antiparticles and color charge (§2)
3. Consequences of the $P = -1$ configuration in spin-2 $tR$-type (§3)

§1–§3 are self-contained within the combinatorial framework of [4]. §4 adopts the same stance as [4] §9 as an interpretive example, showing the correspondence between $M_F$ and charge quantization.

**Referenced papers**:
- [3] Full Spherical Coverage of Central Projection in Discrete Space and the Number-Theoretic Necessity of Five-Dimensional Background Space (W5)
- [4] Set Structure of the Six-Dimensional Hyperrectangle and Its Combinatorial Properties (W8)
- [5] Interactions of Shape-Invariant Waves (W11)

---

## §1 Definition of Signed Area

### 1.1 2-Faces

We use the set signs from [4] §1.3.3.

**Definition 1.1 (2-face)**. Among the six axes $(x, y, z, t, R, Q)$, a set in which exactly two axes $i, j$ have nonzero signs ($\sigma_i, \sigma_j \in \{+1, -1\}$) and all remaining axes have zero signs is called the **2-face** formed by axes $(i, j)$.

### 1.2 Signed Area

By [4] §1.1, each coordinate can take nonzero values $v_i, v_j \neq 0$ (either real-valued or integer-valued).

**Definition 1.2 (Signed area)**. For the 2-face formed by axes $(i, j)$, the **signed area** $M(i, j)$ is defined as:

$$M(i, j) = v_i \cdot v_j$$

$M$ is real-valued and carries two pieces of information: its **sign** (orientation) and its **absolute value** (magnitude of area):

$$\text{sign}(M) = \text{sign}(v_i) \cdot \text{sign}(v_j) \in \{+1, -1\}$$

$$|M| = |v_i| \cdot |v_j|$$

$\text{sign}(M) = +1$ is called **positive orientation** and $\text{sign}(M) = -1$ is called **reverse orientation**. The specific values of $|v_i|$ and $|v_j|$ are not constrained by this definition.

**Remark**. For cancellation of signed areas (e.g., particle-antiparticle pair annihilation) to satisfy $M + M' = 0$, the condition $|M| = |M'|$, i.e., equality of the absolute values of the axes, is required. The antiparticle definition in [5] §5.4 (only spatial axis signs are inverted; absolute values are invariant) automatically satisfies this condition.

**Proposition 1.2 (Signed area on integer lattices)**. By [3], coordinates that realize full spherical coverage of central projection in discrete space take values on the integer lattice $\mathbb{Z}^n$. Then for nonzero coordinates $v_i, v_j \in \mathbb{Z} \setminus \{0\}$:

(i) $M(i,j) = v_i \cdot v_j \in \mathbb{Z} \setminus \{0\}$

(ii) $\text{sign}(M) \in \{+1, -1\}$

*Proof*. (i) The integer ring $\mathbb{Z}$ is closed under multiplication ($a, b \in \mathbb{Z} \implies ab \in \mathbb{Z}$) and is an integral domain ($a \neq 0$ and $b \neq 0 \implies ab \neq 0$). Since $v_i, v_j \in \mathbb{Z} \setminus \{0\}$, we have $M = v_i \cdot v_j \in \mathbb{Z} \setminus \{0\}$. (ii) Since $M \neq 0$, $\text{sign}(M) \in \{+1, -1\}$ is well-defined. $\square$

**Remark**. Definition 1.2 is also applicable to the continuous case where coordinate values are real-valued, but the two-valuedness of $\text{sign}(M)$ is guaranteed solely by the condition $v_i, v_j \neq 0$. On integer lattices, this condition is automatically satisfied by the lattice structure.

### 1.3 Relation to the Sign Product $P$

**Proposition 1.1**. The sign product $P = (-1)^{|\{i \in \{t, R\} : \sigma_i < 0\}|}$ from [4] Proposition 7.1 coincides with the sign of $M(t, R)$ when both the $t$-axis and $R$-axis are nonzero: $P = \text{sign}(M(t, R))$.

*Proof*. When $\text{sign}(v_t), \text{sign}(v_R) \in \{+1, -1\}$:

| $\text{sign}(v_t)$ | $\text{sign}(v_R)$ | Number of negatives | $P$ | $\text{sign}(M)$ |
|:---:|:---:|:---:|:---:|:---:|
| $+$ | $+$ | 0 | $+1$ | $+1$ |
| $+$ | $-$ | 1 | $-1$ | $-1$ |
| $-$ | $+$ | 1 | $-1$ | $-1$ |
| $-$ | $-$ | 2 | $+1$ | $+1$ |

$P = \text{sign}(M)$ holds for all four configurations. $\square$

**Remark**. When only one of $t, R$ is nonzero, $P = (-1)^{0 \text{ or } 1}$, and the 2-face $M(t,R)$ is not defined (a single axis does not constitute a face). In this case, $P$ functions as a generalization of $\text{sign}(M)$.

---

## §2 Fermion Area

### 2.1 Definition

By [4] Definition 2.1, fermion types have exactly $k = 2$ of the four $xyzt$ axes nonzero.

**Definition 2.1 (Fermion area)**. The signed area $M_F = v_i \cdot v_j$ of the 2-face formed by the two nonzero axes $(i, j) \subseteq \{x, y, z, t\}$ of a fermion is called the **fermion area**. Its **orientation** is given by $\text{sign}(M_F) = \text{sign}(v_i) \cdot \text{sign}(v_j) \in \{+1, -1\}$.

### 2.2 Orientation of Fermion Area in [4] Table A

We compute $\text{sign}(M_F)$ for all fermions in [4] §9.9 Table A. The absolute values $|v_i|$ of each axis do not affect the results of this section.

**Charged leptons** (nonzero axes: 1 spatial axis + $t$-axis, $v_t > 0$, $Q = 4$, $c_1 = 1$):

| Particle | Set | Nonzero axes | $\text{sign}(M_F)$ |
|:---|:---|:---:|:---:|
| $e^-$ | $(+, 0, 0, +, 0, 4)$ | $x, t$ | $(+)(+)= +1$ |
| $\mu^-$ | $(0, +, 0, +, 0, 4)$ | $y, t$ | $(+)(+)= +1$ |
| $\tau^-$ | $(0, 0, +, +, 0, 4)$ | $z, t$ | $(+)(+)= +1$ |

**Neutrinos** (nonzero axes: 2 spatial axes, $t = 0$):

| Particle | Set | Nonzero axes | $\text{sign}(M_F)$ |
|:---|:---|:---:|:---:|
| $\nu_e$ | $(+, +, 0, 0, 0, 0)$ | $x, y$ | $(+)(+)= +1$ |
| $\nu_\mu$ | $(0, +, +, 0, 0, 0)$ | $y, z$ | $(+)(+)= +1$ |
| $\nu_\tau$ | $(+, 0, +, 0, 0, 0)$ | $x, z$ | $(+)(+)= +1$ |

**Up-type quarks** (nonzero axes: 1 spatial axis + $t$-axis, $v_t > 0$, $Q \in \{5,6,7\}$, $c_1 = 1$):

| Particle | Set | Nonzero axes | $\text{sign}(M_F)$ |
|:---|:---|:---:|:---:|
| $u$ | $(+, 0, 0, +, 0, Q)$ | $x, t$ | $(+)(+)= +1$ |
| $c$ | $(0, +, 0, +, 0, Q)$ | $y, t$ | $(+)(+)= +1$ |
| $t$ | $(0, 0, +, +, 0, Q)$ | $z, t$ | $(+)(+)= +1$ |

**Down-type quarks** (nonzero axes: 1 spatial axis + $t$-axis, $v_t < 0$, $Q \in \{1,2,3\}$, $c_1 = 0$):

| Particle | Set | Nonzero axes | $\text{sign}(M_F)$ |
|:---|:---|:---:|:---:|
| $d$ | $(+, 0, 0, -, 0, Q)$ | $x, t$ | $(+)(-)= -1$ |
| $s$ | $(0, +, 0, -, 0, Q)$ | $y, t$ | $(+)(-)= -1$ |
| $b$ | $(0, 0, +, -, 0, Q)$ | $z, t$ | $(+)(-)= -1$ |

**Summary**: It is solely the **orientation** (sign) of $M_F$ that determines the charge structure; it does not depend on the absolute value $|M_F| = |v_i| \cdot |v_j|$.

| Fermion type | Type of nonzero axes | $Q$ | $\text{sign}(M_F)$ |
|:---|:---|:---:|:---:|
| Charged lepton | (spatial, $t>0$) | 4 | $+1$ |
| Neutrino | (spatial, spatial) | 0 | $+1$ |
| Up-type quark | (spatial, $t>0$) | 5,6,7 | $+1$ |
| Down-type quark | (spatial, $t<0$) | 1,2,3 | $-1$ |

### 2.3 Fermion Area of Antiparticles

By [5] §5.4, for antiparticles the signs of the spatial axes $(x, y, z)$ are inverted, while the sign of the $t$-axis and the $Q$ value remain unchanged.

**Proposition 2.1**. When the fermion face consists of a (spatial, $t$) combination, the fermion area of the antiparticle has the opposite sign to that of the particle: $M_F^{\text{anti}} = -M_F$.

*Proof*. For the antiparticle, the spatial axis value is inverted to $-v_{\text{spatial}}$, while the $t$-axis remains unchanged. Therefore, $M_F^{\text{anti}} = (-v_{\text{spatial}}) \cdot v_t = -(v_{\text{spatial}} \cdot v_t) = -M_F$. $\square$

**Proposition 2.2**. When the fermion face consists of a (spatial, spatial) combination, the fermion area of the antiparticle is identical to that of the particle: $M_F^{\text{anti}} = M_F$.

*Proof*. For the antiparticle, both spatial axis values are inverted to $-v_{\text{spatial}_1}, -v_{\text{spatial}_2}$. Therefore, $M_F^{\text{anti}} = (-v_{\text{spatial}_1}) \cdot (-v_{\text{spatial}_2}) = v_{\text{spatial}_1} \cdot v_{\text{spatial}_2} = M_F$. $\square$

**Corollary 2.1**. Neutrinos in [4] Table A (spatial $\times$ spatial, $\text{sign}(M_F) = +1$) have antineutrinos that also possess $\text{sign}(M_F) = +1$. For charged leptons and all quarks (spatial $\times$ $t$), the antiparticle has $M_F$ with inverted sign ($\text{sign}(M_F^{\text{anti}}) = -\text{sign}(M_F)$).

### 2.4 Independence of Color Charge and Fermion Area

**Proposition 2.3**. The value of the $Q$-axis does not affect the fermion area $M_F$. Fermions with the same $\sigma_{\text{spatial}}, \sigma_t$ have identical $M_F$ regardless of the value of $Q$.

*Proof*. $\text{sign}(M_F) = \text{sign}(v_i) \cdot \text{sign}(v_j)$ is determined solely by the values of the $xyzt$ axes, and $Q$ does not appear in the definition. $\square$

---

## §3 The $P = -1$ Configuration of Spin-2 $tR$-Type

By [4] Proposition 7.1, spin-2 $tR$-type admits both $P = +1$ and $P = -1$ configurations. In [4] §9.7, the physical correspondence of $P = -1$ was left unresolved.

**Proposition 3.1**. The $P = -1$ configuration of spin-2 $tR$-type is equivalent to $\sigma_t$ and $\sigma_R$ having opposite signs.

*Proof*. By Proposition 1.1, $P = \text{sign}(M(t, R)) = \text{sign}(v_t) \cdot \text{sign}(v_R)$. $P = -1$ $\iff$ $\text{sign}(v_t) \cdot \text{sign}(v_R) = -1$ $\iff$ $v_t$ and $v_R$ have opposite signs. $\square$

**Proposition 3.2**. Applying the sign rule for attraction/repulsion from [5] §4.3 (same relative sign of $t$-axis components implies attraction; opposite signs imply repulsion) to the graviton interpretation of [4] §9.5:

- $P = +1$ ($t, R$ same sign): $t$-component has a single sign $\to$ attraction only
- $P = -1$ ($t, R$ opposite signs): $t$-component and $R$-component have opposite signs $\to$ repulsive along the $R$-axis direction

*Proof*. In [5] §4.3, the sign of interaction between direction-constrained waves is determined by $\text{sign}(k_t^{(1)}) \cdot \text{sign}(k_t^{(2)})$. When the spherical wave of spin-2 $tR$-type follows this sign rule, $P = +1$ implies both axes share the same directionality, yielding attraction only. For $P = -1$, the action along the $R$-direction has the opposite sign to the $t$-direction, acting repulsively with respect to the spatial scale. $\square$

**Remark**. The interpretation of the physical consequences of the $P = -1$ configuration is addressed in §4.3.

---

## §4 Interpretive Example: Derivation of Charge Structure

This section has the same status as [4] §9; it is an interpretive example that maps the combinatorial results of §1–§3 onto the Standard Model.

### 4.1 Third Component of Weak Isospin $I_3$

Using the $\text{sign}(M_F)$ introduced in §2.2 and the $c_1$ bit (weak isospin) from [4] §9.2, we define the third component of weak isospin $I_3$ of the Standard Model as follows.

For quarks ($c_2 c_3 \neq 00$), since $\text{sign}(M_F)$ directly distinguishes up-type ($+1$) from down-type ($-1$), we set $I_3 = \text{sign}(M_F)/2$.

For leptons ($c_2 c_3 = 00$), neutrinos ($Q = 0$, spatial $\times$ spatial) and charged leptons ($Q = 4$, spatial $\times$ $t$) both have $\text{sign}(M_F) = +1$ and cannot be distinguished by $\text{sign}(M_F)$ alone. We define $I_3 = (1 - 2c_1)/2$ using the $c_1$ bit.

| Fermion type | $Q$ | $c_1$ | $\text{sign}(M_F)$ | $I_3$ | SM $I_3$ |
|:---|:---:|:---:|:---:|:---:|:---:|
| Neutrino | $0$ | $0$ | $+1$ | $(1 - 2 \cdot 0)/2 = +1/2$ | $+1/2$ |
| Charged lepton | $4$ | $1$ | $+1$ | $(1 - 2 \cdot 1)/2 = -1/2$ | $-1/2$ |
| Up-type quark | $5,6,7$ | $1$ | $+1$ | $+1/2$ | $+1/2$ |
| Down-type quark | $1,2,3$ | $0$ | $-1$ | $-1/2$ | $-1/2$ |

The value of $I_3$ agrees with the Standard Model for all fermions. For quarks, $\text{sign}(M_F) = 2c_1 - 1$ holds, so $I_3 = \text{sign}(M_F)/2 = (2c_1 - 1)/2$ has the opposite sign relation to $(1 - 2c_1)/2$ for leptons.

### 4.2 Determination of Weak Hypercharge $Y$ from Color Charge Components $c_2 c_3$

The Standard Model charge formula (Gell-Mann–Nishijima relation):

$$Q_{\text{charge}} = I_3 + \frac{Y}{2}$$

Assume that $Y$ depends only on the color charge components $c_2 c_3$: $Y = Y_L$ when $c_2 c_3 = 00$ (leptons), and $Y = Y_Q$ when $c_2 c_3 \neq 00$ (quarks). Then the following two conditions uniquely determine $Y_L, Y_Q$.

**Condition (i): Integer charges for leptons**. For any fermion with $c_2 c_3 = 00$, $Q_{\text{charge}} = I_3 + Y_L/2 \in \mathbb{Z}$ must hold. Since $I_3 \in \{+1/2, -1/2\}$, $Y_L$ must be odd.

**Condition (ii): Integer charges for color-neutral baryons**. The total charge of a three-body color-neutral system ($c_2 c_3 = 01, 10, 11$, one each) with $c_2 c_3 \neq 00$, $Q_{\text{baryon}} = \sum_{a=1}^{3}(I_3^{(a)} + Y_Q / 2) \in \mathbb{Z}$, must hold for any combination of $I_3^{(a)} \in \{+1/2, -1/2\}$.

$$Q_{\text{baryon}} = \sum I_3^{(a)} + \frac{3 Y_Q}{2}$$

Since $\sum I_3^{(a)} \in \{-3/2, -1/2, +1/2, +3/2\}$ (always a half-integer), $3 Y_Q / 2$ must also be a half-integer. That is, $Y_Q = (2n+1)/3$ ($n \in \mathbb{Z}$).

**Determination of the minimal solution**. From condition (i), the odd number with smallest $|Y_L|$ is $Y_L = \pm 1$. Requiring negative charge for charged leptons:

$$Q_{\text{charge}}(e^-) = -\frac{1}{2} + \frac{Y_L}{2} < 0 \implies Y_L < 1$$

Hence $Y_L = -1$. From condition (ii), the solution with smallest $|Y_Q|$ is $n = 0$: $Y_Q = 1/3$.

**Proposition 4.1 (Unique determination of $Y$ from integer charge conditions)**. From the two-valuedness of $I_3 \in \{+1/2, -1/2\}$ and the condition that charges of observable particles (leptons and color-neutral baryons) are integers, the minimal solution for $Y$ is uniquely determined as:

$$Y_L = -1 \quad (c_2 c_3 = 00), \qquad Y_Q = +\frac{1}{3} \quad (c_2 c_3 \neq 00)$$

### 4.3 Consequences

**Consequence 1: Reproduction of the Gell-Mann–Nishijima formula**. Substituting the $I_3$ from §4.1 and $Y$ from Proposition 4.1:

| Fermion | $\text{sign}(M_F)$ | $I_3$ | $Y$ | $Q_{\text{charge}} = I_3 + Y/2$ | SM value |
|:---|:---:|:---:|:---:|:---:|:---:|
| $\nu$ | $+1$ | $+1/2$ | $-1$ | $0$ | $0$ ✓ |
| $e^-$ | $+1$ | $-1/2$ | $-1$ | $-1$ | $-1$ ✓ |
| $u$ | $+1$ | $+1/2$ | $+1/3$ | $+2/3$ | $+2/3$ ✓ |
| $d$ | $-1$ | $-1/2$ | $+1/3$ | $-1/3$ | $-1/3$ ✓ |

This is consistent with the charge map $f_Q^{\text{charge}}$ of [4] Definition 9.1.

**Consequence 2: Automatic satisfaction of anomaly cancellation**. The solution from Proposition 4.1 satisfies:

$$Y_L + 3 Y_Q = -1 + 3 \times \frac{1}{3} = 0$$

This corresponds to the gauge anomaly cancellation condition in the Standard Model. In this framework, it is derived algebraically from the structure of $c_2 c_3$ (one color-neutral state $c_2 c_3 = 00$ + three color charge states $c_2 c_3 \neq 00$) and the integer charge condition.

**Consequence 3: Charge equality of proton and electron**. The charge of the proton $(uud)$:

$$Q_p = \sum I_3^{(a)} + \frac{3 Y_Q}{2} = \left(\frac{1}{2} + \frac{1}{2} - \frac{1}{2}\right) + \frac{1}{2} = +1$$

The charge (absolute value) of the electron:

$$|Q_e| = \left|I_3(e^-) + \frac{Y_L}{2}\right| = \left|-\frac{1}{2} - \frac{1}{2}\right| = 1$$

$|Q_p| = |Q_e| = 1$; the reason that the proton and electron carry charges of equal magnitude is determined algebraically from the color charge structure of $c_2 c_3$ ($1 + 3$ classification) and the two-valuedness of $I_3$ ($\pm 1/2$). There is no need to introduce fractional charges $+2/3, -1/3$ as fundamental quantities.

**Consequence 4: Interpretation of spin-2 $P = -1$**. By Proposition 3.2 of §3, the graviton configuration with $P = -1$ acts repulsively along the $R$-axis direction. Since the $R$-axis in [4] §9.8 is the axis that links the subjective coordinate system to the spatial scale of the six-dimensional structure, the $P = -1$ configuration may correspond to **accelerating expansion of the spatial scale**. This is a candidate for the geometric description of the cosmological constant $\Lambda > 0$ (dark energy). However, whether the $P = -1$ configuration is physically realized (existence of selection rules) remains an open question, and this paper limits itself to stating the consequence.

---

## §5 Conclusion

In this paper, we provided the definition of the signed area $M(\sigma)$ declared in [4] §7 and derived the following three results.

1. The **signed area $M(i,j) = v_i \cdot v_j$** (Definition 1.2) coincides with the restriction of the sign product $P$ of [4] Proposition 7.1 to the case where both $t$ and $R$ are nonzero (Proposition 1.1).
2. The **fermion area $M_F$** (Definition 2.1) is computable for all fermions in [4] Table A, and we demonstrated sign inversion under antiparticles (Proposition 2.1), invariance of neutrino antiparticles (Proposition 2.2), and independence from the $Q$-axis (Proposition 2.3).
3. The **spin-2 $P = -1$** configuration corresponds to opposite signs of $v_t, v_R$, and from the sign rule of [5], it was derived to be repulsive along the $R$-axis direction (Propositions 3.1, 3.2).

In the interpretive example (§4):
- $I_3$ (quarks: $\text{sign}(M_F)/2$; leptons: $(1-2c_1)/2$) agrees with the third component of weak isospin in the Standard Model (§4.1)
- The integer charge condition uniquely determines $Y_L = -1, Y_Q = +1/3$ (Proposition 4.1)
- The anomaly cancellation condition $Y_L + 3Y_Q = 0$ is automatically satisfied (Consequence 2)
- The charge equality of the proton and electron is derived from the $1+3$ structure of $c_2 c_3$ (Consequence 3)

were demonstrated.

---

## References

- [3] N. Kihara, "Full Spherical Coverage of Central Projection in Discrete Space and the Number-Theoretic Necessity of Five-Dimensional Background Space," Zenodo preprint, DOI: 10.5281/zenodo.19624957 (2026).
- [4] N. Kihara, "Set Structure of the Six-Dimensional Hyperrectangle and Its Combinatorial Properties," Zenodo preprint, DOI: 10.5281/zenodo.19721125 (2026).
- [5] N. Kihara, "Interactions of Shape-Invariant Waves — Axial Displacement Transfer, Wave Packet Deformation, and Retroactive Construction of Causality," Zenodo preprint, DOI: 10.5281/zenodo.19721128 (2026).
