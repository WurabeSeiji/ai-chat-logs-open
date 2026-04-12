# Information Transmission Model via Delay Circuits: Summary of 5 Papers

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, Faculty of Engineering Science, graduated)  
**ORCID:** 0009-0004-6753-4020  
**Date:** April 2026

---

## Paper List

| # | Title | DOI |
|---|-------|-----|
| [1] | An Information-Theoretic Framework for Information Transmission | [10.5281/zenodo.19534345](https://doi.org/10.5281/zenodo.19534345) |
| [2] | Symmetries Inherent in the Delay Circuit Model | [10.5281/zenodo.19534347](https://doi.org/10.5281/zenodo.19534347) |
| [3] | Realization of Simple Harmonic Oscillation and Sinusoidal Waves in the Delay Circuit Model | [10.5281/zenodo.19534349](https://doi.org/10.5281/zenodo.19534349) |
| [4] | Realization of Perfectly Elastic Collision in the Delay Circuit Model | [10.5281/zenodo.19534353](https://doi.org/10.5281/zenodo.19534353) |
| [5] | Realization of a Wave Packet Collapse Model in the Delay Circuit Model | [10.5281/zenodo.19534355](https://doi.org/10.5281/zenodo.19534355) |

---

## Purpose of This Research

Starting from only two computational elements — the delay circuit $F_D$ (which transmits information with a one-step delay) and the NOT operation circuit (which inverts information instantaneously) — we deductively construct patterns ranging from oscillation, propagation, and collision to wave packet collapse. The starting point is the minimal definition that "the transmission of information $I$ is accompanied by a delay." No physical claims are made.

---

## Constructed Computational Propositions

### Paper [1]: Basic Structure of Information Transmission

**1. Three Properties of Information $I$**  
Three properties are defined for information $I$: existence of a null state, equivalence judgment, and inversion ($\overline{\overline{I}} = I$).  
(§2)

**2. Delay Circuit $F_D$ and NOT Operation Circuit**  
$F_D$ is defined as a circuit that outputs input information after a one-step delay. The NOT operation circuit is defined as a circuit that inverts input information without delay.  
(§3)

**3. Simple Harmonic Oscillation Pattern**  
By feedback connection of $F_D$ and NOT, the oscillation pattern $x(n) = I \cdot (-1)^{n+1}$ is constructed.  
(§4.1)

**4. Constant-Velocity Propagation Pattern**  
By serial connection of $F_D$ (open string model), the linear propagation pattern with arrival position $k = n$ is constructed.  
(§4.2)

**5. Fundamental String Vibration Pattern**  
By the closed string model ($F_D$ serial + NOT), a fundamental vibration pattern where one oscillation $= 2n$ operations is constructed.  
(§4.3)

---

### Paper [2]: Framework Extension

**6. Inversion Symmetry**  
The operation rule of $F_D$ is invariant under the transformation $I \to \overline{I}$.  
(§3.1)

**7. Superclass and Meta-Information**  
A superclass for grouping information of different classes, and meta-information (counter $n$, operation function $*$) attached to information are defined.  
(§3.2, §3.3)

**8. External Intervention**  
An operation that forcibly rewrites the internal information of a specific $F_D$ from outside the system is defined. From the perspective of the system's interior, this is observed as an event that suddenly appears from outside causality.  
(§3.5)

---

### Paper [3]: Four Types of Wave Patterns

**9. Class $A$ (Square Wave)**  
A square wave is constructed as a class without an operation function $*$.  
(§3.3)

**10. Class $B$ (Transverse Wave / Sinusoidal Wave)**  
By setting $w = w_{\max} \cdot \sin(n/m \cdot 360°)$ as the operation function $*$, a transverse wave with displacement perpendicular to the propagation direction is constructed.  
(§3.4)

**11. Class $C$ (Longitudinal Wave)**  
Using the same $*$ as Class $B$, a longitudinal wave is constructed by interpreting the output $w$ as displacement in the propagation direction.  
(§3.5)

**12. Class $D$ (Rotational Wave)**  
By setting $w = (n/m \cdot 360°) \bmod 360°$ as $*$, a rotational wave is constructed.  
(§3.6)

The definitions of $F_D$ and the NOT operation circuit are not modified at all. Differences in wave types are expressed solely through the choice of class definition and $*$.

---

### Paper [4]: Perfectly Elastic Collision

**13. Class $E$ (Transverse Wave with Position Information)**  
Class $B$ is extended to define Class $E$ with initial position $x_0$, current position $x$, and displacement $dx$.  
(§3)

**14. Non-Interaction Model (Pass-Through)**  
A model is constructed where two instances of Class $E$ propagate in opposite directions and pass through each other without interaction at the same position.  
(§4)

**15. Class $F$ (with Collision Functionality)**  
Class $F$ is defined by adding collision radius $r$ and collision condition $(x_1 - x_2)^2 < \max(r_1^2, r_2^2)$ to Class $E$. Using $w_{\max}$ as an analogue of mass, $dx$ is updated according to the formula for one-dimensional perfectly elastic collision.  
(§5)

**16. Conservation of Momentum and Energy**  
It is verified that the momentum analogue $Q = w_{\max,1} \cdot dx_1 + w_{\max,2} \cdot dx_2$ and the energy analogue $K = w_{\max,1} \cdot dx_1^2 + w_{\max,2} \cdot dx_2^2$ are conserved before and after collision.  
(§5.2)

---

### Paper [5]: Wave Packet Collapse

**17. Class $G$ (Promotion of Wave Packet Radius and Maximum Wave Height to Information)**  
Class $G$ is defined by promoting the meta-information $r$ (wave packet radius) and $w_{\max}$ (maximum wave height) of Class $F$ to information.  
(§2)

**18. Mass Analogue $M = w_{\max} \cdot r$**  
The mass analogue is defined as $M = w_{\max} \cdot r$, and it is shown that $M$ is conserved before and after collision.  
(§2.2)

**19. Wave Packet Collapse Rule**  
A rule is constructed where the wave packet radius collapses as $r^{\text{new}} = \min(r_1, r_2)$ upon collision, and the maximum wave height increases as $w_{\max}^{\text{new}} = M / r^{\text{new}}$ due to conservation of $M$.  
(§2.4)

**20. Wave Packet Collapse Sequence**  
A process is realized where a photon with a wide wave packet ($r = 48$, $w_{\max} = 1$) collides with an electron with a narrow wave packet ($r = 6$), causing the wave packet to collapse to $r = 6$, $w_{\max} = 8$.  
(§3)

---

## Framework Structure

```
Definition of Information I (null state, equivalence, inversion)
│
├── Basic Elements (Paper [1])
│    ├── Delay Circuit F_D (transmits with one-step delay)
│    ├── NOT Operation Circuit (inverts instantaneously)
│    └── Three Basic Patterns (oscillation, propagation, string vibration)
│
├── Framework Extension (Paper [2])
│    ├── Inversion Symmetry
│    ├── Superclass, Meta-Information, Operation Function *
│    └── External Intervention
│
└── Applications (Papers [3]–[5])
     ├── Four Types of Wave Patterns (square, transverse, longitudinal, rotational)
     ├── Perfectly Elastic Collision (momentum and energy conservation)
     └── Wave Packet Collapse (conservation of M = w_max · r)
```

This framework constructs a series of patterns from oscillation to wave packet collapse by changing only the class definition and the choice of operation function $*$, without any modification to the definitions of $F_D$ and the NOT operation circuit. Each pattern takes a form mathematically identical to known physical formulas, but this is a consequence of the consistency of the computational model, not a physical claim.

---

## What This Research Does Not Claim

This research is limited to describing the consistency of the computational model. The following claims are not made:

- That delay circuits represent physical entities or devices
- That the derived patterns constitute derivations or proofs of physical laws
- That the physical mechanism of wave packet collapse is explained
- That the mass analogue $w_{\max}$ is identical to physical mass
- That external intervention corresponds to quantum mechanical observation

These physical identifications require additional assumptions beyond the computational propositions of this research (Paper [5] §4).

---

## References

![QR Paper 1](qr_paper1.png)
[1] Kihara, N. (2026). *An Information-Theoretic Framework for Information Transmission*. DOI: 10.5281/zenodo.19534345.

&nbsp;

![QR Paper 2](qr_paper2.png)
[2] Kihara, N. (2026). *Symmetries Inherent in the Delay Circuit Model*. DOI: 10.5281/zenodo.19534347.

&nbsp;

![QR Paper 3](qr_paper3.png)
[3] Kihara, N. (2026). *Realization of Simple Harmonic Oscillation and Sinusoidal Waves in the Delay Circuit Model*. DOI: 10.5281/zenodo.19534349.

&nbsp;

![QR Paper 4](qr_paper4.png)
[4] Kihara, N. (2026). *Realization of Perfectly Elastic Collision in the Delay Circuit Model*. DOI: 10.5281/zenodo.19534353.

&nbsp;

![QR Paper 5](qr_paper5.png)
[5] Kihara, N. (2026). *Realization of a Wave Packet Collapse Model in the Delay Circuit Model*. DOI: 10.5281/zenodo.19534355.
