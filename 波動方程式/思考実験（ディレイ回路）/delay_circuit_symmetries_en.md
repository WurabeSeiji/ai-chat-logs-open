# Symmetries Inherent in the Delay Circuit Model

**Author:** Noriaki Kihara
**Affiliation:** WF System Co., Ltd. (Osaka University, School of Engineering Science, graduated)
**Date:** April 2026
**Category:** Research Note
**DOI:** [10.5281/zenodo.19534347](https://doi.org/10.5281/zenodo.19534347)

---

## Abstract



---

## 1. Introduction



---

## 2. Preliminaries: Summary of Definitions from the Previous Paper



---

## 3. Organization of Symmetries

### 3.1 Symmetry of $I$ and $\overline{I}$

In the previous paper [1] §2.2 Property 5, the inversion $\overline{I}$ was defined as an operation within the same class. That is, $I$ and $\overline{I}$ belong to the same class.

The operation rule of the delay circuit $F_D$ depends solely on the equivalence check between the input information and the internal information (whether $\mathrm{IN} \neq I$ holds). Replacing $I$ with $\overline{I}$ does not change the equivalence check or the structure of the operation rule.

**Definition (Inversion Symmetry):** The operation rule of the delay circuit $F_D$ is invariant under the transformation $I \to \overline{I}$. This is called inversion symmetry.

This property is a trivial consequence that follows directly from the definition ([1] §2.2) that inversion is an operation within the same class and satisfies $\overline{\overline{I}} = I$.

### 3.2 Definition of the Superclass

The previous paper [1] dealt only with operations within the same class (the delay circuit $F_D$ and the NOT operation circuit). In this section, we define a concept for handling information belonging to different classes in a unified manner.

A tuple $(A, B, C_1, C_2, \ldots, C_m)$ of information belonging to different classes is defined as a **superclass** $\mathcal{I}$. The number of child elements is arbitrary (including zero), and neither the class nor the value of each child element is constrained. When there are zero child elements, this corresponds to $\mathrm{null}$.

The superclass is required to satisfy the following properties.

#### Property 1: Order Independence

The superclass does not depend on the order in which child elements appear.

$$
(A, B) = (B, A)
$$

#### Property 2: Transparency of null

A $\mathrm{null}$ element as a child element of the superclass is equivalent to its absence.

$$
(A, \mathrm{null}, B) = (A, B)
$$

#### Property 3: Existence of Inverse Operations

For any operation on a superclass $\mathcal{I}$, an inverse operation exists.

The operation satisfying Properties 1 through 3 above is denoted $*$.

#### Definition of Meta-Information

The base class of the superclass has **meta-information** that is distinct from the information $I$ itself. Meta-information consists of the following two components:

- **Counter $n$:** Accumulates the number of passes through the delay circuit. Its initial value is $\mathrm{null}$.
- **Operation function $*$:** The operation executed when passing through the delay circuit. Its initial value is $\mathrm{null}$. The operation function $*$ can reference $n$.

Meta-information accompanies information $I$ but is not the value of $I$ itself. Therefore, the definition of information $I$ in [1] §2.2 remains unchanged.

$n$ is an exception to the inverse operation rule: even when an inverse operation is performed, $n$ is preserved.

### 3.3 Extended Definition of $F_D$

The operation rule of $F_D$ defined in [1] §2.2 (copying input information to internal information and passing it to the output) is not changed. When information of the superclass passes through, $F_D$ performs the following operations on the meta-information of the passing information, in addition to the original operation rule:

1. Increment the counter $n$ by $+1$ (if $n = \mathrm{null}$, set it from $0$ to $1$)
2. Execute the operation function $*$ if it is not $\mathrm{null}$

The entity that operates on meta-information is $F_D$; the information $I$ itself does not need to recognize that it "has passed through $F_D$."

The NOT operation circuit does not operate on meta-information. This structurally corroborates the fact stated in [1] §2.3 that "because the NOT operation circuit's operation involves no delay, it is not counted in the delay count $n$."

When an operation generates child classes on its own, $n$ can be re-initialized.

### 3.4 Consistency between the Superclass and the Delay Circuit

The superclass $\mathcal{I}$ satisfies the definition of information $I$ in [1] §2.2 ("any information that can be represented numerically"). That is, a tuple of information $(A, B)$ or $(C_1, C_2, \ldots, C_m)$ belonging to the superclass can each be treated as a single piece of information $I$.

Therefore, the definitions of the delay circuit $F_D$ and the NOT operation circuit ([1] §2.2–§2.4) apply to superclass information without modification. The meta-information operations added in §3.3 are not a modification of the original operation rule of $F_D$ but rather the addition of a new operational layer.

### 3.5 External Intervention from Outside the System

The operation rule of $F_D$ defined in the previous paper [1] and the meta-information operations in §3.3 specify the deterministic behavior inside the system. In this section, we define intervention from outside the system.

**External Intervention:** We allow the operation of forcibly overwriting the internal information of a specific $F_D$ from outside the system. This operation is performed independently of the operation rules inside the system.

External intervention is an exception to Property 3 of §3.2: information that has undergone external intervention is excluded from the scope of inverse operations. That is, external intervention introduces irreversibility into the interior of the system.

From the perspective of an observer inside the system, external intervention is an event that appears suddenly from outside the causal order.

Note: This definition suggests that the present model admits a hierarchical structure. The relationship between an upper-level system and a lower-level system is outside the scope of this paper.

---

## 4. Discussion



---

## 5. Conclusion

Formulating this model as equations and mapping it to a physical framework are outside the scope of this paper.

---

## References

[1] Noriaki Kihara, "Information-Theoretic Organization of Information Transmission," Research Note, April 2026.
