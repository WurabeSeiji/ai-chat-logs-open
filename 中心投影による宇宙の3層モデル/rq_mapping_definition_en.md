# Construction of $(R, Q)$ Mapping Adopting a 6-Dimensional Hypercube

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v1.0  
**Date**: April 2026  
**DOI**: [10.5281/zenodo.19692853](https://doi.org/10.5281/zenodo.19692853)

---

## Abstract

This paper defines an interface connecting the subjective spaces of the three-layer structure model of central projection [M1] with the combinatorial structure of a 6-dimensional hypercube [W8].

This paper adopts two premises and, under those premises, defines terminology and notation. It contains no theorem proofs. The purpose of this paper is to provide a common definitional foundation that subsequent papers may reference.

**Keywords**: 6-dimensional hypercube, central projection, subjective space, metadata, $(R, Q)$ mapping

---

## §1 Adopted Premises

This paper adopts the following two premises. These are not objects of proof in this paper but are axiomatically adopted as the foundation for the definitions constructed herein.

### Premise I (Adoption of the 6-Dimensional Hypercube)

We adopt the 6-dimensional hypercube defined in [W8] (a combinatorial structure with axes $x, y, z, t, R, Q$). Here, $x, y, z, t, R$ are five axes taking real values, and $Q \in \{0, 1, 2, 3, 4, 5, 6, 7\}$ is an 8-valued discrete label axis with 3-bit encoding $Q = 4R + 2G + B$ ($R, G, B \in \{0, 1\}$) ([W8] Definitions 1.1, 1.2).

### Premise II (Metadata-Holding Capability)

A subjective space can hold arbitrary metadata that is numerically representable.

**Remark**: Premise I is an external assumption of this paper and does not claim that the 6-dimensional hypercube is derived from the geometry of central projection (see [W8] §1). Premise II does not restrict in advance the types of metadata that can be held.

---

## §2 Definitions

### Definition 2.1 (Schema and Instance)

The combinatorial structure of the 6-dimensional hypercube adopted in Premise I is called the **schema** $\mathcal{S}$. For each subjective space $S^n(R_1)$ ($R_1 > 0$), the concrete entity obtained by assigning $\mathcal{S}$ is called an **instance**.

### Definition 2.2 ($(R, Q)$ Mapping)

The map that identifies instances is called the $(R, Q)$ mapping and is defined as follows:

$$
\varphi : \mathcal{I} \to \mathbb{R}_{> 0} \times \{0, 1, 2, 3, 4, 5, 6, 7\}
$$

Here, $\mathcal{I}$ is the set of instances, and $\varphi$ assigns to each instance a pair consisting of a curvature radius $R_1 \in \mathbb{R}_{> 0}$ and a discrete label $Q \in \{0, 1, \ldots, 7\}$.

**Remark**: The $R$ component is identified with the curvature radius parameter $R_1$ of [M1] §3.1. The $Q$ component is the value of the label axis of [W8], held as metadata by virtue of Premise II.

### Definition 2.3 (Uniqueness of Schema and Multiplicity of Instances)

The schema $\mathcal{S}$ is unique. No restriction is imposed on the cardinality of the instance set $\mathcal{I}$. It is not excluded that distinct instances may share the same $(R_1, Q)$ value.

**Remark**: The discrimination of multiple instances sharing the same $(R_1, Q)$ is outside the definitional scope of this paper and is deferred to subsequent papers.

---

## §3 Relationship with Prior Works

### §3.1 Relationship with the Three-Layer Structure Model [M1]

The instances in Definition 2.1 of this paper are a concretization of the structure permitted by Condition 3 (non-constraint on multiplicity) of [M1] §2. [M1] sets as a condition that infinitely many subjective spaces may coexist at different curvature radii, and this paper constructs its definitions within the scope of this condition.

### §3.2 Relationship with the 6-Dimensional Hypercube [W8]

This paper adopts the combinatorial structure of [W8] as the schema $\mathcal{S}$, but the 62-state enumeration results and the interpretation examples of §9 of [W8] are not included in the definitions of this paper. The internal structure of [W8] is referenced independently as the content of the schema.

---

## §4 What This Paper Does Not Claim

1. That the 6-dimensional hypercube is isomorphic to the subjective spaces of central projection.
2. That the 6-dimensional hypercube is derived from the geometry of central projection.
3. That the definitions of this paper cannot be constructed with structures other than the 6-dimensional hypercube.
4. The physical existence of instances.
5. Physical interpretations of the $(R, Q)$ mapping (such as the micro black hole interpretation).
6. Methods for discriminating multiple instances sharing the same $(R_1, Q)$.

---

## References

- [P1]: Noriaki Kihara, Geometric Formulation of Four-Dimensional Spacetime via Central Projection, Zenodo, 2026. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)
- [M1]: Noriaki Kihara, Three-Layer Structure Model of Central Projection: Nested Structure of Subjective Spaces via $R$–$R_1$–$R_0$ and Middleware Geometry, Zenodo, 2026. DOI: [10.5281/zenodo.19691713](https://doi.org/10.5281/zenodo.19691713)
- [M2]: Noriaki Kihara, Regularity of Central Projection Between Hyperspherical Shells: A Formulation Without Tangent Hyperplanes and Degeneracy at $R = 0$, Zenodo, 2026. DOI: [10.5281/zenodo.19692192](https://doi.org/10.5281/zenodo.19692192)
- [W8]: Noriaki Kihara, Set Structure and Combinatorial Properties of a 6-Dimensional Hypercube, Zenodo, 2026. DOI: [10.5281/zenodo.19688521](https://doi.org/10.5281/zenodo.19688521)

---

*v1.0*
