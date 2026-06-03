# A Thought Experiment on a Model Particle — From the Phase-Area Floor to the Admissibility Condition for Dimension

**Author**: Noriaki Kihara
**Affiliation**: WF System Co., Ltd. / Osaka University, School of Engineering Science (graduate)
**ORCID**: [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)
**Version**: v1 (stealth draft)
**Date**: June 2026
**License**: CC BY 4.0
**Concept DOI**: [10.5281/zenodo.20528511](https://doi.org/10.5281/zenodo.20528511)
**Version DOI (v1.0)**: [10.5281/zenodo.20528512](https://doi.org/10.5281/zenodo.20528512)

---

## Character of this paper

**This is an observational paper, cast in the form of a thought experiment. It is not a proof paper or a claim paper.**

This paper does not modify standard quantum theory. It proposes no new physical theory. It proves no new mathematical theorem. It changes no observable prediction whatsoever. What this paper does is, through **a thought experiment that posits a single model particle**, to take two pieces of already standardly established mathematics — the Robertson uncertainty floor (the absence of zero-area states in phase space), and the Hurwitz–Frobenius series stating that normed division algebras exist only in dimensions 1, 2, 4, 8 — and **lay them side by side and observe them through the single lens of phase-space area**.

The "particle" treated here is not a specific real particle (photon, electron, etc.) but **an abstract model particle that carries a single internal phase plane**. We avoid naming a concrete particle in order to keep the thought experiment general, and to make explicit that no claim depends on a particular statistics, mass, or spin.

This paper does **not** claim:

- any change to the mathematical predictions of standard quantum theory, quantum field theory, or special relativity;
- a physical derivation that spacetime dimension is 3+1 (we observe only the mathematical admissibility condition on dimension);
- a derivation of the metric signature, Minkowski structure, or the imaginary unit;
- new physical constants, cross sections, decay rates, or new particles;
- the proof of any new mathematical theorem;
- any metaphysical claim about the origin of the universe.

Evaluation and interpretation are left to the reader.

---

## Abstract

We posit a single model particle and follow, as a thought experiment, how it can exist on the phase plane $(q,p)$. The starting points are two known facts. First, the Robertson uncertainty relation forbids a zero-area, point-like state ($\Delta q=\Delta p=0$) in phase space; a minimal symplectic area $\sqrt{\det\Sigma}=\tfrac12$ (in $\hbar=1$ units) remains as a floor. Second, the algebras whose norm is multiplicative under the product (algebras satisfying the sum-of-squares product identity) are limited to the reals, complexes, quaternions, and octonions, of dimensions 1, 2, 4, 8 (Hurwitz).

We lay these two side by side through the lens of phase-area. However, the re-readings (ii), (iii) below and in §6 are **not** derivations from standard theory but **the author's correspondence hypothesis**, connecting the vibrational degrees of freedom to four spacetime axes (made explicit as a three-layer division in the in-text section "On the second half of this paper"). (i) The model particle cannot collapse to a point and necessarily has a minimal spread of area $\tfrac12$. In the isotropic minimum-uncertainty state the per-axis standard deviation is $\tfrac{1}{\sqrt2}$; the notation $\pm\tfrac12$ used here refers not to a per-axis width but to the floor $\tfrac12$ itself (the zero-point half). (ii) If the two axes of the phase plane are nondimensionalized, both into a frequency dimension, as $(\nu_1,\nu_2)$, the area can formally be written as the bare product $\nu_1\cdot\nu_2$ without an explicit conversion constant (this is a choice of natural units and of a map; it does not mean that $\hbar$ loses its physical role, and amplitude and frequency are not a standard canonical conjugate pair). The ratio $k=\nu_2/\nu_1$ behaves as a squeeze parameter, with $k=1$ isotropic and $k\neq1$ anisotropic. (iii) Reading momentum and energy as dependent quantities of the spatial-direction and time-direction frequencies yields a **formal four-component quantity** obtained by lining up the dependent quantities of four frequencies. For this to become a physical four-momentum $(E,\boldsymbol{p})$, a Minkowski metric must be introduced separately.

Finally, we observe the admissibility condition for dimension. The phase-area floor is a continuous constraint and by itself does not fix the dimension to an integer. What pins the dimension to the sharp integer four is the existence of a division algebra — a discrete condition stricter than the floor. The conditions used to sieve are not chosen by us after the fact: they are symmetries appearing in the structure of standard quantum theory — norm preservation (unitarity), divisibility (reversible evolution), non-commutativity of rotation (angular-momentum algebra; excludes dimensions 1, 2), and associativity of composition (operator composition; excludes dimension 8) — **transposed and read onto the composition algebra of the four frequency axes**. Under this reading, the smallest dimension satisfying all of them simultaneously is four (the quaternions $\mathbb{H}$) (shown as a table over dimensions 1–6 in §6.3). We note that this paper retains, as standard, that the state space is the complex field ($\mathbb{C}$); what is linked to the quaternions is not the scalar field of states but the composition algebra of the four spacetime axes. We record the coincidence that two independent constraints — the continuous floor (cannot be a point) and the discrete ring (closes only at four) — converge on the same dimension four. Why the system reaches four and stops there (the stabilization of dimension) we do not prove; in §7 we record it only as a conjecture under the assumption that "dimension number too bears an uncertainty" (the similar prior works [12][13] are context, not grounds).

This paper makes no new claim; it confines itself to a side-by-side observation of these known structures within a single thought experiment. The introduction of a signature (Minkowski metric) is not treated here and is left as future work.

---

## §1 Setting of the thought experiment: a single model particle

### 1.1 What is assumed and what is not

Consider a single model particle. About this particle, the only thing this paper assumes is the following one point — it possesses a vibrating degree of freedom. Mass, charge, spin, statistics, and a background spacetime (externally given coordinate axes or clocks) are not assumed. We also avoid identifying it with a concrete real particle. Naming a photon or an electron would invite questions this paper does not treat — "why that particle?", "because it is a boson?", "because it is massless?" The particle here is purely an abstract unit of vibration carrying a single phase plane.

Starting from this minimal assumption, we trace in order "what it means, in phase space, for a vibrating unit to exist." All the mathematics used at each step is standardly established; this paper merely cites it and lines it up under one lens.

### 1.2 The stage: the phase plane $(q,p)$

We take the stage on which an observer reads off quantities about this particle to be the phase plane spanned by a conjugate pair $(q,p)$. The invariant that carries meaning on the phase plane is neither length nor volume but **area** (of action dimension), because $q$ and $p$ are a conjugate pair and the area-preserving transformations are the symplectic ones. This setting is common with the author's observational paper [10], and this paper proceeds with the thought experiment on that framework. Hereafter, "area" consistently means the symplectic area $\sqrt{\det\Sigma}$ ($\Sigma$ is the covariance matrix).

---

## §2 Cannot become a point — the floor $\tfrac12$ and the internality of the phase plane

### 2.1 What the Robertson floor forbids

The Robertson uncertainty relation (Robertson 1929 [2]) is

$$\Delta A\cdot\Delta B \ge \tfrac12\bigl|\langle[A,B]\rangle\bigr|,$$

and substituting $A=q,\ B=p,\ [q,p]=i$ ($\hbar=1$) gives $\Delta q\cdot\Delta p\ge\tfrac12$. More precisely, the $\mathrm{Sp}(2,\mathbb{R})$-invariant symplectic area

$$\sqrt{\det\Sigma}=\sqrt{(\Delta q)^2(\Delta p)^2-\mathrm{Cov}(q,p)^2}\ \ge\ \tfrac12$$

has a floor (de Gosson's quantum blob [4]).

Let us state this precisely. What is forbidden is a **zero-area state** in phase space (a point-like state satisfying $\Delta q=\Delta p=0$), not the origin as a mean value $(\langle q\rangle,\langle p\rangle)=(0,0)$. The latter is an ordinary, attainable state. Thus our observation is "a zero-area point-like state does not exist," not "the coordinate origin does not exist."

### 2.2 For vibration to exist, an area is required

For the model particle to vibrate means that phase advances on the phase plane. The advance of phase (rotation) requires an **area** to rotate within, and an area is spanned by at least two axes (a conjugate pair). On a single axis, rotation = vibration cannot be defined.

Here we note the first branch point of the thought experiment. If we go looking for the area that supports the particle's vibration **outside** the particle (in a background spacetime), we smuggle in the unassumed structure of a background spacetime. This paper avoids that and regards the phase plane as **internal** to the particle. That each mode of a field carries the structure of a harmonic oscillator (the conjugate pair $q,p$) is a standard fact, and this paper reads the vibration on this internal area. That is, no external partner is needed to turn the vibration; we assume that the area is present inside the model particle from the start.

### 2.3 The floor $\tfrac12$ and its zero-point notation $\pm\tfrac12$

In the isotropic minimum-uncertainty state with zero correlation ($\mathrm{Cov}(q,p)=0$), the per-axis standard deviation is $\Delta q=\Delta p=\tfrac{1}{\sqrt2}$, and their product equals the floor $\Delta q\cdot\Delta p=\tfrac12$.

Let us define the terminology precisely. What this paper writes hereafter as $\pm\tfrac12$ is **not the per-axis standard deviation** (that is $\tfrac{1}{\sqrt2}\approx0.707$, not $\tfrac12$). What $\pm\tfrac12$ denotes is **the area floor $\tfrac12$ itself** — the irreducible half-quantum corresponding to the zero-point $\tfrac12$ of semiclassical quantization $\oint p\,dq=(n+\tfrac12)h$. The $\pm$ is a sign expressing that this half-quantum floor acts in both of the two conjugate directions; it is not a numerical per-axis width. Hence it is not an object to be multiplied as in $0.5\times0.5$, but a notation that attributes a single quantity, the floor, to two directions.

What matters is that this floor $\tfrac12$ is not an "uncertainty of frequency" measured by an external clock or ruler. It is the internal spread — the minimal area of the phase plane — itself, expressing that the model particle cannot collapse to a point.

---

## On the second half of this paper — the boundary between known mathematics and the author's correspondence hypothesis

Because the character of this paper branches from here, we divide the discussion into three layers and state the standing of each at the outset. This is a premise for reading §3–§7.

**Layer one: established mathematics (the object of observation).** The Robertson uncertainty relation used in §2 and the Hurwitz–Frobenius–type dimensional restriction on normed division algebras (1, 2, 4, 8) referenced in §6 are both known, established mathematical facts. This paper does not prove them; it cites and juxtaposes them.

**Layer two: the author's correspondence hypothesis (a re-reading peculiar to this paper).** The parts that connect this mathematics from "the vibrational degree of freedom internal to the model particle" to "four spacetime axes" — taking the two axes of the phase plane as amplitude and frequency (§3), aligning both axes as homogeneous dimensionless quantities (§4), corresponding the four frequencies to a four-momentum (§5), and transposing the symmetries of quantum theory onto the composition algebra of the four axes to select the quaternions (§6.3) — are **not derivations from standard theory**. They are a **working hypothesis (correspondence assumption)** based on the author's structural intuition, whose validity this paper does not prove. In particular, that amplitude and frequency form a standard canonical conjugate pair, and that the composition algebra of the four spacetime axes must be a normed division algebra, are both assumptions imposed by this paper, not requirements of standard quantum theory.

**Layer three: the more speculative part.** §7's "dimension number too bears an uncertainty" and "stabilizes at four" are the author's non-derivational conjectures, going one step further than layer two.

Given the above, the connections to "four frequencies," "quaternions," and "four spacetime axes" that follow should be read not as derivation but as a **correspondence hypothesis**. What this paper records is the structural coincidence that two mathematics of different origin (the phase-area floor and the discrete dimensional series of division algebras) appear to point to the same dimension four under the single lens of the model particle's vibrational degree of freedom. This paper does not call this coincidence a derivation, nor does it derive that real spacetime is four-dimensional. Omitting this stance would obscure the gist of the paper (what is observed and what is posited as a hypothesis), so we state it explicitly here rather than cut it.

---

## §3 The two axes of the phase plane — amplitude and frequency, time as a dependent label

### 3.1 Choice of labels for the conjugate axes

What labels do we attach to the two axes of the phase plane? Conventionally one takes $(q,p)=$ (position, momentum) or (time, frequency). But in this thought experiment we have not assumed a background spacetime (external clock or ruler). Placing external time on a conjugate axis would smuggle in unassumed structure.

So we choose labels that close the plane using only quantities internal to the model particle. We take one axis as **amplitude** (the size of the swing) and the other as **frequency** (the speed at which phase rotates). Amplitude is internal to the particle, and frequency is internal too.

Let us be precise here. Amplitude and frequency are **not** a canonical conjugate pair in the standard sense — a pair satisfying $[q,p]=i$. In a harmonic oscillator the frequency is a parameter of the Hamiltonian and the amplitude is the size of the state; they do not hold a canonical commutation relation (there is the time–bandwidth Gabor uncertainty, but that is different in character from a Robertson-type canonical commutation relation). This paper introduces these two axes formally and observationally as **effective two variables** describing the state of a vibrational mode; it does not claim that (amplitude, frequency) is a canonical pair. Transposing the Robertson floor onto this plane and reading it there is itself part of the correspondence hypothesis described in "On the second half of this paper." Time can be treated as a **dependent label** for reading off, from outside, how far phase has advanced on this plane. The phrasing "how many times per unit time" holds only once the rotation of the plane is projected onto an external reference. Therefore in this paper we **take the principal axes of the plane to be amplitude and frequency and treat time as a dependent quantity** — we posit this choice as an assumption. This is not the exclusive claim that "time cannot be a principal axis," but means that taking frequency as a principal axis is consistent with our framework (we have not excluded an alternative choice that places time on a principal axis).

### 3.2 Anti-correlation of amplitude and frequency

When the amplitude axis and frequency axis span an area $\tfrac12$, squeezing one stretches the other. Pushing toward large frequency (the fast, fine-vibration side) narrows the amplitude spread; pushing toward small frequency (the long-wavelength limit) stretches the amplitude side. The area stays invariant while only the shape stretches. This is a standard area-preserving transformation (squeezing); the same kind of transformation is treated in the author's [10], but the discussion below does not depend on it.

---

## §4 Nondimensionalization and the ratio $k$ — isotropy and differentiation

### 4.1 Aligning both axes to frequency

If the amplitude axis and frequency axis are treated with different dimensions, forming an area requires interposing a conversion constant ($\hbar$). So we nondimensionalize both axes into a frequency dimension and write $(\nu_1,\nu_2)$. This nondimensionalization is for expressing the area using only the phase plane internal to the model particle, a convenience that avoids interposing a dimensional conversion of background-spacetime origin. The area becomes the bare product

$$\text{area}\ \propto\ \nu_1\cdot\nu_2,$$

so that formally it can be written without an explicit conversion constant. The floor is the dimensionless $\tfrac12$. Let us state this precisely too. This is a **choice of map** — adopting natural units $\hbar=1$ and mapping the two axes onto homogeneous dimensionless variables — and does not mean that $\hbar$ loses its physical role in standard physics. Identifying the symplectic area of action dimension with the bare frequency product $\nu_1\nu_2$ properly requires specifying this map and a scale; this paper adopts it only as an observational convenience (this map too belongs to the correspondence hypothesis described in "On the second half of this paper").

### 4.2 The ratio $k$ and isotropy / anisotropy

Define the ratio of the two axes as

$$k\equiv\frac{\nu_2}{\nu_1}.$$

$k=1$ is the isotropic state in which the two axes are on equal footing (a circle on the phase plane), and $k\neq1$ is anisotropic (an ellipse). $k$ is the parameter of an area-preserving squeeze, with $k=1$ corresponding to isotropy and $k\neq1$ to anisotropic distortion — the discussion below closes with this geometry alone. (The correspondence linking this kind of squeeze to a velocity parameter is also treated in the author's [10], but this section does not depend on that correspondence.)

When anisotropy arises, one of the two axes (the side whose ratio narrows) can afterwards be named "spatial-direction frequency (wavelength side)" and the other "time-direction frequency." On this view, the distinction between space and time can be read not as a property inscribed a priori on the two axes but as a label generated by anisotropy. This, however, is this paper's reading; it is not a claim that replaces the determination mechanism of spacetime structure in standard theory (nor is it proven). Under this reading, at the isotropic ($k=1$) stage no privileged axis that could be called "spatial direction" has yet been established.

---

## §5 Four-momentum as a dependent quantity

### 5.1 Momentum is a dependent quantity of frequency

We have $p=\hbar k$ ($k$ the wavenumber = spatial-direction frequency), so momentum can be written as the spatial-direction frequency times a conversion constant (the de Broglie relation). Using this relation, this paper adopts a counting in which the degrees of freedom of the phase plane are counted not by position and momentum separately but by the **number of frequencies**. This does not deny that the standard phase space treats position and momentum as independent coordinates; it is this paper's way of recounting the same degrees of freedom with frequency in the leading role. The same structure by which time was a dependent label is adopted for momentum as well.

### 5.2 Extension to many axes and the four-momentum

We extend the two-axis $(\nu_1,\nu_2)$ phase plane to many axes and consider four frequencies $\nu_1,\dots,\nu_4$ (the reason the dimension is four is observed in §6). A dependent quantity hangs off each of the four. If we (afterwards) assign three of the four to the spatial direction and one to the time direction, the dependent quantities of the former correspond to the three components of momentum $\boldsymbol{p}=(p_x,p_y,p_z)$, and the dependent quantity of the latter to the energy $E$. $E=\hbar\nu$ and $\boldsymbol{p}=\hbar\boldsymbol{k}$ are not separate relations but the four components of the single structure "a dependent quantity attaches to a frequency." What is obtained here, however, is only a **formal four-component quantity** lining up four frequency-like variables. For this to become a physical four-momentum $p_\mu=(E,\boldsymbol{p})$, there must be, behind it, a Minkowski metric $\eta_{\mu\nu}$ and spacetime translation symmetry so that it behaves as a representation of Lorentz transformations. Since this paper excludes the signature structure as in §8.4, all that can be said at this stage is: "if a Minkowski metric and translation symmetry can be introduced separately, this four-component quantity may possibly be put in correspondence with a four-momentum."

Here we record an important observation. **In this paper's construction, up to this stage the four frequencies are treated on equal footing, and we have not yet introduced any information distinguishing which is space and which is time.** The structure in which only one component flips sign (the Minkowski metric) is, in our view, not a property introduced at this symmetric stage but something that enters only after something breaks the symmetry. What breaks that symmetry we do not treat (an open problem in §8.4); this paper treats only the symmetric four. This is not a claim that the Minkowski metric or signature structure is "later-born" in standard theory, but a limitation on how far this thought experiment goes.

---

## §6 Why are there four frequencies — the discrete condition of division algebras

### 6.1 The floor does not fix the dimension to an integer

The phase-area floor $\tfrac12$ is a continuous constraint. It guarantees the impossibility of collapsing to a point, but by itself does not fix the dimension to an integer value. If only the floor were acting, the dimension would be free to fluctuate, undetermined between three and four, four and five. To nail the dimension to the sharp integer four requires a discrete condition stricter than the floor. That is the existence condition of a division algebra.

### 6.2 The admissible dimensions are limited to 1, 2, 4, 8

That the area (norm) spanned by the two axes is preserved as a product under composition of the axes — this is the condition that "composition of vibrations does not break the norm." The algebras satisfying this condition (normed division algebras for which the sum-of-squares product identity holds) are limited to the reals $\mathbb{R}$, complexes $\mathbb{C}$, quaternions $\mathbb{H}$, and octonions $\mathbb{O}$, of dimensions 1, 2, 4, 8 (Hurwitz 1898/1923 [5]). There is no three-dimensional normed division algebra. A qualification is needed here. There are many algebraic structures handling three dimensions (rotations $\mathfrak{so}(3)$, the cross product, Clifford algebras, etc.); it is not that "three dimensions has no algebra." What does not exist is a **normed division algebra with three real basis elements**, and the intuition "three does not close (in this sense)" is supported only by that theorem. Note that even for general real division algebras with associativity dropped, the admissible dimensions are limited to 1, 2, 4, 8, as shown by algebraic topology (Bott–Milnor, Kervaire 1958 [7]).

### 6.3 From here it is not observation but the author's interpretation — the reading that points to four

**Let us state at the outset that the character changes from this section on.** Up to §6.2, we juxtaposed and observed known mathematics — the uncertainty floor, and that normed division algebras are limited to dimensions 1, 2, 4, 8. That is a juxtaposition of known structures no one can dispute. **In contrast, this section (§6.3) and the next (§7) are how the author reads those known structures — interpretation and intuition, not a claim as an observational paper.** When we say below that "four is selected," it does not mean that standard quantum theory or a mathematical theorem requires four. It is a description of the author's view that, placing certain conditions on the composition algebra of the four spacetime axes and reading it, one arrives at four. The reader should take this not as an observational result but as the author's interpretation.

#### The transposition itself is an assumption

We consider which of the admissible 1, 2, 4, 8 to assign to the algebra of the four spacetime axes. The conditions this section places are symmetries appearing in the structure of standard quantum theory, **transposed and read onto the composition algebra of the four axes**. Let us state this precisely. These symmetries (unitarity, reversibility, non-commutativity of angular momentum, associativity of operator composition) are properly **properties of the state space (complex Hilbert space) and the operator algebra**, and require nothing of the spacetime dimension. Standard quantum theory holds perfectly well on a two-dimensional ($\mathbb{C}$) Hilbert space. Therefore the **transposition itself — imposing these on the algebra of the four spacetime axes — is this paper's (unproven) move**, whose justification is not in this paper but belongs to the author's intuition. The full weight of the conclusion rests on this single move — and we take it on in the main text, not in a corner of the fine print.

#### After transposing, the conditions that actually act are essentially two

Listing the transposed conditions, the ones that act independently are the following two.

- **The product is non-commutative**: in a commutative algebra ($\mathbb{R}, \mathbb{C}$) two independent rotation planes cannot stand. This drops dimensions 1, 2.
- **The product is associative**: the composition order can be defined. The octonions $\mathbb{O}$ are non-associative, so this drops dimension 8.

"Norm preservation (unitarity)" and "reversibility (divisibility)" could also be listed as conditions, but since a normed division algebra is at the same time a division algebra, the two select the same dimension set $\{1,2,4,8\}$ and add no independent sieving (Hurwitz [5]). That is, within $\{1,2,4,8\}$ the independent conditions pointing to four are **only the two: non-commutativity and associativity**. As a table:

| Dim | Algebra | Normed division algebra (= divisible) | Product non-commutative | Product associative | All conditions |
|---|---|:---:|:---:|:---:|:---:|
| 1 | $\mathbb{R}$ | ○ | ✗ | ○ | ✗ |
| 2 | $\mathbb{C}$ | ○ | ✗ | ○ | ✗ |
| 3 | — | ✗ | — | — | ✗ |
| **4** | $\mathbb{H}$ | ○ | ○ | ○ | **○** |
| 5 | — | ✗ | — | — | ✗ |
| 6 | — | ✗ | — | — | ✗ |
| 8 | $\mathbb{O}$ | ○ | ○ | ✗ | ✗ |

(The first column collapses "normed = divisible" into one. 3, 5, 6 have no normed division algebra to begin with.) Under this reading, the smallest dimension simultaneously satisfying divisibility, non-commutativity, and associativity is four (the quaternions $\mathbb{H}$) (the associative real division algebras are, up to isomorphism, the three $\mathbb{R},\mathbb{C},\mathbb{H}$: Frobenius [6]). Consistency with the discrete structure (Hurwitz integers, $D_4$ lattice) is in §6.4.

#### Self-criticism regarding circularity

If, above, we motivate "non-commutativity" by the non-commutativity of the angular-momentum algebra $[J_i,J_j]=i\varepsilon_{ijk}J_k$ ($\mathrm{SU}(2)$), that non-commutativity derives from the rotations $\mathrm{SO}(3)$ of three-dimensional space. Then we would be deriving four while presupposing three-dimensional space — **circular**. To avoid this, the non-commutativity we make act here must be purified to **the non-commutativity of the quaternion product**. But the non-commutativity of the angular-momentum algebra and that of the $\mathbb{H}$ product are structurally similar yet **not identical**, and identifying them is itself this paper's (unproven) view. On this point too, this section remains the author's reading rather than a derivation.

#### Relation to standard quantum theory

To repeat, for safety: that superposition of states and probability amplitudes are written with complex numbers (complex Hilbert space) is retained as standard, and what is linked to the quaternions is not the scalar field of states but the composition algebra of the four spacetime axes onto which this paper has transposed the conditions. This section contradicts no prediction of standard quantum theory. The process by which the system expands from one dimension to four, and the reason it does not proceed beyond five dimensions, are treated as the author's conjecture in §7.

### 6.4 The four as discrete integers

The integer version of the quaternions comes in two kinds, and a distinction is needed here. The Lipschitz integers (quaternions all of whose coordinates are integers) correspond naturally to the $\mathbb{Z}^4$ lattice. The Hurwitz integers add to those the half-integer-coordinate points, a more symmetric four-dimensional lattice structure closely related to the $D_4$ lattice. That is, $D_4$, the Hurwitz integers, and $\mathbb{Z}^4$ are mutually related but **not simply the same object** ($D_4$ should not be identified with a mere $\mathbb{Z}^4$). The uniqueness of the densest packing of the $D_4$ lattice is classically known [9]. Whereas the complexes (Gaussian integers) are the two-dimensional discrete version, four appears as the direction in which rotation is non-commutative and associative and in which all directions mesh on the integer lattice. The four frequencies $\nu_1,\dots,\nu_4$ can be read as quantities riding on these discrete integers. The notation $\pm\tfrac12$ acts on each of the four as the minimal continuous slack that keeps an integer from collapsing to a point (an integer $n$ plus the irreducible half-floor, i.e., the $\tfrac12$ of $n+\tfrac12$).

### 6.5 Coincidence with an independent series on the physics side (reference)

For reference, this also meshes with an independent theorem on the physics side. It is known that supersymmetric Yang–Mills theory holds in spacetime dimensions limited to 3, 4, 6, 10, with the reason lying in the existence of normed division algebras in the dimensions 1, 2, 4, 8 two below (Baez–Huerta [8]). This paper neither derives nor claims this correspondence, but records that the view of four frequencies riding on the quaternions sits naturally atop a known series.

---

## §7 A conjecture about why it stabilizes at four

**This section is not a proof but a conjecture.** Up through §6.3 we confirmed the static view (under the correspondence hypothesis) that "transposing the symmetries of standard quantum theory onto the algebra of the four spacetime axes and reading them, the smallest dimension satisfying them is four." What remains is the question deferred at the end of §6.3 and in §8.4 — why the system reaches four and stops there. To this, this section records not a definite answer but a single conjecture, as an observation.

### 7.1 Assumption: dimension number too bears an uncertainty

Here this paper places one assumption explicitly (not a derivation). Position and momentum, as a conjugate pair, bear an uncertainty floor and cannot collapse to a zero-area point-like state — **by analogy** with this, we assume that **dimension number too is a fluctuating quantity and, unless a special symmetry is broken, has fluctuations rather than a fixed integer**.

Let us clearly separate the difference in standing. The §2 floor comes out rigorously, as the Robertson inequality, from the canonical conjugate pair $[q,p]=i$. By contrast, this section's "uncertainty of dimension number" is, insofar as we have not specified what observable is conjugate to the dimension number $D$, not a rigorous Robertson-type floor but an **analogy**. To make it rigorous would require additional structure — a dimension-number operator, its conjugate, a state space carrying the dimension fluctuations, dynamics of inter-dimensional transition — which this paper does not supply. Therefore the §2 floor and this section's "fluctuation of dimension" differ in standing, and the latter is merely the starting point of the author's conjecture.

There exist similar prior works for the view that treats dimension not as a fixed integer but as a dynamical, scale-dependent quantity. It is known that several independent approaches to quantum gravity show a change of effective dimension at short distances (dimensional reduction) [12], and that in causal dynamical triangulations the spectral dimension runs in a scale-dependent way [13]; in particular [12] discusses this reduction by analogy with spontaneous symmetry breaking in field theory. **However, these are context showing that treating dimension as a non-fixed quantity is not outlandish; they are not grounds for this paper's assumption.** The direction is in fact opposite (they handle dimension decreasing at short distances), and this paper does not depend on these results. The above assumption is purely what this paper posits as a premise, and we observe its consequences below.

### 7.2 Conjecture: expansion by fluctuation and stabilization at four

Under this assumption, we place the following conjectures in order. For all of them this paper has no means of quantification and cannot make the probabilities explicit.

First, since there is no canonical phase plane in the usual sense at zero dimension, the §2 Robertson floor cannot be applied directly. We place this as the analogy of §7.1 — if we allow fluctuation of dimension number, "zero dimension" can be read not as a completely fixed nothing but as an indeterminate state containing the possibility of dimension generation, we assume.

Second, by this fluctuation, we assume there exists a probability that the dimension expands from zero to one or more dimensions. Neither the value of the probability nor the mechanism of expansion is supplied by this paper.

Third, we assume the dimension fluctuates without stabilizing and can expand up to four dimensions.

Fourth, once four dimensions is reached, the division-algebra condition seen in §6.3 — the condition the four-axis algebra should satisfy when the symmetries of standard quantum theory are transposed and read — acts. This is a discrete constraint stronger than the continuous uncertainty floor. Four dimensions is the smallest dimension satisfying this constraint, and here the system settles into a stable state.

Fifth, since this stability is stronger than the uncertainty fluctuation, the fluctuation cannot cross the barrier, and re-fluctuation from four dimensions to lower or higher dimensions is forbidden — or so we conjecture. That is, four is selected as the dimension "reachable by fluctuation, and once reached, not exitable by fluctuation."

To repeat, the picture itself in the fourth and fifth points — "stabilizes at four, cannot cross the barrier" — is not proven by this paper but posited as conjecture. In particular we consciously note that the fifth has the largest leap among the conjectures. The §6.3 division-algebra condition is a **static algebraic fact** (only at four do divisibility, non-commutativity, and associativity all hold); "closing" and dynamically "being unable to cross a barrier at a stable point" are different things. To assert the latter requires dynamics such as a potential or an action, which this paper does not possess. We do not quantify the metaphors of a potential valley or a stability barrier either. Evaluation is left to the reader.

### 7.3 Confluence of two constraints (recorded only)

Organizing this conjecture, dimension four appears as the confluence point of two independent constraints. One is the continuous floor (cannot be a point, hence the dimension does not collapse to zero and can fluctuate). The other is the discrete algebraic condition (when the fluctuation reaches four, it closes there for the first time, making four a stable point). The former is uncertainty, the latter division algebras — two mathematics of different origin both pointing to the same dimension four.

This paper confines itself to **recording** this coincidence, and leaves to the reader's evaluation whether to see it as accidental or necessary, or whether it can be reduced to a deeper single root. From within this paper's framework, neither of these two constraints is derived; both are cited as external established mathematics.

---

## §8 Relation to existing research, and what this paper does not claim

### 8.1 Relation to existing research

Each component of this paper is standardly established. The interpretation of uncertainty via phase-space area is Robertson 1929 [2]; the time–bandwidth version is Gabor 1946 [3]; the quantum blob as an area invariant is de Gosson [4]. That the dimensions of normed division algebras are limited to 1, 2, 4, 8 is Hurwitz [5]; that associative real division algebras are limited to $\mathbb{R},\mathbb{C},\mathbb{H}$ is Frobenius [6]; the dimensional restriction on general real division algebras is Bott–Milnor / Kervaire [7]; the correspondence between division algebras and supersymmetry is Baez–Huerta [8]; the $D_4$ lattice packing is [9]. The view of phase-space area as a common currency, and the squeeze = boost correspondence, are the author's observational paper [10]. This paper merely re-lays these out, from the phase-area lens, within the setting of a single model particle's thought experiment.

### 8.2 What this paper does not claim

- Any change to the mathematical predictions of standard quantum theory, quantum field theory, or special relativity (this paper is observationally equivalent to them).
- A physical derivation that spacetime dimension is 3+1. What this paper observes is the mathematical condition of admissible dimension for an algebra whose norm closes as a product, not the determination mechanism of real spacetime dimension. We do not equate the two.
- A derivation of the metric signature, Minkowski structure, or imaginary unit. We treat only up to the stage where the four frequencies are on equal footing (§5.2) and do not introduce the sign flip.
- New physical constants, cross sections, decay rates, or new particles.
- The proof of any new mathematical theorem. The cited theorems (Hurwitz, Frobenius, etc.) are all externally established results.
- Any metaphysical claim about the origin of the universe, singularities, or "design." §7 confines itself to a conjecture about stabilization at four under the assumption that "dimension number too bears an uncertainty" (not a proof); anything beyond is left to the reader.

### 8.3 What this paper records

- That the model particle, because of the phase-area floor $\tfrac12$, cannot take a point-like state and necessarily has a finite spread. That the per-axis standard deviation of the isotropic minimum-uncertainty state is $\tfrac{1}{\sqrt2}$, and that this paper's notation $\pm\tfrac12$ refers not to a per-axis width but to the floor $\tfrac12$ (the zero-point half) itself.
- That taking the two axes of the phase plane as amplitude and frequency and nondimensionalizing both into a frequency dimension makes the area the bare product $\nu_1\cdot\nu_2$. That in this framework one can consistently treat time as a dependent label and take the principal axes of the plane as amplitude and frequency (without excluding a choice that places time on a principal axis).
- That the ratio $k=\nu_2/\nu_1$ separates isotropy ($k=1$) from anisotropy ($k\neq1$), and that the space/time distinction can be read as a label generated by the anisotropy of the two axes (not a claim replacing the spacetime structure of standard theory).
- That reading momentum and energy as dependent quantities of spatial- and time-direction frequencies yields a formal four-component quantity lining up the dependent quantities of four frequencies; that for this to become a physical four-momentum a Minkowski metric must be introduced separately, and that at the four-frequency stage the sign asymmetry is not yet introduced.
- That what fixes the dimension to four is not the continuous floor but the discrete condition of division algebras. That transposing the symmetries appearing in the structure of standard quantum theory (norm preservation, divisibility, non-commutativity of rotation, associativity of composition) onto the composition algebra of the four frequency axes and reading them, the smallest dimension satisfying them all is four. That these symmetries are not chosen by this paper after the fact but are structures really present in standard quantum theory, while transposing them onto the four-spacetime-axis algebra is itself this paper's move, and the state space being the complex field is retained as standard (the table and caveats of §6.3).
- That, assuming "dimension number too bears an uncertainty," four can be selected as a stable point "reachable by fluctuation, and once reached not exitable by fluctuation" (§7, conjecture; not a proof, and the similar prior works [12][13] are not grounds).
- The coincidence that two independent constraints, the continuous floor and the discrete ring, both converge on dimension four (interpretation left to the reader).

### 8.4 Future work (points made explicit as unresolved)

1. **Dynamical mechanism of dimensional expansion**: what §6.3 made explicit is only the static view (under the correspondence hypothesis) that "transposing the symmetries of standard quantum theory onto the algebra of the four spacetime axes and reading them, the smallest dimension satisfying them is four." The process by which the system expands from one dimension to four while accumulating symmetries, and the reason it does not proceed beyond five dimensions, are a separate question from that fact and are not treated here. Whether the expansion is a dynamical evolution or merely a logical ordering of the conditions for existence is left to a separate paper.
2. **Mechanism of introducing the signature (Minkowski metric)**: the mechanism by which one of the four equal frequencies falls to the hyperbolic (time) side is not treated here. What breaks the symmetry (a choice of the observer's projection, or an internal asymmetry) is unresolved and left to a separate paper.
3. **Connection to the real spacetime dimension**: the relation between this paper's "four as an admissible dimension" and the observed spacetime dimension 3+1 (how the split into three space and one time is determined) is, like the signature problem of §5.2, unresolved.
4. **Making a low-dimensional toy model explicit**: writing the extension from the two axes $(\nu_1,\nu_2)$ to four as a concrete action of the quaternions $\mathbb{H}$, and showing explicitly how the floor $\tfrac12$ is distributed over each axis.
5. **Bridge to the eight-dimensional phase space**: in the standard symplectic structure, each spacetime axis $x^\mu$ has a conjugate momentum $p_\mu$, and the phase space is eight-dimensional. Folding momentum into a dependent quantity of frequency as in §5.1 makes this canonical-pair structure invisible. The correspondence between this paper's "counting by the number of frequencies" and the eight-dimensional phase space is not supplied here and is left unresolved.
6. **Integration with the author's framework**: working out in detail, in a separate paper, the consistency between this paper's four frequencies and the author's framework of $\mathbb{Z}^4$ discrete spacetime, the $D_4$ lattice, and the Standard-Model correspondence.

---

## Related work (same author)

This paper is self-contained with external standard references alone. The author's observational papers [10] (a re-reading of symplectic symmetry, Stone's theorem, and Wick rotation with phase-space area as a common currency) and [11] (an observation of structural correspondences between signal/control theory and quantum theory) are in the vicinity of the same theme, but this paper does not depend on their conclusions. The thought experiment here refers only to the phase-area framework of [10] as a premise.

---

## References

[1] W. Heisenberg (1927). *Über den anschaulichen Inhalt der quantentheoretischen Kinematik und Mechanik*. Z. Phys. **43**, 172.
[2] H. P. Robertson (1929). *The uncertainty principle*. Phys. Rev. **34**, 163.
[3] D. Gabor (1946). *Theory of communication*. J. IEE **93**, 429.
[4] M. de Gosson (2013). *Quantum blobs*. Found. Phys. **43**, 440.
[5] A. Hurwitz (1898). *Über die Composition der quadratischen Formen von beliebig vielen Variablen*. Nachr. Ges. Wiss. Göttingen, 309; (1923) Math. Ann. **88**, 1. (Normed division algebras are limited to dimensions 1, 2, 4, 8.)
[6] F. G. Frobenius (1878). *Über lineare Substitutionen und bilineare Formen*. J. Reine Angew. Math. **84**, 1. (Associative real division algebras are limited to $\mathbb{R},\mathbb{C},\mathbb{H}$.)
[7] R. Bott, J. Milnor (1958). *On the parallelizability of the spheres*. Bull. Amer. Math. Soc. **64**, 87; M. Kervaire (1958). *Non-parallelizability of the n-sphere for n>7*. Proc. Natl. Acad. Sci. **44**, 280. (Real division algebras are limited to dimensions 1, 2, 4, 8.)
[8] J. C. Baez, J. Huerta (2010). *Division algebras and supersymmetry I*. In *Superstrings, Geometry, Topology, and C\*-algebras*, Proc. Sympos. Pure Math. **81**, 65; arXiv:0909.0551.
[9] J. H. Conway, N. J. A. Sloane (1999). *Sphere Packings, Lattices and Groups*, 3rd ed. Springer. ($D_4$ lattice.)
[10] N. Kihara (2026). *Conserved Quantities and Uncertainty as Phase-Space Area — A Unified View of Symplectic Symmetry, Wick Rotation, and Stone's Theorem*. Observational paper (same author). Concept DOI: 10.5281/zenodo.20521566.
[11] N. Kihara (2026). *An Observation of Structural Correspondences between Signal/Control Theory and Quantum Mechanics, Quantum Optics, and Open Quantum Systems*. Observational paper (same author). Concept DOI: 10.5281/zenodo.20521598.
[12] S. Carlip (2017). *Dimension and dimensional reduction in quantum gravity*. Class. Quantum Grav. **34**, 193001; arXiv:1705.05417. (The view of dimension as a non-fixed effective quantity, by analogy with spontaneous symmetry breaking. Referenced as similar context for §7's assumption, not as grounds.)
[13] J. Ambjørn, J. Jurkiewicz, R. Loll (2005). *The spectral dimension of the universe is scale dependent*. Phys. Rev. Lett. **95**, 171301; arXiv:hep-th/0505113. (Scale dependence of the spectral dimension. As above, referenced as context, not grounds.)

---

Author: Noriaki Kihara / WF System Co., Ltd. / ORCID [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020) / CC BY 4.0

---

## Revision history

- **v1 (2026-06)**: First edition (stealth draft). As a thought experiment on a single model particle, juxtaposed-observation of: (§1) setting; (§2) the floor $\tfrac12$, the internality of the phase plane, and the absence of point-like states; (§3) the two axes of amplitude and frequency, and time as a dependent label; (§4) nondimensionalization and the ratio $k$, isotropy/anisotropy; (§5) the four-component quantity as a dependent quantity and the equal footing of the four; (§6) the admissibility of dimension four via the discrete condition of division algebras (Hurwitz/Frobenius, with a table transposing the symmetries of standard quantum theory onto the algebra of the four spacetime axes, making explicit that the smallest dimension satisfying all is four; the dynamical mechanism of expansion deferred); (§7) under the assumption that "dimension number too bears an uncertainty," recording as conjecture the expansion by fluctuation from zero to four and the stabilization at four (not a proof; the similar prior works Carlip [12], AJL [13] are context, not grounds), and the record of the confluence of two constraints. The introduction of a signature (Minkowski metric) was explicitly left out as future work. Refers to the framework of the author's observational papers [10][11] as a premise.
- **v1 revision (peer-review-reflected, 2026-06)**: Reflecting three AI reviews (Claude.ai / ChatGPT / Grok). (a) Added "On the second half of this paper" at the §2/§3 boundary, making explicit the three layers — established mathematics (layer one), the author's correspondence hypothesis (layer two, the re-reading onto four spacetime axes in §3–§6.3), and conjecture (layer three, §7) — and declaring the second half a "working hypothesis, not a derivation, not a claim as an observational paper." (b) Separated §2.3's $\pm\tfrac12$ from the per-axis standard deviation ($1/\sqrt2$), redefining it as the zero-point half (the floor $\tfrac12$ itself). (c) Collapsed the §6.3 table's "norm preservation = divisibility" into one column, making explicit that the independently acting conditions are the two, non-commutativity and associativity; demoted "QM requires" to "transpose and read," and self-criticized the circularity of angular-momentum non-commutativity. (d) Made explicit in §3.1 that amplitude and frequency are not a canonical conjugate pair, in §4.1 that $\hbar$'s role does not vanish, and in §5.2 that the four-momentum remains a formal four-component quantity without a Minkowski metric. (e) Qualified §6.2's "three does not close" to "limited to a normed division algebra with three real basis elements," and dissolved the identification of Lipschitz/Hurwitz/$\mathbb{Z}^4$/$D_4$ in §6.4. (f) Separated in §7.1 that the uncertainty of dimension number is an analogy because the conjugate is unspecified, and was self-aware in §7.2 of the difference between algebraic "closing" and dynamical "stability." Retitled from "toward the determination of dimension" to "toward the admissibility condition for dimension."
