# Curvature Oscillation Readout via Area Cross-Terms in Complex Rotating Waves

## ―― Minimal Application to Gravitational Waves

Noriaki Kihara (木原範昭)\
WF System Co., Ltd.\
ORCID: 0009-0004-6753-4020

Type: Research Note / Hypothesis and Observational Paper\
Date: September 1, 2026\
Version DOI: 10.5281/zenodo.22230941\
Concept DOI: 10.5281/zenodo.22230941

------------------------------------------------------------------------

## Abstract

Previous work [1] established that for a geodetic unit cell placed in a positively curved constant-curvature space, the geodetic line length at the boundary is conserved independent of the radius of curvature, while curvature distortion appears for the first time from the area spanned by two directions. This note applies that result to the complex rotating wave

$$z=a+ib$$

Upon squaring, we obtain

$$z^2=(a^2-b^2)+2iab$$

The cross term $(2iab)$ is an area-type quantity arising from the product of two independent components $(a,b)$. Furthermore, setting

$$a=A\cos\theta,\qquad b=A\sin\theta$$

we get

$$2ab=A^2\sin 2\theta$$

This area-type quantity is not static; it oscillates at twice the original phase angle.

In this note, we combine the previous geometric result that curvature effects emerge from two-dimensional area rather than one-dimensional geodetic lines with the area cross-term in the square of a complex wave. We propose a minimal hypothesis: the area oscillation of a complex rotating wave in the imaginary direction couples with curvature and is read out from outside as a real-valued oscillation, with gravitational waves as one candidate manifestation. Using GW150914 as a representative example, if we adopt the Schwarzschild radius as the external curvature scale, the natural frequency scale

$$c/(2\pi r_s) \approx 2.5 \times 10^2 \text{ Hz}$$

falls within the same order as the observed frequency band of 35--250 Hz. Furthermore, $r_s/D \approx 1.5 \times 10^{-20}$, and the difference from the observed peak strain $1.0 \times 10^{-21}$ is a dimensionless coefficient of approximately $6.5 \times 10^{-2}$. However, this latter coefficient is not derived from this note; we limit verification to dimensional consistency rather than amplitude prediction.

------------------------------------------------------------------------

## 1. Motivation—Curvature Distortion Begins from Area, Not Length

Previous work [1] evaluated the distortion of a regular geodetic cell of unit side length placed in a positively curved constant-curvature space $S^d(R)$, separately for side length, angle, area, and volume.

The simplest result is that the geodetic line length at the boundary remains exactly 1 by construction. By contrast, in a geodetic square spanned by two directions, angle and area exhibit curvature dependence. The vertex angle is

$$\theta(R)=\arccos\left[-\tan^2\left(\frac{1}{2R}\right)\right],$$

and the area is

$$A(R) = R^2\left[4\theta(R)-2\pi\right]$$

In the flat limit,

$$A(R) = 1+\frac{1}{6R^2}+O(R^{-4}).$$

For $d$-dimensional volume, we have

$$V_d(R) = 1+\frac{d(d-1)}{12R^2}+O(R^{-4})$$

Since

$$\frac{d(d-1)}{12} = \frac{1}{6}\binom{d}{2}$$

the lowest-order volume distortion can be organized as the sum of area distortions created by independent pairs of two directions [1].

Thus the starting point for this note, derived from previous work, is

$$\boxed{\text{The first geometric gateway to curvature distortion is two-dimensional area, not one-dimensional length}}$$

In this note, we consider what happens when this result is applied to a rotating wave expressed in complex numbers.

------------------------------------------------------------------------

## 2. Squaring a Complex Rotating Wave and the Area Cross-Term

Let the complex wave be

$$z=a+ib$$

where $(a,b)$ are real components. Using amplitude $A$ and phase $\theta$, we can write

$$a=A\cos\theta,\qquad b=A\sin\theta$$

Squaring this wave gives

$$z^2 = (a+ib)^2 = (a^2-b^2)+2iab. \tag{1}$$

The real part is

$$a^2-b^2 = A^2(\cos^2\theta-\sin^2\theta) = A^2\cos2\theta, \tag{2}$$

and the coefficient of the imaginary part is

$$2ab = 2A^2\cos\theta\sin\theta = A^2\sin2\theta. \tag{3}$$

Therefore,

$$\boxed{z^2=A^2(\cos2\theta+i\sin2\theta)} \tag{4}$$

What this note focuses on is equation (3).

Whereas $a^2$ and $b^2$ are squares of single-direction components, $ab$ is the product of two independent components $(a,b)$. Following the framing of reference [1]—"curvature distortion emerges from the area spanned by two directions"—the quantity $2ab$ is naturally selected as the minimal area-type cross-product in the complex rotating wave that can couple with curvature.

Moreover, this is not a static area quantity. From equation (3),

$$2ab=A^2\sin2\theta$$

is a real-valued oscillation that undergoes sign reversal at twice the original phase $\theta$.

This fact permits the following geometric interpretation, distinct from mere notational convenience.

------------------------------------------------------------------------

## 3. Hypothesis of Curvature Readout

A one-dimensional real wave

$$x=A\cos\theta$$

represents oscillation along a single geodetic direction. In this case, according to reference [1], no area cell spanned by two directions arises.

By contrast,

$$z=A(\cos\theta+i\sin\theta)$$

rotates in the complex plane and simultaneously possesses two orthogonal components. Its square necessarily contains the area-type cross-term

$$2iab=iA^2\sin2\theta$$

### 3.1 Geometric Correspondence between the Imaginary Axis and the Curvature Direction

Let us begin with the fundamental equation of a 3-sphere:

$$x^2 + y^2 + z^2 = R^2$$

Rewriting this,

$$x^2 + y^2 + z^2 - R^2 = 0$$

Organizing by the radial direction $r^2 := x^2 + y^2 + z^2$, we obtain

$$r^2 - R^2 = 0$$

Now applying the complex substitution $R' = iR$, we get

$$r^2 - (iR)^2 = 0 \quad \Rightarrow \quad r^2 + R^2 = 0$$

This is the equation of the **light cone** (null surface).

In summary:
- **Real axis** $R$: sphere (defines positive curvature)
- **Imaginary axis** $iR$: light cone (null geodesics)

The fact that the imaginary direction is the curvature direction is not merely a formal re-interpretation but rather a **geometric necessity between sphere coordinates and light cone coordinates**. The light cone is a characteristic expression of where curvature governs spacetime structure, and the generator direction is precisely the imaginary axis $iR$.

![Geometric correspondence between sphere and light cone](sphere_lightcone_diagram.svg)

**Figure 1:** The left diagram shows a three-dimensional sphere $x^2 + y^2 + z^2 = R^2$ in real space (positive curvature). The right diagram shows the two-dimensional light cone $r^2 + (iR)^2 = 0$ after the complex substitution $R' = iR$ (null surface). The imaginary axis $iR$ being the curvature direction becomes geometrically self-evident.

With this correspondence in place, we posit the following hypothesis:

### 3.2 Hypothesis of Curvature Readout

> Hypothesis (Curvature Readout Hypothesis):\
> When a complex rotating wave exists in a space with curvature, the area-type cross-term
> $(2iab)$ spanned by two directions can couple with curvature. When the curvature direction is expressed as a direction not directly readable from the internal real-valued geodetic directions, this imaginary cross-term can be read out from external observation as a real-valued periodic oscillation.

Writing only the formal direction transformation via the imaginary unit,

$$i(2iab)=-2ab=-A^2\sin2\theta. \tag{5}$$

Accordingly, the candidate oscillation that is read out is real-valued.

In this paper, equations (1)--(4) are identities, but we interpret their imaginary direction physically as a curvature direction.

Using the area correction from reference [1], the minimal candidate in the positive-curvature case is conceptually

$$h_\text{curv} \propto (k_s(R)-1)(-2ab), \tag{6}$$

$$k_s(R) = R^2 \left[4\arccos\left(-\tan^2\frac{1}{2R}\right)-2\pi\right]. \tag{7}$$

In the weak curvature limit,

$$k_s(R)-1 = \frac{1}{6R^2}+O(R^{-4})$$

thus

$$h_\text{curv} \propto -\frac{A^2}{6R^2}\sin2\theta +O(R^{-4}). \tag{8}$$

Equations (6)--(8) are candidate couplings; the proportionality constant, normalization, and choice of whether to use $k_s$ or its reciprocal as the measure remain undetermined in this note. What is established here is that this curvature-derived component vanishes in the flat limit $R\to\infty$, and that its lowest-order curvature dependence has the same $R^{-2}$ form as the area distortion from reference [1].

If this interpretation is correct, gravitational waves should not be thought of as "the complex rotating wave itself," but rather

$$\boxed{\text{Complex rotating wave} \to \text{area cross-term} \to \text{coupling with curvature direction} \to \text{real-valued oscillation readable from outside}}$$

Moreover, the fact that the double-angle structure emerges automatically from equation (4) suggests a formal correspondence with the fact that gravitational waves possess a double-angle polarization structure unlike ordinary vector waves. However, this note does not claim that the tensor perturbation $h_{\mu\nu}$ of general relativity or spin-2 representations are derived from equation (4).

------------------------------------------------------------------------

## 4. Non-uniqueness of Observational Scale and the Schwarzschild Radius

Here we encounter one issue.

From the complex wave alone, one cannot uniquely determine the radius of curvature $(R)$. Depending on whether the same local wave is interpreted as a molecular-scale system, a stellar system, or the entire universe, the curvature scale appropriate for observational adoption may differ.

In this note, we do not hide this as a defect but explicitly state it as observational dependence in the readout.

The black hole is an object that minimizes this ambiguity. Given a mass $M$, we can adopt a natural length scale borrowed externally: the Schwarzschild radius

$$r_s=\frac{2GM}{c^2} \tag{9}$$

In this note, we do not regard this as derived from the complex wave. Rather, we use it as an externally specified readout scale given by designating a black hole as the observational target.

The corresponding natural time scale is

$$t_s=\frac{r_s}{c},$$

and the frequency scale is

$$f_s=\frac{1}{2\pi t_s} = \frac{c}{2\pi r_s}. \tag{10}$$

Because the cross-term in equation (3) has a double-angle structure relative to phase, the readout frequency could be $2f_s$ depending on the definition of the internal phase $\theta$. We do not fix this factor 2 in this note; determining the observational mapping is left as a future task.

------------------------------------------------------------------------

## 5. Order-of-Magnitude Check Using GW150914

In the initial analysis by LIGO of GW150914, the source masses of the two black holes were approximately

$$M_1=36_{-4}^{+5}M_\odot,\qquad M_2=29_{-4}^{+4}M_\odot,$$

the luminosity distance was

$$D_L=410_{-180}^{+160} \text{ Mpc},$$

the signal swept upward from approximately 35 Hz to 250 Hz, and the peak strain was

$$h_\text{peak}\simeq1.0\times10^{-21}$$

[2,3].

For an order-of-magnitude check, taking the total mass before merger as

$$M=M_1+M_2\simeq65M_\odot$$

equation (9) gives

$$r_s \simeq 1.92\times10^5 \text{ m} \simeq 192 \text{ km}. \tag{11}$$

Therefore equation (10) becomes

$$f_s = \frac{c}{2\pi r_s} \simeq 2.49\times10^2 \text{ Hz}. \tag{12}$$

This falls within the same order as the 35--250 Hz observed in GW150914, and is immediately consistent as a representative strong-gravity timescale.

Next, for a radiation field with amplitude that decays as $(1/D)$ at large distance, examining the simplest dimensionless ratio constructible from the Schwarzschild radius yields

$$\frac{r_s}{D_L} \simeq \frac{1.92\times10^5}{410\times10^6\times3.0857\times10^{16}} \simeq 1.52\times10^{-20}. \tag{13}$$

The ratio of the observed peak strain to this quantity is

$$C_\text{obs} \equiv \frac{h_\text{peak}}{r_s/D_L} \simeq 6.6\times10^{-2}. \tag{14}$$

That is,

$$h_\text{peak} \simeq 0.066 \cdot \frac{r_s}{D_L}$$

However, the coefficient (0.066) in equation (14) is not derived from the complex wave model of this note. Therefore, this is not a prediction but a dimensionless coefficient obtained by back-calculation from observations. In standard gravitational wave theory as well, the fact that strain from a strong-gravity source is related to the ratio of the source's gravitational radius to distance emerges naturally from dimensional analysis [4]. What this note has verified is the restricted fact that the natural scale proposed for curvature readout is not off by many orders of magnitude compared to a representative real event.

What is more directly interesting in this note is the frequency side. By providing only $r_s$ externally, we obtain

$$f_s\simeq249 \text{Hz}$$

which falls within the observed band of 35--250 Hz. This is not immediately grounds for rejecting the hypothesis.

------------------------------------------------------------------------

## 6. Discussion

The structure of this note is intentionally simple.

First, as a result from reference [1], we have

$$\text{1D length: no curvature distortion},$$

$$\text{2D area: first curvature distortion}$$

Second, squaring a complex rotating wave yields

$$(a+ib)^2=(a^2-b^2)+2iab$$

and the product of two directions $(2ab)$ necessarily appears.

Third,

$$2ab=A^2\sin2\theta$$

so the area-type cross-term is a periodic oscillation of real value.

Combining these three points, we obtain an extremely simple hypothesis:

$$\boxed{\text{As a candidate for the wave quantity with which curvature can first couple, }2iab\text{ is selected}}$$

In this view, the effect vanishes in flat space or for waves along one-dimensional geodetic lines. Only when a complex rotating wave that sweeps out area exists in a curved background does a vibrational component in the curvature direction arise.

This mechanism is not fundamentally limited to black holes. Similar effects can occur at sites of local gravitational curvature. However, in microscopic domains, discretization and quantization of phase and other constraints are dominant, and this effect may be difficult to observe directly. This point remains a hypothesis in this note.

Additionally, the transformation $(z^2)$ under

$$z\rightarrow e^{i\alpha}z$$

is

$$z^2\rightarrow e^{2i\alpha}z^2$$

This double-angle structure formally resembles the double-angle dependence appearing in the polarization of gravitational waves. Whether the double-angle transformation on the complex plane in equation (4) allows us to derive the tensor perturbation $h_{\mu\nu}$ of general relativity and spin-2 representations is beyond the scope of this paper and left for future work. Establishing precise correspondence requires theoretical examination of spacetime tensor structure, gauge freedom, two polarizations, propagation equations, and energy flux.

------------------------------------------------------------------------

## 7. Scope of Claims and Conclusion

The claims made in this note are limited to the following scope.

1. In reference [1], the geodetic line length itself is not distorted by curvature, and curvature distortion appears for the first time from the area spanned by two directions.
2. The square of a complex rotating wave $(z=a+ib)$ necessarily contains an area-type cross-term $(2iab)$ arising from the product of two directions.
3. If $(a=A\cos\theta, b=A\sin\theta)$, its coefficient is $2ab=A^2\sin2\theta$, a double-angle oscillation with real value.
4. Reading this cross-term as a coupling component with the curvature direction yields a minimal hypothesis: the area oscillation in the imaginary direction is read out from the outside as a real-valued curvature oscillation.
5. For black holes, adopting the Schwarzschild radius as the external readout scale, GW150914 gives $$c/(2\pi r_s)\sim2.5\times10^2 \text{ Hz}$$, which falls within the same order as the observed gravitational wave frequency band.
6. $(r_s/D_L\sim1.5\times10^{-20})$ also shows no enormous discrepancy with the observed strain $(10^{-21})$, but the proportionality constant is not derived in this note, so it should not be regarded as an amplitude prediction.

What this note presents is the **coupling mechanism between complex rotating waves and curvature** derived from

$$\boxed{\text{the geometric result from reference [1] that "curvature distortion begins from area"} + \text{the mathematical identity that "the square of a complex wave contains oscillating area cross-terms"}}$$

Specifically, we propose that when a complex rotating wave exists in curved space, the area cross-term in its imaginary direction couples with curvature and is read out from outside as an observable real-valued oscillation. That the frequency scale of GW150914 (≈250 Hz) and amplitude scale (≈$10^{-21}$) fall within the same order as the natural curvature scale constructed from the Schwarzschild radius suggests that this hypothesis is not merely a formal re-interpretation but may bear on the essence of the gravitational wave phenomenon.

Complete verification requires deriving the coupling rule in equation (6), amplitude normalization, polarization tensor, propagation rules, and energy flux from a single wave model. This note positions itself as a minimal research note preceding such investigation.

------------------------------------------------------------------------

## References

[1] Noriaki Kihara, "Distortions of Geodetic Unit Cells in Positively Curved Constant-Curvature Space—Exact Evaluation of Side Length, Angle, Area, and Volume," v1.4, Zenodo, Version DOI: 10.5281/zenodo.20684135; Concept DOI: 10.5281/zenodo.20680269 (2026).

[2] B. P. Abbott et al. (LIGO Scientific Collaboration and Virgo Collaboration), "Observation of Gravitational Waves from a Binary Black Hole Merger," Physical Review Letters 116, 061102 (2016). DOI: 10.1103/PhysRevLett.116.061102.

[3] B. P. Abbott et al. (LIGO Scientific Collaboration and Virgo Collaboration), "Properties of the Binary Black Hole Merger GW150914," Physical Review Letters 116, 241102 (2016). DOI: 10.1103/PhysRevLett.116.241102.

[4] LIGO Scientific Collaboration and Virgo Collaboration, "The basic physics of the binary black hole merger GW150914," Annalen der Physik 529, 1600209 (2017). DOI: 10.1002/andp.201600209.

------------------------------------------------------------------------
