# Deriving the Form of the Born Distribution from the Reproducing-Kernel Property of a Localized Odd-Harmonic Wave — Reducing the Remaining Postulates to the Squaring Rule and Randomness

**Noriaki Kihara**
(Derived from peer-review dialogue notes. This paper is not a claim of established physics, but an organization of the exact consequences that follow once a single premise is added to a mathematical observation.)

Version v0.1 (2026-06-27)
DOI (Version): 10.5281/zenodo.20965527
DOI (Concept): 10.5281/zenodo.20965526
Zenodo: https://zenodo.org/records/20965527

---

## Abstract

In quantum mechanics, the Born rule — that the observed probability is given by the square of the base wave, $|\psi|^2$ — is posited as an independent postulate in standard theory. Here we regard the **constant-amplitude odd-harmonic sum** (isolated peak wave) $S_N$ on the half-wavelength phase interval $\varphi\in[-\pi/2,\pi/2]$ as the localized kernel of observation, and we adopt the following two as **premises**: (i) observation is the convolution of the localized kernel with **phase differences**, and the probability is the **square** of the observed amplitude (the amplitude-to-probability squaring rule); (ii) phase differences finer than the band are not observed (finite-$N$ truncation). Then, because $S_N$ is the **truncated reproducing kernel** of the odd-harmonic basis on this interval, the observed amplitude **reproduces without distortion** any band-limited base wave $\psi_{\rm base}$ (exact condition $\mathrm{supp}(\hat\psi_{\rm base})\subseteq\{n\ \text{odd}:|n|\le N\}$), and the observed distribution **coincides exactly with $|\psi_{\rm base}|^2$** (after normalization). When the base wave is complex, the shifted kernel $S_N(\varphi-\varphi_0)$ is the reproducing kernel of the complex odd-harmonic basis, so the squaring becomes a genuine modulus $|Z|^2=Z\bar Z\ (\ne Z^2)$, reaching a square that includes the real–imaginary cross structure. We emphasize that this is not the Born rule of configuration-space quantum mechanics itself, but a correspondence on the half-wavelength, odd-harmonic band model.

What this paper **derives is the "form" of the distribution** (faithful reproduction by the reproducing-kernel property); the **squaring rule itself, and the probability interpretation, remain postulates**. We therefore do not claim to have "derived the Born rule" or "solved the measurement problem." Whereas standard theory posits the Born rule $|\psi|^2$ wholesale, this paper derives the form of the distribution from the reproducing-kernel property of the localized kernel, reducing the remaining postulates to the two points "squaring rule" and "the existence of randomness." The origin of randomness (why a single shot lands at one point and scatters across trials) is unsolved, and we do not touch it here.

---

## 1. Introduction

The Born rule $P(x)\propto|\psi(x)|^2$ is at the core of the predictive power of quantum mechanics, yet in the standard formulation it is required as a basic postulate. Attempts to derive "why the square" and "why the form $|\psi|^2$" include Gleason's theorem [3], which derives the squaring rule from the measure theory of Hilbert space; Zurek's approach [4], which derives it from environment-induced invariance (envariance); decision-theoretic derivations [7]; and recent work on the operational redundancy of the measurement postulates [8]. On the other hand, de Broglie–Bohm theory [5] and its quantum-equilibrium account [6] assume the initial distribution $|\psi|^2$ and then discuss dynamical relaxation.

This paper **differs in purpose** from these. We **do not attempt** to derive the squaring rule itself. Instead, we ask what the **form** of the observed distribution becomes once we add the single premise "observation = phase-difference convolution + squaring" to the single localized wave (the constant-amplitude odd-harmonic isolated peak wave) established in the starting paper [2]. The conclusion is that this form coincides **exactly** with $|\psi_{\rm base}|^2$ by the reproducing-kernel property. Thereby the "form" part of the Born squared structure ceases to be a free postulate, and the remaining postulates are reduced to the two points "the amplitude-to-probability squaring rule" and "the existence of randomness."

The contribution of this paper is modest: it is neither a claim of a new law of physics nor a claim of new mathematics. The reproducing-kernel identity we use is a classical fact of Fourier analysis (reproducing property of the Dirichlet kernel / RKHS); the novelty is not there. The contribution is the **observation** that, once one sets up the **mapping** ⟨the physical localized wave of [2] $=$ its reproducing kernel $=$ the model of observation⟩, the form of the Born distribution follows exactly under a single premise. We avoid exaggeration and make explicit at each point what is derived and what is presupposed (especially why it is not circular) (§3.1, §8, §9).

---

## 2. Setup: the localized kernel (isolated peak wave)

### 2.1 The isolated peak wave

The object of the starting paper [2] is the constant-amplitude odd-harmonic sum on the half-wavelength interval $\varphi\in[-\pi/2,\pi/2]$:

$$
S_N(\varphi)=\sum_{m=0}^{(N-1)/2}\cos\big((2m+1)\varphi\big)=\frac{\sin((N+1)\varphi)}{2\sin\varphi}
$$

($N$ odd). The central main peak is $S_N(0)=(N+1)/2$, the localization width is $\sim 1/(N+1)$, and the peak-normalized form converges to the universal $\sin u/u$ ($u=(N+1)\varphi$) (Fig. 1). These are rigorously established in [2] and are taken as the starting point here.

![Fig. 1 Convergence of the localized kernel S_N](born_localization_kernel.png)

**Fig. 1.** The localized kernel (constant-amplitude odd-harmonic isolated peak wave) $S_N(\varphi)=\sin((N+1)\varphi)/(2\sin\varphi)$. Normalizing by the peak $(N+1)/2$: (left) the main-peak width narrows as $\sim1/(N+1)$; (right) in the scaled variable $u=(N+1)\varphi$ it converges to the universal form $\sin u/u$. $N=9,99,999$.

### 2.2 Orthogonality of the basis and the reproducing-kernel property (the key)

On the interval $[-\pi/2,\pi/2]$ the odd-harmonic cosines are orthogonal:

$$
\int_{-\pi/2}^{\pi/2}\cos((2m+1)\varphi)\cos((2m'+1)\varphi)\,d\varphi=\frac{\pi}{2}\,\delta_{mm'}.
$$

Hence $S_N$ is the **truncated reproducing kernel** of this basis, truncated with all coefficients equal to 1. That all coefficients are equal (= all amplitudes 1) is the condition for being a reproducing kernel, and the "constant amplitude" property of [2] plays an essential role here. Without constant amplitude it is not a reproducing kernel, and the faithful reproduction below does not hold.

**This identity itself is a classical fact.** $S_N$ is a Dirichlet-type kernel adapted to the odd-harmonic band, and the statement "a band-limited function is faithfully reproduced by convolution with the reproducing kernel of the same band" is nothing but a standard fact of Fourier analysis (reproducing property of the Dirichlet kernel, band-limited sampling, reproducing kernel Hilbert spaces, RKHS) [10,11]. We do not claim this as a new discovery. Rather, the novelty of this paper is not the identity but the **mapping** ⟨this kernel $=$ the physical localized wave of [2] $=$ the model of observation⟩: the observation kernel was not artificially chosen to be the identity map; the point is the correspondence that the localized wave determined independently from the physics of [2] **happens** to coincide with this reproducing kernel (circularity is treated head-on in §3.1). Note that, as shown in §4–5, this kernel has **no edge effect**: reproduction is **exact pointwise over the entire interval** (including the boundary $\varphi_0=\pm\pi/2$), without the edge degradation (Gibbs-like blurring) attendant on ordinary windowed/truncated kernels. This is a consequence of each complex mode's reproduction integral closing exactly independently of $\varphi_0$ (§5.1).

In complex notation, by $\cos(n\varphi)=(e^{in\varphi}+e^{-in\varphi})/2$,

$$
S_N(\psi)=\frac12\!\!\sum_{\substack{n\ \text{odd}\\|n|\le N}}\!\! e^{in\psi},
$$

so $S_N$ is also the truncated kernel of the **complex odd-harmonic basis** $\{e^{in\varphi}:n\ \text{odd}\}$. This basis is complete on the interval: $\{e^{in\varphi}:n\ \text{odd}\}=e^{i\varphi}\cdot\{e^{i2m\varphi}:m\in\mathbb Z\}$, where the latter is complete on an interval of length $\pi$, and $e^{i\varphi}$ is unimodular and hence preserves completeness (both even and odd functions are representable). This completeness guarantees the complex extension of §5.

---

## 3. Observation model (explicit statement of the premises adopted)

We model observation as follows. This is a **premise (model choice), not a derivation**, and we say so explicitly.

**(i) Only phase differences are observed ⇒ convolution by the localized kernel + squaring.**
Observation scans the localized kernel $S_N$ across the base state (scan = convolution, depending only on the difference $\varphi_0-\varphi$) and takes the absolute square of the observed amplitude:

$$
P_N(\varphi_0)=\Big|\,(S_N*\psi_{\rm base})(\varphi_0)\,\Big|^2
=\Big|\int_{-\pi/2}^{\pi/2}S_N(\varphi_0-\varphi)\,\psi_{\rm base}(\varphi)\,d\varphi\Big|^2 .
$$

The operational motivation for this modeling is matched filtering / heterodyne detection that reads off the overlap with a localized reference wave. That the Ramsey metrology of atomic clocks is a "frequency = readout of the relative phase (beat) against a reference oscillator" [9] gives an operational grounding for "observation = phase difference / reference comparison." However, **representing measurement by this localized kernel is itself a premise** (the choice of taking the measurement basis / POVM in the family of localized waves). We do not pass this off as a derivation.

**(ii) Unobservable phase differences are ignored ⇒ finite-$N$ truncation.**
Taking phase differences finer than the band to contribute nothing to observation, we truncate the kernel at finite $N$. If $N$ covers the bandwidth of the base wave, reproduction becomes exact, as shown below.

**On the squaring rule (the core limitation of this paper).** The "squaring" in (i) above is a premise. When the convolution result $Z$ is real, $\bar Z=Z$, so $|Z|^2=Z^2$, and the power of the square (the exponent 2) is not derived independently. What this paper rigorously **derives** is the **form** part: that if one squares, the result is not distorted into something else but is **exactly** $|\psi_{\rm base}|^2$. The squaring rule (amplitude → probability) and the probability interpretation remain postulates. To read it as a probability density one assumes $\int_{-\pi/2}^{\pi/2}|\psi_{\rm base}|^2\,d\varphi<\infty$ (trivial by band-limitedness).

### 3.1 Reply to the circularity objection

The strongest anticipated objection is: "since you set observation to be convolution with the reproducing kernel (= the identity map), you merely assumed $|\psi|^2$ to obtain $|\psi|^2$." We answer this head-on. The key is to separate the **degree of freedom that could distort** from **what is accepted as a postulate**.

- **The form of the kernel is a degree of freedom; it was not chosen to be the identity map.** A general kernel (an odd-harmonic sum with non-uniform amplitudes $\sum_m a_m\cos((2m+1)\varphi)$, or another window) reproduces the base wave **with distortion** (the mode coefficient $c_m$ turns into $a_m c_m$, and if $a_m\not\equiv$ const the form changes). It is undistorted **only** when all amplitudes are equal, and that constant-amplitude kernel is precisely the localized wave that [2] derived independently from physics. Therefore the "identity map" is not an assumption but **a consequence of the fact that the physical localized wave of [2] coincides with the reproducing kernel.** Moreover, the negative control of §6 (the envelope reading is non-integrable as $1/\sin^2$) independently supports that "projection + squaring" is the unique operationalization of the form.
- **The only thing left as a postulate is the power of the square.** The argument above guarantees "the form is undistorted," but it does **not** derive the **squaring rule (exponent 2)** that turns amplitude into probability. Making the division of labor explicit: §4–§6 fix the **form** (the content of $|\cdot|$, namely $\psi_{\rm base}$, and the choice of projection vs. envelope), but the **power 2** is what was declared as a postulate in §3.

In short, if there is any circularity, it is **only in the squaring rule**, and that has been declared a postulate from the outset. What this paper newly **derives** is "that the form is undistorted," and this is a nontrivial consequence (it fails unless the amplitudes are equal) for the real degree of freedom that is the choice of kernel.

![Fig. 2 Mechanism: an equal-amplitude kernel reproduces without distortion, a non-equal-amplitude kernel distorts](born_mechanism.png)

**Fig. 2 (mechanism).** Why reproduction occurs, and why equal amplitude is essential. With the same base wave $\psi_{\rm base}=\cos\varphi+0.5\cos3\varphi-0.3\cos5\varphi$ ($N=9$): (top) the equal-amplitude kernel $S_N$ (= reproducing kernel) samples the base wave at each $\varphi_0$ ((1)), reproduces it without distortion ((2)), and squares to coincide with $|\psi_{\rm base}|^2$ ((3)). (bottom) the non-equal-amplitude kernel $\tilde S_N$ with amplitudes tapered to $1/(2m+1)$, by the same procedure, reproduces **with distortion** ((2′)) and **does not coincide** with $|\psi_{\rm base}|^2$ ((3′)). Numerically, too, the reproduction deviation is $\sim10^{-16}$ (machine precision) for equal amplitude vs. $0.53$ for non-equal amplitude. A visualization of the reason it is not circular (§3.1), against the real degree of freedom that is the choice of kernel.

---

## 4. Main result: faithful reproduction of the form by the reproducing-kernel property (real, even base)

### 4.1 The kernel identity (exact by symbolic computation)

Take the lowest odd harmonic $\cos\varphi$ as the base wave. The kernel identity

$$
\int_{-\pi/2}^{\pi/2}\cos\big((2m+1)(\varphi_0-\varphi)\big)\cos\varphi\,d\varphi=\frac{\pi}{2}\cos\varphi_0\cdot\delta_{m,0}
$$

holds. All $m\ge1$ vanish by orthogonality, and only $m=0$ (the fundamental) survives. By symbolic integration with `sympy` we verified $m=0,\dots,4$, obtaining $\pi\cos\varphi_0/2$ for $m=0$ and zero otherwise (Appendix A).

### 4.2 Projection + squaring = the square of the base wave (exact for all $N$)

From the identity, every term of $S_N$ with $m\ge1$ vanishes, and

$$
\boxed{\;(S_N*\cos)(\varphi_0)=\frac{\pi}{2}\cos\varphi_0
\quad\Longrightarrow\quad
\big|(S_N*\cos)(\varphi_0)\big|^2=\frac{\pi^2}{4}\cos^2\varphi_0\;}
$$

holds **exactly for every odd $N\ge1$**. No $N\to\infty$ limit is needed, and moreover it is **exact pointwise over the whole range of $\varphi_0$ (including the boundary)** (no edge effect; since, as in §5.1, the reproduction integral closes independently of $\varphi_0$, there is none of the edge degradation seen in ordinary truncated kernels). Numerical verification gives, for $N=1,9,31,99$, a normalized deviation $\sim10^{-16}$ (machine precision, Table 1). Fig. 3(a) shows the points for $N=1,9,99$ all lying on the target curve $\cos^2\varphi_0$.

### 4.3 General (real, even) band-limited base wave

For $\psi_{\rm base}=\sum_m c_m\cos((2m+1)\varphi)$ ($c_m\in\mathbb R$), provided $N$ covers its highest degree (exact condition $\mathrm{supp}(\hat\psi_{\rm base})\subseteq\{n\ \text{odd}:|n|\le N\}$),

$$
\big|(S_N*\psi_{\rm base})(\varphi_0)\big|^2=\frac{\pi^2}{4}\,|\psi_{\rm base}(\varphi_0)|^2 .
$$

Numerical verification ($\psi_{\rm base}=\cos\varphi+0.5\cos3\varphi-0.3\cos5\varphi$, highest degree 5): for $N\ge5$ satisfying the band condition, the deviation is $\sim10^{-15}$; for $N<5$ the 5th mode is out of band and is dropped, giving a mismatch (Fig. 3(b), Table 1). The correspondence — premise (ii) "ignore unobservable phase differences" = truncation, exact once it covers the band — is confirmed numerically.

---

## 5. Complex extension: the genuine Born modulus $|Z|^2=Z\bar Z\ne Z^2$

The base waves of §4 are restricted to even and real because the odd harmonics are **cosines**, in which case $|Z|^2=Z^2$ and the complex structure intrinsic to the Born rule (the real–imaginary cross terms) does not appear. This is a degenerate case. This section extends to complex base waves exactly, using the fact that the shifted kernel $S_N(\varphi-\varphi_0)$ is the reproducing kernel of the **complex** odd-harmonic basis.

### 5.1 Mode-by-mode reproduction (the core)

For a complex mode $\psi(\varphi)=e^{ik\varphi}$ ($k$ odd, $|k|\le N$), using the representation of §2.2 and the orthogonality relation

$$
\int_{-\pi/2}^{\pi/2}e^{i(k-n)\varphi}\,d\varphi=\pi\,\delta_{n,k}\qquad(k-n\ \text{is even, so}\ \sin((k-n)\pi/2)=0)
$$

we get

$$
(S_N*e^{ik\cdot})(\varphi_0)=\frac12\sum_{n}e^{-in\varphi_0}\!\int_{-\pi/2}^{\pi/2}\!e^{i(k-n)\varphi}d\varphi
=\frac{\pi}{2}\,e^{ik\varphi_0}.
$$

By linearity, for a band-limited complex wave $\psi_{\rm base}=\sum_{|k|\le N}c_k e^{ik\varphi}$ ($c_k\in\mathbb C$),

$$
\boxed{\;(S_N*\psi_{\rm base})(\varphi_0)=\frac{\pi}{2}\,\psi_{\rm base}(\varphi_0)
\quad\Longrightarrow\quad
\big|(S_N*\psi_{\rm base})(\varphi_0)\big|^2=\frac{\pi^2}{4}\,|\psi_{\rm base}(\varphi_0)|^2
=\frac{\pi^2}{4}\big(\mathrm{Re}^2+\mathrm{Im}^2\big)\;}
$$

holds exactly.

### 5.2 Consequence

If $\psi_{\rm base}$ is complex, the convolution result is also complex, and the squaring becomes a genuine modulus $|Z|^2=Z\bar Z\ne Z^2$. Numerical verification ($\psi_c=e^{i\varphi}+(0.5-0.3i)e^{i3\varphi}+(0.2+0.4i)e^{-i5\varphi}$): for $N\ge5$, $|conv|^2/\mathrm{norm}$ agrees with $|\psi_c|^2=\mathrm{Re}^2+\mathrm{Im}^2$ to a deviation $\sim10^{-15}$ (Fig. 3(d), Table 1). Furthermore, for a pure complex mode $\psi=e^{i\varphi}$, $|Z|^2/\mathrm{norm}=1$ (constant $=\cos^2+\sin^2$), whereas $Z^2$ oscillates as $e^{2i\varphi_0}$ ($\max|\mathrm{Im}(Z^2)|=1$); we verified this and thereby made explicit that the squaring here is $|\cdot|^2$, not $(\cdot)^2$.

Hence there is no need to narrow the scope to "even, real." The method holds for (complex) $L^2$ base waves on the half-wavelength interval, in the odd-harmonic expansion. We do not use the vague phrase "arbitrary base wave"; we state the scope in this form.

---

## 6. Negative control: the envelope reading cannot be a probability

The reading "average the fast oscillation $\sin^2((N+1)\varphi)$ (→ $1/2$) and take the envelope as the observed distribution" gives an observation weight $\propto 1/((N+1)\sin\varphi)^2$, which diverges as $1/\varphi^2$ at the center $\varphi=0$ and **cannot even be normalized** (cannot be a probability distribution). Its shape is also $1/\sin^2$, different from $\cos^2$ (Fig. 3(c)). This confirms that the correct operationalization of the two premises is not the envelope average but **projection + squaring**, reinforcing that the choice of our observation model is not arbitrary.

![Fig. 3 Main result (4 panels)](born_from_localization.png)

**Fig. 3.** Convolution by the localized kernel with phase differences + squaring reproduces the Born distribution $|\psi_{\rm base}|^2$ of the base wave. (a) for base $=\cos$, projection + squaring agrees exactly with $\cos^2$ for all $N$ ($=1,9,99$). (b) general real base (modes 1,3,5) is exact for $N\ge5$, the 5th mode dropped for $N=3$. (c) negative control = envelope $1/\sin^2$ diverges at the center, non-normalizable, different from $\cos^2$. (d) for a complex base, the genuine modulus $|Z|^2=\mathrm{Re}^2+\mathrm{Im}^2$ ($N=5,31$ agree, $N=3$ dropped). Horizontal axis in % of the half-wavelength.

---

## 7. Making the mechanism transparent: the role of the sine components

The shifted kernel decomposes as $S_N(\varphi_0-\varphi)=\sum_m\big[\cos\cdot\cos+\sin\cdot\sin\big]$, i.e. the reproducing kernel (the $\cos\cos$ part) plus a $\sin\cdot\sin$ part.

- **Even, real base (§4):** against an even $\psi_{\rm base}$ the sine harmonics are odd × even = odd, integrating to zero over the symmetric interval. Only the cosine harmonics contribute to reproduction. This is the true nature of the "cancellation," and in this degenerate case $|Z|^2=Z^2$.
- **General, complex base (§5):** the sine harmonics do **not** vanish; they reproduce the **odd / imaginary parts.** The mode-by-mode reproduction $(S_N*e^{ik\varphi})(\varphi_0)=(\pi/2)e^{ik\varphi_0}$ transparently shows the roles of both the cosine and sine components.

Numerically too, with a wave that is real but has an odd component, $\psi=\cos\varphi+0.7\sin\varphi$, we verified that the sine harmonics correctly reproduce the odd part $0.7\sin\varphi$ and are not cancelled (deviation $\sim10^{-16}$ for $N=1,9$, Table 1).

---

## 8. Numerical verification

`born_from_localization.py` includes symbolic computation (`sympy`), numerical convolution, the negative control, and the complex extension, reproducing all results of this paper to machine precision.

**Table 1. Maximum normalized deviation (machine-precision agreement in all tests)**

| Test | Base wave | Condition | $\max$ deviation |
|:--|:--|:--|:--|
| §4.2 projection + squaring | $\cos\varphi$ (real, even) | all $N=1,9,31,99$ | $\sim3\times10^{-16}$ |
| §4.3 general (real) | $\cos+0.5\cos3-0.3\cos5$ | $N\ge5$ | $\sim10^{-15}$ |
| §4.3 dropout | same | $N=3$ (5th mode missing) | $O(1)$ mismatch |
| §5 complex | $e^{i\varphi}+(0.5{-}0.3i)e^{i3\varphi}+(0.2{+}0.4i)e^{-i5\varphi}$ | $N\ge5$ | $\sim2\times10^{-15}$ |
| §5.2 $|Z|^2\ne Z^2$ | $e^{i\varphi}$ | $|Z|^2{=}1$ const, $Z^2$ oscillates | — |
| §7 odd component | $\cos\varphi+0.7\sin\varphi$ | $N=1,9$ | $\sim9\times10^{-16}$ |
| §6 negative control | envelope $1/\sin^2$ | diverges at center | non-normalizable (breaks) |

---

## 9. Discussion — what is derived and what is presupposed

We separate the roles of each element:

- **Constant-amplitude odd harmonics (the property of the starting point [2])** = truncated reproducing kernel ⇒ the convolution **reproduces the base wave as is.** Without constant amplitude it is not a reproducing kernel and does not reproduce. The "all amplitudes 1" of [2] becomes necessary here.
- **"Only phase differences are observable" (premise (i))** = convolution onto the difference + squaring ($Z\bar Z$) ⇒ makes the distribution $|\psi|^2$ (Born) rather than $\psi$.
- **"Ignore unobservable phase differences" (premise (ii))** = finite-$N$ truncation ⇒ exact once it covers the band.

**What is derived (a mathematical fact, with no extra assumptions):** the **form** $|\psi_{\rm base}|^2$ of the observed distribution comes out exactly from the reproducing-kernel property. For a complex base it comes out as a genuine modulus including the real–imaginary cross terms.

**What remains as a premise:** (a) the **squaring rule** of amplitude → probability (the exponent 2 is not derived independently); (b) the **measurement model** representing observation by the convolution of the localized kernel (the choice of POVM / measurement basis); (c) the **existence of randomness** (why a single shot lands at one point and scatters across trials). In particular (c) is an unsolved problem that this paper does not touch.

Contrast with prior work: Gleason [3], Zurek [4], and the decision-theoretic derivations [7] attempt to **derive** the squaring rule itself from deep assumptions (measure-theoretic noncontextuality / envariance / rational preferences). This paper does not attempt that; it **only** shows that, leaving the squaring rule as a postulate, the form of the distribution coincides with $|\psi_{\rm base}|^2$ by the reproducing-kernel property. Compared with de Broglie–Bohm [5] / quantum equilibrium [6], which assume the initial distribution $|\psi|^2$, the premise left here (phase-difference observation = self-comparison) is operationally grounded in atomic-clock metrology [9], but it remains a postulate nonetheless.

Hence the precise claim is:

> Adopting as a **postulate** that "observation is the phase-difference convolution by the localized kernel, and the probability is the square of that amplitude," the reproducing-kernel property of the constant-amplitude odd-harmonic localized wave makes the observed distribution coincide **exactly** with $|\psi_{\rm base}|^2$ of a (complex) band-limited base wave on the half-wavelength interval. The **form** of the Born square is derived from the localized kernel, and the remaining postulates are reduced to "the squaring rule" and "the existence of randomness."

---

## 10. Conclusion

From a single localized kernel (the constant-amplitude odd-harmonic isolated peak wave [2]) and the single premise "observation = phase-difference convolution + squaring," we have shown, by the reproducing-kernel property, that the form of the observed distribution coincides exactly with the (post-normalization) distribution $|\psi_{\rm base}|^2$ of the base wave. For a complex base the squaring becomes a genuine modulus $|Z|^2=Z\bar Z\ne Z^2$, corresponding to the complex structure of the Born rule on the half-wavelength, odd-harmonic band model. No other framework is needed; it is self-contained with a single localized kernel and two premises.

This paper has neither "derived the Born rule" nor "solved the measurement problem." It is a **reduction of scope**: of the Born squared structure, it derives the "form" from the reproducing-kernel property of the localized kernel and narrows the remaining postulates to "the squaring rule" and "the existence of randomness." Whereas standard theory posits the Born rule wholesale, the postulates here are one step lighter. The origin of randomness remains unsolved.

---

## Appendix A: `sympy` output of the kernel identity

```
I_m(phi0) = int_{-pi/2}^{pi/2} cos((2m+1)(phi0-phi)) cos(phi) dphi
  m=0:  pi*cos(phi0)/2
  m=1:  0
  m=2:  0
  m=3:  0
  m=4:  0
```

Only $m=0$ gives $(\pi/2)\cos\varphi_0$; all $m\ge1$ are zero (orthogonality). Hence $(S_N*\cos)(\varphi_0)=(\pi/2)\cos\varphi_0$ for all odd $N$.

## Appendix B: Reproduction

```bash
python3 born_from_localization.py     # symbolic + numeric + complex + negative control (machine precision)
python3 born_fig2.py                   # Fig. 1 (born_localization_kernel.png/svg, convergence of the kernel)
python3 born_mechanism.py              # Fig. 2 (born_mechanism.png/svg, mechanism: equal vs non-equal amplitude)
python3 born_fig.py                    # Fig. 3 (born_from_localization.png/svg, 4 panels)
```

## Figures (list)

In the text, Fig. 1 is placed just after §2.1, Fig. 2 just after §3.1, and Fig. 3 just after §6 (numbered in order of appearance). Files (PNG = display, SVG = vector master):

- **Fig. 1** (property of the kernel): `born_localization_kernel.png` / `.svg`. The localized kernel $S_N$, peak-normalized, (left) narrows as $\sim1/(N+1)$, (right) converges to the universal form $\sin u/u$ in $u=(N+1)\varphi$. Consistent with the figure of the starting paper [2].
- **Fig. 2** (mechanism): `born_mechanism.png` / `.svg`. The equal-amplitude kernel samples → reproduces without distortion → squares to coincide with $|\psi_{\rm base}|^2$ (top); the non-equal-amplitude kernel $1/(2m+1)$ distorts and mismatches by the same procedure (bottom). A visualization of "undistorted only for equal amplitude" (§3.1). Reproduction deviation: equal amplitude $\sim10^{-16}$ / non-equal $0.53$.
- **Fig. 3** (verification of reproduction): `born_from_localization.png` / `.svg`. (a) for base $=\cos$, projection + squaring agrees exactly with $\cos^2$ for all $N$; (b) general real base (modes 1,3,5) exact for $N\ge5$, dropped for $N=3$; (c) negative control = envelope $1/\sin^2$ diverges at center, non-normalizable; (d) for a complex base, the genuine modulus $|Z|^2=\mathrm{Re}^2+\mathrm{Im}^2$ ($N=5,31$ agree, $N=3$ dropped).

## References

[1] M. Born, "Zur Quantenmechanik der Stoßvorgänge," *Z. Phys.* **37**, 863 (1926).
[2] N. Kihara, "An Observation on the Isolated Peak Wave of a Constant-Amplitude Odd-Harmonic Sum on a Half-Wavelength Phase Interval and Its Localization," Zenodo, v0.4 (2026-06-25), Concept DOI: 10.5281/zenodo.20833096.
[3] A. M. Gleason, "Measures on the closed subspaces of a Hilbert space," *J. Math. Mech.* **6**, 885 (1957).
[4] W. H. Zurek, "Probabilities from entanglement, Born's rule $p_k=|\psi_k|^2$ from envariance," *Phys. Rev. A* **71**, 052105 (2005).
[5] D. Bohm, "A Suggested Interpretation of the Quantum Theory in Terms of 'Hidden' Variables. I, II," *Phys. Rev.* **85**, 166, 180 (1952).
[6] A. Valentini and H. Westman, "Dynamical origin of quantum probabilities," *Proc. R. Soc. A* **461**, 253 (2005).
[7] D. Deutsch, "Quantum theory of probability and decisions," *Proc. R. Soc. A* **455**, 3129 (1999); D. Wallace, *The Emergent Multiverse* (Oxford, 2012).
[8] L. Masanes, T. Galley, M. Müller, "The measurement postulates of quantum mechanics are operationally redundant," *Nat. Commun.* **10**, 1361 (2019).
[9] N. F. Ramsey, "A Molecular Beam Resonance Method with Separated Oscillating Fields," *Phys. Rev.* **78**, 695 (1950) (Ramsey / heterodyne phase readout as reference metrology).
[10] E. M. Stein and R. Shakarchi, *Fourier Analysis: An Introduction* (Princeton Univ. Press, 2003) (Dirichlet kernel and reproduction/sampling of band-limited functions); Y. Katznelson, *An Introduction to Harmonic Analysis*, 3rd ed. (Cambridge Univ. Press, 2004).
[11] N. Aronszajn, "Theory of reproducing kernels," *Trans. Amer. Math. Soc.* **68**, 337 (1950) (reproducing kernel Hilbert spaces, RKHS).

---

*(This paper derives from peer-review dialogue notes (Noriaki Kihara × Iris, 2026-06-27). It is not a claim of established physics but the consequence of a mathematical observation plus an explicitly stated premise. The paper is self-contained with a single localized kernel and two premises, and states explicitly that the origin of randomness is unsolved. The verification code `born_from_localization.py` reproduces all results to machine precision.)*
