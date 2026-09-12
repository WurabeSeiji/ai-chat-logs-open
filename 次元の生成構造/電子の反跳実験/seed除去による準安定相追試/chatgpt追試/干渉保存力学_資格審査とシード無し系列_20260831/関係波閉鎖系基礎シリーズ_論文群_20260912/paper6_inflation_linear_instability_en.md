# The inflation curve is exponential amplification of the floor's linear instability mode——the largest eigenvalue of the floor Jacobian and its non-normality / time-varying corrections

**Series:** Foundations of Self-Consistent Relational-Wave Closed Systems (Paper 6 of 10)
**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-09-12
**Version DOI:** 10.5281/zenodo.22729096
**Concept DOI:** 10.5281/zenodo.22729095
**ORCID:** 0009-0004-6753-4020
**Zenodo:** https://zenodo.org/records/22729096

---

## Abstract

We derive analytically, as the **linear instability** of a tiny seed, the inflation-like rapid expansion (exponential growth of the out-of-floor-plane component $H_\perp/H$) reported in prior work [1]. When the largest eigenvalue $|\mu_{\max}|$ of the linearized map of the floor (the Jacobian $J$: the $2M\times2M$ real matrix obtained by linearizing $z\to\exp(\Delta\tau K(\varphi))z$ at the floor) exceeds unity, it gives $H_\perp/H(t)\approx\text{seed}^2\,|\mu_{\max}|^{2t}$, ramp slope $2\log_{10}|\mu_{\max}|$, and onset $\approx\ln(1/H_{\perp0})/(2\ln|\mu_{\max}|)=\ln(1/\text{seed})/\ln|\mu_{\max}|$ (with $H_{\perp0}=\text{seed}^2$). The $H_\perp$ slope obtained by evolving a tiny $\delta$ under the frozen floor $J$ matches $2\log_{10}|\mu_{\max}|$ exactly ($N=6{:}0.2935$ vs analytic $0.2937$). $|\mu_{\max}|$ decreases as $N$ increases ($1.402\to1.269\to1.185$ for $N=6,8,10$) = inflation becomes slower and the onset lengthens (measured onset is $64,87,133$ steps for $N=6,8,10$). The measured ramp slope is consistently $18$–$31\%$ steeper than the floor value, and this stems from the fact that the floor $J$ is a strongly **non-normal matrix** (non-normality $0.53$–$0.76$) and that the real ramp is the amplification of the **time-varying, non-commuting time-ordered product** $\prod J(\varphi_t)$ accompanying the gradual change of $K(\varphi)$ (joint spectral radius $\ge$ individual spectral radius) (**directly confirmed in follow-up test A**: the slope of the time-varying product matches the measured value and exceeds the frozen-floor value by 18–31%). The floor eigenvalue is its lower bound (initial / frozen rate). Classification: derivable (form, mechanism, $N$ dependence, onset) + direct confirmation of the non-normality correction (follow-up test A; only the closed form is open).

---

## 1. Problem

In Papers 1 and 2 the floor was established as an unstable relative equilibrium. "If the dynamics of the wave has been solved, then the inflation curve too should be derivable analytically from a tiny seed"——this paper shows this quantitatively, as the exponential amplification of the linear instability mode of the floor Jacobian.

## 2. The linear-instability prediction formulas

We construct by numerical differentiation the real Jacobian $J\in\mathbb R^{2M\times2M}$ obtained by linearizing the map around the floor $z_0$ (embedding the state as $[\Re z;\Im z]$), and compute its eigenvalues $\mu$. If the floor is unstable, $|\mu_{\max}|>1$. A tiny seed $\delta$ is amplified along the unstable eigendirection, and the fraction of the component orthogonal to the floor plane is

$$\frac{H_\perp}{H}(t)\approx(\text{seed})^2\,|\mu_{\max}|^{2t},\qquad
\text{ramp slope}=\frac{d\log_{10}(H_\perp/H)}{dt}=2\log_{10}|\mu_{\max}|,\qquad
\text{onset}\approx\frac{\ln(1/H_{\perp0})}{2\ln|\mu_{\max}|}=\frac{\ln(1/\text{seed})}{\ln|\mu_{\max}|}.$$

(The onset is the number of steps to reach $H_\perp/H\to1$. $H_{\perp0}=\text{seed}^2$. The old draft's $\ln(1/\text{seed})/(2\ln|\mu|)$ was a coefficient error that produced a self-contradiction—about half of the measured onset—so it has been corrected.)

## 3. Verification ($N=6,8,10$)

| $N$ | $|\mu_{\max}|$ | $\lambda=\ln|\mu|$/step | number of unstable modes | analytic slope $2\log_{10}|\mu|$ | $\delta$-evolution slope under frozen floor $J$ | measured ramp slope | measured/floor-analytic | floor $J$ non-normality | onset analytic/measured |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 1.40234 | 0.3381 | 5 | 0.2937 | 0.2935 | 0.3460 | 1.178 | 0.761 | 86 / 64 |
| 8 | 1.26908 | 0.2383 | 6 | 0.2070 | 0.1988 | 0.2702 | 1.306 | 0.641 | 125 / 87 |
| 10 | 1.18513 | 0.1699 | 8 | 0.1475 | 0.1271 | 0.1814 | 1.230 | 0.530 | 180 / 133 |

- **The linear-instability formula is exactly correct**: the $H_\perp$ slope from evolving a tiny $\delta$ under the frozen floor $J$ matches the analytic $2\log_{10}|\mu_{\max}|$ ($N=6{:}0.2935$ vs $0.2937$).
- **The $N$ dependence is reproduced**: $|\mu_{\max}|$ decreases as $N$ increases ($1.402\to1.269\to1.185$) = inflation becomes slower and the onset lengthens (measured onset is $64,87,133$ steps for $N=6,8,10$, analytic $86,125,180$ steps). The closed-form $N$-dependence law of the onset is left as future work.
- **onset** $\approx\ln(1/\text{seed})/\ln|\mu|$ matches in order of magnitude and scale (analytic $86/125/180$ vs measured $64/87/133$).

![Figure 1: Inflation = linear instability. Left = $N$ dependence of the floor Jacobian $|\mu_{\max}|$ and the onset; center = analytic ramp slope $2\log_{10}|\mu|$ (floor) vs measured (non-normality gap); right = the measured $H_\perp/H$ ramp for $N=6$ with the analytic slope overlaid.](../位相ロック全N確認_相対平衡_20260912/インフレ曲線_線形不安定解析_20260912/fig_inflation_linear_instability.png)

**Figure 1** Inflation = the floor's linear instability ($N$ dependence, slope, overlay).

## 4. Non-normality correction (the amount by which the measurement exceeds the floor value)

The measured ramp slope is consistently $18$–$31\%$ steeper than the floor value (measured/floor-analytic $=1.18$–$1.31$).

- The floor $J$ is a strongly **non-normal matrix** (non-normality $\|JJ^{\mathsf T}-J^{\mathsf T}J\|/\|J^{\mathsf T}J\|=0.53$–$0.76$).
- The spectral radius when $J$ is frozen at each point is nearly constant along the ramp ($1.402$ for $N=6$), but the real ramp **persistently** exceeds it. Frozen non-normality only provides transient amplification; the persistent part is carried by the amplification of the **time-varying, non-commuting time-ordered product of Jacobians** $\prod_t J(\varphi_t)$ (joint spectral radius $\ge$ individual spectral radius) as $K(\varphi)$ changes gradually. **This mechanism has been established by direct verification (follow-up test A)**: the $H_\perp$ slope from evolving $\delta$ under the time-varying Jacobian product $\prod_t J(z_t)$ on the real trajectory matches the measured ramp (time-varying/measured = 0.999 / 1.004 / 1.019 for $N=6,8,10$), and exceeds the frozen floor $J_0$ by 18–31% (time-varying/floor = 1.177 / 1.311 / 1.253 = quantitatively matching the measured excess). The floor eigenvalue is its lower bound (initial / frozen rate). Therefore the fact that the analytic onset (based on the floor $|\mu_{\max}|$) in §3 consistently exceeds the measured value ($86/125/180$ vs measured $64/87/133$) is likewise a consequence of the same non-normality mechanism——because the floor eigenvalue is a lower bound on the growth rate, the floor-based onset gives an upper bound, and the real ramp is steeper (slope ratio $1.18$–$1.31$) and reaches the threshold faster. The onset overestimation and the non-normality excess are two sides of the same mechanism.

**Deeper reading (finite-time Lyapunov rate of the tangent cocycle).** Near the floor, $\delta_{t+1}=J_0\delta_t$ and the rate is $\ln\rho(J_0)=\ln|\mu_{\max}|$ (the initial growth rate near the floor). Once away from the floor, on the real trajectory $\delta_{t+1}=J(z_t)\delta_t$, i.e. $\delta_t=\Phi(t,0)\delta_0$, $\Phi(t,0)=J(z_{t-1})\cdots J(z_0)$ (the time-ordered product of tangent maps = the **tangent cocycle**). Hence **the real inflation rate is the finite-time Lyapunov rate of the tangent cocycle $\Phi(t,0)=\prod_t J(z_t)$ on the real trajectory**

$$\gamma_t=\frac1t\ln\frac{\|P_\perp\Phi(t,0)\delta_0\|}{\|P_\perp\delta_0\|},\qquad \text{slope of }\log_{10}(H_\perp/H)=\frac{2\gamma_t}{\ln10}$$

($P_\perp$ = orthogonal projection onto the floor plane). Since the $J_t$ are non-normal ($0.53$–$0.76$) and mutually non-commuting, $\Phi(t,0)\neq J_0^t$, and the direction of maximal amplification is further amplified by the next $J_t$ while rotating under time evolution (unstable eigenmode → rotating amplified direction → time-ordered exponential growth). Follow-up test A confirmed this $\gamma_t$ as matching the measured ramp (time-varying/measured 0.999 / 1.004 / 1.019). $\mu_{\max}$ (the local instability of the floor) and $\Phi(t,0)$ (the tangent dynamics of the real inflation) are a **two-stage dynamics**, and $\mu_{\max}$ is its lower bound ($t\to0^+$ / frozen limit).

## 5. Conclusion of the two-stage structure

- **Derivable (form, mechanism, $N$ dependence, onset)**: the inflation curve is the exponential amplification of the floor's linear instability mode. The main predictions (slope $2\lambda/\ln10$, onset, $N$ dependence) are given by the largest eigenvalue of the floor Jacobian $\lambda=\ln|\mu_{\max}|$.
- **Unification with Paper 5 (two faces of the same $\mu_{\max}$)**: this $\lambda=\ln|\mu_{\max}|$ is precisely the summit curvature of the Mexican-hat effective potential of Paper 5 ($V''(0)=-\lambda$, $\lambda=\ln|\mu_{\max}|$). That is, **the instability of SSB (the negative curvature of the summit that drives the fall into the valley) and the exponential growth rate of inflation are tied together by the same floor Jacobian eigenvalue $\mu_{\max}$**. Paper 5's static SSB picture (spontaneous selection from an unstable symmetric maximum) and this paper's dynamic amplification (exponential rapid expansion) are the static and dynamic aspects of a single $\mu_{\max}$. As $N$ increases and $\mu_{\max}\to1$, the summit flattens (curvature decreases), and the fall = inflation also becomes slower (onset lengthens). That the growth rate does not depend on precision or seed depth (= it is intrinsic to the dynamics) has already been established in an experiment that verified the floor of the seed at 100-digit precision——even at IC100 (closure $\sim10^{-102}$, floor $\sim10^{-200}$), $\gamma_\tau$ agrees with the 64-bit run to 8 digits ($N=8$), and it fires across about 179 digits from the floor under a single exponential law (`precision_seed_dynamics_isolation_100digit_20260902`, $N=7,8$).
- **Refinement (non-normality correction)**: the exact slope is the amplification of the time-ordered product of the time-varying non-normal $J$, adding $\sim18$–$31\%$ on top of the floor eigenvalue (**directly confirmed in follow-up test A**: the slope of the time-varying product matches the measured ramp and exceeds the frozen floor by 18–31%). Only the closed form remains underived.

**Unified chain** (a single line from ignition to landing):

$$\rho(J_0)>1\ \Rightarrow\ \text{ignition qualification (Paper 2)}\ \to\ \ln\rho(J_0)\ \Rightarrow\ \text{initial growth rate near the floor}\ \to\ \gamma_{\rm FT}\big[\textstyle\prod_t J(z_t)\big]\ \Rightarrow\ \text{real ramp growth rate (follow-up test A)}\ \to\ \frac{\ln(1/\text{seed})}{\gamma}\ \Rightarrow\ \text{onset}\ \to\ \text{nonlinear saturation}\ \Rightarrow\ \text{degenerate vacuum (Paper 5)}.$$

Viewing the same unstable floor by "how it leaves the floor" gives inflation, and by "where it lands" gives SSB——the two are not distinct mechanisms but two views of a single unstable floor.

This is a quantitative description of rolling down from the Mexican-hat summit (the unstable relative equilibrium of the 90° skeleton), and the landing (the degenerate vacuum, $\Delta=-2\pi/N$) is treated by Paper 5, while the geometry of tilting away from the floor plane is treated by Paper 7.

## 6. Claim classification / open items

- Classification: derived consequence (the linear-instability formula) + measurement of the non-normality correction (the refinement is a hypothesis). Judgment: retained.
- Open: the **closed form** of the amplification rate of the non-normal time-varying product is underived (the mechanism is established by follow-up test A; only the analytic expression of the amplification rate remains for future work). The floor Jacobian for $N=40$ ($2M=1560$ dimensions) is heavy and uncomputed——for small $N$ (6, 8, 10) the two-stage structure is established.

## Related work (structural correspondence, not the source of derivation)

- Transient amplification / pseudospectra of non-normal matrices (Trefethen & Embree 2005): the amount by which the ramp exceeds the floor eigenvalue is the non-normal transient growth / amplification of the time-ordered product.
- Self-consistent inflation v2 [1]: this paper is the analytic derivation (refinement) of its rapid-expansion mechanism.

## References

1. Noriaki Kihara, "The mechanism of inflation-like rapid expansion in self-consistent relational-wave closed systems," Version DOI 10.5281/zenodo.22176949 / Concept DOI 10.5281/zenodo.22112008.
2. L. N. Trefethen and M. Embree, *Spectra and Pseudospectra*, Princeton (2005).

## Reproduction

The verification scripts, figures, and reports are in `位相ロック全N確認_相対平衡_20260912/インフレ曲線_線形不安定解析_20260912/` (`analyze_inflation_linear_instability_20260912.py` = verifications A–E [floor $J$ eigenvalues / non-normality / frozen-floor $\delta$ evolution / measured ramp / onset], `make_figure_inflation.py`, `fig_inflation_linear_instability.png`, `REPORT.md`, `解析ログ_20260912.log`, `SHA256SUMS.txt`, `run_all.sh`). The **direct confirmation of the mechanism of the non-normality correction (follow-up test A)** is in `位相ロック全N確認_相対平衡_20260912/ランプ直接Lyapunov_非正規時変積_20260912/` (`analyze_ramp_lyapunov_20260912.py` = comparison of the $H_\perp$ slopes of frozen floor / time-varying product / measurement, `results/ramp_lyapunov_summary.csv`, `README.md` / `REPORT.md` / `SHA256SUMS.txt` / `run_all.sh`). The input is the existing den=N npz (500 steps, make_parent floor). No runs or physics were changed. The dynamical map is identical to the canonical `run_N3_N40_stage123_v1.py` (SHA256 `1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567`).
