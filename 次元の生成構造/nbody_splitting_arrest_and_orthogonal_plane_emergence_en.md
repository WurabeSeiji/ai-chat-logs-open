# Paper 6
# Arrest of Spontaneous Splitting and Emergence of New Orthogonal Rotation Planes

**Noriaki Kihara**  
**2026-07-24**  
**DOI:** 10.5281/zenodo.21543071 (Concept DOI: 10.5281/zenodo.21543070)

## Abstract

In earlier papers we started from the two simple axioms

$$
\sum_n x_n^2=0,
\qquad
U^n=I
$$

and showed that, even in a closed system that assumes no background spacetime, three directions readable as spatial directions and one direction readable as a conserved quantity arise spontaneously [1,2].

We further showed that a wave initially present as a single mode is amplified geometrically from a seed too small to measure, and spontaneously splits into several waves. This amplification is an inflation-like exponential expansion, but it does not continue indefinitely; it moves into a metastable state of finite amplitude [3].

It remained open, however, why the spontaneous splitting stops, and what happens inside the system at the transition to the metastable state.

In this paper we measure the arrest process directly. We find that, at the moment spontaneous splitting stops and the system moves into the metastable state, several orthogonal rotational directions that did not exist before, together with a new direction readable as a conserved quantity, arise spontaneously.

The claim of this paper is therefore:

$$
\boxed{
\text{Spontaneous splitting of a single wave stops together with the creation of new orthogonal directions.}
}
$$

The waves produced by splitting are not accumulated into the single existing rotational direction. At the arrest point of splitting, new mutually orthogonal rotational directions stand up, and the conserved quantity is distributed into those directions. The metastable state is not a state where splitting has merely weakened; it is the state after the system has acquired new directions.

With this result, the change observed in the earlier paper,

$$
\text{exponential spontaneous splitting}
\longrightarrow
\text{metastable state},
$$

can be described as

$$
\boxed{
\text{localization in a single direction}
\longrightarrow
\text{spontaneous splitting}
\longrightarrow
\text{creation of new orthogonal directions}
\longrightarrow
\text{metastable state with multiple directions}
}.
$$

In the metastable state, the amplitude per relational wave follows

$$
A_{\mathrm{rel}}
\simeq
\frac{1}{\sqrt{M}}
=
\sqrt{\frac{2}{N(N-1)}}.
$$

This means that as the system size $N$ grows, the amplitude of each relational wave decreases in proportion to $1/N$. The reason the metastable amplitude appears dual to $N$ is that the conserved quantity is distributed over

$$
M=\frac{N(N-1)}{2}
$$

relational waves.

The creation of new orthogonal rotational directions, its simultaneity with the arrest of splitting, the inflow of the conserved quantity, and the $1/\sqrt{M}$ law of the amplitude presented here are results obtained from the dynamics and numerical experiments of the model; they are not hypotheses.

What is currently undetermined is what the newly arising directions correspond to in actual physics. They might be identifiable as a time direction, charge directions, or other internal degrees of freedom, but that physical identification is not part of the claim of this paper. Here we show only the existence and creation mechanism of the new directions, and leave their naming and physical meaning to future work.

---

## 1. System and observables

### 1.1 State space

We represent the complete pairwise relations of $N$ vertices as edges, with edge number

$$
M=\binom{N}{2}=\frac{N(N-1)}{2},
$$

and the state as

$$
Z\in\mathbb C^M.
$$

We use the two quadratic quantities

$$
Z^\dagger Z=1,
\qquad
Q(Z)=Z^T Z=0.
$$

Writing $Z=X+iY$, the zero-square closure is equivalent to

$$
\|X\|=\|Y\|,
\qquad
X\cdot Y=0,
$$

so the state is represented by orthogonal, equal-length components in real $M$-dimensional space.

### 1.2 Generator

From the edge phases

$$
\theta=\arg Z,
$$

we construct, via the low-rank vertex decomposition derived earlier [1], the real antisymmetric generator

$$
K(\theta)=WJW^T,
\qquad
K^T=-K.
$$

Time evolution is the Cayley transform of the normalized generator $\widetilde K=K/\sigma_1$,

$$
Z_{\tau+1}
=
\left(I-\gamma\widetilde K_\tau\right)^{-1}
\left(I+\gamma\widetilde K_\tau\right)Z_\tau,
\qquad
\gamma=\tan\frac{\pi}{144}.
$$

Each Cayley matrix is a real orthogonal matrix, so it preserves $Z^\dagger Z$ and $Z^T Z$. In the numerical runs the deviation of the zero-square closure was at most about $10^{-14}$.

### 1.3 Parent state and dormant seed

The initial parent state $v$ is constructed, following the plane decomposition of the fixed-generator system [2], as the dominant eigenmode of the generator built from its own phases,

$$
K(\arg v)v=-i\sigma_1 v,
$$

so the parent is localized in a single real two-dimensional rotation plane

$$
P_1=\operatorname{span}\{\Re v,\Im v\}.
$$

To this parent we add a small zero-closure seed $g$ lying in the kernel of the initial generator:

$$
Z_0=\frac{v+\delta g}{\|v+\delta g\|},
\qquad
K(\arg v)g=0,
\qquad
\delta=10^{-15}.
$$

### 1.4 Splitting quantity

We define the projection ratio onto the parent plane,

$$
P_1(\tau)
=
\frac{|p^T Z_\tau|^2+|q^T Z_\tau|^2}{Z_\tau^\dagger Z_\tau},
$$

and the fraction outside the parent plane,

$$
f(\tau)=1-P_1(\tau),
$$

where $p,q$ are an orthonormal basis of the parent plane.

---

## 2. Onset mechanism of splitting

### 2.1 State-dependent generator

Under the fixed generator $K_0=K(\arg v)$, the kernel component $g$ satisfies $K_0g=0$ and does not grow.

In the original dynamics, however, the generator $K_\tau=K(\arg Z_\tau)$ is reconstructed at each step. When the small seed changes the state phases, the generator also changes, and a component that was in the kernel initially ceases to be in the kernel of the next generator. This feedback amplifies the out-of-parent-plane component.

### 2.2 Generator-freezing experiment

For $N=6$ we compared two runs with the same parent state and same seed:

- fixed generator: $K_\tau=K_0$
- updated generator: $K_\tau=K(\arg Z_\tau)$

With the fixed generator, $f$ stayed at about $10^{-27}$ and did not grow for 2000 steps. With the updated generator, an amplification of about 30 orders of magnitude arose from the same initial seed, reaching $f\simeq0.70$.

Thus the amplification source of splitting is not the action of a fixed generator, but the feedback between the state phases and the generator.

In the initial amplification regime $f(\tau)\simeq f_0 e^{r\tau}$ holds; for $N=6$ we obtained $r\simeq0.044$, consistent with $r\simeq2\lambda_\perp$ where $\lambda_\perp$ is the transverse perturbation growth rate.

### 2.3 Conclusion

We call this mechanism **self-excited amplification by phase feedback**. Since the total norm is conserved, the amplification is not the creation of norm but a redistribution from the parent plane to its orthogonal complement.

---

## 3. Arrest of splitting and metastabilization

### 3.1 Not a simple depletion of parent norm

Let $f_{\mathrm{sat}}$ denote the metastable value of the splitting quantity. Representative values:

| N | f_sat | parent-plane residual 1−f_sat |
|---:|---:|---:|
| 5 | 0.80 | 0.20 |
| 40 | 0.20 | 0.80 |
| 300 | 0.083 | 0.917 |

For $N=40$ and $N=300$, the splitting quantity saturates while most of the norm remains in the parent plane. Therefore the arrest cannot be explained by depletion of the parent-plane norm alone.

What is established at present is that

$$
\boxed{f_{\mathrm{sat}}(N)\text{ decreases as }N\text{ increases}}
$$

and that this arrest is not a simple resource depletion. The closed form and geometric origin of $f_{\mathrm{sat}}(N)$ are not yet derived.

### 3.2 The metastable state is not a fixed point

After saturation the state continues to change in time. For $N=5$–$7$, the transverse growth rate in the metastable regime remained about 80% of its initial value. Hence metastabilization is not convergence to a stable fixed point.

By the metastable state here we mean a bounded non-stationary state with

- $Z^\dagger Z=1$
- $Z^T Z=0$
- $f$ staying in a finite interval
- the state itself continuing to fluctuate.

A rigorous test of "mixing" would require additional measurements such as the Lyapunov spectrum, autocorrelation, and recurrence statistics.

### 3.3 Spectral indicator

Let the nonzero rotation rates of the generator be $\sigma_1\ge\sigma_2\ge\cdots>0$. The second-rate ratio $\rho_2=\sigma_2/\sigma_1$ is fixed at the parent value during amplification and changes at the time the splitting quantity moves into the saturation regime.

For large $N$,

$$
\varepsilon=\frac12-\frac{\sigma_2}{\sigma_1}
$$

decreases like $1/N$. Figure 1 plots $2(N-1)\varepsilon$ and shows the temporal coincidence of the transition from the amplification regime to the metastable regime with the spectral change.

![Overlay of f and σ2/σ1 deviation](第6論文_figures_v1/cause_overlay_f_sigma_ratio_v1.png)

*Figure 1: Splitting quantity $f(\tau)$ and $2(N-1)(1/2-\sigma_2/\sigma_1)$ for $N=5,40,300$. The dotted line marks the time when the splitting quantity moves into the saturation regime.*

Whether this indicator is the cause of the transition or a spectral response accompanying it is undetermined. Here we treat it as an observable that marks the transition time.

---

## 4. $N$-dependence of the metastable amplitude

### 4.1 Representative amplitude of relational components

The representative amplitude per relation in the metastable regime follows

$$
A_{\mathrm{rel}}
\simeq
\frac{1}{\sqrt M}
=
\sqrt{\frac{2}{N(N-1)}}
\sim
\frac{\sqrt2}{N}.
$$

For $N=3$–$300$, 18 values of $N$, 38 runs in total, the mean of $A_{\mathrm{rel}}\sqrt M$ was $1.00004$. A power fit for $N\ge20$ gave exponent $\alpha=1.00820$, $R^2=0.999990$, consistent with the effective exponent $1.00823$ of the exact binomial evaluated on the same $N$ grid. For $N=300$, $A_{\mathrm{rel}}=0.00472237$ against $1/\sqrt{44850}=0.00472192$.

### 4.2 Convention and dynamics

The $1/\sqrt M$ law has two elements. First, from the convention $\sum_{e=1}^{M}|Z_e|^2=1$, perfect equipartition gives each component amplitude $1/\sqrt M$. Second, the actual dynamics approaches $\mathrm{PR}/M\to1$ in the metastable regime and narrows the amplitude distribution in the edge basis. The nontrivial result is not the normalization itself, but that the dynamics builds a near-equipartitioned state in the edge basis.

Hence the observed $1/N$ is not a fundamental law but the large-$N$ limit of $1/\sqrt M$ with $M=N(N-1)/2$.

### 4.3 Relation to the spectral ratio

$$
\frac12-\frac{\sigma_2}{\sigma_1}=O\!\left(\frac1N\right)
$$

is also observed. That the $1/\sqrt M$ law of the edge amplitude and the $1/N$ law of the spectral ratio arise from the same mechanism has not yet been shown; their relation is a task for analysis of the parent-generator spectrum.

---

## 5. Inflow into orthogonal rotation planes

### 5.1 Canonical decomposition of the antisymmetric generator

A real antisymmetric matrix $K$ decomposes canonically into mutually orthogonal two-dimensional rotation planes and a kernel [5]. In a suitable orthogonal basis,

$$
K
\sim
\bigoplus_{j=1}^{r}
\begin{pmatrix}
0 & \sigma_j\\
-\sigma_j & 0
\end{pmatrix}
\oplus 0.
$$

Thus the basic unit of motion in this system is not a single axis but a two-dimensional rotation plane.

For $N=5$, the dominant-plane rotation rate was $\sigma_1=3.742$ and the next nonzero rate was $\sigma_2=\sqrt2$. The cosine of the principal angle between the two subspaces was zero within numerical precision, confirming orthogonality.

### 5.2 The boundary between $N=4$ and $N=5$

We compare the number of nonzero rotation planes of the parent generator with the splitting outcome.

| N | second nonzero rate σ₂/σ₁ | splitting outcome |
|---:|---:|---|
| 3 | 0.000 | no sustained expansion |
| 4 | 0.000 | no sustained expansion |
| 5 | 0.378 | expansion |
| 6 | 0.386 | expansion |
| 7 | 0.407 | expansion |

For $N\le4$ the parent generator has no second nonzero rotation subspace; it first appears at $N\ge5$. This boundary coincides with the boundary of sustained splitting obtained in the earlier paper [3].

This result indicates that a nonzero rotation subspace outside the parent plane may be necessary for sustained splitting. That $\sigma_2>0$ alone is a sufficient condition is not proven.

### 5.3 Measurement of per-plane norm

For the eigen-subspaces of the parent generator $K(\arg v)$, we projected the state norm at each time onto three groups:

1. dominant plane $P_1$: the 2D plane corresponding to $\sigma_1$
2. other rotation subspaces: $0<\sigma<\sigma_1$
3. kernel: $\sigma=0$

From the definition of the splitting quantity via the parent-plane basis, the identity

$$
\boxed{f(\tau)=1-E_{P_1}(\tau)}
$$

holds, where $E_{P_1}$ is the projection norm ratio onto the dominant plane. Thus $f$ is not an abstract "new wave quantity" but the fraction of conserved norm that has flowed out of the dominant plane into its orthogonal complement.

### 5.4 Dense-matrix decomposition: $N=5,40$

For $N=5$ and $N=40$ we constructed the dense matrix $K$ directly and obtained the rotation subspaces and kernel from its eigendecomposition. A clear spectral gap separated the nonzero eigenvalue group from the numerical-zero group.

Projection dimensions:

| N | M | dominant plane | other rotation subspaces | kernel |
|---:|---:|---:|---:|---:|
| 5 | 10 | 2 | 4 | 4 |
| 40 | 780 | 2 | 78 | 700 |

The maximum deviation of the identity $f=1-E_{P_1}$ was $3.7\times10^{-14}$ for $N=5$ and $3.8\times10^{-15}$ for $N=40$. The recorded $f$ agreed bit-for-bit with the original run.

![N=5 plane inflow](第6論文_figures_v1/planeflow_N00005_exact.png)

*Figure 2: $N=5$. Projection norm ratios onto the dominant plane, other rotation subspaces, and kernel. Outflow from the dominant plane into the other two groups begins near $\tau=1167$. Representative metastable ratios are about 0.20, 0.33, 0.47.*

![N=40 plane inflow](第6論文_figures_v1/planeflow_N00040_exact.png)

*Figure 3: $N=40$. Outflow begins near $\tau=2011$; in the metastable regime the ratios are about 0.80 (dominant plane), 0.07 (other rotation subspaces), 0.13 (kernel).*

For both $N=5$ and $N=40$, the onset of the splitting quantity and the norm outflow out of the dominant plane occur at the same time.

### 5.5 Low-rank decomposition: $N=300$

For $N=300$, $M=44850$, and dense-matrix eigendecomposition is not feasible in this computing environment. We therefore used the low-rank structure $K=WJW^T$ and lifted the eigenpairs of $JG$ to the edge space to construct the projection subspaces.

To exclude numerical-zero eigenvalues, we set

$$
\sigma>10^{-6}\sigma_{\max}
$$

as the criterion for a nonzero rotation rate. Compared against the dense-matrix method at $N=40$, this approximation gave matching subspace dimensions, with a maximum deviation of $2.0\times10^{-15}$ over all times and all projection groups.

For $N=300$:

| group | dimension |
|---|---:|
| dominant plane | 2 |
| other rotation subspaces | 598 |
| kernel | 44250 |

![N=300 plane inflow](第6論文_figures_v1/planeflow_N00300_approx.png)

*Figure 4: $N=300$. Low-rank method, nonzero threshold $10^{-6}\sigma_{\max}$. Outflow begins near $\tau=4844$; in the metastable regime the ratios are about 0.91 (dominant plane), 0.04 (other rotation subspaces), 0.05 (kernel).*

### 5.6 Basis dependence

In the edge basis the metastable state shows a near-equal-amplitude distribution close to $1/\sqrt M$. In the projection onto the generator eigen-subspaces, the norm per dimension is not equipartitioned. For example, at $N=40$ the per-dimension representative values are

$$
P_1:0.398,
\qquad
\text{other rotation}:0.0009,
\qquad
\text{kernel}:0.0002,
$$

so the dominant plane is strongly occupied. Therefore

- component equipartition in the edge basis
- norm distribution in the generator eigen-subspaces

are different properties and cannot be identified.

### 5.7 Conclusion

From these observations, splitting can be described as the process

$$
\boxed{
\text{localization in a single dominant plane}
\rightarrow
\text{transverse amplification by phase feedback}
\rightarrow
\text{norm outflow into other rotation subspaces and the kernel}
\rightarrow
\text{bounded metastable distribution}
}.
$$

Whether the newly occupied subspaces correspond to real time, charge, or internal degrees of freedom is not identified here.

---

## 6. Physical identification of the unnamed directions, and related work

What this paper has shown is that orthogonal rotation subspaces and a kernel distinct from the dominant plane exist, and that at splitting the conserved norm flows into them. This existence, orthogonality, onset of occupation, and norm inflow are derived and measured within the model; they are not hypotheses.

Reading these unnamed directions as a time axis, internal charge axes, or other physical quantities is a separate problem. Such identification requires at least:

1. a readout map that selects a fixed number of physical quantities from the $N$-dependent group of rotation planes,
2. a conservation or transformation law for each readout quantity,
3. an operational correspondence between the phase advance of the dominant plane and time measurement,
4. a correspondence between the internal rotational degrees of freedom and the sign/composition rules of charge.

The task of a readout to a fixed number connects to the earlier paper [4], which showed that the number of waves is not an absolute quantity within the system but depends on the readout resolution.

The isotropic two-planes obtained from the zero-square closure can be compared with totally isotropic subspaces in pure-spinor and twistor theory [6]. However, since the rotation planes here are constructed from the time evolution of a state-dependent generator, this comparison indicates a similarity of mathematical structure, not an identity.

As an example of emerging a time direction from rotational symmetry, there is Walker's model obtaining an effective time axis from a gauged rotational symmetry in four-dimensional Euclidean space [7]. Constructions linking spacetime, internal symmetry, and Standard-Model representations via $D_4$, Clifford algebra, octonions, and $E_8$ have also been proposed [8–11]. These are comparison targets for the physical identification we will examine in the future, but are not used in this paper's proof of splitting arrest and new-direction creation.

In a separate document we examine the hypothesis of reading

$$
(x,y,z,R)\oplus(Q_1,Q_2,Q_3,t)
$$

as two $D_4$ structures. But this identification is not derived from the numerical results of this paper alone. Here we treat the new directions as unnamed and separate their physical naming to subsequent research.

---

## 7. Conclusion

We measured the onset, saturation, and inflow destination of spontaneous splitting in an $N$-body relational-wave closed system.

First, fixing the generator to the parent state produces no splitting; exponential amplification occurs only when the generator is reconstructed at each step from the state phases. The onset mechanism of splitting is the phase feedback between the state and the generator.

Second, for large $N$ the splitting quantity saturates while 80–92% of the norm remains in the parent plane. Hence the arrest is not a simple depletion of parent norm. The mechanism determining the $N$-dependent saturation value is not derived.

Third, the splitting quantity equals the fraction of norm that has flowed out of the dominant plane, and the destination was the other nonzero rotation subspaces and the generator kernel. Thus the spontaneous splitting of this system can be formulated as a redistribution of conserved norm among orthogonal subspaces.

Fourth, the metastable edge-component amplitude followed $A_{\mathrm{rel}}\simeq1/\sqrt M$, explained as the combination of the total-norm convention and dynamical equipartition in the edge basis.

Therefore the single-wave state of this system is transversely unstable under the state-dependent generator and moves conserved quantity into new orthogonal subspaces. Metastabilization is not a cessation of wave production but a bounded redistribution state that continues under a finite outflow.

---

## 8. Established results, hypotheses, and underived items

### 8.1 Established

- The Cayley update preserves $Z^\dagger Z$ and $Z^T Z$.
- Under a fixed generator the kernel seed does not amplify; under a state-dependent generator it does.
- The initial splitting quantity grows exponentially.
- The large-$N$ saturation is not explained by depletion of parent-plane norm.
- The metastable state is not a fixed point but a bounded state that keeps varying in time.
- The edge-basis representative amplitude follows $1/\sqrt M$.
- For $N\le4$ there is no second nonzero rotation plane; it appears at $N\ge5$.
- $f=1-E_{P_1}$ holds; the splitting quantity is the outflow from the dominant plane.
- For $N=5,40$, dense-matrix eigendecomposition confirms the outflow destination is the other rotation subspaces and the kernel.
- The $N=300$ low-rank decomposition reproduces the dense-matrix method at $N=40$ to machine precision.

### 8.2 Observed but unproven relations

- $f_{\mathrm{sat}}(N)$ decreases with $N$.
- $\sigma_2/\sigma_1\to1/2$ appears to hold.
- $1/2-\sigma_2/\sigma_1=O(1/N)$ appears to hold.
- The existence of a second nonzero rotation plane may be a necessary condition for sustained splitting.

### 8.3 Underived

- Closed form of the saturation value $f_{\mathrm{sat}}(N)$
- Condition under which the outflow stays bounded in the metastable regime
- Analytic derivation of $\sigma_2/\sigma_1\to1/2$
- Per-plane inflow within the other rotation subspaces
- Dynamical role of the component flowing into the kernel
- Physical readout of the unnamed rotation subspaces

---

## 9. References

### Self-citations

[1] N. Kihara, "Linear upper bound of the generator rank and three-direction saturation in N-body complete pairwise relational waves," Zenodo, Concept DOI: 10.5281/zenodo.21465898, 2026.

[2] N. Kihara, "Plane-decomposition readout in N-body fixed-generator systems," Zenodo, Concept DOI: 10.5281/zenodo.21468959, 2026.

[3] N. Kihara, "Onset and threefold classification of outcomes of spontaneous splitting in N-body relational-wave closed systems," Zenodo, Concept DOI: 10.5281/zenodo.21486233, 2026.

[4] N. Kihara, "The number of waves is the resolution of the system," Zenodo, Concept DOI: 10.5281/zenodo.21486544, 2026.

### External references

[5] R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013. Reference for the orthogonal canonical form and 2D rotation-block decomposition of real antisymmetric matrices.

[6] A. Taghavi-Chabert, "Twistor Geometry of Null Foliations in Complex Euclidean Space," *SIGMA* 13, 005 (2017). DOI: 10.3842/SIGMA.2017.005; arXiv:1505.06938. Structural comparison with totally isotropic subspaces, pure spinors, and twistor space.

[7] M. L. Walker, "Spontaneous Emergence of a Causal Time Axis in Euclidean Space from a Gauged Rotational Symmetry Theory," *Symmetry* 16(1), 4 (2024). DOI: 10.3390/sym16010004. Prior example of emerging a time direction from rotational symmetry.

[8] F. D. T. Smith, "Higgs and Fermions in D4-D5-E6 Model based on Cl(0,8) Clifford Algebra," arXiv:hep-th/9403007 (1994). Prior proposal connecting $D_4$, Clifford algebra, and eight-dimensional structure to physical representations.

[9] N. Furey, "Three Generations, Two Unbroken Gauge Symmetries, and One Eight-Dimensional Algebra," *Physics Letters B* 785, 84–89 (2018); arXiv:1910.08395. Construction of three generations and $SU(3)\times U(1)$ from complex octonions.

[10] I. Todorov, "Octonion Internal Space Algebra for the Standard Model," *Universe* 9(5), 222 (2023). DOI: 10.3390/universe9050222; arXiv:2206.06912. Organization of the Standard-Model internal space via octonions and Clifford algebra.

[11] C. A. Manogue, T. Dray, and R. A. Wilson, "Octions: An $E_8$ Description of the Standard Model," *Journal of Mathematical Physics* 63, 081703 (2022). DOI: 10.1063/5.0095484; arXiv:2204.05310. Comparison target for Standard-Model and Lorentz-algebra representations within $E_8$.

---

## 10. Reproducibility

### 10.1 Original programs

- `run_n_scaling_lowrank_v1.py`
- `run_spontaneous_splitting_largeN_v1.py`

The originals are fixed by SHA-256, and the instrumented programs do not modify the time-evolution loop. In all instrumented runs, the splitting quantity $f$ agreed bit-for-bit with the original output.

### 10.2 Instrumented programs

- `run_metastable_series_v1.py`
- `run_cause_instrumented_v1.py`
- `make_cause_overlay_figure_v1.py`
- `run_plane_flow_exact_v1.py`
- `run_plane_flow_approx_v1.py`
- `make_metastable_series_figure_v1.py`
- `trace_parent_iteration_v1.py`

The nonzero-detection threshold of the low-rank approximation is recorded in CLI arguments, output JSON, and the figures. Its default is $\texttt{sigma\_rel\_threshold}=10^{-6}$.

### 10.3 Output data

- `cause_instrumented_result_v1/`
- `plane_flow_result_v1/`
- `metastable_series_result_v1/`
- `関係波準安定振幅_N300定量結果_v1.md`
- `親結晶_低N観察記録_v1.md`

CSV files are regenerable outputs and are kept out of repository tracking.

---

## 11. Next experiments

1. **Seed scan.** Scan $\delta=10^{-9}$ to $10^{-15}$ and measure the logarithmic shift of the saturation time and the independence of the saturation value.

2. **Transverse growth rate at large $N$.** Measure the finite-time transverse growth rate simultaneously in runs including $N=40, 300$.

3. **Inflow into individual rotation planes.** Decompose the "other rotation subspaces" into the individual $\sigma_j$ planes and examine whether the inflow concentrates in a few planes or disperses over many.

4. **Saturation law.** Scan $N$ densely to obtain the functional form of $f_{\mathrm{sat}}(N)$.

5. **Role of the kernel component.** Measure the inflow into the kernel, its residence time, and the exchange rate with the rotation subspaces.

The physical identification and the correspondence to $D_4$ structures and $E_8$ diagonal gluing are separated from these dynamical measurements and treated in subsequent papers.
