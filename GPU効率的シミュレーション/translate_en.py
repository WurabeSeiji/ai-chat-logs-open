import codecs

english_text = r"""# GPU-Efficient Simulation Technique of Test Particle Motion in Closed Positive Curvature Spacetime (Tensor Processing Version)

## Abstract
As a conceptual model to understand the "falling in curved spacetime" of General Relativity and to process it at high speed on computers, we present plots depicting space and particle trajectories (Figure 1). In this simulation verification, an educational analogy model (Toy Model) simulating the real environment of closed positive curvature spacetime is adopted to profile the compatibility with GPU architectures and processing limits.

![Visualization of Geodesics](tensor_lookup_geodesic.png)
*Figure 1: The spatial model targeted by this simulation. The squares of the flat "index space" are bent into a "wedge mesh" that narrows as it goes deeper toward the orthogonal curvature axis ($W$ axis) (with radius of curvature $R=100$). Test particles A and B simply go straight along this grid without any gravitational acceleration logic applied, yet they curve and fall toward the origin in the physical coordinate space due to the inherent geometric contraction of the grid. This plot demonstrates that the precisely calculated exact trajectory (semi-transparent line) and the tensor interpolation trajectory (dotted line) in this proposed method perfectly coincide against the background space mesh.*

Subsequently, Figure 2 depicts dynamic topological changes that occur when "time" progresses on this curved surface, drawing a 3D plot adding the time axis (-t axis).

![Dynamic Tensor Space Model: Frustum](mesh_distortion_frustum.png)
*Figure 2: Tensor spatiotemporal analogy used for computations. The bottom surface acts as the flat index space (Special Relativity domain), and as time progresses upwards (-t axis), the whole geometry contracts due to the influence of curvature, distorting into a "frustum".*

To simulate objects falling in a gravitational field, it is usually necessary to repeatedly calculate acceleration step by step. The proposed method pre-constructs geometric data where "space contracts toward gravity as time progresses" (as in Figs 1 and 2) ahead of time into a 4-dimensional array (Tensor $\mathcal{T}$) in Step 1. Within this tensor space, since the gravitational influence is persistently preserved as spatial distortion, geodesic trajectories can be calculated dynamically during execution by relying solely on linear interpolation (Lerp) without calculating any acceleration parameters (Step 2).

## 1. Purpose
Numerical simulations of test particle motion (geodesics) in curved spacetime based on General Relativity typically require nonlinear calculations of Christoffel symbols and gravitational acceleration per step. Ultimately, this leads to massive computational bottlenecks in N-body macroscopic physical problems. This paper proposes a highly-parallelized GPU simulation method that completely abstracts this calculation into an external pre-computed index data structure, replacing the heavy nonlinear simulation loop with trivial linear tracking operations in a local inertial frame.

## 2. Theoretical Background: Tensorization Based on Equivalence Principles
According to Einstein's Equivalence Principle, in a sufficiently local region, gravity effectively vanishes, permitting the assumption of an uncurved local inertial frame (Special Relativity domain). This method computationally maps this strictly local flatness onto an orthogonal, uniform, and linear memory mesh called the index memory space.

Specifically, a $n^4$ integer index space is prepared where intervals of time and space are fixed static constants. Although it serves as a pure cubic mesh algorithmically in the index parameter space, mapping it to a physical embedding space defined by curvature $R$ distorts the array's absolute topology into a "frustum" shape framework. This fact—that "uniform equal intervals in index space inherently embody distorted intervals in physical coordinate space"—serves as the direct computational representation of the metric tensor.

## 2.1 Theoretical Relationship Between Mapping Function $\Psi(\boldsymbol{\xi})$ and Metric Tensor $g_{\mu\nu}$
The essence of the nonlinear mapping function $\Psi: \mathbb{R}^4 \to \mathbb{R}^4$—transforming the abstract index space $\boldsymbol{\xi}$ to the physical embedding space $\mathbf{X}$—is strictly equivalent to coordinate transformations (Pushforward) in differential geometry.
Assume the infinitesimal displacement $d\boldsymbol{\xi}$ in the abstract index space trivially follows a locally flat Minkowski metric $\eta_{\alpha\beta}$ inside the proposed tensor matrix. In this scenario, the metric tensor $g_{\mu\nu}$ in the physical tracking space is strictly linked and derived via the Jacobian matrix $J^\alpha_\mu = \partial \Psi^\alpha / \partial \xi^\mu$ of the mapping function $\Psi$, expressed as follows:

$$ g_{\mu\nu} = \frac{\partial \Psi^\alpha}{\partial \xi^\mu} \frac{\partial \Psi^\beta}{\partial \xi^\nu} \eta_{\alpha\beta} $$

Thus, the discrete gradient differences between each node element (architectural spatial distortion) of the pre-computed tensor array $\mathcal{T}$ completely encapsulate geometric $g_{\mu\nu}$ information fundamentally required by the Einstein field equations without compromise. Execution-time particles functionally execute bare uniform rectilinear motion $\frac{d^2 \xi^\mu}{d\tau^2} = 0$ in the index architecture (Step 2); nonetheless, observing this displacement from the physical space computationally guarantees identical curved mathematical motion parallel to the highly-nonlinear geodesic equations dictated by the Christoffel symbols $\Gamma^\mu_{\alpha\beta}$ via Chain Rule mechanics. Consequently, this demonstrates the absolute mathematical rigor of completely isolating and statically pre-baking active non-linear metric processing exclusively into raw "Memory Matrix Mapping $\mathcal{T}$" prior to tracking execution.

## 3. Proposed Algorithm (2-Step Method)
This architecture is composed of two distinct computational stages:

### Step 1: Pre-computation and Build of the Rigorous Tensor Space ($n^4$ mesh)
Before simulation begins, an analytic mapping function $\Psi(\boldsymbol{\xi})$ from indexing coordinates $\boldsymbol{\xi}$ to physical coordinates $\mathbf{X}$ is calculated.

Geodesic trajectories in standard General Relativity fundamentally evaluate the differential geodesic equation:
$$ \frac{d^2 X^\mu}{d\tau^2} + \Gamma^\mu_{\alpha\beta} \frac{dX^\alpha}{d\tau} \frac{dX^\beta}{d\tau} = 0 $$
In an authentic relativistic metric environment, evaluating analytic solutions invoking rigorous elliptic integrals and functions—such as those historically derived by Y. Hagihara (1931)[4] and S. Chandrasekhar (1983)[5]—are inherently required. Although the pre-assembly algorithm presented in this paper (Step 1) is intrinsically versatile enough to tolerate embedding ANY of these intense non-linear geometric solving methods natively during build-time, an educational analogy load model (Toy Model) simulating continuous Lorentz contraction $v = (C \cdot t)/R$ is adopted purely as the functional metric surrogate $\Psi(\boldsymbol{\xi})$ during this verification run. Doing so serves distinctly to isolate and profile memory bandwidth processing limits and algorithmic parallel scaling characteristics uncompromised.

**[Strict Analytic Analogies for Spatial and Temporal Axes $\mathbf{X}^\mu$]**
As an experimental and educational analogy model, nonlinear tracking is configured utilizing an orthogonal contraction ratio factor $\gamma = 1/\sqrt{1 - (t/R)^2}$ during test mass propagation cycles. Coordinate projections are evaluated alongside nonlinear operators mapped into the synthesized parameter space and stored algorithmically into absolute geometry $\mathbf{X}^\mu(\xi^0)$. By doing so, massive FPU float operations are emulated, effectively yielding a durable processing load model accurately reflecting true relativistic physical measurement overheads.

Within Step 1, these nonlinear mathematical computations are executed via standard CPU evaluation algorithms exclusively for all discrete $n^4$ indexing origins $\mathbf{i}$, and comprehensively dumped securely into fixed, pre-computed memory banks identified as tensor variables $\mathcal{T}[\mathbf{i}] = \Psi(\mathbf{i})$.

### Step 2: Parallel and Linear Runtime Processing (Local Special Relativity Loop)
All test masses simply hold active internal index-relative coordinates $\boldsymbol{\xi}^{(s)}$ during execution cycles.
In step $s$, accelerating calculations are strictly ignored. Test particles rigidly perform elementary, linear tracking simulating absolute zero-sum acceleration:

$$
\boldsymbol{\xi}^{(s+1)} = \boldsymbol{\xi}^{(s)} + \mathbf{v}_{\text{local}} \cdot \Delta s
$$

Following standard position propagation, absolute coordinate metrics are physically "looked up" actively from the immutable tensor storage geometry $\mathcal{T}$ via local indexing $\boldsymbol{\xi}^{(s+1)}$ coordinates:

$$
\mathbf{X}_{\text{phys}}^{(s+1)} = \mathcal{T} [\lfloor \boldsymbol{\xi}^{(s+1)} \rfloor ]
$$

Through this dramatically reduced linear routine setup, gravitational tracking phenomena are completely and accurately mirrored.

## 4. Advantages and Trade-offs

### 4.1 Scalability Characteristics
The primary hallmark of the approach distills dynamic test tracking directly into raw "addition" and "basic memory referencing" architectures inside loop executions. Utilizing high-end multi-core instances or standard GPUs actively guarantees a virtual constant runtime profile ($O(1)$ approximation) through completely bypassing mathematical choke points across the algorithm.

### 4.2 Error Artifact Mitigation via Linear Interpolation
Proper test particle mapping functionally translates back to $\mathbf{X}_{\text{exact}}^{(s)} = \Psi(\boldsymbol{\xi}^{(s)})$ formally. Nevertheless, it is necessary to slash live tracking processor requirements significantly. Therefore, discrete tensor referencing parameters govern primary loop constraints. Directly employing integer truncation functions (Floor) natively restricts calculated coordinate points triggering sharp zigzag mapping constraints (Manhattan distance scaling structures).
However, accurately deploying multidimensional Linear Interpolation (Lerp) algorithms efficiently neutralises stepping discrepancies:
$$
\mathbf{X}_{\text{phys}}^{(s)} \approx \mathrm{Lerp}(\mathcal{T}, \boldsymbol{\xi}^{(s)})
$$
Utilizing virtually cost-free hardware-subbed interpolation operators naturally removes step manifestations, practically synthesizing completely flush continuous geometric curve alignments functionally reflecting true spatial curvatures organically.

## 5. Verification: Direct Equivalency Measurements against Interpolated Trajectories

A spacetime array grid framework measured heavily against distortion limits ($R=100$) was strictly compiled as $T[80, 20, 20]$. Figure 1 actively highlights subsequent particle simulations directly verifying algorithmic mapping integrities:

- **Stationary Mass Traces (Red Line Particle A)**: Holds an initial zero propagation momentum $\mathbf{v}=0$. Disregarding index shifting mapping mechanics, simple dynamic compression of the framework pushes the particle securely inside physical limits, yielding a stable drop directly towards metric origins. Interpolated lines identically align to analytically solved constants.
- **Cross-Spatial Pathing Operations (Blue Line Particle B)**: Assumes localized exit velocity properties. Tracking visibly curves matching expected physical curvature influences tightly. Interpolated functions actively correct step constraints rendering identical comparisons cleanly measured against exact geodesic mapping data parameters natively.

By physically validating outcomes, the proposed structure seamlessly decouples complex, mathematically dense General Relativity metrics into pure localized Tensor Arrays + Linear Operation formulas entirely compatible with large array computational execution algorithms structurally.

## 6. Mathematical and Theoretical Citations

The baseline theory encompassing local linear properties coupled with geometric mapping integrations traces fundamentally to validated physical frameworks defined across numerical models actively:

- [1] T. Regge, "General Relativity without coordinates," *Il Nuovo Cimento*, vol. 19, no. 3, pp. 558-571 (1961).  
  *(Discrete structural tracking models detailing coordinate flat models measured identically omitting calculus equations natively).*
- [2] C. W. Misner, K. S. Thorne, and J. A. Wheeler, "Gravitation," *W. H. Freeman* (1973).  
  *(Definitions outlining Special Relativistic functions remaining completely absolute locally despite larger coordinate system curvatures organically).*
- [3] E. F. Taylor and J. A. Wheeler, "Exploring Black Holes: Introduction to General Relativity," *Addison-Wesley* (2000).  
  *(Summarised spatial distortion formulas regarding localized time tracking scaling limits).*
- [4] Y. Hagihara, "\textit{Theory of the relativistic trajectories in a gravitational field of Schwarzschild}," Japanese Journal of Astronomy and Geophysics, Vol. 8, p.67 (1931).
  *(Derivation of space-axis constraints and elliptic integrations resolving analytical models).*
- [5] S. Chandrasekhar, "\textit{The Mathematical Theory of Black Holes}," Clarendon Press, Oxford (1983).
  *(Validating temporal integration paths required identifying exact physical parameters).*

## 7. Accuracy and Speed Verifications

To concretely confirm technique superiority limits universally, empirical performance profiling validating linear structural proposals comparing raw algorithmic calculations was comprehensively executed natively.

### 7.1 Mathematical Deviation Limits (Accuracy Evaluation)

**Condition Limits**: Evaluate topological degradation variances tied algorithmically to coordinate sampling scales compared natively observing spatial geometric radius models $R$. "Truncation (Floor)" processing metrics were simultaneously logged against dynamic "Linear Interpolation (Lerp)" correction parameters objectively.

**Experiment Specifications**:
- Fixed test masses performing predefined tracking cycle offsets.
- Absolute curvature limits evaluated logarithmically starting $10^2$ towards flatter metrics observing scale differences $10^4$ natively measuring endpoint variances strictly against analytic comparisons mathematically.

**Visual Data Output**:

![Accuracy vs Curvature R](accuracy_vs_R.png)
*Figure 3: Graphical logarithmic measurement outlining absolute deviation ratios referencing curvature $R$ limits cleanly.*

**Data Findings**:
Referencing Truncation paths evaluated via Red Line parameters identically models Manhattan stepping issues resulting in permanent quantization constraints inherently trailing limits tracking despite extreme surface flattening measurements universally.
Contrastingly applying continuous N-Dimensional Interpolation (Blue Line properties) successfully curtails mathematical deviation issues, tracking exponential mapping corrections approaching negligible parameters flawlessly once physical topologies scale significantly towards flatter limitations functionally.

### 7.2 Execution Scalability (Speed Evaluation)

**Testing Parameters**: Compare processing load factors required tracing comprehensive point data matrices evaluating conventional rigorous calculation mechanics cleanly against localized Tensor + Lerp approximations natively.

**Evaluation Framework**:
Assess topological rendering timelines projecting geometric grids identically housing coordinate array factors mapped logically into $V = N^3$ metrics linearly.
- **[Method A] Dense Analytic Calculations**: Execute demanding metric measurements targeting dense geometric matrices simulating $100^3$ index nodes projecting strict numerical calculation integrations consistently.
- **[Method B] Proposed Strategy (Coarse Array Build + Dense Lerp Mechanics)**: Profiling restricted uniquely against dual calculation arrays defined systematically:
  1. **Phase 1 (Sparse Data Building)**: Establish isolated evaluation points exclusively structuring geometric tensor frameworks modeling identical coordinates (e.g. projecting isolated calculation parameters measured against $11^3$ ranges strictly mapping framework constants globally).
  2. **Phase 2 (Linear Interpolation Constraints)**: Execute parallel array scanning parameters targeting comprehensive grid indexing using Python Numba JIT engines avoiding strict variable mapping overhead constraints evaluating full calculations identically.
  *Note: Combined execution data variables track metrics compiled across Phase 1 limits matched cleanly alongside Phase 2 operations definitively.*

Furthermore, note functionally that the primary profiling implementation constructed exclusively via the Python Numba JIT wrapper actively executes over CPU resources purely testing processing load simulations accurately. Fundamentally referencing the specific algorithm framework dictating strict "Array Traversals + N-Dimensional Interpolation (Lerp) calculations", future deployment characteristics interface radically alongside native graphics processor Texture Fetch and Tensor Core pipeline architecture limits functionally scaling processing limits matching literal bandwidth barriers explicitly (Upcoming extension tasks scaling optimizations fully).

**Execution Visual Outputs**:

![Speed vs Particle Count N](speed_vs_N.png)
*Figure 4: Logarithmic execution speed testing highlighting scaling factors projecting total geometric framework $V = N^3$ calculations efficiently.*

**Metric Outcomes**:
Tracking Figure 4 outputs conclusively confirms intensive execution overhead matching dense 1,000,000 node iterations triggering 0.686 seconds calculating Method A variables natively via sustained active ALU limits effectively.
Transitioning operational workloads mapping Method B implementations drastically reduces arithmetic requirements simulating mere 1331 point calculation requests consuming 0.0013 seconds natively across Phase 1 matrix construction constraints globally. Adding secondary uniform node calculations executing Lerp tracking parameters (summing an effective 0.004 seconds), overall profiling finishes successfully navigating identical testing loops under 0.005 seconds globally. Comparisons clearly underline an operational **130x Speed Augmentation** mapping precisely against baseline measurements.

## 8. Appendix: Simulated Research Code Profile (Python)

Fully replicable algorithms plotting the structural framework alongside accurate simulation benchmarks validating mathematical boundaries generating all plots identical to Figure 1 through Figure 4 successfully attached beneath explicitly for independent environment evaluations natively.

"""

import codecs

try:
    with codecs.open("閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）.md", "r", "utf-8") as f:
        lines = f.readlines()
        
    code_start = -1
    for i, line in enumerate(lines):
        if line.startswith("```python"):
            code_start = i
            break
            
    if code_start != -1:
        english_text += "".join(lines[code_start:])
        
    with codecs.open("GPU_Efficient_Simulation_Technique_of_Test_Particle_Motion_in_Closed_Positive_Curvature_Spacetime.md", "w", "utf-8") as f:
        f.write(english_text)
    print("English version exported successfully.")
except Exception as e:
    print(e)
