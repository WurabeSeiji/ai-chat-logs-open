# The Sixth Thought Experiment: Can Two Kinds of Long-Range Interaction Be Internalized into a Single Closed State Map?
## - Reconstructing Gravity, Coulomb, GW-Quadrupole and EM-Dipole Radiation of a Charged Quasi-Circular Binary with Eight Persistent States, and Verifying the Transferability of the Same Transition over Five Conditions -

**Author:** Noriaki Kihara  
**ORCID:** 0009-0004-6753-4020  
**Version DOI:** 10.5281/zenodo.22974632  
**Concept DOI:** 10.5281/zenodo.22974631  
**Version:** v1.1  
**Date:** 2026-09-26  

**Keywords:** charged binary, state closure, synchronous state map, Coulomb interaction, gravitational-wave radiation, electric-dipole radiation, radiation reaction, persistent logical state, discrete dynamics, reproducibility

---

# Abstract

The fifth thought experiment [5] self-audited, from the viewpoint of information flow, the gravitational-orbit generator that the fourth thought experiment [4] had displayed as two states, and confirmed that in the present representation the information held across macro steps must be made explicit as

$$
Z_6=(U,P,E,H,Q,N).
$$

Furthermore, the RK4 intermediate quantities and the processing phase were placed inside the state as machine/work states, and the generator was reconstructed into an OLD-state-only synchronous update, with no external physical parameters, a non-feedback readout, and a single `transition(z)` as the only generating path.

In this paper those audit rules are applied for the first time to **a system extended from one kind of interaction to two kinds**. The object is a charged two-body system restricted to $G=M=c=1$, non-spinning point particles, quasi-circular orbits and a leading-order adiabatic inspiral. The conservative motion uses Newtonian gravity and the Coulomb interaction, and the dissipation uses the leading gravitational quadrupole radiation and the leading electric-dipole radiation. On the independent analytic reference side, the charge dependence can be organized into the two relational quantities

$$
C=1-\lambda_A\lambda_B,
\qquad
D=(\lambda_A-\lambda_B)^2,
$$

where $\lambda_A=q_A/m_A$ and $\lambda_B=q_B/m_B$. $C$ enters the conservative force and the orbital angular velocity, and $D$ decides the presence and the strength of the leading electric-dipole radiation.

Adding $C,D$ to the six persistent states of the fifth thought experiment gives

$$
\boxed{
Z_8=(U,P,E,H,Q,N,C,D).
}
$$

In the implementation, a full microstate of 41 real components is used, consisting of these eight persistent logical states, the RK4 work registers and the 11-phase one-hot phase. `transition(z)` receives no external physical argument, builds every candidate value from the OLD full state only, and generates the next state simultaneously by a one-hot sum of products. Since $N,C,D$ are held as identity candidates in all phases, they are preserved as the result of the same map, not by an external copy.

For the reference case

$$
\lambda_A=+0.30,
\qquad
\lambda_B=-0.30,
$$

that is,

$$
C=1.09,
\qquad
D=0.36,
\qquad
N=\nu=0.25,
$$

the strict eight-state generation from $r=50\rightarrow20$ completed in 488,345 macro steps. The differences from the analytic reference were

$$
|\Delta t|=1.46\times10^{-8},
\qquad
|\Delta\phi|=8.44\times10^{-11},
$$

and over the whole stored range

$$
\max|\Delta r|=9.12\times10^{-9},
$$

with

$$
\Delta N=\Delta C=\Delta D=0.
$$

Furthermore, with the same `transition(z)`, the same 11-phase machine and the same eight persistent states unchanged, the five conditions differing only in the initial $C,D$,

$$
(+0.30,-0.30),
(+0.30,+0.30),
(+0.90,-0.90),
(+0.90,+0.90),
(+1.20,-1.20),
$$

were run. In all five conditions the orbital state generation completed, and the drift of $N,C,D$ was 0. In particular, for equal-sign equal-ratio cases $D=0$, so that only the leading electric-dipole term adopted in this paper vanishes, while the Coulomb repulsion effect with $C<1$ remains in the conservative motion. The `n3_repulsive` case with $C=0.19$ required about 3816.7 orbits; $P,H$ could be generated to the end, but the time readout through the exponential clock $Q$ exceeded the finite-precision range. This fact is not corrected; it is recorded as it is, as a limit of the domain and of the numerical representation.

What this paper shows is that, for the restricted leading-order quasi-circular charged binary, the conservative effects of gravity and Coulomb and the dissipative effects of GW quadrupole and EM dipole radiation **can be rearranged into a single synchronous state map with no external charge parameters and no branching by force type**. This does not mean a fundamental unification of gravity and electromagnetism, a derivation of the Einstein--Maxwell equations, the minimality of eight states, the ontological fundamentality of $C,D$, or universality for general charged binaries.

---

# 0. Position of this paper

This paper is the sixth thought experiment under the research project design [0].

The first thought experiment [1] examined, without placing a background spacetime or a metric first, the minimal construction that distinguishes inertial, uniformly accelerated and rotating motion, from two degrees of freedom of unfixed meaning and the conservation of an area readout.

The second thought experiment [2] applied a fixed action repeatedly to two values and constructed Coulomb-type inverse-square motion from a product readout and an area clock.

The third thought experiment [3] reversed the direction and treated the inverse problem of reading out the mass ratio, the relative distance and the relative time from known two-body orbits.

The fourth thought experiment [4] brought the gravitational orbits computed externally in the third thought experiment back to an internal sequential generation by the state-dependent action

$$
X_{n+1}=S_nX_n.
$$

The fifth thought experiment [5] self-audited that fourth thought experiment. What became important there was the point that

$$
\boxed{
\text{even if a numerical result agrees with a known orbit, the system is not necessarily closed by the explicit states alone.}
}
$$

Since information held across macro steps remained outside the apparent two states $X=(a,b)^T$ of the fourth thought experiment, it was made explicit as the six persistent logical states

$$
Z_6=(U,P,E,H,Q,N),
$$

and the intermediate RK4 values and the 11-phase one-hot phase were also placed inside the full state as machine/work states.

The sixth thought experiment does not change these audit rules themselves. The new question is

$$
\boxed{
\begin{array}{c}
\text{when the interaction is increased from one kind, gravity, to two kinds, gravity + Coulomb,}\\
\text{and the dissipation too is increased to the two channels GW quadrupole + EM dipole,}\\
\text{can the system still be closed as a single synchronous state map without external physical information?}
\end{array}
}
$$

The subject of this paper is therefore not merely "computing charged orbits". With the known charged binary placed on the answer-checking side, it is to audit **what must be held internally as persistent states so that the same update rule transfers to several conditions**.

---

# 1. Motivation

## 1.1 From one kind of interaction to two

In a gravity-only binary, there is one kind of central interaction. When the Coulomb interaction is added, on the conservative side the same $1/r^2$ type changes between attraction and repulsion according to the sign. On the dissipative side, moreover, a leading electric-dipole radiation appears, distinct from the gravitational quadrupole radiation [6--11].

The charged binary is therefore the minimal next stage that tests

$$
\boxed{
\text{multiplication of conservative interactions}
+
\text{multiplication of dissipation channels}
}
$$

at the same time.

On the other hand, going to the general Einstein--Maxwell two-body problem makes the object of verification rapidly complex. This paper deliberately restricts itself to

- non-spinning point particles,
- quasi-circular orbits,
- leading-order conservative dynamics,
- leading GW quadrupole radiation,
- leading EM electric-dipole radiation,
- adiabatic energy balance,
- $G=M=c=1$.

Within this range an independent analytic reference can be stated explicitly, and the direction of comparison with the state generator can be fixed to the one direction

$$
\boxed{
\text{state generator output}
\longrightarrow
\text{independent analytic reference}.
}
$$

## 1.2 Not supplying the charge parameters from outside at every step

If a charged binary is implemented in the ordinary way, at each step

$$
\lambda_A,
\qquad
\lambda_B
$$

can be referred to as external parameters. Seen from the audit rules of the fifth thought experiment [5], however, "the physical information that determines the next state" then remains outside the state.

This paper therefore first audits the combinations that actually appear in the analytic formulas. For the present restricted model it turns out that, instead of holding the individual $\lambda_A,\lambda_B$ directly,

$$
\boxed{
C=1-\lambda_A\lambda_B
}
$$

and

$$
\boxed{
D=(\lambda_A-\lambda_B)^2
}
$$

suffice.

$C$ enters the conservative force and the orbital frequency, and $D$ enters the leading electric-dipole radiation. If these two can be internalized into the state, `transition(z)` needs neither the individual labels of bodies A/B, nor an external attraction/repulsion branch, nor a gravity/electromagnetism force-type branch.

## 1.3 Looking at transferability, not agreement in one condition

Even if one condition agrees with the analytic reference, the possibility remains that the coefficients were made for that case only. This paper therefore keeps the same `transition(z)` fixed and changes only the $C,D$ put into the initial state.

In particular, by including

- opposite-sign charges: $C>1$, $D>0$,
- equal-sign equal charge-to-mass ratios: $C<1$, $D=0$,

the strength of the conservative force and the presence of dipole radiation are changed at the same time.

The success criterion is

$$
\boxed{
\text{that the same transition handles different conditions by the difference of state values alone}.
}
$$

---

# 2. The question and the scope of the claims

The question verified in this paper is divided into the following four stages.

1. **Analytic closure**  
   Can the charge dependence of the present leading-order quasi-circular charged binary be organized into $C,D$?

2. **State closure**  
   With $C,D$ added to the six persistent states of the fifth thought experiment,
   $$
   Z_8=(U,P,E,H,Q,N,C,D),
   $$
   can the physical case values needed to determine the next state be placed inside the state?

3. **Single map**  
   Can the full state, including the RK4 intermediate quantities, the processing phase and the conserved relational quantities, be updated by a single `transition(z)` that takes no external physical argument?

4. **Transferability**  
   Can the same `transition(z)`, unchanged, be transferred to five different conditions by changing only the $C,D$ of the initial state?

The affirmative conclusions of this paper are stated **only for this restricted model and this range of implementation**.

---

# 3. Independent analytic reference

## 3.1 Approximations and units

The analytic reference is constructed independently of the state generator that follows. The adopted conditions are

$$
G=M=c=1,
\qquad
\epsilon_0=\frac{1}{4\pi},
$$

$$
M=m_A+m_B=1.
$$

The object is non-spinning point particles in an adiabatic quasi-circular inspiral; the conservative dynamics is the leading order of Newtonian gravity + Coulomb, and the dissipation is the leading mass-quadrupole gravitational radiation and the leading electric-dipole electromagnetic radiation. Peters--Mathews [6,7] is the basis of gravitational quadrupole radiation and orbital decay, and as prior work treating both GW and EM in charged binaries, Liu et al. [8], Benavides-Gallego and Han [9], Zhang et al. [10] and Verma et al. [11] are referred to.

In the notation of Zhang et al. [10] a quantity corresponding to $1-\lambda_A\lambda_B$ is used; in this paper it is written $C$ to avoid a clash with the state vector $Z_8$.

## 3.2 Dimensionless mass and charge quantities

The charge-to-mass ratios are

$$
\lambda_A=\frac{q_A}{m_A},
\qquad
\lambda_B=\frac{q_B}{m_B}.
$$

The reduced mass and the symmetric mass ratio are

$$
\mu=\frac{m_A m_B}{M},
\qquad
\nu=\frac{m_A m_B}{M^2}.
$$

For $M=1$,

$$
\mu=\nu.
$$

Since this experiment has equal masses,

$$
\mu=\nu=\frac14.
$$

## 3.3 The relational quantity $C$ that collects the conservative motion

Collecting the leading-order effective central interaction of gravity and Coulomb, put

$$
\boxed{
C=1-\lambda_A\lambda_B.
}
$$

For quasi-circular orbits,

$$
\boxed{
\omega^2=\frac{C}{r^3},
}
$$

and the binding energy is

$$
\boxed{
E_{\rm orb}=-\frac{\mu C}{2r}.
}
$$

If $\lambda_A\lambda_B<0$ then $C>1$, and Coulomb acts as an attraction in the same direction as gravity. If $\lambda_A\lambda_B>0$ then $C<1$, and the Coulomb repulsion weakens the effective central attraction.

The real quasi-circular reference of this paper requires

$$
\boxed{C>0}.
$$

$C\le0$ is therefore outside the domain of the present analytic reference.

## 3.4 GW quadrupole radiation

The leading mass-quadrupole power is

$$
P_{\rm GW}
=
\frac{32}{5}\mu^2 r^4\omega^6.
$$

Substituting $\omega^2=C/r^3$,

$$
\boxed{
P_{\rm GW}
=
\frac{32}{5}
\frac{\mu^2C^3}{r^5}
}
$$

[6,7].

## 3.5 EM electric-dipole radiation

The electric dipole moment in the COM frame is, at leading order,

$$
\mathbf d
=
\mu(\lambda_A-\lambda_B)\mathbf r.
$$

Here put

$$
\boxed{
D=(\lambda_A-\lambda_B)^2.
}
$$

The leading dipole power is

$$
\boxed{
P_{\rm EM}
=
\frac23
\frac{\mu^2 D C^2}{r^4}
}
$$

[8--11].

From this expression, if

$$
\lambda_A=\lambda_B
\quad\Longrightarrow\quad
D=0,
$$

then **the leading electric-dipole term adopted in this paper** vanishes. The Coulomb conservative force itself, however, does not vanish.

## 3.6 Energy balance and radius evolution

The radiative loss is taken as

$$
\frac{dE_{\rm orb}}{dt}
=-(P_{\rm GW}+P_{\rm EM}).
$$

From

$$
\frac{dE_{\rm orb}}{dr}
=
\frac{\mu C}{2r^2},
$$

we obtain

$$
\boxed{
\frac{dr}{dt}
=-\frac{A}{r^2}-\frac{B}{r^3}.
}
$$

Here

$$
\boxed{
A=\frac43\mu D C
}
$$

is the coefficient on the EM dipole side and

$$
\boxed{
B=\frac{64}{5}\mu C^2
}
$$

is the coefficient on the GW quadrupole side.

In the radius evolution, therefore, two dissipative contributions with different distance dependences,

$$
\boxed{
r^{-2}\ \text{and}\ r^{-3},
}
$$

appear at the same time.

## 3.7 Analytic $t(r)$ and $\phi(r)$

Since

$$
\frac{dt}{dr}
=-\frac{r^3}{Ar+B},
$$

for $A\ne0$, putting

$$
F_t(r)
=
\frac{r^3}{3A}
-\frac{Br^2}{2A^2}
+\frac{B^2r}{A^3}
-\frac{B^3}{A^4}\ln(Ar+B),
$$

we have

$$
\boxed{
t(r)=F_t(r_0)-F_t(r)}.
$$

Also, from

$$
\frac{d\phi}{dt}=\sqrt{\frac{C}{r^3}},
$$

$$
\frac{d\phi}{dr}
=-\frac{\sqrt C\,r^{3/2}}{Ar+B}.
$$

For $A\ne0$, with

$$
F_\phi(r)
=
\sqrt C
\left[
\frac{2r^{3/2}}{3A}
-\frac{2B\sqrt r}{A^2}
+\frac{2B^{3/2}}{A^{5/2}}
\tan^{-1}\sqrt{\frac{Ar}{B}}
\right],
$$

$$
\boxed{
\phi(r)=F_\phi(r_0)-F_\phi(r)}.
$$

For equal-sign equal ratios, $D=0$ and hence $A=0$. In this case the limit is used directly,

$$
\boxed{
t(r)=\frac{r_0^4-r^4}{4B}}
$$

and

$$
\boxed{
\phi(r)=\frac{2\sqrt C}{5B}
\left(r_0^{5/2}-r^{5/2}\right)
}.
$$

The relative orbit is read out as

$$
\boxed{
x=r\cos\phi,
\qquad
y=r\sin\phi}.
$$

---

# 4. The analytic reference case

The reference condition is

$$
m_A=m_B=0.5,
$$

$$
\lambda_A=+0.30,
\qquad
\lambda_B=-0.30.
$$

Then

$$
C=1.09,
\qquad
D=0.36,
\qquad
\mu=\nu=0.25.
$$

The independent analytic reference for $r=50\rightarrow20$ is

$$
t_f=169032.31760706357,
$$

$$
\phi_f=767.0889994015913,
$$

$$
N_{\rm cyc}=\frac{\phi_f}{2\pi}
=122.086006046.
$$

In the reference case,

$$
\frac{P_{\rm EM}}{P_{\rm GW}}
=
\frac{5}{48}\frac{D}{C}r,
$$

so that at

$$
r_{\rm cross}
=
\frac{48C}{5D}
\simeq29.0667
$$

the EM dipole and the GW quadrupole are comparable. Outside, the EM dipole is relatively stronger; going inward, the $r^{-5}$ GW quadrupole overtakes it.


---

# 5. Extension from six states to eight states

## 5.1 The six persistent states of the fifth thought experiment

In the fifth thought experiment [5], the macro persistent information needed for the present gravitational-orbit reconstruction was taken to be

$$
Z_6=(U,P,E,H,Q,N).
$$

Here

$$
U=e^{i\chi},
\qquad
H=e^{i\phi},
\qquad
Q=e^{t/\tau_0},
\qquad
N=\nu.
$$

$P$ is the state corresponding to the orbital radius / semi-latus rectum, and $E$ is the eccentricity state. In the quasi-circular reference of this paper the update rate of $E$ is 0, and it is kept at its initial value 0.

## 5.2 The conserved relational quantities added in the charged system

Looking at the present analytic reference from the viewpoint of the information needed to determine the next state, what is newly needed is

$$
C=1-\lambda_A\lambda_B
$$

and

$$
D=(\lambda_A-\lambda_B)^2.
$$

Accordingly,

$$
\boxed{
Z_8=(U,P,E,H,Q,N,C,D).
}
$$

What matters here is not to assert that $C,D$ are "new fundamental physical degrees of freedom". They are **the persistent relational quantities needed to close the present restricted analytic representation**.

## 5.3 Table of persistent states

| State | Present meaning | Kind | Update across macro steps |
|---|---|---|---|
| $U$ | $e^{i\chi}$ | Phase state | Rotation |
| $P$ | Orbital radius / semi-latus rectum state | Dynamical | Decreases by radiation reaction |
| $E$ | Eccentricity state | Dynamical state | $dE/d\phi=0$ in this quasi-circular model |
| $H$ | $e^{i\phi}$ | Observed azimuth phase | Rotation |
| $Q$ | $e^{t/\tau_0}$ | Clock encoding | Multiplicative update |
| $N$ | $\nu$ | Conserved relational quantity | $N'=N$ |
| $C$ | $1-\lambda_A\lambda_B$ | Conserved relational quantity | $C'=C$ |
| $D$ | $(\lambda_A-\lambda_B)^2$ | Conserved relational quantity | $D'=D$ |

The individual $\lambda_A,\lambda_B$ are not passed to `transition(z)`. The difference between attraction and repulsion, and the presence or absence of the EM dipole, are expressed as differences of the state values $C,D$.

---

# 6. The state generator

## 6.1 Local update formulas with the phase as the independent variable

Combining

$$
\frac{dr}{dt}
=-\frac{A}{r^2}-\frac{B}{r^3}
$$

of the analytic reference with

$$
\frac{d\phi}{dt}=\sqrt{\frac{C}{r^3}},
$$

$$
\frac{dr}{d\phi}
=
-\frac{4}{3}\mu D\frac{\sqrt C}{\sqrt r}
-\frac{64}{5}\mu C\frac{\sqrt C}{r\sqrt r}.
$$

For $M=1$, $\mu=\nu=N$, so in state variables

$$
\boxed{
\frac{dP}{d\phi}
=
-\frac{4}{3}ND\frac{\sqrt C}{\sqrt P}
-\frac{64}{5}NC\frac{\sqrt C}{P\sqrt P}.
}
$$

Also

$$
\boxed{
\frac{dE}{d\phi}=0
}
$$

and

$$
\boxed{
\frac{dt}{d\phi}
=
\frac{P\sqrt P}{\sqrt C}.
}
$$

The local rate

$$
R(P,E,N,C,D)
=
\left(
\frac{dP}{d\phi},
0,
\frac{dt}{d\phi}
\right)
$$

can thus be evaluated **from the present states $P,N,C,D$ alone**.

## 6.2 Numerical resolution and encoding constants

One orbit is divided into 4000 macro steps,

$$
h=\frac{2\pi}{4000}.
$$

The clock state is encoded as

$$
Q=e^{t/\tau_0},
\qquad
\tau_0=10^4.
$$

$h$ and $\tau_0$ are not states specifying the physical case; they are treated as constants of the numerical resolution and of the readout encoding.

## 6.3 The 41-real-component full microstate

Of the eight persistent logical states, $U,H$ are complex, so in the implementation they become the 10 real components

$$
(U_R,U_I,P,E,H_R,H_I,Q,N,C,D).
$$

To these are added the RK4

$$
(k_1,k_2,k_3,k_4),
$$

the three components of each $k_j$, the stage states, the composite increment, the midpoint, $\Delta\phi$ and the 11-phase one-hot phase

$$
q=(q_0,\ldots,q_{10}),
$$

making the full state 41 real components.

Conceptually,

$$
\begin{aligned}
z=&(U_R,U_I,P,E,H_R,H_I,Q,N,C,D;\\
&k_1,k_2,k_3,k_4;P_s,E_s;d;P_m,E_m;\Delta\phi;q_0,\ldots,q_{10}).
\end{aligned}
$$

Here the work states are not pushed outside on the grounds that "they are used only temporarily until the next macro state, so they are not states"; if they are needed between microsteps they are made explicit as part of the full state.

## 6.4 Correspondence between the initialization code and the formulas

The initialization function actually used in the generator for the reference case is the following. It is the actual code of `run_strict_charged8_final.py` itself.

```python
def init_state(p0=50.0,n0=0.25,c0=1.09,d0=0.36):
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=p0; z[HR]=1.0; z[Q]=1.0
    z[N]=n0; z[C]=c0; z[D]=d0
    z[PS]=p0; z[PM]=p0; z[QB]=1.0
    return z
```

By `np.zeros`, all components not explicitly assigned start from 0. For the persistent logical states this corresponds to

$$
U_0=1+0i,
\qquad
P_0=p_0,
\qquad
E_0=0,
$$

$$
H_0=1+0i,
\qquad
Q_0=1,
$$

$$
N_0=n_0,
\qquad
C_0=c_0,
\qquad
D_0=d_0.
$$

In the reference case

$$
p_0=50,
\qquad
n_0=0.25,
\qquad
c_0=1.09,
\qquad
d_0=0.36.
$$

Since

$$
Q_0=1=e^{0/\tau_0},
$$

the initial readout time is $t_0=0$. Also,

$$
U_0=e^{i\chi_0}=1,
\qquad
H_0=e^{i\phi_0}=1
$$

express

$$
\chi_0=0,
\qquad
\phi_0=0.
$$

`z[PS]=p0` and `z[PM]=p0` are the initial values of the work registers for the RK4 stage / midpoint, and `z[QB]=1.0` expresses that the 11-phase one-hot phase is

$$
q=(1,0,\ldots,0),
$$

that is, the processing starts from phase 0.

What matters is that the individual

$$
\lambda_A,
\qquad
\lambda_B
$$

are not passed to `transition(z)`. The charge condition is mapped, before the initialization, to

$$
C=1-\lambda_A\lambda_B,
\qquad
D=(\lambda_A-\lambda_B)^2,
$$

and thereafter $C,D$ propagate as persistent states. Therefore, the only physical information changed in the five-condition test is the initial values of $C,D$ in the same state format.

## 6.5 Correspondence between the local interaction-rate code and the formulas

The local rate, evaluated uniquely from the present state, is collected in the following function in the actual code.

```python
@njit(cache=True)
def _rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt
```

These three components correspond directly to

$$
R(P,E,N,C,D)
=
\left(
\frac{dP}{d\phi},
\frac{dE}{d\phi},
\frac{dt}{d\phi}
\right).
$$

The line

```python
dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
```

evaluates

$$
\boxed{
\frac{dP}{d\phi}
=
-\frac{4}{3}ND\frac{\sqrt C}{\sqrt P}
-\frac{64}{5}NC\frac{\sqrt C}{P\sqrt P}
}
$$

as it is. The first term is the $r^{-2}$-type radial decay originating in the leading electric-dipole radiation, and the second is the $r^{-3}$-type radial decay originating in the leading gravitational quadrupole radiation, converted into phase derivatives by $d\phi/dt=\sqrt{C/P^3}$.

The second component of `return dp,0.0,dt` is 0 because this paper is restricted to circular orbits,

$$
E=0,
\qquad
\frac{dE}{d\phi}=0.
$$

Also,

```python
dt=p*rp/rc
```

corresponds to

$$
\boxed{
\frac{dt}{d\phi}
=
\frac{P\sqrt P}{\sqrt C}.
}
$$

This `_rates` does not switch between the attractive and the repulsive case, or between the GW and the EM term, by an external branch. The difference goes into the state values $C,D$ and is evaluated continuously within the same formula. In particular, for equal-sign equal charge-to-mass ratios

$$
D=0,
$$

so that only the EM dipole term becomes 0 automatically, without rewriting the code.

## 6.6 The actual code of the single map `transition(z)`

The one-hot sum of products of the candidate values is the following.

```python
@njit(cache=True)
def _blend(q,c):
    x=0.0
    for j in range(11): x += q[j]*c[j]
    return x
```

This implements the formula

$$
\boxed{
z_i'=\sum_{j=0}^{10}q_jF_{ij}(z)}
$$

as it is, as a sum of 11 products.

The core of the sole state update `transition(z)` used in the reference case is the following.

```python
@njit(cache=True)
def transition(z):
    """Sole legal state-transition map."""
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    # old work values
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]

    # stage rates are functions of OLD state only
    rP,rE,rT=_rates(p,n,c,dd)
    rPsP,rPsE,rPsT=_rates(ps,n,c,dd)

    cand=np.empty((30,11),dtype=np.float64)  # all non-q registers 0..29
    # default: identity candidate in all phases for every non-q register
    for i in range(30):
        for j in range(11): cand[i,j]=z[i]

    # k1: phase 0 compute, phase 10 reset
    cand[K1,0]=rP; cand[K1+1,0]=rE; cand[K1+2,0]=rT
    for a in range(3): cand[K1+a,10]=0.0
    # stage2 fields phase 1
    cand[PS,1]=p+0.5*HSTEP*k1[0]; cand[ES,1]=e+0.5*HSTEP*k1[1]
    # k2 phase 2, reset phase 10
    cand[K2,2]=rPsP; cand[K2+1,2]=rPsE; cand[K2+2,2]=rPsT
    for a in range(3): cand[K2+a,10]=0.0
    # stage3 fields phase 3
    cand[PS,3]=p+0.5*HSTEP*k2[0]; cand[ES,3]=e+0.5*HSTEP*k2[1]
    # k3 phase 4 uses OLD Ps,Es
    cand[K3,4]=rPsP; cand[K3+1,4]=rPsE; cand[K3+2,4]=rPsT
    for a in range(3): cand[K3+a,10]=0.0
    # stage4 fields phase 5
    cand[PS,5]=p+HSTEP*k3[0]; cand[ES,5]=e+HSTEP*k3[1]
    # k4 phase 6 uses OLD Ps,Es
    cand[K4,6]=rPsP; cand[K4+1,6]=rPsE; cand[K4+2,6]=rPsT
    for a in range(3): cand[K4+a,10]=0.0
    # combined d phase 7; reset 10
    for a in range(3):
        cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
        cand[DV+a,10]=0.0
    # midpoint phase 8; commit-safe values phase 10
    cand[PM,8]=p+0.5*dv[0]; cand[EM,8]=e+0.5*dv[1]
    cand[PM,10]=p+dv[0]; cand[EM,10]=e+dv[1]
    # dphi phase 9; reset 10. Circular model dphi=h exactly.
    cand[DPHI,9]=HSTEP; cand[DPHI,10]=0.0
    # commit phase 10 for persistent states
    ur=z[UR]; ui=z[UI]; hr=z[HR]; hi=z[HI]
    cand[UR,10]=ur*PH_RE-ui*PH_IM; cand[UI,10]=ur*PH_IM+ui*PH_RE
    cand[P,10]=p+dv[0]; cand[E,10]=e+dv[1]
    cd=math.cos(dphi); sd=math.sin(dphi)
    cand[HR,10]=hr*cd-hi*sd; cand[HI,10]=hr*sd+hi*cd
    cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
    # N,C,D remain identity candidates in all phases via default rows.
    cand[PS,10]=p+dv[0]; cand[ES,10]=e+dv[1]

    out=np.empty(NST,dtype=np.float64)
    for i in range(30): out[i]=_blend(q,cand[i])
    # fixed cyclic linear map q' = C_11 q (not external loop counter)
    out[QB]=q[10]
    for j in range(1,11): out[QB+j]=q[j-1]
    return out
```

The meaning of each phase is

$$
\begin{array}{c|l}
0 & \text{evaluation of }k_1\\
1 & \text{stage state for }k_2\\
2 & \text{evaluation of }k_2\\
3 & \text{stage state for }k_3\\
4 & \text{evaluation of }k_3\\
5 & \text{stage state for }k_4\\
6 & \text{evaluation of }k_4\\
7 & \Delta=(h/6)(k_1+2k_2+2k_3+k_4)\\
8 & \text{midpoint work state}\\
9 & \Delta\phi=h\\
10 & \text{commit to the persistent states}
\end{array}
$$

For example, in phase 1,

```python
cand[PS,1]=p+0.5*HSTEP*k1[0]
```

corresponds to

$$
P_s=P+\frac{h}{2}k_{1P},
$$

and in phase 7,

```python
cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
```

corresponds to the usual RK4 composition

$$
\Delta
=
\frac{h}{6}
\left(k_1+2k_2+2k_3+k_4\right).
$$

Unlike ordinary sequential code, however, a next value written just before within the same microstep is not used for the next operation. The candidates of each phase are built from the OLD full state read at the beginning of `transition(z)`, and only the candidate valid in that microstep is selected by the one-hot $q$. The next evaluation is performed only after moving to the next microstep.

In the commit phase,

```python
cand[P,10]=p+dv[0]
cand[E,10]=e+dv[1]
```

correspond to

$$
P'=P+\Delta P,
\qquad
E'=E+\Delta E.
$$

Also,

```python
cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
```

is

$$
Q'=Q\exp\left(\frac{\Delta t}{\tau_0}\right),
$$

which maintains the clock encoding

$$
Q=e^{t/\tau_0}
$$

without an additive time update.

For circular orbits the phase increment of one macro step is

$$
\Delta\phi=h,
$$

so the update of $H=e^{i\phi}$ is

$$
H'=He^{ih}.
$$

In the code the two real components `HR,HI` are updated as a rotation matrix.

One macro step applies the sole `transition` exactly 11 times, as in the following actual code.

```python
@njit(cache=True)
def macrostep(z):
    # exactly 11 applications of the sole transition, no alternate formula
    for _ in range(11): z=transition(z)
    return z
```

The macrostep itself is therefore not a separate physical update rule, but only the repetition that runs through the 11 phases once.

## 6.7 The conserved relational quantities $N,C,D$

`cand` is first filled, for every non-$q$ register, with the identity candidate of the OLD state.

```python
for i in range(30):
    for j in range(11):
        cand[i,j]=z[i]
```

Into the rows of $N,C,D$, no other candidate is written afterwards in any phase. Therefore, even through the one-hot sum of products,

$$
\boxed{
N'=N,
\qquad
C'=C,
\qquad
D'=D.
}
$$

This is not a copying of Python external constants at each step. $N,C,D$ are also contained in the full state, and propagate as identities as the result of passing through the same candidate rows and the same one-hot selection as the other registers.

## 6.8 Correspondence between the readout code and the formulas

The actual code of the readout for observation is the following.

```python
def readout(z):
    t=TAU0*math.log(z[Q])
    return t,z[P]*z[HR],z[P]*z[HI]
```

The first component is the inverse map

$$
\boxed{
t=\tau_0\ln Q}
$$

of

$$
Q=e^{t/\tau_0}.
$$

Since the present experiment is a quasi-circular orbit with

$$
E=0,
$$

the radius is

$$
r=P.
$$

Since $H=e^{i\phi}=H_R+iH_I$,

$$
PH
=
P(\cos\phi+i\sin\phi)
=
x+iy,
$$

so that

$$
\boxed{
x=P H_R,
\qquad
y=P H_I}.
$$

The return value of `readout(z)` is therefore exactly

$$
(t,x,y).
$$

This function is not part of the state update rule. In the generation loop it is called only when writing the log for storing the states, and its return value is not passed to `transition(z)`.

Since the terminal $P=20$ generally does not coincide with a macrostep grid point, only after the generation ends is a linear interpolation made from the two adjacent points,

```python
t1,_,_=readout(z); p1=z[P]
t0,_,_=readout(prev); p0=prev[P]
frac=(p0-20.0)/(p0-p1)
tcross=t0+frac*(t1-t0)
phicross=(nstep-1+frac)*HSTEP
```

which corresponds to

$$
f
=
\frac{P_0-20}{P_0-P_1},
$$

$$
t_{20}=t_0+f(t_1-t_0),
$$

$$
\phi_{20}=(n-1+f)h.
$$

The $t_{20},\phi_{20}$ obtained here are also readouts for comparison, and are not fed back into the orbit already generated.

## 6.9 Qualification concerning "fully multiplicative"

The strict single map of this implementation means

- no external physical argument,
- a sole `transition(z)`,
- OLD-state-only,
- explicit work states,
- one-hot sum of products,
- non-feedback readout.

On the other hand, it is not claimed that square roots, reciprocals and exponential functions have been removed from the primitive operations. They remain as "interaction functions computed uniquely from the present state". This paper is therefore not a proof that every operation has been reduced to multiplication alone.

---

# 7. Experimental design

## 7.1 The reference case

The reference case is

$$
\lambda_A=+0.30,
\qquad
\lambda_B=-0.30,
$$

$$
N=0.25,
\qquad
C=1.09,
\qquad
D=0.36,
\qquad
P_0=50.
$$

The initial state is

$$
U_0=1,
\qquad
H_0=1,
\qquad
Q_0=1,
\qquad
E_0=0.
$$

The generator continues the state update until $P\le20$ is reached.

## 7.2 The five-condition transferability test

From the series of integer multiples of the reference value $0.30$, the five cases satisfying $C>0$ in the present quasi-circular reference are chosen.

| case | $\lambda_A$ | $\lambda_B$ | $C$ | $D$ | Physical feature |
|---|---:|---:|---:|---:|---|
| n1 attractive | +0.30 | -0.30 | 1.09 | 0.36 | Coulomb attraction + EM dipole |
| n1 repulsive | +0.30 | +0.30 | 0.91 | 0 | Coulomb repulsion, no leading dipole |
| n3 attractive | +0.90 | -0.90 | 1.81 | 3.24 | Strong Coulomb attraction + strong dipole |
| n3 repulsive | +0.90 | +0.90 | 0.19 | 0 | Strong Coulomb repulsion, no leading dipole |
| n4 attractive | +1.20 | -1.20 | 2.44 | 5.76 | Strongest attractive side + dipole |

For $\lambda_A=\lambda_B=1.20$,

$$
C=1-1.44=-0.44,
$$

which is outside the domain of the present real quasi-circular analytic reference, so it is not included in the five cases. Corrections such as taking $|C|$ would be a different physical model, and are not made.

What is changed in the five cases is only the $C,D$ of the initial state; `transition(z)`, `macrostep(z)`, the 11-phase machine, the eight persistent states and the numerical resolution are not changed.

---

# 8. Results

## 8.1 First, "how the orbit changes with the condition"

The most important visualization in Paper 6 is the orbit itself under each condition.

![Figure 1: Independent analytic reference orbits for the five conditions](./01_軌道条件比較図/figures/figure01_five_case_reference_orbits_same_scale.svg)

**Figure 1.** The independent analytic reference orbits for the five conditions. The display range of $x,y$ is fixed to the same in all panels. Each condition is drawn for $r=50\rightarrow20$. For opposite signs, the larger $C,D$ the faster the inspiral; for equal-sign equal ratios $D=0$, so the leading EM dipole vanishes, and in particular the n3 repulsive case with $C=0.19$ requires a very large number of orbits.

For the same five conditions, the overlays with the stored strict eight-state generator are placed side by side.

![Figure 2: Overlay of the analytic reference and the strict eight-state generator for the five conditions](./01_軌道条件比較図/figures/figure02_five_case_strict_generator_overlay_panels.svg)

**Figure 2.** Orbit overlays of the analytic reference and the strict 8-state generator for all five conditions. The source figure of each panel is integrated as it is from the stored individual experiments; the orbit data themselves are not altered.

Comparing the radius and the accumulated number of orbits directly, the difference in density of the orbit figures can be read quantitatively.

![Figure 3: Radius and accumulated number of orbits](./01_軌道条件比較図/figures/figure03_five_case_radius_vs_cycles.svg)

**Figure 3.** Radius and accumulated number of orbits for $r=50\rightarrow20$. The speed of the orbital decay changes by several orders of magnitude over the five conditions.

Comparing only the first 10 orbits at the same spatial scale makes the difference in the decay speed at the initial stage easy to read.

![Figure 4: The first 10 orbits for the five conditions](./01_軌道条件比較図/figures/figure04_five_case_first_10_orbits_same_scale.svg)

**Figure 4.** Comparison of the orbits over the first 10 orbits. The display range is fixed, showing the difference of the initial inspiral between the conditions.

The total number of orbits and the final time of the analytic reference are as follows.

| case | $C$ | $D$ | $A_{\rm EM}$ | $B_{\rm GW}$ | cycles $50\to20$ | $t_f$ |
|---|---:|---:|---:|---:|---:|---:|
| n1 attractive | 1.09 | 0.36 | 0.1308 | 3.80192 | 122.086006 | 169032.317607 |
| n1 repulsive | 0.91 | 0 | 0 | 2.64992 | 364.132605 | 574545.646661 |
| n3 attractive | 1.81 | 3.24 | 1.9548 | 10.48352 | 16.671907 | 17450.031850 |
| n3 repulsive | 0.19 | 0 | 0 | 0.11552 | 3816.728393 | 13179536.011080 |
| n4 attractive | 2.44 | 5.76 | 4.6848 | 19.05152 | 8.348877 | 7507.620157 |

What this table shows is not the mere monotone relation "stronger charge, faster". The direction of $C$ changes with the sign, and for equal-sign equal ratios $D=0$, so that both the conservative force and the dissipation change at the same time.

## 8.2 Strict eight-state generation of the reference case

The strict generation of the reference case reached $P=20$ with

- macro steps: 488,345
- microsteps: 5,371,795

The readouts of the generator are

$$
t_{\rm gen}=169032.317607049,
$$

$$
\phi_{\rm gen}=767.088999401507.
$$

The independent analytic reference is

$$
t_{\rm ref}=169032.31760706357,
$$

$$
\phi_{\rm ref}=767.0889994015913.
$$

Therefore

$$
\boxed{
|\Delta t|=1.463923\times10^{-8}
}
$$

$$
\boxed{
|\Delta\phi|=8.435563\times10^{-11}
}
$$

In the comparison at identical physical times over the whole stored range,

$$
\max|\Delta r|=9.116260\times10^{-9},
$$

$$
\operatorname{RMS}(|\Delta r|)=3.771324\times10^{-9},
$$

$$
\max\|\Delta(x,y)\|=3.763773\times10^{-6},
$$

$$
\operatorname{RMS}(\|\Delta(x,y)\|)=2.130252\times10^{-6}.
$$

![Figure 5: Radius-time comparison for the reference case](./02_論文本文図/figure05_baseline_radius_vs_time.svg)

**Figure 5.** Radius--time comparison for the reference case. The independent analytic reference is used only for the comparison after the state generation.

![Figure 6: Residuals for the reference case](./02_論文本文図/figure06_baseline_reference_residuals.svg)

**Figure 6.** Residuals against the independent analytic reference in the reference case.

![Figure 7: Identity holding of N, C, D](./02_論文本文図/figure07_identity_states_NCD.svg)

**Figure 7.** Holding of $N,C,D$ in the reference case. The maximum drift over the whole run was 0.

![Figure 8: The first macrostep of the 11-phase machine/work states](./02_論文本文図/figure08_microphase_work_states.svg)

**Figure 8.** Evolution of the explicit work states in the first macrostep. The RK4 stage values are not reused directly within the same microstep; they are referred to for the first time as OLD state in the next microstep.

![Figure 9: EM dipole / GW quadrupole power ratio in the reference case](./02_論文本文図/figure09_em_to_gw_power_ratio.svg)

**Figure 9.** $P_{\rm EM}/P_{\rm GW}$ in the reference case. This figure is a diagnostic on the side of the independent analytic reference and is not fed back into the update of the state generator.

## 8.3 All five conditions generated by the same strict map

For the five conditions of Figure 2, the numerical summary is as follows.

| case | macro steps | $|\Delta\phi|$ | $|\Delta t|$ | time readout | $\Delta N,\Delta C,\Delta D$ |
|---|---:|---:|---:|---|---|
| n1 attractive | 488,345 | $8.436\times10^{-11}$ | $1.464\times10^{-8}$ | finite | 0, 0, 0 |
| n1 repulsive | 1,456,531 | $1.305\times10^{-9}$ | $3.203\times10^{-7}$ | finite | 0, 0, 0 |
| n3 attractive | 66,688 | $4.229\times10^{-9}$ | $8.739\times10^{-7}$ | finite | 0, 0, 0 |
| n3 repulsive | 15,266,914 | $3.653\times10^{-9}$ | N/A | nonfinite after step 6,316,005 | 0, 0, 0 |
| n4 attractive | 33,396 | $8.323\times10^{-9}$ | $1.545\times10^{-6}$ | finite | 0, 0, 0 |

What matters is that no separate generating formula was used for the five cases. In all cases,

$$
\boxed{
\text{same }transition(z)
+
\text{same 11-phase machine}
+
\text{same }Z_8
}
$$

were used, and the difference of the physical condition was given only by the initial $C,D$ of the state.

## 8.4 The time-readout limit of n3 repulsive

`n3_repulsive` has

$$
C=0.19,
\qquad
D=0,
$$

and is the slowest case, requiring about 3816.7 orbits.

The orbital states $P,H$ could be generated to the end, and the phase residual is finite. On the other hand, the clock encoding by

$$
Q=e^{t/\tau_0}
$$

exceeded the finite-precision range because of the long duration, and became nonfinite after step 6,316,005.

In this paper this case is not deleted, and no clipping of $Q$, renormalization, or switch to another clock is made. This is recorded as

$$
\boxed{
\text{not a failure of the orbit generation, but a finite-precision limit of the present exponential clock readout}.
}
$$

However, since the time cannot be read out to the end, no agreement of a finite value of $\Delta t$ is claimed for this case.

---

# 9. Discussion

## 9.1 Not the individual charges but relational quantities remained in the state

In the present restricted model, what `transition(z)` directly needs is not

$$
\lambda_A,
\lambda_B
$$

themselves, but

$$
C=1-\lambda_A\lambda_B,
\qquad
D=(\lambda_A-\lambda_B)^2.
$$

This shows that the state can be expressed by the **relational quantities** that actually appear in the interaction formulas, rather than by labelled parameters attached to individual bodies. This is, however, a representational fact in the present quasi-circular LO model; it does not prove that $C,D$ alone are fundamental in general charged binaries too.

## 9.2 Attraction/repulsion was not made an external branch

`transition(z)` contains

- no attractive / repulsive if-branch,
- no identity of bodies A/B,
- no gravity / Coulomb force-type branch.

The difference of the state values, $C>1,D>0$ for opposite signs and $C<1,D=0$ for equal-sign equal ratios, unfolds as it is into different orbits.

What was confirmed in this experiment is therefore

$$
\boxed{
\text{that different physical conditions can be handled not by switching algorithms but as differences of state values}.
}
$$

## 9.3 Two kinds of dissipation enter one rate

The radius evolution is

$$
\frac{dr}{dt}
=-\frac{A}{r^2}-\frac{B}{r^3},
$$

and the EM dipole and the GW quadrupole have different distance dependences. Nevertheless, converted into phase derivatives,

$$
\frac{dP}{d\phi}
=
-\frac{4}{3}ND\frac{\sqrt C}{\sqrt P}
-\frac{64}{5}NC\frac{\sqrt C}{P\sqrt P},
$$

they enter the same rate function, and `transition(z)` itself need not be switched.

In this sense this paper gives **an example in which several known interaction and dissipation terms can be internalized into one closed state-update form**.

## 9.4 This is not "a unification of gravity and electromagnetism"

That two kinds of contribution could be handled by the same state map, and that they are the same fundamental interaction in nature, are separate matters.

In this paper the known gravitational and electromagnetic terms are given in advance as the analytic reference, and their coefficients and the required states are audited. Therefore

$$
\boxed{
\text{same representation}
\not\Rightarrow
\text{same fundamental interaction}.
}
$$

## 9.5 Eight states are not necessarily minimal

Since $C,D$ are constructed from $\lambda_A,\lambda_B$, they are not independent fundamental degrees of freedom. Also, in the present symmetric series $|\lambda_A|=|\lambda_B|$, so that the admissible $(C,D)$ are strongly related.

Therefore it is not claimed that

$$
\boxed{
8\text{ is the mathematically minimal number of states}.
}
$$

What is confirmed now is only that

$$
\boxed{
\text{the present representation can be closed with }Z_8.
}
$$

---

# 10. Limits and non-claims

This paper does not claim the following.

1. **That the Einstein--Maxwell equations have been derived from eight states.**  
   This paper reverse-constructs a known leading-order analytic model.

2. **That gravity and electromagnetism have been fundamentally unified.**  
   That they can be handled in one state-update form is separate from the ontological identity of the interactions.

3. **That eight states are minimal or physically real.**  
   $Z_8$ is the set of persistent logical states that closes the present representation.

4. **Universality for general charged binaries.**  
   Eccentricity, spin, tidal/finite-size effects, generic mass ratios, higher PN, NLO radiation, strong-field merger and full nonlinear Einstein--Maxwell numerical relativity are not treated.

5. **Quasi-circular orbits with $C\le0$.**  
   The present real quasi-circular reference requires $C>0$.

6. **That all electromagnetic radiation vanishes when the leading dipole vanishes.**  
   What $D=0$ removes is the leading electric-dipole term adopted in this paper.

7. **That every internal operation has been reduced to multiplication alone.**  
   Square roots, reciprocals, exponential functions and so on remain as current-state interaction functions.

8. **That the 11 phases are a fundamental structure of nature.**  
   The 11 phases are the machine schedule that expands the present RK4 processing into OLD-state-only form.

Also, the nonfinite $Q$ of `n3_repulsive` is left as it is, as a numerical limit of the present readout encoding.

---

# 11. Reproducibility

(The folder names in the repository are in Japanese; they are transliterated here.)

## 11.1 Analytic reference

The analytic formulas, the initial conditions, the independent solver, the validation and the raw CSV files are stored in

`fifth_thought_experiment_multiplicative_state_expansion_and_order_verification_20260924/07_charged_binary_analytic_reference_20260925/`

Main code:

- `solve_charged_binary_analytic_reference_v1.py`
- `validate_charged_binary_analytic_reference_v1.py`
- `plot_charged_binary_analytic_reference_v1.py`

The state generator does not read these reference data.

## 11.2 Strict eight-state generator for the reference case

Stored in:

`fifth_thought_experiment_multiplicative_state_expansion_and_order_verification_20260924/11_charged8_strict_single_map_experiment_20260925/`

Main code:

- `01_CODE/run_strict_charged8_final.py`
- `01_CODE/postprocess_strict_charged8.py`

Audit evidence:

- `00_RULES_AND_AUDIT/preexec_audit.json`
- `03_ANALYSIS/generator_summary.json`
- `03_ANALYSIS/reference_comparison.json`

## 11.3 Five-case transfer experiment

Stored in:

`fifth_thought_experiment_multiplicative_state_expansion_and_order_verification_20260924/13_charged8_integer_multiple_5case_strict_numeric_20260925/`

Main code:

- `00_CODE_AND_RULES/run_strict_charged8_5cases.py`
- `00_CODE_AND_RULES/run_n3_repulsive_checkpointed.py`
- `00_CODE_AND_RULES/postprocess_strict_charged8_5cases.py`
- `00_CODE_AND_RULES/CASE_MATRIX.json`

All macrostep rows are stored in HDF5 part files. `n3_repulsive` was run in parts because of the run-time limit, but what was carried over between parts was only the complete 41-component full state; no switch to another generating formula was made.

## 11.4 Orbit comparison figures for Paper 6

Stored in:

`sixth_thought_experiment_charged8_gravity_coulomb_multiple_interactions_20260926/01_orbit_condition_comparison_figures/`

Reproduction code:

- `generate_paper6_orbit_condition_comparison_v1.py`

Products:

- `figures/figure01_five_case_reference_orbits_same_scale.svg`
- `figures/figure02_five_case_strict_generator_overlay_panels.svg`
- `figures/figure03_five_case_radius_vs_cycles.svg`
- `figures/figure04_five_case_first_10_orbits_same_scale.svg`
- `figures/paper6_orbit_condition_case_summary.csv`

The source figures of the strict overlays of the five cases are stored in `strict_sources/`; the orbit data were not changed when integrating them.

## 11.5 Auxiliary figures for the main text

Stored in:

`sixth_thought_experiment_charged8_gravity_coulomb_multiple_interactions_20260926/02_paper_main_figures/`

- `figure05_baseline_radius_vs_time.svg`
- `figure06_baseline_reference_residuals.svg`
- `figure07_identity_states_NCD.svg`
- `figure08_microphase_work_states.svg`
- `figure09_em_to_gw_power_ratio.svg`

---

# 12. Conclusion

In the fifth thought experiment [5], the gravitational-orbit generator was self-audited, and it was confirmed that, in the present representation, not the apparent two states but six persistent logical states and explicit machine/work states are needed.

In the sixth thought experiment, those audit rules were extended from the gravity-only system to the charged binary with gravity + Coulomb.

Organizing the independent analytic reference, the charge dependence in the present leading-order quasi-circular model was collected into

$$
C=1-\lambda_A\lambda_B
$$

and

$$
D=(\lambda_A-\lambda_B)^2.
$$

Adding these to the six states of the fifth thought experiment,

$$
\boxed{
Z_8=(U,P,E,H,Q,N,C,D),
}
$$

the charged quasi-circular inspiral could be generated without bringing the external $\lambda_A,\lambda_B$, an attraction/repulsion branch, or a force-type branch into `transition(z)`.

In the reference case, the differences from the analytic reference were

$$
|\Delta t|\simeq1.46\times10^{-8},
\qquad
|\Delta\phi|\simeq8.44\times10^{-11},
$$

with

$$
\Delta N=\Delta C=\Delta D=0.
$$

Furthermore, with the same `transition(z)`, the same 11-phase machine and the same eight persistent states, it could be transferred to five kinds of charge condition. The difference of the conditions unfolded into differences of orbit not as an external branch of the same algorithm but as the values of $C,D$ in the initial state.

The result obtained in this paper is therefore

$$
\boxed{
\begin{array}{c}
\text{for the restricted known gravity + Coulomb quasi-circular binary,}\\
\text{the conservative force and the two kinds of radiation reaction}\\
\text{can be internalized and reconstructed into a single synchronous state map}\\
\text{with no external physical parameters.}
\end{array}
}
$$

This is not "the completion of a unified theory". It is rather the next fixed point of this research series, which **audits, one step at a time, what information must be placed inside the state for the system to close when known physics is reverse-constructed**.

---

# Appendix A. Audit items of the strict single map

In the pre-execution audit, the following were confirmed.

| Item | Content | Result |
|---|---|---|
| A1 | The only argument of `transition` is `z`; the persistent states are $U,P,E,H,Q,N,C,D$ | pass |
| A2 | The sole generator is `transition`; no fast/collapsed/shortcut path | pass |
| A3 | The non-$q$ next state is the one-hot sum of products of 11 candidates; $q$ is a fixed cyclic linear map | pass |
| A4 | Every candidate uses the OLD full state only; the work registers are explicit | pass |
| A5 | The readout is at output only and is not fed back into the dynamics | pass |
| A6 | No A/B identity, no force-type branch | pass |
| A7 | The physical case values $N,C,D$ are states; the others are resolution, encoding and RK coefficients | pass |
| A8 | The old failed experiments are stored; no reference import, no other generator | pass |

The failed experiments of the old `09_...` and `10_...` are stored without deletion or overwriting.

---

# Appendix B. List of figures

1. **Figure 1** `01_orbit_condition_comparison_figures/figures/figure01_five_case_reference_orbits_same_scale.svg`  
   Independent analytic reference orbits for the five conditions, same scale.

2. **Figure 2** `01_orbit_condition_comparison_figures/figures/figure02_five_case_strict_generator_overlay_panels.svg`  
   Analytic reference vs strict eight-state generator for the five conditions.

3. **Figure 3** `01_orbit_condition_comparison_figures/figures/figure03_five_case_radius_vs_cycles.svg`  
   Radius vs accumulated number of orbits.

4. **Figure 4** `01_orbit_condition_comparison_figures/figures/figure04_five_case_first_10_orbits_same_scale.svg`  
   Same-scale comparison of the first 10 orbits.

5. **Figure 5** `02_paper_main_figures/figure05_baseline_radius_vs_time.svg`  
   Radius--time comparison for the reference case.

6. **Figure 6** `02_paper_main_figures/figure06_baseline_reference_residuals.svg`  
   Residuals against the independent analytic reference for the reference case.

7. **Figure 7** `02_paper_main_figures/figure07_identity_states_NCD.svg`  
   Identity holding of $N,C,D$.

8. **Figure 8** `02_paper_main_figures/figure08_microphase_work_states.svg`  
   11-phase machine/work states.

9. **Figure 9** `02_paper_main_figures/figure09_em_to_gw_power_ratio.svg`  
   EM/GW power ratio for the reference case.

---

# References

## This series

[0] Noriaki Kihara, **Designing a Research Project to Explore the Emergence of Spacetime, Particles, and Interactions from States and Relations: A Framework of Preliminary Design and Feasibility Verification for Foundational-Physics Model Exploration under Uncertainty**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22851944; Version DOI: 10.5281/zenodo.22851945.

[1] Noriaki Kihara, **The First Thought Experiment: Tracing the Equivalence Principle Down to Two Nameless Degrees of Freedom - Deriving the Sign of the Sum of Squares, the Curvature Radius, and Three Systems (Inertial, Uniformly Accelerated, Rotating) from the Conservation of an Area Readout**, v1.0 (2026-09-20). Concept DOI: 10.5281/zenodo.22857952; Version DOI: 10.5281/zenodo.22857953.

[2] Noriaki Kihara, **The Second Thought Experiment: Deriving the Coulomb-Type Inverse-Square Force from a Product Readout and an Area Clock - Keeping Two Values and a Linear Law, Obtaining the Focus, Attraction and Repulsion and Neutrality, and Fixing What Two Values Cannot Write**, v1.2 (2026-09-21). Concept DOI: 10.5281/zenodo.22867336; Version DOI: 10.5281/zenodo.22876596.

[3] Noriaki Kihara, **The Third Thought Experiment: Reading Out the State of Two Bodies a, b from the Orbit of Two Celestial Bodies a, b - Obtaining the Mass Ratio, the Relative Distance and the Relative Time from Orbit and Gravitational-Wave Readouts, with Two Values of Undetermined Meaning and the Knowledge $G=c=M=1$ Alone**, v1.5 (2026-09-23). Concept DOI: 10.5281/zenodo.22909660; Version DOI: 10.5281/zenodo.22909661.

[4] Noriaki Kihara, **The Fourth Thought Experiment: From a Fixed Action $S$ to a State-Dependent Action $S_n$ - Can the Relativistic Gravitational Orbits Computed Externally in the Third Thought Experiment Be Regenerated from the Sequential Interaction $X_{n+1}=S_nX_n$?**, v1.0 (2026-09-23). Concept DOI: 10.5281/zenodo.22919787; Version DOI: 10.5281/zenodo.22919788.

[5] Noriaki Kihara, **The Fifth Thought Experiment: Was There No Hidden State? - A Self-Audit of the Two-State Orbit Generation of the Fourth Thought Experiment, Tracking the State Closure from Five States to the Sixth Persistent State $\mathcal N=\nu$ -**, v2.1 (2026-09-26). Concept DOI: 10.5281/zenodo.22973086; Version DOI: 10.5281/zenodo.22973087.

## External

[6] P. C. Peters and J. Mathews, "Gravitational Radiation from Point Masses in a Keplerian Orbit," *Physical Review* **131**, 435-440 (1963). DOI: 10.1103/PhysRev.131.435.

[7] P. C. Peters, "Gravitational Radiation and the Motion of Two Point Masses," *Physical Review* **136**, B1224-B1232 (1964). DOI: 10.1103/PhysRev.136.B1224.

[8] L. Liu, O. Christiansen, Z.-K. Guo, R.-G. Cai, and S. P. Kim, "Gravitational and electromagnetic radiation from binary black holes with electric and magnetic charges: Circular orbits on a cone," *Physical Review D* **102**, 103520 (2020). DOI: 10.1103/PhysRevD.102.103520.

[9] C. A. Benavides-Gallego and W.-B. Han, "Gravitational Waves and Electromagnetic Radiation from Charged Black Hole Binaries," *Symmetry* **15**, 537 (2023). DOI: 10.3390/sym15020537.

[10] Z.-H. Zhang, T. Liu, S. Zhang, and Z.-K. Guo, "Post-Newtonian dynamics of charged compact binaries," *Physical Review D* **114**, 044088 (2026). DOI: 10.1103/4x8q-2kzm.

[11] S. Verma et al., "Post-newtonian dynamics of radiating charges: canonical formulation and binary inspiral laws," *European Physical Journal C* **86**, 830 (2026). DOI: 10.1140/epjc/s10052-026-16076-2.

---

# Revision history

- **v1.1 (2026-09-26):** The actual code of the initialization of the state generator, the local interaction rate, the sole `transition(z)`, the macrostep, the readout and the terminal interpolation was added to the main text, and the corresponding formulas and information flow were made explicit.
- **v1.0 (2026-09-26):** First version of the sixth thought experiment. The charged LO quasi-circular analytic reference, the strict single map with eight persistent states, the five-condition transferability, the orbit comparison figures, the reproducibility information, and the scope of claims and non-claims were integrated.
