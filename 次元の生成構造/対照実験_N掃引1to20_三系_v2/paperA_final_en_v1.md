# Space, Matter and the Clock Do Not Share a Common Birth Condition — What a Seed Must Be to Make a Particle, and the Lower Bound of Resolution

**Author:** Noriaki Kihara (WF System Co., Ltd.)　**Date:** 2026-08-10　**Version:** v1

**Version DOI:** [10.5281/zenodo.21874482](https://doi.org/10.5281/zenodo.21874482)　**Concept DOI:** [10.5281/zenodo.21874481](https://doi.org/10.5281/zenodo.21874481)

---

## What kind of paper this is

This paper does not propose a new dynamics. It imports, read-only and without any
modification, the dynamics established in the previous paper, and changes only
**where the seed is placed, how strong it is, its phase, the resolution, and the
number of divisions of the period**, then measures what happens.

Every number quoted here was recomputed independently from the stored run data and
checked mechanically against the existing records (Appendix E). Where an experiment
was run after its prediction had been fixed in a document, this is stated in the text.

---

## Premise

In the previous paper, 62 wave shapes were identified as candidates corresponding to
the particles of the Standard Model. In the course of that work it was found that,
in the inflation experiment, waves appearing at odd harmonics behave in a fermionic
way and waves at even harmonics behave in a bosonic way; and that when the reflection
(mixing) ratio of the interaction is around 0.7, a state mixing the fermionic and the
bosonic appears. The mechanism behind that mixing was left open.

In particular, geometric-series expansion of the inflationary kind occurred **without**
any seed, whereas a seed was **indispensable** for particle-like waves to appear. Only
one kind of seed, of comparatively large amplitude, had been tried, so which seed
produces which particle-like wave was also left open.

The influence of the resolution N was likewise left unanalysed.

This paper addresses those questions and obtains definite results.

---

## The results, sorted into five levels of confidence

So that the reader does not confuse them, the findings are classified up front.

**(1) Strong measured regularities**

- The onset of expansion is governed by the logarithm of the seed's **complex sum**
  (R² = 0.99987)
- The time at which particle-like waves finish appearing is governed by a power of the
  **power placed on odd bands** (R² = 0.996, 14 points)
- With a seed on even bands only, the odd bands stay **exactly zero**
  (42000 updates × 8 strengths)
- The classification by greatest common divisor **follows the number of divisions of
  the period** (278 pairs, zero failures)

**(2) Understood down to the mechanism**

- The onset depends on the complex sum because the first half of one update looks only
  at that sum and applies the same rotation to every location
- Nothing happens with an even-band-only seed because the coefficient of the interaction
  is built from the ratio of odd-band to even-band power; if the odd bands are zero,
  the coefficient is zero. That the odd bands never appear also follows, as parity
  conservation, from the three-wave form of the interaction

**(3) Empirical only**

- The exponent −1.073 itself (why it is close to −1 is not explained)
- That selectivity toward the targeted location peaks at intermediate strength
- The differences of regime with resolution
- Whether the first half of an update reads the **sum** or the **sum of squares** has
  not been separated experimentally (the update procedure says it is the sum)

**(4) Recorded as anomalies, interpretation withheld**

- At resolution 4 the clock runs 9.9 times faster
- At resolution 4 only, the readout crosses the value 0.6972
- At resolutions 8 and 10 the initial state could not be built
- At resolutions 13 and 20 the competition between planes reappears

**(5) Correspondence hypotheses to real physics (not tested here)**

- That odd and even bands correspond to fermions and bosons
- That the birth of a third direction corresponds to the birth of real space
- That the clock the system carries corresponds to physical time
- That a seed corresponds to a source of particle production

**None of (5) is tested in this paper**; all of it is inherited from the previous paper.
What is measured here is only which quantities become readable, and under which
conditions, inside the model.

---

## How the measurements were conducted

What is most often doubted in independent research is whether the conditions were
chosen after seeing the results. The following procedure was used.

- **Some experiments had their predictions fixed in a document before being run.**
  The four conditions of the phase-cancelled seed (section 3) and the test of changing
  the number of divisions of the period (section 8) were run after the predicted values
  had been written down, and hit as predicted
- **The dynamics is imported read-only, and checks its own digest before and after each
  run.** A mismatch stops the run on the spot
- **Every number in the claims was recomputed from the stored data and checked against
  the existing records** (Appendix E)
- **All 462 per-run record figures were inspected.** As a result, six of the claims in
  the text were found to need correction, and were corrected. In particular, the
  "window of strength" in section 2 was found by this inspection to be an artefact of
  the observation time
- Digests of 126 runs, 1174 MB of stored data and 13 programs are listed in the appendix

---

## The experimental system

Only what is needed to read the claims is given here.

### Points, relations, resolution

The system consists of **N points**. Every pair of points carries a **relation**, so
there are N(N−1)/2 relations. This N is called the **resolution**. The main experiment
uses N = 12, hence 66 relations.

### Waves and "locations"

**Each relation carries one wave, but that wave is not a single complex number.**
It is a set of values arranged along two periodic directions. The first direction is
divided into 16 parts of one turn, the second into 8, so there are 128 values in all.

A wave along a periodic direction can be classified by how many times its phase turns
in one lap. That count is called the **winding number**. Write k for the winding number
of the first direction and m for that of the second. The 128 values can be separated
into components indexed by the pair (k, m). In this paper a **location** means such a
pair. There are 128 locations.

**A location is not an address that the wave remembers.** Because the wave is a
distribution rather than a single complex number, its components simply carry indices.

Locations sharing the same k are grouped into a **band**; there are 16 bands. The
previous paper found that **waves on odd bands behave in a fermionic way and waves on
even bands in a bosonic way**, and this paper inherits that correspondence unchanged.
![Figure 1　The system. The map of 128 locations carried by each relation, and the relation between the background oscillation (black), the seed (red) and the location it targets (green frame). Dark columns are odd bands.](figA01_lattice_and_seeds_v1.png)

**Figure 1　The system. The map of 128 locations carried by each relation, and the relation between the background oscillation (black), the seed (red) and the location it targets (green frame). Dark columns are odd bands.**



### The background oscillation and the seed

Before a run starts, a wave of large amplitude is placed at location (k = 2, m = 0).
This is the **background oscillation**. It is common to every condition, including the
"no seed" case below. **It is because of this background oscillation that
inflation-like expansion occurs even without a seed.**

On top of that, adding a small amplitude at specified locations is called **seeding**.
The magnitude added is the **seed strength**, and the set of locations is the
**seed placement**.

### Updates, and the count kept from outside

The system repeats discrete updates. **There are no seconds in this paper; all time is
measured in numbers of updates.** That count is a scale we keep from outside; it is not
inside the system.

One update has two halves. The first forms, for each relation, the **sum of all 128
values taken as complex numbers**, decides a rotation from the argument of that sum,
and **applies the same rotation to all 128 locations**. The second half is the
interaction in which three waves meet and make another.
![Figure 2　One representative run (mixed 8 locations, strength 0.01, resolution 12, 42000 updates). Space, matter and the clock are born in this order, and the plane finally settles on one.](figA03_full_run_example_v1.png)

**Figure 2　One representative run (mixed 8 locations, strength 0.01, resolution 12, 42000 updates). Space, matter and the clock are born in this order, and the plane finally settles on one.**



### The clock — the time the system has of its own

As the system develops, it reaches a state in which **the whole turns by a constant
angle without changing shape**. Once that happens, one can read from the state itself
how many degrees it turned per update. **That is the tick of time for the system, and
in this paper it is called the clock.**

Two conditions must hold for it to be readable: the overlap with the immediately
preceding state must exceed a floor, and the amplitude of the component carrying the
rotation must exceed a floor. The first update at which both hold is when the
**clock is born**.

**The count of updates kept from outside and this clock are different things.**
"The clock is born at the 2nd update" means that after two updates on the outside
scale, a tick of time has appeared inside the system.

### The condensate

At the location where the background oscillation sits, the phases of the gathered
values align so that **squaring them and adding gives zero** (they cancel). This state
is called the **condensate**. Particle-like waves are born on top of this lump. Here we
say the condensate exists when the residue of that cancellation is below 10⁻⁵.

### The third direction

The distribution of the background oscillation initially lies within a single plane.
When a component growing outside that plane appears, a **third direction** stands up.
The first update at which the fraction outside exceeds 5% is when **space is born**.

But "having grown outside" does not fix the direction. If two or more planes compete
as the plane to grow out of, there is no telling which to take as reference. We measure
separately the degree to which one plane is singled out, and when that exceeds a fixed
value we say **the third direction is determined**.
![Figure 3　The whole main experiment. 5 seed placements × 8 strengths = 40 conditions, run with identical dynamics, resolution and number of updates. A red dash means not born.](figA02_experiment_matrix_v1.png)

**Figure 3　The whole main experiment. 5 seed placements × 8 strengths = 40 conditions, run with identical dynamics, resolution and number of updates. A red dash means not born.**



---

## 1) The form of the seed gives partial control over which kind of particle is born

Keeping the winding number k of the first direction fixed, we compared **two seeds
differing only in the winding number m of the second direction**, placed with the same
strength. One targets a partner that is an uncharged wave (neutrino type); the other
targets a partner that is an electron-type wave.

The number of updates at which inflation starts, and the number at which the clock is
born, **agreed exactly for 7 of the 8 strengths**. They differed by a single update only
at the weakest strength, where the run is indistinguishable from the no-seed case.

**Even so, the locations occupied by the resulting waves split into 16 and 64.**
For the electron type the pattern is a checkerboard: when the first winding number is
even so is the second, and when the first is odd so is the second.

That is, **"when it is born" and "what is born" are decided separately**. The former
is insensitive to the seed placement; the latter is decided by it.
![Figure 4　Locations occupied, for five seed placements (strength 0.1, after 42000 updates). Red frames mark the seed. The neutrino type occupies 16, the electron type a checkerboard of 64, and the even-band-only seed just 4.](figA06_support_structure_v1.png)

**Figure 4　Locations occupied, for five seed placements (strength 0.1, after 42000 updates). Red frames mark the seed. The neutrino type occupies 16, the electron type a checkerboard of 64, and the even-band-only seed just 4.**



### Being born as intended happens only at intermediate strength

Dividing the power at the targeted location by the power at non-targeted locations
measures how concentrated the result is at the target.

| Seed strength | 10⁻⁸ | 10⁻⁴ | **10⁻³** | 10⁻² | 0.03 and above |
|---|---:|---:|---:|---:|---:|
| Concentration | ≈ 0 | 0.009 | **10⁴** | 2×10³ | 0.08 |

**If the seed is too strong the wave spreads over all 128 locations and the aim stops
working.** The aim works best near strength 10⁻³, where concentration at the target
reaches ten thousand-fold; increasing the strength drops it by six orders of magnitude.

This concentration also varies by a factor of 120 with resolution (57 at N = 3, 0.47 at
N = 20).
![Figure 5　Concentration at the targeted location peaks near strength 10⁻³ and falls by six orders of magnitude as the seed is strengthened.](figA18_selectivity_window_v1.png)

**Figure 5　Concentration at the targeted location peaks near strength 10⁻³ and falls by six orders of magnitude as the seed is strengthened.**


![Figure 6　The map of 128 locations when concentrated (strength 10⁻³) and when fully spread (strength 0.1), for the same electron-type seed.](figA06b_ledger_pair_v1.png)

**Figure 6　The map of 128 locations when concentrated (strength 10⁻³) and when fully spread (strength 0.1), for the same electron-type seed.**



---

## 2) The seed amplitude appears to have a window, but that appearance is created by the length of the observation

If one looks only as far as 42000 updates, it appears that too weak a seed yields
inflation and three directions but no particle-like wave, that a stronger seed yields
particle-like waves, and that strengthening further has no additional effect.

However, when the run was extended to 300000 updates, **at strength 10⁻², where nothing
appeared to happen within 42000 updates, particle-like waves rose all at once at update
262,751**. It was not that it had no effect; it took time to take effect.

| Seed strength | 10⁻³ | **10⁻²** | 0.03 | 0.044 | 0.1 |
|---|---:|---:|---:|---:|---:|
| Updates until they finish appearing | not reached in 42000 | **262,751** | 25,760 | 13,391 | 2,352 |

Within the measured range no boundary strength appears; the data connect smoothly
through the following empirical law.

```
(updates until they finish appearing) ∝ (power placed on odd bands)^(−1.073)
goodness of fit R² = 0.996 (14 points)
```

Extrapolating this law, the number of updates diverges as the power placed on odd bands
approaches zero. **Therefore, seen within any finite observation time, an apparent lower
bound must arise.** Indeed a seed placed on no odd band at all (power exactly zero)
never finishes, at any strength.

However, **it has not been shown that every non-zero seed necessarily finishes in the
infinite-time limit.** What has been shown is that what looked like a threshold at
42000 updates collapsed at 300000, and that the 14 measurable points lie on this law.

Within the measured range, then, the strength of the seed does not decide whether
particle-like waves can form, but how long they take to form.

The amount produced differs qualitatively before and after.

- Before, the amount on odd bands is **exactly what was put in** (the seed power itself)
- After, **odd and even bands are roughly half and half**

The lower bound "too weak a seed yields no particle-like wave" is likewise likely an
artefact of the observation time. A long run is needed to check this (not done).
![Figure 7　Updates until the particle-like waves finish appearing. Stars mark points that only the 300000-update runs could reach. Points above the 42000 line do not fail to finish; they finish outside that line. Right: the same data against the power placed on odd bands.](figA16_saturation_power_law_v1.png)

**Figure 7　Updates until the particle-like waves finish appearing. Stars mark points that only the 300000-update runs could reach. Points above the 42000 line do not fail to finish; they finish outside that line. Right: the same data against the power placed on odd bands.**


![Figure 8　Space, matter and the clock do not share a birth condition. All four combinations occur.](figA05_three_births_v1.png)

**Figure 8　Space, matter and the clock do not share a birth condition. All four combinations occur.**



---

## 3) The strength of the seed changes the run-up to the onset of inflation

The run-up shortens in proportion to the logarithm of the seed strength.

```
run-up = 9.892 − 48.611 × ln(magnitude of the seed summed as complex numbers)
goodness of fit R² = 0.99987
```

What matters here is **not the total power of the seed but the magnitude obtained by
adding the seed as complex numbers**. Organised by total power instead, the scatter
worsens by a factor of 7.7.

This follows from the update procedure. The first half of an update **looks only at the
sum of the 128 values taken as complex numbers** and applies the same rotation to every
location. Two initial states with the same sum are therefore indistinguishable to it.

This was tested directly. Without changing the seed's power or placement at all, the
**phases were arranged to cancel**, so that only the sum becomes zero. It was then
predicted, **before measuring**, that the result would agree — down to the number of
updates — with the case having no seed on even bands. **All four strengths hit exactly.**

| Seed strength | ordinary seed | phase-cancelled seed | no seed on even bands |
|---:|---:|---:|---:|
| 0.01 | 74 | **116** | 116 |
| 0.03 | 29 | **43** | 43 |
| 0.044 | 22 | **33** | 33 |
| 0.1 | 12 | **17** | 17 |
![Figure 9　The run-up. Left: organised by the seed summed as complex numbers. Right: the same 15 points organised by total power.](figA10_runup_log_law_v1.png)

**Figure 9　The run-up. Left: organised by the seed summed as complex numbers. Right: the same 15 points organised by total power.**


![Figure 10　The phase-cancellation control. Without changing power or placement, the result agrees update for update with the case having no seed on even bands.](figA11_phase_cancellation_v1.png)

**Figure 10　The phase-cancellation control. Without changing power or placement, the result agrees update for update with the case having no seed on even bands.**



### This control experiment alone cannot separate the sum from the sum of squares

The cancellation used a triple of equal amplitudes with phases 0, 2π/3 and 4π/3.
That arrangement makes the complex sum zero, but **it also makes the sum of squares
zero** (both of order 10⁻¹⁵).

**This experiment alone therefore cannot tell whether the first half of an update reads
the "sum" or the "sum of squares".** The paper says it is the sum because the update
procedure itself uses only the argument of the sum. The experiment is consistent with
that, but it does not on its own single out the sum.

To separate them one needs an arrangement whose sum is zero but whose sum of squares is
not (for example a pair with phase difference π). That was not done here.

### The complex sum matters only for the onset

The "updates until the waves finish appearing" of section 2 is instead decided by total
power. Two different rules of decision coexist in the same system.

| What is observed | Decided by | Quantity that matters |
|---|---|---|
| Onset of inflation | first half of an update | the **sum** taken as complex numbers |
| Time when particle-like waves finish appearing | second half | the total **power** |

The first half looks only at the sum. The second is proportional to a coefficient built
from the ratio of odd-band to even-band power. Hence the quantity that matters depends
on which phenomenon is observed.

A seed placed only on even bands satisfies both at once. Its odd-band power is zero, so
the finishing time is infinite (and indeed it never finishes at any strength), while its
onset lies on the straight line of the complex sum.
![Figure 11　Two rules of decision in one system. Left: the onset, decided by the first half of an update. Right: the finishing time, decided by the second half.](figA17_two_laws_v1.png)

**Figure 11　Two rules of decision in one system. Left: the onset, decided by the first half of an update. Right: the finishing time, decided by the second half.**


![Figure 12　Updates to onset for all 40 conditions. Re-plotted against the seed summed as complex numbers (right), the five placements fall on one curve.](figA09_tau_space_all_conditions_v1.png)

**Figure 12　Updates to onset for all 40 conditions. Re-plotted against the seed summed as complex numbers (right), the five placements fall on one curve.**



---

## 4) There is a minimum resolution: at N = 1, 2 neither inflation nor three directions appeared

At N = 1 there is not a single relation between points, and at N = 2 only one, so the
system cannot be built at all.

### There are singular resolutions: at N = 3, 4 the third direction is not determined

**At N = 3 and N = 4 the planes compete and no single plane is selected, so the third
direction is not determined.** Three planes compete at N = 3 and two at N = 4. Moreover
this competition **never settles; it keeps fluctuating over the whole run.** There is no
steady state.

N = 5 and 6 show yet another aspect: **quiet stretches alternate with violently
fluctuating ones.** Only from N = 7 upward does a single plane dominate.

**At N = 8 and N = 10 the initial state itself could not be built.** This, however, is
because the implementation stops after one failed attempt; it has not been confirmed
that it is structurally impossible (unverified).

**Seeding acts to resolve the competition.** Among the 16 resolutions that could be
built, the number showing competition falls from 8 with no seed, to 4 with a seed at one
location, to 2 with seeds at eight locations. **N = 3 and N = 4 alone cannot be resolved
by any amount of seeding.**

There is no relation of the kind "resolutions that compete have longer run-ups". The two
are separate properties.
![Figure 13　Competition between planes, by resolution and seed. Yellow marks N = 3, 4 (competing under every condition); red marks resolutions where the initial state could not be built. Right: the more seed locations, the fewer resolutions compete.](figA12_resolution_regimes_v1.png)

**Figure 13　Competition between planes, by resolution and seed. Yellow marks N = 3, 4 (competing under every condition); red marks resolutions where the initial state could not be built. Right: the more seed locations, the fewer resolutions compete.**


![Figure 14　The judgements for four conditions over resolutions 1–20. In the no-seed rows only space is marked, with matter and the clock absent at every resolution.](figA13_nsweep_birth_matrix_v1.png)

**Figure 14　The judgements for four conditions over resolutions 1–20. In the no-seed rows only space is marked, with matter and the clock absent at every resolution.**



### Three anomalies that appear only at N = 4

1. Two planes compete and the third direction is not determined
2. **The clock ticks 9.9 times faster than at other resolutions**
3. **It is the only resolution at which the readout crosses the special value 0.6972
   reported in the previous paper**

Examining all 126 runs, the only ones to exceed that value are **the two runs with a
single-location seed at N = 4**, reaching a maximum of 0.8019. The neutrino type and the
electron type agree to 12 digits on that value. The closest approach is 0.0000068
(0.001% in relative terms).

**But it only crosses; it does not stay.** The readout lies within ±0.001 of that value
for **only 9 of the 42000 updates**. At N = 4 the system keeps fluctuating violently
after the waves have finished appearing, and the fluctuation happens to cross the value.
At N = 12 the maximum is 0.6786, short of that value by 2.7%.
![Figure 15　The anomaly at resolution 4. The number of competing planes never settles (left), and that fluctuation crosses the special value 0.6972 (right).](figA04b_N4_anomaly_v1.png)

**Figure 15　The anomaly at resolution 4. The number of competing planes never settles (left), and that fluctuation crosses the special value 0.6972 (right).**



### From N ≥ 12 up to N = 20, tested, there was broadly no large change of behaviour

The run-up stays within 206–277 updates with no monotone trend.

**There are two exceptions.** Competition between planes reappears at N = 13 and N = 20.
In particular **at N = 20 the third direction is not determined**.

### After the waves finish appearing, the plane stops being determined again

"A single plane dominates from N = 7 upward" refers to the state before the
particle-like waves finish appearing. Afterwards the number of competing planes keeps
fluctuating at every resolution.

---

## 5) There are not four things that "are born" but six, and a strong seed destroys the last two

What the experiment judges is the following six.

Whether the system could be built / whether **space is born** (a third direction stands
up) / whether **matter is born** (waves appear on odd bands) / whether **the clock is
born** (the system can read its own time) / whether **the third direction is
determined** / whether **the condensate exists**.

As the seed is strengthened, **the first four increase but the last two are lost**.

| Seed strength | space | matter | clock | third direction determined | condensate |
|---:|:--:|:--:|:--:|:--:|:--:|
| 10⁻¹⁵, 10⁻⁸ | ○ | ○ | **×** | ○ | ○ |
| 10⁻⁴, 10⁻³ | ○ | ○ | ○ | ○ | ○ |
| 10⁻² | ○ | ○ | ○ | ○ | **×** (with many seed locations) |
| 0.03 and above | ○ | ○ | ○ | **×** | **×** |

**A strong seed produces matter and a clock at the cost of the certainty of the third
direction and of the condensate.** Things do not simply come into being; there is a
trade.

Furthermore, extending the run to 300000 updates, **the condensate is eventually
destroyed even at strength 10⁻²**. That it looked intact at 42000 updates only means it
had not yet been destroyed.
![Figure 16　The six criteria. Above the thick line the first four increase with seed strength; below it the last two are lost.](figA19_six_criteria_tradeoff_v1.png)

**Figure 16　The six criteria. Above the thick line the first four increase with seed strength; below it the last two are lost.**



---

## 6) With a seed on even bands only, the interaction itself never occurs even once

With a seed on even bands only, **no wave whatsoever appears on the odd bands**, over
the whole of 42000 updates and all 8 strengths. More than that, **only 2 of the 16 bands
ever carry a wave** — the background oscillation and the seed itself — and **not even
the partner band that this seed was aiming at stands up.**

The reason is this. The coefficient that sets the strength of the interaction is built
from the ratio of odd-band to even-band power. If the odd bands are zero, that
coefficient is zero. If the coefficient is zero, the whole interaction expression is
zero. **With a seed on even bands only, the interaction never occurs and the system runs
on the first half of the update alone.**

Why no wave appears on the odd bands follows from the update procedure itself. The
interaction has the form of three waves meeting, so whether the outgoing wave lies on an
even or an odd band is decided by the sum of the parities of the three incoming ones.
Gather only even ones and only even ones come out.

In exchange, this seed **preserves the certainty of the third direction and the
condensate at every strength**. Nothing is destroyed because no interaction occurs. It
is the clearest control for the trade described in section 5.
![Figure 17　Growth by band. Left: the neutrino type, where the targeted partner band stands up. Middle: the even-band-only seed, where only 2 bands ever carry a wave. Right: the coefficient setting the interaction strength, identically zero for the even-band-only seed.](figA04_band_evolution_sumrule_v1.png)

**Figure 17　Growth by band. Left: the neutrino type, where the targeted partner band stands up. Middle: the even-band-only seed, where only 2 bands ever carry a wave. Right: the coefficient setting the interaction strength, identically zero for the even-band-only seed.**


![Figure 18　With a seed on even bands only, the odd bands stay exactly zero over the whole of 42000 updates and all 8 strengths.](figA07_parity_selection_rule_v1.png)

**Figure 18　With a seed on even bands only, the odd bands stay exactly zero over the whole of 42000 updates and all 8 strengths.**



---

## 7) Time appears twice in this system: the order of updates is not time

Two things that could be called time appear here. One is the **count of updates kept
from outside**, which is merely a serial number we attach in order to run the
experiment; it is not inside the system. The other is the **clock the system reads from
its own state**, which exists only once the whole turns by a constant angle without
changing shape.

**These two must not be conflated.** No quantity named time is fed into the dynamics of
this system; what is fed in is only the update procedure. Nevertheless, as the system
develops, a tick of rotation becomes readable from the relations among states.
**Time is not something given, but something that became readable.**

### Having a clock, its being constant, and its rate are three different things

The number of updates at which the system becomes able to read its own time, and the
number at which the rate it reads settles to a constant value, are different quantities.

- With seeds at eight locations and 4000 updates, **the rate does not settle for N ≥ 15**
  (the time has become readable)
- Applying the phase cancellation of section 3, **the number of updates to settle falls
  from 15638 to 18, a factor of 870** (at strength 0.1, from 96 to 6)
- The rate of the clock itself varies with condition. Of 132 runs, **67 deviated from the
  reference by more than 5%**. The largest is the factor 9.9 at N = 4

**That the clock has been born does not mean the system keeps a constant time.**
![Figure 19　Updates at which the clock is born (left) and the peak of the mixing ratio (right). In the grey region on the left, no clock is born within the 42000 updates observed.](figA08_amplitude_window_v1.png)

**Figure 19　Updates at which the clock is born (left) and the peak of the mixing ratio (right). In the grey region on the left, no clock is born within the 42000 updates observed.**



---

## 8) The classification by greatest common divisor found in the previous paper survives a change in the number of divisions of the period

The previous paper found that, inside the background sea, the properties of a wave depend
not on the winding number of the second direction itself but on **the greatest common
divisor of that winding number and the "number of divisions of the period"**.

By the number of divisions of the period is meant, as in the description of the system
above, into how many parts one lap of the second direction is divided. A winding number
has meaning only as a remainder modulo that number, so changing it should change which
winding numbers fall into the same class. The previous paper registered this as a
question to be tested. This paper tested it.

Re-running at the same number of divisions, 16, **all 62 values agreed bit for bit** with
the stored results of the previous paper. Divisions of 8, 12 and 32 were then added, and
**all 278 pairs came out as predicted, with zero failures**.

The decisive case is the pair of winding numbers 1 and 3. If the greatest common divisor
governs, this pair should fall in the same class at divisions 8, 16 and 32, and in
**different classes at 12 alone**. Measurement gives **agreement to 14 digits** at 8, 16
and 32, and a **60% difference at 12 alone**.
![Figure 20　The test of changing the number of divisions of the period. Points of the same colour line up at the same height. Winding numbers 1 and 3 agree to 14 digits at divisions 8, 16 and 32, and differ by 60% at 12 alone.](figA14_divisor_class_register_order_v1.png)

**Figure 20　The test of changing the number of divisions of the period. Points of the same colour line up at the same height. Winding numbers 1 and 3 agree to 14 digits at divisions 8, 16 and 32, and differ by 60% at 12 alone.**



---

## Open questions

1. Whether the failure to build the initial state at N = 8 and N = 10 is structural or
   an artefact of the implementation cutting off after one attempt
2. Whether "too weak a seed yields no particle-like wave" is likewise an artefact of the
   observation time
3. The long run at strength 10⁻² with many seed locations. The saturation law (Figure 7)
   is drawn from 14 points; with those 2 added it becomes 16. The exponent may shift, but
   the conclusion that a finite observation time must produce an apparent lower bound
   does not
4. Whether the classification by greatest common divisor also holds in the many-body
   system (section 8 was verified in the two-wave system)
5. A test changing the number of divisions of the first direction (16 in this experiment)
6. Measurement of the coefficient the interaction actually uses (the relation to the
   value 0.6972 of the previous paper is decided there)
7. Matching the 62 kinds of the previous paper against which location in this experiment
   corresponds to which particle

The following go beyond the scope of this paper and are deferred to a sequel.

8. **An orthogonal experiment separating the run-up law from the saturation law.**
   The phase-cancellation experiment here fixed power and placement and varied the
   complex sum. Its orthogonal counterpart would **fix the sum and vary only the power**;
   the onset should be unchanged and only the finishing time should move. The same
   experiment would settle the "sum versus sum of squares" question of section 3, by
   including an arrangement whose sum is zero but whose sum of squares is not
9. **Falsification of the saturation law by extrapolation.** Freeze the exponent at the
   value from the existing 14 points, predict the finishing time for unmeasured strengths
   in advance, and go for the hit without refitting
10. **A formal audit of how far the system closes on relations alone.** For the three
    readouts used here (the plane and the third direction; bands and mixing ratio; the
    clock criterion), classify whether the inputs close on sets of relations alone, or
    whether quantities attached to individual points have crept in.

    **This must be taken in two stages.** Even if the readouts are confirmed to close on
    relations alone, that shows only that **the readout maps close on relations**, not
    that **the system itself runs on relations alone**. For the latter, the same audit is
    required of the construction of the initial state, the first half of the update, and
    the interaction.

    That is, one must check in turn, for each stage — building the state, updating,
    interacting, reading out — whether any quantity attached to a point itself has slipped
    in. This paper has not verified even the last of those stages.

---

# Appendix　Full list of runs, figures and programs

This appendix was generated automatically by `make_paperA_appendix_v1.py`, which scans
the folder. Nothing in it was transcribed by hand. The Japanese edition carries the
identical tables; only the headings are translated here.

## Appendix A　The runs

126 stored run files (NPZ), 1174 MB in total. Runs differing in seed placement, seed
strength, number of updates or resolution are counted as different runs.

- 40 conditions of the main experiment: 5 seed placements × 8 strengths
  (resolution 12, 42000 updates), all present
- Resolution sweep N = 1…20 for 4 conditions (no seed, neutrino type, electron type,
  mixed), 4000 updates
- Reproduction controls, phase-cancelled arms and their replicate, resolution-4 runs,
  and long runs at 300000 updates

## Appendix B　The figures

20 figures in the body; 462 per-run record figures; 41 contact sheets used to inspect
all of them. Naming is
`fig_<seed placement>[_T<updates>][_d<strength>][_rep-<tag>]_<family>[_N<resolution>]_v2.png`,
so any figure cited in the text is identified uniquely by its file name.

| Family | Content | Count |
|---|---|---:|
| 4panel | space / number of competing planes / clock / cancellation residue | 126 |
| mix | power per band and mixing ratio | 126 |
| ledger | map of the 128 locations, and growth at the targeted one | 78 |
| summary | per-resolution summary | 66 |
| birth_matrix | the six criteria | 66 |

## Appendix C　The programs

Thirteen programs, plus five imported read-only. SHA-256 digests of every one are listed
in the Japanese edition and in `論文A_凍結マニフェスト_v1.json`.

| Role | File |
|---|---|
| Main experiment | `run_nsweep_three_series_v2.py`, `run_tb_nsweep_1to20_v1.py` |
| Additional seed placements | `run_missing_seed_sweeps_T42000_v1.py` |
| Phase-cancelled seed | `run_phase_balanced_mixed_v1.py`, `run_phase_balanced_mixed_grid_v1.py` |
| Test of the number of divisions | `run_divisor_class_register_order_v1.py` |
| Recomputation and cross-check of every claimed number | `aggregate_paperA_claims_v1.py` |
| Figures | `make_paperA_figures_v2.py`, `make_paperA_figures_audit_v1.py` |
| Closest approach to 0.6972 over all runs | `probe_alpha_root_closest_approach_v1.py` |
| Contact sheets for full visual inspection | `make_contact_sheets_v1.py` |
| This appendix | `make_paperA_appendix_v1.py` |
| Long runs and resolution-4 runs | `run_stage4_longtime_orchestrator_v1.py` |

The dynamics itself (interaction, plane and third-direction readout, band and mixing
readout, clock criterion) and the two-wave map used in section 8 are imported read-only
and never modified. Each run verifies their digests before and after.

## Appendix D　Pre-registrations, reports and records

Pre-registration documents (predictions fixed before running), result reports, the
full-figure audit record, and the JSON files holding every number.

## Appendix E　Results of the cross-check

Every number claimed was recomputed independently from the stored data and checked
mechanically against the existing records.

| Check | Result |
|---|---|
| Onset and clock counts recomputed from stored data vs. the records | 40 conditions, 0 mismatches |
| Number of locations occupied | agrees |
| Odd bands exactly zero for the even-band-only seed | maximum 0.0 |
| Run-up law | 9.8922 − 48.6108 × ln(complex sum), R² = 0.999868 |
| Same points organised by total power | R² = 0.992202 |
| The four phase-cancelled conditions | all as predicted |
| Reproduction of the previous paper (16 divisions) | 62 values, largest discrepancy 0.0e+00 |
| Test of changing the number of divisions | 278 pairs, 0 failures |
| Law for the finishing time | ∝ (odd-band power)^(−1.073), R² = 0.9960 (14 points) |
