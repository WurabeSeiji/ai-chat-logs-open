# I Put Gravity and the Coulomb Force into One and the Same State Update - After First Flushing Out Every "Hidden State": the Fifth and Sixth Thought Experiments

At the end of the last article, I wrote this.

> The next problem is whether the rule for making the action, assembled here from a known solution and a known expression, can be derived from a smaller number of relational quantities.

I intended to start on that right away.

But there was one thing to do first.

Inside the implementation that I was calling "an orbit made from internal states and interactions alone", was there really no information brought in from outside?

I decided to doubt that point first.

## First, I audited my own previous result

In the fourth thought experiment, by letting the action itself change with the state, a family of orbits including the periastron advance and the shrinking by gravitational waves could be regenerated with one and the same rule.

On the surface, the orbit was made from two internal states.

But when I reread the program from the viewpoint of "which information is remembered until the next computation", it was not two.

The phase of the orbit, the size of the orbit, the flattening, the azimuth, the clock.

And on top of that, the mass-ratio information that decides the reaction of the gravitational waves was still being given from outside.

In other words, the numbers agreeing and the system being closed by the states alone were not the same thing.

That is the fifth thought experiment.

## I stopped saying "it is an intermediate calculation, so it is not a state"

In ordinary numerical computation, intermediate values are just working variables.

But if a value is used in the next computation, information is being remembered in the meantime.

Then, under the rules of this research, it must not be hidden outside.

I decided to make the intermediate values explicit as states too, for as long as they are needed.

In the previous orbit computation, several stages of calculation happen inside one update.

I split this into 11 small operations, and even put the information

"which operation am I doing now"

inside the state.

Furthermore, the mass-ratio information that had been passed from outside was moved inside, as the sixth persistent state.

Even with this change, the existing orbit output did not change.

What matters here is not that I found out "six is the true number of states in nature".

Rather the opposite.

Do not hide, on the grounds of convenience, the information needed to decide the next state.

Having carried that audit method through to the end, once, is the result.

## So this time I did the experiment in the reverse direction

From here it is the sixth thought experiment.

Until now, this research has been thinking in the direction: can motion be produced from as few states and relations as possible?

This time I went somewhat in reverse.

First, following the currently known theory, I make the "motion that is the answer".

The object is a two-body system in which gravity and the Coulomb force act at the same time.

And as the causes of the orbit shrinking, I put in both

radiation by gravitational waves, and

dipole radiation by electromagnetic waves.

However, I am not trying to solve the whole general charged two-body problem.

I restricted it to a range where an analytic reference can be built clearly: no spin, quasi-circular orbits, and a low-order approximation.

And that known orbit is placed only on the answer-checking side.

The regenerating side does not look at the answer while generating.

## I removed the branches "this is gravity" and "this is electromagnetism"

What I wanted to test here was not the computation of a charged binary itself.

The question is more limited.

Can motion in which gravity and the Coulomb force act at the same time, and in which the orbit shrinks by two kinds of radiation, be written as a single closed state update?

Inside the generator I did not place the branches

"this is particle A",

"this is particle B",

"here I compute gravity",

"here I compute the Coulomb force",

"this is the attractive case",

"this is the repulsive case".

Instead, I put into the state only two relational quantities decided by the charge condition.

One is a quantity expressing the strength of the conservative interaction of gravity and the Coulomb force combined.

The other is decided by the difference of the charge-to-mass ratios of the two bodies, and expresses the strength of the electromagnetic dipole radiation.

I added these two to the previous six states.

The persistent states became eight.

Here too, however, I am not saying "eight is the minimal number of states in nature". Only that, in order to close this approximate model, the present representation makes eight explicit.

## Once it starts moving, it does not look at the external orbit solution

The orbit for answer checking is built independently from the known theory and stored.

The state generator does not read it.

What is used during generation is the present state only.

The operation that makes the next state from one state was fixed to one and the same transition.

Including the intermediate calculations, it reads only the previous state and makes the next state all at once.

The position and time that are read out are not returned to the generator either.

Here I used, as they were, the audit rules made in the fifth thought experiment.

## This is how the orbit came out

Since the shrinking of the orbit is hard to see in an ordinary plot, where the lines overlap, I drew the first 10 orbits with the radius on a logarithmic scale.

[Image position: the first 10 orbits made by the strict eight-state generator, with the radial direction on a logarithmic scale so that the slight orbital decay per orbit is visible; drawn from the stored trajectory only
figure07_first_10_orbits_log_radial_scale (from 11_.../04_FIGURES_SVG)]

This figure makes it easy to see that the orbit goes a little further inward with each revolution.

But with this figure alone, the possibility remains that it was tuned for one condition only.

So, as last time, I moved the same update rule to other conditions.

## I changed the conditions, but not the update rule

I changed the charge condition to five kinds.

A condition with opposite-sign charges, where the attraction becomes stronger.

A condition with equal-sign charges, where a Coulomb repulsion enters.

A condition with strong electromagnetic dipole radiation.

A condition where the leading electromagnetic dipole radiation vanishes because the charge-to-mass ratios are equal.

Five cases with rather different properties.

Even so, what I changed was only the two relational quantities put into the state at the start.

The transition was not changed.

The 11-stage internal processing was not changed.

The kinds of states were not changed.

The numerical resolution was not changed.

Even so, in all five cases the orbital states could be generated to the end.

## The differences between the orbits are very large

The number of orbits needed to shrink from the same radius 50 to 20 changed greatly with the condition.

In the fast condition, about 8.35 orbits.

In another condition, about 16.7 orbits.

In the reference condition, about 122 orbits.

On the repulsive side, about 364 orbits.

In the slowest condition, about 3817 orbits.

With the same update rule, only the difference of the relational quantities put into the initial state unfolded into orbits this different.

[Image position: comparison of the number of orbits needed to shrink from radius 50 to 20 for the five charge conditions; the same state update rule, only the relational quantities of the initial state changed
figure03_five_case_radius_vs_cycles.png]

What I want to look at here is not the numbers of orbits themselves.

Whether different interaction conditions could be handled within the same state update, rather than by a separate program for each case.

That is the criterion this time.

## I also checked against the reference orbit

For the reference case, the orbit made by the state generator was compared, after the generation had finished, with the independent analytic reference.

The difference in the final time was about 1.5 x 10^-8.

The difference in the accumulated phase was about 8.4 x 10^-11.

The positional orbits overlap as well.

Also, the three conserved relational quantities placed inside did not drift during the computation.

Here too, however, what matters is not only that the numbers agreed well.

That the answer-checking data were not returned to the generator.

That the update formula was not switched from case to case.

That the necessary physical conditions were not passed as external arguments every time.

I attach importance to the agreement having been obtained while keeping these three.

## One failure, left as it is

Of the five cases, the slowest orbit took about 3817 revolutions.

The state of the orbit itself could be generated to the end.

But because the time was stored as an exponential, over a long duration it exceeded the limit of the numerical representation.

I did not switch the clock alone to another method, or renormalize it midway.

For that case I recorded, as it is,

"the orbital state could be generated, but the present way of reading out the clock reached the limit of finite precision".

This is a part I value quite highly in this research.

Do not erase what did not work by changing the rules afterwards.

## I have not "completely discarded" background spacetime and particles either

This is easily misunderstood, so let me write it clearly.

This time, to make the orbit on the answer-checking side, I used the known theory of gravity and electromagnetism.

Therefore, it is not that

"gravity and the Coulomb force were derived from zero, without general relativity or electromagnetism".

Rather, I tested

whether, placing known physics as the answer for the moment, the side that regenerates that motion can make the next state from the present state alone, without referring to an external orbit solution, the labels of particles A/B, or a branch by force type such as gravity / electromagnetism.

This is quite the reverse of the previous direction, "can physics be derived from a few relations?".

First, decompose known physics into a form closed by states alone.

After that, consider which states are really necessary, and whether they can be reduced further.

This time I took that order.

## What I am not saying

Let me state clearly, this time as well, what these results do not say.

Gravity and electromagnetism have not been unified.

The Einstein-Maxwell equations have not been derived from eight states.

There is no proof that eight states is the minimal number of states in nature.

Nor am I claiming that the two charge relational quantities are fundamental degrees of freedom of nature.

General eccentric orbits, spin, and mergers in strong gravitational fields have not been reproduced either.

What was treated this time is a limited model using quasi-circular orbits, leading-order conservative dynamics, gravitational quadrupole radiation and electromagnetic dipole radiation.

Nor have square roots, reciprocals and exponential functions all been derived from multiplication alone.

The "closed state update" of this time means that the necessary physical information was not hidden outside, and was written as one synchronous update.

## What was found

What was found in the fifth thought experiment was

that even if the numbers are right, if necessary information is hidden outside, one cannot say "generated from states alone".

So I flushed out the information being held to decide the next state, and brought it back into explicit states, including the intermediate calculations.

In the sixth thought experiment, without changing those audit rules, I increased the interactions from one kind to two.

Gravity and the Coulomb force.

Gravitational-wave radiation and electromagnetic dipole radiation.

Even so, with one and the same state update rule, it could be moved to five kinds of orbit with rather different properties.

Joined into one line, it comes out like this.

Last time, it turned out that letting the action change with the state allows the move to a family of gravitational orbits.

This time, I first flushed out whether any information was hidden in that implementation, and on that basis tested whether the known dynamics including gravity and the Coulomb force could be put into the same closed state update.

This is not yet the stage of "having derived physics".

Rather, it is the stage of decomposing known physics once more into states and relations, and fixing, one by one,

what is really necessary,

what was merely convenient external information,

and how far it can be carried by the same update rule.

Next: of these eight states, how many are really independent?

And can the form of the interaction itself, given this time from the known theory, be derived from fewer relations?

I return there.

## Paper information

This article is based on the two papers of the fifth and sixth thought experiments.

### The fifth thought experiment

"The Fifth Thought Experiment: Was There No Hidden State?
- A Self-Audit of the Two-State Orbit Generation of the Fourth Thought Experiment, Tracking the State Closure from Five States to the Sixth Persistent State N = nu -"

Version: v2.1

Concept DOI (always to the latest version)
https://doi.org/10.5281/zenodo.22973086

Version DOI (fixed to v2.1)
https://doi.org/10.5281/zenodo.22973087

English PDF (direct download)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E4%BA%94%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E4%B9%97%E6%B3%95%E7%8A%B6%E6%85%8B%E5%86%85%E9%83%A8%E5%B1%95%E9%96%8B%E3%81%A8%E8%BF%91%E4%BC%BC%E6%AC%A1%E6%95%B0%E6%A4%9C%E8%A8%BC_20260924/thought_experiment_05_hidden_state_audit_en_v2.1.pdf

### The sixth thought experiment

"The Sixth Thought Experiment: Can Two Kinds of Long-Range Interaction Be Internalized into a Single Closed State Map?
- Reconstructing Gravity, Coulomb, GW-Quadrupole and EM-Dipole Radiation of a Charged Quasi-Circular Binary with Eight Persistent States, and Verifying the Transferability of the Same Transition over Five Conditions -"

Version: v1.1

Concept DOI (always to the latest version)
https://doi.org/10.5281/zenodo.22974631

Version DOI (fixed to v1.1)
https://doi.org/10.5281/zenodo.22974632

English PDF (direct download)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%85%AD%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E8%8D%B7%E9%9B%BB8%E7%8A%B6%E6%85%8B_%E9%87%8D%E5%8A%9B%E3%82%AF%E3%83%BC%E3%83%AD%E3%83%B3%E8%A4%87%E6%95%B0%E7%9B%B8%E4%BA%92%E4%BD%9C%E7%94%A8_20260926/thought_experiment_06_charged8_state_closure_en_v1.1.pdf

### The previous article

"Repeating the Same Action Was Not Enough - When I Let the Action Itself Change with the State, Nine Gravitational Orbits All Fitted, the Fourth Thought Experiment"
https://note.com/kiharanoriaki/n/n1a56045ab5c3

### The article before that

"I Tried to Read Out Two Bodies, and the Observer Appeared as a Third State - reading the mass ratio, the distance and the time from an orbit and a gravitational wave, the third thought experiment"
https://note.com/kiharanoriaki/n/n7ab90a782c0f

### The second thought experiment

"In this model, attraction and repulsion came out of one and the same rule - deriving a Coulomb-type force while keeping two numbers, the second thought experiment"
https://note.com/kiharanoriaki/n/n1d87bdbfa230

### The design document of the research project

"Design the Search Before Searching for the Answer - rebuilding the search for the foundations of physics as a single project"
https://note.com/kiharanoriaki/n/n4c66305cb6cb

#TheoreticalPhysics #Gravity #Electromagnetism #CoulombForce #GeneralRelativity #DiscreteDynamics #Simulation #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
