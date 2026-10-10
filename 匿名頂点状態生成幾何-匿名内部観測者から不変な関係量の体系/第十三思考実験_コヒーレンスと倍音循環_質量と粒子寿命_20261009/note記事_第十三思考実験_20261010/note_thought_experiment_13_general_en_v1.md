# Why Does a Particle That Has Acquired Mass Remain Instead of Vanishing? - The Harmonics That Seemed to Flow Out Were Coming Back: the Thirteenth Thought Experiment

There is a word that often comes up in explanations of quantum mechanics.

Decoherence.

A quantum superposition, through its involvement with the surrounding environment, comes to look like the probabilities of the everyday world. That is the usual explanation.

Recently, reading Dr. Masahiro Hotta's explanation of the many-worlds interpretation, wave-packet collapse and decoherence, a naive question came to me.

Can wave-packet collapse and decoherence really be separated as two different phenomena?

And in the first place, for whom is coherence a property?

This is the story of where that question led, which turned out to be an unexpected place.

## For whom is coherence a property?

My clue was my own experience with laser holography.

In holography, the slightest vibration suddenly disturbs the interference image.

But if the phase only shifts by a fixed amount, the fringes do not disappear; they just move sideways.

The image disappears when the phase wobbles during the exposure and gets averaged out.

Even with the same light, whether interference can be read depends on what you overlap it with and on the time window in which you record it.

If so, then coherence is not a label attached to the wave itself, but a relation between the partner wave and the conditions of readout.

That is what I thought.

## Align the harmonics and you get a sharp peak

Now consider a wave made of harmonics.

On top of a fundamental wave, we add waves with 2, 3, 4, ... times its frequency, lining up their crests.

If we add N of them, the combined wave becomes a sharp peak at one point and almost vanishes elsewhere.

The width of the peak shrinks roughly as 1/N.

[Insert image]

Figure:
`fig00_harmonic_peak_sharpening.png`

Caption:
"N harmonics added with their crests aligned. Grey is a single tone (N = 1), the same height everywhere. Blue is 4 harmonics, orange is 16. The more harmonics, the thinner and sharper the peak. The horizontal axis is one period of the fundamental."

This is a localized wave, a lump of waves that looks like a particle.

In an earlier article, I wrote that when such a localized wave is sent through a double slit, it interferes as a wave and yet shows up at one point like a particle.

The important point is that this peak is narrow not only in space but also in time.

The technique for making extremely short light pulses with lasers (mode locking) also creates a sharp peak in time by aligning the phases of many frequencies.

## Sharp waves rarely meet

From here the first answer comes into view.

If the partner of a sharp wave is a single tone, the partner has the same height everywhere, so interference can be read wherever and whenever they overlap.

But if the partner is also a wave with a sharp peak, then if the two peaks are off by even a little, the overlap almost disappears.

The window in which they can meet becomes narrow, about 1/N, in both space and time.

The same wave interferes readably with a single tone, and hardly at all with another sharp wave.

So we cannot talk about coherence without saying what the partner is.

**Coherence is a relation between the partner wave and the readout window**

This is my answer to the first question.

## Even if the background is a single tone, harmonics move into it

Next I thought about what, concretely, the environment is.

In explanations of decoherence, the environment is the rest of the system that we have decided not to observe, and it is removed in the calculation.

But to think about what is actually happening, we have to decide what kind of wave the environment is.

In the simplest case, let the environment be a single-tone wave with aligned phase.

What happens if a particle-like wave with many harmonics collides with it?

Since the single-tone wave has no harmonics, it seems that nothing would happen.

However, if there is scattering that exchanges components between the two waves, the harmonics on the particle side move over to the background.

The fact that the background was initially a single tone does not prevent the harmonics from moving.

Here a new question arose.

Even if the particle acquires something like mass through its involvement with the background, if it keeps pouring its harmonics into the background each time, won't it eventually lose its sharpness and disappear?

Why does a particle that has acquired mass remain, instead of dissipating and vanishing?

## In the July experiment, the harmonics that flowed out came back

Faced with this question, I remembered a numerical experiment I ran in July 2026.

It is an experiment in which two waves, A and B, collide again and again inside a closed box.

A starts as a spread-out simple wave. B is a sharp wave with aligned harmonics.

At each collision, part is reflected and part moves to the partner. I set the fraction that is reflected to R = 0.70.

[Insert image]

Figure:
`fig01_exchange_R070_waveform_evolution.png`

Caption:
"The shapes of the two waves A (blue) and B (orange) after repeated collisions. From top: collisions 0, 1, 2, 3, 5, 10, 20 and 42. At first B is sharp and A is broad, but with more collisions the sharpness moves to A and then back to B. At collision 42 the two have almost the same shape. Each line is scaled so that its maximum is 1."

The sharpness was not lost one way.

It moved from B to A, and back to B. It kept going back and forth.

Moreover, this calculation contains no term that lets energy escape to the outside (no radiation term).

Even if, seen from the particle side, the harmonics seemed to decrease, they had not flowed out and vanished; they had moved to the partner and were on their way back.

## This time, I checked again with the same program

For this paper, I copied the original program of that July experiment without changing a single character and ran it again.

Both the data and the figures matched the July results completely. The figures are identical pixel for pixel.

Then I recorded the state of the waves at every collision up to the 128th.

I found four things.

First. There is no dissipation. The total strength of each wave stayed at 1 no matter how many times they collided.

Second. The height of the peak changes. With the strength unchanged, when the sharpness moves to the partner, the wave spreads out and its peak drops. The peak height changed by a factor of about 30, from about 0.004 to 0.125. The figure above scales each line to a maximum of 1, so this change cannot be seen there.

Third. When the average order of the harmonics (the effective order) of A and B are added, the sum was always exactly 33. A has order 1 and B has an average order of 32, and the two only mix in a fixed proportion, so the sum does not change. The fraction that moves to the partner in the first collision is 1 − R = 0.30, which is the same as the fraction transmitted.

Fourth. The back and forth can be read as an angle. When the phase of the exchange reaches 90 degrees the two waves are split half and half, at 180 degrees they are completely swapped, and at 360 degrees they are back where they started. At collision 42 the two had almost the same shape because it is close to the half-and-half point.

At R = 0.70, however, they never return exactly. The angle advanced per collision is a ratio that does not divide 360 degrees. Up to collision 128, the closest approach to a return is collision 103, where 1.6 parts in 10,000 remained on the partner side.

[Insert image]

Figure:
`fig02_normalization_trace_R070.png`

Caption:
"The record up to collision 128. Top: the amount removed by the re-scaling of strength that the original program performs after each collision. It is only the size of rounding error, so the re-scaling does nothing. Middle: the effective orders of A and B. The sum of the two lines is always 33. Solid lines are the original program, dashed lines the calculation without re-scaling; they overlap exactly. Bottom: the absolute height of the peak. With the strength unchanged, it goes up and down by a factor of about 30."

## When R is set to a value that divides evenly, the state itself returned

So what value of R would make the waves return exactly?

In fact, the value R = 0.70 has a history.

In the experiment of July 13, when I varied the reflected fraction R coarsely from 0 to 1, the difference between the two waves was smallest at 0.7, and the exchange looked cleanest there. So I had been using 0.70.

Later, in papers of July 15 and 18, I examined the neighbourhood of 0.7 very finely. There was one deep dip around 0.6972, and when I went finer still, its centre converged exactly to the value R = cos²(23π/124) ≈ 0.6972.

This value is one of the values for which the angle advanced per collision is a ratio that divides 360 degrees exactly (a finite-order root). What matters is only that it is a ratio of integers that divides evenly; the numbers 23 and 124 themselves have no special meaning.

This time, keeping the calculation procedure of the original program as it was, I replaced only R with this value and recorded up to collision 248.

- At collision 62, the shapes of A and B were completely swapped.
- At collision 124, the state of the two waves itself returned to the initial state. The deviation from the initial state was about 5 parts in 10 to the 15th power, only the size of the computer's calculation error.
- At R = 0.70, collision 103 only brought the shapes almost back. This time, not only the shapes but the state itself, including the phases of the waves, came back.

I also ran the same calculation twice and confirmed that every file of data and figures it produced matched, without a single bit of difference.

[Insert image]

Figure:
`fig04_R12423_180deg_steps_052-072.png`

Caption:
"Shapes at collisions 52 to 72 with R = cos²(23π/124). At collision 62 in the middle, the shapes of A (blue) and B (orange) are completely swapped."

[Insert image]

Figure:
`fig05_R12423_360deg_steps_114-134.png`

Caption:
"Shapes at collisions 114 to 134 in the same run. At collision 124 in the middle, they return to the initial shapes: A a broad simple wave and B a sharp wave. The shapes before and after it repeat the sequence of collisions 0 to 10."

But the return here happened because I set R to that value. The period 124 is fixed from the start by the value of R that was put in. It does not show that the period arose from within the interaction.

## It is not stable because it has a period

When you hear that something goes back and forth and returns, you are tempted to say it is stable because it has a period.

Written as a formula, this is the condition that after n collisions it returns: U to the power n = 1.

But if you put this condition in at the start, of course it returns.

That does not explain why it is stable.

In the July experiment, and also in this time's exact return at collision 124, the reflected fraction R was a value I gave by hand. Whether it returns was decided from the start by that value.

So I reverse the order.

First there is an interaction. When exchanges are repeated many times, ask which states remain and which fade away.

If a period, or phases that line up at integer ratios (phase locking), appear, they are results.

**The period is not the cause of stability but a result of what remains stably**

This is the order I held most important in this thought experiment.

## Does the internal back and forth look like a clock running slow?

Suppose a particle exchanges harmonics with the background while continuing to appear from outside as the same sharp peak.

Even if nothing is lost overall, inside the particle the exchange with the background keeps going round and round.

Here I asked one question.

Could what looked like acquiring mass be a slowing of the clock caused by this internal back and forth?

The idea of linking matter to an internal period goes back to de Broglie a hundred years ago. De Broglie determined the internal period from the mass. My question goes the other way: can we read a mass-like property from the internal back and forth?

Moreover, in relativity, the slowing of clocks and the shortening of lengths in the direction of motion happen together.

Looking only at time would be one-sided.

I want to read, from the same state, the period in the time direction and the width of the peak in the space direction separately, and check whether they change together by the same rule.

But if the formulas of relativity are put into the calculation, that is not a derivation. I will not put them in at the start, and compare afterwards.

## Why do particles with lifetimes exist?

There is one more important point.

In the world there are not only particles that are stable forever, but also particles that break down after a certain time.

A completely closed back and forth alone does not produce a lifetime.

If there is a path by which a little leaks out of the loop of exchange, the state may be kept for a long time but eventually breaks down.

Can stable particles, long-lived but decaying particles, and quickly decaying particles be handled in the same framework? This too is a question still to be checked.

## What is the background?

Up to now I have kept writing "background".

But as I thought about it, I noticed one more strange thing.

In my way of thinking, there are no space coordinates or time coordinates at the starting point. There is only the overall state and the interaction that changes it.

Then where is the background wave? When did the collision happen? And how many background waves are there in the first place?

None of these can be decided in advance.

Two waves can be said to be separate only when their responses to the interaction differ.

So the background in this thought experiment is not a container with a place, a time and a number of waves, but the nameless part of the overall state that has not yet been distinguished by the interaction.

Thinking this way, one more question remains.

If the background and the particle are not separate things, where do the kinds of particles, such as electrons and photons, come from?

Perhaps the kinds of particles are not fixed in advance, but appear distinguishable as a result of reading out the quantities conserved by the interaction.

I have carried this question over to the next paper.

## What I am not saying

Let me make clear what this thought experiment does not say.

- I am not saying that I have derived mass, the relation between time and space in relativity, quantum energy levels or the lifetimes of real particles.
- I have not yet shown that a period (U to the power n = 1) appears as a result of the interaction. In the July experiment R was a value given by hand, and the period was fixed from the start by that value.
- There is not yet any result in which phases lining up at integer ratios, discrete energies and lifetimes until breakdown came out of the same calculation.
- What I checked this time is that the July experiment was reproduced with the same program and its contents recorded (no dissipation, strength conserved, peak height changing, the sum kept at 33), and that when R is set to a value that divides evenly, the state itself returns at collision 124.

## What I found

I found three things.

First. Coherence is not a property of the wave itself but a relation between the partner wave and the readout window. A sharp wave with aligned harmonics has a narrow window for meeting a partner, in both space and time.

Second. Even if, seen from the particle side, the harmonics seem to flow out, in the exchange between two closed waves the harmonics that flowed out come back. This time, redoing the calculation with the original program, I confirmed that there is no dissipation, that only the peak height changes while the strength is conserved, and that the sum of the effective orders is kept at exactly 33. When R is set to a value that divides evenly, the state itself returns to the start at collision 124.

Third. It is not stable because it has a period; first there is an interaction, and a period appears in what remains stably. In this order, the next task is to ask whether mass, the slowing of clocks, the shortening of lengths, lifetimes and the kinds of particles can be read from a single interaction.

If a single-tone background and a sharp harmonic wave are made to exchange without any special condition on the period, and only the states that remain stably are picked out, can their period, the slowing of their clock, the width of their peak and their lifetime all be read at the same time?

That is what I want to check by calculation.

## About the paper

This article is based on the following paper.

### The thirteenth thought experiment

"The Thirteenth Thought Experiment: Why Does a Particle That Has Acquired Mass Not Dissipate and Vanish?
- Coherence as a Relation, Harmonic Exchange of Localized Waves, Internal Clocks and Particle Lifetimes Considered from a Single Interaction -"

Version: v2.0 (2026-10-10)

Concept DOI (always to the latest version)
https://doi.org/10.5281/zenodo.23266250

Version DOI (fixed to v2.0)
https://doi.org/10.5281/zenodo.23280070

English PDF (Zenodo)
https://zenodo.org/records/23280070/files/thought_experiment_13_coherence_harmonic_exchange_mass_lifetime_en_v2.0.pdf

English PDF (direct download from GitHub)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%8D%81%E4%B8%89%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E3%82%B3%E3%83%92%E3%83%BC%E3%83%AC%E3%83%B3%E3%82%B9%E3%81%A8%E5%80%8D%E9%9F%B3%E5%BE%AA%E7%92%B0_%E8%B3%AA%E9%87%8F%E3%81%A8%E7%B2%92%E5%AD%90%E5%AF%BF%E5%91%BD_20261009/thought_experiment_13_coherence_harmonic_exchange_mass_lifetime_en_v2.0.pdf

Folder of the paper (figures, and the programs and data of the re-run of the July experiment and of the record at R = cos²(23π/124); GitHub)
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E5%8C%BF%E5%90%8D%E9%A0%82%E7%82%B9%E7%8A%B6%E6%85%8B%E7%94%9F%E6%88%90%E5%B9%BE%E4%BD%95-%E5%8C%BF%E5%90%8D%E5%86%85%E9%83%A8%E8%A6%B3%E6%B8%AC%E8%80%85%E3%81%8B%E3%82%89%E4%B8%8D%E5%A4%89%E3%81%AA%E9%96%A2%E4%BF%82%E9%87%8F%E3%81%AE%E4%BD%93%E7%B3%BB/%E7%AC%AC%E5%8D%81%E4%B8%89%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_%E3%82%B3%E3%83%92%E3%83%BC%E3%83%AC%E3%83%B3%E3%82%B9%E3%81%A8%E5%80%8D%E9%9F%B3%E5%BE%AA%E7%92%B0_%E8%B3%AA%E9%87%8F%E3%81%A8%E7%B2%92%E5%AD%90%E5%AF%BF%E5%91%BD_20261009

The Japanese full text and the package of the re-run of the July experiment (zip) are included in the same Zenodo record. The package of the record at R = cos²(23π/124) is in the GitHub folder above. The Japanese version is the original; the English version is a translation.

### The Japanese version of this article

https://note.com/kiharanoriaki/n/n4274f6d4ada7

### The article on localized waves and the double slit

"The Double-Slit Puzzle Without Observation or "Wave-Packet Collapse"?"
https://note.com/kiharanoriaki/n/n701e9d57d7bb

### The article on the July experiment

"Does a Wave Packet Collapse by Observation, or Gather Through Interaction?"
https://note.com/kiharanoriaki/n/n7c59e0a3b13d

### The previous article

"Where Does 1/137 Come From? - Rereading Time, Mass, Charge and Symmetry from the Light-Cone Equation: the Twelfth Thought Experiment"
https://note.com/kiharanoriaki/n/n936f889736eb

### The design document of the research project

"Design the Search Before Searching for the Answer - rebuilding the search for the foundations of physics as a single project"
https://note.com/kiharanoriaki/n/n4c66305cb6cb

#TheoreticalPhysics #QuantumMechanics #Decoherence #Coherence #WavePacketCollapse #Harmonics #Mass #ParticleLifetime #Relativity #ThoughtExperiment #IndependentResearch #OpenScience #Preprint
