# Hiding a Coin with a Cup: Can It Be Written with Nothing but Overlapping Waves? — Reproducing Professor Hotta's "Coin and Cup" Without Placing Spacetime

In an article on note, "The Linear Algebra of Quantum Mechanics Is Merely a Notation for Probability" (in Japanese), Professor Masahiro Hotta places a single coin, marked "+" on the front and "−" on the back, inside an opaque cup. While the cup covers it, the coin stays front-up with probability q and flips over with probability 1 − q. Then the cup is removed, and one looks at whether the coin shows front or back.

From this classical stochastic process, the article builds, step by step, the density matrix, the Bloch representation, CPTP maps, projection matrices, the Born rule, and state collapse, and shows that the linear algebra used in quantum mechanics can be introduced as a notation for probability.

Professor Hotta's article (in Japanese):
https://note.com/quantumuniverse/n/n897cf1fca3e9

I wondered whether this thought experiment could be reproduced even without assuming a background spacetime.

The point was to realize two states — the coin hidden by the cup, and the coin visible after the cup is lifted — using only the interference of waves, with no background spacetime assumed.

I recorded that thought experiment in the order in which it was carried out, arranged it as a paper, and published it on Zenodo (Concept DOI: https://doi.org/10.5281/zenodo.24562867). The paper contains hypotheses, and the verification is still to come. This article introduces its contents in order.

This paper does not object to the conclusion of Professor Hotta's article, and it does not use the derivations in the article. What I borrowed is only the setting of a cup and a coin, which anyone can picture.

## 1. These Are the Two States I Want to Reproduce

[Figure 1 here: note_fig1_schematic_covered_en.png]

Figure 1. The cup placed over the coin. The coin is hidden and cannot be seen.

[Figure 2 here: note_fig2_schematic_beside_en.png]

Figure 2. The cup placed beside the coin. The coin is visible.

The table is drawn in black, the coin in red, and the cup in green. The colors are the same in the figures that follow.

Ordinarily we think that a place called the table exists first, and that the coin and the cup are put on it. Here, no place is put first. Both the coin and the cup are represented as waves. And I want the difference between Figure 1 and Figure 2 to come out of nothing but the relation between the two waves.

## 2. The First Method: Make Both the Coin and the Cup a Peak

At first I used the method I have used before. On a circle S¹ with radius of curvature R, harmonics are added together to make a wave that rises at only one place. Harmonics are the waves that go around the circle exactly once, twice, three times, and so on. When they are added, a peak forms at one place.

I took a small peak as the coin and a large peak as the cup.

[Figure 3 here: note_fig3_peak_cup_over_coin_en.png]

Figure 3. The coin (red, small peak) and the cup (green, large peak) placed on the circle. The cup is placed over the coin.

What is used to read the relation between the two waves is the cross term. If two quantities a and b are added and squared,

　(a + b)² = a² + b² + 2ab

Here a² is the term of the coin alone, and b² is the term of the cup alone. Only the last one, 2ab, tells us about the partner. This is the cross term: the product of the two waves. It stands only where both are present.

At first sight, this worked. A wave of the cross term stands only when the coin and the cup overlap.

[Figure 4 here: note_fig4_peak_cross_overlap_en.png]

Figure 4. The cup placed over the coin. A peak stands in the cross term (blue).

When they do not overlap, it stays zero.

[Figure 5 here: note_fig5_peak_cross_apart_en.png]

Figure 5. The cup shifted away from the coin. The cross term (blue) stays almost zero.

## 3. But Then I Was Stuck

The state in which the coin and the cup overlap can be told apart from the state in which they do not.

However, I could not tell apart the state in which the coin is not hidden (that is, visible) from the state in which the coin is hidden (that is, not visible).

All the cross term says is that the two overlap. It does not say that the cup covers the coin. A product has the same value when the order is exchanged. Even when I swapped the sizes of the coin and the cup entirely, the cross term returned the same value. Which one hides and which one is hidden is written nowhere.

At first I thought as follows. The system is a circle S¹ inside two dimensions, so there is no front and back. That must be why the relation of hiding cannot be expressed. So I considered moving to a three-dimensional spherical shell S³ inside four dimensions. I also considered placing a separate observer to decide what is in front and what is behind.

All of these go in the direction of adding something.

## 4. After a Night's Sleep

When I woke up the next morning, a different solution had come to me.

**Instead of a wave that represents the cup, why not make a wave out of the places where the cup is not?**

Where the cup covers, the height of the wave is zero. Everywhere else, it has a constant height. As a wave shape, this is what is called a well-shaped waveform.

It is more complicated than a simple peak of harmonics, but it can be realized by thinking of it as a wave with a DC component (a constant height applied uniformly around the whole circle). Any wave shape can be made by adding harmonics. The shape of a well is one of them.

The coin remains the same peak as before. The multiplication is also the same as before. The only thing changed is the wave shape of the cup.

## 5. Hidden When Placed Over, Visible When Shifted

[Figure 6 here: note_fig6_well_hidden_en.png]

Figure 6. The cup (green, well-shaped) placed over the coin (red). No peak stands in the cross term (blue). The coin is hidden.

Where the coin is, the height of the cup's wave is zero. Because it is a product, the cross term is also zero. In the calculation, the cross term at the center of the coin was about 3% of its value when the coin is not covered. It is not exactly zero because the harmonics are cut off at a finite number; it gets smaller as more harmonics are added.

[Figure 7 here: note_fig7_well_partly_visible_en.png]

Figure 7. The cup shifted a little. The rim of the cup lies over the coin. A peak of the cross term stands only on the side that is not covered. The coin is partly visible.

[Figure 8 here: note_fig8_well_visible_en.png]

Figure 8. The cup shifted off the coin. The peak of the coin stands as it is in the cross term. The coin is visible.

Hidden when placed over. Partly visible when the rim lies over it. Visible when shifted away.

The difference between Figure 1 and Figure 2 came out of two waves and their product alone. There was no need to add a three-dimensional spherical shell, an observer, or a special operation whose result depends on the order of multiplication. The circle was enough.

## 6. Why This Works

In the first method, the coin and the cup were both the same kind of peak. If both are peaks, the difference between the one that hides and the one that is hidden is not in the wave shapes. So it does not come out of the multiplication either.

The moment the cup was made well-shaped, the difference in roles entered the shape of the wave.

I noticed one more thing. Outside the well, the well-shaped wave has a constant height all around the circle. This plays the role of the place where the coin can be seen. In other words, without preparing a separate background called the table, the wave of the cup itself serves that role too.

A cup is a place that makes the coin not visible; everywhere else, the coin is visible as it is. The well-shaped wave is this statement put directly into a shape.

## 7. Not "Is It There or Not" but "Is It Visible or Not"

When I thought of the well-shaped cup, I reconsidered the question itself.

What had been on my mind from the start was the knowledge before the cup is lifted. Before lifting, we do not know how many coins are inside, whether front or back, or even whether there is anything at all. This does not mean that there are zero coins. It is a state in which there is no information in the first place. Borrowing a word from programming, I decided to call it null.

Null cannot be read directly, because there is no way to confirm that nothing is there.

So I stopped judging whether the coin is there or not. Instead, I judge whether the coin is visible or not.

- Visible: the height of the cross term is larger than the limit of what can be read (the resolution)
- Not visible: the height of the cross term is at or below the limit of what can be read

In doing so, I gave up one thing: distinguishing "not visible" from "not there".

In fact, the cross term at the place of the coin when the cup is placed over it, and the cross term at a place where there was never a coin, were about the same size. Both are below the limit of what can be read. Whether a coin is hidden under the cup, or there is no coin there, cannot be told from the cross term.

What we learn when the cup is placed over the coin is only that the coin is no longer visible. Whether there is a coin inside the cup, we do not know.

This is exactly our own state before the cup is lifted.

## 8. Being Able to Move the Cup

There is one more thing I learned in this thought experiment.

To put the cup down and to lift it, one must be able to move the cup while it remains the same cup.

On a circle, shifting a peak or a well takes nothing more than rotating the position of each harmonic. Neither the shape nor the size changes.

But when I tried to do the same on a straight line with ends, it could not be done. On a line with both ends fixed, a peak made by adding waves is pulled toward an end, and when one tries to move it to another place, both its shape and its size change. The same cup cannot be carried to another place.

For the very familiar act of putting a cup down and lifting it to hold, the system that carries it must be finite, closed, and without ends. In the paper, I go on from here to consider how far the shape of the system is narrowed down. That part is a hypothesis, and its verification is still to come.

Then why can a cup be moved on a real table? In my view, the table is actually curved and closed. The curvature is so gentle that it merely looks flat.

## 9. What Was Checked and What Remains to Be Checked

What was checked

- With the product of a peak and a peak, the cross term stands only when they overlap, and it has the same value when the two are exchanged
- When the cup is made well-shaped, the coin is hidden when the cup is placed over it, partly visible when the rim lies over it, and visible when the cup is shifted away
- On a line with both ends fixed, a peak of the same shape cannot be made at another place

The wave shapes and the numbers were computed and checked with the attached programs. The paper has 17 figures.

What remains to be checked

- Whether the act of putting a cup down and lifting it may be read as a wave
- Whether it may be said that the condition for observation to hold is that the cross term can be read
- Whether this experiment on the circle can also be drawn on the three-dimensional spherical shell S³

The parts used are all already known: the superposition of harmonics (Fourier series), the well-shaped waveform, and multiplication. There is no new theorem. The content this time is that, when known parts are put side by side, the cup and the coin can be written in this way.

## 10. Written in the Order in Which I Thought

This paper is written not from the conclusion but in the order in which I thought. The place where I got stuck and the idea that came after a night's sleep are written as they happened.

The thought experiment was carried out as a dialogue with AI (ChatGPT and Claude). Along the way, misreadings occurred many times. They are also recorded as they are, at the places where they occurred, because readers are likely to misread at the same places.

I am deeply grateful to Professor Hotta for explaining how the formalism of quantum mechanics comes about with the setting of a cup and a coin, which anyone can picture. The contents of this article and of the paper are my own attempt made by borrowing that setting, and any errors are entirely my responsibility.

I would be grateful to be told of any errors, and for any advice.

### Paper

"The Coin and Cup Thought Experiment ― Reverse-Deriving the Requirement Specification of a System from the Conditions under Which an Action Holds ―"

Concept DOI (always the latest version):
https://doi.org/10.5281/zenodo.24562867

Version DOI (v1.0):
https://doi.org/10.5281/zenodo.24562868

English PDF (Zenodo, direct download):
https://zenodo.org/records/24562868/files/coin_and_cup_thought_experiment_en.pdf

Paper folder (manuscripts in Japanese and English, figures and the programs that make them, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E5%A0%80%E7%94%B0%E5%85%88%E7%94%9F%E3%81%AEnote%E3%81%A8%E3%82%B3%E3%82%A4%E3%83%B3%E3%81%A8%E3%82%B3%E3%83%83%E3%83%97%E3%81%AE%E6%80%9D%E8%80%83%E5%AE%9F%E9%A8%93_20261124

### Previous Article

"Can Physics Be Read Out from a Single Multiplication, Without Placing Spacetime or Force Fields? — The Forms of Two Equations Were Inside the Product of Octonions"
https://note.com/kiharanoriaki/n/ndf5e50ae911e

#TheoreticalPhysics #QuantumMechanics #Observation #ThoughtExperiment #CoinAndCup #Waves #Interference #FourierSeries #Spacetime #IndependentResearch #OpenScience #Zenodo
