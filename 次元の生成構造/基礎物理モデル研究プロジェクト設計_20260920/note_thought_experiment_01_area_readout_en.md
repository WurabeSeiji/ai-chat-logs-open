# Two Numbers for State and Relation, and Elliptic and Hyperbolic Orbits Appeared

-- A report on the first thought experiment carried out under the design document

Noriaki Kihara

## This continues the previous article

In the previous article I wrote about how I rebuilt my search for the foundations of physics as a single project. The design document starts from nothing but nameless "states" and the "relations" between them, with no spacetime and no particles given in advance.

This article reports the first thought experiment carried out under that design.

The result first. "State" and "relation", with neither space nor time defined, are represented by two numbers a and b, drawn as coordinates on a two-dimensional XY plane. On that plane I placed the simplest possible rule of interaction. Elliptic orbits and hyperbolic orbits appeared. I also confirmed that the phase along the orbit becomes discrete, meaning that it advances in jumps.

## The starting point is two numbers

The reason for two is simple. To express "state" and "relation" as distinct things, one number is not enough. At least two are needed.

Write the two numbers as a and b and draw them as coordinates on a two-dimensional XY plane. The pair (a, b) is one point on the plane. Which of the two is the state and which is the relation is left undecided, because nothing inside the system can decide it. The two numbers have no names, and swapping them counts as the same thing.

This plane is not space. It is just a picture of two numbers side by side. There is no time yet, either.

## The rule that nothing can be read alone

The design document had one more principle. A state cannot be read by itself. Reading anything out requires an interaction.

Put into mathematics, this becomes a strong constraint. Think of a point on the plane as an arrow from the origin. Nothing can be read from one arrow. Something can be read only when there are two. There is essentially only one way of reading that satisfies this condition, and it is the area of the parallelogram spanned by the two arrows.

The length of an arrow, or the angle between two arrows, cannot be read at this stage. Measuring lengths and angles needs a ruler, and nobody has supplied one yet. Everything starts from the fact that only area can be read.

## Then the simplest rule

Next comes the rule of interaction. At every step the two numbers are mixed in fixed proportions. This is a linear interaction, the simplest rule one can think of.

And, as the design document requires, the quantity that can be read must be conserved. Here that means the area spanned by the arrow at one step and the arrow at the next step does not change from step to step.

With these two in place, the content of the rule collapses into a single number. Mixing two numbers should have four degrees of freedom, yet only one of them means anything when seen from inside.

## Elliptic and hyperbolic orbits appeared

The size of that one number decides how the point moves. There are three kinds of motion, shown side by side in Figure 1.

====================
[Insert Figure 1 here]
File: fig01_area_readout_three_classes_en_v1.png
Caption: Figure 1. Three kinds of motion from one rule. Top left: an elliptic orbit. Top right: hyperbolic orbits. Bottom left: the boundary case, a straight line in equal steps. Bottom right: the area readout (orange line) does not change from step to step.
====================

In one range the point keeps circling the center and traces an ellipse (Figure 1, top left). The area swept per step is constant (Figure 1, bottom right), and the acceleration always points to the center. Constant swept area is the same form as Kepler's second law, in which a planet sweeps equal areas in equal times around the Sun. A point pulled toward a center, tracing an ellipse while keeping the swept area constant, is the same picture as an orbit in a gravitational field.

In another range the point recedes exponentially and traces a hyperbola (Figure 1, top right). This is the same shape as the path drawn in spacetime by a body that keeps a constant acceleration in relativity.

On the boundary between the two ranges the point moves in a straight line in equal steps (Figure 1, bottom left). This is the form of inertial motion, with no acceleration.

Neither space nor time was given, and yet the basic forms of motion, which are rotation, uniform acceleration, and inertia, all came out of one rule. Complex numbers were not put in at the start either. On the elliptic side, something that does the job of the imaginary unit appears from the interaction itself.

In this system the first thing that can be read about motion is acceleration, not velocity. Velocity changes sign depending on which way the direction of time is taken, so it cannot be decided from inside. Acceleration can be read regardless of direction. And there is only one acceleration. Whether you call it "the orbit is curving" or "the point is being pulled to the center", it is the same readout. This has the same form as Einstein's equivalence principle, that acceleration and gravity cannot be told apart from inside. The two settings in which Einstein considered the equivalence principle, a rotating system and a uniformly accelerated system, correspond to the elliptic side and the hyperbolic side of this one rule.

The quantity that fixes the size of an orbit turns out to be exactly the radius of curvature of the orbit. The ruler that measures it was not supplied from outside. The interaction decides the ruler.

## The phase also turned out to be discrete

The motion along the ellipse does not slide continuously. It advances by a fixed angle at each step. Each of the dots on the ellipse in Figure 1 (top left) is the position at one step. Like a clock hand jumping from mark to mark, the position along the orbit, that is, the phase, moves in jumps.

How many steps it takes to come back to the starting position is decided by that angle. If the angle is an awkward value, the point never comes back exactly.

The design document has another idea that matters here. Every readout has a limit to how finely it can distinguish things. If there is such a limit, then there is a first step count at which the point has come back so close to the start that it cannot be told apart from it. So an integer appears.

This integer cannot take arbitrary values. It takes only special numbers decided by the angle. If the angle is set by the golden ratio, for example, the numbers that appear are 3, 5, 8, 13, 21, and so on, the Fibonacci numbers. The larger the orbit, the larger the integer (Figure 2). The integer that says after how many steps the orbit closes turned out to be the size of the orbit measured in units of the finest distinguishable difference.

====================
[Insert Figure 2 here]
File: fig02_area_readout_closure_staircase_en_v1.png
Caption: Figure 2. Number of steps needed to come back to the starting position. Horizontal axis: the size of the orbit divided by the finest distinguishable difference. Vertical axis: number of steps. The four lines are four different step angles, and the blue one is the golden-ratio case. In every case the number rises like a staircase and takes only special integers.
====================

There was no need to decide in advance that the orbit closes after so many steps. Given a continuous quantity and a finite fineness, the discrete integers come out on their own.

## The failed check was left as it is

In the previous article I wrote down prohibitions for myself. This time I put them into practice.

For every numerical check, the condition "if this happens, it fails" was declared before the run. Of fifteen checks, one actually failed. That result was left in the paper as it is, and I looked for the cause. It was loss of significance in the arithmetic. From that cause I made a new prediction, declared it in advance as well, tested it, and it passed.

For drafting, for checking the formulas and the numerics, and for review, I used several AIs in separate roles. Along the way an error was found in one proposition. Both the side that wrote it and the side that checked it had missed it, and the third, reviewing side found it. This history is also recorded in the work log of the paper. Preserve not only what was discovered, but why the computation was performed and where it went wrong. This became the first working example of the discipline I wrote about last time.

## What was found

Represent state and relation, with neither space nor time defined, by two numbers. Draw them as a point on a two-dimensional plane. Place a linear interaction and conservation on it. That alone produced elliptic and hyperbolic orbits, and a discrete phase.

The mathematics used is classical, a structure well known in geometry. What is new is the route. Starting from "nothing can be read alone" and "there are no names", there was nowhere else to arrive but this structure.

On to the next thought experiment.

## The paper

- Title: The First Thought Experiment: Tracing the Equivalence Principle Down to Two Nameless Degrees of Freedom -- Deriving the Sign of the Sum of Squares, the Curvature Radius, and Three Systems (Inertial, Uniformly Accelerated, Rotating) from the Conservation of an Area Readout (v1.0, Japanese and English versions)
- Concept DOI (always resolves to the latest version): https://doi.org/10.5281/zenodo.22857952
- Version DOI (v1.0, fixed): https://doi.org/10.5281/zenodo.22857953
- Paper PDF (English, direct download): https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/次元の生成構造/基礎物理モデル研究プロジェクト設計_20260920/thought_experiment_01_equivalence_principle_area_readout_en_v1.pdf
- Previous article (Japanese): https://note.com/kiharanoriaki/n/nbd87af574beb

---

#physics #gravity #equivalenceprinciple #relativity #thoughtexperiment #keplerslaws #fibonacci #spacetime #theoreticalphysics #independentresearch #openscience #preprint
