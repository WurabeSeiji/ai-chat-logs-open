# Is the "Color" of Quarks Inside Eight Numbers? — I Tried Writing a One-Line Formula from a Public Lecture in Octonions

On November 18, 2026, I attended the lecture "How Far Can the World of Elementary Particles Be Reproduced on a Computer?" by Professor Tetsuya Onogi, given in "Science Night 2026," a public lecture series of the Graduate School of Science, The University of Osaka.

It was a talk for a general audience and easy to follow, and yet the difficult concepts were not blurred; each was explained carefully, one by one. Even for me, without professional training in physics, it was very easy to understand.

What impressed me most was that quantum chromodynamics (QCD), the theory of the strong force, can be written in just one line.

I began studying theoretical physics on my own this year, and I was in the middle of studying the octonions, a system of numbers made of eight real numbers. So I considered whether the formula I had learned in the lecture could be written in octonions as it stands. I organized the record of that attempt into the form of a paper and published it on Zenodo (Concept DOI: https://doi.org/10.5281/zenodo.24369465).

In this article, I introduce its content step by step.

## 1. The One-Line Formula from the Lecture

Protons and neutrons are made of three particles called quarks. What binds quarks together is the strong force, and the particles that carry this force are called gluons. The theory that describes this world is QCD.

QCD is written in the following single line.

　L = −(1/4g²) F·F + Σ q̄ (iγD − m) q

The first term is the gluon-only term, and the second term is the term for quarks and their interaction. I understood that the only quantities that can be adjusted in this formula are the strength of the force g and the quark masses m.

Being able to write a formula, however, is different from being able to solve it. The lecture introduced lattice QCD, a method in which spacetime is approximated by a lattice, quarks are placed on the lattice sites, gluons are placed on the lines connecting the sites, and the calculation is done on a computer.

## 2. Three Colors, Eight Gluons

Quarks have three kinds of states called "color." This has nothing to do with visible colors; the name was given because there are three of them.

The gluons that carry the force come in eight kinds. The set of transformations that mix the colors is called SU(3), and the number of its independent transformations is exactly eight.

Three and eight. These two numbers form the skeleton of QCD.

## 3. The Octonions Split into 1 + 3

A real number is one number, a complex number is a pair of two numbers, and a quaternion is a set of four numbers. Next come the octonions, sets of eight numbers. It is known that there are only these four systems of numbers in which multiplication is possible and magnitude is preserved under multiplication.

The octonions have seven units whose square is −1. Choose one of them and let it play the role of the i of the complex numbers. Then the eight numbers split as follows.

　eight numbers = one complex number + three complex numbers

This is the splitting 2 + 6 = 8. The part "three complex numbers" corresponds to the three color components of a quark. The part "one complex number" is the component without color.

## 4. The Transformations That Preserve the Multiplication Rule Form SU(3)

The transformations that do not break the multiplication rule of the octonions have 14 independent directions in total. If we collect only those that do not move the unit chosen above, eight remain.

These eight coincide exactly with SU(3), the transformations that mix the colors. In other words, **the number eight of the gluons appears inside the octonions**.

One point needs care: the eight gluons are not the eight numbers of the octonions themselves. The gluons are on the side of the transformations that act on the octonions.

This correspondence is an already known construction, shown by Günaydin and Gürsey in 1973. I learned of it afterwards and re-checked it with my own multiplication table.

## 5. Writing the Formula in Octonions

With these preparations, the formula from the lecture can be re-read as follows.

- Quark: a field with three colors → a field taking values in the "three complex numbers" part of the octonions
- Gluon: the SU(3) gauge field → a field taking values in the transformations of the octonions that do not move the chosen unit
- Sum over colors → the "one complex number" component of an octonion product

The form of the formula does not change. Its content is also the same as standard QCD. What changes is that the three of color and the eight of the gluons appear as properties of the octonions, not as numbers given from outside.

Interestingly, when two octonions are multiplied, both the inner product over colors and the antisymmetric symbol that groups the three colors appear in the result. The former is the form that makes a colorless combination from a quark and an antiquark, and the latter is the form that makes a colorless combination from three quarks.

## 6. Writing the Gluon-Only Term in Components Gives Seven Conditions

Next, I wrote down the first term of the formula, the gluon-only part, in the components of the octonions.

Apart from the chosen unit, six directions remain. The rotations among six directions are as many as the ways of choosing two directions, which is 15. The gluon field is not these 15 moving freely; it satisfies seven conditions determined by the multiplication table.

　15 − 7 = 8

The eight that remain are the eight gluons.

These seven conditions are the table of signs of the multiplication copied as it is. Geometrically, they can be read as follows.

- The six directions split into three planes. The operation of multiplying by the chosen unit rotates each of the three planes by 90°
- The eight gluons consist of two rotations within the planes and six rotations that turn two planes in a linked way

## 7. The Formula on the Lattice Can Also Be Written

The formula on the lattice introduced in the lecture (the Wilson action) can be written in octonions in the same way. The variable placed on a line connecting lattice sites is an element of SU(3), so it is read as a transformation of the octonions and rewritten as a sum over the eight components. Here too, the content is the same as the standard formula.

## 8. What Came Out of the Octonions and What Was Still Placed from Outside

After the rewriting, I sorted what appeared as a property of the octonions from what I placed from outside exactly as in standard QCD.

What appeared from the octonion side

- The three components of color
- SU(3), the transformations that mix the colors
- The eight gluons
- The two forms that make colorless combinations
- The component without color

What was still placed from outside

- Spacetime, and derivatives with respect to spacetime
- The sign that separates time from space
- Spinors, which describe how electrons and quarks turn
- The kinds of quarks and their masses
- The strength of the force
- The lattice, and the way of summing as a quantum theory

Only the part concerning color appeared from the octonion side. Everything concerning spacetime is placed from outside. The choice of one unit out of the seven is also made by hand. I do not yet understand what selects that one.

## 9. This Is a Reference for Checking Answers

This paper contains no new theorem. It does not derive QCD from the octonions. It only rewrites the content of standard QCD in the language of octonions.

I wrote it nevertheless because it is needed for the goal of my research. I am studying whether spacetime, particles, and interactions can be read out solely from the relations among states, without placing a background spacetime from the start. The Standard Model is not something to be replaced; it is the answer that must be reproduced.

As I go on with the study in octonions, I will need something to compare against when checking whether what corresponds to the strong force has come out. This rewriting is that **reference for checking answers**. I summarized the points of comparison in five items.

- Whether the three components of color and the component without color appear in the 1 + 3 splitting
- Whether the transformations that carry the interaction appear as eight transformations satisfying the seven conditions
- Whether the colorless combinations appear in the form of the inner product and the antisymmetric symbol
- Whether the gluon-only term appears in the form written down here
- Whether the adjustable quantities are only the strength of the force and the quark masses

If points of disagreement appear, I will record them as well.

The formulas in the paper were checked numerically in 25 items with the attached program.

I have only just begun studying this subject. I would be grateful if you could point out any errors or offer advice.

Finally, I would like to express my sincere gratitude to Professor Onogi for a lecture that was both easy to follow and sincere.

### Paper

"QCD Expressed in Octonions ― Writing the Lagrangian of Quantum Chromodynamics As It Stands in Octonions, as a Reference for Checking Answers ―"

Concept DOI (always the latest version):
https://doi.org/10.5281/zenodo.24369465

Version DOI (v1.0):
https://doi.org/10.5281/zenodo.24369466

English PDF (Zenodo, direct download):
https://zenodo.org/records/24369466/files/qcd_expressed_in_octonions_en.pdf

Paper folder (manuscripts and verification program, GitHub):
https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/QCD%E3%82%92%E5%85%AB%E5%85%83%E6%95%B0%E3%81%A7%E8%A1%A8%E7%8F%BE%E3%81%97%E3%81%9F%E3%82%89_20261119

### Previous Article

"Which Axes Are Spacetime and Which Are Internal? — Two Radii Carried by All Waves and the Covariance of Conservation Laws"
https://note.com/kiharanoriaki/n/nbd90fac6a6f2

#TheoreticalPhysics #ParticlePhysics #QuantumChromodynamics #QCD #Quarks #Gluons #Octonions #LatticeQCD #IndependentResearch #OpenScience #Zenodo
