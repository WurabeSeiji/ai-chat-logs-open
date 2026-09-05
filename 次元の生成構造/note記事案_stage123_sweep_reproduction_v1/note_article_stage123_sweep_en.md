# The Dynamics That Drives Inflation Is Made of Three Parts

In my July article "Why Directions Stop at Three, and Waves Grow on Their Own," I reported a strange phenomenon. A seed too small to be measured — an energy of about one part in a trillion, times another trillion, times another million — amplified itself exponentially with nothing added from outside, showing a rise reminiscent of the inflation of the early universe.

In the follow-up "Inflation Was Already Written into the Toppling Wave," I studied how this amplification arises, revising and retracting several claims as I closed in on the core of the mechanism.

This is the sequel. There was a single question.

What, exactly, is the dynamics that drives this amplification made of?

## What I Did

Change the program one line at a time, and run it again. I did only that, exhaustively.

The dynamics decomposes into three parts:

- Part 1: the basic device that rotates the whole collection of waves, and the way time is stepped
- Part 2: a normalization that, when building the interaction, discards the magnitude of each wave and uses only its direction
- Part 3: an operation that extracts only a particular half of the interaction, making the rotation a real rotation

I then reran, from exactly the same initial state, dynamics with each part removed one at a time. The initial state is made fully deterministic from a random-seed formula, and I verified every time that it is bit-for-bit identical to the initial values used in the July computation.

## What I Found

- Amplification occurs only when Part 2 and Part 3 are both present.
- Remove Part 2 alone, and the amplification is replaced by a fake at the very first step: noise far larger than the seed is injected at step one, and the seed's amplification can no longer be observed.
- Remove Part 3 alone, and the same happens: noise jumps in at step one, and the system breaks in a different way.
- I also tried the shortcut of doing Part 2 just once at the beginning. This breaks the most violently of all. The normalization is a rule applied at every step, not a preprocessing.
- Running the full three-part dynamics across system sizes N=3 to 40 — 228 runs in total — the first step stayed at the seed scale in every run, and in most of them the same shape of amplification curve as in July appeared. This phenomenon is not a coincidence of one particular system size.

## Seeing It in Figures

Here are the amplification curves for all 38 systems, N=3 to 40. The vertical axis is the fraction of energy outside the initial set of waves (which I call the parent); the bottom is the seed scale, the top is saturation. In every panel a straight uphill line — exponential amplification that keeps growing at a constant rate — appears, starting from the seed.

![Amplification curves for all systems N=3..40](../電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/results/fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png)

The next three figures look at the same runs through a different window. Each of the 780 waves is plotted on the complex plane as an arrow with a direction and a magnitude.

At the start. All waves are neatly aligned into four bundles. This is the shape of the parent, corresponding to the starting point of the amplification curve (the state with only the seed).

![Complex-plane figure at the start](../電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/fig_complex_plane_step0_N3_N40_stage123.png)

At the end. The four bundles are gone, and the waves become a ring scattered in all directions. This is the state after the amplification curve has finished climbing and saturated. What I want you to notice is that the waves on the ring have almost exactly the same magnitude. In the measurements, not a single wave in any system had a magnitude below half the median. The directions scatter apart, while the magnitudes line up neatly. This is the hallmark of the final state of this dynamics.

![Complex-plane figure at the end](../電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/fig_complex_plane_final_N3_N40_stage123.png)

This is the end-state ring at maximum magnification. Note that the window of this figure is not the origin, but a tiny region — about one ten-thousandth of the ring radius — centered on one cluster on the ring. The tick labels may make the points look scattered at small values, but every point shown here has the same magnitude as the ring. What to read here is that the waves gathered in the same direction do not fall exactly onto a single point but form a loose bundle of very slight width. This is the fine organization of the final state that the amplification curve does not show.

![Zoom of the end state](../電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/fig_complex_plane_final_zoom_N3_N40_stage123.png)

## The Real Protagonist Here Is the Method of Verification

This paper contains no physical interpretation. What it records is only what is needed for anyone to reach the same result.

- That the rerun of the July computation is bit-for-bit identical to the stored values in every row
- That the regenerated initial values are bit-for-bit identical to the July initial values
- That the starting state of all 228 runs is bit-for-bit identical to those initial values

Only when this three-stage cross-check passes does the comparison of adding and removing parts become meaningful. Every file carries a SHA256 fingerprint ledger, and there is a correspondence table linking equations to program line numbers.

## About the Paper

The paper is for researchers and is quite hard. It is split into an overview and three chapters. It centers on the equation–program correspondence, the verification procedure, and the ledger of all data; do not expect it to be an enjoyable read.

The paper package (English versions, programs, and all data bundled; about 710MB) is on Zenodo:
https://doi.org/10.5281/zenodo.22317636

The English PDFs can also be downloaded directly here:

Overview (structure of the whole and the inventory of contents)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/paper_overview/overview_stage123_sweep_en.pdf

Chapter 1: Generation of the initial data
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/paper_ch1_static_parents/ch1_static_parents_en.pdf

Chapter 2: Definition of the dynamics and the sweep body (main chapter)
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/paper_ch2_sweep_dynamics/ch2_sweep_stage123_en.pdf

Chapter 3: Complex-plane readout figures
https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/paper_ch3_complex_plane/ch3_complex_plane_en.pdf

The original July paper (Onset and Threefold Classification of Outcomes of Spontaneous Splitting)
https://doi.org/10.5281/zenodo.21486234
