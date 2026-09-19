# What Remains When You Shave Off the Assumptions

-- A thought experiment that reads loops, flux, and condensates from nothing but sums of squared complex numbers

Noriaki Kihara

## It starts from a formula that usually ends with "trivially true"

Square a complex number z, and you still have a complex number. Add complex numbers together, and you still have a complex number. So the sum of squares

Σ zₙ² = K

is just an ordinary sum of complex numbers. Moreover, every complex number K has a square root, so we can write K = C². In other words,

Σ zₙ² = C²

looks special but says nothing new. In a mathematics textbook the story would end with one word: trivial.

What the paper does is a thought experiment: without breaking this triviality in any way, how much physical structure can be read out of it? Not a single new equation appears. What appears is only an order of reading.

## The first step is a rereading

Give the squared quantity a name: Zₙ := zₙ². Then the formula above becomes

Σ Zₙ = K,

a perfectly ordinary sum of complex numbers. There is no need to treat C² as special "because it is a square." In the world of squared quantities, K and C² are just complex numbers.

This small rereading is the core of the whole paper.

[Insert figure 1 here: fig1_complex_closure_loop_residual_en.png -- the overall picture: the trivial identity reread as a sum of complex relations and a loop residual]

When the sum K is zero, drawing the arrows head to tail brings you back to the starting point. That is a closed loop. When it is not zero, the loop fails to close and a leftover K remains. We call it the loop residual. Close, or leave a remainder. That is all.

## With three, a loop appears for the first time

Suppose there are three complex relations whose sum is zero:

W₁ + W₂ + W₃ = 0.

Then any one of them is automatically determined by the other two. And connecting the three arrows draws a closed triangle.

What matters here is the order. We did not assume a triangle first and then derive a closing condition. The algebraic fact that the sum is zero comes first; the triangle can be read afterwards.

[Insert figure 2 here: fig2_three_body_triangular_loop_en.png -- if three complex vectors sum to zero, a closed triangle can be read afterwards]

## With four, there is a trap

With three, there are three points and three relations -- the numbers happen to coincide. With four, a complete network in which everyone connects to everyone needs six relations.

It is tempting to declare "there are four of them, so surely they are all connected" -- but that is an assumption. Doing the actual computation, with only five relations -- four outer edges plus one diagonal -- two triangular loops can close simultaneously. So the closing of loops does not imply a complete network.

In the discipline of this paper, an underivable assumption gets deleted. The complete network was discarded right here.

[Insert figure 3 here: fig3_four_vertex_triangular_loops_en.png -- with 4 vertices and 5 relations the triangular loops still close; complete connectivity is an extra assumption]

## Glue them together, and the inside disappears

Here comes the best part.

Glue two triangles together along one shared edge. Follow each loop in a consistent orientation, and the shared edge appears once with a plus sign and once with a minus sign. Add them, and they cancel exactly to zero. The internal relation disappears, and only the relations along the outer rim remain.

Reading this as a flow per face -- a flux -- we get the structure: shared internal parts cancel, and only the external boundary remains. It looks exactly like Stokes' theorem from calculus, but the theorem was not assumed; this is the result of following the additive fact that shared terms appear with opposite signs.

And when many faces are glued until the outer shell closes completely, we decided to call that body a condensate -- because from the outside it can be read as a single object.

[Insert figure 4 here: fig4_internal_cancellation_condensate_en.png -- shared-edge cancellation, the closed polyhedron, and the hierarchical growth of condensates]

For a closed polyhedron, adding up all the face residuals always gives zero. This is not an assumption: it follows identically from the structure alone -- every edge belongs to two faces with opposite orientations. And yet, even though the total is zero, each individual face residual need not be zero.

In other words, a structure is possible in which nonzero exchange keeps going on inside while the whole remains closed. From the outside, a quiet single lump; inside, exchange continues. Borrowing the vocabulary of gauge theory, a reading close to the exchange of virtual particles becomes available (though the paper states explicitly that this is an interpretation candidate, not a derivation).

Even better: the face residual is itself a complex number. So the quantity left on the outside of a condensate can be used as a relation at the next level, and the same story repeats. Relations, faces, condensates, relations between condensates. The type stays complex all the way. Climbing the hierarchy requires no new kind of number.

## The discipline of this paper: delete, don't add

When building a model, one usually keeps adding assumptions. This paper goes the other way: every time a hidden assumption was found, it was deleted.

- Complete network -> deleted (counterexample found)
- Comparison by ratio -> removed from the foundation (it smuggles in a comparison standard)
- A same-or-different detector -> removed from the foundation
- Presupposing closure at zero -> deleted
- Presupposing the triangle -> deleted (readable after the sum is zero)

[Insert figure 5 here: fig5_thought_experiment_flowchart_en.png -- at each stage, check what is being assumed, shave it off, and keep only the readouts that survive]

The question left at the end is this: not what was introduced, but what never had to be introduced.

## Something I should say honestly

After finishing v1.1 of this paper, I read a recent paper by another research group (in quantum information and holography, an example in which entanglement measurable only among three or more parts keeps growing even after the part-by-part information has saturated) -- and it stopped me cold.

My own paper also quietly carried an assumption: that pairwise relations and loop residuals suffice to read all the information of a many-body system. Nobody had pointed it out, but in v1.2 I withdrew that assumption and moved it to the list of open problems. The discipline of deleting assumptions once found was applied to the paper itself.

So I do not claim this framework is the final description of the world. The only claim is that it can serve as a minimal foundation for clarifying which assumptions are truly indispensable.

## The paper

The paper is published on Zenodo in both Japanese and English (figures and the reproduction scripts for the figures included).

- Title: Reading Relations, Loops, Flux, and Condensates from the Closure of Complex Numbers -- A Minimal-Assumption Thought Experiment That Shaves Off Axiom Candidates (v1.3)
- Concept DOI (always resolves to the latest version): https://doi.org/10.5281/zenodo.22841916
- Version DOI (v1.3, fixed): https://doi.org/10.5281/zenodo.22841917
- Zenodo record: https://zenodo.org/records/22841917
- Paper PDF (English, direct download): https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/平方数を基本量とした場合の検討/paper_complex_square_closure_flux_condensate_readout_en_v1_3.pdf
- Paper PDF (Japanese, direct download): https://github.com/WurabeSeiji/ai-chat-logs-open/raw/main/平方数を基本量とした場合の検討/paper_complex_square_closure_flux_condensate_readout_ja_v1_3.pdf

---

#physics #mathematics #complexnumbers #geometry #gaugetheory #thoughtexperiment #theoreticalphysics #mathematicalphysics #science #researchnotes #preprint #openscience
