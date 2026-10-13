import Glitchy

-- These commands report every axiom used transitively by the main results.
-- Expected: only Lean's standard propext, Classical.choice, and Quot.sound.
#print axioms Glitchy.general_k_glitches
#print axioms Glitchy.general_k_sharp_bounds
#print axioms Glitchy.general_k_neighbourhood
#print axioms Glitchy.general_k_no_surviving_insertions
#print axioms Glitchy.general_k_few_edits
#print axioms Glitchy.residue_mass
#print axioms Glitchy.Edited.displacement
#print axioms Glitchy.Edited.inserted_card_le
#print axioms Glitchy.modelI_bound
#print axioms Glitchy.modelI_sharp
#print axioms Glitchy.example_counts_correct
#print axioms Glitchy.example_fourier_mass
#print axioms Glitchy.example_bounds
#print axioms Glitchy.two_glitch_constant
#print axioms Glitchy.one_outcome_cancels
#print axioms Glitchy.insertion_budget_bound
#print axioms Glitchy.retained_normalized
#print axioms Glitchy.retained_factor_tendsto
#print axioms Glitchy.no_uniform_exact_bound
