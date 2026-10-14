import Glitchy

-- These commands report every axiom used transitively by the main results.
-- Expected: only Lean's standard propext, Classical.choice, and Quot.sound.
--
-- `#print axioms foo` lists every axiom that the proof of `foo` depends on, following all
-- the lemmas it uses, including Mathlib's. It is how you check that nothing was assumed
-- without proof. The three standard axioms are:
--   * `propext`          two logically equivalent propositions are equal;
--   * `Classical.choice` the axiom of choice (it also gives the law of excluded middle);
--   * `Quot.sound`       quotient types behave as expected (it also gives function
--                        extensionality).
-- Almost all of Mathlib relies on these three. An unfinished proof (`sorry`) would show up
-- as `sorryAx`.
-- Run this file with:  lake env lean Audit.lean
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

-- Section 4.2: actual edit histories through rounded-window QFT mass.
#print axioms Glitchy.general_k_nondividing
#print axioms Glitchy.general_k_nondividing_universal
#print axioms Glitchy.general_k_nondividing_canonical
#print axioms Glitchy.cyclic_matching_moment
#print axioms Glitchy.cyclicProbability_eq_fourier
#print axioms Glitchy.cyclicProbability_sum
#print axioms Glitchy.Edited.goodPairs_bound
#print axioms Glitchy.Edited.matching_bound
#print axioms Glitchy.Edited.normalised_matching
#print axioms Glitchy.fejerSquare_expansion
#print axioms Glitchy.fejerSquare_geometric
#print axioms Glitchy.fejer_expectation_lower
#print axioms Glitchy.fejer_mass_of_matching
#print axioms Glitchy.roundedWindows_separation
#print axioms Glitchy.roundedWindows_full_of_lt
#print axioms Glitchy.canonicalPeaks_nearest
#print axioms Glitchy.nondividing_tie_regression
#print axioms Glitchy.nondividing_mixed_regression
#print axioms Glitchy.nondividing_long_regression
