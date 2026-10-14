import Mathlib

/-!
Literal edits of an infinite tagged sequence (the report's Section 4.1, Step 2).
An original entry's tag is its original index t ≥ 0; insertion number A (counting from 0)
receives tag -(A+1), so the j-th insertion gets tag -j. Output symbols are assigned
separately, so inserted symbols may coincide with each other or with originals.
-/

-- This is the first file in the suggested reading order (see the guide in `Glitchy.lean`).
-- Lean ideas introduced here:
--   * `def` with `if ... then ... else ...`
--   * `inductive`: defining a property by listing the only rules that can establish it
--   * proof by `induction` over those rules
--   * the everyday tactics `intro`, `simp`, `split_ifs`, `push_cast`, `omega`

-- Everything until `end Glitchy` gets the prefix `Glitchy.`, so `insertAt` below is really
-- `Glitchy.insertAt`. This keeps our names from clashing with Mathlib's.
namespace Glitchy

-- A sequence is a function `τ : ℕ → ℤ`: position `x` holds the tag `τ x` (the report writes
-- τ(x)). The tag records where an entry came from (its provenance), not the symbol it shows.
--
-- `insertAt τ p A` is the sequence after inserting one new entry at position `p`:
--   * positions before `p` are unchanged;
--   * position `p` gets the fresh tag `-(A+1)`, where A = number of earlier insertions;
--   * positions after `p` move right by one, so position `x` now holds the old `τ (x - 1)`.
-- Its arguments are `τ`, `p`, `A`, and a position `x`; the result is an integer tag.
-- `(A : ℤ)` turns the natural number A into an integer so that it can be negated.
def insertAt (τ : ℕ → ℤ) (p A : ℕ) (x : ℕ) : ℤ :=
  if x < p then τ x else if x = p then -((A : ℤ) + 1) else τ (x - 1)

-- `deleteAt τ p` removes the entry at position `p`; everything after it moves left by one.
def deleteAt (τ : ℕ → ℤ) (p : ℕ) (x : ℕ) : ℤ :=
  if x < p then τ x else τ (x + 1)

-- `Edited A D τ` says: τ can be reached from the unedited sequence by exactly A insertions
-- and D deletions, in some order, at any positions.
--
-- `inductive Edited : ℕ → ℕ → (ℕ → ℤ) → Prop` declares a new property (`Prop`) of two
-- naturals and a sequence. Each `|` line is a rule (a "constructor"), and these rules are
-- the ONLY ways to prove an `Edited` statement:
--   * `original`: the unedited sequence x ↦ x (each tag is the original position) is
--                 `Edited 0 0`. `(x : ℤ)` casts x into ℤ.
--   * `insert`:   if `h` proves `Edited A D τ`, then inserting at any `p` gives
--                 `Edited (A + 1) D`. The new tag is -(A+1), so the tags are -1, -2, -3, ...
--   * `delete`:   if `h` proves `Edited A D τ`, then deleting at any `p` gives
--                 `Edited A (D + 1)`.
-- `{A D τ}` are implicit arguments: Lean infers them from the type of `h`.
-- Because these are the only rules, a fact about every edited sequence can be proved by
-- `induction`: check it for `original`, then show `insert` and `delete` preserve it.
inductive Edited : ℕ → ℕ → (ℕ → ℤ) → Prop
  | original : Edited 0 0 (fun x => (x : ℤ))
  | insert {A D τ} (h : Edited A D τ) (p : ℕ) :
      Edited (A + 1) D (insertAt τ p A)
  | delete {A D τ} (h : Edited A D τ) (p : ℕ) :
      Edited A (D + 1) (deleteAt τ p)

-- Statement: every tag is at least -A. Inserted tags are -1, ..., -A, and originals are ≥ 0.
-- Naming it `Edited.tag_lower` enables "dot notation": for `h : Edited A D τ`, writing
-- `h.tag_lower x` means `Edited.tag_lower h x`. Later proofs in this file use that.
/-- No tag is below minus the number of insertions. -/
theorem Edited.tag_lower {A D τ} (h : Edited A D τ) :
    ∀ x, -(A : ℤ) ≤ τ x := by
  -- Induct on how `h` was built: one case per constructor of `Edited`.
  induction h with
  -- Base case: the tag is x ≥ 0 = -0. `intro x` names the variable from `∀ x`;
  -- `simp` closes the goal. A `;` separates two tactics written on one line.
  | original => intro x; simp
  -- `@insert A D τ h p ih` names the pieces of the `insert` rule. The `@` lets us name its
  -- implicit arguments too. Here A, D, τ describe the sequence *before* this insertion, and
  -- `ih` (the induction hypothesis) is the claim for that sequence: `∀ x, -(A : ℤ) ≤ τ x`.
  -- The goal is `∀ x, -((A + 1 : ℕ) : ℤ) ≤ insertAt τ p A x`.
  | @insert A D τ h p ih =>
    intro x
    -- Replace `insertAt` by its definition (the if-then-else).
    simp only [insertAt]
    -- Split the `if`s into three goals: x < p, x = p, and x > p.
    split_ifs
    -- `·` works on one goal. `have := ih x` adds the IH at x as a fact named `this`.
    -- `push_cast` rewrites `((A + 1 : ℕ) : ℤ)` as `(A : ℤ) + 1`, and `omega` (an automatic
    -- solver for linear arithmetic over ℕ and ℤ) finishes.
    · have := ih x; push_cast; omega
    -- x = p: the tag is exactly -(A+1).
    · push_cast; omega
    -- x > p: the tag is the old τ (x - 1) ≥ -A > -(A+1).
    · have := ih (x - 1); push_cast; omega
  | @delete A D τ h p ih =>
    intro x
    simp only [deleteAt]
    -- Deleting only moves old tags around, so `ih` proves both cases.
    -- `t₁ <;> t₂` runs `t₂` on every goal that `t₁` produces.
    split_ifs <;> apply ih

-- Statement: if position x holds an original entry (tag ≥ 0), its displacement
-- s = x - τ x (current position minus original index) lies in [-D, A]: the report's
-- Eq. (shift-range). Each earlier insertion pushes the entry right by one, and each earlier
-- deletion pulls it left by one.
-- `(x : ℤ)` casts the position into ℤ so the difference may be negative.
/-- Surviving original entries have displacement in [-D,A]. -/
theorem Edited.displacement {A D τ} (h : Edited A D τ) :
    ∀ x, 0 ≤ τ x → -(D : ℤ) ≤ (x : ℤ) - τ x ∧ (x : ℤ) - τ x ≤ A := by
  induction h with
  -- Base case: the displacement is x - x = 0. `intro x hx` names both x and the hypothesis
  -- `hx : 0 ≤ τ x`.
  | original => intro x hx; simp
  | @insert A D τ h p ih =>
    intro x hx
    -- Unfold `insertAt` in the hypothesis `hx` and in the goal. `⊢` stands for the goal.
    simp only [insertAt] at hx ⊢
    -- Split into cases. `with hxp heq` names the case facts, e.g. `hxp : x < p`.
    split_ifs at hx ⊢ with hxp heq
    -- x < p: nothing moved, so the IH applies directly.
    · have := ih x hx; push_cast; omega
    -- x = p: this is the new entry. Its tag -(A+1) is negative, which contradicts `hx`.
    · omega
    -- x > p: the entry came from position x - 1 and moved right by one.
    · have := ih (x - 1) hx
      -- Subtraction on ℕ is truncated (0 - 1 = 0), so the cast of x - 1 equals
      -- (x : ℤ) - 1 only when x ≥ 1. Record that fact.
      have : 0 < x := by omega
      push_cast
      omega
  | @delete A D τ h p ih =>
    intro x hx
    simp only [deleteAt] at hx ⊢
    split_ifs at hx ⊢ with hxp
    -- x < p: unchanged.
    · have := ih x hx; push_cast; omega
    -- x ≥ p: the entry came from position x + 1 and moved left by one.
    -- `push_cast at *` tidies the casts in every hypothesis and in the goal.
    · have := ih (x + 1) hx; push_cast at *; omega

-- Statement: different positions carry different tags. `Function.Injective τ` means
-- `τ x = τ y → x = y`. This is what makes provenance well defined.
/-- Tags stay distinct, including when an inserted entry is later deleted. -/
theorem Edited.injective {A D τ} (h : Edited A D τ) : Function.Injective τ := by
  induction h with
  -- Base case: `hxy` says (x : ℤ) = (y : ℤ). `simpa using hxy` simplifies it to x = y,
  -- which is the goal.
  | original => intro x y hxy; simpa using hxy
  | @insert A D τ h p ih =>
    intro x y hxy
    -- Collect lower bounds on the old tags so `omega` can see that the new tag -(A+1)
    -- differs from every old tag, since those are all ≥ -A.
    have lx := h.tag_lower x
    have ly := h.tag_lower y
    have lx' := h.tag_lower (x - 1)
    have ly' := h.tag_lower (y - 1)
    simp only [insertAt] at hxy
    -- Each of x and y is before, at, or after p: 3 × 3 = 9 cases.
    split_ifs at hxy with hx hx' hy hy' hy hy'
    -- In every case (`all_goals`), first try "the old tags were distinct (`ih`), then
    -- omega". If that fails, because `hxy` does not compare two old tags, try omega alone.
    -- `first | t₁ | t₂` runs t₁, and runs t₂ only if t₁ fails.
    all_goals first
    | have := ih hxy; omega
    | omega
  | @delete A D τ h p ih =>
    intro x y hxy
    simp only [deleteAt] at hxy
    -- 2 × 2 = 4 cases. In each, `ih` says the two old positions are equal, and omega
    -- turns that into x = y (or a contradiction).
    split_ifs at hxy with hx hy hy
    all_goals have := ih hxy; omega

-- After `open Finset`, `univ`, `card_le_card`, ... can be written without `Finset.` in front.
open Finset

-- Statement: at most A of the first N positions hold an inserted (negative) tag.
-- Reading the statement from the inside out:
--   `(univ : Finset (Fin N))`              the positions {0, ..., N-1};
--   `.filter (fun x : Fin N => τ x < 0)`   keep the positions whose tag is negative;
--   `.card`                                count them.
-- (`τ x` with `x : Fin N` silently converts x into a natural number first.)
-- Proof idea: x ↦ τ x is injective and lands in {-A, ..., -1}, which has A elements.
/-- At most A inserted entries can survive in any prefix of length N. -/
theorem Edited.inserted_card_le {A D τ} (h : Edited A D τ) (N : ℕ) :
    ((univ : Finset (Fin N)).filter (fun x : Fin N => τ x < 0)).card ≤ A := by
  -- `let` gives the set of positions a short name.
  let S := (univ : Finset (Fin N)).filter (fun x : Fin N => τ x < 0)
  -- On the first N positions, x ↦ τ x is still injective: it is `τ` composed with the
  -- injective map `Fin.val : Fin N → ℕ`.
  have hi : Function.Injective (fun x : Fin N => τ x) :=
    h.injective.comp Fin.val_injective
  -- The tags of S lie in the integer interval [-A, -1]. (`Icc a b` = {a, ..., b}.)
  have hs : S.image (fun x : Fin N => τ x) ⊆ Icc (-(A : ℤ)) (-1) := by
    -- To prove X ⊆ Y, take an element `t` with `ht : t ∈ X`.
    intro t ht
    -- `mem_image : t ∈ S.image f ↔ ∃ x ∈ S, f x = t`, and `.mp` uses the → direction.
    -- `obtain ⟨x, hx, rfl⟩` unpacks the ∃ into x and `hx : x ∈ S`. The `rfl` pattern
    -- uses the equation f x = t to replace t by `τ x` everywhere.
    obtain ⟨x, hx, rfl⟩ := mem_image.mp ht
    -- `mem_filter : x ∈ S.filter P ↔ x ∈ S ∧ P x`; `.2` takes the second half,
    -- so `hn : τ x < 0`.
    have hn := (mem_filter.mp hx).2
    -- Membership in [-A, -1] is the pair of inequalities -A ≤ τ x and τ x ≤ -1.
    -- `⟨_, _⟩` builds that pair. `by omega` turns τ x < 0 into τ x ≤ -1.
    exact mem_Icc.mpr ⟨h.tag_lower x, by omega⟩
  -- A subset is no larger: (S.image _).card ≤ (Icc (-A) (-1)).card.
  have hc := card_le_card hs
  -- An injective image has the same size as S. `at hc` rewrites inside hc, not the goal.
  rw [card_image_of_injective _ hi] at hc
  -- `Int.card_Icc` computes the size of [-A, -1], which is A. `simpa ... using hc`
  -- simplifies hc and the goal until they match.
  simpa [Int.card_Icc] using hc

end Glitchy
