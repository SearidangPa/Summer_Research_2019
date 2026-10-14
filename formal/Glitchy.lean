-- This is the root file of the library. `import Glitchy` (as `Audit.lean` does) loads all of it.
-- `import Glitchy.Edits` loads the file `Glitchy/Edits.lean` and everything it proves.
import Glitchy.Edits
import Glitchy.Boxes
import Glitchy.Fourier
import Glitchy.Counting
import Glitchy.Literal
import Glitchy.ModelI
import Glitchy.Examples
import Glitchy.Nondividing
import Glitchy.NondividingExamples

/-! Kernel-checkable general-k total-mass arguments in Sections 4.1 and 4.2. -/

/-!
# Reading these files if you are new to Lean

## Suggested reading order

1. `Glitchy/Edits.lean`    – what an edited sequence is (insertions and deletions)
2. `Glitchy/Counting.lean` – pure inequalities (Cauchy–Schwarz plus algebra)
3. `Glitchy/Boxes.lean`    – how many (symbol, residue) boxes the first N entries occupy
4. `Glitchy/Fourier.lean`  – characters mod r, orthogonality, and Parseval's identity
5. `Glitchy/Literal.lean`  – joins the pieces into the main theorems
6. `Glitchy/ModelI.lean`   – the simpler Model I
7. `Glitchy/Examples.lean` – the worked example from the figure, checked by computation

Section 4.2 continues with:

8. `Glitchy/Matching.lean` – circular Fourier moments and the QFT reindexing bridge
9. `Glitchy/PairCounting.lean` – matching pairs counted by induction on literal edits
10. `Glitchy/Fejer.lean` – expectation inequalities and the two trigonometric regimes
11. `Glitchy/FejerKernel.lean` – squared character sums, expansion, and concentration
12. `Glitchy/PeakGeometry.lean` – arbitrary rounded peaks and circular window separation
13. `Glitchy/Nondividing.lean` – the complete nondivisible theorem and corollaries
14. `Glitchy/NondividingExamples.lean` – rounding ties, mixed edits, and boundary regressions

Each Lean idea is explained in detail the first time it appears in this order and
more briefly after that.

## How to explore interactively (Neovim)

Install the lean.nvim plugin (https://github.com/Julian/lean.nvim). It starts Lean's
language server and adds an *infoview* window. On Neovim 0.12 or later:

    vim.pack.add { "https://github.com/Julian/lean.nvim" }
    vim.g.lean_config = { mappings = true }

lean.nvim does not install Lean itself; install `elan` as described in `README.md`.

The infoview opens beside any `.lean` file. Put the cursor at the end of a proof line and
the infoview shows the current **goal**: the facts you have (above `⊢`) and what is left
to prove (after `⊢`). Watching the goal change line by line is the fastest way to learn.

* `K` shows the type and documentation of the name under the cursor.
* `<C-]>` (or your go-to-definition key) jumps to its definition, including into Mathlib.
  `<C-o>` jumps back.
* `<LocalLeader>i` toggles the infoview, `<LocalLeader><Tab>` moves into it and back,
  and `<LocalLeader>r` restarts Lean for the current file.
* `<LocalLeader>\` on a symbol such as ℝ shows how to type it.

Typing symbols in insert mode: `\N` ℕ, `\Z` ℤ, `\R` ℝ, `\C` ℂ, `\to` →, `\le` ≤, `\ne` ≠,
`\sum` ∑, `\and` ∧, `\forall` ∀, `\exists` ∃, `\<` ⟨, `\>` ⟩, `\l` ←, `\sub` ⊆, `\mid` ∣.

## Comments

* `-- ...` runs to the end of the line.
* `/- ... -/` is a block comment.
* `/-- ... -/` documents the declaration right after it (press `K` on its name to see it).
* `/-! ... -/` documents the whole file or section.

## Declarations

* `def f (x : ℕ) : ℤ := body` defines a function `f` from ℕ to ℤ.
* `theorem name (hypotheses) : statement := proof` states and proves a result.
* `:= by` followed by tactics builds the proof step by step (see the cheat sheet below).
  Without `by`, the proof is written directly as an expression ("term mode").

## Propositions are types; proofs are values

A statement such as `0 < N` is a type, and a proof of it is a value of that type.
A hypothesis `(hN : 0 < N)` is an ordinary argument, so a theorem is a function from
proofs of its hypotheses to a proof of its conclusion. That is why theorems are
*applied* like functions: `foo N hN` uses theorem `foo` with `N` and the proof `hN`.
The same arrow `→` means both "function from" and "implies".
By convention, hypothesis names start with `h` (`hN`, `hrN`, `hsum`, ...).

## Brackets around arguments

* `(x : ℕ)` — explicit: you pass it yourself.
* `{x : ℕ}` — implicit: Lean works it out from the other arguments.
* `[Fintype B]` — instance: Lean finds it automatically. Here it means "B is finite".

## Notation used throughout

* `ℕ` (naturals, includes 0), `ℤ`, `ℝ`, `ℂ`. `ℕ → ℤ` is the type of functions from ℕ to ℤ,
  so an infinite sequence of integers is a function `s : ℕ → ℤ`.
* `fun x => e` is the function x ↦ e. Application is juxtaposition: `f x`, `f x y`.
* `(x : ℤ)` for `x : ℕ` casts x into ℤ; Lean displays such a cast as `↑x`.
* `Fin N` is {0, 1, ..., N-1}; for `x : Fin N`, `x.val` is the underlying natural number.
* `ZMod r` is the integers mod r.
* `Finset X` is a finite set of elements of X. `univ` is the set of *all* elements of a
  finite type. `S.filter P` keeps the elements satisfying P, `S.image f` is f(S), and
  `S.card` is the number of elements.
* `∑ b, f b` sums over every b; `∑ x ∈ S, f x` sums over the finite set S.
* `∧` and, `∨` or, `¬` not, `↔` iff, `∀ x, P x` for all, `∃ x, P x` there exists.
* `h.1` and `h.2` are the two halves of `h : A ∧ B`.
* `⟨a, b⟩` builds a pair, a proof of `A ∧ B`, or a proof of `∃ x, P x` (witness first).
* For `h : A ↔ B`, `h.mp : A → B` and `h.mpr : B → A`.
* `h₁.trans h₂` chains `a ≤ b` and `b ≤ c` into `a ≤ c` (also works for `=` and `⊆`).
* `Finset.card_le_card` is the lemma `card_le_card` in namespace `Finset`; after
  `open Finset` it can be written `card_le_card`.

## Tactic cheat sheet

* `intro x` — goal `∀ x, P x` (or `A → B`): name the variable (or hypothesis) and move it
  above the line.
* `exact e` — finish the goal with `e`, which must prove exactly the goal.
* `apply e` — use a lemma whose conclusion matches the goal; its hypotheses become the
  new goals.
* `refine e` — like `exact`, but each `?_` hole in `e` becomes a new goal.
* `have h : P := proof` — prove an intermediate fact `h`. `have := ...` names it `this`.
* `let x := e` — a local abbreviation.
* `rw [eq]` — rewrite left side → right side of `eq` in the goal. `rw [← eq]` goes
  right → left, and `rw [eq] at h` rewrites inside hypothesis `h`.
* `simp` — simplify with Mathlib's standard rewrite rules. `simp [l₁, l₂]` adds rules,
  `simp only [l₁, l₂]` uses only those, and `simpa [..] using e` simplifies both the goal
  and `e`, then closes the goal with `e`.
* `simp_rw [l]` — rewrite even inside `∑` and `fun` (where `rw` cannot reach).
* `unfold f` — replace `f` by its definition.
* `constructor` — split a goal `A ∧ B` into goals `A` and `B`.
* `·` — focus on the next goal (the indented lines below it prove that goal).
* `obtain ⟨x, hx⟩ := h` — unpack `h : ∃ x, P x` (or `A ∧ B`) into its pieces.
* `by_cases h : P` — split into the case `P` and the case `¬P`.
* `by_contra h` — assume the negation of the goal and prove `False`.
* `split_ifs` — split an `if ... then ... else ...` into one case per branch.
* `induction h with | c₁ ... => ... | c₂ ... => ...` — one case per constructor.
* `calc` — a chain `a = b`, `_ ≤ c`, ...; each `_` stands for the previous right side.
* `congr 1` — reduce `f a = f b` to `a = b`.
* `ext x` — prove two sets or functions are equal by checking each element `x`.
* `omega` — linear arithmetic over ℕ and ℤ (handles ℕ's truncated subtraction).
* `linarith` — linear arithmetic over ℝ using the hypotheses.
* `nlinarith` — `linarith` that also multiplies hypotheses together (some nonlinear goals).
* `positivity` — proves goals of the form `0 < e`, `0 ≤ e`, or `e ≠ 0`.
* `norm_num` — arithmetic with explicit numbers.
* `ring` — identities in commutative rings. `ring_nf` puts both sides in a normal form.
* `field_simp` — clear denominators (needs to know they are nonzero).
* `push_cast`, `norm_cast`, `exact_mod_cast` — move casts such as ℕ → ℝ around.
* `decide` — evaluate a finite, computable statement; Lean's kernel re-checks the result.
* `t₁ <;> t₂` — run `t₂` on every goal that `t₁` produces.
-/
