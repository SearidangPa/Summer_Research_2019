import Mathlib

/-!
Literal edits of an infinite tagged sequence. Original tags are nonnegative
integers; insertion number A receives tag -(A+1). Output symbols are assigned
separately, so inserted symbols may coincide with each other or with originals.
-/
namespace Glitchy

def insertAt (s : ℕ → ℤ) (p A : ℕ) (x : ℕ) : ℤ :=
  if x < p then s x else if x = p then -((A : ℤ) + 1) else s (x - 1)

def deleteAt (s : ℕ → ℤ) (p : ℕ) (x : ℕ) : ℤ :=
  if x < p then s x else s (x + 1)

inductive Edited : ℕ → ℕ → (ℕ → ℤ) → Prop
  | original : Edited 0 0 (fun x => (x : ℤ))
  | insert {A D s} (h : Edited A D s) (p : ℕ) :
      Edited (A + 1) D (insertAt s p A)
  | delete {A D s} (h : Edited A D s) (p : ℕ) :
      Edited A (D + 1) (deleteAt s p)

/-- No tag is below minus the number of insertions. -/
theorem Edited.tag_lower {A D s} (h : Edited A D s) :
    ∀ x, -(A : ℤ) ≤ s x := by
  induction h with
  | original => intro x; simp
  | @insert A D s h p ih =>
    intro x
    simp only [insertAt]
    split_ifs
    · have := ih x; push_cast; omega
    · push_cast; omega
    · have := ih (x - 1); push_cast; omega
  | @delete A D s h p ih =>
    intro x
    simp only [deleteAt]
    split_ifs <;> apply ih

/-- Surviving original entries have displacement in [-D,A]. -/
theorem Edited.displacement {A D s} (h : Edited A D s) :
    ∀ x, 0 ≤ s x → -(D : ℤ) ≤ (x : ℤ) - s x ∧ (x : ℤ) - s x ≤ A := by
  induction h with
  | original => intro x hx; simp
  | @insert A D s h p ih =>
    intro x hx
    simp only [insertAt] at hx ⊢
    split_ifs at hx ⊢ with hxp heq
    · have := ih x hx; push_cast; omega
    · omega
    · have := ih (x - 1) hx
      have : 0 < x := by omega
      push_cast
      omega
  | @delete A D s h p ih =>
    intro x hx
    simp only [deleteAt] at hx ⊢
    split_ifs at hx ⊢ with hxp
    · have := ih x hx; push_cast; omega
    · have := ih (x + 1) hx; push_cast at *; omega

/-- Tags stay distinct, including when an inserted entry is later deleted. -/
theorem Edited.injective {A D s} (h : Edited A D s) : Function.Injective s := by
  induction h with
  | original => intro x y hxy; simpa using hxy
  | @insert A D s h p ih =>
    intro x y hxy
    have lx := h.tag_lower x
    have ly := h.tag_lower y
    have lx' := h.tag_lower (x - 1)
    have ly' := h.tag_lower (y - 1)
    simp only [insertAt] at hxy
    split_ifs at hxy with hx hx' hy hy' hy hy'
    all_goals first
    | have := ih hxy; omega
    | omega
  | @delete A D s h p ih =>
    intro x y hxy
    simp only [deleteAt] at hxy
    split_ifs at hxy with hx hy hy
    all_goals have := ih hxy; omega

open Finset

/-- At most A inserted entries can survive in any prefix of length N. -/
theorem Edited.inserted_card_le {A D s} (h : Edited A D s) (N : ℕ) :
    ((univ : Finset (Fin N)).filter (fun x : Fin N => s x < 0)).card ≤ A := by
  let S := (univ : Finset (Fin N)).filter (fun x : Fin N => s x < 0)
  have hi : Function.Injective (fun x : Fin N => s x) :=
    h.injective.comp Fin.val_injective
  have hs : S.image (fun x : Fin N => s x) ⊆ Icc (-(A : ℤ)) (-1) := by
    intro t ht
    obtain ⟨x, hx, rfl⟩ := mem_image.mp ht
    have hn := (mem_filter.mp hx).2
    exact mem_Icc.mpr ⟨h.tag_lower x, by omega⟩
  have hc := card_le_card hs
  rw [card_image_of_injective _ hi] at hc
  simpa [Int.card_Icc] using hc

end Glitchy
