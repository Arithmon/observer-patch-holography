import Mathlib

/-!
# Component-local source programs

The model's algebra, response, overlap and entropy construction is analytic.
This file proves the all-program information obstruction, including local
adaptation, retained records, arbitrary local randomness tapes and stopping
rules. `State` includes every local controller and checkpoint field. The
one-step interface excludes an oracle that reads another component to pick
the receiver's move; such an oracle is a communication channel.
-/

set_option autoImplicit false

namespace OPH.SourceSelectionLocality

variable {Site State Output : Type*}

/-- A component may use its own complete state and the public step number.
Any independent randomness tape can be included in its initial state. -/
def localRun (step : ℕ → Site → State → State) (initial : Site → State) :
    ℕ → Site → State
  | 0 => initial
  | n + 1 => fun i => step n i (localRun step initial n i)

theorem localRun_eq (step : ℕ → Site → State → State)
    (x y : Site → State) (receiver : Site) (h : x receiver = y receiver)
    (n : ℕ) : localRun step x n receiver = localRun step y n receiver := by
  induction n with
  | zero => exact h
  | succ n ih => simp only [localRun, ih]

theorem whole_transcript_eq (step : ℕ → Site → State → State)
    (x y : Site → State) (receiver : Site) (h : x receiver = y receiver) :
    (fun n => localRun step x n receiver) =
      (fun n => localRun step y n receiver) := by
  funext n
  exact localRun_eq step x y receiver h n

/-- Even a function of the entire infinite local transcript cannot distinguish
the interventions. Measurability is unnecessary for this equality theorem. -/
theorem any_transcript_decoder_eq (step : ℕ → Site → State → State)
    (x y : Site → State) (receiver : Site) (h : x receiver = y receiver)
    (decode : (ℕ → State) → Output) :
    decode (fun n => localRun step x n receiver) =
      decode (fun n => localRun step y n receiver) := by
  rw [whole_transcript_eq step x y receiver h]

/-- A local stopping rule, including one that never stops, cannot provide
different answers in locally indistinguishable source worlds. -/
theorem stopping_read_eq (step : ℕ → Site → State → State)
    (x y : Site → State) (receiver : Site) (h : x receiver = y receiver)
    (read : (ℕ → State) → Option (ℕ × Output)) :
    read (fun n => localRun step x n receiver) =
      read (fun n => localRun step y n receiver) :=
  any_transcript_decoder_eq step x y receiver h read

theorem cannot_read_both (step : ℕ → Site → State → State)
    (x y : Site → State) (receiver : Site) (h : x receiver = y receiver)
    (decode : (ℕ → State) → Output) (a b : Output) (hab : a ≠ b) :
    ¬ (decode (fun n => localRun step x n receiver) = a ∧
       decode (fun n => localRun step y n receiver) = b) := by
  rintro ⟨ha, hb⟩
  apply hab
  rw [← ha, ← hb]
  exact any_transcript_decoder_eq step x y receiver h decode

/-- A separately supplied cross-component copy is a concrete deletion control:
it breaks the factorization used above and really communicates. -/
def copyAcross (source receiver : Site) [DecidableEq Site]
    (x : Site → State) : Site → State :=
  Function.update x receiver (x source)

theorem copyAcross_reads [DecidableEq Site] (source receiver : Site)
    (x : Site → State) : copyAcross source receiver x receiver = x source := by
  simp [copyAcross]

def population (n : ℕ) : ℕ := 10 * 4 ^ n + 2

theorem population_zero : population 0 = 12 := by norm_num [population]

theorem population_recurrence (n : ℕ) :
    population (n + 1) = 4 * population n - 6 := by
  dsimp [population]
  rw [pow_succ]
  omega

theorem population_successor_not_cube (n q : ℕ) : population (n + 1) ≠ q ^ 3 := by
  intro h
  have hm := congrArg (fun x : ℕ => x % 4) h
  have hq : q % 4 = 0 ∨ q % 4 = 1 ∨ q % 4 = 2 ∨ q % 4 = 3 := by omega
  rcases hq with hq | hq | hq | hq <;>
    norm_num [population, pow_succ, Nat.add_mod, Nat.mul_mod, hq] at hm

theorem population_not_cube (n q : ℕ) : population n ≠ q ^ 3 := by
  cases n with
  | zero =>
    rw [population_zero]
    by_cases hq : q ≤ 2
    · interval_cases q <;> norm_num
    · have hq3 : 3 ≤ q := by omega
      have hp := Nat.pow_le_pow_left hq3 3
      norm_num at hp
      omega
  | succ n => exact population_successor_not_cube n q

/-- Support of a supplied operational classical channel, not just an
algebraic restriction. A decoder must be right on every possible output. -/
def Separates (C : State → Output → Prop) : Prop :=
  ∀ x x' y, C x y → C x' y → x = x'

def Decodes (C : State → Output → Prop) (d : Output → State) : Prop :=
  ∀ x y, C x y → d y = x

theorem decoder_implies_separation (C : State → Output → Prop)
    (d : Output → State) (h : Decodes C d) : Separates C := by
  intro x x' y hxy hx'y
  exact (h x y hxy).symm.trans (h x' y hx'y)

theorem separation_implies_decoder [Nonempty State]
    (C : State → Output → Prop) (h : Separates C) : ∃ d, Decodes C d := by
  classical
  let d : Output → State := fun y =>
    if hy : ∃ x, C x y then Classical.choose hy else Classical.choice inferInstance
  refine ⟨d, ?_⟩
  intro x y hxy
  have hy : ∃ z, C z y := ⟨x, hxy⟩
  dsimp [d]
  rw [dif_pos hy]
  exact h _ _ y (Classical.choose_spec hy) hxy

theorem exact_decoder_iff [Nonempty State] (C : State → Output → Prop) :
    (∃ d, Decodes C d) ↔ Separates C := by
  constructor
  · rintro ⟨d, hd⟩
    exact decoder_implies_separation C d hd
  · exact separation_implies_decoder C

theorem collision_excludes_decoder (C : State → Output → Prop)
    (x x' : State) (y : Output) (hne : x ≠ x') (hxy : C x y) (hx'y : C x' y) :
    ¬ ∃ d, Decodes C d := by
  rintro ⟨d, hd⟩
  exact hne (decoder_implies_separation C d hd x x' y hxy hx'y)

theorem separated_channels_compose {Next : Type*}
    (C : State → Output → Prop) (B : Output → Next → Prop)
    (hC : Separates C) (hB : Separates B) :
    Separates (fun x z => ∃ y, C x y ∧ B y z) := by
  rintro x x' z ⟨y, hxy, hyz⟩ ⟨y', hx'y', hy'z⟩
  have hyy' : y = y' := hB y y' z hyz hy'z
  subst y'
  exact hC x x' y hxy hx'y'

end OPH.SourceSelectionLocality
