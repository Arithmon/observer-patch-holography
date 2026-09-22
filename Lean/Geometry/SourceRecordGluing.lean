import Geometry.SourceNetCausalCone
import Mathlib

/-!
# The proposed record-gluing law and its finite consequences

The extra source law is stated in `code/source_selection_model/RECORD_GLUING.md`.
It is a hypothesis about physical source implementations, not a new Lean axiom.
The results below prove the geometry of its primitive process theory and the
all-word value/intervention/resource consequences of primitive lowering.
They do not assert that the current simulator satisfies that source law.
-/

set_option autoImplicit false

namespace OPH.SourceRecordGluing

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

structure Leg (E : Type*) where
  displacement : E
  duration : ℝ

def Admissible (c : ℝ) (leg : Leg E) : Prop :=
  0 ≤ leg.duration ∧ ‖leg.displacement‖ ≤ c * leg.duration

/-- The actual primitive grammar: a scaled unit flight or a local wait.
The cone inequality is not used to define these primitives. -/
def Primitive (c : ℝ) (leg : Leg E) : Prop :=
  (∃ (u : E) (a : ℝ), ‖u‖ = 1 ∧ 0 ≤ a ∧
    leg.displacement = a • u ∧ leg.duration = a / c) ∨
  (leg.displacement = 0 ∧ 0 ≤ leg.duration)

theorem primitive_admissible {c : ℝ} (hc : 0 < c) (leg : Leg E)
    (h : Primitive c leg) : Admissible c leg := by
  rcases h with ⟨u, a, hu, ha, hd, ht⟩ | ⟨hd, ht⟩
  · constructor
    · rw [ht]
      exact div_nonneg ha hc.le
    · rw [hd, ht, norm_smul, Real.norm_of_nonneg ha, hu, mul_one,
        mul_div_cancel₀ _ hc.ne']
  · constructor
    · exact ht
    · rw [hd, norm_zero]
      exact mul_nonneg hc.le ht

def displacement (word : List (Leg E)) : E := (word.map Leg.displacement).sum
def duration (word : List (Leg E)) : ℝ := (word.map Leg.duration).sum

theorem rotated_scaled_flight (R : E ≃ₗᵢ[ℝ] E) (e : E) (he : ‖e‖ = 1)
    {scale tau : ℝ} (hs : 0 ≤ scale) (ht : 0 < tau) :
    ‖scale • R e‖ = (1 / tau) * (scale * tau) := by
  rw [norm_smul, Real.norm_of_nonneg hs, R.norm_map, he]
  field_simp

omit [NormedSpace ℝ E] in
theorem word_duration_nonnegative {c : ℝ} (word : List (Leg E))
    (h : ∀ leg ∈ word, Admissible c leg) : 0 ≤ duration word := by
  induction word with
  | nil => simp [duration]
  | cons leg rest ih =>
    have hl := (h leg (by simp)).1
    have hr := ih (fun x hx => h x (by simp [hx]))
    simpa [duration] using add_nonneg hl hr

omit [NormedSpace ℝ E] in
theorem word_speed_bound {c : ℝ} (word : List (Leg E))
    (h : ∀ leg ∈ word, Admissible c leg) :
    ‖displacement word‖ ≤ c * duration word := by
  induction word with
  | nil => simp [displacement, duration]
  | cons leg rest ih =>
    have hl := (h leg (by simp)).2
    have hr := ih (fun x hx => h x (by simp [hx]))
    calc
      ‖displacement (leg :: rest)‖ = ‖leg.displacement + displacement rest‖ := by
        simp [displacement]
      _ ≤ ‖leg.displacement‖ + ‖displacement rest‖ := norm_add_le _ _
      _ ≤ c * leg.duration + c * duration rest := add_le_add hl hr
      _ = c * duration (leg :: rest) := by simp [duration]; ring

/-- Normalizing a nonzero displacement constructs its unit flight control. -/
theorem normalized_control (v : E) (hv : v ≠ 0) :
    ‖(1 / ‖v‖) • v‖ = 1 ∧ ‖v‖ • ((1 / ‖v‖) • v) = v := by
  have hn : ‖v‖ ≠ 0 := norm_ne_zero_iff.mpr hv
  constructor
  · rw [norm_smul, Real.norm_of_nonneg (by positivity)]
    exact one_div_mul_cancel hn
  · rw [smul_smul]
    simp [hn]

theorem flight_wait_realizes {c T : ℝ} (hc : 0 < c) (v : E)
    (hv : ‖v‖ ≤ c * T) :
    ∃ word : List (Leg E), word.length ≤ 2 ∧
      (∀ leg ∈ word, Primitive c leg) ∧
      displacement word = v ∧ duration word = T := by
  let flight : Leg E := ⟨v, ‖v‖ / c⟩
  let wait : Leg E := ⟨0, T - ‖v‖ / c⟩
  have hf : Primitive c flight := by
    by_cases hz : v = 0
    · right
      simp [flight, hz]
    · left
      obtain ⟨hu, hd⟩ := normalized_control v hz
      exact ⟨(1 / ‖v‖) • v, ‖v‖, hu, norm_nonneg v, hd.symm, rfl⟩
  have hw : Primitive c wait := by
    right
    exact ⟨rfl, sub_nonneg.mpr ((div_le_iff₀ hc).mpr (by simpa [mul_comm] using hv))⟩
  refine ⟨[flight, wait], by simp, ?_, ?_, ?_⟩
  · intro leg hl
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hl
    rcases hl with rfl | rfl
    · exact hf
    · exact hw
  · simp [displacement, flight, wait]
  · simp [duration, flight, wait]

theorem record_cone_iff {c T : ℝ} (hc : 0 < c) (v : E) :
    (∃ word : List (Leg E), (∀ leg ∈ word, Primitive c leg) ∧
       displacement word = v ∧ duration word = T) ↔
      0 ≤ T ∧ ‖v‖ ≤ c * T := by
  constructor
  · rintro ⟨word, hw, rfl, rfl⟩
    have ha := fun leg hleg => primitive_admissible hc leg (hw leg hleg)
    exact ⟨word_duration_nonnegative word ha, word_speed_bound word ha⟩
  · rintro ⟨hT, hv⟩
    obtain ⟨word, _, hw, hd, ht⟩ := flight_wait_realizes hc v hv
    exact ⟨word, hw, hd, ht⟩

/-- Access is defined by primitive execution, not by a metric-ball menu. -/
def BudgetRead (c T : ℝ) (x y : E) : Prop :=
  ∃ word : List (Leg E), (∀ leg ∈ word, Primitive c leg) ∧
    displacement word = y - x ∧ duration word = T

theorem budget_read_iff {c T : ℝ} (hc : 0 < c) (hT : 0 ≤ T) (x y : E) :
    BudgetRead c T x y ↔ ‖y - x‖ ≤ c * T := by
  exact (record_cone_iff hc (y - x)).trans (and_iff_right hT)

def GeneratedPath (S : Set E) (c T : ℝ) (k : ℕ) (x y : E) : Prop :=
  ∃ p : ℕ → E, p 0 = x ∧ p k = y ∧
    (∀ j ≤ k, p j ∈ S) ∧ (∀ j < k, BudgetRead c T (p j) (p (j + 1)))

theorem generated_path_iff {S : Set E} {c T : ℝ} (hc : 0 < c) (hT : 0 ≤ T)
    (k : ℕ) (x y : E) :
    GeneratedPath S c T k x y ↔
      OPH.SourceNetCausalCone.Reachable S (c * T) k x y := by
  constructor <;> rintro ⟨p, h0, hk, hS, hp⟩ <;>
    refine ⟨p, h0, hk, hS, ?_⟩
  · exact fun j hj => (budget_read_iff hc hT _ _).mp (hp j hj)
  · exact fun j hj => (budget_read_iff hc hT _ _).mpr (hp j hj)

/-- The covering theorem is joined to actual primitive executions.
The all-neighbour menu is a conclusion of the flight/wait grammar. -/
theorem transport_generates_covering_cone {Ω S : Set E} {h c T : ℝ}
    (hΩ : Convex ℝ Ω) (hS : S ⊆ Ω)
    (hcover : OPH.SourceNetCausalCone.Covers Ω S h)
    (hh : 0 ≤ h) (hc : 0 < c) (hT : 0 ≤ T)
    {x y : E} (hx : x ∈ S) (hy : y ∈ S) {k : ℕ} (hk : 0 < k) :
    (‖y - x‖ ≤ (k : ℝ) * (c * T - 2 * h) → GeneratedPath S c T k x y) ∧
    (GeneratedPath S c T k x y → ‖y - x‖ ≤ (k : ℝ) * (c * T)) := by
  constructor
  · intro hd
    exact (generated_path_iff hc hT k x y).mpr
      (OPH.SourceNetCausalCone.covering_constructs_path hΩ hS hcover hh hx hy hk hd)
  · intro hp
    exact OPH.SourceNetCausalCone.reachable_outer ((generated_path_iff hc hT k x y).mp hp)

omit [NormedSpace ℝ E] in
/-- A native implementation with nonnegative extra delay cannot outrun
the ideal primitive dependency path. -/
theorem positive_overhead_preserves_cone {c T overhead : ℝ}
    (hc : 0 ≤ c) (he : 0 ≤ overhead) (v : E) (hv : ‖v‖ ≤ c * T) :
    ‖v‖ ≤ c * (T + overhead) := by nlinarith

omit [NormedSpace ℝ E] in
theorem strict_slack (c : ℝ) (hc : 0 < c) (T : ℝ) (v : E)
    (hv : ‖v‖ < c * T) :
    0 < T - ‖v‖ / c ∧ ∀ overhead,
      overhead < T - ‖v‖ / c → ‖v‖ / c + overhead < T := by
  constructor
  · exact sub_pos.mpr ((div_lt_iff₀ hc).mpr (by simpa [mul_comm] using hv))
  · intro overhead h
    linarith

/-- Homogeneity and constant latency on the response orbit force a radial
latency. The source Lie theorem supplies the transitive SO(3) orbit. -/
theorem invariant_latency_is_euclidean (F : E → ℝ) (e : E)
    (hzero : F 0 = 0)
    (hscale : ∀ (a : ℝ), 0 ≤ a → ∀ v, F (a • v) = a * F v)
    (horbit : ∀ u, ‖u‖ = 1 → F u = F e) (v : E) :
    F v = ‖v‖ * F e := by
  by_cases hv : v = 0
  · subst v
    simp [hzero]
  · obtain ⟨hu, hmul⟩ := normalized_control v hv
    calc
      F v = F (‖v‖ • ((1 / ‖v‖) • v)) := congrArg F hmul.symm
      _ = ‖v‖ * F ((1 / ‖v‖) • v) := hscale _ (norm_nonneg v) _
      _ = ‖v‖ * F e := by rw [horbit _ hu]

section Lowering

variable {Instruction NativeInstruction AbstractState NativeState : Type*}

def run {I S : Type*} (step : I → S → S) : List I → S → S
  | [], state => state
  | instruction :: rest, state => run step rest (step instruction state)

theorem run_append {I S : Type*} (step : I → S → S) (a b : List I) (s : S) :
    run step (a ++ b) s = run step b (run step a s) := by
  induction a generalizing s with
  | nil => rfl
  | cons i a ih => exact ih (step i s)

/-- Only primitive correctness is assumed. All-word correctness, including
every intervention in the input state, is a conclusion. -/
theorem primitive_lowering_extends
    (abstractStep : Instruction → AbstractState → AbstractState)
    (nativeStep : NativeInstruction → NativeState → NativeState)
    (lower : Instruction → List NativeInstruction) (project : NativeState → AbstractState)
    (primitive : ∀ i s, project (run nativeStep (lower i) s) = abstractStep i (project s))
    (program : List Instruction) (s : NativeState) :
    project (run nativeStep (program.flatMap lower) s) =
      run abstractStep program (project s) := by
  induction program generalizing s with
  | nil => rfl
  | cons i program ih =>
    simp only [List.flatMap_cons, run_append, run]
    rw [ih, primitive]

theorem lowering_retains_event_count (lower : Instruction → List NativeInstruction)
    (program : List Instruction) :
    (program.flatMap lower).length = (program.map (fun i => (lower i).length)).sum := by
  induction program with
  | nil => rfl
  | cons i program ih => simp [ih]

end Lowering

/-- Finite source preparations use q^3 independently addressable records. -/
theorem grid_population (q : ℕ) : Fintype.card (Fin 3 → Fin q) = q ^ 3 := by
  simp

omit [NormedSpace ℝ E] in
/-- A uniform perturbation of the event assignment changes the causal
margin by at most 2*spatialError + 2*c*timeError. -/
theorem causal_margin_stability {c spaceError timeError t s t' s' : ℝ}
    (hc : 0 ≤ c) (x y x' y' : E)
    (hx : ‖x' - x‖ ≤ spaceError) (hy : ‖y' - y‖ ≤ spaceError)
    (ht : |t' - t| ≤ timeError) (hs : |s' - s| ≤ timeError) :
    |(c * (s' - t') - ‖y' - x'‖) - (c * (s - t) - ‖y - x‖)| ≤
      2 * spaceError + 2 * c * timeError := by
  have hd : |‖y' - x'‖ - ‖y - x‖| ≤ 2 * spaceError := by
    calc
      _ ≤ ‖(y' - x') - (y - x)‖ := abs_norm_sub_norm_le _ _
      _ = ‖(y' - y) - (x' - x)‖ := by congr 1; abel
      _ ≤ ‖y' - y‖ + ‖x' - x‖ := norm_sub_le _ _
      _ ≤ _ := by linarith
  have htime : |(s' - t') - (s - t)| ≤ 2 * timeError := by
    calc
      _ = |(s' - s) - (t' - t)| := by congr 1; ring
      _ ≤ |s' - s| + |t' - t| := abs_sub _ _
      _ ≤ _ := by linarith
  calc
    _ = |c * ((s' - t') - (s - t)) - (‖y' - x'‖ - ‖y - x‖)| := by congr 1; ring
    _ ≤ |c * ((s' - t') - (s - t))| + |‖y' - x'‖ - ‖y - x‖| := abs_sub _ _
    _ = c * |(s' - t') - (s - t)| + |‖y' - x'‖ - ‖y - x‖| := by rw [abs_mul, abs_of_nonneg hc]
    _ ≤ _ := by nlinarith

end OPH.SourceRecordGluing
