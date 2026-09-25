import NativeRepairFourCycle

/-!
# Frozen geometry and unselected geometric source attachments

A repair which changes only the load component cannot change a geometry-only
readout, regardless of the schedule or preparation. An unobserved volume
attachment is not identified by the same retained load history. The four-cycle
control shows that matching the diagonal of a covariance matrix does not select
its separated-position entries. The actual twelve-port witness is separately
verified by exact finite enumeration; no physical curvature is formalized here.
-/
namespace OPH.NativeGeometricSourceIdentification

variable {Geometry Load Move : Type*}

/-- A load-only repair preserves the supplied geometric component. -/
def frozenStep (repair : Load → Move → Load) (state : Geometry × Load)
    (move : Move) : Geometry × Load := (state.1, repair state.2 move)

theorem frozen_step_geometry (repair : Load → Move → Load)
    (state : Geometry × Load) (move : Move) :
    (frozenStep repair state move).1 = state.1 := rfl

/-- Arbitrary finite schedules, including null attempts, preserve geometry. -/
theorem frozen_history_geometry (repair : Load → Move → Load)
    (schedule : List Move) (state : Geometry × Load) :
    (schedule.foldl (frozenStep repair) state).1 = state.1 := by
  induction schedule generalizing state with
  | nil => rfl
  | cons move rest ih =>
    simp only [List.foldl_cons]
    rw [ih]
    rfl

theorem geometry_readout_history (repair : Load → Move → Load)
    (schedule : List Move) (state : Geometry × Load) (volume : Geometry → ℝ) :
    volume (schedule.foldl (frozenStep repair) state).1 = volume state.1 := by
  rw [frozen_history_geometry]

/-- A log-volume contrast with its own fixed reference is exactly zero. -/
noncomputable def logVolumeContrast (volume reference : ℝ) : ℝ :=
  (Real.log volume - Real.log reference) / 3

theorem frozen_volume_contrast (repair : Load → Move → Load)
    (schedule : List Move) (state : Geometry × Load) (volume : Geometry → ℝ) :
    logVolumeContrast (volume (schedule.foldl (frozenStep repair) state).1)
      (volume state.1) = 0 := by
  rw [geometry_readout_history repair schedule state volume]
  simp [logVolumeContrast]

/-- No preparation or cross-position weighting changes a zero field's covariance. -/
theorem zero_cross_covariance {Sample Port : Type*} [Fintype Sample]
    (weight : Sample → ℝ) (q : Sample → Port → ℝ)
    (hq : ∀ s i, q s i = 0) (i j : Port) :
    ∑ s, weight s * q s i * q s j = 0 := by
  simp [hq]

/-- Coarse-graining also preserves the zero geometric source. -/
theorem zero_coarse_readout {Port : Type*} [Fintype Port]
    (weight q : Port → ℝ) (hq : ∀ i, q i = 0) :
    ∑ i, weight i * q i = 0 := by
  simp [hq]

/-- Positive volume attachments can encode any supplied load gain. -/
theorem attached_volume_positive (gain reading : ℝ) :
    0 < Real.exp (3 * gain * reading) := Real.exp_pos _

theorem attached_volume_contrast (gain reading : ℝ) :
    logVolumeContrast (Real.exp (3 * gain * reading)) 1 = gain * reading := by
  simp [logVolumeContrast, Real.log_exp]
  ring

/-- Two retained histories that agree cannot identify differing attached fields. -/
theorem no_factor_through_identical_record {World Record : Type*}
    (record : World → Record) (target : World → ℝ) (a b : World)
    (same : record a = record b) (different : target a ≠ target b) :
    ¬ ∃ decode : Record → ℝ, ∀ w, decode (record w) = target w := by
  rintro ⟨decode, h⟩
  apply different
  rw [← h a, ← h b, same]

/-- The gain is an unobserved attachment, not a new native history measurement. -/
theorem load_history_does_not_select_volume_gain (reading : ℝ) (hne : reading ≠ 0) :
    ¬ ∃ decoded : ℝ, ∀ gain : ℝ,
      decoded = logVolumeContrast (Real.exp (3 * gain * reading)) 1 := by
  rintro ⟨decoded, h⟩
  have hzero := h 0
  have hone := h 1
  rw [attached_volume_contrast] at hzero hone
  simp only [zero_mul, one_mul] at hzero hone
  exact hne (hone.symm.trans hzero)

open Matrix OPH.Thermodynamics.NativeRepairFourCycle

/-- Covariances normalized to unit trace; their diagonal entries agree. -/
def normalizedLoad : Matrix (Fin 4) (Fin 4) ℚ := (1/3 : ℚ) • Pi

def normalizedDrive : Matrix (Fin 4) (Fin 4) ℚ := (1/24 : ℚ) • (L * L)

theorem same_diagonal_distinct_separation :
    (∀ i, normalizedLoad i i = 1/4 ∧ normalizedDrive i i = 1/4) ∧
    normalizedLoad 0 2 = -1/12 ∧ normalizedDrive 0 2 = 1/12 := by
  constructor
  · intro i
    fin_cases i <;> norm_num [normalizedLoad, normalizedDrive, Pi, L,
      Matrix.mul_apply, Fin.sum_univ_succ]
  · norm_num [normalizedLoad, normalizedDrive, Pi, L,
      Matrix.mul_apply, Fin.sum_univ_succ, Matrix.cons_val_two]

theorem normalized_source_shapes_distinct : normalizedLoad ≠ normalizedDrive := by
  intro h
  have hentry := congrArg (fun M : Matrix (Fin 4) (Fin 4) ℚ => M 0 2) h
  change normalizedLoad 0 2 = normalizedDrive 0 2 at hentry
  have hc := same_diagonal_distinct_separation
  rw [hc.2.1, hc.2.2] at hentry
  norm_num at hentry

#print axioms frozen_step_geometry
#print axioms frozen_history_geometry
#print axioms geometry_readout_history
#print axioms frozen_volume_contrast
#print axioms zero_cross_covariance
#print axioms zero_coarse_readout
#print axioms attached_volume_positive
#print axioms attached_volume_contrast
#print axioms no_factor_through_identical_record
#print axioms load_history_does_not_select_volume_gain
#print axioms same_diagonal_distinct_separation
#print axioms normalized_source_shapes_distinct
end OPH.NativeGeometricSourceIdentification
