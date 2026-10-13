import SourceLaw.HistoryRegister

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The two-point flux rule and its Euclidean time face

A two-point flux discretization gives a seam between two cells the conductance
`face measure / centre distance` and a cell the capacity `cell volume`. The scalar packet's
spatial conductances `w_u w_v / h` and capacities `w_a w_b w_c` have this form; the golden
carrier files check it against the packet data.

Applying the same rule to the time faces of the space-time product mesh, whose face measure is
the spatial cell volume `M_i`, gives the conductance `M_i / Δτ`; this is the register's time-seam
conductance `J_i / Δτ` with `J = M` (`flux_time_face`, true by definition). With space-time cell
volumes `M_i Δτ` as capacities every time seam has the
rate `2/Δτ²` (`flux_spacetime_time_rate`) and every spatial seam keeps its original rate
(`flux_spacetime_spatial_rate`). This is a candidate principle (isotropy of the Euclidean
discretization, units with `c = 1`), not a derivation of the inertia.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

variable {𝕜 : Type*} [Field 𝕜]

/-- Two-point flux conductance: face measure over centre distance. -/
def fluxCond (face dist : 𝕜) : 𝕜 := face / dist

section Rates

variable {ι : Type*}

/-- With the time coefficient set to `J = M`, the register's time-seam conductance is the
two-point-flux value `M_i / Δτ` (face measure `M_i`, distance `Δτ`); this holds by definition. -/
theorem flux_time_face (M : ι → 𝕜) (dτ : 𝕜) {E : Type*} (cE : E → 𝕜) (T : ℕ)
    (i : ι) (t : Fin T) :
    regCond cE M dτ T (.inr (i, t)) = fluxCond (M i) dτ := rfl

/-- With space-time cell volumes `M_i Δτ` as capacities, a time seam has rate `2/Δτ²`. -/
theorem flux_spacetime_time_rate {T : ℕ} (M : ι → 𝕜) (dτ : 𝕜) (i : ι) (a b : Fin (T + 1))
    (hM : M i ≠ 0) (hd : dτ ≠ 0) :
    seamRate (fun p : ι × Fin (T + 1) => M p.1 * dτ) (fluxCond (M i) dτ) (i, a) (i, b)
      = 2 / dτ ^ 2 := by
  unfold seamRate fluxCond
  field_simp
  ring

/-- With space-time cell volumes, a spatial seam (face measure `Δτ · face`) keeps its rate. -/
theorem flux_spacetime_spatial_rate {T : ℕ} (M : ι → 𝕜) (face dist dτ : 𝕜) (i j : ι)
    (a : Fin (T + 1)) (hd : dτ ≠ 0) (hi : M i ≠ 0) (hj : M j ≠ 0) :
    seamRate (fun p : ι × Fin (T + 1) => M p.1 * dτ) (fluxCond (dτ * face) dist) (i, a) (j, a)
      = seamRate M (fluxCond face dist) i j := by
  unfold seamRate fluxCond
  field_simp

end Rates

#print axioms flux_time_face
#print axioms flux_spacetime_time_rate
#print axioms flux_spacetime_spatial_rate

end OPH.SourceLaw
