import ObserverPatchHolography.EinsteinBranch.EdgeCenterTiltCocycle

/-!
# The survival cocycle does not select a tilt

Every real exponent realizes the existing positive continuous cocycle.
Nonnegative exponents also obey survival at positive times. These are
counterexamples to parameter selection by the cocycle law alone; no
measured spectral index or cosmological source is identified.
-/
namespace OPH.EinsteinBranch.EdgeCenterTiltFreedom

noncomputable def exponentialCocycle (theta : ℝ) :
    OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle where
  u := fun t => Real.exp (-(theta * t))
  pos := fun _ => Real.exp_pos _
  mul := by
    intro s t
    rw [mul_add, neg_add, Real.exp_add]
  cont := by fun_prop

theorem arbitrary_tilt (theta : ℝ) :
    (exponentialCocycle theta).tilt = theta := by
  apply OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle.tilt_eq_of_exp
  intro s
  rfl

theorem cocycle_does_not_select_tilt (theta : ℝ) :
    ∃ S : OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle,
      S.tilt ≠ theta := by
  refine ⟨exponentialCocycle (theta + 1), ?_⟩
  rw [arbitrary_tilt]
  linarith

theorem nonnegative_tilt_survival (theta t : ℝ) (htheta : 0 ≤ theta) (ht : 0 ≤ t) :
    0 < (exponentialCocycle theta).u t ∧ (exponentialCocycle theta).u t ≤ 1 := by
  constructor
  · exact Real.exp_pos _
  · change Real.exp (-(theta * t)) ≤ 1
    apply Real.exp_le_one_iff.mpr
    exact neg_nonpos.mpr (mul_nonneg htheta ht)

#print axioms arbitrary_tilt
#print axioms cocycle_does_not_select_tilt
#print axioms nonnegative_tilt_survival
end OPH.EinsteinBranch.EdgeCenterTiltFreedom
