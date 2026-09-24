import Mathlib.Tactic

set_option autoImplicit false

namespace OPH.M1VacuumFidelity

/-!
Finite algebra reductions for the analytic full-spectrum vacuum theorem.
The Gaussian state interpretation, spectral counting, time averages and
all-level asymptotic conclusions are proved in the accompanying derivation.
-/

/-- A symmetric symplectic mode has the stated nonnegative one-step defect. -/
theorem one_step_excitation (a b c : ℝ) (h : a ^ 2 - b * c = 1) :
    (2 * a ^ 2 + b ^ 2 + c ^ 2 - 2) / 4 = (b + c) ^ 2 / 4 := by
  nlinarith

/-- The second-order canonical matrix has determinant one. -/
theorem second_order_symplectic (z : ℝ) :
    (1 - z ^ 2 / 2) ^ 2 - z * (-z + z ^ 3 / 4) = 1 := by
  ring

/-- Exact produced particle number after one second-order step. -/
theorem second_order_excitation (z : ℝ) :
    (2 * (1 - z ^ 2 / 2) ^ 2 + z ^ 2 + (-z + z ^ 3 / 4) ^ 2 - 2) / 4 =
      z ^ 6 / 64 := by
  ring

/-- The invariant positive form differs from the original vacuum form. -/
theorem invariant_form (a b c x y : ℝ) (h : a ^ 2 - b * c = 1) :
    -c * (a * x + b * y) ^ 2 + b * (c * x + a * y) ^ 2 =
      -c * x ^ 2 + b * y ^ 2 := by
  calc
    _ = (a ^ 2 - b * c) * (-c * x ^ 2 + b * y ^ 2) := by ring
    _ = _ := by rw [h]; ring

/-- The signed composition retains the intended first-order duration. -/
theorem composition_duration (r w : ℝ) (h : (2 - r) * w = 1) :
    2 * w + (-r * w) = 1 := by
  nlinarith

/-- The cubic defect cancels for the declared real cubic root. -/
theorem cubic_cancellation (r w : ℝ) (h : r ^ 3 = 2) :
    2 * w ^ 3 + (-r * w) ^ 3 = 0 := by
  calc
    _ = (2 - r ^ 3) * w ^ 3 := by ring
    _ = 0 := by rw [h]; ring

/-- A bounded spectral maximum and its mean force a bulk fraction. -/
theorem bulk_fraction (N K M : ℝ) (hM : 0 < M)
    (h : N * M ≤ K * (2 * M) + (N - K) * (M / 2)) : N / 3 ≤ K := by
  nlinarith

/-- Fourth-order energy error on the explicit sparse schedule vanishes. -/
theorem repair_energy_scale (x : ℝ) (hx : x ≠ 0) :
    x ^ 96 * (1 / x ^ 20) ^ 8 * (x ^ 7) ^ 9 = 1 / x := by
  field_simp

/-- The global vacuum number error has the stronger stated decay. -/
theorem repair_number_scale (x : ℝ) (hx : x ≠ 0) :
    x ^ 96 * (1 / x ^ 20) ^ 8 * (x ^ 7) ^ 8 = 1 / x ^ 8 := by
  field_simp

/-- The second-order program fails on the same explicit schedule. -/
theorem second_order_energy_scale (x : ℝ) (hx : x ≠ 0) :
    x ^ 96 * (1 / x ^ 20) ^ 4 * (x ^ 7) ^ 5 = x ^ 51 := by
  field_simp

end OPH.M1VacuumFidelity
