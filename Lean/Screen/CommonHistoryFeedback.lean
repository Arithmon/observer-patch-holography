import Mathlib

/-! Exact classical record restoration and its finite accumulation factor.
The action stability, operator norm, and continuous comparison remain
analytical hypotheses; these lemmas do not construct physical controls. -/

namespace OPH.Screen.CommonHistoryFeedback

theorem restoration_error (x η : ℝ) : 2 * (x / 2 + η) - x = 2 * η := by
  ring

theorem restoration_error_abs (x η : ℝ) :
    |2 * (x / 2 + η) - x| = 2 * |η| := by
  rw [restoration_error, abs_mul]
  norm_num

theorem bounded_restoration (x η ε : ℝ) (h : |η| ≤ ε) :
    |2 * (x / 2 + η) - x| ≤ 2 * ε := by
  rw [restoration_error_abs]
  linarith

theorem triangular_accumulation (n : ℕ) (ε : ℝ) :
    2 * ε * (∑ k ∈ Finset.range n, ((k : ℝ) + 1)) =
      ε * (n : ℝ) * ((n : ℝ) + 1) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ]
    push_cast
    nlinarith

end OPH.Screen.CommonHistoryFeedback
