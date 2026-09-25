import Geometry.SourceNetLayeredOrder

/-!
# Uniform read counts and finite conformal windows

The first results use the actual layered read-chain theorem: equal endpoint
layers imply equal elapsed uniform-tick clocks, independently of spatial
motion. Such a reading cannot also equal its strictly smaller Lorentz factor.

The finite-window control evaluates the exact polynomial obtained by
integrating the supplied conformal profile `σ(t)=1+t` over a vertical
four-dimensional diamond, with the common spatial-volume constant omitted.
The polynomial evaluation and clock mismatch are kernel-checked here; the
geometric integral identifying that polynomial is established analytically
in the paper. No FLRW dynamics or physical clock is inferred.
-/
namespace OPH.SourceNetClockSeparation
open OPH.SourceNetLayeredOrder

variable {S : Type*} {step : S → S → Prop}

/-- Every one-layer read contributes the same declared tick. -/
def uniformReadClock (tick : ℝ) (n : ℕ) : ℝ := tick * n

/-- The clock depends only on the layer separation of the endpoints. -/
theorem uniform_read_clock_endpoints (tick : ℝ) (c : ℕ → ℕ × S) (n : ℕ)
    (hc : IsReadChain step c n) :
    uniformReadClock tick n = tick * ((c n).1 - (c 0).1 : ℕ) := by
  rw [readChain_layer c n hc]
  simp [uniformReadClock]

/-- Spatially distinct read chains between the same layers read the same time. -/
theorem uniform_read_clock_motion_independent (tick : ℝ)
    (c d : ℕ → ℕ × S) (n m : ℕ)
    (hc : IsReadChain step c n) (hd : IsReadChain step d m)
    (hstart : (c 0).1 = (d 0).1) (hend : (c n).1 = (d m).1) :
    uniformReadClock tick n = uniformReadClock tick m := by
  rw [uniform_read_clock_endpoints tick c n hc,
    uniform_read_clock_endpoints tick d m hd, hstart, hend]

/-- Positive uniform layer time is not its nonzero-boost proper-time value. -/
theorem uniform_read_clock_not_lorentz (tick beta : ℝ) (n : ℕ)
    (htick : 0 < tick) (hn : 0 < n) (hbeta : 0 < beta) (hsub : beta < 1) :
    uniformReadClock tick n ≠
      uniformReadClock tick n * Real.sqrt (1 - beta ^ 2) := by
  have htime : 0 < uniformReadClock tick n := by
    unfold uniformReadClock
    exact mul_pos htick (by exact_mod_cast hn)
  have harg : 0 ≤ 1 - beta ^ 2 := by nlinarith
  have hs := Real.sq_sqrt harg
  have hsmall : Real.sqrt (1 - beta ^ 2) < 1 := by
    have hnonneg := Real.sqrt_nonneg (1 - beta ^ 2)
    nlinarith
  exact ne_of_gt (mul_lt_of_lt_one_right htime hsmall)

/-- Rational part of the continuum diamond mass for `σ(t)=1+t`. -/
def radiationDiamondMass (T : ℚ) : ℚ :=
  T ^ 4 * (99 * T ^ 4 + 672 * T ^ 3 + 1792 * T ^ 2 + 2240 * T + 1120) / 35840

/-- The vertical proper duration for this supplied conformal profile. -/
def radiationProperDuration (T : ℚ) : ℚ := T + T ^ 2 / 2

theorem radiation_window_ratio :
    radiationDiamondMass 1 / radiationDiamondMass (1 / 2) = 1516288 / 44451 ∧
    radiationProperDuration 1 / radiationProperDuration (1 / 2) = 12 / 5 := by
  norm_num [radiationDiamondMass, radiationProperDuration]

/-- Taking the fourth root of the mass ratio cannot equal the proper-time ratio. -/
theorem radiation_window_clock_mismatch :
    radiationDiamondMass 1 / radiationDiamondMass (1 / 2) ≠
      (radiationProperDuration 1 / radiationProperDuration (1 / 2)) ^ 4 := by
  norm_num [radiationDiamondMass, radiationProperDuration]

#print axioms uniform_read_clock_endpoints
#print axioms uniform_read_clock_motion_independent
#print axioms uniform_read_clock_not_lorentz
#print axioms radiation_window_ratio
#print axioms radiation_window_clock_mismatch
end OPH.SourceNetClockSeparation
