import Mathlib

/-!
# Pointwise temporal-filter comparison

This is the rational integrand inequality used by the common deterministic
filter comparison. Integration against a Fourier weight and the stochastic
covariance representation are separate analytical steps, not proved here.
-/
namespace OPH.Thermodynamics.TemporalFilterMonotonicity

theorem filter_integrand_monotone (a b omega : ℝ)
    (ha : 0 < a) (hab : a ≤ b) :
    a^2 / (a^2 + omega^2) ≤ b^2 / (b^2 + omega^2) := by
  have hb : 0 < b := lt_of_lt_of_le ha hab
  have hdena : 0 < a^2 + omega^2 := by positivity
  have hdenb : 0 < b^2 + omega^2 := by positivity
  apply (div_le_div_iff₀ hdena hdenb).mpr
  have hsq : a^2 ≤ b^2 := by nlinarith
  have hm := mul_nonneg (sub_nonneg.mpr hsq) (sq_nonneg omega)
  nlinarith

#print axioms filter_integrand_monotone
end OPH.Thermodynamics.TemporalFilterMonotonicity
