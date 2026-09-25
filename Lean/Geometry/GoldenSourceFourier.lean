import Geometry.GoldenSourceAssignment
import Mathlib

/-!
# Finite Fourier factorization of a Cartesian golden orbit

The phasor uses the golden step before taking fractional parts: integer
Fourier modes have the same phase on the unit torus. The finite three-axis
sum factorizes and each one-axis sum obeys the geometric-series identity.
No fixed-mode asymptotic, infinite-volume structure factor, or random source
law is formalized by these finite algebraic statements.
-/
namespace OPH.GoldenSourceFourier
open scoped BigOperators

noncomputable def goldenPhasor (n : ℤ) : ℂ :=
  Complex.exp (2 * Real.pi * Complex.I * (n : ℂ) * (Real.goldenRatio : ℂ))

/-- Integer Fourier phases are unchanged by taking the orbit coordinate modulo one. -/
theorem integer_phase_fract (n : ℤ) (x : ℝ) :
    Complex.exp (2 * Real.pi * Complex.I * (n : ℂ) * ((Int.fract x : ℝ) : ℂ)) =
      Complex.exp (2 * Real.pi * Complex.I * (n : ℂ) * (x : ℂ)) := by
  have hx : ((Int.fract x : ℝ) : ℂ) = (x : ℂ) - ((⌊x⌋ : ℤ) : ℂ) := by
    simp only [Int.fract, Complex.ofReal_sub, Complex.ofReal_intCast]
  rw [hx]
  have he : 2 * (Real.pi : ℂ) * Complex.I * (n : ℂ) * ((x : ℂ) - ((⌊x⌋ : ℤ) : ℂ)) =
      2 * Real.pi * Complex.I * (n : ℂ) * (x : ℂ) +
        ((-(n * ⌊x⌋) : ℤ) : ℂ) * (2 * Real.pi * Complex.I) := by
    push_cast
    ring
  rw [he, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one]

/-- The phase at the actual fractional golden orbit point is its phasor power. -/
theorem orbit_phase (n : ℤ) (b : ℕ) :
    Complex.exp (2 * Real.pi * Complex.I * (n : ℂ) *
      (OPH.GoldenSourceAssignment.orbit b : ℂ)) = goldenPhasor n ^ b := by
  rw [OPH.GoldenSourceAssignment.orbit, integer_phase_fract]
  unfold goldenPhasor
  rw [← Complex.exp_nat_mul]
  congr 1
  push_cast
  ring

noncomputable def axisCoefficient (q : ℕ) (n : ℤ) : ℂ :=
  (∑ b ∈ Finset.range q, goldenPhasor n ^ b) / q

noncomputable def cubeCoefficient (q : ℕ) (nx ny nz : ℤ) : ℂ :=
  (∑ bx ∈ Finset.range q, (∑ iy ∈ Finset.range q, (∑ bz ∈ Finset.range q,
    goldenPhasor nx ^ bx * goldenPhasor ny ^ iy * goldenPhasor nz ^ bz))) / (q : ℂ)^3

/-- The actual finite Cartesian sum is a product of its one-axis means. -/
theorem cube_coefficient_factorization (q : ℕ) (nx ny nz : ℤ) :
    cubeCoefficient q nx ny nz =
      axisCoefficient q nx * axisCoefficient q ny * axisCoefficient q nz := by
  have hsum : (∑ bx ∈ Finset.range q, (∑ iy ∈ Finset.range q,
      (∑ bz ∈ Finset.range q, goldenPhasor nx ^ bx * goldenPhasor ny ^ iy * goldenPhasor nz ^ bz))) =
      (∑ bx ∈ Finset.range q, goldenPhasor nx ^ bx) *
      (∑ iy ∈ Finset.range q, goldenPhasor ny ^ iy) *
      (∑ bz ∈ Finset.range q, goldenPhasor nz ^ bz) := by
    simp only [← Finset.mul_sum, ← Finset.sum_mul]
  unfold cubeCoefficient axisCoefficient
  rw [hsum]
  ring

/-- Division-free geometric identity, also valid at the zero Fourier mode. -/
theorem axis_coefficient_geometric_identity (q : ℕ) (n : ℤ) :
    (q : ℂ) * axisCoefficient q n * (goldenPhasor n - 1) = goldenPhasor n ^ q - 1 := by
  by_cases hq : q = 0
  · subst q
    simp [axisCoefficient]
  · have hq' : (q : ℂ) ≠ 0 := by exact_mod_cast hq
    unfold axisCoefficient
    rw [mul_div_cancel₀ _ hq']
    exact geom_sum_mul (goldenPhasor n) q

/-- The geometric quotient requires the nonresonance denominator explicitly. -/
theorem axis_coefficient_geometric_quotient (q : ℕ) (n : ℤ)
    (hn : goldenPhasor n ≠ 1) :
    axisCoefficient q n = (goldenPhasor n ^ q - 1) /
      ((q : ℂ) * (goldenPhasor n - 1)) := by
  by_cases hq : q = 0
  · subst q
    simp [axisCoefficient]
  · have hq' : (q : ℂ) ≠ 0 := by exact_mod_cast hq
    apply (eq_div_iff (mul_ne_zero hq' (sub_ne_zero.mpr hn))).mpr
    have h := axis_coefficient_geometric_identity q n
    simpa only [mul_assoc, mul_left_comm, mul_comm] using h

/-- A zero mode contributes exactly one when the population is nonempty. -/
theorem axis_coefficient_zero (q : ℕ) (hq : 0 < q) : axisCoefficient q 0 = 1 := by
  have hq' : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  simp [axisCoefficient, goldenPhasor, hq']

#print axioms integer_phase_fract
#print axioms orbit_phase
#print axioms cube_coefficient_factorization
#print axioms axis_coefficient_geometric_identity
#print axioms axis_coefficient_geometric_quotient
#print axioms axis_coefficient_zero
end OPH.GoldenSourceFourier
