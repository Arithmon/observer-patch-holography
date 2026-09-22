import FiniteConditionalRepair

set_option autoImplicit false

open Filter Topology

namespace OPH.Thermodynamics

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- Entropy of the existing normalized Gibbs family, with the existing
`klTerm` convention at zero. -/
noncomputable def gibbsEntropy (E : Ω → ℝ) (β : ℝ) : ℝ :=
  shannon (gibbs E β)

/-- Ground-level degeneracy on the finite carrier. -/
noncomputable def groundDegeneracy (E : Ω → ℝ) (E0 : ℝ) : ℕ :=
  (Finset.univ.filter (fun x => E x = E0)).card

theorem gibbsEntropy_eq_sum (E : Ω → ℝ) (β : ℝ) :
    gibbsEntropy E β = ∑ x, Real.negMulLog (gibbs E β x) := by
  unfold gibbsEntropy shannon
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro x _
  by_cases h : gibbs E β x = 0
  · simp [h]
  · simp [klTerm, h, Real.negMulLog, neg_mul]

#print axioms gibbsEntropy_eq_sum

theorem gibbsEntropy_normalization [Nonempty Ω] (E : Ω → ℝ) (β : ℝ) :
    gibbsEntropy E β = Real.log (partitionZ E β) +
      β * ∑ x, gibbs E β x * E x := by
  have hsum : ∑ x, gibbs E β x = 1 := by
    simp only [gibbs, ← Finset.sum_div]
    exact div_self (partitionZ_pos E β).ne'
  rw [gibbsEntropy_eq_sum]
  have ht : ∀ x, Real.negMulLog (gibbs E β x) =
      gibbs E β x * Real.log (partitionZ E β) + β * (gibbs E β x * E x) := by
    intro x
    have hl : Real.log (gibbs E β x) = -β * E x - Real.log (partitionZ E β) := by
      rw [gibbs, gibbsWeight, Real.log_div (Real.exp_pos _).ne'
        (partitionZ_pos E β).ne', Real.log_exp]
    rw [Real.negMulLog, hl]
    ring
  simp_rw [ht]
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.mul_sum, hsum, one_mul]

#print axioms gibbsEntropy_normalization

/-- The log-sum bound consumes the same entropy convention. -/
theorem gibbsEntropy_le_log_card [Nonempty Ω] (E : Ω → ℝ) (β : ℝ) :
    gibbsEntropy E β ≤ Real.log (Fintype.card Ω) := by
  have hsum : ∑ x, gibbs E β x = 1 := by
    simp only [gibbs, ← Finset.sum_div]
    exact div_self (partitionZ_pos E β).ne'
  have h := logSum (gibbs E β) (fun _ => (1 : ℝ))
    (fun x => (div_pos (Real.exp_pos _) (partitionZ_pos E β)).le)
    (fun _ => zero_le_one) (by simp)
  rw [hsum] at h
  have h' : -Real.log (Fintype.card Ω : ℝ) ≤
      ∑ x, klTerm (gibbs E β x) 1 := by
    simpa [klTerm, Real.log_inv] using h
  simpa only [gibbsEntropy, shannon, neg_neg] using neg_le_neg h'

#print axioms gibbsEntropy_le_log_card

/-- Subtracting a common ground energy leaves normalized weights unchanged. -/
theorem gibbs_shift (E : Ω → ℝ) (E0 β : ℝ) (x : Ω) :
    gibbs E β x = gibbs (fun y => E y - E0) β x := by
  have hw : ∀ y, gibbsWeight (fun z => E z - E0) β y =
      gibbsWeight E β y * Real.exp (β * E0) := by
    intro y
    unfold gibbsWeight
    rw [show -β * (E y - E0) = -β * E y + β * E0 by ring, Real.exp_add]
  have hz : partitionZ (fun y => E y - E0) β =
      partitionZ E β * Real.exp (β * E0) := by
    simp only [partitionZ, hw, Finset.sum_mul]
  rw [gibbs, gibbs, hw, hz, mul_div_mul_right _ _ (Real.exp_pos _).ne']

#print axioms gibbs_shift

/-- Pointwise low-temperature convergence, under exactly the surface's
positive gap and inhabited ground-sector hypotheses. -/
theorem gibbs_tendsto_ground (E : Ω → ℝ) (E0 Δ : ℝ) (hΔ : 0 < Δ)
    (hground : ∀ x, E x = E0 ∨ E0 + Δ ≤ E x)
    (hne : (Finset.univ.filter (fun x => E x = E0)).Nonempty) (x : Ω) :
    Tendsto (fun β => gibbs E β x) atTop
      (𝓝 (if E x = E0 then (groundDegeneracy E E0 : ℝ)⁻¹ else 0)) := by
  have hg : (groundDegeneracy E E0 : ℝ) ≠ 0 := by
    exact_mod_cast (Finset.card_pos.mpr hne).ne'
  have hw : ∀ y, Tendsto (fun β => gibbsWeight (fun z => E z - E0) β y)
      atTop (𝓝 (if E y = E0 then (1 : ℝ) else 0)) := by
    intro y
    by_cases hy : E y = E0
    · simp [gibbsWeight, hy]
    · have hpos : 0 < E y - E0 := by
        rcases hground y with h | h
        · exact (hy h).elim
        · linarith
      rw [if_neg hy]
      simpa only [gibbsWeight, Function.comp_def, neg_mul] using
        Real.tendsto_exp_neg_atTop_nhds_zero.comp
          (Filter.Tendsto.atTop_mul_const hpos tendsto_id)
  have hz : Tendsto (fun β => partitionZ (fun y => E y - E0) β) atTop
      (𝓝 (groundDegeneracy E E0 : ℝ)) := by
    simpa [partitionZ, groundDegeneracy, Finset.sum_boole] using
      tendsto_finset_sum Finset.univ (fun y _ => hw y)
  have h := (hw x).div hz hg
  have heq : (fun β => gibbs E β x) =
      (fun β => gibbsWeight (fun y => E y - E0) β x /
        partitionZ (fun y => E y - E0) β) :=
    funext (fun β => gibbs_shift E E0 β x)
  rw [heq]
  simpa only [Pi.div_apply, ite_div, one_div, zero_div] using h

#print axioms gibbs_tendsto_ground

/-- Residual Gibbs entropy is precisely the logarithm of ground degeneracy. -/
theorem gibbsEntropy_tendsto (E : Ω → ℝ) (E0 Δ : ℝ) (hΔ : 0 < Δ)
    (hground : ∀ x, E x = E0 ∨ E0 + Δ ≤ E x)
    (hne : (Finset.univ.filter (fun x => E x = E0)).Nonempty) :
    Tendsto (gibbsEntropy E) atTop (𝓝 (Real.log (groundDegeneracy E E0))) := by
  have hg : (groundDegeneracy E E0 : ℝ) ≠ 0 := by
    exact_mod_cast (Finset.card_pos.mpr hne).ne'
  have h := tendsto_finset_sum Finset.univ (fun x _ =>
    Real.continuous_negMulLog.continuousAt.tendsto.comp
      (gibbs_tendsto_ground E E0 Δ hΔ hground hne x))
  have hs : (∑ x : Ω, Real.negMulLog
      (if E x = E0 then (groundDegeneracy E E0 : ℝ)⁻¹ else 0)) =
      Real.log (groundDegeneracy E E0) := by
    simp only [apply_ite, Real.negMulLog_zero, ← Finset.sum_filter]
    rw [Finset.sum_const, nsmul_eq_mul]
    change (groundDegeneracy E E0 : ℝ) *
      Real.negMulLog (groundDegeneracy E E0 : ℝ)⁻¹ = _
    rw [Real.negMulLog, Real.log_inv]
    field_simp
  simpa only [Function.comp_def, ← gibbsEntropy_eq_sum, hs] using h

#print axioms gibbsEntropy_tendsto

/-- The zero-entropy third law holds exactly for a unique ground state. -/
theorem gibbsEntropy_tendsto_zero_iff (E : Ω → ℝ) (E0 Δ : ℝ) (hΔ : 0 < Δ)
    (hground : ∀ x, E x = E0 ∨ E0 + Δ ≤ E x)
    (hne : (Finset.univ.filter (fun x => E x = E0)).Nonempty) :
    Tendsto (gibbsEntropy E) atTop (𝓝 0) ↔ groundDegeneracy E E0 = 1 := by
  have hlim := gibbsEntropy_tendsto E E0 Δ hΔ hground hne
  have hg : (0 : ℝ) < groundDegeneracy E E0 := by
    exact_mod_cast Finset.card_pos.mpr hne
  constructor
  · intro h
    have heq : Real.log (groundDegeneracy E E0) = 0 := tendsto_nhds_unique hlim h
    have hc : (groundDegeneracy E E0 : ℝ) = 1 := by
      rcases Real.log_eq_zero.mp heq with h | h | h <;> linarith
    exact_mod_cast hc
  · intro h
    simpa [h] using hlim

#print axioms gibbsEntropy_tendsto_zero_iff

/-- Two ground states and one excited state: genuinely two energy levels. -/
def doubleGroundEnergy (x : Fin 3) : ℝ := if x = 2 then 1 else 0

theorem doubleGround_degeneracy : groundDegeneracy doubleGroundEnergy 0 = 2 := by
  have hu : (Finset.univ : Finset (Fin 3)) = {0, 1, 2} := by decide
  norm_num [groundDegeneracy, hu, doubleGroundEnergy, Fin.ext_iff, Finset.filter_insert, Finset.filter_singleton]

#print axioms doubleGround_degeneracy

theorem doubleGround_entropy_limit :
    Tendsto (gibbsEntropy doubleGroundEnergy) atTop (𝓝 (Real.log 2)) ∧
      Real.log (2 : ℝ) ≠ 0 ∧
      ¬ Tendsto (gibbsEntropy doubleGroundEnergy) atTop (𝓝 0) := by
  have hground : ∀ x, doubleGroundEnergy x = 0 ∨ 0 + 1 ≤ doubleGroundEnergy x := by
    intro x
    by_cases hx : x = 2 <;> simp [doubleGroundEnergy, hx]
  have hne : (Finset.univ.filter (fun x => doubleGroundEnergy x = 0)).Nonempty := by
    exact ⟨0, by simp [doubleGroundEnergy, Fin.ext_iff]⟩
  have hlim := gibbsEntropy_tendsto doubleGroundEnergy 0 1 (by norm_num) hground hne
  have hiff := gibbsEntropy_tendsto_zero_iff doubleGroundEnergy 0 1 (by norm_num) hground hne
  rw [doubleGround_degeneracy] at hlim hiff
  refine ⟨by simpa using hlim, (Real.log_pos (by norm_num : (1 : ℝ) < 2)).ne', ?_⟩
  simpa using hiff

#print axioms doubleGround_entropy_limit

end OPH.Thermodynamics
