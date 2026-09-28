import WhitneyCertifiedConsumers

set_option autoImplicit false
set_option maxHeartbeats 2000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneyOmittedCellCounterexample

open OPH.ConeCochainBridge
open OPH.WhitneySourceGeometry
open OPH.WhitneySourceAssembly
open OPH.WhitneyCertifiedConsumers
open OPH.WhitneyActualCertificates

noncomputable section

/-- A concrete actual edge field supported on the first radial edge. -/
def omittedWitness : ConeOne :=
  (fun p => if p = 0 then 1 else 0, 0)

theorem omittedWitness_ne_zero : omittedWitness ≠ 0 := by
  intro h
  have hp := congrArg (fun x : ConeOne => x.1 (0 : Fin 12)) h
  norm_num [omittedWitness] at hp

theorem first_cell_curl_witness :
    (localCurlReal.mulVec (edgeRestriction 0 omittedWitness)) 0 = 1 := by
  have he0 : edgeRestriction 0 omittedWitness 0 = 1 := by
    rw [edgeRestriction_apply]
    norm_num [edgePositive, signR]
    have hi : edgeIndex 0 0 = coneOneIndexEquiv (Sum.inl 0) := by decide +kernel
    rw [hi, coneOneCoordinate_inl]
    simp [omittedWitness]
  have he1 : edgeRestriction 0 omittedWitness 1 = 0 := by
    rw [edgeRestriction_apply]
    norm_num [edgePositive, signR]
    have hi : edgeIndex 0 1 = coneOneIndexEquiv (Sum.inl 1) := by decide +kernel
    rw [hi, coneOneCoordinate_inl]
    simp [omittedWitness]
  have he3 : edgeRestriction 0 omittedWitness 3 = 0 := by
    rw [edgeRestriction_apply]
    norm_num [edgePositive, signR]
    have hi : edgeIndex 0 3 = coneOneIndexEquiv (Sum.inr 0) := by decide +kernel
    rw [hi, coneOneCoordinate_inr]
    simp [omittedWitness]
  have hrow : ∀ e : Fin 6,
      localCurlReal 0 e = ![(1 : ℝ), -1, 0, 1, 0, 0] e := by
    intro e
    fin_cases e <;>
      norm_num [localCurlReal, localCurl, q, OPH.WhitneyFiniteCertificate.eval_mk]
  simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_succ]
  simp_rw [hrow]
  simp [he0, he1, he3]

theorem first_cell_curl_ne_zero :
    localCurlReal.mulVec (edgeRestriction 0 omittedWitness) ≠ 0 := by
  intro h
  have h0 := congrFun h 0
  rw [first_cell_curl_witness] at h0
  exact one_ne_zero h0

/-- The malformed assembly obtained by dropping the first actual tetrahedron. -/
def omittedCellStiffness : LinearMap.BilinForm ℝ ConeOne :=
  ∑ t ∈ (Finset.univ.erase (0 : Fin 20)),
    sourceVolume t • localStiffness.compl₁₂ (edgeRestriction t) (edgeRestriction t)

@[simp] theorem omittedCellStiffness_apply (x y : ConeOne) :
    omittedCellStiffness x y =
      ∑ t ∈ (Finset.univ.erase (0 : Fin 20)),
        sourceVolume t * localStiffness (edgeRestriction t x) (edgeRestriction t y) := by
  simp [omittedCellStiffness, LinearMap.compl₁₂_apply]

theorem complete_eq_omitted_add_first (x : ConeOne) :
    sourceStiffness x x = omittedCellStiffness x x +
      sourceVolume 0 * localStiffness (edgeRestriction 0 x) (edgeRestriction 0 x) := by
  rw [sourceStiffness, assemble_apply, omittedCellStiffness_apply]
  exact (Finset.sum_erase_add
    (s := (Finset.univ : Finset (Fin 20)))
    (f := fun t => sourceVolume t * localStiffness (edgeRestriction t x)
      (edgeRestriction t x)) (Finset.mem_univ (0 : Fin 20))).symm

theorem omitted_cell_changes_actual_stiffness :
    omittedCellStiffness omittedWitness omittedWitness <
      sourceStiffness omittedWitness omittedWitness := by
  have hlocal :
      0 < localStiffness (edgeRestriction 0 omittedWitness)
        (edgeRestriction 0 omittedWitness) := by
    rw [localStiffness_apply]
    exact (realForm_positive faceMass_posDef)
      (localCurlReal.mulVec (edgeRestriction 0 omittedWitness))
      first_cell_curl_ne_zero
  rw [complete_eq_omitted_add_first]
  have hv := volume_positive
  exact lt_add_of_pos_right _ (mul_pos hv hlocal)

#print axioms omitted_cell_changes_actual_stiffness

end

end OPH.WhitneyOmittedCellCounterexample
