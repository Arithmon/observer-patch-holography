import WhitneySourceAssembly

set_option autoImplicit false
set_option maxHeartbeats 180000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneySourceAssembly

open OPH.ConeCochainBridge
open OPH.WhitneySourceGeometry

theorem localCurlReal_mulVec_3 (z : Fin 6 → ℝ) :
    localCurlReal.mulVec z 3 = z 3 - z 4 + z 5 := by
  simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
    OPH.WhitneyFiniteCertificate.eval_mk, OPH.WhitneySourceGeometry.q]
  ring

def signZ (b : Bool) : ℤ := if b then 1 else -1

theorem signZ_cast (b : Bool) : (signZ b : ℝ) = signR b := by
  cases b <;> norm_num [signZ, signR]

def boundaryCoefficient (t : Fin 20) (l : Fin 6) (e : Fin 30) : ℤ :=
  if edgeIndex t l = coneOneIndexEquiv (Sum.inr e) then signZ (edgePositive t l) else 0

def outerCoefficient (t : Fin 20) (e : Fin 30) : ℤ :=
  boundaryCoefficient t 3 e - boundaryCoefficient t 4 e + boundaryCoefficient t 5 e

theorem outerCoefficient_source : ∀ t e,
    OPH.LocalFaceMaxwellAction.faceIncidenceZ t e = outerCoefficient t e := by
  decide +kernel

theorem outer_face_binding : ∀ t,
    faceIndex t 3 = coneTwoIndexEquiv (Sum.inr t) ∧ facePositive t 3 = true := by
  decide +kernel

theorem boundary_slot (t : Fin 20) (l : Fin 6) (x : ConeOne)
    (h : ∃ e : Fin 30, edgeIndex t l = coneOneIndexEquiv (Sum.inr e)) :
    edgeRestriction t x l =
      ∑ e : Fin 30, (boundaryCoefficient t l e : ℝ) * x.2 e := by
  obtain ⟨e, he⟩ := h
  rw [edgeRestriction_apply, he]
  rw [Finset.sum_eq_single e]
  · simp [boundaryCoefficient, he, coneOneCoordinate, coneOneIndexEquiv, signZ_cast]
  · intro e' _ hne
    have hneq : edgeIndex t l ≠ coneOneIndexEquiv (Sum.inr e') := by
      rw [he]
      intro h
      exact hne (Sum.inr.inj (coneOneIndexEquiv.injective h.symm))
    simp [boundaryCoefficient, hneq]
  · simp

theorem boundary_slot_3 : ∀ t, ∃ e : Fin 30,
    edgeIndex t 3 = coneOneIndexEquiv (Sum.inr e) := by decide +kernel
theorem boundary_slot_4 : ∀ t, ∃ e : Fin 30,
    edgeIndex t 4 = coneOneIndexEquiv (Sum.inr e) := by decide +kernel
theorem boundary_slot_5 : ∀ t, ∃ e : Fin 30,
    edgeIndex t 5 = coneOneIndexEquiv (Sum.inr e) := by decide +kernel

theorem faceRestriction_coneCurl_3 (t : Fin 20) (x : ConeOne) :
    faceRestriction t (coneCurl x) 3 =
      localCurlReal.mulVec (edgeRestriction t x) 3 := by
  rw [faceRestriction_apply, localCurlReal_mulVec_3,
    boundary_slot t 3 x (boundary_slot_3 t),
    boundary_slot t 4 x (boundary_slot_4 t),
    boundary_slot t 5 x (boundary_slot_5 t)]
  obtain ⟨hf, hp⟩ := outer_face_binding t
  rw [hf, hp]
  rw [coneTwoCoordinate_inr]
  simp only [signR, if_true, one_mul]
  change OPH.LocalFaceMaxwellAction.faceCurvature x.2 t = _
  rw [OPH.LocalFaceMaxwellAction.faceCurvature_apply]
  simp only [OPH.LocalFaceMaxwellAction.faceIncidenceR]
  calc
    (∑ e : Fin 30,
        (OPH.LocalFaceMaxwellAction.faceIncidenceZ t e : ℝ) * x.2 e) =
      ∑ e : Fin 30, ((boundaryCoefficient t 3 e : ℝ) -
        boundaryCoefficient t 4 e + boundaryCoefficient t 5 e) * x.2 e := by
          apply Finset.sum_congr rfl
          intro e _
          rw [outerCoefficient_source t e]
          simp [outerCoefficient]
    _ = (∑ e : Fin 30, (boundaryCoefficient t 3 e : ℝ) * x.2 e) -
        (∑ e : Fin 30, (boundaryCoefficient t 4 e : ℝ) * x.2 e) +
        (∑ e : Fin 30, (boundaryCoefficient t 5 e : ℝ) * x.2 e) := by
          rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
          apply Finset.sum_congr rfl
          intro e _
          push_cast
          ring

#print axioms faceRestriction_coneCurl_3

end OPH.WhitneySourceAssembly
