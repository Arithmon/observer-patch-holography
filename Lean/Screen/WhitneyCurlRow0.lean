import WhitneySourceAssembly

set_option autoImplicit false
set_option maxHeartbeats 300000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneySourceAssembly

open OPH.ConeCochainBridge
open OPH.WhitneySourceGeometry

theorem row0_binding (t : Fin 20) :
    ∃ (a b : Fin 12) (e : Fin 30),
      edgeIndex t 0 = coneOneIndexEquiv (Sum.inl a) ∧
      edgeIndex t 1 = coneOneIndexEquiv (Sum.inl b) ∧
      edgePositive t 0 = true ∧ edgePositive t 1 = true ∧
      edgeIndex t 3 = coneOneIndexEquiv (Sum.inr e) ∧
      faceIndex t 0 = coneTwoIndexEquiv (Sum.inl e) ∧
      facePositive t 0 = edgePositive t 3 ∧
      globalEdgeVertices (edgeIndex t 3) =
        orientedPair (edgePositive t 3) (a.succ, b.succ) := by
  fin_cases t <;> decide +kernel

theorem localCurlReal_mulVec_0 (z : Fin 6 → ℝ) :
    localCurlReal.mulVec z 0 = z 0 - z 1 + z 3 := by
  simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
    OPH.WhitneyFiniteCertificate.eval_mk, OPH.WhitneySourceGeometry.q]
  ring

theorem faceRestriction_coneCurl_0 (t : Fin 20) (x : ConeOne) :
    faceRestriction t (coneCurl x) 0 =
      localCurlReal.mulVec (edgeRestriction t x) 0 := by
  obtain ⟨a, b, e, h0, h1, hp0, hp1, h3, hf, hsign, hend⟩ := row0_binding t
  rw [faceRestriction_apply, localCurlReal_mulVec_0,
    edgeRestriction_apply, edgeRestriction_apply, edgeRestriction_apply]
  rw [hsign, hf]
  rw [h0, h1, h3, hp0, hp1]
  have hend' := hend
  rw [h3] at hend'
  by_cases hs : edgePositive t 3 = true
  · simp [hs, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
      edgeRestriction, coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hs, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring
  · have hsfalse : edgePositive t 3 = false := by
      exact Bool.eq_false_of_not_eq_true hs
    simp [hsfalse, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
      edgeRestriction, coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hsfalse, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring

#print axioms faceRestriction_coneCurl_0

end OPH.WhitneySourceAssembly
