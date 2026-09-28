import WhitneySourceAssembly

set_option autoImplicit false
set_option maxHeartbeats 300000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneySourceAssembly

open OPH.ConeCochainBridge
open OPH.WhitneySourceGeometry

theorem row1_binding (t : Fin 20) :
    ∃ (a b : Fin 12) (e : Fin 30),
      edgeIndex t 0 = coneOneIndexEquiv (Sum.inl a) ∧
      edgeIndex t 2 = coneOneIndexEquiv (Sum.inl b) ∧
      edgePositive t 0 = true ∧ edgePositive t 2 = true ∧
      edgeIndex t 4 = coneOneIndexEquiv (Sum.inr e) ∧
      faceIndex t 1 = coneTwoIndexEquiv (Sum.inl e) ∧
      facePositive t 1 = edgePositive t 4 ∧
      globalEdgeVertices (edgeIndex t 4) =
        orientedPair (edgePositive t 4) (a.succ, b.succ) := by
  fin_cases t <;> decide +kernel

theorem localCurlReal_mulVec_1 (z : Fin 6 → ℝ) :
    localCurlReal.mulVec z 1 = z 0 - z 2 + z 4 := by
  simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
    OPH.WhitneyFiniteCertificate.eval_mk, OPH.WhitneySourceGeometry.q]
  ring

theorem faceRestriction_coneCurl_1 (t : Fin 20) (x : ConeOne) :
    faceRestriction t (coneCurl x) 1 =
      localCurlReal.mulVec (edgeRestriction t x) 1 := by
  obtain ⟨a, b, e, h0, h2, hp0, hp2, h4, hf, hsign, hend⟩ := row1_binding t
  rw [faceRestriction_apply, localCurlReal_mulVec_1,
    edgeRestriction_apply, edgeRestriction_apply, edgeRestriction_apply]
  rw [hsign, hf, h0, h2, h4, hp0, hp2]
  have hend' := hend
  rw [h4] at hend'
  by_cases hs : edgePositive t 4 = true
  · simp [hs, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hs, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring
  · have hsfalse : edgePositive t 4 = false := Bool.eq_false_of_not_eq_true hs
    simp [hsfalse, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hsfalse, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring

#print axioms faceRestriction_coneCurl_1

end OPH.WhitneySourceAssembly
