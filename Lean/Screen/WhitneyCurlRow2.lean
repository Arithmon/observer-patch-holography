import WhitneySourceAssembly

set_option autoImplicit false
set_option maxHeartbeats 300000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneySourceAssembly

open OPH.ConeCochainBridge
open OPH.WhitneySourceGeometry

theorem row2_binding (t : Fin 20) :
    ∃ (a b : Fin 12) (e : Fin 30),
      edgeIndex t 1 = coneOneIndexEquiv (Sum.inl a) ∧
      edgeIndex t 2 = coneOneIndexEquiv (Sum.inl b) ∧
      edgePositive t 1 = true ∧ edgePositive t 2 = true ∧
      edgeIndex t 5 = coneOneIndexEquiv (Sum.inr e) ∧
      faceIndex t 2 = coneTwoIndexEquiv (Sum.inl e) ∧
      facePositive t 2 = edgePositive t 5 ∧
      globalEdgeVertices (edgeIndex t 5) =
        orientedPair (edgePositive t 5) (a.succ, b.succ) := by
  fin_cases t <;> decide +kernel

theorem localCurlReal_mulVec_2 (z : Fin 6 → ℝ) :
    localCurlReal.mulVec z 2 = z 1 - z 2 + z 5 := by
  simp [localCurlReal, localCurl, Matrix.mulVec, dotProduct, Fin.sum_univ_succ,
    OPH.WhitneyFiniteCertificate.eval_mk, OPH.WhitneySourceGeometry.q]
  ring

theorem faceRestriction_coneCurl_2 (t : Fin 20) (x : ConeOne) :
    faceRestriction t (coneCurl x) 2 =
      localCurlReal.mulVec (edgeRestriction t x) 2 := by
  obtain ⟨a, b, e, h1, h2, hp1, hp2, h5, hf, hsign, hend⟩ := row2_binding t
  rw [faceRestriction_apply, localCurlReal_mulVec_2,
    edgeRestriction_apply, edgeRestriction_apply, edgeRestriction_apply]
  rw [hsign, hf, h1, h2, h5, hp1, hp2]
  have hend' := hend
  rw [h5] at hend'
  by_cases hs : edgePositive t 5 = true
  · simp [hs, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hs, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring
  · have hsfalse : edgePositive t 5 = false := Bool.eq_false_of_not_eq_true hs
    simp [hsfalse, globalEdgeVertices, orientedPair, coneOneIndexEquiv] at hend'
    rcases hend' with ⟨hl, hr⟩
    simp [coneOneCoordinate, coneTwoCoordinate, coneOneIndexEquiv,
      coneTwoIndexEquiv, signR, hsfalse, coneCurl,
      OPH.DiscreteCoulombGreen.realCoboundary_apply, hl, hr]
    ring

#print axioms faceRestriction_coneCurl_2

end OPH.WhitneySourceAssembly
