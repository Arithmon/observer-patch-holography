import WhitneyActualCertificates
import WhitneyConeMass
import SeamCurrentEdge30Moment

set_option autoImplicit false
set_option maxHeartbeats 300000
set_option maxRecDepth 8192

open scoped BigOperators Matrix goldenRatio

namespace OPH.WhitneySourceGeometry

open OPH.WhitneyFiniteCertificate
open OPH.WhitneyActualCertificates
open OPH.LocalFaceMaxwellAction
open OPH.SeamCurrentCarrierQuotient

/-! The finite tables below are bound back to the current source definitions
by the `*_source_bound` theorems.  They are indexing devices, not independent
geometry authorities. -/

def tetraVertex : Fin 20 → Fin 4 → Fin 13 := ![
  ![0,1,2,3], ![0,1,4,2], ![0,1,3,5], ![0,1,7,4], ![0,1,5,7],
  ![0,2,6,3], ![0,2,4,8], ![0,2,8,6], ![0,3,9,5], ![0,3,6,9],
  ![0,4,7,10], ![0,4,10,8], ![0,5,11,7], ![0,5,9,11], ![0,6,8,12],
  ![0,6,12,9], ![0,7,11,10], ![0,8,10,12], ![0,9,12,11], ![0,10,11,12]
]

def edgeIndex : Fin 20 → Fin 6 → Fin 42 := ![
  ![0,1,2,12,13,17], ![0,3,1,14,12,18], ![0,2,4,13,15,21],
  ![0,6,3,16,14,24], ![0,4,6,15,16,27], ![1,5,2,19,17,22],
  ![1,3,7,18,20,25], ![1,7,5,20,19,30], ![2,8,4,23,21,28],
  ![2,5,8,22,23,31], ![3,6,9,24,26,33], ![3,9,7,26,25,35],
  ![4,10,6,29,27,34], ![4,8,10,28,29,37], ![5,7,11,30,32,36],
  ![5,11,8,32,31,38], ![6,10,9,34,33,39], ![7,9,11,35,36,40],
  ![8,11,10,38,37,41], ![9,10,11,39,40,41]
]

def edgePositive : Fin 20 → Fin 6 → Bool := ![
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,false], ![true,true,true,true,true,true],
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,false], ![true,true,true,true,true,true],
  ![true,true,true,true,true,true], ![true,true,true,true,true,false],
  ![true,true,true,true,true,false], ![true,true,true,true,true,true],
  ![true,true,true,true,true,false], ![true,true,true,true,true,true]
]

def faceIndex : Fin 20 → Fin 4 → Fin 50 := ![
  ![0,1,5,30], ![2,0,6,31], ![1,3,9,32], ![4,2,12,33],
  ![3,4,15,34], ![7,5,10,35], ![6,8,13,36], ![8,7,18,37],
  ![11,9,16,38], ![10,11,19,39], ![12,14,21,40], ![14,13,23,41],
  ![17,15,22,42], ![16,17,25,43], ![18,20,24,44], ![20,19,26,45],
  ![22,21,27,46], ![23,24,28,47], ![26,25,29,48], ![27,28,29,49]
]

def facePositive : Fin 20 → Fin 4 → Bool := ![
  ![true,true,true,true], ![true,true,false,true], ![true,true,true,true],
  ![true,true,false,true], ![true,true,true,true], ![true,true,false,true],
  ![true,true,true,true], ![true,true,false,true], ![true,true,false,true],
  ![true,true,true,true], ![true,true,true,true], ![true,true,false,true],
  ![true,true,false,true], ![true,true,true,true], ![true,true,true,true],
  ![true,true,false,true], ![true,true,false,true], ![true,true,true,true],
  ![true,true,false,true], ![true,true,true,true]
]

def localEdge : Fin 6 → Fin 4 × Fin 4 :=
  ![(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def localFace : Fin 4 → Fin 4 × Fin 4 × Fin 4 :=
  ![(0,1,2),(0,1,3),(0,2,3),(1,2,3)]

def globalEdgeVertices (e : Fin 42) : Fin 13 × Fin 13 :=
  match finSumFinEquiv.symm e with
  | Sum.inl p => (0, p.succ)
  | Sum.inr s => ((seamLeft s).succ, (seamRight s).succ)

def localEdgeVertices (t : Fin 20) (e : Fin 6) : Fin 13 × Fin 13 :=
  (tetraVertex t (localEdge e).1, tetraVertex t (localEdge e).2)

def orientedPair (positive : Bool) (p : Fin 13 × Fin 13) : Fin 13 × Fin 13 :=
  if positive then p else (p.2, p.1)

theorem edge_table_source_bound :
    ∀ t e, globalEdgeVertices (edgeIndex t e) =
      orientedPair (edgePositive t e) (localEdgeVertices t e) := by
  decide +kernel

theorem tetra_table_source_bound : ∀ t,
    tetraVertex t 0 = 0 ∧
    tetraVertex t 1 = (faceVertices t).1.succ ∧
    tetraVertex t 2 = (faceVertices t).2.1.succ ∧
    tetraVertex t 3 = (faceVertices t).2.2.succ := by
  decide +kernel

def globalFaceVertices (f : Fin 50) : Fin 13 × Fin 13 × Fin 13 :=
  match finSumFinEquiv.symm f with
  | Sum.inl e => (0, (seamLeft e).succ, (seamRight e).succ)
  | Sum.inr t =>
      ((faceVertices t).1.succ, (faceVertices t).2.1.succ,
        (faceVertices t).2.2.succ)

def localFaceVertices (t : Fin 20) (f : Fin 4) : Fin 13 × Fin 13 × Fin 13 :=
  (tetraVertex t (localFace f).1,
    tetraVertex t (localFace f).2.1,
    tetraVertex t (localFace f).2.2)

def orientedTriple (positive : Bool) (p : Fin 13 × Fin 13 × Fin 13) :
    Fin 13 × Fin 13 × Fin 13 :=
  if positive then p else (p.1, p.2.2, p.2.1)

theorem face_table_source_bound :
    ∀ t f, globalFaceVertices (faceIndex t f) =
      orientedTriple (facePositive t f) (localFaceVertices t f) := by
  decide +kernel

theorem edgeIndex_surjective : Function.Surjective (fun p : Fin 20 × Fin 6 =>
    edgeIndex p.1 p.2) := by
  decide +kernel

theorem faceIndex_surjective : Function.Surjective (fun p : Fin 20 × Fin 4 =>
    faceIndex p.1 p.2) := by
  decide +kernel

def q (a b : ℚ) : Q5 := QuadraticAlgebra.mk a b
def phiQ5 : Q5 := q (1/2) (1/2)

def vertexQ5 : Fin 13 → Fin 3 → Q5 := ![
  ![q 0 0,q 0 0,q 0 0],
  ![q 0 0,q (-1) 0,-phiQ5], ![q (-1) 0,-phiQ5,q 0 0],
  ![-phiQ5,q 0 0,q (-1) 0], ![q 1 0,-phiQ5,q 0 0],
  ![q 0 0,q 1 0,-phiQ5], ![-phiQ5,q 0 0,q 1 0],
  ![phiQ5,q 0 0,q (-1) 0], ![q 0 0,q (-1) 0,phiQ5],
  ![q (-1) 0,phiQ5,q 0 0], ![phiQ5,q 0 0,q 1 0],
  ![q 1 0,phiQ5,q 0 0], ![q 0 0,q 1 0,phiQ5]
]

theorem vertexQ5_source_bound (p : Fin 12) (d : Fin 3) :
    eval (vertexQ5 p.succ d) = OPH.SeamCurrentEdge30Moment.portVector p d := by
  fin_cases p <;> fin_cases d <;>
    simp [vertexQ5, q, phiQ5, OPH.SeamCurrentEdge30Moment.portVector,
      Real.goldenRatio]
  all_goals ring

def ray (t : Fin 20) (i : Fin 3) (d : Fin 3) : Q5 :=
  vertexQ5 (tetraVertex t i.succ) d - vertexQ5 (tetraVertex t 0) d

def dotQ5 (x y : Fin 3 → Q5) : Q5 := ∑ d, x d * y d

def rayGram (t : Fin 20) : Matrix (Fin 3) (Fin 3) Q5 :=
  fun i j => dotQ5 (ray t i) (ray t j)

def expectedRayGram : Matrix (Fin 3) (Fin 3) Q5 :=
  fun i j => (if i = j then q 2 0 else 0) + phiQ5

theorem source_tetrahedra_congruent : ∀ t, rayGram t = expectedRayGram := by
  decide +kernel

def signedDeterminant (t : Fin 20) : Q5 :=
  let r0 := ray t 0
  let r1 := ray t 1
  let r2 := ray t 2
  r0 0 * (r1 1 * r2 2 - r1 2 * r2 1) -
    r1 0 * (r0 1 * r2 2 - r0 2 * r2 1) +
    r2 0 * (r0 1 * r1 2 - r0 2 * r1 1)

theorem source_tetrahedra_signedDeterminant :
    ∀ t, signedDeterminant t = q 3 1 := by
  decide +kernel

def volumeQ5 : Q5 := q 3 1 * q (1/6) 0

theorem volumeQ5_eq_determinant_div_six (t : Fin 20) :
    volumeQ5 = signedDeterminant t / (6 : Q5) := by
  rw [source_tetrahedra_signedDeterminant]
  decide +kernel

theorem volume_positive : 0 < eval volumeQ5 := by
  apply (positive_iff_eval_pos volumeQ5).mp
  decide +kernel

/-- Gram matrix of the four barycentric gradients, normalized before the
positive tetrahedron volume is restored. -/
def gradientGram : Matrix (Fin 4) (Fin 4) Q5 := !![
  q (21/2) (-9/2), q (-7/2) (3/2), q (-7/2) (3/2), q (-7/2) (3/2);
  q (-7/2) (3/2), q (3/2) (-1/2), q 1 (-1/2), q 1 (-1/2);
  q (-7/2) (3/2), q 1 (-1/2), q (3/2) (-1/2), q 1 (-1/2);
  q (-7/2) (3/2), q 1 (-1/2), q 1 (-1/2), q (3/2) (-1/2)
]

theorem gradientGram_rows_sum_zero : ∀ i, (∑ j, gradientGram i j) = 0 := by
  decide +kernel

theorem gradientGram_inverts_source_rays : ∀ t i j,
    (∑ k : Fin 3, rayGram t i k * gradientGram k.succ j.succ) =
      if i = j then 1 else 0 := by
  decide +kernel

def moment (i j : Fin 4) : Q5 := if i = j then q (1/10) 0 else q (1/20) 0

/-- The four positive degree-two tetrahedron quadrature points used by the
independent source verifier, represented exactly in `Q5`. -/
def quadratureA : Q5 := q (1/4) (3/20)
def quadratureB : Q5 := q (1/4) (-1/20)

def quadratureBarycentric (r i : Fin 4) : Q5 :=
  if r = i then quadratureA else quadratureB

/-- Normalized quadrature moment.  Multiplication by `volumeQ5` restores the
physical tetrahedron weight `volume/4` at each point. -/
def quadratureMoment (i j : Fin 4) : Q5 :=
  q (1/4) 0 * ∑ r : Fin 4, quadratureBarycentric r i * quadratureBarycentric r j

theorem four_point_quadrature_moment_exact :
    ∀ i j, quadratureMoment i j = moment i j := by
  decide +kernel

theorem four_point_quadrature_weighted_exact : ∀ i j,
    q (1/4) 0 * volumeQ5 *
        ∑ r : Fin 4, quadratureBarycentric r i * quadratureBarycentric r j =
      volumeQ5 * moment i j := by
  intro i j
  rw [← four_point_quadrature_moment_exact i j]
  simp only [quadratureMoment]
  ring

def derivedEdgeMass : Matrix (Fin 6) (Fin 6) Q5 := fun e f =>
  let i := (localEdge e).1
  let j := (localEdge e).2
  let k := (localEdge f).1
  let l := (localEdge f).2
  moment i k * gradientGram j l - moment i l * gradientGram j k -
    moment j k * gradientGram i l + moment j l * gradientGram i k

theorem edgeMass_from_source_geometry : derivedEdgeMass = edgeMass := by
  decide +kernel

def quadratureEdgeMass : Matrix (Fin 6) (Fin 6) Q5 := fun e f =>
  let i := (localEdge e).1
  let j := (localEdge e).2
  let k := (localEdge f).1
  let l := (localEdge f).2
  quadratureMoment i k * gradientGram j l -
    quadratureMoment i l * gradientGram j k -
    quadratureMoment j k * gradientGram i l +
    quadratureMoment j l * gradientGram i k

theorem edgeMass_from_four_point_quadrature : quadratureEdgeMass = edgeMass := by
  rw [← edgeMass_from_source_geometry]
  funext e f
  simp only [quadratureEdgeMass, derivedEdgeMass]
  rw [four_point_quadrature_moment_exact, four_point_quadrature_moment_exact,
    four_point_quadrature_moment_exact, four_point_quadrature_moment_exact]

def cycleFace (f : Fin 4) : Fin 3 → Fin 4 × Fin 4 × Fin 4 := fun r =>
  let p := localFace f
  match r.1 with
  | 0 => p
  | 1 => (p.2.1, p.2.2, p.1)
  | _ => (p.2.2, p.1, p.2.1)

def crossDot (i j k l : Fin 4) : Q5 :=
  gradientGram i k * gradientGram j l - gradientGram i l * gradientGram j k

def derivedFaceMass : Matrix (Fin 4) (Fin 4) Q5 := fun f h =>
  4 * ∑ r : Fin 3, ∑ s : Fin 3,
    moment (cycleFace f r).1 (cycleFace h s).1 *
      crossDot (cycleFace f r).2.1 (cycleFace f r).2.2
        (cycleFace h s).2.1 (cycleFace h s).2.2

theorem faceMass_from_source_geometry : derivedFaceMass = faceMass := by
  decide +kernel

def quadratureFaceMass : Matrix (Fin 4) (Fin 4) Q5 := fun f h =>
  4 * ∑ r : Fin 3, ∑ s : Fin 3,
    quadratureMoment (cycleFace f r).1 (cycleFace h s).1 *
      crossDot (cycleFace f r).2.1 (cycleFace f r).2.2
        (cycleFace h s).2.1 (cycleFace h s).2.2

theorem faceMass_from_four_point_quadrature : quadratureFaceMass = faceMass := by
  rw [← faceMass_from_source_geometry]
  funext f h
  simp only [quadratureFaceMass, derivedFaceMass]
  congr 1
  apply Finset.sum_congr rfl
  intro r _
  apply Finset.sum_congr rfl
  intro s _
  rw [four_point_quadrature_moment_exact]

def localCurl : Matrix (Fin 4) (Fin 6) Q5 := !![
  q 1 0, q (-1) 0, q 0 0, q 1 0, q 0 0, q 0 0;
  q 1 0, q 0 0, q (-1) 0, q 0 0, q 1 0, q 0 0;
  q 0 0, q 1 0, q (-1) 0, q 0 0, q 0 0, q 1 0;
  q 0 0, q 0 0, q 0 0, q 1 0, q (-1) 0, q 1 0
]

def derivedStiffness : Matrix (Fin 6) (Fin 6) Q5 :=
  localCurl.transpose * faceMass * localCurl

theorem stability24_source_identity :
    stability24 = fun i j => 24 * edgeMass i j - derivedStiffness i j := by
  decide +kernel

theorem stability48_source_identity :
    stability48 = fun i j => 48 * edgeMass i j - derivedStiffness i j := by
  decide +kernel

#print axioms edge_table_source_bound
#print axioms face_table_source_bound
#print axioms source_tetrahedra_congruent
#print axioms source_tetrahedra_signedDeterminant
#print axioms edgeMass_from_source_geometry
#print axioms faceMass_from_source_geometry
#print axioms four_point_quadrature_moment_exact
#print axioms edgeMass_from_four_point_quadrature
#print axioms faceMass_from_four_point_quadrature
#print axioms stability24_source_identity

end OPH.WhitneySourceGeometry
