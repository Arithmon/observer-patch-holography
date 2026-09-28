import WhitneyAlgebraicLDL
import Mathlib.LinearAlgebra.Matrix.Reindex

set_option autoImplicit false
set_option maxHeartbeats 4000000
set_option maxRecDepth 8192

open scoped Matrix

namespace OPH.WhitneyCertificatePipeline

open OPH.WhitneyFiniteCertificate
open OPH.WhitneyAlgebraicLDL

/-- The trusted checker exposed by the executable certification pipeline. -/
def checkLDL {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (certificate : Certificate n) : Bool :=
  ldlCheck A certificate

theorem checkLDL_sound {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    {certificate : Certificate n} (h : checkLDL A certificate = true) :
    (evalMatrix A).PosDef :=
  ldl_checked_sound h

/-- Mathematical producer used for the completeness statement.  The command-line
producer is intentionally untrusted and its output is accepted only after
`checkLDL` is reduced by the kernel. -/
noncomputable def produceLDL {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (hA : (evalMatrix A).PosDef) : Certificate n :=
  Classical.choose (posDef_ldl_complete A hA)

theorem produceLDL_complete {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (hA : (evalMatrix A).PosDef) :
    checkLDL A (produceLDL A hA) = true :=
  Classical.choose_spec (posDef_ldl_complete A hA)

theorem certificate_exists_iff_posDef {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) :
    (∃ certificate, checkLDL A certificate = true) ↔ (evalMatrix A).PosDef := by
  simpa [checkLDL] using (posDef_iff_exists_ldl_certificate A).symm

/-- Exact congruence data for an arbitrary lawful signed reexpression.  A
signed permutation is the sparse special case, while the definition also
supports any exact invertible change of coordinates. -/
structure LawfulSignedReexpression {n : Nat}
    (source target : Matrix (Fin n) (Fin n) Q5) where
  forward : Matrix (Fin n) (Fin n) Q5
  backward : Matrix (Fin n) (Fin n) Q5
  forward_backward : forward * backward = 1
  backward_forward : backward * forward = 1
  target_eq : target = forward.transpose * source * forward

def LawfulSignedReexpression.identity {n : Nat}
    (source : Matrix (Fin n) (Fin n) Q5) :
    LawfulSignedReexpression source source where
  forward := 1
  backward := 1
  forward_backward := by simp
  backward_forward := by simp
  target_eq := by simp

def LawfulSignedReexpression.trans {n : Nat}
    {source middle target : Matrix (Fin n) (Fin n) Q5}
    (first : LawfulSignedReexpression source middle)
    (second : LawfulSignedReexpression middle target) :
    LawfulSignedReexpression source target := by
  rcases first with ⟨firstForward, firstBackward, hfb, hbf, hMiddle⟩
  rcases second with ⟨secondForward, secondBackward, hsf, hsb, hTarget⟩
  refine {
    forward := firstForward * secondForward
    backward := secondBackward * firstBackward
    forward_backward := ?_
    backward_forward := ?_
    target_eq := ?_ }
  · calc
      (firstForward * secondForward) * (secondBackward * firstBackward) =
          firstForward * (secondForward * secondBackward) * firstBackward := by
            noncomm_ring
      _ = 1 := by rw [hsf]; simpa using hfb
  · calc
      (secondBackward * firstBackward) * (firstForward * secondForward) =
          secondBackward * (firstBackward * firstForward) * secondForward := by
            noncomm_ring
      _ = 1 := by rw [hbf]; simpa using hsb
  · calc
      target = secondForward.transpose * middle * secondForward := hTarget
      _ = secondForward.transpose *
          (firstForward.transpose * source * firstForward) * secondForward := by
            rw [hMiddle]
      _ = (firstForward * secondForward).transpose * source *
          (firstForward * secondForward) := by
            rw [Matrix.transpose_mul]
            noncomm_ring

theorem LawfulSignedReexpression.forward_injective {n : Nat}
    {source target : Matrix (Fin n) (Fin n) Q5}
    (transport : LawfulSignedReexpression source target) :
    Function.Injective (evalMatrix transport.forward).mulVec := by
  rcases transport with ⟨forward, backward, hfb, hbf, htarget⟩
  intro x y hxy
  change evalMatrix forward *ᵥ x = evalMatrix forward *ᵥ y at hxy
  have h := congrArg (fun z => evalMatrix backward *ᵥ z) hxy
  dsimp at h
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
    ← evalMatrix_mul, hbf] at h
  simpa [evalMatrix] using h

theorem LawfulSignedReexpression.backward_injective {n : Nat}
    {source target : Matrix (Fin n) (Fin n) Q5}
    (transport : LawfulSignedReexpression source target) :
    Function.Injective (evalMatrix transport.backward).mulVec := by
  rcases transport with ⟨forward, backward, hfb, hbf, htarget⟩
  intro x y hxy
  change evalMatrix backward *ᵥ x = evalMatrix backward *ᵥ y at hxy
  have h := congrArg (fun z => evalMatrix forward *ᵥ z) hxy
  dsimp at h
  rw [Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
    ← evalMatrix_mul, hfb] at h
  simpa [evalMatrix] using h

theorem LawfulSignedReexpression.source_eq {n : Nat}
    {source target : Matrix (Fin n) (Fin n) Q5}
    (transport : LawfulSignedReexpression source target) :
    source = transport.backward.transpose * target * transport.backward := by
  rcases transport with ⟨forward, backward, hfb, hbf, htarget⟩
  calc
    source = (forward * backward).transpose * source *
        (forward * backward) := by rw [hfb]; simp
    _ = backward.transpose * (forward.transpose * source * forward) * backward := by
          rw [Matrix.transpose_mul]
          noncomm_ring
    _ = backward.transpose * target * backward := by
          rw [← htarget]

theorem LawfulSignedReexpression.posDef_iff {n : Nat}
    {source target : Matrix (Fin n) (Fin n) Q5}
    (transport : LawfulSignedReexpression source target) :
    (evalMatrix target).PosDef ↔ (evalMatrix source).PosDef := by
  constructor
  · intro hTarget
    rw [transport.source_eq, evalMatrix_mul, evalMatrix_mul,
      evalMatrix_transpose]
    exact hTarget.conjTranspose_mul_mul_same transport.backward_injective
  · intro hSource
    rw [transport.target_eq, evalMatrix_mul, evalMatrix_mul,
      evalMatrix_transpose]
    exact hSource.conjTranspose_mul_mul_same transport.forward_injective

theorem certificate_exists_transport_iff {n : Nat}
    {source target : Matrix (Fin n) (Fin n) Q5}
    (transport : LawfulSignedReexpression source target) :
    (∃ certificate, checkLDL target certificate = true) ↔
      ∃ certificate, checkLDL source certificate = true := by
  rw [certificate_exists_iff_posDef, certificate_exists_iff_posDef,
    transport.posDef_iff]

/-- Matrix reindexing is the exact coordinate action underlying signed source
reexpressions after signs have been absorbed into the linear equivalence. -/
def transportMatrix {n m : Type*} (e : n ≃ m) (A : Matrix n n Q5) : Matrix m m Q5 :=
  Matrix.reindex e e A

@[simp] theorem transport_identity {n : Type*} (A : Matrix n n Q5) :
    transportMatrix (Equiv.refl n) A = A := by
  exact Matrix.reindex_refl_refl A

theorem transport_composition {n m k : Type*} (e : n ≃ m) (f : m ≃ k)
    (A : Matrix n n Q5) :
    transportMatrix f (transportMatrix e A) = transportMatrix (e.trans f) A := by
  rfl

/-- A source transport is lawful when it supplies a two-sided coordinate
equivalence and the transported exact matrix.  D, C, orientation, assembly and
consumer compatibility are supplied by the source-level specialization. -/
structure LawfulReexpression {n m : Type*} where
  coordinate : n ≃ m

def LawfulReexpression.mapMatrix {n m : Type*} (t : LawfulReexpression (n := n) (m := m))
    (A : Matrix n n Q5) : Matrix m m Q5 :=
  transportMatrix t.coordinate A

@[simp] theorem LawfulReexpression.mapMatrix_identity {n : Type*}
    (A : Matrix n n Q5) :
    (LawfulReexpression.mk (Equiv.refl n)).mapMatrix A = A :=
  transport_identity A

theorem LawfulReexpression.mapMatrix_composition {n m k : Type*}
    (first : LawfulReexpression (n := n) (m := m))
    (second : LawfulReexpression (n := m) (m := k)) (A : Matrix n n Q5) :
    second.mapMatrix (first.mapMatrix A) =
      (LawfulReexpression.mk (first.coordinate.trans second.coordinate)).mapMatrix A :=
  transport_composition first.coordinate second.coordinate A

#print axioms checkLDL_sound
#print axioms produceLDL_complete
#print axioms certificate_exists_iff_posDef
#print axioms certificate_exists_transport_iff
#print axioms transport_identity
#print axioms transport_composition

end OPH.WhitneyCertificatePipeline
