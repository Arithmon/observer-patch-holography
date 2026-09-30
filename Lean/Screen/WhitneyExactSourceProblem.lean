import WhitneyCertificatePipeline
import SeamCurrentEdge30Moment

set_option autoImplicit false
set_option maxHeartbeats 1000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix goldenRatio

namespace OPH.WhitneyExactSourceProblem

open OPH.WhitneyFiniteCertificate

def q (a b : ℚ) : Q5 := QuadraticAlgebra.mk a b
def phiQ5 : Q5 := q (1 / 2) (1 / 2)

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

def gradientGram : Matrix (Fin 4) (Fin 4) Q5 := !![
  q (21/2) (-9/2), q (-7/2) (3/2), q (-7/2) (3/2), q (-7/2) (3/2);
  q (-7/2) (3/2), q (3/2) (-1/2), q 1 (-1/2), q 1 (-1/2);
  q (-7/2) (3/2), q 1 (-1/2), q (3/2) (-1/2), q 1 (-1/2);
  q (-7/2) (3/2), q 1 (-1/2), q 1 (-1/2), q (3/2) (-1/2)
]

def localEdge : Fin 6 → Fin 4 × Fin 4 :=
  ![(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def localGradient : Matrix (Fin 6) (Fin 4) Q5 := fun edge vertex =>
  if vertex = (localEdge edge).1 then q (-1) 0
  else if vertex = (localEdge edge).2 then q 1 0
  else q 0 0

def localFace : Fin 4 → Fin 4 × Fin 4 × Fin 4 :=
  ![(0,1,2),(0,1,3),(0,2,3),(1,2,3)]

def moment (i j : Fin 4) : Q5 := if i = j then q (1/10) 0 else q (1/20) 0

def edgeMass : Matrix (Fin 6) (Fin 6) Q5 := fun e f =>
  let i := (localEdge e).1
  let j := (localEdge e).2
  let k := (localEdge f).1
  let l := (localEdge f).2
  moment i k * gradientGram j l - moment i l * gradientGram j k -
    moment j k * gradientGram i l + moment j l * gradientGram i k

def cycleFace (f : Fin 4) : Fin 3 → Fin 4 × Fin 4 × Fin 4 := fun r =>
  let p := localFace f
  match r.1 with
  | 0 => p
  | 1 => (p.2.1, p.2.2, p.1)
  | _ => (p.2.2, p.1, p.2.1)

def crossDot (i j k l : Fin 4) : Q5 :=
  gradientGram i k * gradientGram j l - gradientGram i l * gradientGram j k

def faceMass : Matrix (Fin 4) (Fin 4) Q5 := fun f h =>
  4 * ∑ r : Fin 3, ∑ s : Fin 3,
    moment (cycleFace f r).1 (cycleFace h s).1 *
      crossDot (cycleFace f r).2.1 (cycleFace f r).2.2
        (cycleFace h s).2.1 (cycleFace h s).2.2

def localCurl : Matrix (Fin 4) (Fin 6) Q5 := !![
  q 1 0, q (-1) 0, q 0 0, q 1 0, q 0 0, q 0 0;
  q 1 0, q 0 0, q (-1) 0, q 0 0, q 1 0, q 0 0;
  q 0 0, q 1 0, q (-1) 0, q 0 0, q 0 0, q 1 0;
  q 0 0, q 0 0, q 0 0, q 1 0, q (-1) 0, q 1 0
]

def stiffness : Matrix (Fin 6) (Fin 6) Q5 :=
  localCurl.transpose * faceMass * localCurl

theorem local_chain : localCurl * localGradient = 0 := by decide +kernel

def stability24 : Matrix (Fin 6) (Fin 6) Q5 := fun i j =>
  24 * edgeMass i j - stiffness i j

def stability48 : Matrix (Fin 6) (Fin 6) Q5 := fun i j =>
  48 * edgeMass i j - stiffness i j

/-- Finite exact source data carried by the executable problem.  Keeping D, C,
and both mass matrices in one value prevents a signed coordinate change from
silently updating only one part of the assembly. -/
structure ExactWhitneySource where
  gradient : Matrix (Fin 6) (Fin 4) Q5
  curl : Matrix (Fin 4) (Fin 6) Q5
  edgeMass : Matrix (Fin 6) (Fin 6) Q5
  faceMass : Matrix (Fin 4) (Fin 4) Q5
  chain : curl * gradient = 0

def source : ExactWhitneySource where
  gradient := localGradient
  curl := localCurl
  edgeMass := edgeMass
  faceMass := faceMass
  chain := local_chain

def ExactWhitneySource.stiffness (problem : ExactWhitneySource) :
    Matrix (Fin 6) (Fin 6) Q5 :=
  problem.curl.transpose * problem.faceMass * problem.curl

def ExactWhitneySource.stability (scale : Q5) (problem : ExactWhitneySource) :
    Matrix (Fin 6) (Fin 6) Q5 := fun i j =>
  scale * problem.edgeMass i j - problem.stiffness i j

/-- The nonidentity edge signed permutation emitted by
`WhitneySourceProblemCLI`.  Rows are target coordinates and columns are source
coordinates. -/
def runtimeEdgeForward : Matrix (Fin 6) (Fin 6) Q5 := !![
  q 0 0, q (-1) 0, q 0 0, q 0 0, q 0 0, q 0 0;
  q 1 0, q 0 0, q 0 0, q 0 0, q 0 0, q 0 0;
  q 0 0, q 0 0, q 1 0, q 0 0, q 0 0, q 0 0;
  q 0 0, q 0 0, q 0 0, q 1 0, q 0 0, q 0 0;
  q 0 0, q 0 0, q 0 0, q 0 0, q 1 0, q 0 0;
  q 0 0, q 0 0, q 0 0, q 0 0, q 0 0, q 1 0
]

/-- The nonidentity face signed permutation emitted by the runtime source
problem. -/
def runtimeFaceForward : Matrix (Fin 4) (Fin 4) Q5 := !![
  q 0 0, q (-1) 0, q 0 0, q 0 0;
  q 1 0, q 0 0, q 0 0, q 0 0;
  q 0 0, q 0 0, q 1 0, q 0 0;
  q 0 0, q 0 0, q 0 0, q 1 0
]

def runtimeTarget : ExactWhitneySource where
  gradient := runtimeEdgeForward * localGradient
  curl := runtimeFaceForward * localCurl * runtimeEdgeForward.transpose
  edgeMass := runtimeEdgeForward * edgeMass * runtimeEdgeForward.transpose
  faceMass := runtimeFaceForward * faceMass * runtimeFaceForward.transpose
  chain := by decide +kernel

/-- Exact source-level contract for a simultaneous reexpression of D, C and
the two mass matrices. -/
structure LawfulExactReexpression (source target : ExactWhitneySource) where
  edgeForward : Matrix (Fin 6) (Fin 6) Q5
  edgeBackward : Matrix (Fin 6) (Fin 6) Q5
  faceForward : Matrix (Fin 4) (Fin 4) Q5
  faceBackward : Matrix (Fin 4) (Fin 4) Q5
  edge_forward_backward : edgeForward * edgeBackward = 1
  edge_backward_forward : edgeBackward * edgeForward = 1
  face_forward_backward : faceForward * faceBackward = 1
  face_backward_forward : faceBackward * faceForward = 1
  gradient_transport : target.gradient = edgeForward * source.gradient
  curl_transport : target.curl = faceForward * source.curl * edgeBackward
  edgeMass_transport :
    target.edgeMass = edgeForward * source.edgeMass * edgeForward.transpose
  faceMass_transport :
    target.faceMass = faceForward * source.faceMass * faceForward.transpose

def LawfulExactReexpression.identity (problem : ExactWhitneySource) :
    LawfulExactReexpression problem problem where
  edgeForward := 1
  edgeBackward := 1
  faceForward := 1
  faceBackward := 1
  edge_forward_backward := by simp
  edge_backward_forward := by simp
  face_forward_backward := by simp
  face_backward_forward := by simp
  gradient_transport := by simp
  curl_transport := by simp
  edgeMass_transport := by simp
  faceMass_transport := by simp

def LawfulExactReexpression.trans {first middle last : ExactWhitneySource}
    (left : LawfulExactReexpression first middle)
    (right : LawfulExactReexpression middle last) :
    LawfulExactReexpression first last where
  edgeForward := right.edgeForward * left.edgeForward
  edgeBackward := left.edgeBackward * right.edgeBackward
  faceForward := right.faceForward * left.faceForward
  faceBackward := left.faceBackward * right.faceBackward
  edge_forward_backward := by
    calc
      (right.edgeForward * left.edgeForward) *
          (left.edgeBackward * right.edgeBackward) =
        right.edgeForward * (left.edgeForward * left.edgeBackward) *
          right.edgeBackward := by simp only [Matrix.mul_assoc]
      _ = 1 := by rw [left.edge_forward_backward]; simpa using right.edge_forward_backward
  edge_backward_forward := by
    calc
      (left.edgeBackward * right.edgeBackward) *
          (right.edgeForward * left.edgeForward) =
        left.edgeBackward * (right.edgeBackward * right.edgeForward) *
          left.edgeForward := by simp only [Matrix.mul_assoc]
      _ = 1 := by rw [right.edge_backward_forward]; simpa using left.edge_backward_forward
  face_forward_backward := by
    calc
      (right.faceForward * left.faceForward) *
          (left.faceBackward * right.faceBackward) =
        right.faceForward * (left.faceForward * left.faceBackward) *
          right.faceBackward := by simp only [Matrix.mul_assoc]
      _ = 1 := by rw [left.face_forward_backward]; simpa using right.face_forward_backward
  face_backward_forward := by
    calc
      (left.faceBackward * right.faceBackward) *
          (right.faceForward * left.faceForward) =
        left.faceBackward * (right.faceBackward * right.faceForward) *
          left.faceForward := by simp only [Matrix.mul_assoc]
      _ = 1 := by rw [right.face_backward_forward]; simpa using left.face_backward_forward
  gradient_transport := by
    rw [right.gradient_transport, left.gradient_transport, Matrix.mul_assoc]
  curl_transport := by
    rw [right.curl_transport, left.curl_transport]
    simp only [Matrix.mul_assoc]
  edgeMass_transport := by
    rw [right.edgeMass_transport, left.edgeMass_transport,
      Matrix.transpose_mul]
    simp only [Matrix.mul_assoc]
  faceMass_transport := by
    rw [right.faceMass_transport, left.faceMass_transport,
      Matrix.transpose_mul]
    simp only [Matrix.mul_assoc]

def runtimeLawfulReexpression : LawfulExactReexpression source runtimeTarget where
  edgeForward := runtimeEdgeForward
  edgeBackward := runtimeEdgeForward.transpose
  faceForward := runtimeFaceForward
  faceBackward := runtimeFaceForward.transpose
  edge_forward_backward := by decide +kernel
  edge_backward_forward := by decide +kernel
  face_forward_backward := by decide +kernel
  face_backward_forward := by decide +kernel
  gradient_transport := rfl
  curl_transport := rfl
  edgeMass_transport := rfl
  faceMass_transport := rfl

theorem runtimeTarget_stability24_transport :
    runtimeTarget.stability (q 24 0) =
      runtimeEdgeForward * source.stability (q 24 0) *
        runtimeEdgeForward.transpose := by
  decide +kernel

#print axioms vertexQ5_source_bound
#print axioms local_chain
#print axioms LawfulExactReexpression.trans
#print axioms runtimeTarget_stability24_transport

end OPH.WhitneyExactSourceProblem
