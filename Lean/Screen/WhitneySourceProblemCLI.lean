import WhitneyCertificateCLI
import WhitneyExactSourceProblem

set_option autoImplicit false

namespace OPH.WhitneySourceProblemCLI

open Lean
open OPH.WhitneyFiniteCertificate
open OPH.WhitneyCertificateCLI
open OPH.WhitneyExactSourceProblem

def q5Rows {m n : Nat} (matrix : Matrix (Fin m) (Fin n) Q5) :
    Array (Array Q5) :=
  Array.ofFn fun i => Array.ofFn fun j => matrix i j

def wireRows (rows : Array (Array Q5)) :
    Array (Array WireQ5) :=
  rows.map fun row => row.map q5Wire

def wireMatrix {n : Nat} (matrix : Matrix (Fin n) (Fin n) Q5) :
    Array (Array WireQ5) := wireRows (q5Rows matrix)

def edgeReindex : Array WireSignedIndex :=
  #[⟨1, -1⟩, ⟨0, 1⟩, ⟨2, 1⟩, ⟨3, 1⟩, ⟨4, 1⟩, ⟨5, 1⟩]

def faceReindex : Array WireSignedIndex :=
  #[⟨1, -1⟩, ⟨0, 1⟩, ⟨2, 1⟩, ⟨3, 1⟩]

def sourceAssembly : WireAssembly :=
  let gradient := q5Rows localGradient
  let curl := q5Rows localCurl
  let edge := q5Rows edgeMass
  let face := q5Rows faceMass
  { gradient := wireRows gradient
    curl := wireRows curl
    edgeMass := wireRows edge
    faceMass := wireRows face
    edgeReindex := edgeReindex
    faceReindex := faceReindex
    transportedGradient := wireRows (transportGradient gradient edgeReindex)
    transportedCurl := wireRows (transportCurl curl faceReindex edgeReindex)
    transportedEdgeMass := wireRows (transportSquare edge edgeReindex)
    transportedFaceMass := wireRows (transportSquare face faceReindex) }

def sourceProblem (sourceSha256 : String) : WireProblem :=
  ⟨"oph.whitney.exact-problem.v2", sourceSha256, sourceAssembly, #[
    ⟨"edge_mass", wireMatrix edgeMass, none⟩,
    ⟨"face_mass", wireMatrix faceMass, none⟩,
    ⟨"stability_24", wireMatrix stability24, none⟩,
    ⟨"stability_48", wireMatrix stability48, none⟩]⟩

def emitSourceProblem (sourceSha256 : String) (output : System.FilePath) : IO Unit :=
  IO.FS.writeFile output (Json.pretty (toJson (sourceProblem sourceSha256)) ++ "\n")

end OPH.WhitneySourceProblemCLI
