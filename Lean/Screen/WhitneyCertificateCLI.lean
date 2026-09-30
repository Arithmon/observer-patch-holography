import WhitneyCertificatePipeline
import Lean.Data.Json.FromToJson
import Lean.Data.Json.Printer

set_option autoImplicit false
set_option maxRecDepth 8192

namespace OPH.WhitneyCertificateCLI

open Lean
open OPH.WhitneyFiniteCertificate
open OPH.WhitneyAlgebraicLDL
open OPH.WhitneyCertificatePipeline

structure WireQ5 where
  a : String
  b : String
  deriving DecidableEq, FromJson, ToJson

structure WireCertificate where
  lower : Array (Array WireQ5)
  diagonal : Array WireQ5
  deriving DecidableEq, FromJson, ToJson

structure WireTarget where
  name : String
  matrix : Array (Array WireQ5)
  certificate : Option WireCertificate := none
  deriving DecidableEq, FromJson, ToJson

structure WireSignedIndex where
  source : Nat
  sign : Int
  deriving DecidableEq, FromJson, ToJson

structure WireAssembly where
  gradient : Array (Array WireQ5)
  curl : Array (Array WireQ5)
  edgeMass : Array (Array WireQ5)
  faceMass : Array (Array WireQ5)
  edgeReindex : Array WireSignedIndex
  faceReindex : Array WireSignedIndex
  transportedGradient : Array (Array WireQ5)
  transportedCurl : Array (Array WireQ5)
  transportedEdgeMass : Array (Array WireQ5)
  transportedFaceMass : Array (Array WireQ5)
  deriving DecidableEq, FromJson, ToJson

structure WireProblem where
  schema : String
  sourceSha256 : String
  assembly : WireAssembly
  targets : Array WireTarget
  deriving DecidableEq, FromJson, ToJson

structure WireCertifiedTarget where
  name : String
  matrix : Array (Array WireQ5)
  certificate : WireCertificate
  deriving DecidableEq, FromJson, ToJson

structure WireResult where
  schema : String
  sourceSha256 : String
  targets : Array WireCertifiedTarget
  deriving DecidableEq, FromJson, ToJson

def parseRat (value : String) : Except String Rat := do
  match value.splitOn "/" with
  | [numerator] =>
      let some n := numerator.toInt? | throw s!"invalid rational numerator: {value}"
      pure (Rat.ofInt n)
  | [numerator, denominator] =>
      let some n := numerator.toInt? | throw s!"invalid rational numerator: {value}"
      let some d := denominator.toInt? | throw s!"invalid rational denominator: {value}"
      if d ≤ 0 then throw s!"nonpositive rational denominator: {value}"
      pure (Rat.divInt n d)
  | _ => throw s!"invalid rational syntax: {value}"

def WireQ5.parse (value : WireQ5) : Except String Q5 := do
  let a ← parseRat value.a
  let b ← parseRat value.b
  if toString a != value.a || toString b != value.b then
    throw "noncanonical exact scalar"
  pure (QuadraticAlgebra.mk a b)

def q5Wire (value : Q5) : WireQ5 :=
  ⟨toString value.re, toString value.im⟩

def encodeExact (problem : WireProblem) : ByteArray :=
  (Json.compress (toJson problem) ++ "\n").toUTF8

def decodeExact (bytes : ByteArray) : Except String WireProblem := do
  let some text := String.fromUTF8? bytes | throw "exact problem is not UTF-8"
  let json ← Json.parse text
  let problem ← fromJson? json
  if encodeExact problem != bytes then
    throw "noncanonical, duplicate, or unsupported exact-problem fields"
  pure problem

def parseMatrix (rows : Array (Array WireQ5)) : Except String (Array (Array Q5)) := do
  let n := rows.size
  if n = 0 then throw "empty exact matrix"
  let mut result := #[]
  for row in rows do
    if row.size != n then throw "exact matrix is not square"
    result := result.push (← row.mapM WireQ5.parse)
  pure result

def parseRectMatrix (rows : Array (Array WireQ5)) (height width : Nat) :
    Except String (Array (Array Q5)) := do
  if rows.size != height then throw "exact rectangular matrix row count"
  let mut result := #[]
  for row in rows do
    if row.size != width then throw "exact rectangular matrix column count"
    result := result.push (← row.mapM WireQ5.parse)
  pure result

def validateSignedMap (name : String) (mapping : Array WireSignedIndex) (n : Nat) :
    Except String Unit := do
  if mapping.size != n then throw s!"{name} signed reindex dimension"
  let mut seen := Array.replicate n false
  for entry in mapping do
    if entry.source >= n then throw s!"{name} signed reindex range"
    if entry.sign != 1 && entry.sign != -1 then throw s!"{name} signed reindex sign"
    if seen[entry.source]! then throw s!"{name} signed reindex is not injective"
    seen := seen.set! entry.source true

def signed (left right : Int) (value : Q5) : Q5 :=
  if left = right then value else -value

def transportSquare (matrix : Array (Array Q5))
    (mapping : Array WireSignedIndex) : Array (Array Q5) :=
  mapping.map fun i => mapping.map fun j =>
    signed i.sign j.sign ((matrix[i.source]!)[j.source]!)

def transportGradient (gradient : Array (Array Q5))
    (edges : Array WireSignedIndex) : Array (Array Q5) :=
  edges.map fun edge => (gradient[edge.source]!).map fun value =>
    if edge.sign = 1 then value else -value

def transportCurl (curl : Array (Array Q5)) (faces edges : Array WireSignedIndex) :
    Array (Array Q5) :=
  faces.map fun face => edges.map fun edge =>
    signed face.sign edge.sign ((curl[face.source]!)[edge.source]!)

def chainEntry (curl gradient : Array (Array Q5)) (face vertex : Nat) : Q5 :=
  (Array.range gradient.size).foldl
    (fun total edge => total + (curl[face]!)[edge]! * (gradient[edge]!)[vertex]!) 0

def stiffnessEntry (curl faceMass : Array (Array Q5)) (i j : Nat) : Q5 :=
  (Array.range faceMass.size).foldl (fun total a =>
    total + (Array.range faceMass.size).foldl (fun inner b =>
      inner + (curl[a]!)[i]! * (faceMass[a]!)[b]! * (curl[b]!)[j]!) 0) 0

def stabilityMatrix (bound : Nat) (edgeMass curl faceMass : Array (Array Q5)) :
    Array (Array Q5) :=
  (Array.range edgeMass.size).map fun i =>
    (Array.range edgeMass.size).map fun j =>
      (QuadraticAlgebra.mk (bound : Rat) 0) * (edgeMass[i]!)[j]! -
        stiffnessEntry curl faceMass i j

def requiredTarget (problem : WireProblem) (name : String) : Except String WireTarget := do
  let hits := problem.targets.filter fun target => target.name = name
  if hits.size != 1 then throw s!"required exact target census: {name}"
  let some target := hits[0]? | throw s!"missing exact target: {name}"
  pure target

def validateAssembly (problem : WireProblem) : Except String Unit := do
  let wire := problem.assembly
  let edgeMass ← parseMatrix wire.edgeMass
  let faceMass ← parseMatrix wire.faceMass
  let nEdge := edgeMass.size
  let nFace := faceMass.size
  if nEdge = 0 || nFace = 0 || wire.gradient.isEmpty then
    throw "empty source assembly"
  let nVertex := wire.gradient[0]!.size
  if nVertex = 0 then throw "empty source vertex space"
  let gradient ← parseRectMatrix wire.gradient nEdge nVertex
  let curl ← parseRectMatrix wire.curl nFace nEdge
  validateSignedMap "edge" wire.edgeReindex nEdge
  validateSignedMap "face" wire.faceReindex nFace
  let gradient' ← parseRectMatrix wire.transportedGradient nEdge nVertex
  let curl' ← parseRectMatrix wire.transportedCurl nFace nEdge
  let edgeMass' ← parseMatrix wire.transportedEdgeMass
  let faceMass' ← parseMatrix wire.transportedFaceMass
  if edgeMass'.size != nEdge || faceMass'.size != nFace then
    throw "transported source mass dimension"
  for face in [:nFace] do
    for vertex in [:nVertex] do
      if chainEntry curl gradient face vertex != 0 then throw "source curl-gradient law"
      if chainEntry curl' gradient' face vertex != 0 then
        throw "transported curl-gradient law"
  if gradient' != transportGradient gradient wire.edgeReindex then
    throw "signed gradient transport"
  if curl' != transportCurl curl wire.faceReindex wire.edgeReindex then
    throw "signed curl transport"
  if edgeMass' != transportSquare edgeMass wire.edgeReindex then
    throw "signed edge-mass transport"
  if faceMass' != transportSquare faceMass wire.faceReindex then
    throw "signed face-mass transport"
  let expected := #[
    ("edge_mass", edgeMass), ("face_mass", faceMass),
    ("stability_24", stabilityMatrix 24 edgeMass curl faceMass),
    ("stability_48", stabilityMatrix 48 edgeMass curl faceMass)]
  for item in expected do
    let target ← requiredTarget problem item.1
    if (← parseMatrix target.matrix) != item.2 then
      throw s!"source assembly target mismatch: {item.1}"

def matrixOfRows (rows : Array (Array Q5)) : Matrix (Fin rows.size) (Fin rows.size) Q5 :=
  fun i j => (rows[i.1]!)[j.1]!

structure RuntimeCertificate where
  lower : Array (Array Q5)
  diagonal : Array Q5

def RuntimeCertificate.asCertificate (certificate : RuntimeCertificate)
    (n : Nat) : Certificate n :=
  ⟨fun i j => (certificate.lower[i.1]!)[j.1]!, fun i => certificate.diagonal[i.1]!⟩

def produceRuntimeLDL (matrix : Array (Array Q5)) : Except String RuntimeCertificate := do
  let n := matrix.size
  if n = 0 then throw "empty exact matrix"
  for row in matrix do
    if row.size != n then throw "exact matrix is not square"
  let mut lower := Array.replicate n (Array.replicate n 0)
  let mut diagonal := Array.replicate n 0
  for i in [:n] do
    lower := lower.set! i ((lower[i]!).set! i 1)
    let mut pivot := (matrix[i]!)[i]!
    for k in [:i] do
      pivot := pivot - (lower[i]!)[k]! * (lower[i]!)[k]! * diagonal[k]!
    if positive pivot then
      diagonal := diagonal.set! i pivot
    else
      throw s!"nonpositive exact LDL pivot at {i}"
    for j in [i + 1:n] do
      let mut numerator := (matrix[j]!)[i]!
      for k in [:i] do
        numerator := numerator - (lower[j]!)[k]! * (lower[i]!)[k]! * diagonal[k]!
      lower := lower.set! j ((lower[j]!).set! i (numerator / pivot))
  pure ⟨lower, diagonal⟩

def RuntimeCertificate.toWire (certificate : RuntimeCertificate) : WireCertificate :=
  ⟨certificate.lower.map (fun row => row.map q5Wire), certificate.diagonal.map q5Wire⟩

def WireCertificate.parse (certificate : WireCertificate) (n : Nat) :
    Except String RuntimeCertificate := do
  if certificate.lower.size != n || certificate.diagonal.size != n then
    throw "certificate dimension mismatch"
  let lower ← parseMatrix certificate.lower
  let diagonal ← certificate.diagonal.mapM WireQ5.parse
  pure ⟨lower, diagonal⟩

def certifyTarget (target : WireTarget) : Except String WireCertifiedTarget := do
  let rows ← parseMatrix target.matrix
  let certificate ← produceRuntimeLDL rows
  let checked := checkLDL (matrixOfRows rows)
    (certificate.asCertificate rows.size)
  if checked then
    pure ⟨target.name, target.matrix, certificate.toWire⟩
  else
    throw s!"generated certificate failed trusted checker: {target.name}"

abbrev ScalarData := Int × Nat × Int × Nat

def scalarData (z : Q5) : ScalarData := (z.re.num, z.re.den, z.im.num, z.im.den)

def targetData {n : Nat} (matrix : Matrix (Fin n) (Fin n) Q5)
    (certificate : Certificate n) : Array (Array (Array ScalarData)) :=
  #[Array.ofFn (fun i => Array.ofFn (fun j => scalarData (matrix i j))),
    Array.ofFn (fun i => Array.ofFn (fun j => scalarData (certificate.factor i j))),
    #[Array.ofFn (fun i => scalarData (certificate.diagonal i))]]

def identifier (value : String) : String :=
  value.foldl (fun out c => if c.isAlphanum then out.push c else out.push '_') ""

def ratLiteral (value : Rat) : String :=
  if value.den = 1 then s!"({value.num} : Rat)"
  else s!"(({value.num} : Rat) / {value.den})"

def q5Literal (value : Q5) : String :=
  s!"QuadraticAlgebra.mk {ratLiteral value.re} {ratLiteral value.im}"

def arrayLiteral (values : Array String) : String :=
  "#[" ++ String.intercalate ", " values.toList ++ "]"

def matrixLiteral (rows : Array (Array Q5)) : String :=
  arrayLiteral (rows.map (fun row => arrayLiteral (row.map q5Literal)))

def certificateLiteral (certificate : RuntimeCertificate) : String :=
  s!"⟨matrixOfRows {matrixLiteral certificate.lower}, fun i => ({arrayLiteral (certificate.diagonal.map q5Literal)})[i.1]!⟩"

def renderWitness (sourceSha256 : String)
    (targets : Array (WireTarget × Array (Array Q5) × RuntimeCertificate)) : String :=
  let header := "import WhitneyCertificateCLI\n\nset_option autoImplicit false\n" ++
    "set_option maxHeartbeats 8000000\nset_option maxRecDepth 32768\n\n" ++
    "open OPH.WhitneyFiniteCertificate\nopen OPH.WhitneyAlgebraicLDL\n" ++
    "open OPH.WhitneyCertificatePipeline\nopen OPH.WhitneyCertificateCLI\n\n" ++
    "namespace OPH.WhitneyGeneratedCertificate\n\n" ++
    s!"def sourceSha256 : String := {repr sourceSha256}\n\n"
  let body := targets.foldl (fun output item =>
    let target := item.1
    let rows := item.2.1
    let certificate := item.2.2
    let name := identifier target.name
    output ++ s!"def {name}Matrix : Matrix (Fin {rows.size}) (Fin {rows.size}) Q5 :=\n" ++
      s!"  matrixOfRows {matrixLiteral rows}\n" ++
      s!"def {name}Certificate : Certificate {rows.size} :=\n" ++
      s!"  {certificateLiteral certificate}\n" ++
      s!"theorem {name}_checked : checkLDL {name}Matrix {name}Certificate = true := by\n" ++
      "  decide +kernel\n" ++
      s!"theorem {name}_posDef : (evalMatrix {name}Matrix).PosDef :=\n" ++
      s!"  checkLDL_sound {name}_checked\n" ++
      s!"#print axioms {name}_checked\n#print axioms {name}_posDef\n" ++
      s!"def {name}Data : Array (Array (Array ScalarData)) := " ++
      arrayLiteral (#[rows, certificate.lower, #[certificate.diagonal]].map (fun matrix =>
        arrayLiteral (matrix.map (fun row => arrayLiteral (row.map (fun z =>
          s!"({z.re.num}, {z.re.den}, {z.im.num}, {z.im.den})")))))) ++ "\n" ++
      s!"theorem {name}Data_bound : {name}Data = targetData {name}Matrix {name}Certificate := by\n" ++
      "  decide +kernel\n" ++ s!"#print axioms {name}Data_bound\n\n") ""
  header ++ body ++ "end OPH.WhitneyGeneratedCertificate\n"

def run (input output witness : System.FilePath) : IO Unit := do
  let bytes ← IO.FS.readBinFile input
  let problem ← IO.ofExcept (decodeExact bytes)
  if problem.schema != "oph.whitney.exact-problem.v2" then
    throw (IO.userError s!"unsupported exact problem schema: {problem.schema}")
  if problem.targets.isEmpty then throw (IO.userError "exact problem has no targets")
  let roundTrip ← IO.ofExcept (decodeExact (encodeExact problem))
  if roundTrip != problem then throw (IO.userError "decode_encode regression")
  IO.ofExcept (validateAssembly problem)
  let mut certified : Array WireCertifiedTarget := #[]
  let mut witnessInputs : Array (WireTarget × Array (Array Q5) × RuntimeCertificate) := #[]
  for target in problem.targets do
    let rows ← IO.ofExcept (parseMatrix target.matrix)
    let certificate ← IO.ofExcept (match target.certificate with
      | some supplied => supplied.parse rows.size
      | none => produceRuntimeLDL rows)
    let checked := checkLDL (matrixOfRows rows) (certificate.asCertificate rows.size)
    if !checked then throw (IO.userError s!"trusted checker rejected generated target: {target.name}")
    certified := certified.push ⟨target.name, target.matrix, certificate.toWire⟩
    witnessInputs := witnessInputs.push (target, rows, certificate)
  let result : WireResult :=
    ⟨"oph.whitney.exact-certificate.v1", problem.sourceSha256, certified⟩
  IO.FS.writeFile output (Json.pretty (toJson result) ++ "\n")
  IO.FS.writeFile witness (renderWitness problem.sourceSha256 witnessInputs)

def main (args : List String) : IO UInt32 := do
  match args with
  | [input, output, witness] =>
      try
        run input output witness
        pure 0
      catch error =>
        IO.eprintln error.toString
        pure 1
  | _ =>
      IO.eprintln "usage: whitneyCert <problem.json> <certificate.json> <witness.lean>"
      pure 2

end OPH.WhitneyCertificateCLI
