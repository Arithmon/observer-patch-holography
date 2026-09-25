import PortFrameGram

namespace OPH.GaloisPortFrames

/-!
# Two incidence-derived port Gram tables

Both branches are defined directly from adjacency and antipodes. The positive
specialization is compared with `PortFrameGram.g5` after the two-table
definition. Integer pairs represent `a + b sqrt(5)`; the finite identities
are checked by kernel reduction. The exact spherical-degree calculation is
in the independent Python certificate, not formalized here. No support
attachment or physical selection is asserted.
-/

/-- Coefficients in `Z[sqrt(5)]`. -/
abbrev Z5 := ℤ × ℤ

/-- Field conjugation on the integer coefficient ring. -/
def sigma (x : Z5) : Z5 := (x.1, -x.2)

/-- Multiplication in the quadratic coefficient ring. -/
def mul (x y : Z5) : Z5 :=
  (x.1 * y.1 + 5 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)

/-- Sum of coefficient pairs. -/
def sumZ (xs : List Z5) : Z5 :=
  xs.foldr (fun x s => (x.1 + s.1, x.2 + s.2)) (0, 0)

/-- `false` denotes the positive branch; `true` its conjugate. -/
def branchSign (conjugate : Bool) : ℤ := if conjugate then -1 else 1

/-- Five times the normalized Gram table, independently for each branch. -/
def gram5 (conjugate : Bool) (i j : Fin 12) : Z5 :=
  if i = j then (5, 0)
  else if OPH.PortFrameGram.adj i j then (0, branchSign conjugate)
  else if j = OPH.PortFrameGram.antipode i then (-5, 0)
  else (0, -branchSign conjugate)

/-- Matrix multiplication of the scaled Gram table with itself. -/
def square (c : Bool) (i j : Fin 12) : Z5 :=
  sumZ ((List.finRange 12).map fun k => mul (gram5 c i k) (gram5 c k j))

/-- Left multiplication by the committed adjacency. -/
def adjacencyProduct (c : Bool) (i j : Fin 12) : Z5 :=
  sumZ ((List.finRange 12).map fun k =>
    if OPH.PortFrameGram.adj i k then gram5 c k j else (0, 0))

/-- Left multiplication by `L = 5 I - A`. -/
def laplacianProduct (c : Bool) (i j : Fin 12) : Z5 :=
  (5 * (gram5 c i j).1 - (adjacencyProduct c i j).1,
   5 * (gram5 c i j).2 - (adjacencyProduct c i j).2)

/-- The declared positive table is the positive specialization. -/
theorem positive_specialization :
    ∀ i j : Fin 12, gram5 false i j = OPH.PortFrameGram.g5 i j := by decide

/-- Conjugation exchanges the complete two-table family. -/
theorem galois_exchange :
    ∀ c : Bool, ∀ i j : Fin 12, sigma (gram5 c i j) = gram5 (!c) i j := by decide

/-- Each branch is symmetric. -/
theorem gram_symmetric :
    ∀ c : Bool, ∀ i j : Fin 12, gram5 c i j = gram5 c j i := by decide

/-- Each normalized table has trace twelve. -/
theorem trace_twelve :
    ∀ c : Bool, sumZ ((List.finRange 12).map fun i => gram5 c i i) = (60, 0) := by
  decide

set_option maxRecDepth 8192 in
set_option maxHeartbeats 4000000 in
/-- `(5G)^2 = 20(5G)`, hence `G^2 = 4G`, for both branches. -/
theorem gram_square :
    ∀ c : Bool, ∀ i j : Fin 12,
      square c i j = (20 * (gram5 c i j).1, 20 * (gram5 c i j).2) := by decide

set_option maxHeartbeats 4000000 in
/-- The two column spaces have adjacency eigenvalues `+sqrt(5)` and `-sqrt(5)`. -/
theorem adjacency_eigenvalues :
    ∀ c : Bool, ∀ i j : Fin 12,
      adjacencyProduct c i j = mul (0, branchSign c) (gram5 c i j) := by decide

set_option maxHeartbeats 4000000 in
/-- The associated Laplacian eigenvalues are `5-sqrt(5)` and `5+sqrt(5)`. -/
theorem laplacian_eigenvalues :
    ∀ c : Bool, ∀ i j : Fin 12,
      laplacianProduct c i j = mul (5, -branchSign c) (gram5 c i j) := by decide

/-- A one-step rotation of the link at port zero. -/
def orderFive : Fin 12 → Fin 12 := ![0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11]

/-- The selected combinatorial action has exact order five. -/
theorem order_five :
    (∀ i : Fin 12, orderFive (orderFive (orderFive (orderFive (orderFive i)))) = i) ∧
      orderFive 1 ≠ 1 := by decide

/-- The action preserves the committed adjacency. -/
theorem orderFive_adjacency :
    ∀ i j : Fin 12, OPH.PortFrameGram.adj (orderFive i) (orderFive j) =
      OPH.PortFrameGram.adj i j := by decide

/-- Four times the exact carrier rotation for each branch. -/
def rotation4 (c : Bool) : Fin 3 → Fin 3 → Z5 :=
  ![![(-1, branchSign c), (1, branchSign c), (-2, 0)],
    ![(-1, -branchSign c), (2, 0), (-1, branchSign c)],
    ![(2, 0), (-1, branchSign c), (1, branchSign c)]]

/-- Twice the cyclic coordinates of each port, in the two embeddings. -/
def vector2 (c : Bool) : Fin 12 → Fin 3 → Z5 :=
  let p : Z5 := (1, branchSign c)
  let n : Z5 := (-1, -branchSign c)
  ![![(0,0),(2,0),p], ![(2,0),p,(0,0)], ![p,(0,0),(2,0)],
    ![(-2,0),p,(0,0)], ![(0,0),(-2,0),p], ![p,(0,0),(-2,0)],
    ![n,(0,0),(2,0)], ![(0,0),(2,0),n], ![(2,0),n,(0,0)],
    ![n,(0,0),(-2,0)], ![(-2,0),n,(0,0)], ![(0,0),(-2,0),n]]

/-- The matrices implement the same order-five port action in both frames. -/
theorem rotation_realizes_action :
    ∀ c : Bool, ∀ i : Fin 12, ∀ k : Fin 3,
      sumZ ((List.finRange 3).map fun j => mul (rotation4 c k j) (vector2 c i j)) =
        (4 * (vector2 c (orderFive i) k).1, 4 * (vector2 c (orderFive i) k).2) := by
  decide

/-- The trace pair is `4 phi = 2+2 sqrt(5)` and its conjugate. -/
theorem orderFive_trace_pair :
    ∀ c : Bool,
      sumZ ((List.finRange 3).map fun i => rotation4 c i i) =
        (2, 2 * branchSign c) := by decide

/-- Both displayed matrices are orthogonal after division by four. -/
theorem rotation_orthogonal :
    ∀ c : Bool, ∀ i j : Fin 3,
      sumZ ((List.finRange 3).map fun k => mul (rotation4 c i k) (rotation4 c j k)) =
        (if i = j then (16, 0) else (0, 0)) := by decide

#print axioms positive_specialization
#print axioms galois_exchange
#print axioms gram_symmetric
#print axioms trace_twelve
#print axioms gram_square
#print axioms adjacency_eigenvalues
#print axioms laplacian_eigenvalues
#print axioms order_five
#print axioms orderFive_adjacency
#print axioms rotation_realizes_action
#print axioms orderFive_trace_pair
#print axioms rotation_orthogonal

end OPH.GaloisPortFrames
