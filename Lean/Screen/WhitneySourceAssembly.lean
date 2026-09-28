import WhitneySourceGeometry

set_option autoImplicit false
set_option maxHeartbeats 2000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneySourceAssembly

open OPH.WhitneyFiniteCertificate
open OPH.WhitneyActualCertificates
open OPH.WhitneySourceGeometry
open OPH.ConeCochainBridge

noncomputable section

def Positive {V : Type*} [AddCommGroup V] [Module ℝ V]
    (B : LinearMap.BilinForm ℝ V) : Prop :=
  ∀ x : V, x ≠ 0 → 0 < B x x

def JointKernelTrivial {T V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (r : T → V →ₗ[ℝ] W) : Prop :=
  ∀ x : V, (∀ t : T, r t x = 0) → x = 0

def assemble {T V W : Type*} [Fintype T]
    [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W]
    (w : T → ℝ) (B : T → LinearMap.BilinForm ℝ W)
    (r : T → V →ₗ[ℝ] W) : LinearMap.BilinForm ℝ V :=
  ∑ t, w t • (B t).compl₁₂ (r t) (r t)

@[simp] theorem assemble_apply {T V W : Type*} [Fintype T]
    [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W]
    (w : T → ℝ) (B : T → LinearMap.BilinForm ℝ W)
    (r : T → V →ₗ[ℝ] W) (x y : V) :
    assemble w B r x y = ∑ t, w t * B t (r t x) (r t y) := by
  simp [assemble, LinearMap.compl₁₂_apply]

theorem Positive.nonneg {V : Type*} [AddCommGroup V] [Module ℝ V]
    {B : LinearMap.BilinForm ℝ V} (hB : Positive B) (x : V) : 0 ≤ B x x := by
  by_cases hx : x = 0
  · simp [hx]
  · exact le_of_lt (hB x hx)

theorem weighted_sum_eq_zero_iff {T : Type*} [Fintype T]
    (w f : T → ℝ) (hw : ∀ t, 0 < w t) (hf : ∀ t, 0 ≤ f t) :
    (∑ t, w t * f t) = 0 ↔ ∀ t, f t = 0 := by
  have hnonneg : ∀ t ∈ (Finset.univ : Finset T), 0 ≤ w t * f t :=
    fun t _ => mul_nonneg (le_of_lt (hw t)) (hf t)
  rw [Finset.sum_eq_zero_iff_of_nonneg hnonneg]
  constructor
  · intro h t
    exact (mul_eq_zero.mp (h t (Finset.mem_univ t))).resolve_left (ne_of_gt (hw t))
  · intro h t _
    rw [h t, mul_zero]

theorem assemble_eq_zero_iff {T V W : Type*} [Fintype T]
    [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W]
    (w : T → ℝ) (B : T → LinearMap.BilinForm ℝ W)
    (r : T → V →ₗ[ℝ] W) (hw : ∀ t, 0 < w t)
    (hB : ∀ t, Positive (B t)) (x : V) :
    assemble w B r x x = 0 ↔ ∀ t, r t x = 0 := by
  rw [assemble_apply, weighted_sum_eq_zero_iff w
    (fun t => B t (r t x) (r t x)) hw (fun t => (hB t).nonneg (r t x))]
  constructor
  · intro h t
    by_contra hne
    exact (ne_of_gt (hB t (r t x) hne)) (h t)
  · intro h t
    simp [h t]

theorem assemble_positive_iff {T V W : Type*} [Fintype T]
    [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W]
    (w : T → ℝ) (B : T → LinearMap.BilinForm ℝ W)
    (r : T → V →ₗ[ℝ] W) (hw : ∀ t, 0 < w t)
    (hB : ∀ t, Positive (B t)) :
    Positive (assemble w B r) ↔ JointKernelTrivial r := by
  constructor
  · intro hp x hx
    by_contra hne
    exact (ne_of_gt (hp x hne)) ((assemble_eq_zero_iff w B r hw hB x).mpr hx)
  · intro hk x hx
    have hn : 0 ≤ assemble w B r x x := by
      rw [assemble_apply]
      exact Finset.sum_nonneg (fun t _ =>
        mul_nonneg (le_of_lt (hw t)) ((hB t).nonneg (r t x)))
    exact lt_of_le_of_ne hn (fun heq =>
      hx (hk x ((assemble_eq_zero_iff w B r hw hB x).mp heq.symm)))

theorem assemble_symm {T V W : Type*} [Fintype T]
    [AddCommGroup V] [Module ℝ V] [AddCommGroup W] [Module ℝ W]
    (w : T → ℝ) (B : T → LinearMap.BilinForm ℝ W)
    (r : T → V →ₗ[ℝ] W) (hB : ∀ t, (B t).IsSymm) :
    (assemble w B r).IsSymm := by
  constructor
  intro x y
  simp only [assemble_apply]
  apply Finset.sum_congr rfl
  intro t _
  rw [(hB t).eq]

def realForm {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    LinearMap.BilinForm ℝ (Fin n → ℝ) := Matrix.toBilin' (evalMatrix A)

theorem realForm_positive {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    (hA : (evalMatrix A).PosDef) : Positive (realForm A) := by
  intro x hx
  rw [realForm, Matrix.toBilin'_apply']
  simpa using hA.dotProduct_mulVec_pos hx

theorem realForm_symm {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    (hA : (evalMatrix A).PosDef) : (realForm A).IsSymm := by
  constructor
  intro x y
  simp only [realForm, Matrix.toBilin'_apply]
  have hsymm : (evalMatrix A).IsSymm := by
    apply Matrix.IsSymm.ext
    intro i j
    have h := congrFun (congrFun hA.1.eq i) j
    simpa [Matrix.conjTranspose_apply] using h
  calc
    (∑ i, ∑ j, x i * evalMatrix A i j * y j) =
        ∑ j, ∑ i, y j * evalMatrix A j i * x i := by
          rw [Finset.sum_comm]
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          rw [hsymm.apply]
          ring
    _ = ∑ i, ∑ j, y i * evalMatrix A i j * x j := rfl

def signR (positive : Bool) : ℝ := if positive then 1 else -1

theorem signR_ne_zero (b : Bool) : signR b ≠ 0 := by
  cases b <;> norm_num [signR]

def coneOneIndexEquiv : Fin 12 ⊕ Fin 30 ≃ Fin 42 := finSumFinEquiv
def coneTwoIndexEquiv : Fin 30 ⊕ Fin 20 ≃ Fin 50 := finSumFinEquiv

def coneOneCoordinate (e : Fin 42) : ConeOne →ₗ[ℝ] ℝ where
  toFun x := match coneOneIndexEquiv.symm e with
    | Sum.inl p => x.1 p
    | Sum.inr s => x.2 s
  map_add' x y := by
    generalize h : coneOneIndexEquiv.symm e = z
    cases z <;> simp [h]
  map_smul' a x := by
    generalize h : coneOneIndexEquiv.symm e = z
    cases z <;> simp [h]

def coneTwoCoordinate (f : Fin 50) : ConeTwo →ₗ[ℝ] ℝ where
  toFun x := match coneTwoIndexEquiv.symm f with
    | Sum.inl e => x.1 e
    | Sum.inr t => x.2 t
  map_add' x y := by
    generalize h : coneTwoIndexEquiv.symm f = z
    cases z <;> simp [h]
  map_smul' a x := by
    generalize h : coneTwoIndexEquiv.symm f = z
    cases z <;> simp [h]

@[simp] theorem coneOneCoordinate_inl (p : Fin 12) (x : ConeOne) :
    coneOneCoordinate (coneOneIndexEquiv (Sum.inl p)) x = x.1 p := by
  simp [coneOneCoordinate, coneOneIndexEquiv]

@[simp] theorem coneOneCoordinate_inr (e : Fin 30) (x : ConeOne) :
    coneOneCoordinate (coneOneIndexEquiv (Sum.inr e)) x = x.2 e := by
  simp [coneOneCoordinate, coneOneIndexEquiv]

@[simp] theorem coneTwoCoordinate_inl (e : Fin 30) (x : ConeTwo) :
    coneTwoCoordinate (coneTwoIndexEquiv (Sum.inl e)) x = x.1 e := by
  simp [coneTwoCoordinate, coneTwoIndexEquiv]

@[simp] theorem coneTwoCoordinate_inr (t : Fin 20) (x : ConeTwo) :
    coneTwoCoordinate (coneTwoIndexEquiv (Sum.inr t)) x = x.2 t := by
  simp [coneTwoCoordinate, coneTwoIndexEquiv]

def edgeRestriction (t : Fin 20) : ConeOne →ₗ[ℝ] (Fin 6 → ℝ) where
  toFun x e := signR (edgePositive t e) * coneOneCoordinate (edgeIndex t e) x
  map_add' x y := by
    funext e
    simp [mul_add]
  map_smul' a x := by
    funext e
    simp only [map_smul, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

@[simp] theorem edgeRestriction_apply (t : Fin 20) (x : ConeOne) (e : Fin 6) :
    edgeRestriction t x e =
      signR (edgePositive t e) * coneOneCoordinate (edgeIndex t e) x := rfl

def faceRestriction (t : Fin 20) : ConeTwo →ₗ[ℝ] (Fin 4 → ℝ) where
  toFun x f := signR (facePositive t f) * coneTwoCoordinate (faceIndex t f) x
  map_add' x y := by
    funext f
    simp [mul_add]
  map_smul' a x := by
    funext f
    simp only [map_smul, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

@[simp] theorem faceRestriction_apply (t : Fin 20) (x : ConeTwo) (f : Fin 4) :
    faceRestriction t x f =
      signR (facePositive t f) * coneTwoCoordinate (faceIndex t f) x := rfl

theorem edgeRestriction_joint : JointKernelTrivial edgeRestriction := by
  intro x hx
  apply Prod.ext
  · funext p
    have hcoord : coneOneCoordinate (coneOneIndexEquiv (Sum.inl p)) x = 0 := by
      obtain ⟨⟨t, e⟩, hindex⟩ := edgeIndex_surjective
        (coneOneIndexEquiv (Sum.inl p))
      change edgeIndex t e = coneOneIndexEquiv (Sum.inl p) at hindex
      have hz := congrFun (hx t) e
      change signR (edgePositive t e) * coneOneCoordinate (edgeIndex t e) x = 0 at hz
      rw [hindex] at hz
      exact (mul_eq_zero.mp hz).resolve_left (signR_ne_zero _)
    simpa [coneOneCoordinate] using hcoord
  · funext e
    have hcoord : coneOneCoordinate (coneOneIndexEquiv (Sum.inr e)) x = 0 := by
      obtain ⟨⟨t, l⟩, hindex⟩ := edgeIndex_surjective
        (coneOneIndexEquiv (Sum.inr e))
      change edgeIndex t l = coneOneIndexEquiv (Sum.inr e) at hindex
      have hz := congrFun (hx t) l
      change signR (edgePositive t l) * coneOneCoordinate (edgeIndex t l) x = 0 at hz
      rw [hindex] at hz
      exact (mul_eq_zero.mp hz).resolve_left (signR_ne_zero _)
    simpa [coneOneCoordinate] using hcoord

theorem faceRestriction_joint : JointKernelTrivial faceRestriction := by
  intro x hx
  apply Prod.ext
  · funext e
    have hcoord : coneTwoCoordinate (coneTwoIndexEquiv (Sum.inl e)) x = 0 := by
      obtain ⟨⟨t, f⟩, hindex⟩ := faceIndex_surjective
        (coneTwoIndexEquiv (Sum.inl e))
      change faceIndex t f = coneTwoIndexEquiv (Sum.inl e) at hindex
      have hz := congrFun (hx t) f
      change signR (facePositive t f) * coneTwoCoordinate (faceIndex t f) x = 0 at hz
      rw [hindex] at hz
      exact (mul_eq_zero.mp hz).resolve_left (signR_ne_zero _)
    simpa [coneTwoCoordinate] using hcoord
  · funext t
    have hcoord : coneTwoCoordinate (coneTwoIndexEquiv (Sum.inr t)) x = 0 := by
      obtain ⟨⟨c, f⟩, hindex⟩ := faceIndex_surjective
        (coneTwoIndexEquiv (Sum.inr t))
      change faceIndex c f = coneTwoIndexEquiv (Sum.inr t) at hindex
      have hz := congrFun (hx c) f
      change signR (facePositive c f) * coneTwoCoordinate (faceIndex c f) x = 0 at hz
      rw [hindex] at hz
      exact (mul_eq_zero.mp hz).resolve_left (signR_ne_zero _)
    simpa [coneTwoCoordinate] using hcoord

def sourceVolume (_ : Fin 20) : ℝ := eval volumeQ5
def localEdgeMass (_ : Fin 20) : LinearMap.BilinForm ℝ (Fin 6 → ℝ) :=
  realForm edgeMass
def localFaceMass (_ : Fin 20) : LinearMap.BilinForm ℝ (Fin 4 → ℝ) :=
  realForm faceMass

def sourceMassOne : LinearMap.BilinForm ℝ ConeOne :=
  assemble sourceVolume localEdgeMass edgeRestriction

def sourceMassTwo : LinearMap.BilinForm ℝ ConeTwo :=
  assemble sourceVolume localFaceMass faceRestriction

theorem sourceMassOne_positive : Positive sourceMassOne := by
  exact (assemble_positive_iff sourceVolume localEdgeMass edgeRestriction
    (fun _ => volume_positive) (fun _ => realForm_positive edgeMass_posDef)).mpr
    edgeRestriction_joint

theorem sourceMassTwo_positive : Positive sourceMassTwo := by
  exact (assemble_positive_iff sourceVolume localFaceMass faceRestriction
    (fun _ => volume_positive) (fun _ => realForm_positive faceMass_posDef)).mpr
    faceRestriction_joint

theorem sourceMassOne_symm : sourceMassOne.IsSymm :=
  assemble_symm sourceVolume localEdgeMass edgeRestriction
    (fun _ => realForm_symm edgeMass_posDef)

theorem sourceMassTwo_symm : sourceMassTwo.IsSymm :=
  assemble_symm sourceVolume localFaceMass faceRestriction
    (fun _ => realForm_symm faceMass_posDef)

def localCurlReal : Matrix (Fin 4) (Fin 6) ℝ := localCurl.map eval

#print axioms assemble_positive_iff
#print axioms edgeRestriction_joint
#print axioms faceRestriction_joint
#print axioms sourceMassOne_positive
#print axioms sourceMassTwo_positive

end

end OPH.WhitneySourceAssembly
