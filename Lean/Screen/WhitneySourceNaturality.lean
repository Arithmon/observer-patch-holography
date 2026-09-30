import WhitneyCertifiedConsumers
import WhitneyFrameMorphism
import WhitneyExactSourceProblem

set_option autoImplicit false
set_option maxHeartbeats 2000000
set_option maxRecDepth 8192

open scoped BigOperators

namespace OPH.WhitneySourceNaturality

open OPH.ConeCochainBridge
open OPH.WhitneyConeModes
open OPH.WhitneySourceAssembly
open OPH.WhitneyCertifiedConsumers
open OPH.WhitneyConeMass
open OPH.WhitneyQuantumBridge
open OPH.WhitneyFrameMorphism
open OPH.WhitneyFiniteCertificate

noncomputable section

def pullback {V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (e : V ≃ₗ[ℝ] W)
    (B : LinearMap.BilinForm ℝ W) : LinearMap.BilinForm ℝ V :=
  B.compl₁₂ e.toLinearMap e.toLinearMap

@[simp] theorem pullback_apply {V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (e : V ≃ₗ[ℝ] W)
    (B : LinearMap.BilinForm ℝ W) (x y : V) :
    pullback e B x y = B (e x) (e y) := rfl

@[simp] theorem pullback_identity {V : Type*} [AddCommGroup V] [Module ℝ V]
    (B : LinearMap.BilinForm ℝ V) :
    pullback (LinearEquiv.refl ℝ V) B = B := by
  apply LinearMap.ext
  intro x
  apply LinearMap.ext
  intro y
  rfl

theorem pullback_composition {U V W : Type*}
    [AddCommGroup U] [Module ℝ U] [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (e : U ≃ₗ[ℝ] V) (f : V ≃ₗ[ℝ] W)
    (B : LinearMap.BilinForm ℝ W) :
    pullback (e.trans f) B = pullback e (pullback f B) := by
  apply LinearMap.ext
  intro x
  apply LinearMap.ext
  intro y
  rfl

theorem pullback_positive_iff {V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (e : V ≃ₗ[ℝ] W)
    (B : LinearMap.BilinForm ℝ W) :
    Positive (pullback e B) ↔ Positive B := by
  constructor
  · intro h y hy
    let x := e.symm y
    have hx : x ≠ 0 := by
      intro hz
      apply hy
      simpa [x] using congrArg e hz
    simpa [x] using h x hx
  · intro h x hx
    apply h (e x)
    simpa using e.injective.ne hx

/-- Public API for a lawful signed source reexpression.  The edge and face
coordinate maps must move the two mass forms, curl, and the mass-weighted
constraint together.  Signs and permutations are represented by the two
linear equivalences. -/
structure LawfulWhitneyReexpression
    (Edge Face : Type*) [AddCommGroup Edge] [Module ℝ Edge]
    [AddCommGroup Face] [Module ℝ Face] where
  edge : ConeOne ≃ₗ[ℝ] Edge
  face : ConeTwo ≃ₗ[ℝ] Face
  massOne : LinearMap.BilinForm ℝ Edge
  massTwo : LinearMap.BilinForm ℝ Face
  gradient : ConeZero →ₗ[ℝ] Edge
  curl : Edge →ₗ[ℝ] Face
  constraint : Edge →ₗ[ℝ] Module.Dual ℝ ConeZero
  massOne_symm : massOne.IsSymm
  massTwo_symm : massTwo.IsSymm
  massOne_transport : ∀ x y, massOne (edge x) (edge y) = sourceMassOne x y
  massTwo_transport : ∀ x y, massTwo (face x) (face y) = sourceMassTwo x y
  gradient_transport : ∀ x, gradient x = edge (coneGradient x)
  curl_transport : ∀ x, curl (edge x) = face (coneCurl x)
  constraint_transport : ∀ x,
    constraint (edge x) = coneConstraintWithMass sourceMassOne x

namespace LawfulWhitneyReexpression

variable {Edge Face : Type*} [AddCommGroup Edge] [Module ℝ Edge]
  [AddCommGroup Face] [Module ℝ Face]

/-- Any pair of invertible signed coordinate changes induces a lawful source
reexpression.  D/gradient, C/curl, both mass forms, and the constraint are
defined from the same two equivalences, so a caller cannot move only one. -/
def ofEquivs (edge : ConeOne ≃ₗ[ℝ] Edge) (face : ConeTwo ≃ₗ[ℝ] Face) :
    LawfulWhitneyReexpression Edge Face where
  edge := edge
  face := face
  massOne := pullback edge.symm sourceMassOne
  massTwo := pullback face.symm sourceMassTwo
  gradient := edge.toLinearMap.comp coneGradientLinear
  curl := face.toLinearMap.comp (coneCurlLinear.comp edge.symm.toLinearMap)
  constraint := (coneConstraintWithMass sourceMassOne).comp edge.symm.toLinearMap
  massOne_symm := by
    constructor
    intro x y
    exact sourceMassOne_symm.eq (edge.symm x) (edge.symm y)
  massTwo_symm := by
    constructor
    intro x y
    exact sourceMassTwo_symm.eq (face.symm x) (face.symm y)
  massOne_transport _ _ := by simp
  massTwo_transport _ _ := by simp
  gradient_transport _ := rfl
  curl_transport _ := by simp
  constraint_transport _ := by simp

def identity : LawfulWhitneyReexpression ConeOne ConeTwo :=
  ofEquivs (LinearEquiv.refl ℝ ConeOne) (LinearEquiv.refl ℝ ConeTwo)

@[simp] theorem ofEquivs_edge_identity : identity.edge = LinearEquiv.refl ℝ ConeOne := rfl

@[simp] theorem ofEquivs_face_identity : identity.face = LinearEquiv.refl ℝ ConeTwo := rfl

theorem ofEquivs_edge_composition {Edge₂ : Type*} [AddCommGroup Edge₂] [Module ℝ Edge₂]
    (edge₁ : ConeOne ≃ₗ[ℝ] Edge) (edge₂ : Edge ≃ₗ[ℝ] Edge₂)
    (face₁ : ConeTwo ≃ₗ[ℝ] Face) (face₂ : Face ≃ₗ[ℝ] Face) :
    (ofEquivs (edge₁.trans edge₂) (face₁.trans face₂)).edge = edge₁.trans edge₂ := rfl

theorem ofEquivs_face_composition {Face₂ : Type*} [AddCommGroup Face₂] [Module ℝ Face₂]
    (edge₁ : ConeOne ≃ₗ[ℝ] Edge) (edge₂ : Edge ≃ₗ[ℝ] Edge)
    (face₁ : ConeTwo ≃ₗ[ℝ] Face) (face₂ : Face ≃ₗ[ℝ] Face₂) :
    (ofEquivs (edge₁.trans edge₂) (face₁.trans face₂)).face = face₁.trans face₂ := rfl

def stiffness (transport : LawfulWhitneyReexpression Edge Face) :
    LinearMap.BilinForm ℝ Edge :=
  transport.massTwo.compl₁₂ transport.curl transport.curl

theorem stiffness_symm (transport : LawfulWhitneyReexpression Edge Face) :
    transport.stiffness.IsSymm := by
  constructor
  intro x y
  exact transport.massTwo_symm.eq (transport.curl x) (transport.curl y)

def frameIsometry (transport : LawfulWhitneyReexpression Edge Face) :
    FrameIsometry sourceMassOne (coneStiffnessWithMass sourceMassTwo)
      (coneConstraintWithMass sourceMassOne)
      transport.massOne transport.stiffness transport.constraint where
  toLinearEquiv := transport.edge
  mass_preserving := transport.massOne_transport
  stiffness_preserving x y := by
    simp only [stiffness, LinearMap.compl₁₂_apply]
    rw [transport.curl_transport, transport.curl_transport,
      transport.massTwo_transport]
    rfl
  constrained_iff x := by
    rw [transport.constraint_transport]

def frame (transport : LawfulWhitneyReexpression Edge Face) :
    PositiveNormalFrame (ι := Fin 30) transport.massOne transport.stiffness
      transport.constraint :=
  FrameIsometry.mapFrame transport.frameIsometry sourcePositiveNormalFrame

theorem certified_source_to_classical_consumer
    (transport : LawfulWhitneyReexpression Edge Face) (q velocity : Fin 30 → ℝ) :
    (∑ x, (velocity x ^ 2 - sourcePositiveNormalFrame.omega x ^ 2 * q x ^ 2) / 2) =
      ∑ y, ((FrameIsometry.transition transport.frameIsometry
        sourcePositiveNormalFrame transport.frame velocity y) ^ 2 -
        transport.frame.omega y ^ 2 *
          (FrameIsometry.transition transport.frameIsometry
            sourcePositiveNormalFrame transport.frame q y) ^ 2) / 2 :=
  FrameIsometry.classicalAction_naturality transport.frameIsometry
    sourcePositiveNormalFrame transport.frame transport.massOne_symm q velocity

theorem certified_source_to_quantum_consumer
    (transport : LawfulWhitneyReexpression Edge Face) (hbar : ℝ)
    (p : ParticlePolynomial (Fin 30)) :
    quantumHamiltonian hbar transport.frame.omega
        (FrameIsometry.polynomialMap transport.frameIsometry
          sourcePositiveNormalFrame transport.frame p) =
      FrameIsometry.polynomialMap transport.frameIsometry
        sourcePositiveNormalFrame transport.frame
        (quantumHamiltonian hbar sourcePositiveNormalFrame.omega p) :=
  FrameIsometry.quantumHamiltonian_naturality transport.frameIsometry
    sourcePositiveNormalFrame transport.frame transport.massOne_symm
    transport.stiffness_symm hbar p

theorem certified_source_to_gradient_consumer
    (transport : LawfulWhitneyReexpression Edge Face) (u : ConeZero) :
    transport.curl (transport.gradient u) = 0 := by
  rw [transport.gradient_transport, transport.curl_transport]
  simp [coneCurl_gradient]

theorem certified_source_to_dynamics_consumer
    (transport : LawfulWhitneyReexpression Edge Face) (y : Edge) (hy : y ≠ 0) :
    transport.massOne y y - transport.stiffness y y / 48 >
      transport.massOne y y / 2 := by
  let x := transport.edge.symm y
  have hx : x ≠ 0 := by
    intro h
    apply hy
    simpa [x] using congrArg transport.edge h
  have h := source_half_mass_gap_positive x hx
  have hm := transport.massOne_transport x x
  have hs : transport.stiffness (transport.edge x) (transport.edge x) =
      coneStiffnessWithMass sourceMassTwo x x := by
    simp only [stiffness, LinearMap.compl₁₂_apply]
    rw [transport.curl_transport, transport.massTwo_transport]
    rfl
  have he : transport.edge x = y := transport.edge.apply_symm_apply y
  rw [he] at hm hs
  rw [← hm, ← hs] at h
  exact h

end LawfulWhitneyReexpression

/-- The current Lean carrier orders radial faces first, while the Python source
geometry orders the twenty outer faces before the thirty radial faces. -/
abbrev PythonConeTwo := (Fin 20 → ℝ) × (Fin 30 → ℝ)

def faceBlockSwap : ConeTwo ≃ₗ[ℝ] PythonConeTwo where
  toFun x := (x.2, x.1)
  invFun x := (x.2, x.1)
  left_inv _ := rfl
  right_inv _ := rfl
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

@[simp] theorem faceBlockSwap_apply (x : ConeTwo) :
    faceBlockSwap x = (x.2, x.1) := rfl

@[simp] theorem faceBlockSwap_symm_apply (x : PythonConeTwo) :
    faceBlockSwap.symm x = (x.2, x.1) := rfl

def pythonFaceRestriction (t : Fin 20) : PythonConeTwo →ₗ[ℝ] (Fin 4 → ℝ) :=
  (faceRestriction t).comp faceBlockSwap.symm.toLinearMap

def pythonMassTwo : LinearMap.BilinForm ℝ PythonConeTwo :=
  assemble sourceVolume localFaceMass pythonFaceRestriction

theorem pythonMassTwo_eq_pullback :
    pythonMassTwo = pullback faceBlockSwap.symm sourceMassTwo := by
  apply LinearMap.ext
  intro x
  apply LinearMap.ext
  intro y
  simp [pythonMassTwo, sourceMassTwo, pythonFaceRestriction, assemble_apply]

theorem sourceMassTwo_roundTrip :
    pullback faceBlockSwap pythonMassTwo = sourceMassTwo := by
  rw [pythonMassTwo_eq_pullback, ← pullback_composition]
  convert pullback_identity sourceMassTwo using 1

theorem pythonMassTwo_positive : Positive pythonMassTwo := by
  rw [pythonMassTwo_eq_pullback, pullback_positive_iff]
  exact sourceMassTwo_positive

def pythonConeCurl : ConeOne →ₗ[ℝ] PythonConeTwo :=
  faceBlockSwap.toLinearMap.comp coneCurlLinear

@[simp] theorem pythonConeCurl_apply (x : ConeOne) :
    pythonConeCurl x = faceBlockSwap (coneCurl x) := rfl

def pythonStiffness : LinearMap.BilinForm ℝ ConeOne :=
  pythonMassTwo.compl₁₂ pythonConeCurl pythonConeCurl

theorem pythonStiffness_eq_source :
    pythonStiffness = coneStiffnessWithMass sourceMassTwo := by
  apply LinearMap.ext
  intro x
  apply LinearMap.ext
  intro y
  simp only [pythonStiffness, LinearMap.compl₁₂_apply, pythonConeCurl_apply,
    coneStiffnessWithMass, LinearMap.mk₂_apply]
  rw [pythonMassTwo_eq_pullback]
  rfl

theorem pythonStiffness_symm : pythonStiffness.IsSymm := by
  rw [pythonStiffness_eq_source]
  exact coneStiffnessWithMass_symm sourceMassTwo sourceMassTwo_symm

theorem pythonMassTwo_symm : pythonMassTwo.IsSymm := by
  rw [pythonMassTwo_eq_pullback]
  constructor
  intro x y
  exact sourceMassTwo_symm.eq (faceBlockSwap.symm x) (faceBlockSwap.symm y)

def sourcePythonLawfulReexpression :
    LawfulWhitneyReexpression ConeOne PythonConeTwo where
  edge := LinearEquiv.refl ℝ ConeOne
  face := faceBlockSwap
  massOne := sourceMassOne
  massTwo := pythonMassTwo
  gradient := coneGradientLinear
  curl := pythonConeCurl
  constraint := coneConstraintWithMass sourceMassOne
  massOne_symm := sourceMassOne_symm
  massTwo_symm := pythonMassTwo_symm
  massOne_transport _ _ := rfl
  massTwo_transport x y := by
    rw [pythonMassTwo_eq_pullback]
    rfl
  gradient_transport _ := rfl
  curl_transport _ := rfl
  constraint_transport _ := rfl

def runtimeStability24Reexpression :
    OPH.WhitneyCertificatePipeline.LawfulSignedReexpression
      (OPH.WhitneyExactSourceProblem.source.stability
        (OPH.WhitneyExactSourceProblem.q 24 0))
      (OPH.WhitneyExactSourceProblem.runtimeTarget.stability
        (OPH.WhitneyExactSourceProblem.q 24 0)) where
  forward := OPH.WhitneyExactSourceProblem.runtimeEdgeForward.transpose
  backward := OPH.WhitneyExactSourceProblem.runtimeEdgeForward
  forward_backward := by decide +kernel
  backward_forward := by decide +kernel
  target_eq := by
    rw [OPH.WhitneyExactSourceProblem.runtimeTarget_stability24_transport]
    rfl

/-- A checked local stability certificate packaged with two fixed, independently
proved reexpressions: the local signed runtime permutation and the global
Python face-order convention.  These have different carriers and are not one
coordinate policy.  No theorem here identifies the local permutation with the
global face-order map; the global consumer results follow from the separately
certified source masses. -/
structure CertifiedLawfulWhitneyReexpression where
  sourceCertificate : Certificate 6
  sourceCertificate_checked :
    OPH.WhitneyCertificatePipeline.checkLDL
      (OPH.WhitneyExactSourceProblem.source.stability
        (OPH.WhitneyExactSourceProblem.q 24 0)) sourceCertificate = true

namespace CertifiedLawfulWhitneyReexpression

def sourceTransport (_bundle : CertifiedLawfulWhitneyReexpression) :
    LawfulWhitneyReexpression ConeOne PythonConeTwo :=
  sourcePythonLawfulReexpression

def exactSourceTransport (_bundle : CertifiedLawfulWhitneyReexpression) :
    OPH.WhitneyExactSourceProblem.LawfulExactReexpression
      OPH.WhitneyExactSourceProblem.source
      OPH.WhitneyExactSourceProblem.runtimeTarget :=
  OPH.WhitneyExactSourceProblem.runtimeLawfulReexpression

def certificateTransport (_bundle : CertifiedLawfulWhitneyReexpression) :
    OPH.WhitneyCertificatePipeline.LawfulSignedReexpression
      (OPH.WhitneyExactSourceProblem.source.stability
        (OPH.WhitneyExactSourceProblem.q 24 0))
      (OPH.WhitneyExactSourceProblem.runtimeTarget.stability
        (OPH.WhitneyExactSourceProblem.q 24 0)) :=
  runtimeStability24Reexpression

theorem target_certificate_exists
    (bundle : CertifiedLawfulWhitneyReexpression) :
    ∃ certificate, OPH.WhitneyCertificatePipeline.checkLDL
      (OPH.WhitneyExactSourceProblem.runtimeTarget.stability
        (OPH.WhitneyExactSourceProblem.q 24 0)) certificate = true :=
  (OPH.WhitneyCertificatePipeline.certificate_exists_transport_iff
    bundle.certificateTransport).2
      ⟨bundle.sourceCertificate, bundle.sourceCertificate_checked⟩

theorem certified_classical_consumer
    (bundle : CertifiedLawfulWhitneyReexpression) (q velocity : Fin 30 → ℝ) :
    (∑ x, (velocity x ^ 2 - sourcePositiveNormalFrame.omega x ^ 2 * q x ^ 2) / 2) =
      ∑ y, ((FrameIsometry.transition bundle.sourceTransport.frameIsometry
        sourcePositiveNormalFrame bundle.sourceTransport.frame velocity y) ^ 2 -
        bundle.sourceTransport.frame.omega y ^ 2 *
          (FrameIsometry.transition bundle.sourceTransport.frameIsometry
            sourcePositiveNormalFrame bundle.sourceTransport.frame q y) ^ 2) / 2 :=
  bundle.sourceTransport.certified_source_to_classical_consumer q velocity

theorem certified_quantum_consumer
    (bundle : CertifiedLawfulWhitneyReexpression) (hbar : ℝ)
    (p : ParticlePolynomial (Fin 30)) :
    quantumHamiltonian hbar bundle.sourceTransport.frame.omega
        (FrameIsometry.polynomialMap bundle.sourceTransport.frameIsometry
          sourcePositiveNormalFrame bundle.sourceTransport.frame p) =
      FrameIsometry.polynomialMap bundle.sourceTransport.frameIsometry
        sourcePositiveNormalFrame bundle.sourceTransport.frame
        (quantumHamiltonian hbar sourcePositiveNormalFrame.omega p) :=
  bundle.sourceTransport.certified_source_to_quantum_consumer hbar p

theorem certified_gradient_consumer
    (bundle : CertifiedLawfulWhitneyReexpression) (u : ConeZero) :
    bundle.sourceTransport.curl (bundle.sourceTransport.gradient u) = 0 :=
  bundle.sourceTransport.certified_source_to_gradient_consumer u

theorem certified_dynamics_consumer
    (bundle : CertifiedLawfulWhitneyReexpression) (y : ConeOne) (hy : y ≠ 0) :
    bundle.sourceTransport.massOne y y -
        bundle.sourceTransport.stiffness y y / 48 >
      bundle.sourceTransport.massOne y y / 2 :=
  bundle.sourceTransport.certified_source_to_dynamics_consumer y hy

end CertifiedLawfulWhitneyReexpression

def actualCertifiedSourceReexpression :
    CertifiedLawfulWhitneyReexpression where
  sourceCertificate := OPH.WhitneyActualCertificates.stability24Certificate
  sourceCertificate_checked := by decide +kernel

theorem actual_runtime_target_certificate_exists :
    ∃ certificate, OPH.WhitneyCertificatePipeline.checkLDL
      (OPH.WhitneyExactSourceProblem.runtimeTarget.stability
        (OPH.WhitneyExactSourceProblem.q 24 0)) certificate = true :=
  actualCertifiedSourceReexpression.target_certificate_exists

/-- The actual source-coordinate reexpression as a #916 physical-sector
isomorphism.  M₁ and the constraint stay on the edge carrier; M₂ and curl
move together through the nonidentity face block permutation. -/
def sourcePythonFrameIsometry :
    FrameIsometry sourceMassOne (coneStiffnessWithMass sourceMassTwo)
      (coneConstraintWithMass sourceMassOne)
      sourceMassOne pythonStiffness (coneConstraintWithMass sourceMassOne) := by
  simpa [LawfulWhitneyReexpression.stiffness, sourcePythonLawfulReexpression,
    pythonStiffness] using sourcePythonLawfulReexpression.frameIsometry

def pythonSourceFrame :
    PositiveNormalFrame (ι := Fin 30) sourceMassOne pythonStiffness
      (coneConstraintWithMass sourceMassOne) :=
  FrameIsometry.mapFrame sourcePythonFrameIsometry sourcePositiveNormalFrame

theorem source_python_classical_naturality (q velocity : Fin 30 → ℝ) :
    (∑ x, (velocity x ^ 2 - sourcePositiveNormalFrame.omega x ^ 2 * q x ^ 2) / 2) =
      ∑ y, ((FrameIsometry.transition sourcePythonFrameIsometry
        sourcePositiveNormalFrame pythonSourceFrame velocity y) ^ 2 -
        pythonSourceFrame.omega y ^ 2 *
          (FrameIsometry.transition sourcePythonFrameIsometry
            sourcePositiveNormalFrame pythonSourceFrame q y) ^ 2) / 2 :=
  FrameIsometry.classicalAction_naturality sourcePythonFrameIsometry
    sourcePositiveNormalFrame pythonSourceFrame sourceMassOne_symm q velocity

theorem source_python_quantum_naturality (hbar : ℝ)
    (p : ParticlePolynomial (Fin 30)) :
    quantumHamiltonian hbar pythonSourceFrame.omega
        (FrameIsometry.polynomialMap sourcePythonFrameIsometry
          sourcePositiveNormalFrame pythonSourceFrame p) =
      FrameIsometry.polynomialMap sourcePythonFrameIsometry
        sourcePositiveNormalFrame pythonSourceFrame
        (quantumHamiltonian hbar sourcePositiveNormalFrame.omega p) :=
  FrameIsometry.quantumHamiltonian_naturality sourcePythonFrameIsometry
    sourcePositiveNormalFrame pythonSourceFrame sourceMassOne_symm
    pythonStiffness_symm hbar p

#print axioms pullback_composition
#print axioms sourceMassTwo_roundTrip
#print axioms pythonStiffness_eq_source
#print axioms LawfulWhitneyReexpression.certified_source_to_classical_consumer
#print axioms LawfulWhitneyReexpression.certified_source_to_quantum_consumer
#print axioms LawfulWhitneyReexpression.certified_source_to_gradient_consumer
#print axioms LawfulWhitneyReexpression.certified_source_to_dynamics_consumer
#print axioms CertifiedLawfulWhitneyReexpression.target_certificate_exists
#print axioms CertifiedLawfulWhitneyReexpression.certified_classical_consumer
#print axioms CertifiedLawfulWhitneyReexpression.certified_quantum_consumer
#print axioms CertifiedLawfulWhitneyReexpression.certified_dynamics_consumer
#print axioms actual_runtime_target_certificate_exists
#print axioms source_python_classical_naturality
#print axioms source_python_quantum_naturality

end

end OPH.WhitneySourceNaturality
