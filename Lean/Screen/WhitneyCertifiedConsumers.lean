import WhitneyCurlRow0
import WhitneyCurlRow1
import WhitneyCurlRow2
import WhitneyCurlRow3
import WhitneyConeMassNaturality
import WhitneyMaxwellDynamics

set_option autoImplicit false
set_option maxHeartbeats 2000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix

namespace OPH.WhitneyCertifiedConsumers

open OPH.ConeCochainBridge
open OPH.WhitneyFiniteCertificate
open OPH.WhitneyActualCertificates
open OPH.WhitneySourceGeometry
open OPH.WhitneySourceAssembly
open OPH.WhitneyConeMass
open OPH.WhitneyConeMassNaturality
open OPH.WhitneyQuantumBridge
open OPH.WhitneyFrameNaturality

noncomputable section

theorem faceRestriction_coneCurl (t : Fin 20) (x : ConeOne) :
    faceRestriction t (coneCurl x) =
      localCurlReal.mulVec (edgeRestriction t x) := by
  funext f
  fin_cases f
  · exact faceRestriction_coneCurl_0 t x
  · exact faceRestriction_coneCurl_1 t x
  · exact faceRestriction_coneCurl_2 t x
  · exact faceRestriction_coneCurl_3 t x

def localStiffness : LinearMap.BilinForm ℝ (Fin 6 → ℝ) :=
  realForm derivedStiffness

def sourceStiffness : LinearMap.BilinForm ℝ ConeOne :=
  assemble sourceVolume (fun _ => localStiffness) edgeRestriction

theorem eval_derivedStiffness :
    evalMatrix derivedStiffness =
      localCurlReal.transpose * evalMatrix faceMass * localCurlReal := by
  simp [derivedStiffness, localCurlReal, evalMatrix, Matrix.map_mul,
    Matrix.transpose_map]

theorem localStiffness_apply (x y : Fin 6 → ℝ) :
    localStiffness x y =
      realForm faceMass (localCurlReal.mulVec x) (localCurlReal.mulVec y) := by
  rw [localStiffness, realForm, realForm, eval_derivedStiffness,
    Matrix.toBilin'_apply', Matrix.toBilin'_apply']
  calc
    x ⬝ᵥ (localCurlReal.transpose * evalMatrix faceMass * localCurlReal) *ᵥ y =
        x ⬝ᵥ localCurlReal.transpose *ᵥ
          (evalMatrix faceMass *ᵥ (localCurlReal *ᵥ y)) := by
            rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
    _ = (localCurlReal *ᵥ x) ⬝ᵥ
          (evalMatrix faceMass *ᵥ (localCurlReal *ᵥ y)) := by
            rw [Matrix.dotProduct_mulVec, Matrix.vecMul_transpose]

theorem sourceStiffness_eq_coneStiffness :
    sourceStiffness = coneStiffnessWithMass sourceMassTwo := by
  apply LinearMap.ext
  intro x
  apply LinearMap.ext
  intro y
  change sourceStiffness x y = sourceMassTwo (coneCurl x) (coneCurl y)
  rw [sourceStiffness, sourceMassTwo, assemble_apply, assemble_apply]
  apply Finset.sum_congr rfl
  intro t _
  rw [localStiffness_apply, ← faceRestriction_coneCurl,
    ← faceRestriction_coneCurl]
  rfl

def localGap24 : LinearMap.BilinForm ℝ (Fin 6 → ℝ) := realForm stability24

def sourceGap24 : LinearMap.BilinForm ℝ ConeOne :=
  assemble sourceVolume (fun _ => localGap24) edgeRestriction

theorem localGap24_identity (x : Fin 6 → ℝ) :
    localGap24 x x = 24 * realForm edgeMass x x - localStiffness x x := by
  rw [localGap24, localStiffness, realForm, realForm, realForm,
    Matrix.toBilin'_apply, Matrix.toBilin'_apply, Matrix.toBilin'_apply,
    stability24_source_identity]
  simp only [evalMatrix, Matrix.map_apply, map_sub, map_mul]
  simp_rw [mul_sub, sub_mul, Finset.sum_sub_distrib]
  rw [show eval (24 : Q5) = (24 : ℝ) by simpa using map_natCast eval 24]
  rw [Finset.mul_sum]
  simp_rw [Finset.mul_sum]
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  ring

theorem localStiffness_le (x : Fin 6 → ℝ) :
    localStiffness x x ≤ 24 * realForm edgeMass x x := by
  have hnonneg := (realForm_positive stability24_posDef).nonneg x
  rw [← localGap24, localGap24_identity] at hnonneg
  linarith

theorem sourceGap24_positive : Positive sourceGap24 := by
  exact (assemble_positive_iff sourceVolume (fun _ => localGap24) edgeRestriction
    (fun _ => volume_positive) (fun _ => realForm_positive stability24_posDef)).mpr
    edgeRestriction_joint

theorem sourceGap24_identity (x : ConeOne) :
    sourceGap24 x x = 24 * sourceMassOne x x - sourceStiffness x x := by
  rw [sourceGap24, sourceMassOne, sourceStiffness, assemble_apply, assemble_apply,
    assemble_apply]
  rw [Finset.mul_sum]
  simp_rw [localGap24_identity]
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro t _
  simp only [localEdgeMass]
  ring

theorem source_stiffness_bound (x : ConeOne) :
    sourceStiffness x x ≤ 24 * sourceMassOne x x := by
  rw [sourceStiffness, sourceMassOne, assemble_apply, assemble_apply]
  calc
    (∑ t, sourceVolume t * localStiffness (edgeRestriction t x)
      (edgeRestriction t x)) ≤
      ∑ t, sourceVolume t * (24 * realForm edgeMass (edgeRestriction t x)
        (edgeRestriction t x)) := by
          apply Finset.sum_le_sum
          intro t _
          exact mul_le_mul_of_nonneg_left (localStiffness_le (edgeRestriction t x))
            (le_of_lt volume_positive)
    _ = 24 * ∑ t, sourceVolume t * realForm edgeMass (edgeRestriction t x)
        (edgeRestriction t x) := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro t _
          ring

theorem source_cone_stiffness_bound (x : ConeOne) :
    coneStiffnessWithMass sourceMassTwo x x ≤ 24 * sourceMassOne x x := by
  rw [← sourceStiffness_eq_coneStiffness]
  exact source_stiffness_bound x

theorem source_stiffness_bound_strict (x : ConeOne) (hx : x ≠ 0) :
    sourceStiffness x x < 24 * sourceMassOne x x := by
  have hgap := sourceGap24_positive x hx
  rw [sourceGap24_identity] at hgap
  linarith

theorem source_cone_stiffness_bound_strict (x : ConeOne) (hx : x ≠ 0) :
    coneStiffnessWithMass sourceMassTwo x x < 24 * sourceMassOne x x := by
  rw [← sourceStiffness_eq_coneStiffness]
  exact source_stiffness_bound_strict x hx

theorem source_half_mass_gap_positive (x : ConeOne) (hx : x ≠ 0) :
    sourceMassOne x x -
        coneStiffnessWithMass sourceMassTwo x x / 48 > sourceMassOne x x / 2 := by
  have hstrict := source_cone_stiffness_bound_strict x hx
  linarith

def sourcePositiveNormalFrame :
    PositiveNormalFrame (ι := Fin 30) sourceMassOne
      (coneStiffnessWithMass sourceMassTwo) (coneConstraintWithMass sourceMassOne) :=
  conePositiveNormalFrameWithMass sourceMassOne sourceMassTwo
    sourceMassOne_symm sourceMassTwo_symm sourceMassOne_positive sourceMassTwo_positive

theorem source_same_action_normal_modes (q velocity : Fin 30 → ℝ) :
    (sourceMassOne (reconstruct sourcePositiveNormalFrame.vector velocity)
        (reconstruct sourcePositiveNormalFrame.vector velocity) -
      coneStiffnessWithMass sourceMassTwo
        (reconstruct sourcePositiveNormalFrame.vector q)
        (reconstruct sourcePositiveNormalFrame.vector q)) / 2 =
      ∑ x, (velocity x ^ 2 - sourcePositiveNormalFrame.omega x ^ 2 * q x ^ 2) / 2 := by
  simpa [sourcePositiveNormalFrame] using
    cone_same_action_with_mass sourceMassOne sourceMassTwo sourceMassOne_symm
      sourceMassTwo_symm sourceMassOne_positive sourceMassTwo_positive q velocity

theorem source_same_modes_quantized_energy (hbar : ℝ) (hh : 0 < hbar)
    (p : ParticlePolynomial (Fin 30)) :
    (∑ x, (1 / 2 : ℂ) •
      (momentumOperator x (sourcePositiveNormalFrame.omega x)
          (positionScale hbar (sourcePositiveNormalFrame.omega x))
          (momentumOperator x (sourcePositiveNormalFrame.omega x)
            (positionScale hbar (sourcePositiveNormalFrame.omega x)) p) +
        ((sourcePositiveNormalFrame.omega x ^ 2 : ℝ) : ℂ) •
          positionOperator x (positionScale hbar (sourcePositiveNormalFrame.omega x))
            (positionOperator x (positionScale hbar (sourcePositiveNormalFrame.omega x)) p))) =
      quantumHamiltonian hbar sourcePositiveNormalFrame.omega p := by
  simpa [sourcePositiveNormalFrame] using
    cone_same_quantized_energy_with_mass sourceMassOne sourceMassTwo sourceMassOne_symm
      sourceMassTwo_symm sourceMassOne_positive sourceMassTwo_positive hbar hh p

theorem source_corrected_kinetic_positive (v : ConeOne) (hv : v ≠ 0) :
    0 < sourceMassOne v v - (1 / 2 : ℝ) ^ 2 / 12 *
      coneStiffnessWithMass sourceMassTwo v v := by
  exact OPH.WhitneyMaxwellDynamics.corrected_kinetic_positive
    (1 / 2 : ℝ) 24 sourceMassOne (coneStiffnessWithMass sourceMassTwo)
    sourceMassOne_positive source_cone_stiffness_bound (by norm_num) v hv

theorem source_sign_quantum_naturality (hbar : ℝ)
    (p : ParticlePolynomial (Fin 30)) :
    quantumHamiltonian hbar (signFlipFrame sourcePositiveNormalFrame).omega
        (polynomialTransport sourcePositiveNormalFrame
          (signFlipFrame sourcePositiveNormalFrame) p) =
      polynomialTransport sourcePositiveNormalFrame
        (signFlipFrame sourcePositiveNormalFrame)
        (quantumHamiltonian hbar sourcePositiveNormalFrame.omega p) := by
  simpa [sourcePositiveNormalFrame] using
    actual_cone_sign_quantum_naturality sourceMassOne sourceMassTwo
      sourceMassOne_symm sourceMassTwo_symm sourceMassOne_positive sourceMassTwo_positive
      hbar p

#print axioms sourceStiffness_eq_coneStiffness
#print axioms source_stiffness_bound
#print axioms source_stiffness_bound_strict
#print axioms source_half_mass_gap_positive
#print axioms source_same_action_normal_modes
#print axioms source_same_modes_quantized_energy
#print axioms source_corrected_kinetic_positive
#print axioms source_sign_quantum_naturality

end

end OPH.WhitneyCertifiedConsumers
