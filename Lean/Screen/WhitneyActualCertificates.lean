import WhitneyGeneratedCertificate

set_option autoImplicit false
set_option maxHeartbeats 8000000
set_option maxRecDepth 8192

namespace OPH.WhitneyActualCertificates

open OPH.WhitneyFiniteCertificate
open OPH.WhitneyAlgebraicLDL
def edgeMass : Matrix (Fin 6) (Fin 6) Q5 :=
  OPH.WhitneyGeneratedCertificate.edge_massMatrix
def edgeMassCertificate : Certificate 6 :=
  OPH.WhitneyGeneratedCertificate.edge_massCertificate

def faceMass : Matrix (Fin 4) (Fin 4) Q5 :=
  OPH.WhitneyGeneratedCertificate.face_massMatrix
def faceMassCertificate : Certificate 4 :=
  OPH.WhitneyGeneratedCertificate.face_massCertificate

def stability24 : Matrix (Fin 6) (Fin 6) Q5 :=
  OPH.WhitneyGeneratedCertificate.stability_24Matrix
def stability24Certificate : Certificate 6 :=
  OPH.WhitneyGeneratedCertificate.stability_24Certificate

def stability48 : Matrix (Fin 6) (Fin 6) Q5 :=
  OPH.WhitneyGeneratedCertificate.stability_48Matrix
def stability48Certificate : Certificate 6 :=
  OPH.WhitneyGeneratedCertificate.stability_48Certificate

theorem edgeMass_checked : ldlCheck edgeMass edgeMassCertificate = true := by
  exact OPH.WhitneyGeneratedCertificate.edge_mass_checked

theorem faceMass_checked : ldlCheck faceMass faceMassCertificate = true := by
  exact OPH.WhitneyGeneratedCertificate.face_mass_checked

theorem stability24_checked : ldlCheck stability24 stability24Certificate = true := by
  exact OPH.WhitneyGeneratedCertificate.stability_24_checked

theorem stability48_checked : ldlCheck stability48 stability48Certificate = true := by
  exact OPH.WhitneyGeneratedCertificate.stability_48_checked

theorem edgeMass_posDef : (evalMatrix edgeMass).PosDef :=
  ldl_checked_sound edgeMass_checked

theorem faceMass_posDef : (evalMatrix faceMass).PosDef :=
  ldl_checked_sound faceMass_checked

theorem stability24_posDef : (evalMatrix stability24).PosDef :=
  ldl_checked_sound stability24_checked

theorem stability48_posDef : (evalMatrix stability48).PosDef :=
  ldl_checked_sound stability48_checked

def tamper00 {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    Matrix (Fin n) (Fin n) Q5 :=
  fun i j => if i.1 = 0 ∧ j.1 = 0 then A i j + 1 else A i j

theorem edgeMass_tamper_rejected :
    ldlCheck (tamper00 edgeMass) edgeMassCertificate = false := by
  decide +kernel

theorem stability24_tamper_rejected :
    ldlCheck (tamper00 stability24) stability24Certificate = false := by
  decide +kernel

#print axioms edgeMass_checked
#print axioms faceMass_checked
#print axioms stability24_checked
#print axioms stability48_checked
#print axioms edgeMass_posDef
#print axioms faceMass_posDef
#print axioms stability24_posDef
#print axioms edgeMass_tamper_rejected
#print axioms stability24_tamper_rejected

end OPH.WhitneyActualCertificates
