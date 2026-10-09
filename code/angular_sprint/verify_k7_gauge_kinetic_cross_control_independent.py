#!/usr/bin/env python3
"""Independent custody, typing and exact-arithmetic verifier for the K7/OPH packet."""
from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path

SCHEMA="arithmon.oph.k7_gauge_kinetic_cross.v1"
OPH="0660c94573c8f79c1b86955fbd1d7d927a0e4fa2"
K7="210480b2ba8a13c98e6a6f146fbfd965e4904c00"
K7L="2553ed170d04db8777f84fb2f9ef2b66b386f7fb"
EXPECTED=[
 ("Arithmon/K7",K7,"publications/ERRATUM_v3.5.md","2bcd161a6684f93f83e3541cbd0682738fa21b70","ERRATUM_v3.5.md"),
 ("Arithmon/K7",K7,"publications/papers/markdown/k7_framework_3_5_main.md","3a0895fdaa9624eaed9d083e45566e75b3b1036b","k7_framework_3_5_main.md"),
 ("Arithmon/K7",K7,"docs/openwave-candidate/honest_ledger.md","39efc4eace29269f99a37fc94662fd757574efeb","honest_ledger.md"),
 ("Arithmon/K7-Lean",K7L,"GIFT/Relations/GaugeSector.lean","664e6e620de8f11d628a2c0beb3c5c90d6c190f2","GaugeSector.lean")]
OPH_PINS={"Lean/Screen/KineticFormDichotomy.lean":"671578153bb70bf590e380fb7fd035062a341155","Lean/Screen/RGRepresentationFrontier.lean":"c35eb3bc1f8e5f8bd4633bdea07911049bb985d0","code/angular_sprint/kinetic_form_selection_certificate.py":"9cadce4736c666dc48a242ccfb46011d00b24189","code/angular_sprint/runtime/kinetic_form_selection_receipt.json":"7469c729fd54cc7f558830c03fbae1bf510e5f32"}

def blob(raw): return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def require(ok,msg):
 if not ok: raise ValueError(msg)
def determinant(a,b,c):
 return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def validate(doc, root):
 require(doc.get("schema")==SCHEMA,"schema mismatch")
 pins=doc["pins"]
 require(pins["oph_commit"]==OPH,"OPH commit pin changed")
 require(pins["k7_commit"]==K7 and pins["k7_fetched_main_head"]==K7 and pins["k7_lean_commit"]==K7L and pins["k7_lean_fetched_main_head"]==K7L,"external K7 commit/head pin changed")
 require(pins["oph_sources"]==OPH_PINS,"OPH path/blob provenance changed")
 for path,expected in OPH_PINS.items():
  require(blob((root/path).read_bytes())==expected,f"OPH frozen bytes changed: {path}")
 sources=pins["k7_sources"]
 require(len(sources)==len(EXPECTED),"K7 provenance source omitted/added")
 for s,(repo,commit,path,gitblob,name) in zip(sources,EXPECTED):
  require((s["repository"],s["commit"],s["path"],s["git_blob"],s["vendored_path"])==(repo,commit,path,gitblob,f"code/angular_sprint/external/k7_gauge/{name}"),f"K7 source pin changed: {path}")
  raw=(root/s["vendored_path"]).read_bytes()
  require(hashlib.sha256(raw).hexdigest()==s["sha256"] and len(raw)==s["bytes"],f"K7 vendored bytes changed: {path}")
  require(blob(raw)==gitblob,f"K7 git blob mismatch: {path}")

 # Re-derive OPH columns and cofactors from the declared family formula.
 ng,nh=3,1
 k=(F(10,3),F(2),F(2))
 b=(F(20,9)*ng+F(1,6)*nh,-F(22,3)+F(4,3)*ng+F(1,6)*nh,-F(11)+F(4,3)*ng)
 cof=(k[1]*b[2]-k[2]*b[1],k[2]*b[0]-k[0]*b[2],k[0]*b[1]-k[1]*b[0])
 oc=doc["oph_contract"]
 require(tuple(map(F,oc["kinetic_column"]))==k and tuple(map(F,oc["beta_column"]))==b,"OPH k/b changed")
 require((ng,nh)==(oc["nG"],oc["nH"]),"OPH matter completion changed")
 require(tuple(map(F,oc["cofactors"]))==cof and tuple(map(int,oc["integer_zero_locus"]))==(69,-333,218),"OPH determinant derivation changed")
 require(oc["vector_role"]=="alpha_inverse=(alpha_Y^-1,alpha_2^-1,alpha_3^-1); issue-639 column sealed and not read","OPH vector role changed")
 require(oc["u1_normalization"]=="declared hypercharge Y" and oc["matter_branch_selected_physical"] is False,"OPH typing/selection boundary changed")

 candidates=doc["k7_candidates"]
 t1=candidates["type1_structural_bundle"]
 require(t1["alpha_em_inverse"]["exact"]=="267489/1952" and t1["alpha_em_inverse"]["scale"]=="alpha_em^-1(0)","low-energy alpha relation was relabeled")
 require(t1["alpha_em_inverse"]["public_data_used_in_construction"] is False,"public value consumed as input")
 require(t1["sin2_theta_W"]["exact"]=="3/13" and t1["alpha_s"]["exact"]=="sqrt(2)/12","Type-I source values drift")
 require(t1["assembled_vector"]["classification"]=="MIXED_SCALE_VECTOR" and t1["assembled_vector"]["same_object_contract"] is False,"mixed Type-I vector misclassified")
 require(t1["assembled_vector"]["values"]=={"alpha_Y_inverse":"(10/13)*(267489/1952)","alpha_2_inverse":"(3/13)*(267489/1952)","alpha_3_inverse":"6*sqrt(2)"},"mixed vector component order/value changed")
 bt=candidates["b_test_holonomy_bundle"]
 require(bt["inputs"]["alpha_em_inverse"]=="91*sqrt(2)" and bt["inputs"]["sin2_theta_W"]=="3/13" and bt["inputs"]["alpha_s"]=="sqrt(2)/12","B-test inputs mismatch")
 deriv=bt["normalization_derivation"]
 require(F(91,1)/(F(7)+F(70,3))==3 and deriv["common_factor"]=="alpha_em_inverse/(7+70/3)=3*sqrt(2)","B-test common factor not reconstructed")
 require(deriv["GUT_ratio"]=="14:7:2" and deriv["alpha_em_inverse_relation"]=="alpha_2_inverse + alpha_Y_inverse" and deriv["alpha_Y_inverse_from_alpha_1_GUT_inverse"]=="(5/3)*alpha_1_GUT_inverse","B-test normalization derivation drift")
 require(deriv["GUT_inverse_vector"]==["42*sqrt(2)","21*sqrt(2)","6*sqrt(2)"] and deriv["Y_inverse_vector"]==["70*sqrt(2)","21*sqrt(2)","6*sqrt(2)"],"B-test converted vector mismatch")
 require(bt["values"]["alpha_Y_inverse"]=="70*sqrt(2)" and bt["scale"].startswith("CONFLICTED:") and bt["field_content"].startswith("MSSM"),"B-test vector typing/scale conflict drift")
 require(bt["same_object_contract"] is False and bt["scheme"]=="UNKNOWN","B-test scheme/contract was promoted")
 rg=candidates["rge_bundle"]
 require(rg["found"] is True and rg["classification"]=="CONDITIONAL_NUMERICAL_RGE_OUTPUT__NO_EXACT_SOURCE_CLEAN_REPRODUCER","paper RGE candidate omitted/mistyped")
 require(rg["evidence_scope"]==["vendored K7 publications/ERRATUM_v3.5.md","vendored K7 publications/papers/markdown/k7_framework_3_5_main.md","vendored K7 docs/openwave-candidate/honest_ledger.md","vendored K7-Lean GIFT/Relations/GaugeSector.lean"],"RGE evidence scope broadened beyond pinned surfaces")
 require("pinned K7 source surfaces consumed by this packet" in rg["record"] and rg["reproducible_artifact"]=="NOT ESTABLISHED BY PINNED SOURCE SURFACES","unsupported repository-wide absence claim")
 require(rg["displayed_outputs"]["alpha_em_inverse_MZ"]=="131.19" and rg["displayed_outputs"]["sin2_theta_W_MZ"]=="0.2377" and rg["displayed_outputs"]["alpha_s_split_spectrum"]=="0.1224","paper RGE outputs drift")
 require(rg["scale"]=="M_Z per §5.3 output table" and rg["same_object_contract"] is False and rg["exact_vector"].startswith("NOT AVAILABLE"),"RGE output promoted to an exact OPH vector")

 controls=doc["controls"]
 m=controls["mixed_scale_probe"]
 require(m["classification"]=="NONPHYSICAL_MIXED_SCALE_CONTROL" and m["eligible_for_framework_cross_verdict"] is False,"mixed-scale probe promoted")
 require(m["xY"]==["(10/13)*(267489/1952)","(3/13)*(267489/1952)","6*sqrt(2)"],"mixed vector component/order changed")
 a0=F(267489,1952)
 mixed_r=F(69)*F(10,13)*a0-F(333)*F(3,13)*a0
 mixed_s=F(218)*6
 require(mixed_r==F(-82654101,25376) and mixed_s==F(1308),"mixed probe not independently reconstructed")
 require(m["exact_D"]=={"rational_part":str(mixed_r),"sqrt2_part":str(mixed_s)},"mixed control exact value changed")
 require(mixed_r!=0 or mixed_s!=0,"near-zero treated as zero")
 bb=controls["b_test_branch_mismatch_probe"]
 require(bb["classification"]=="RGE_BRANCH_MISMATCH_CONTROL" and bb["eligible_for_framework_cross_verdict"] is False,"B-test mismatch promoted")
 require(bb["exact_D"]=={"rational_part":"0","sqrt2_part":"-855"},"B-test control exact value changed")
 require((F(69*70-333*21+218*6))==F(-855),"B-test determinant not independently reconstructed")
 # A complete first-row GUT rescaling multiplies the determinant by 3/5.
 x=(F(70),F(21),F(6)); kg=(F(2),F(2),F(2)); bg=(F(41,6),F(-19,6),F(-7))
 # x and the corresponding k,b components all rescale on row 1.
 yg=(F(3,5)*x[0],x[1],x[2]); kyg=(F(3,5)*kg[0],kg[1],kg[2]); byg=(F(3,5)*bg[0],bg[1],bg[2])
 require(determinant(yg,kyg,byg)==F(3,5)*determinant(x,kg,bg),"row covariance does not hold")
 cov=controls["normalization_covariance"]
 require(cov["zero_status_preserved"] is True and cov["coherent_first_row_factor"]=="3/5","normalization covariance missing")
 require(cov["all_x_k_b_first_entries_transformed"] is True and cov["hostile_partial_row_transform_rejected"] is True,"partial normalization flags drift")
 xzero=(F(333),F(69),F(0)); kz=(F(10,3),F(2),F(2)); bz=(F(41,6),F(-19,6),F(-7))
 ydet=determinant(xzero,kz,bz)
 gutdet=determinant((F(3,5)*xzero[0],xzero[1],xzero[2]),(F(3,5)*kz[0],kz[1],kz[2]),(F(3,5)*bz[0],bz[1],bz[2]))
 require(ydet==0 and gutdet==F(3,5)*ydet==0,"Y/GUT exact zero covariance failed")
 partials=cov["partial_row_transform_determinants"]
 only_x=determinant((F(3,5)*xzero[0],xzero[1],xzero[2]),kz,bz)
 only_k=determinant(xzero,(F(3,5)*kz[0],kz[1],kz[2]),bz)
 only_b=determinant(xzero,kz,(F(3,5)*bz[0],bz[1],bz[2]))
 require((only_x,only_k,only_b)==(F(5106,5),F(-644),F(-1886,5)),"hostile partial-row arithmetic changed")
 require(tuple(map(F,(partials["x_only"],partials["k_only"],partials["b_only"])))==(only_x,only_k,only_b) and all(v!=0 for v in (only_x,only_k,only_b)),"partial row conversion incorrectly preserves zero")
 require(doc["comparison"]["same_object_contract"] is False and doc["comparison"]["admissible_candidate"] is None,"no-admissible-vector gate changed")
 require(doc["primary_verdict"]=="NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES","primary verdict drift")
 trust=doc["trust_boundary"]
 for field in ("public_measurement_consumed_as_cross_input","public_measurement_values_parsed_by_cross_control","oph_sealed_comparison_opened","k7_modified","oph_matter_branch_selected","physical_cross_framework_identity_claimed"):
  require(trust[field] is False,f"trust boundary violated: {field}")
 require(trust["vendored_main_paper_contains_measurement_comparison_text"] is True,"measurement-text disclosure missing")
 require(trust["off_plane_does_not_falsify_oph_globally"] is True and trust["oph_port_response_identified_with_k7"] is False and trust["matter_trace_physically_selected"] is False,"OPH branch firewall missing")
 return True

def main():
 p=argparse.ArgumentParser();p.add_argument("--receipt",required=True);a=p.parse_args()
 path=Path(a.receipt).resolve(); root=path.parents[3]
 doc=json.loads(path.read_text());validate(doc,root);print("PASS independent K7/OPH custody and exact contract audit")
if __name__=="__main__":main()
