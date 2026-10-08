import copy
import importlib.util
import json
import shutil
import hashlib
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
RECEIPT=HERE/"runtime"/"k7_gauge_kinetic_cross_control.json"
spec=importlib.util.spec_from_file_location("cross_verifier",HERE/"verify_k7_gauge_kinetic_cross_control_independent.py")
ver=importlib.util.module_from_spec(spec);spec.loader.exec_module(ver)
pspec=importlib.util.spec_from_file_location("cross_producer",HERE/"k7_gauge_kinetic_cross_control.py")
prod=importlib.util.module_from_spec(pspec);pspec.loader.exec_module(prod)

def receipt(): return json.loads(RECEIPT.read_text())
def rejected(mutator):
 d=copy.deepcopy(receipt());mutator(d)
 try: ver.validate(d,ROOT)
 except (ValueError,KeyError,TypeError): return
 raise AssertionError("hostile mutation passed")

def test_clean_packet_passes(): assert ver.validate(receipt(),ROOT)
def test_exact_quadratic_field_operations():
 a=prod.quadratic(1,1); sq=prod.q2mul(a,a)
 assert sq==prod.quadratic(3,2)
 assert prod.q2sub(sq,prod.quadratic(3,2))==prod.quadratic()
 assert prod.q2zero(prod.quadratic()) and not prod.q2zero(sq)
 assert prod.q2scale(a, __import__("fractions").Fraction(3,5))==prod.quadratic(__import__("fractions").Fraction(3,5),__import__("fractions").Fraction(3,5))
def test_mixed_scale_is_control_only():
 d=receipt(); c=d["controls"]["mixed_scale_probe"]
 assert c["classification"]=="NONPHYSICAL_MIXED_SCALE_CONTROL" and not c["eligible_for_framework_cross_verdict"]
def test_b_test_is_control_only():
 d=receipt(); c=d["controls"]["b_test_branch_mismatch_probe"]
 assert c["classification"]=="RGE_BRANCH_MISMATCH_CONTROL" and not c["eligible_for_framework_cross_verdict"]
def test_b_test_scale_conflict_is_preserved():
 assert receipt()["k7_candidates"]["b_test_holonomy_bundle"]["scale"].startswith("CONFLICTED:")
def test_type_three_output_is_recorded_but_not_promoted():
 r=receipt()["k7_candidates"]["rge_bundle"]
 assert r["found"] and r["displayed_outputs"]["alpha_em_inverse_MZ"]=="131.19"
 assert r["exact_vector"].startswith("NOT AVAILABLE") and not r["same_object_contract"]
def test_normalization_covariance_recorded(): assert receipt()["controls"]["normalization_covariance"]["zero_status_preserved"]
def test_no_admissible_vector_is_not_an_off_plane_claim():
 d=receipt(); assert d["comparison"]["admissible_candidate"] is None and not d["comparison"]["same_object_contract"]

# Hostile controls 1–25. Each alters a pinned premise, exact value, type, or
# trust boundary. Together these ensure no cosmetic receipt edit can close it.
def test_hostile_alpha_zero_relabelled_mz(): rejected(lambda d:d["k7_candidates"]["type1_structural_bundle"]["alpha_em_inverse"].update(scale="M_Z"))
def test_hostile_mixed_scale_promotion(): rejected(lambda d:d["controls"]["mixed_scale_probe"].update(eligible_for_framework_cross_verdict=True))
def test_hostile_untransformed_gut_u1(): rejected(lambda d:d["oph_contract"].update(u1_normalization="GUT inserted into Y row"))
def test_hostile_transform_x_only(): rejected(lambda d:d["controls"]["normalization_covariance"].update(all_x_k_b_first_entries_transformed=False))
def test_hostile_transform_k_only(): rejected(lambda d:d["controls"]["normalization_covariance"].update(hostile_partial_row_transform_rejected=False))
def test_hostile_swap_y_and_weak(): rejected(lambda d:d["oph_contract"].update(vector_role="swapped alpha_2 and alpha_Y"))
def test_hostile_alpha_s_instead_inverse(): rejected(lambda d:d["controls"]["mixed_scale_probe"].update(xY=["a","b","sqrt(2)/12"]))
def test_hostile_float_sqrt2(): rejected(lambda d:d["controls"]["mixed_scale_probe"]["exact_D"].update(sqrt2_part="1849.??"))
def test_hostile_beta_mutation(): rejected(lambda d:d["oph_contract"]["beta_column"].__setitem__(0,"7"))
def test_hostile_kinetic_mutation(): rejected(lambda d:d["oph_contract"]["kinetic_column"].__setitem__(0,"4"))
def test_hostile_ng_nh_mutation(): rejected(lambda d:d["oph_contract"].update(nG=2))
def test_hostile_mssm_beta_substitution(): rejected(lambda d:d["oph_contract"].update(beta_column=["33/5","1","-3"]))
def test_hostile_btest_branch_promotion(): rejected(lambda d:d["controls"]["b_test_branch_mismatch_probe"].update(eligible_for_framework_cross_verdict=True))
def test_hostile_public_measurement_input(): rejected(lambda d:d["trust_boundary"].update(public_measurement_read=True))
def test_hostile_sealed_column_read(): rejected(lambda d:d["trust_boundary"].update(oph_sealed_comparison_opened=True))
def test_hostile_mixed_verdict_promotion(): rejected(lambda d:d.update(primary_verdict="K7_GAUGE_VECTOR_EXACTLY_OFF_OPH_MATTER_TRACE_RG_PLANE"))
def test_hostile_exact_near_zero_claim(): rejected(lambda d:d["controls"]["mixed_scale_probe"]["exact_D"].update(rational_part="0",sqrt2_part="0"))
def test_hostile_external_bytes_mutation_and_rehashed_json(tmp_path):
 mini=tmp_path/"repo"
 for rel in ver.OPH_PINS:
  dest=mini/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
 d=receipt()
 for s in d["pins"]["k7_sources"]:
  dest=mini/s["vendored_path"];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/s["vendored_path"],dest)
 target=mini/d["pins"]["k7_sources"][-1]["vendored_path"]
 target.write_bytes(target.read_bytes()+b"\\n-- mutation")
 d["pins"]["k7_sources"][-1]["sha256"]=hashlib.sha256(target.read_bytes()).hexdigest()
 try: ver.validate(d,mini)
 except ValueError: return
 raise AssertionError("mutated source with updated JSON digest passed pinned blob check")
def test_hostile_external_sha_rehash_does_not_change_commit_pin(): rejected(lambda d:d["pins"].update(k7_commit="f"*40))
def test_hostile_erratum_omitted(): rejected(lambda d:d["pins"]["k7_sources"].pop(0))
def test_hostile_lean_arithmetic_does_not_fix_scale(): rejected(lambda d:d["k7_candidates"]["type1_structural_bundle"]["alpha_em_inverse"].update(scale="M_Z"))
def test_hostile_port_branch_identified_with_k7(): rejected(lambda d:d["trust_boundary"].update(oph_port_response_identified_with_k7=True))
def test_hostile_matter_branch_selected(): rejected(lambda d:d["oph_contract"].update(matter_branch_selected_physical=True))
def test_hostile_off_plane_global_falsification(): rejected(lambda d:d["trust_boundary"].update(off_plane_does_not_falsify_oph_globally=False))
