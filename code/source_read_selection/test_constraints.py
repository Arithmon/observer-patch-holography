from copy import deepcopy
from fractions import Fraction as F
import ast

import pytest

from . import codec,constraints,verify,verify_constraints
from .test_selection import cli


@pytest.fixture(scope="module")
def packet():
    return codec.load_artifact(codec.HERE/"controls.json")["projection_controls"]


def test_exact_projection_and_same_reference_cover(packet):
    codec.equal(constraints.build(),packet,"projection producer drift")
    result = verify_constraints.verify(packet)
    assert result["correlated_failure_mass"]=="7/16"
    assert result["full_selected_failure_mass"]=="1/6"
    assert result["coarse_selected_failure_mass"]=="0"
    assert result["nonconvex_log_score_ratio"]=="343/256"


def test_projection_verifier_is_independent(packet,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("called constrained producer")
    monkeypatch.setattr(constraints,"build",forbidden)
    monkeypatch.setattr(constraints,"group_projection",forbidden)
    verify_constraints.verify(packet)
    tree = ast.parse((codec.HERE/"verify_constraints.py").read_text())
    assert all("constraints" not in ast.unparse(node) for node in ast.walk(tree)
               if isinstance(node,(ast.Import,ast.ImportFrom)))


MUTATIONS = [
    lambda p:p.update(scope="A3_forces_M1"),
    lambda p:p.update(extra=True),
    lambda p:p.pop("successful_face"),
    lambda p:p["correlated_words"].update(selected=["0","3/4","1/8","1/8"]),
    lambda p:p["correlated_words"].update(selected=["3/8","3/8","1/8","1/8"]),
    lambda p:p["correlated_words"].update(reference=["1/4"]*4),
    lambda p:p["correlated_words"].update(groups=[0,1,0,1]),
    lambda p:p["correlated_words"].update(targets=["1/2","1/2"]),
    lambda p:p["correlated_words"].update(ratios=["1","1"]),
    lambda p:p["correlated_words"].update(failure_mass="0"),
    lambda p:p["correlated_words"].update(cut_avoidance_mass="0"),
    lambda p:p["correlated_words"].update(product_determinant="0"),
    lambda p:p["correlated_words"]["word_order"].reverse(),
    lambda p:p["correlated_words"]["receiver_coefficients"].__setitem__(2,"1/4"),
    lambda p:p["same_reference_cover_comparison"].update(middle_mass="1/4"),
    lambda p:p["same_reference_cover_comparison"].update(coarse_selected=["1/4","1/2","1/4"]),
    lambda p:p["same_reference_cover_comparison"].update(coarse_reference=["3/4","1/4"]),
    lambda p:p["same_reference_cover_comparison"].update(coarse_matrix=[[1,0,0],[0,1,1]]),
    lambda p:p["same_reference_cover_comparison"].update(full_matrix=[[1,1,1]]),
    lambda p:p["same_reference_cover_comparison"].update(constraint_and_coarse_determinant=0),
    lambda p:p["same_reference_cover_comparison"].update(failure_atom=1),
    lambda p:p["same_reference_cover_comparison"]["coarse_failure_reconstruction"].update(offset="0"),
    lambda p:p["same_reference_cover_comparison"]["coarse_failure_reconstruction"].update(coefficients=["0","1"]),
    lambda p:p["same_reference_cover_comparison"]["full_projection"].update(selected=["0","1/2","1/2"]),
    lambda p:p["same_reference_cover_comparison"]["full_projection"].update(extra=True),
    lambda p:p["nonconvex_control"].update(convex=True),
    lambda p:p["nonconvex_control"].update(competitor_score_ratio="1/2"),
    lambda p:p["nonconvex_control"].update(selected=["1/4"]*4),
    lambda p:p["nonconvex_control"]["feasible"].pop(),
    lambda p:p["successful_face"].update(status="source_derived"),
    lambda p:p["successful_face"].update(excluded_atoms=[]),
]


@pytest.mark.parametrize("mutation",MUTATIONS)
def test_constrained_resealed_forgery_fails(packet,mutation,tmp_path):
    changed = deepcopy(packet)
    mutation(changed)
    with pytest.raises((ValueError,KeyError,TypeError)):
        verify_constraints.verify(changed)
    artifact = codec.load_artifact(codec.HERE/"controls.json")
    artifact["projection_controls"] = changed
    controls = tmp_path/"controls.json"
    controls.write_bytes(codec.canonical(artifact))
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt["controls_sha256"] = codec.digest(artifact)
    path = tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(receipt))
    result = cli(controls,path)
    assert result.returncode!=0 and "verification failed" in result.stderr


@pytest.mark.parametrize("bad",[None,True,1,1.0,"1/0","+1","0/1","2/4","NaN",""])
def test_invalid_projection_rational_rejected(bad,packet):
    changed = deepcopy(packet["correlated_words"])
    changed["reference"][0] = bad
    with pytest.raises(ValueError):
        verify_constraints.certify_projection(changed)


@pytest.mark.parametrize("field,value",[
    ("reference",["0","1/2","1/4","1/4"]),
    ("reference",["-1/8","5/8","1/4","1/4"]),
    ("selected",["3/16","9/16","1/8","1/4"]),
    ("groups",[False,0,1,1]),
    ("groups",[0,0,0,0]),
    ("ratios",["3/2"]),
    ("targets",["1","0"]),
])
def test_projection_hypotheses_not_silently_relaxed(packet,field,value):
    changed = deepcopy(packet["correlated_words"])
    changed[field] = value
    with pytest.raises(ValueError):
        verify_constraints.certify_projection(changed)


def test_feasible_but_nonoptimal_candidate_is_rejected(packet):
    candidate = deepcopy(packet["correlated_words"])
    candidate["selected"] = ["1/4","1/2","1/8","1/8"]
    # Correct normalization and both moment constraints; forged dual ratio.
    candidate["ratios"] = ["2","1/2"]
    with pytest.raises(ValueError,match="KL moment certificate"):
        verify_constraints.certify_projection(candidate)


def test_coarse_cover_injectivity_is_not_history_support_completeness(packet):
    cover = packet["same_reference_cover_comparison"]
    coarse = cover["coarse_matrix"]
    selected = list(map(F,cover["coarse_selected"]))
    witness = list(map(F,cover["feasible_positive_witness"]))
    assert verify_constraints.determinant([[0,1,0],*coarse])==-1
    assert all(v>0 for v in verify_constraints.image(coarse,selected))
    assert selected[0]==0 and witness[0]>0
    # Equal columns at history atoms 0 and 1 rule out a linear positive
    # representation of the history-0 indicator by this coarse cover.
    assert [row[0] for row in coarse]==[row[1] for row in coarse]
    assert [1,0,0][0] != [1,0,0][1]
