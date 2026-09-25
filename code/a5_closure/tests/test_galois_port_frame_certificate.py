"""Hostile controls for the explicitly frozen two-frame packet."""
from copy import deepcopy
from pathlib import Path
import importlib.util
import sys
import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
import galois_port_frame_certificate as producer
import verify_galois_port_frame_independent as verifier


@pytest.fixture(scope="module")
def packet():
    return producer.build()


def test_independent_replay_and_committed_reference(packet):
    assert verifier.verify(packet)["verified"]["minus"]["degree"]==7
    assert packet==verifier.load(verifier.REFERENCE)


MUTATIONS = [
    "sqrt_sign_fixed", "branch_omitted", "geometric_adjacency", "orientation_reversed",
    "face_removed", "face_duplicated", "chord_relabelled", "wrong_generator_order",
    "laplacian_sign", "selected_gram_source", "unsupported_attachment", "downstream_target",
    "degree_forged", "galois_coordinate_omitted", "rotation_trace_forged", "regular_value_boundary",
]


@pytest.mark.parametrize("mutation",MUTATIONS)
def test_required_mutations_fail_closed(packet,mutation):
    p=deepcopy(packet)
    minus=p["family"]["minus"]
    if mutation=="sqrt_sign_fixed":
        minus["gram"]=deepcopy(p["family"]["plus"]["gram"])
    elif mutation=="branch_omitted":
        del p["family"]["minus"]
    elif mutation=="geometric_adjacency":
        p["distances"]=[[{0:0,1:2,2:1,3:3}[d] for d in row] for row in p["distances"]]
    elif mutation=="orientation_reversed":
        p["faces"]=[[a,c,b] for a,b,c in p["faces"]]
    elif mutation=="face_removed":
        p["faces"].pop()
    elif mutation=="face_duplicated":
        p["faces"][-1]=p["faces"][0]
    elif mutation=="chord_relabelled":
        minus["chord_distance_order"]=[1,2,3]
    elif mutation=="wrong_generator_order":
        p["generator"]=list(range(12))
    elif mutation=="laplacian_sign":
        minus["laplacian_eigenvalue"]=["5","-1"]
    elif mutation=="selected_gram_source":
        p["mathematical_inputs"]=["Lean/Screen/PortFrameGram.lean#g5"]
    elif mutation=="unsupported_attachment":
        p["support_attachment"]={"interface":"PORT-GRAM-SUPPORT-ATTACHMENT","status":"PROVED","source_theorem":"unbound"}
    elif mutation=="downstream_target":
        p["numerical_target"]="27^phi"
    elif mutation=="degree_forged":
        minus["degree_controls"][0]["degree"]=1
    elif mutation=="galois_coordinate_omitted":
        minus["vectors"]=deepcopy(p["family"]["plus"]["vectors"])
    elif mutation=="rotation_trace_forged":
        minus["rotation_trace"]=["1/2","1/2"]
    elif mutation=="regular_value_boundary":
        minus["degree_controls"][0]["direction"]=[0,0,1]
    with pytest.raises(ValueError):
        verifier.verify(p)


def test_producer_dependency_mutation():
    text=(BASE/"galois_port_frame_certificate.py").read_text()
    verifier.source_contract(text)
    with pytest.raises(ValueError,match="import boundary"):
        verifier.source_contract(text+"\nfrom PortFrameGram import g5\n")
    with pytest.raises(ValueError,match="downstream target"):
        verifier.source_contract(text+'\nnumerical_target = "27**phi"\n')


def test_no_producer_import_in_independent_verifier():
    tree=verifier.ast.parse((BASE/"verify_galois_port_frame_independent.py").read_text())
    imports=[node for node in verifier.ast.walk(tree) if isinstance(node,(verifier.ast.Import,verifier.ast.ImportFrom))]
    assert all("galois_port_frame_certificate" not in verifier.ast.unparse(node) for node in imports)


def test_source_response_exchange():
    import galois_source_response_control as control
    result=control.verify()
    assert result["bounded_nonselection"]
    assert result["plus"]==result["minus"]
    assert not result["literal_block_swap_equals_galois"]


def test_source_control_rejects_unchanged_quadrupole():
    import galois_source_response_control as control
    p=control.c.strict_load(control.MODEL/"response.json")
    changed=control.conjugate_serialized(p)
    for g,old in zip(changed["generators"],p["generators"]):
        for i in range(3):
            for j in range(3):
                g[0][i][j][2:]=old[0][i][j][2:]
    with pytest.raises(ValueError,match="quadrupole"):
        control.replay(changed,True)
