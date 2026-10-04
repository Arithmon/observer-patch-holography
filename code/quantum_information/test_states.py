"""Independent identities, domain failures and retained audit counterexamples."""

import itertools
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import logm

from quantum_information import (
    conditional_mutual_information, density_matrix, direct_sum_state,
    faithful_log, mutual_information, partial_trace, probabilities,
    relative_entropy, shannon_entropy, von_neumann_entropy,
)


def test_noncommuting_entropy_matches_independent_matrix_log():
    rho = np.array([[.7, .1j], [-.1j, .3]])
    tau = np.array([[.4, .12], [.12, .6]])
    assert np.linalg.norm(rho@tau-tau@rho) > .01
    expected = np.trace(rho@(logm(rho)-logm(tau))).real
    assert relative_entropy(rho,tau) == pytest.approx(expected,abs=1e-14)
    u = np.array([[1,1j],[1j,1]])/np.sqrt(2)
    assert relative_entropy(u@rho@u.conj().T,u@tau@u.conj().T) == pytest.approx(expected,abs=1e-14)


@pytest.mark.parametrize("leak", (1.,.5,1e-20))
def test_support_escape_is_infinite_even_for_small_positive_mass(leak):
    assert np.isinf(relative_entropy(np.diag([1-leak,leak]),np.diag([1.,0.])))


def test_contained_singular_support_and_zero_probability_terms():
    rho, tau = np.diag([.75,.25,0]), np.diag([.5,.5,0])
    expected = .75*np.log(1.5)+.25*np.log(.5)
    assert relative_entropy(rho,tau) == pytest.approx(expected,abs=1e-14)
    assert relative_entropy(rho,rho) == 0
    assert shannon_entropy([1,0]) == 0
    with pytest.raises(ValueError,match="faithful"):
        faithful_log(tau)


def test_small_positive_eigenvalues_are_not_deleted_or_floored():
    rho = np.diag([1.,1e-20])  # trace rounds to one at this precision
    assert von_neumann_entropy(rho) == pytest.approx(-1e-20*np.log(1e-20),rel=1e-14,abs=0)
    assert relative_entropy(np.diag([0.,1.]),np.diag([1.,1e-310])) == pytest.approx(-np.log(1e-310))


@pytest.mark.parametrize("bad", (
    np.diag([1.1,-.1]), np.eye(2), np.zeros((0,0)), np.ones((2,3)),
    np.array([[.5,.2],[0,.5]]), np.diag([np.nan,1]), np.diag([np.inf,0]),
    np.array([[True,False],[False,False]]), [1,0],
))
def test_invalid_density_is_rejected_by_all_entropy_entry_points(bad):
    for function in (density_matrix,von_neumann_entropy,faithful_log):
        with pytest.raises(ValueError):
            function(bad)
    with pytest.raises(ValueError):
        relative_entropy(bad,np.eye(2)/2)
    with pytest.raises(ValueError):
        relative_entropy(np.eye(2)/2,bad)


@pytest.mark.parametrize("bad", ([],[-.1,1.1],[.2,.2],[np.nan,0],[True,False],[1+0j,0]))
def test_invalid_classical_weights_are_rejected(bad):
    with pytest.raises(ValueError):
        probabilities(bad)
    with pytest.raises(ValueError):
        shannon_entropy(bad)


def test_partial_trace_matches_explicit_index_contraction_on_every_subset():
    dims = [2,3,2]
    rng = np.random.default_rng(17)
    x = rng.normal(size=(12,12))+1j*rng.normal(size=(12,12))
    rho = x@x.conj().T
    rho /= np.trace(rho)
    indices = list(itertools.product(*(range(d) for d in dims)))
    for mask in range(8):
        keep = [i for i in range(3) if mask & (1<<i)]
        out_dims = [dims[i] for i in keep]
        size = int(np.prod(out_dims))
        expected = np.zeros((size,size),complex)
        for i, left in enumerate(indices):
            for j, right in enumerate(indices):
                if all(left[t]==right[t] for t in range(3) if t not in keep):
                    a = np.ravel_multi_index(tuple(left[t] for t in keep),out_dims) if keep else 0
                    b = np.ravel_multi_index(tuple(right[t] for t in keep),out_dims) if keep else 0
                    expected[a,b] += rho[i,j]
        assert np.allclose(partial_trace(rho,dims,list(reversed(keep))),expected,rtol=0,atol=1e-14)


@pytest.mark.parametrize("dims,keep", (([2,2],[0,0]),([2,2],[2]),([2,2],[True]),
                                        ([2,0],[0]),([2,2.0],[0]),([True,4],[0]),([2,3],[0])))
def test_partial_trace_rejects_invalid_tensor_domains(dims,keep):
    with pytest.raises(ValueError):
        partial_trace(np.eye(4)/4,dims,keep)


def test_information_partitions_cannot_overlap():
    with pytest.raises(ValueError,match="disjoint"):
        mutual_information(np.eye(4)/4,[2,2],[0],[0,1])
    with pytest.raises(ValueError,match="disjoint"):
        conditional_mutual_information(np.eye(4)/4,[2,2],[0],[1],[0])


def test_direct_sum_entropy_and_relative_entropy_chain_rules():
    states = [np.array([[.7,.1],[.1,.3]]), np.diag([.2,.3,.5])]
    refs = [np.array([[.4,.05j],[-.05j,.6]]), np.diag([.5,.25,.25])]
    p,q = np.array([.3,.7]), np.array([.6,.4])
    rho,tau = direct_sum_state(p,states),direct_sum_state(q,refs)
    expected_s = -np.dot(p,np.log(p))-sum(w*np.trace(a@logm(a)).real for w,a in zip(p,states))
    expected_d = np.dot(p,np.log(p/q))+sum(
        w*np.trace(a@(logm(a)-logm(b))).real for w,a,b in zip(p,states,refs))
    assert von_neumann_entropy(rho) == pytest.approx(expected_s,abs=1e-14)
    assert relative_entropy(rho,tau) == pytest.approx(expected_d,abs=1e-14)
    # Unequal weights carry their own KL term: omitting it is not a simplification.
    assert np.dot(p,np.log(p/q)) > .1


def test_a3_flagged_reduction_does_not_assume_a_global_extension():
    # AB, BC and CA are all perfectly anticorrelated with uniform singletons.
    # No common classical ABC state can satisfy all three pair constraints.
    assert not any(a!=b and b!=c and c!=a for a,b,c in itertools.product((0,1),repeat=3))
    pair = np.diag([0.,.5,.5,0.])
    weights = np.array([1.,2.,3.])
    normalized = weights/weights.sum()
    rho = direct_sum_state(normalized,[pair]*3)
    tau = direct_sum_state(normalized,[np.eye(4)/4]*3)
    assert weights.sum()*relative_entropy(rho,tau) == pytest.approx(6*np.log(2),abs=1e-14)


def test_discarding_a_retained_label_destroys_markovity():
    # Exact classical A=J=C, each value equally likely. J is the collar label.
    rho = np.diag([.5,0,0,0,0,0,0,.5])
    assert abs(conditional_mutual_information(rho,[2,2,2],[0],[1],[2])) < 1e-14
    erased = partial_trace(rho,[2,2,2],[0,2])
    assert mutual_information(erased,[2,2],[0],[1]) == pytest.approx(np.log(2),abs=1e-14)


def test_invalid_states_and_empty_alignment_fail_under_python_optimization(tmp_path):
    code = """
import numpy as np
from quantum_information import relative_entropy
from collar_alignment.msa_characterizations import is_ec_aligned
for check in (lambda: relative_entropy(np.diag([1.1,-.1]),np.eye(2)/2),
              lambda: is_ec_aligned([])):
    try:
        check()
    except ValueError:
        continue
    raise SystemExit('invalid evidence accepted')
"""
    env = dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    result = subprocess.run([sys.executable,'-O','-c',code],cwd=tmp_path,env=env,
                            text=True,capture_output=True,timeout=30)
    assert result.returncode == 0, result.stderr
