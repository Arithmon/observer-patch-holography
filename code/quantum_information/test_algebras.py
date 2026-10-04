"""Index contractions, Gram identities and counterexamples independent of constructors."""

import itertools

import numpy as np
import pytest
from scipy.linalg import eigvalsh

from quantum_information.algebras import (
    FiniteAlgebra, algebra_intersection, embed_operator, matrix_units,
    separation_modulus, tensor_factor_algebra, tensor_separation_modulus,
)
from quantum_information import partial_trace
from geometry import null_net_standardness as standard


def test_subalgebra_region_cannot_silently_be_replaced_by_a_prefix():
    rho = np.diag([.5,0,.5,0])  # mixed A, pure B
    omega = standard.gns_vector(rho)
    annihilator = np.kron(np.eye(2),np.diag([0.,1.]))
    assert np.linalg.norm(annihilator) > 1
    assert np.linalg.norm(np.kron(annihilator,np.eye(4)) @ omega) == 0
    assert standard.subalgebra_separating(rho,[2,2],[0])
    assert not standard.subalgebra_separating(rho,[2,2],[1])
    # The complementary example excludes a fix that simply rejects B.
    reversed_rho = np.diag([.5,.5,0,0])
    assert not standard.subalgebra_separating(reversed_rho,[2,2],[0])
    assert standard.subalgebra_separating(reversed_rho,[2,2],[1])


def test_embedding_matches_independent_indices_on_every_subset():
    dims = [2,3,2]
    states = list(itertools.product(*(range(d) for d in dims)))
    rng = np.random.default_rng(14)
    for mask in range(8):
        region = [i for i in range(3) if mask & (1 << i)]
        rdims = [dims[i] for i in region]
        size = int(np.prod(rdims))
        small = rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
        expected = np.zeros((12,12),complex)
        for i, left in enumerate(states):
            for j, right in enumerate(states):
                if all(left[k] == right[k] for k in range(3) if k not in region):
                    a = np.ravel_multi_index(tuple(left[k] for k in region),rdims) if region else 0
                    b = np.ravel_multi_index(tuple(right[k] for k in region),rdims) if region else 0
                    expected[i,j] = small[a,b]
        assert np.array_equal(embed_operator(small,dims,list(reversed(region))),expected)
        # Independent pairing identity catches normalized versus ordinary trace.
        x = rng.normal(size=(12,12))+1j*rng.normal(size=(12,12))
        rho = x@x.conj().T
        rho /= np.trace(rho)
        assert np.trace(rho@expected) == pytest.approx(np.trace(partial_trace(rho,dims,region)@small))


@pytest.mark.parametrize("dims,region", (([2,2],[0,0]),([2,2],[2]),([2,2],[-1]),
                                        ([2,2],[True]),([2,0],[0]),([True,2],[0])))
def test_invalid_region_declarations_are_rejected(dims,region):
    with pytest.raises(ValueError):
        embed_operator(np.eye(2),dims,region)
    with pytest.raises(ValueError):
        standard.subalgebra_separating(np.eye(4)/4,dims,region)


def test_tensor_separation_formula_matches_dense_action_on_all_regions():
    dims = [2,2,2]
    rng = np.random.default_rng(43)
    x = rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho = x@x.conj().T + np.eye(8)
    rho /= np.trace(rho)
    omega = standard.gns_vector(rho)
    for mask in range(8):
        region = [i for i in range(3) if mask & (1 << i)]
        small_size = 2**len(region)
        # Construct physical matrix units by bit indices, independently of embed_operator.
        basis = []
        for i in range(small_size):
            for j in range(small_size):
                a = np.zeros((8,8))
                for row in range(8):
                    for col in range(8):
                        r = tuple((row >> (2-k)) & 1 for k in range(3))
                        c = tuple((col >> (2-k)) & 1 for k in range(3))
                        ri = sum(r[k] << (len(region)-1-t) for t,k in enumerate(region))
                        ci = sum(c[k] << (len(region)-1-t) for t,k in enumerate(region))
                        if ri == i and ci == j and all(r[k] == c[k] for k in range(3) if k not in region):
                            a[row,col] = 1
                basis.append(np.kron(a,np.eye(8)))
        action = np.column_stack([a@omega for a in basis])
        expected = np.linalg.svd(action,compute_uv=False)[-1]/np.linalg.norm(basis[0])
        assert tensor_separation_modulus(rho,dims,region) == pytest.approx(expected,rel=1e-12)


def test_general_separation_is_invariant_under_basis_scaling_and_shear():
    identity, z = np.eye(2), np.diag([1.,-1.])
    omega = np.ones(2)/np.sqrt(2)
    for basis in ([identity,z], [identity,identity+z], [1e-80*identity,1e80*z]):
        assert standard.separating_modulus(basis,omega) == pytest.approx(1/np.sqrt(2))
    basis = [identity,(1+.4j)*identity+2*z]
    gram = np.array([[np.vdot(a,b) for b in basis] for a in basis])
    action = np.column_stack([a@omega for a in basis])
    expected = np.sqrt(eigvalsh(action.conj().T@action,gram)[0])
    assert separation_modulus(basis,omega) == pytest.approx(expected)
    with pytest.raises(ValueError,match="dependent"):
        separation_modulus([identity,2*identity],omega)


def test_standardness_does_not_allocate_the_full_left_regular_representation(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("dense GNS representation allocated")
    monkeypatch.setattr(standard,"gns_vector",forbidden)
    monkeypatch.setattr(standard,"left_action",forbidden)
    rho = np.eye(64)/64
    assert standard.full_algebra_standard(rho) == (True,True)
    assert standard.subalgebra_separating(rho,[4,4,4],[0,2])


@pytest.mark.parametrize("bad", ([],[np.eye(2),np.eye(2)], [np.diag([1.,-1.])],
    [np.eye(2),np.array([[0.,1.],[0.,0.]])],
    [np.eye(2),np.array([[0.,1.],[1.,0.]]),np.diag([1.,-1.])],
    [np.eye(2),np.full((2,2),np.nan)], [np.eye(2,dtype=bool)]))
def test_invalid_operator_algebras_are_rejected(bad):
    with pytest.raises(ValueError):
        FiniteAlgebra(bad)


def test_intersection_is_the_shared_algebra_and_not_a_tensor_union():
    left = tensor_factor_algebra([2,2],[0])
    right = tensor_factor_algebra([2,2],[1])
    full = FiniteAlgebra(matrix_units(4))
    common = algebra_intersection([left,right])
    assert common.dimension == 1
    assert algebra_intersection([left,full]).dimension == 4
    assert algebra_intersection([left,left]).dimension == 4
    x = np.arange(16).reshape(4,4)
    assert np.allclose(common.project(x),np.trace(x)*np.eye(4)/4,atol=1e-13)
    with pytest.raises(ValueError):
        algebra_intersection([])


def test_nonfaithful_gns_space_is_a_proper_cyclic_subspace_of_hs():
    rho = np.diag([1.,0.])
    omega = standard.gns_vector(rho)
    actions = [np.kron(a,np.eye(2))@omega for a in matrix_units(2)]
    assert np.linalg.matrix_rank(np.column_stack(actions)) == 2
    assert omega.size == 4
    assert standard.full_algebra_standard(rho) == (False,False)


def test_positive_subnormal_separation_modulus_is_not_lost_in_division():
    tiny = np.nextafter(0.,1.)
    modulus = tensor_separation_modulus(np.diag([1.,tiny]),[2],[0])
    assert modulus > 0
    assert modulus == pytest.approx(np.sqrt(tiny)/np.sqrt(2),rel=1e-14,abs=0)
    assert standard.full_algebra_standard(np.diag([1.,tiny])) == (True,True)
    assert standard.subalgebra_separating(np.diag([1.,tiny]),[2],[0])


def test_rank_ambiguity_is_not_mistaken_for_faithfulness():
    h = np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2
    rho = h@np.diag([.25,.5,.25,0])@h.T
    try:
        result = standard.full_algebra_standard(rho)
    except ValueError as error:
        assert "support is numerically unresolved" in str(error)
    else:
        assert result == (False,False)
