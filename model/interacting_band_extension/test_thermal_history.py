"""Independent full-polynomial, heat, record and hostile-input controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from .thermal_history import (cell, evolve, finite_gibbs_seed, generator, history,
                              preparation_bracket)


def oracle(d):
    """Padded oscillator algebra, independent of exact-form construction."""
    n = d+6
    ladder = np.diag(np.sqrt(np.arange(1, n)), 1)
    q, p = (ladder+ladder.T)/2, 1j*(ladder.T-ladder)
    q2, q4, p2 = (matrix[:d, :d] for matrix in (q@q, q@q@q@q, p@p))
    h1, h2 = p2/2+q2/2+q4/16, p2/2+2*q2+q4/16
    return np.kron(np.eye(d), h1)+np.kron(h2, np.eye(d))+np.kron(q2, q2)/4, np.kron(q[:d, :d], np.eye(d))


def test_seed_encloses_independent_canonical_heat_and_reference():
    seed, budget = finite_gibbs_seed(4, F(2))
    h, _ = oracle(4)
    target = expm(-h)/np.sqrt(np.trace(expm(-2*h)))
    error = np.linalg.norm(seed-target)
    assert error < float(budget['canonical_seed_vector_error_upper']) < 1e-8
    assert np.linalg.norm(seed-seed.T) < 1e-15
    # The reference is canonical, not an arbitrary energy eigenvector frame.
    assert np.linalg.norm(seed-np.diag(np.diag(seed))) > .01


def test_generators_include_all_escaped_rows_and_original_mode_order():
    d = 3
    for kind, coefficient in (('hamiltonian', F(1, 16)), ('coordinate', F(1, 4))):
        operation = generator(d, kind, coefficient)
        h, q = oracle(operation.extended)
        target = (-1j*float(coefficient)*h if kind == 'hamiltonian' else 1j*float(coefficient)*q)
        state = np.eye(d*d, dtype=complex)
        real, imag = operation.apply_interval(state)
        expected = target[:, operation.retained]
        # Independent sqrt arithmetic may differ by a few ulps.
        assert np.max(np.abs((real[0]+real[1])/2+1j*(imag[0]+imag[1])/2-expected)) < 2e-14
        assert np.linalg.norm(expected[np.setdiff1d(np.arange(len(target)),operation.retained)]) > .1
        assert np.allclose(operation.apply_float(state),expected[operation.retained],atol=1e-14)


def test_two_read_joint_work_reference_endpoint_against_direct_oracle():
    d = 3
    h, q = oracle(d)
    initial = expm(-h)/np.sqrt(np.trace(expm(-2*h)))
    kick, control, evolution = expm(-1j*q/4), expm(1j*q/4), expm(-1j*h/16)
    plus, minus = expm(1j*q/2), expm(-1j*q/2)
    kraus = (((1+1j)*plus+(1-1j)*minus)/np.sqrt(8),
             ((1-1j)*plus+(1+1j)*minus)/np.sqrt(8))
    first = evolution@kick@initial
    expected = tuple(kraus[last]@evolution@(control if bit else np.eye(d*d))@kraus[bit]@first
                     for bit in (0,1) for last in (0,1))
    outputs, budget = history(initial.astype(complex), -1, occupation=d, cells=2)
    assert max(np.linalg.norm(a-b) for a,b in zip(outputs,expected,strict=True)) < 2e-12
    assert budget['all_record_branches'] == 4
    assert budget['retained_reference_columns'] == d*d
    assert budget['occupation_residual'] > .02
    # A record law alone can miss work/reference changes. Test all amplitudes.
    wrong = tuple(kraus[last]@evolution@kraus[bit]@first for bit in (0,1) for last in (0,1))
    assert max(np.linalg.norm(a-b) for a,b in zip(outputs,wrong,strict=True)) > .02


def test_full_residual_bounds_larger_space_and_detects_trash_path(monkeypatch):
    d, n = 3, 9
    initial = np.zeros((d*d, 2), complex); initial[-1,0]=1
    endpoint, budget = evolve(initial,'hamiltonian',F(1,32),occupation=d,cells=4)
    h, _ = oracle(n)
    ids = [i+n*j for j in range(d) for i in range(d)]
    lifted = np.zeros((n*n,2),complex);lifted[ids]=initial
    exact = expm(-1j*h/32)@lifted
    lifted[ids]=endpoint
    assert np.linalg.norm(exact-lifted) <= sum(budget.values())
    operation=generator(d,'hamiltonian',F(1,32))
    monkeypatch.setattr(type(operation),'apply_float',lambda self,value: np.zeros_like(value))
    _, corrupt = cell(initial,operation)
    assert corrupt['finite_residual'] > .1


def test_partition_certificate_refines_and_free_target_is_enclosed():
    cold = preparation_bracket(8,12,F(4))
    assert float(cold['canonical_preparation_trace_distance_upper']) < .002
    free = preparation_bracket(4,0,F(1),F(0))
    exact = np.exp(-1.5)/((1-np.exp(-1))*(1-np.exp(-2)))
    assert float(free['infinite_partition_lower']) < exact <= float(free['infinite_partition_upper'])*(1+1e-15)
    assert not free['eigensolver_accuracy_assumed']


@pytest.mark.parametrize('bad',[True,2.0,0,-1,33])
def test_bad_occupation_fails_closed(bad):
    with pytest.raises(ValueError): generator(bad,'hamiltonian',F(1))


@pytest.mark.parametrize('value',[np.nan,np.inf,1e308])
def test_nonfinite_and_overflowing_amplitudes_fail_closed(value):
    state=np.zeros((4,1),complex);state[0,0]=value
    with np.errstate(over='ignore',invalid='ignore'):
        with pytest.raises((ValueError,ArithmeticError)):
            history(state,occupation=2)


def test_boolean_masked_and_nonexact_parameters_fail_closed():
    for state in (np.zeros((4,1),bool), np.ma.array(np.ones((4,1)),mask=True), [[1],[0],[0],[0]]):
        with pytest.raises(ValueError): history(state,occupation=2)
    for kind, coefficient in (('trash',F(1)),('coordinate',True),('coordinate',.25)):
        with pytest.raises(ValueError): generator(2,kind,coefficient)
    with pytest.raises(ValueError): history(np.eye(4,dtype=complex),sign=True,occupation=2)
