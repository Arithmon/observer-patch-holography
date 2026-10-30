"""Independent matrix, complete-export and same-reference composition controls."""
import numpy as np
import pytest
from scipy.linalg import expm
from . import recovery as r
from ..field_thermal import thermal
from ..full_fock_field.field import FieldGeometry
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation


@pytest.mark.parametrize('bits,columns', [(2,2),(3,4)])
def test_general_complex_isometry_on_every_input_and_reference(bits, columns):
    rng=np.random.default_rng(71+bits)
    raw=rng.normal(size=(1<<bits,columns))+1j*rng.normal(size=(1<<bits,columns))
    target=np.linalg.qr(raw)[0]
    word=r.compile_isometry(target)
    initial=np.zeros_like(target);initial[:columns]=np.eye(columns)
    run=r.execute_complete(word,bits,initial)
    np.testing.assert_allclose(run.success,target,atol=2e-10,rtol=0)
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,np.eye(columns),atol=2e-10)


def test_lazy_failure_family_matches_each_actual_sender_reset_Kraus_branch(monkeypatch):
    generator=np.zeros((6,6),complex)
    generator[0,2]=generator[2,0]=.4
    generator[1,3]=generator[3,1]=.7
    full=expm(1j*generator)
    monkeypatch.setattr(compiler,'native_pulse',lambda *_:full.copy())
    initial=np.array([[1,2j],[3j,4]],complex)/6
    tape=[compiler.pulse('X',.2,[0]),compiler.pulse('Z',.3,[0])]
    expected,events=thermal.execute_with_exports(tape,1,initial)
    run=r.execute_complete(tape,1,initial)
    np.testing.assert_allclose(run.success,expected,atol=1e-15)
    for i in range(2):
        for a,b in zip(run.failure(i),events[i],strict=True):np.testing.assert_allclose(a,b,atol=1e-15)
    np.testing.assert_allclose(run.failure_gram,sum(x.conj().T@x for e in events for x in e),atol=1e-15)
    # A later change to source lookup cannot rewrite already captured outcomes.
    monkeypatch.setattr(compiler,'native_pulse',lambda *_:np.eye(6))
    np.testing.assert_allclose(run.failure(0),events[0],atol=1e-15)
    assert run.failure_probability>.01


def test_nonadjacent_pair_exports_keep_incoherent_branches_and_immutable_snapshots(monkeypatch):
    rng = np.random.default_rng(407)
    bits = 3
    initial = rng.normal(size=(8, 3)) + 1j*rng.normal(size=(8, 3))
    raw = rng.normal(size=(6, 6)) + 1j*rng.normal(size=(6, 6))
    full = expm(.2j*(raw+raw.conj().T))
    tape = [compiler.pulse('J', .2, [2, 0]),
            compiler.pulse('X', .4, [1]), compiler.pulse('J', .3, [0, 1])]
    monkeypatch.setattr(compiler, 'native_pulse', lambda *_: full)
    run = r.execute_complete(tape, bits, initial)
    state = initial.copy()
    gram = np.zeros((3, 3), complex)

    # Direct bit-index contractions do not reuse either execution helper.
    for event_index, event in enumerate(tape):
        targets = event['targets']
        remainder = [i for i in reversed(range(bits)) if i not in targets]
        d = 1 << len(targets)
        arranged = np.zeros((d, 1 << len(remainder), 3), complex)
        addresses = []
        for mask, row in enumerate(state):
            local = sum(((mask >> bit) & 1) << (len(targets)-1-j)
                        for j, bit in enumerate(targets))
            outside = sum(((mask >> bit) & 1) << (len(remainder)-1-j)
                          for j, bit in enumerate(remainder))
            addresses.append((local, outside))
            arranged[local, outside] = row
        output = np.einsum('ij,jkr->ikr', full[:, :d], arranged)
        for actual, expected in zip(run.failure(event_index), output[d:], strict=True):
            np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=0)
            gram += expected.conj().T@expected
        state = np.array([output[local, outside] for local, outside in addresses])

    np.testing.assert_allclose(run.success, state, atol=1e-13, rtol=0)
    np.testing.assert_allclose(run.failure_gram, gram, atol=1e-13, rtol=0)
    np.testing.assert_allclose(state.conj().T@state+gram, initial.conj().T@initial,
                               atol=1e-13, rtol=0)
    branches = run.failure(0)
    coherent = sum(branches)
    assert np.linalg.norm(coherent.conj().T@coherent-
                          sum(x.conj().T@x for x in branches)) > 1

    # The retained event must survive mutation of every caller-owned object.
    initial[:] = 999
    full[:] = 0
    tape[0]['targets'][:] = [0, 2]
    for actual, expected in zip(run.failure(0), branches, strict=True):
        np.testing.assert_array_equal(actual, expected)


def test_complex_reference_preparation_and_each_global_pointer_mode():
    rng = np.random.default_rng(303)
    raw = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    rho = raw@raw.conj().T + np.eye(4)
    rho /= np.trace(rho)
    assert np.linalg.norm(rho.imag) > .01
    service = r.RecoveryService(rho, thermal.car_restriction(2, (1,)))
    for pointer in (0, 1):
        run, resources = r.reference_preparation(service, pointer_mode=pointer)
        rows = run.success.reshape(2, -1, 4)
        for bit, branch in enumerate(rows):
            projector = np.diag([int(((mask >> pointer) & 1) == bit) for mask in range(4)])
            # Trace all retained purifier/environment axes independently.
            density = branch.T@branch.conj()
            np.testing.assert_allclose(density, projector@rho@projector, atol=2e-10, rtol=0)
        assert abs(np.vdot(run.success, run.success).real+run.failure_probability-1) < 2e-10
        assert resources['native_pulses'] == len(run.word)
        assert resources['proper_code_load_events'] == len(run.word)
        assert resources['proper_code_export_events'] == len(run.word)


@pytest.fixture(scope='module')
def dimer():
    return r.field_service()


@pytest.fixture(scope='module')
def executed(dimer):
    return dimer.execute()


def test_complete_dimer_recovery_Choi_matches_independent_CAR_formula(dimer,executed):
    rho=thermal.ThermalField().target()['reference']
    a=[np.array(k) for k in thermal.car_restriction(4,(0,1))]
    local=sum(k@rho@k.conj().T for k in a)
    from scipy.linalg import fractional_matrix_power
    root=fractional_matrix_power(rho,.5)
    inverse=fractional_matrix_power(local,-.5)
    kraus=[root@k.conj().T@inverse for k in a]
    actual=executed.success.reshape(4,16,4)
    for wanted,found in zip(kraus,actual,strict=True):np.testing.assert_allclose(found,wanted,atol=2e-9,rtol=0)
    np.testing.assert_allclose(sum(k@local@k.conj().T for k in actual),rho,atol=2e-9,rtol=0)
    np.testing.assert_allclose(executed.success.conj().T@executed.success+executed.failure_gram,np.eye(4),atol=2e-9)
    assert len(executed.word)==len(dimer.word)
    assert executed.failure_weights.shape==(len(dimer.word),4)


def test_compiled_recovery_implements_observation_effects_without_retraction(dimer,executed):
    matrices=executed.success.reshape(4,16,4)
    for at in range(4):
        n=np.diag([(mask>>at)&1 for mask in range(16)])
        effect=sum(k.conj().T@n@k for k in matrices)
        np.testing.assert_allclose(effect,dimer.observation.observe(n),atol=2e-9)
    n=np.diag([mask&1 for mask in range(16)])
    regional=np.diag([mask&1 for mask in range(4)])
    assert np.linalg.norm(dimer.observation.observe(n)-regional)>1e-5
    assert np.linalg.norm(dimer.observation.observe(np.eye(16))-np.eye(4))<1e-12


def test_source_prepared_reference_and_both_pointer_records(dimer):
    run,resources=r.reference_preparation(dimer)
    # Trace every retained environment and purifier, keep both recorded outcomes.
    rows=run.success.reshape(2,-1,16)
    actual=[x.T@x.conj() for x in rows]
    rho=dimer.observation.source_reference
    p=[np.diag([int((mask&1)==bit) for mask in range(16)]) for bit in (0,1)]
    for found,projector in zip(actual,p,strict=True):
        np.testing.assert_allclose(found,projector@rho@projector,atol=2e-9)
    assert abs(sum(np.trace(x).real for x in actual)+run.failure_probability-1)<2e-9
    assert resources['all_pointer_and_export_outputs_retained']
    assert resources['native_pulses']>len(dimer.word)


def test_nested_source_words_preserve_both_environments_and_direct_Petz_channel(dimer):
    local=dimer.observation.target_reference
    inner=r.RecoveryService(local,thermal.car_restriction(2,(0,)))
    word,bits=r.nested_word(dimer,inner)
    initial=np.zeros((1<<bits,2),complex);initial[:2]=np.eye(2)
    run=r.execute_complete(word,bits,initial)
    kraus=run.success.reshape(-1,16,2)
    direct=ReferenceObservation(dimer.observation.source_reference,thermal.car_restriction(4,(0,)))
    for i in range(2):
        for j in range(2):
            x=np.zeros((2,2),complex);x[i,j]=1
            np.testing.assert_allclose(sum(k@x@k.conj().T for k in kraus),direct.recover(x),atol=2e-9)
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,np.eye(2),atol=2e-9)
    wrong=r.RecoveryService(np.eye(4)/4,thermal.car_restriction(2,(0,)))
    with pytest.raises(ValueError,match='intermediate reference'):r.nested_word(dimer,wrong)


@pytest.mark.parametrize('damage',['nonfinite','nonunitary'])
def test_bad_full_primitive_fails_before_output(monkeypatch,damage):
    full=np.eye(6,dtype=complex);full[0,0]=np.nan if damage=='nonfinite' else 2
    monkeypatch.setattr(compiler,'native_pulse',lambda *_:full)
    with pytest.raises((ValueError,ArithmeticError)):
        r.execute_complete([compiler.pulse('X',.2,[0])],1,np.eye(2))


def test_wrong_compiled_word_is_rejected_on_the_complete_input_domain():
    service=r.field_service(thermal.ThermalField(FieldGeometry((1,1,1))),(0,))
    service.word=service.word+[compiler.pulse('X',.5,[0])]
    with pytest.raises(ArithmeticError,match='complete channel'):service.execute()


@pytest.mark.parametrize('value',[np.ones((3,2)),np.eye(8)[:,:5],np.eye(4)*2,np.full((4,2),np.nan),np.ma.array(np.eye(4))])
def test_invalid_isometry_domains_fail_closed(value):
    with pytest.raises((ValueError,ArithmeticError)):r.compile_isometry(value)
