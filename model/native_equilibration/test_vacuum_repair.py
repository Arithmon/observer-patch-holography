"""Source-ground compatibility, full dilation and independent work controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import eigh
import sympy as s

from . import vacuum_repair as repair


def test_actual_native_geometry_records_make_a_joint_uniform_bit():
    laws = repair.native_geometry()
    assert [x["law"] for x in laws] == ["C", "E"]
    for law in laws:
        assert len(law["branches"]) == 4
        assert len({tuple(x["label"]) for x in law["branches"]}) == 4
        for first in (0, 1):
            assert sorted(row["final"] for row in law["branches"] if row["initial"] == first) == [0, 1]


def test_hamiltonian_is_true_compression_not_truncated_coordinate_square():
    cutoff, band, r = 4, 2, F(9, 10)
    omega = math.sqrt(14)
    a = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    q, p = (a+a.T)/math.sqrt(2*omega), 1j*math.sqrt(omega/2)*(a.T-a)
    eye = np.eye(cutoff)
    kron = lambda a, b, c: np.kron(np.kron(a, b), c)
    wrong = float(r)/2*(kron(p@p, eye, eye)+kron(eye, p@p, eye)+kron(eye, eye, p@p))
    wrong += (8+6*float(r)**3)/(2*float(r))*(kron(q@q, eye, eye)+kron(eye, q@q, eye)+kron(eye, eye, q@q))
    wrong -= 4/float(r)*(kron(q, q, eye)+kron(eye, q, q))
    actual = repair.hamiltonian(cutoff, band, r)
    assert np.linalg.norm(wrong-actual) > 10
    # Directly project a larger physical oscillator algebra, whose q,p can
    # reach one additional level. This retains the lost top diagonal.
    large = cutoff+1
    b = np.diag(np.sqrt(np.arange(1, large)), 1)
    qq, pp = (b+b.T)/math.sqrt(2*omega), 1j*math.sqrt(omega/2)*(b.T-b)
    take = np.arange(cutoff)
    Q, P = (qq@qq)[np.ix_(take, take)], (pp@pp)[np.ix_(take, take)]
    direct = float(r)/2*(kron(P, eye, eye)+kron(eye, P, eye)+kron(eye, eye, P))
    direct += (8+6*float(r)**3)/(2*float(r))*(kron(Q, eye, eye)+kron(eye, Q, eye)+kron(eye, eye, Q))
    direct -= 4/float(r)*(kron(q, q, eye)+kron(eye, q, q))
    np.testing.assert_allclose(actual, direct, atol=2e-14, rtol=0)


@pytest.mark.parametrize("band", range(4))
def test_fixed_algebra_compatibility_tangent_and_entropy_response(band):
    r, step = F(1), F(1, 10000)
    response = repair.response(4, band, r)
    left = repair.ground(4, band, r-step)
    right = repair.ground(4, band, r+step)
    ds = (repair.entropy(right[2], 4)-repair.entropy(left[2], 4))/(2*float(step))
    de = (right[1]-left[1])/(2*float(step))
    assert ds == pytest.approx(response["entropy_derivative"], abs=3e-8)
    assert de == pytest.approx(response["Hellmann_Feynman_force"], abs=3e-7)


def test_literal_native_preparation_has_a_variational_gaussian_fidelity_bound():
    vector, row = repair.prepare(8, 0, F(11, 10))
    assert row["native_pulses"] == 12291
    certificate = row["certificate"]
    assert 0 <= certificate["infidelity_upper"] < F(1, 100000)
    assert not certificate["normalization_or_probability_repair_performed"]
    # Completely unrelated excited number data are rejected by the
    # independent energy gate rather than accepted as a prepared vacuum.
    excited = np.zeros(512, complex); excited[-1] = 1
    bad = repair.energy_certificate(excited, 8, 0, F(11, 10))
    assert bad["infidelity_upper"] > 1
    with pytest.raises(ValueError, match="norm"):
        repair.energy_certificate(2*vector, 8, 0, F(11, 10))


def test_exact_binary_rayleigh_enclosure_covers_independent_dense_expectation():
    vector, row = repair.prepare(4, 2, F(9, 10))
    n = float(np.vdot(vector, vector).real)
    energy = float((vector.conj()@repair.hamiltonian(4, 2, F(9, 10))@vector).real/n)
    low, high = row["certificate"]["Rayleigh_energy_interval"]
    # The direct double evaluation has its own rounding; the exact
    # rational enclosure is much narrower than one double ulp.
    assert float(low)-2e-14 <= energy <= float(high)+2e-14


def test_literal_source_prepared_entropy_tangent_uses_same_algebra_and_geometry():
    step = F(1, 10000)
    left, _ = repair.prepare(8, 2, F(1)-step)
    right, _ = repair.prepare(8, 2, F(1)+step)
    slope = (repair.entropy(right, 8)-repair.entropy(left, 8))/(2*float(step))
    assert slope == pytest.approx(repair.response(8, 2, F(1))["entropy_derivative"], abs=3e-8)


def test_native_reset_dilation_covers_every_input_and_untouched_reference():
    dilation = repair.reset_dilation(F(9, 10))
    V = dilation["source_isometry"]
    assert V.shape == (64, 8) and dilation["all_input_columns"] == 8
    assert dilation["ordinary_swap_all_column_error"] < 1e-12
    assert dilation["all_input_isometry_error"] < 1e-10
    assert dilation["fresh_blank_qubits"] == dilation["discarded_environment_qubits"] == 3
    # Input maximally entangled with an untouched reference. Field is
    # reset; its original information remains in environment/reference.
    wave = (V/math.sqrt(8)).reshape(8, 8, 8)  # env,field,reference
    reference = np.einsum("efr,efs->rs", wave, wave.conj())
    np.testing.assert_allclose(reference, np.eye(8)/8, atol=2e-12, rtol=0)


def test_joint_repair_checks_all_geometry_field_inputs_not_only_stationary_input():
    actual, target, report = repair.joint_reset_control()
    assert actual.shape == target.shape == (512, 16)
    assert report["all_joint_input_columns"] == 16
    assert report["complete_native_record_branches"] == 4
    assert report["source_feedback_reads_only_committed_native_geometry_record"]
    np.testing.assert_allclose(actual.conj().T@actual, np.eye(16), atol=2e-11, rtol=0)
    # Keep record/environment and trace them only for this operational
    # joint-state comparison. Every matrix unit is checked independently.
    expected = np.zeros((16, 16), complex)
    for bit, r in enumerate((F(9, 10), F(11, 10))):
        state = repair.ground(2, 2, r)[2]
        expected[8*bit:8*(bit+1), 8*bit:8*(bit+1)] = np.outer(state, state)/2
    ideal_channel = np.einsum("ab,ij->abij", expected, np.eye(16))
    for isometry in (actual, target):
        V = isometry.reshape(4, 2, 8, 8, 16)  # native record,geom,env,field,input
        channel = np.einsum("rgefi,rhekj->gfhkij", V, V.conj()).reshape(16, 16, 16, 16)
        np.testing.assert_allclose(channel, ideal_channel, atol=2e-11, rtol=0)


def test_passivity_cycle_work_matches_independent_covariance_energy():
    values = repair.work()
    assert values["total_mean_quench_work_diagnostic"] > 0
    total = 0.
    rm, rp = .9, 1.1
    radial = (8-4*math.sqrt(2), 8., 8+4*math.sqrt(2))
    for lam, multiplicity in repair.BANDS:
        for k in radial:
            wm, wp = math.sqrt(k+float(lam)*rm**3), math.sqrt(k+float(lam)*rp**3)
            Xm, Pm = rm/(2*wm), wm/(2*rm)
            Xp, Pp = rp/(2*wp), wp/(2*rp)
            Em_to_p = rp*Pm/2+(k+float(lam)*rp**3)*Xm/(2*rp)
            Ep_to_m = rm*Pp/2+(k+float(lam)*rm**3)*Xp/(2*rm)
            total += multiplicity*(Em_to_p-wp/2+Ep_to_m-wm/2)/4
    assert total == pytest.approx(values["total_mean_quench_work_diagnostic"], abs=2e-14)
    assert repair.work(contrast=F())["total_mean_quench_work_diagnostic"] == 0


def test_literal_source_ray_work_interval_is_independent_of_ideal_gaussian_work():
    prepared = {(r, band): repair.prepare(4, band, r)[0]
                for r in (F(9, 10), F(11, 10)) for band in range(4)}
    low, high = repair.source_work_interval(prepared, cutoff=4)
    assert low > F(1, 1000) and low < high
    direct = 0.
    for band, (_, multiplicity) in enumerate(repair.BANDS):
        minus, plus = prepared[F(9, 10), band], prepared[F(11, 10), band]
        Hm, Hp = repair.hamiltonian(4, band, F(9, 10)), repair.hamiltonian(4, band, F(11, 10))
        direct += multiplicity*float((minus.conj()@(Hp-Hm)@minus/np.vdot(minus, minus)
                         +plus.conj()@(Hm-Hp)@plus/np.vdot(plus, plus)).real)/4
    assert float(low)-2e-13 <= direct <= float(high)+2e-13


def test_retained_complete_preparation_family_is_independently_certified():
    family = {(r, band): repair.ground(8, band, r)[2]
              for r in (F(9, 10), F(11, 10)) for band in range(4)}
    bounds, work = repair.certify_prepared(family)
    assert len(bounds) == 2 and work[0] > F(1, 1000)
    corrupt = dict(family)
    excited = np.zeros(512); excited[-1] = 1
    corrupt[F(11, 10), 0] = excited
    with pytest.raises(ValueError, match="fidelity"):
        repair.certify_prepared(corrupt)
    corrupt = dict(family); corrupt.pop((F(9, 10), 3))
    with pytest.raises(ValueError, match="complete"):
        repair.certify_prepared(corrupt)


@pytest.mark.parametrize("bad", [True, 1., F(0), F(3), 1+0j])
def test_original_geometry_domain(bad):
    with pytest.raises(ValueError):
        repair.hamiltonian(4, 2, bad)


@pytest.mark.parametrize("cutoff,band", [(True, 2), (3, 2), (16, 2), (4, False), (4, 4)])
def test_fixed_source_code_and_band_domain(cutoff, band):
    with pytest.raises(ValueError):
        repair.hamiltonian(cutoff, band, F(1))
