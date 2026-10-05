import copy
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import pytest
import sympy as sp

from dispersion_feasibility import bounds, check, kinematics as k

REPORT = Path(__file__).with_name("report.json")


def test_exact_support_moments():
    # Exact polynomial identities; three sampled directions would be inadequate.
    phi = (1+sp.sqrt(5))/2
    vectors = []
    for axis in range(3):
        for sign in (-1, 1):
            vectors.append(tuple(sp.Integer(sign if i == axis else 0) for i in range(3)))
    for signs in itertools.product((-1, 1), repeat=3):
        v = (sp.Rational(signs[0],2), signs[1]*(phi-1)/2, signs[2]*phi/2)
        vectors.extend((v, v[1:]+v[:1], v[2:]+v[:2]))
    x = sp.symbols("x y z")
    for v in vectors:
        assert sp.simplify(sum(t*t for t in v)) == 1
    for degree, factor in ((2,10), (4,6)):
        difference = sp.Poly(sp.expand(sum(sum(v[i]*x[i] for i in range(3))**degree for v in vectors)
                                      - factor*sum(t*t for t in x)**(degree//2)), *x)
        assert all(sp.simplify(c) == 0 for c in difference.coeffs())


def test_budget_and_nominal_inputs():
    b = bounds.verify_budget()
    assert b["d_lo"] < k.D < b["d_hi"]
    assert 500000 < k.MASS_EV < 600000
    assert b["small_leg_gap"] > 4*b["feasible_gap"]


def test_report_independent_replay():
    assert check.check(check.load(REPORT))


@pytest.mark.parametrize("variant", k.VARIANTS)
@pytest.mark.parametrize("momentum", ("1e17", "3e17", "1e20"))
def test_precision_refinement(momentum, variant):
    first = k.witness(momentum,variant,(1,2,3),70)
    second = k.witness(momentum,variant,(1,2,3),100)
    with mp.workdps(110):
        assert abs(mp.mpf(first["collinear_soft_energy_eV"])-mp.mpf(second["collinear_soft_energy_eV"])) < mp.mpf("1e-32")


def test_independent_variational_minimum():
    # Direct minimization of the leading convex y objective via its derivative.
    with mp.workdps(70):
        m,d = k.number(k.MASS_EV),k.number(k.D)
        for value in k.MOMENTA:
            momentum = mp.mpf(value)
            epsilon,x = k.leading(momentum,"photon_electron_positron")
            y = x*(1-x)
            derivative = lambda t: -m*m/t**2+3*d*momentum**4
            if y < mp.mpf("0.25"):
                assert abs(derivative(y)/(3*d*momentum**4)) < mp.mpf("1e-60")
            else:
                assert derivative(y) <= 0
            assert 2*m*m/y**3 > 0
            # Full one-dimensional search with 200 logarithmic shares is an
            # independent numerical control, not the global proof.
            for i in range(201):
                trial = mp.mpf("0.25")*mp.power(10,-mp.mpf(i)/20)
                assert (m*m/trial+3*d*momentum**4*trial)/(4*momentum) >= epsilon*(1-mp.mpf("1e-60"))


def test_equal_share_shortcut_fails():
    with mp.workdps(60):
        momentum = mp.mpf("1e19")
        m,d = k.number(k.MASS_EV),k.number(k.D)
        actual,_ = k.leading(momentum,"photon_electron_positron")
        equal = m*m/momentum+3*d*momentum**3/16
        assert equal/actual > 990


def test_energy_conserving_stationary_maximum_is_rejected():
    # A fake symmetric "threshold" can conserve energy and be stationary.
    # Above the bifurcation it is a maximum in the sharing direction.
    packet = check.load(REPORT)
    row = next(r for r in packet["witnesses"] if
               r["variant"] == "photon_electron_positron" and
               r["hard_momentum_eV"] == "10000000000000000000.0" and
               r["direction"] == [1, 0, 0])
    with mp.workdps(90):
        momentum, m = mp.mpf("1e19"), k.number(k.MASS_EV)
        a, d = mp.sqrt(k.number(k.A2)), k.number(k.D)
        dots = k.projections((1, 0, 0), 90)
        hard_correction = k.shell(momentum, 0, a, dots)[1]
        def residual(s):
            return (2*k.shell((momentum-s)/2, m, a, dots)[1]
                    - hard_correction - k.shell(s, 0, a, dots)[1] - 2*s)
        soft = mp.findroot(residual, m*m/momentum+3*d*momentum**3/16)
        assert abs(residual(soft)) < mp.mpf("1e-60")
        row.update(collinear_soft_energy_eV=mp.nstr(k.shell(soft,0,a,dots)[0],35),
                   outgoing_small_fraction="0.5", energy_residual_eV="0")
    with pytest.raises(ValueError, match="stationary point is not a minimum"):
        check.check(packet)


def test_photon_only_hard_emission_and_shared_shell_stability():
    with mp.workdps(80):
        a,m = mp.sqrt(k.number(k.A2)),k.number(k.MASS_EV)
        initial,emitted = mp.mpf("3e17"),mp.mpf("2e17")
        for n in k.DIRECTIONS:
            dots = k.projections(n,80)
            photon = k.shell(emitted,0,a,dots)[0]
            li_gap = mp.sqrt(initial**2+m*m)-mp.sqrt((initial-emitted)**2+m*m)-photon
            shared_gap = k.shell(initial,m,a,dots)[0]-k.shell(initial-emitted,m,a,dots)[0]-photon
            assert li_gap > mp.mpf("7e-7")
            assert shared_gap < 0


def test_chord_cocycle_and_massive_triangle():
    with mp.workdps(70):
        support = check.direct_orbit()
        p,q = (mp.mpf(".13"),mp.mpf("-.21"),mp.mpf(".05")),(mp.mpf("-.08"),mp.mpf(".07"),mp.mpf(".12"))
        def features(v):
            return [mp.expm1(1j*sum(w[i]*v[i] for i in range(3)))/mp.sqrt(10) for w in support]
        fp,fq,fs = features(p),features(q),features(tuple(p[i]+q[i] for i in range(3)))
        for w,left,right,total in zip(support,fp,fq,fs):
            phase = mp.exp(1j*sum(w[i]*p[i] for i in range(3)))
            assert abs(total-left-phase*right) < mp.mpf("1e-65")
        norm = lambda v: mp.sqrt(sum(abs(t)**2 for t in v))
        assert norm(fs) <= norm(fp)+norm(fq)
        assert mp.sqrt(1+norm(fs)**2) < mp.sqrt(1+norm(fp)**2)+norm(fq)


@pytest.mark.parametrize("mutation", (
    lambda p: p.update(scope="empirically excluded"),
    lambda p: p.update(global_threshold_absolute_error_eV="1e-99"),
    lambda p: p["witnesses"].pop(),
    lambda p: p["witnesses"].__setitem__(1,copy.deepcopy(p["witnesses"][0])),
    lambda p: p["witnesses"][0].update(leading_soft_eV="0"),
    lambda p: p["witnesses"][0].update(collinear_soft_energy_eV="0"),
    lambda p: p["witnesses"][0].update(collinear_soft_energy_eV="NaN"),
    lambda p: p["witnesses"][0].update(collinear_soft_energy_eV="inf"),
    lambda p: p["witnesses"][0].update(outgoing_small_fraction="1/2"),
    lambda p: p["witnesses"][0].update(energy_residual_eV="1/0"),
    lambda p: p["witnesses"][0].update(energy_residual_eV="1e-1001"),
    lambda p: p["witnesses"][0].update(leading_soft_eV="1e100000000000000000000"),
    lambda p: p["witnesses"][0].update(outgoing_small_fraction="0.1"),
    lambda p: p["witnesses"][0].update(direction=[True,0,0]),
    lambda p: p["witnesses"][0].update(variant="universal_all_fields"),
    lambda p: p["witnesses"][0].update(energy_residual_eV="1"),
    lambda p: p["witnesses"][0].update(hard_momentum_eV="1e21"),
    lambda p: p["witnesses"][0].update(unreported_field=1),
))
def test_hostile_report(mutation):
    packet = check.load(REPORT)
    mutation(packet)
    with pytest.raises(ValueError):
        check.check(packet)


@pytest.mark.parametrize("value", ("nan","inf","0","1e16","1e21",True))
def test_invalid_momentum(value):
    with pytest.raises(ValueError):
        k.witness(value,"photon_only")


@pytest.mark.parametrize("precision", (True,70.0,20,121))
def test_invalid_precision(precision):
    with pytest.raises(ValueError):
        k.witness("1e19","photon_only",precision=precision)


def test_loader_rejects_duplicate_and_oversize(tmp_path):
    path = tmp_path/"bad.json"
    path.write_text('{"a":1,"a":2}',encoding="utf-8")
    with pytest.raises(ValueError,match="duplicate"):
        check.load(path)
    path.write_bytes(b" "*50001)
    with pytest.raises(ValueError,match="oversized"):
        check.load(path)


def test_optimized_cli_without_producer(tmp_path):
    # Isolate the checker: its numerical route cannot call the producer.
    import shutil
    package = tmp_path/"dispersion_feasibility"
    package.mkdir()
    for name in ("__init__.py","check.py","bounds.py"):
        shutil.copyfile(Path(__file__).with_name(name),package/name)
    environment = dict(os.environ,PYTHONPATH=str(tmp_path))
    process = subprocess.run([sys.executable,"-O","-m","dispersion_feasibility.check",str(REPORT)],
                             cwd=tmp_path,env=environment,capture_output=True,text=True,timeout=30)
    assert process.returncode == 0, process.stderr
    bad = check.load(REPORT)
    bad["witnesses"][0]["collinear_soft_energy_eV"] = "1"
    path = tmp_path/"bad.json"
    path.write_text(json.dumps(bad),encoding="utf-8")
    process = subprocess.run([sys.executable,"-O","-m","dispersion_feasibility.check",str(path)],
                             cwd=tmp_path,env=environment,capture_output=True,text=True,timeout=30)
    assert process.returncode != 0
