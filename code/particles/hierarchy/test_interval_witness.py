"""Regression tests for the independent hierarchy interval witness.

The witness must certify the bracket and the negative derivative enclosure on
the declared interval, and every antecedent-perturbation control must fail
closed. The stored receipt must agree with the checked-in claims.
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp
from mpmath.libmp import to_rational
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "computations"))

import hierarchy_interval_witness as w  # noqa: E402

iv = w.iv


def endpoints(text):
    return tuple(Fraction(x.strip()) for x in text.strip("[]").split(","))


def exact_binary_endpoints(interval):
    return tuple(Fraction(*to_rational(x)) for x in interval._mpi_)


@pytest.fixture(scope="module")
def replay():
    # Complete declared 64-iteration formula and all eight derivative cells.
    # Low scalar precision must not affect the interval or exported bounds.
    previous = iv.dps
    try:
        iv.dps = 8
        with mp.workdps(8):
            result = w.produce()
            assert iv.dps == 8
            return result
    finally:
        iv.dps = previous


def test_bracket_sign_change_fast():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    f_lo = w.phi(iv.mpf(w.I_LO), w.P_STR, b1, b2, b3)
    f_hi = w.phi(iv.mpf(w.I_HI), w.P_STR, b1, b2, b3)
    assert w.sign_of(f_lo) == 1
    assert w.sign_of(f_hi) == -1


def test_derivative_enclosure_negative_one_piece():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    _, cells = w.subdivision(2)
    res = w.phi(w.Dual.var(cells[0]), w.P_STR, b1, b2, b3)
    assert res.d.b < 0
    assert mp.mpf(res.d.a) > -12 and mp.mpf(res.d.b) < -10


def test_perturbed_beta_breaks_bracket():
    b1, b2, b3 = iv.mpf(33) / 5, iv.mpf(1), iv.mpf(-3)
    f_lo = w.phi(iv.mpf(w.I_LO), w.P_STR, b1, b2 + iv.mpf('0.05'), b3)
    f_hi = w.phi(iv.mpf(w.I_HI), w.P_STR, b1, b2 + iv.mpf('0.05'), b3)
    assert w.sign_of(f_lo) * w.sign_of(f_hi) != -1


def test_receipt_matches_claims(replay):
    receipt = json.loads((HERE / "certificates" / "independent_interval_witness_receipt.json").read_text(encoding="utf-8"))
    assert receipt == replay
    assert receipt["bracket_sign_change"] is True
    assert receipt["derivative_strictly_negative"] is True
    assert receipt["unique_root_certified"] is True
    assert receipt["controls_all_fail_closed"] is True
    lo, hi = (mp.mpf(x) for x in receipt["derivative_enclosure_union"].strip("[]").split(","))
    assert lo < mp.mpf("-10.995768") and hi > mp.mpf("-10.985284"), (
        "independent enclosure must contain the published enclosure")
    assert hi < 0


@pytest.mark.parametrize("value", ["0.1", "-0.1", "1e-90", "-1e90"])
def test_serialization_contains_exact_binary_endpoints_at_low_scalar_precision(value):
    x = iv.mpf(value)
    exact_lo, exact_hi = exact_binary_endpoints(x)
    with mp.workdps(5):
        lo, hi = endpoints(w.interval_str(x))
    assert lo <= exact_lo <= exact_hi <= hi
    assert lo <= Fraction(value) <= hi


def test_pi_is_enclosed_and_the_legacy_serializer_is_falsified():
    x = iv.pi
    exact_lo, exact_hi = exact_binary_endpoints(x)
    with mp.workdps(15):
        legacy = (mp.nstr(mp.mpf(x.a), 25), mp.nstr(mp.mpf(x.b), 25))
        lo, hi = endpoints(w.interval_str(x))
    assert lo <= exact_lo < exact_hi <= hi
    assert Fraction(legacy[1]) < exact_lo


@pytest.mark.parametrize("count", [1, 2, 3, 8, 17])
def test_exact_partition_covers_the_entire_requested_interval(count):
    with mp.workdps(5):
        knots, cells = w.subdivision(count)
    assert knots[0] == Fraction(w.I_LO) and knots[-1] == Fraction(w.I_HI)
    assert len(knots) == count + 1 and len(cells) == count
    for i, cell in enumerate(cells):
        lo, hi = exact_binary_endpoints(cell)
        assert lo <= knots[i] < knots[i + 1] <= hi
        if i:
            assert exact_binary_endpoints(cells[i - 1])[1] >= lo


def test_legacy_partition_omits_the_lower_strip_but_exact_partition_covers_it():
    with mp.workdps(15):
        legacy_start = Fraction(str(mp.mpf(w.I_LO)))
        _, cells = w.subdivision()
    assert legacy_start - Fraction(w.I_LO) == Fraction(3, 500_000_000_000_000_000)
    assert exact_binary_endpoints(cells[0])[0] <= Fraction(w.I_LO) < legacy_start


@pytest.mark.parametrize("count", [0, -1, 8.0, True])
def test_invalid_partition_count_fails_closed(count):
    with pytest.raises(ValueError):
        w.subdivision(count)


@pytest.mark.parametrize("value", [iv.mpf("inf"), iv.mpf("-inf"), iv.mpf("nan")])
def test_nonfinite_interval_cannot_be_exported_as_a_certificate(value):
    with pytest.raises((ValueError, OverflowError)):
        w.interval_str(value)


def test_literal_derivative_cells_and_union_certify_the_claimed_domain(replay):
    part = replay["derivative_subdivision"]
    knots = list(map(Fraction, part["exact_rational_knots"]))
    assert knots[0] == Fraction(w.I_LO) and knots[-1] == Fraction(w.I_HI)
    union_lo, union_hi = endpoints(replay["derivative_enclosure_union"])
    assert union_hi < 0
    for i, (cell, derivative) in enumerate(zip(part["outward_cells"], part["derivative_enclosures"], strict=True)):
        lo, hi = endpoints(cell)
        dlo, dhi = endpoints(derivative)
        assert lo <= knots[i] < knots[i + 1] <= hi
        assert union_lo <= dlo <= dhi <= union_hi < 0
    assert replay["formula_scope"]["mu_iterations"] == 64
    assert replay["formula_scope"]["infinite_sum_or_fixed_point_limit_certified_here"] is False
    assert replay["historical_correction"]["prior_receipt_sha256"] == "ae6262c86b592796b1574f3c58ae4ca44df1516fa169c2b1c6ef6f735a5b2291"


def test_interval_precision_restored_when_production_fails(monkeypatch):
    def failure():
        assert iv.dps == 60
        raise RuntimeError("controlled failure")
    monkeypatch.setattr(w, "_produce_at_working_precision", failure)
    previous = iv.dps
    try:
        iv.dps = 17
        with pytest.raises(RuntimeError, match="controlled failure"):
            w.produce()
        assert iv.dps == 17
    finally:
        iv.dps = previous


@pytest.fixture(scope='module')
def implicit_verifier():
    from types import ModuleType
    path=HERE/'validators/verify_implicit_scale_attachment.py'
    module=ModuleType('independent_implicit_scale_verifier')
    module.__file__=str(path)
    exec(compile(path.read_bytes(),str(path),'exec'),module.__dict__)
    return module


def test_true_implicit_scale_has_independent_fresh_certificate(implicit_verifier):
    previous=iv.dps
    try:
        iv.dps=8
        result=implicit_verifier.verify(implicit_verifier.load())
        assert iv.dps==8
    finally:
        iv.dps=previous
    assert result['true_implicit_scale'] is True
    assert result['unique_hierarchy_root'] is True
    assert result['infinite_representation_sums'] is False


@pytest.mark.parametrize('mutation',['physical','infinite','pixel','count','contraction','pin','pin_float','keys'])
def test_implicit_scope_and_antecedents_fail_closed(implicit_verifier,mutation):
    r=implicit_verifier.load()
    if mutation=='physical': r['scope']['physical_hierarchy_attachment']=True
    elif mutation=='infinite': r['scope']['infinite_representation_sums']=True
    elif mutation=='pixel': r['pixel_decimal']='1.730968209403959324879279847782648941'
    elif mutation=='count': r['finite_representation_cutoffs']['su2_max_n']=128.0
    elif mutation=='contraction': r['contraction_bound']='1'
    elif mutation=='pin': r['source_pins'][next(iter(r['source_pins']))]['sha256']='0'*64
    elif mutation=='pin_float':
        entry=r['source_pins'][next(iter(r['source_pins']))]
        entry['bytes']=float(entry['bytes'])
    else: r['ignored_claim']=True
    with pytest.raises(ValueError): implicit_verifier.verify(r)


@pytest.mark.parametrize('mutation',['narrow_image','fake_derivative'])
def test_implicit_numerical_falsifiers(implicit_verifier,mutation):
    r=implicit_verifier.load()
    if mutation=='narrow_image': r['scale_image_outer']=['7.5e-18','7.5e-18']
    else: r['implicit_residual_derivative_outer']=['-11.001','-10.999']
    with pytest.raises(ValueError): implicit_verifier.verify(r)


def test_implicit_duplicate_keys_rejected(implicit_verifier,tmp_path):
    p=tmp_path/'duplicate.json';p.write_text('{"schema":1,"schema":2}',encoding='utf-8')
    with pytest.raises(ValueError,match='duplicate'): implicit_verifier.load(p)


def test_implicit_scale_global_rectangle_has_positive_denominators(implicit_verifier):
    a=iv.mpf(implicit_verifier.A_DEC);y=iv.mpf(['7e-18','8e-18'])
    image,dy,da,_=implicit_verifier.ingredients(a,y)
    assert y.a<image.a and image.b<y.b
    assert dy.a>0 and dy.b<iv.mpf('0.01') and da.a>0
    # A much lower positive scale hits a forbidden strong-running denominator.
    with pytest.raises(ValueError,match='denominator'):
        implicit_verifier.ingredients(a,iv.mpf('1e-50'))


@pytest.mark.parametrize('token',['1.0','1e9999','NaN'])
def test_implicit_loader_rejects_undeclared_numeric_tokens(implicit_verifier,tmp_path,token):
    p=tmp_path/'number.json';p.write_text('{"value":'+token+'}',encoding='utf-8')
    with pytest.raises(ValueError): implicit_verifier.load(p)
