"""Original-input controls for the official posterior accounting."""
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import math
from pathlib import Path

import pytest

from test_official_desi_dr2_chain_audit import mod


def summarize(rows, cls=mod.Accumulator):
    acc = cls()
    for row in rows:
        acc.add(*row)
    return acc.summary()


def decimal_value(value):
    if isinstance(value, Fraction):
        return Decimal(value.numerator) / Decimal(value.denominator)
    return Decimal(value)


def pairwise_control(rows):
    """Independent centered pair-difference identity, from original scalars."""
    with localcontext() as ctx:
        ctx.prec = 180
        rows = [tuple(map(decimal_value, row)) for row in rows]
        weight = sum(row[0] for row in rows)
        mean = [sum(w * row[k] for w, *row in rows) / weight for k in range(2)]
        cov = [[sum(a[0] * b[0] * (a[i+1]-b[i+1]) * (a[j+1]-b[j+1])
                    for a, b in itertools.combinations(rows, 2)) / weight**2
                for j in range(2)] for i in range(2)]
        return dict(w0_mean=float(mean[0]), wa_mean=float(mean[1]),
                    w0_std=float(cov[0][0].sqrt()), wa_std=float(cov[1][1].sqrt()),
                    w0_wa_covariance=float(cov[0][1]),
                    w0_wa_correlation=float(cov[0][1] / (cov[0][0]*cov[1][1]).sqrt())
                    if cov[0][0] and cov[1][1] else None,
                    weight_concentration_ess_not_autocorrelation_corrected=
                    float(weight**2 / sum(r[0]**2 for r in rows)))


@pytest.mark.parametrize('scale', [1., 2.**600, 2.**-600])
@pytest.mark.parametrize('offset,delta', [(0., 1.), (-1., 2.**-30), (2.**50, .5)])
def test_weighted_components_match_pairwise_original_input_control(scale, offset, delta):
    rows = [(scale, offset-delta, -delta), (2*scale, offset, 2*delta),
            (3*scale, offset+delta, -delta)]
    observed = summarize(rows)
    expected = pairwise_control(rows)
    for key, value in expected.items():
        assert observed[key] == pytest.approx(value, rel=1e-12, abs=0), key


@pytest.mark.parametrize('wa,expected', [(-2.**-60, 0.), (0., 1.), (2.**-60, 1.)])
def test_cpl_boundary_keeps_both_sides(wa, expected):
    assert summarize([(1., -1., wa)])['posterior_mass_w_ge_minus_one_for_0_le_z_le_2'] == expected


def test_chain_decimal_boundary_is_not_rounded_to_binary64(tmp_path):
    path = tmp_path / 'decimal.txt'
    path.write_text('# weight w wa\n1 -1.0000000000000000001 0\n')
    result = mod.read_chain(path).summary()
    assert result['raw_rows_in_monotone_subset'] == 0


def test_positive_complement_is_not_subtracted_from_rounded_one():
    result = summarize([(2.**60, -1., 0.), (1., -2., 0.)])
    assert result['posterior_mass_capacity_loss_somewhere_for_0_le_z_le_2'] == pytest.approx(
        float(Fraction(1, 2**60 + 1)), rel=1e-15, abs=0)
    assert result['w0_std'] > 0


def test_weighted_quantile_preserves_tiny_last_atom():
    assert mod.weighted_quantiles([(0., 2.**60), (1., 1.)], (1.,)) == {'1.000': 1.}


def test_weighted_quantile_exact_half_threshold():
    assert mod.weighted_quantiles([(0., 2.**60), (1., 1.), (2., 2.**60)], (.5,)) == {'0.500': 1.}


@pytest.mark.parametrize('shift', [0., 2.**40])
def test_permutation_and_merge_preserve_all_accounting(shift):
    rows = [(2.**60, shift-1., -2.**-30), (1., shift+1., 2.**-30), (3., shift, 0.)]
    expected = summarize(rows)
    for perm in itertools.permutations(rows):
        acc = mod.Accumulator()
        for row in perm:
            part = mod.Accumulator()
            part.add(*row)
            acc.merge(part)
        assert acc.summary() == expected
        assert summarize(perm) == expected


def test_base_lcdm_narrow_posterior_preserves_spread():
    rows = [(1., 68.-2.**-25, .7), (1., 68.+2.**-25, .7)]
    result = summarize(rows, mod.BaseLCDMAccumulator)
    assert result['H0_km_s_Mpc']['weighted_std'] == 2.**-25
    assert result['OmegaLambda']['weighted_std'] == 0
    with localcontext() as ctx:
        ctx.prec = 180
        # Independent SI formula and two-point spread, no producer conversion.
        values = [3*Decimal(o)*(Decimal(h)*1000/Decimal('3.0856775814913673e22') /
                  Decimal(299792458))**2*Decimal('1.616255e-35')**2 for _, h, o in rows]
        expected = float(abs(values[1]-values[0])/2)
    assert result['Lambda_lP2']['weighted_std'] == pytest.approx(expected, rel=1e-12, abs=0)


def fixture_audit(tmp_path, rows):
    path = tmp_path / 'toy_chain.1.txt'
    path.write_text('# weight w wa\n' + ''.join(f'{w} {x} {y}\n' for w,x,y in rows))
    return mod.audit_dataset(tmp_path, dict(slug='toy', model='toy', directory='toy',
                                          sha256=[mod.sha256(path)]))


@pytest.mark.parametrize('center', [0, 10, 40])
def test_gaussian_tail_survives_through_public_dataset_reader(tmp_path, center):
    # Independent four-corner distribution: covariance I, displacement (0,center).
    result = fixture_audit(tmp_path, [(1, -1+x, center+y) for x in [-1,1] for y in [-1,1]])
    diagnostic = result['fixed_capacity_point_gaussian_diagnostic']
    assert diagnostic['mahalanobis_squared'] == center**2
    if center <= 10:
        assert diagnostic['chi2_2dof_survival'] == pytest.approx(math.exp(-center**2/2), rel=1e-12, abs=0)
    else:
        assert diagnostic['chi2_2dof_survival'] is None
        assert diagnostic['log_chi2_2dof_survival'] == -800.
    # Check inverse via a different implementation at sufficient precision.
    import mpmath
    ctx = mpmath.mp.clone(); ctx.dps = 80
    sigma = diagnostic['two_sided_normal_sigma_equivalent']
    assert float(ctx.log(ctx.erfc(ctx.mpf(sigma)/ctx.sqrt(2)))) == pytest.approx(-center**2/2, rel=1e-12, abs=1e-15)


def test_singular_gaussian_keeps_valid_posterior_evidence(tmp_path):
    result = fixture_audit(tmp_path, [(1, -1, 0), (1, 0, 1)])
    assert result['combined']['raw_rows'] == 2
    diagnostic = result['fixed_capacity_point_gaussian_diagnostic']
    assert diagnostic['status'] == 'unavailable_singular_covariance'
    assert diagnostic['mahalanobis_squared'] is None


@pytest.mark.parametrize('header', ['weight w wa w', 'weight w wa weight'])
def test_ambiguous_chain_columns_are_rejected(tmp_path, header):
    path = tmp_path / 'duplicate.txt'
    path.write_text(f'# {header}\n1 -2 0 1\n')
    with pytest.raises(ValueError, match='duplicate'):
        mod.read_chain(path)


@pytest.mark.parametrize('bad', [True, float('nan'), float('inf'), 1j])
def test_invalid_scalar_does_not_mutate_accumulator(bad):
    acc = mod.Accumulator(); acc.add(1., -1., 0.)
    before = acc.summary()
    with pytest.raises(ValueError):
        acc.add(1., bad, 0.)
    assert acc.summary() == before


def test_unrepresentable_positive_summary_field_is_explicitly_refused():
    acc = mod.Accumulator()
    acc.add(1., -1., 0.)
    acc.add(2.**-1074, -2., 0.)
    # The nonzero excluded mass is resolved, but covariance can vanish exactly;
    # with a second nonzero coordinate the nonzero covariance cannot be encoded.
    acc = mod.Accumulator()
    acc.add(1., -1., 0.)
    acc.add(2.**-1074, -1.5, 1.)
    with pytest.raises(ValueError, match='binary64'):
        acc.summary()
