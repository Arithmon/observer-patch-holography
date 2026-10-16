"""Independent controls for the source-transfer/modular-code correspondence."""
from fractions import Fraction as F
from decimal import Decimal, localcontext
import copy
import math

import numpy as np
import pytest

from . import construction as c


def test_transfer_reference_is_faithful_and_source_spectrum_is_not_palm():
    c.spectrum_audit()
    p = np.array(c.reduced()['matrix'], float)
    pi = np.array(c.reduced()['stationary'], float)
    symmetric = np.diag(np.sqrt(pi))@p@np.diag(1/np.sqrt(pi))
    np.testing.assert_allclose(symmetric, symmetric.T, atol=1e-16)
    expected = np.array([1, (25+math.sqrt(3))/30, (25-math.sqrt(3))/30])
    np.testing.assert_allclose(np.linalg.eigvalsh(symmetric), sorted(expected), atol=1e-15)
    assert np.linalg.eigvalsh(symmetric)[0] > 2/3
    assert abs(np.trace(symmetric)-8/3) < 1e-15
    # The nonzero off-diagonal reference is a transfer-space object, not
    # the diagonal native Palm preparation/reference in configuration space.
    assert symmetric[0, 1] == pytest.approx(1/15)


def test_log_and_complex_phase_enclosures_cover_high_precision_independent_values():
    with localcontext() as ctx:
        ctx.prec = 70
        root = Decimal(3).sqrt()
        for j, value in enumerate((Decimal(1), (25+root)/30, (25-root)/30)):
            energy = -value.ln()
            lo, hi = c.energies()[j]
            assert Decimal(lo.numerator)/Decimal(lo.denominator) <= energy
            assert energy <= Decimal(hi.numerator)/Decimal(hi.denominator)
            assert hi-lo < F(1, 10**22)
    for s in c.SETTINGS:
        for j, (real, imag) in enumerate(c.phase_targets(s)):
            angle = 0 if j == 3 else s*float(sum(c.energies()[j])/2)
            assert abs(float(sum(real)/2)-math.cos(angle)) < 1e-15
            assert abs(float(sum(imag)/2)-math.sin(angle)) < 1e-15


def test_raw_source_comparison_includes_padding_and_unrelated_mode():
    raw = c.compiler.execute(c.diagonal_word(4), 2)
    assert c.certify(raw, 'phase', s=4) < c.ERROR_GATE
    for index in (2, 3):
        bad = raw.copy(); bad[index, 0] += 2**-14
        assert c.certify(bad, 'phase', s=4) > F(1, 10**5)
    assert c.certify(raw.conj(), 'phase', s=4) > F(1, 10)


@pytest.mark.parametrize('j', [1, 2])
def test_both_quadratures_read_the_derived_modular_phase(j):
    energy = float(sum(c.energies()[j])/2)
    for axis, function in (('x', math.cos), ('y', math.sin)):
        result = c.factor_execute(c.program(j, (8,), quadrature=axis))
        expected = (1+function(8*energy))/2
        assert len(result.branches) == 4
        assert result.branches[0].probability == pytest.approx(expected, abs=5e-13)
        assert result.branches[j].probability == pytest.approx(1-expected, abs=5e-13)
        assert all(b.probability >= 0 for b in result.branches)
        assert result.branches[3].probability > 0  # Raw padding leakage retained.


def test_continuing_state_is_not_reprepared_and_checkpoint_is_bound_to_prefix():
    p = c.program(2, (4, 1)); result = c.factor_execute(p)
    prefix = c.factor_execute(p, stop_after=5)
    restored = c.factor_execute(p, checkpoint=copy.deepcopy(prefix.checkpoint))
    for a, b in zip(result.branches, restored.branches, strict=True):
        np.testing.assert_array_equal(a.density, b.density)
        assert a.resources == b.resources
    second_plus = sum(b.probability for b in result.branches if b.record[-1] == 0)
    fresh_plus = c.factor_execute(c.program(2, (1,))).branches[0].probability
    assert abs(second_plus-fresh_plus) > .1
    changed = copy.deepcopy(prefix.checkpoint)
    changed['branches'][0][1][0][0] += .001
    with pytest.raises(ValueError, match='complete source prefix'):
        c.factor_execute(p, checkpoint=changed)


def test_generic_density_runtime_agrees_on_complete_unconditioned_terminal_state():
    p = c.program(1, (4, 1)); factors = c.factor_execute(p)
    events = []
    for event in p.events:
        original = event.instruments[0]
        instrument = c.runtime.Instrument(kraus=(tuple(k for group in original.kraus for k in group),),
                                         resources=original.resources, binding={'all_outputs_forgotten': True})
        events.append(c.runtime.Event(event.name, instruments=(instrument,)))
    unrecorded = c.runtime.Program('independent_density_check', 4, p.initial_density, tuple(events))
    dense = c.runtime.execute(unrecorded).branches[0].density
    np.testing.assert_allclose(sum(b.density for b in factors.branches), dense, atol=1e-14, rtol=0)


def test_raw_map_gain_is_charged_and_resources_do_not_add_peak_widths():
    eps = c.ERROR_GATE
    assert c._sequence_allowance(7) > 7*eps
    assert c._sequence_allowance(7) < F(1, 10**10)
    r = c.factor_execute(c.program(1, (4, 1)))
    assert r.resources['four_outcome_reads'] == 2
    assert r.resources['native_pulses'] == 430
    assert 'maximum_work_qubits' not in r.resources


@pytest.mark.parametrize('bad', [True, -1, 2, 4., '4', None])
def test_modular_parameter_domain_is_explicit(bad):
    with pytest.raises(ValueError):
        c.diagonal_word(bad)


def test_bad_source_arrays_are_rejected_before_subclass_dispatch():
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'copy'):
                raise AssertionError('untrusted array dispatch')
            return super().__getattribute__(name)
    with pytest.raises(ValueError, match='base ndarray'):
        c.certify(np.eye(4).view(Hostile), 'phase')
    with pytest.raises(ValueError):
        c.certify(np.full((4, 4), np.nan), 'phase')
