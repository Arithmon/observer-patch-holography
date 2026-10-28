"""Exact join of source/bit stationarity to the actual C/E FIFO process."""
from fractions import Fraction as F

import pytest

from . import observer, observer_reference as ref, stationary_reference
from tools import strict_json


@pytest.fixture(scope="module")
def stationary_joint():
    return stationary_reference.run()["stationary_joint"]


@pytest.mark.parametrize("c", (F(0), F(1)))
def test_full_record_lift_is_normalized_stationary_and_has_right_projection(stationary_joint, c):
    blocks = ref.lift(stationary_joint, c)
    expected = {(tuple(row["state"]), row["memory"]): F(row["probability"])
                for row in stationary_joint if F(row["probability"])}
    assert ref.source_bit_projection(blocks) == expected
    assert observer.trace(observer.summed_density(blocks)) == 1
    assert len(blocks) > 2  # Representatives were actually replaced by full records.
    certificate = ref.certify(blocks, c)
    assert certificate["stationarity_residual"] == "0"
    assert certificate["faithful_on_declared_direct_sum_of_positive_supports"]
    assert not certificate["faithful_on_full_M78_for_each_record"]
    assert all(memory[0].law == ("C" if c else "E") for memory in blocks)


def test_bit_representatives_are_not_mislabeled_as_full_stationary_state(stationary_joint):
    seed = ref.seed_blocks(stationary_joint)
    assert len(seed) == 2
    assert seed != ref.lift(stationary_joint, F(1))


def test_certificate_rejects_nonstationary_and_unnormalized_inputs(stationary_joint):
    seed = ref.seed_blocks(stationary_joint)
    with pytest.raises(ArithmeticError, match="not stationary"):
        ref.certify(seed, F(1))
    doubled = {memory: {key: 2*value for key, value in density.items()}
               for memory, density in seed.items()}
    with pytest.raises(ArithmeticError, match="not normalized"):
        ref.certify(doubled, F(1))


def test_full_record_reference_receipt_replays():
    assert ref.report() == strict_json.load(ref.HERE/"observer_reference_results.json")
