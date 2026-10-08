"""Missingness must survive original nested integer containers."""
from fractions import Fraction
import importlib.util
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("native_original_input_producer", HERE / "build.py")
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)


@pytest.mark.parametrize("container", [list, tuple])
@pytest.mark.parametrize("hidden", [999, -(2**62)])
@pytest.mark.parametrize("consumer", ["matrix", "left_covariance", "right_covariance"])
def test_nested_masked_rows_cannot_supply_a_hidden_sample(container, hidden, consumer):
    # At a04620a7 the outer container passes the mask guard; NumPy then
    # strips the inner mask. For hidden=999, self covariance would be 249001.
    missing = container((np.array([1]), np.ma.array([hidden], mask=[True])))
    observed = [[1], [3]]
    with pytest.raises(ValueError, match="masked"):
        if consumer == "matrix":
            producer.matrix(missing)
        elif consumer == "left_covariance":
            producer.centered_cov(missing, observed)
        else:
            producer.centered_cov(observed, missing)


@pytest.mark.parametrize("placement", ["outer", "row", "scalar", "object_scalar"])
def test_missingness_is_rejected_before_any_array_conversion(placement):
    if placement == "outer":
        value = np.ma.array([[1], [2]], mask=[[False], [True]])
    elif placement == "row":
        value = [np.ma.array([1], mask=False), np.ma.array([2], mask=True)]
    elif placement == "scalar":
        value = [[np.ma.masked], [2]]
    else:
        value = np.empty((2, 1), dtype=object)
        value[0, 0], value[1, 0] = np.ma.array(1, mask=True), 2
    with pytest.raises(ValueError):
        producer.matrix(value)


@pytest.mark.parametrize("kind", ["signed_rows", "unsigned_rows", "object_rows", "plain_tuple"])
def test_complete_nested_integer_populations_preserve_their_exact_variance(kind):
    if kind == "signed_rows":
        values = [np.array([2**62+1], dtype=np.int64), np.array([2**62+3], dtype=np.int64)]
    elif kind == "unsigned_rows":
        values = [np.array([2**64-1], dtype=np.uint64), np.array([2**64-3], dtype=np.uint64)]
    elif kind == "object_rows":
        values = [np.array([2**150+1], dtype=object), np.array([2**150+3], dtype=object)]
    else:
        values = ((2**150+1,), (2**150+3,))
    assert producer.centered_cov(values, values).tolist() == [[Fraction(1)]]
    assert producer.matrix(values).tolist() == [[Fraction(int(row[0]))] for row in values]
