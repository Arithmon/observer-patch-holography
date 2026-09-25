from fractions import Fraction as F
import importlib.util
from pathlib import Path
import pytest
p=Path(__file__).with_name('certify.py');s=importlib.util.spec_from_file_location('common_filter_certify',p)
c=importlib.util.module_from_spec(s);s.loader.exec_module(c)

@pytest.mark.parametrize('coefficients,rates',[([1],[1]),([1,2],[1,4]),([1,-4],[1,4]),([2,-5,1],[1,2,9])])
def test_fourier_and_time_domain_covariance_agree(coefficients,rates):
    vals=[]
    for lam in [.01,.1,1,10,100]:
        exact=c.variance(lam,coefficients,rates)
        assert float(exact)==pytest.approx(c.fourier_variance(lam,coefficients,rates),rel=2e-9,abs=1e-11)
        vals.append(F(str(lam))*exact)
    assert all(b>=a for a,b in zip(vals,vals[1:]))


def test_required_red_shape_violates_monotone_constraint():
    # k^-3 means lambda^-3/2, with lambda=k².
    lam1,lam2=F(1),F(4)
    target1,target2=F(1),F(1,8)
    assert lam2*target2 < lam1*target1


def test_leak_formula_and_signed_memory_dc():
    assert c.variance(4,[1],[3])==F(1,15)
    assert sum(c/g for c,g in zip([F(1),F(-4)],[F(1),F(4)]))==0
    assert c.variance(4,[1,-4],[1,4])>0


def test_invalid_rates_rejected():
    with pytest.raises(ValueError):c.variance(0,[1],[1])
    with pytest.raises(ValueError):c.variance(1,[1],[-1])
