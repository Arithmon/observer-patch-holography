"""Derived-statistic corruption must be rejected without a large graph replay."""
import copy
import json
from pathlib import Path
import pytest
from oph_exact.verify_source_net_manifold_sampled_independent import check_region

HERE=Path(__file__).resolve().parents[2]
DATA=json.loads((HERE/'source_net_manifold_sampled_q55_q89_2026-09-25.json').read_text())

def failures(row):
    result=[]
    check_region(row,3,DATA['continuum_references'],DATA['spectrum_grid'],result,'mutation control')
    return result

def test_original_derived_statistics_pass():
    assert not failures(DATA['levels'][0]['families'][0]['regions'][0])

@pytest.mark.parametrize('mutation',['normalization','dimension','cdf'])
def test_changed_derived_statistics_fail(mutation):
    row=copy.deepcopy(DATA['levels'][0]['families'][0]['regions'][0])
    if mutation=='normalization':row['chains']['3']['count_over_falling_factorial']*=2
    if mutation=='dimension':row['chains']['3']['inverted_dimension']=9
    if mutation=='cdf':row['interval_spectrum']['cdf_on_grid'][0]=1
    assert failures(row)
