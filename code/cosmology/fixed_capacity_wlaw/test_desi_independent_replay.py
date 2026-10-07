"""Challenge independent replay using small non-degenerate chain fixtures."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

from test_official_desi_dr2_chain_audit import mod

PATH = Path(__file__).with_name('verify_official_desi_dr2_independent.py')
SPEC = importlib.util.spec_from_file_location('desi_independent', PATH)
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


@pytest.fixture
def packet(tmp_path, monkeypatch):
    pins = {}
    for name, spec in mod.DOWNLOAD_SPECS.items():
        base = spec is mod.BASE_LCDM_DATASET
        group = []
        hashes = []
        for i in range(1,5):
            filename = f"{spec['slug']}_chain.{i}.txt"
            path = tmp_path / filename
            header = '# weight H0 omegal\n' if base else '# weight w wa\n'
            rows = [(i,60,'.6'), (2,60,'.8'), (3,80,'.6'), (4,80,'.8')] if base else [
                (i,-2,-1), (2,-2,1), (3,0,-1), (4,0,1)]
            path.write_text(header+''.join(f'{w} {x} {y}\n' for w,x,y in rows))
            digest = mod.sha256(path)
            hashes.append(digest); group.append((filename,digest))
        monkeypatch.setitem(spec, 'sha256', hashes)
        pins['base_lcdm_capacity_display' if base else name] = group
    monkeypatch.setattr(verifier, 'PINS', pins)
    # JSON round-trip exercises exactly the serialized public interface.
    receipt = json.loads(json.dumps(mod.build_receipt(tmp_path), allow_nan=False))
    return receipt, tmp_path


def test_independent_replay_accepts_complete_valid_packet_without_producer(packet, monkeypatch):
    receipt, data_dir = packet
    monkeypatch.setitem(sys.modules, 'official_desi_dr2_chain_audit', None)
    assert verifier.replay(receipt, data_dir) > 250


@pytest.mark.parametrize('field', ['raw_rows','expanded_posterior_weight',
    'weight_concentration_ess_not_autocorrelation_corrected','w0_mean','w0_std','wa_mean','wa_std',
    'w0_wa_covariance','w0_wa_correlation','raw_rows_in_monotone_subset',
    'posterior_mass_w_ge_minus_one_for_0_le_z_le_2',
    'posterior_mass_capacity_loss_somewhere_for_0_le_z_le_2',
    'posterior_mass_w0_gt_minus_one','posterior_mass_wa_nonnegative'])
def test_verifier_rejects_each_cpl_field_mutation(packet, field):
    receipt, data_dir = packet
    summary = receipt['datasets']['DESI_DR2_BAO+CMB']['combined']
    summary[field] += 1
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


@pytest.mark.parametrize('value', [0, -1e-300, None, True, float('nan'), float('inf')])
def test_verifier_cannot_hide_small_lambda_component_behind_absolute_floor(packet, value):
    receipt, data_dir = packet
    receipt['base_lcdm_capacity_display']['combined']['Lambda_lP2']['weighted_std'] = value
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


@pytest.mark.parametrize('field', ['mahalanobis_squared','chi2_2dof_survival',
                                   'log_chi2_2dof_survival','two_sided_normal_sigma_equivalent'])
def test_verifier_rejects_each_gaussian_field_mutation(packet, field):
    receipt, data_dir = packet
    receipt['datasets']['DESI_DR2_BAO+CMB']['fixed_capacity_point_gaussian_diagnostic'][field] = 0
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


@pytest.mark.parametrize('mutation', ['missing_dataset','missing_chain','missing_field','extra_field',
    'changed_input','changed_input_and_pin','wrong_range','wrong_tail','unavailable','bool_count','float_count','confirmation','units'])
def test_verifier_rejects_incomplete_changed_or_mislabelled_evidence(packet, mutation):
    receipt, data_dir = packet
    dataset = receipt['datasets']['DESI_DR2_BAO+CMB']
    if mutation == 'missing_dataset':
        receipt['datasets'].pop('DESI_DR2_BAO+CMB')
    elif mutation == 'missing_chain':
        dataset['chains'].pop()
    elif mutation == 'missing_field':
        dataset['combined'].pop('w0_std')
    elif mutation == 'extra_field':
        dataset['combined']['extra'] = 1
    elif mutation in ('changed_input','changed_input_and_pin'):
        path = data_dir/dataset['chains'][0]['file']
        path.write_text(path.read_text().replace('-2','-3'))
        if mutation == 'changed_input_and_pin':
            dataset['chains'][0]['sha256'] = mod.sha256(path)
    elif mutation == 'wrong_range':
        dataset['chain_range_for_monotone_subset_mass'] = [0,1]
    elif mutation == 'wrong_tail':
        dataset['rare_tail_resolution']['classification'] = 'resolved_in_no_chains'
    elif mutation == 'unavailable':
        dataset['fixed_capacity_point_gaussian_diagnostic']['status'] = 'unavailable'
    elif mutation == 'float_count':
        dataset['combined']['raw_rows'] = float(dataset['combined']['raw_rows'])
    elif mutation == 'confirmation':
        receipt['epistemic_status']['oph_confirmation'] = True
    elif mutation == 'units':
        receipt['base_lcdm_capacity_display']['constants']['Mpc_in_m'] = 1.
    else:
        dataset['combined']['raw_rows'] = True
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)

@pytest.mark.parametrize('path,value', [
    (('subset_definition','redshift_range'), '0 <= z <= 200'),
    (('subset_definition','monotone_capacity_condition'), 'w(a)>=0 at every epoch'),
    (('base_lcdm_capacity_display','sample_level_formula'), 'Lambda*l_P^2 = omegal*H0'),
    (('base_lcdm_capacity_display','model_scope'), 'model-independent physical measurement'),
    (('base_lcdm_capacity_display','classification'), 'OPH confirmed'),
    (('datasets','DESI_DR2_BAO+CMB','source_directory'), 'https://example.invalid/fiction'),
    (('datasets','DESI_DR2_BAO+CMB','rare_tail_resolution','warning'), 'Frequentist OPH exclusion'),
    (('arithmetic','chain_columns'), 'rounded to integers'),
    (('producer','script_sha256'), '0'*64),
    (('source','official_sha256_manifest_sha256'), '0'*64),
    (('epistemic_status','notes'), 'All physical premises have been derived'),
])
def test_verifier_binds_scientific_context_as_well_as_numbers(packet, path, value):
    receipt, data_dir = packet
    target = receipt
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


@pytest.mark.parametrize('path', [(), ('producer',), ('source',), ('arithmetic',),
    ('subset_definition',), ('base_lcdm_capacity_display',),
    ('datasets','DESI_DR2_BAO+CMB'),
    ('datasets','DESI_DR2_BAO+CMB','rare_tail_resolution'),
    ('datasets','DESI_DR2_BAO+CMB','fixed_capacity_point_gaussian_diagnostic')])
def test_verifier_rejects_unverified_extra_contract_fields(packet, path):
    receipt, data_dir = packet
    target = receipt
    for key in path:
        target = target[key]
    target['frozen_prediction_confirmed'] = True
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


def test_chain_identity_must_be_an_integer_not_boolean(packet):
    receipt, data_dir = packet
    receipt['datasets']['DESI_DR2_BAO+CMB']['chains'][0]['chain'] = True
    with pytest.raises(ValueError):
        verifier.replay(receipt, data_dir)


@pytest.mark.parametrize('duplicate', ['schema', 'datasets', 'producer'])
def test_public_cli_rejects_duplicate_json_keys_even_if_values_agree(packet, monkeypatch, duplicate):
    receipt, data_dir = packet
    path = data_dir / 'duplicate.json'
    payload = json.dumps(receipt)
    payload = payload[:-1] + ',' + json.dumps(duplicate) + ':' + json.dumps(receipt[duplicate]) + '}'
    path.write_text(payload, encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['verify', '--data-dir', str(data_dir), '--receipt', str(path)])
    with pytest.raises(ValueError, match='duplicate JSON key'):
        verifier.main()


def test_public_cli_accepts_complete_packet_without_importing_producer(packet, monkeypatch, capsys):
    receipt, data_dir = packet
    path = data_dir / 'receipt.json'
    path.write_text(json.dumps(receipt), encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['verify', '--data-dir', str(data_dir), '--receipt', str(path)])
    monkeypatch.setitem(sys.modules, 'official_desi_dr2_chain_audit', None)
    monkeypatch.setitem(sys.modules, 'official_desi', None)
    verifier.main()
    assert capsys.readouterr().out.startswith('PASS: 20 pinned chains, 25 summaries')


@pytest.mark.parametrize('constant', ['NaN', 'Infinity', '-Infinity'])
def test_public_cli_refuses_nonstandard_json_numbers(packet, monkeypatch, constant):
    receipt, data_dir = packet
    path = data_dir / 'nonfinite.json'
    payload = json.dumps(receipt)
    payload = payload[:-1] + ',"nonfinite":' + constant + '}'
    path.write_text(payload, encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['verify', '--data-dir', str(data_dir), '--receipt', str(path)])
    with pytest.raises(ValueError, match='non-finite JSON constant'):
        verifier.main()
