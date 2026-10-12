"""Replay observer Hilbert, modular-clock, region and occupation attachments."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import experiment, operational, modular_clock, refinement, spatial, history_repair, slow_geometry
from .scalar_coupling import ideal_hamiltonian


def regional_receipt(history):
    d = history['dimension']
    vectors = operational.terminal_vectors(history)
    work = sum(np.outer(v, v.conj()) for v in vectors.values())
    regions = []
    for name in spatial.REGIONS:
        rho = spatial.regional_density(d, work, name)
        regions.append({'region': name, 'matrix_dimension': len(rho),
                        'state_trace': float(np.trace(rho).real),
                        'purity': float(np.trace(rho@rho).real)})
    rho0, _ = operational.scalar_state(history)
    marginal_error = float(np.linalg.norm(rho0-spatial.regional_density(d, work, 'scalar_site0')))
    graph = spatial.interaction_graph(d)
    disconnected = spatial.interaction_graph(d, coupling=0.)
    checks = {
        'same_terminal_history_in_gns_and_regional_maps': marginal_error < 1e-12,
        'every_region_retains_total_history_weight': all(abs(row['state_trace']-1) < 1e-9 for row in regions),
        'specified_interactions_are_present': graph['edges'] == [
            ['scalar_site0', 'scalar_site1'], ['scalar_site0', 'gauge_matter_link']],
        'removing_charge_coupling_removes_its_edge': disconnected['edges'] == [
            ['scalar_site0', 'scalar_site1']],
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'dimension': d, 'access_factor_order': list(spatial.FACTOR_ORDER),
            'regions': regions, 'interaction_graph': graph,
            'zero_coupling_control_graph': disconnected,
            'same_state_marginal_error': marginal_error,
            'checks': checks, 'all_pass': all(checks.values())}


def cutoff_map_receipt():
    v = spatial.occupation_embedding()
    coarse, fine = ideal_hamiltonian(2), ideal_hamiltonian(4)
    compression_error = float(np.linalg.norm(v.conj().T@fine@v-coarse))
    blank = np.eye(len(coarse))[:, 0]
    leakage = float(np.linalg.norm((fine@v-v@coarse)@blank))
    if compression_error > 1e-12 or abs(leakage-float(F(9, 64))) > 1e-12:
        raise AssertionError('occupation map no longer has the stated compression/leakage')
    return {'coarse_dimension': 2, 'fine_dimension': 4,
            'hamiltonian_compression_error': compression_error,
            'blank_dynamical_intertwining_defect': leakage,
            'exact_blank_defect': '9/64',
            'is_exact_generator_natural_refinement': False}


def run(steps=512):
    histories = {d: experiment.history(d, steps) for d in (2, 4)}
    observers = [operational.observation_receipt(history) for history in histories.values()]
    regions = [regional_receipt(history) for history in histories.values()]
    clock = modular_clock.receipt(histories[2])
    repaired_histories = [history_repair.audit(d, steps) for d in (2, 4)]
    geometry = slow_geometry.report()
    return {
        'schema': 'oph.common_source.observer_bridge.v1',
        'source': 'reexecuted complete common-source history with native current 0 and feedback enabled',
        'observers': observers, 'regions': regions, 'modular_clock': clock,
        'schedule_independent_history_repair': repaired_histories,
        'native_slow_geometry': geometry,
        'occupation_map_control': cutoff_map_receipt(),
        'interacting_occupation_refinement': refinement.report(),
        'all_pass': all(row['all_pass'] for row in observers+regions+[clock, geometry]+repaired_histories),
        'scope': {
            'quantum_algebra_and_access_supplied': True,
            'regional_localization_supplied': True,
            'modular_clock_probe_and_schedule_supplied': True,
            'new_modular_clock_has_occupation_limit_certificate': False,
            'physical_clock_calibration_established': False,
            'spatial_continuum_established': False,
            'gravity_or_empirical_agreement_established': False,
            'full_lean_physical_tower_instantiated': False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, help='create a receipt, refusing overwrite')
    group.add_argument('--check', type=Path, help='reexecute and compare a retained receipt')
    args = parser.parse_args()
    result = run()
    if args.write:
        serialized = canonical(result)
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(serialized)
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'all_pass': result['all_pass'],
                      'observer_gns_dimensions': [row['gns_hilbert_dimension'] for row in result['observers']],
                      'modular_clock_gap': result['modular_clock']['clock_parameters']['modular_frequency'],
                      'occupation_refinement': result['interacting_occupation_refinement']['status'],
                      'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()
