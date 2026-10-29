"""Actual source joins and countercontrols across the four provisional lanes."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from model.local_repair_family import local_bias, modular_clock
from model.moving_defect import transport
from model.particle_obstruction import classification
from model.record_action import action
from model.observer.source_discrimination import circuit

ROOT = Path(__file__).resolve().parents[2]


def test_defect_parent_and_family_use_same_fixed_endpoint_source():
    p = modular_clock.parameters()
    assert transport.parameters() == p
    family = classification.family_certificate()
    assert family['clock_parameters'] == p.description()
    h = action.kinetic_parent(tuple(map(F, p.port_weights)))
    # Independently lift the actual residual port permutation to all66 pairs.
    g = next(g for g, _ in classification.fixed_reference_stabilizer(p)
             if g != tuple(range(12)))
    states = local_bias.basis(); index = {x: i for i, x in enumerate(states)}
    permutation = [index[tuple(x[g.index(j)] for j in range(12))] for x in states]
    assert np.max(abs(h[np.ix_(permutation, permutation)]-h)) < 1e-14
    pi = np.array(list(map(float, local_bias.reference(p))))
    assert np.array_equal(pi[permutation], pi)
    # Shared residual symmetry still does not identify kinetic/modular flows.
    assert np.linalg.norm(h*(pi[None, :]-pi[:, None])) > 1e-3


def test_finite_record_point_estimate_does_not_certify_exact_symmetry():
    receipt = json.loads((ROOT/'model/local_repair_family/record_clock_results.json').read_text())
    counts = {tuple(b['seam']): (b['trials'], b['coin1_count']) for b in receipt['acquisition']['blocks']}
    estimate = action.infer_parent(counts)
    assert estimate['status'] == 'accept'
    weights = estimate['relative_weights']
    preserving = [g for g, _ in classification.rotation_representation()
                  if all(weights[g[i]] == weights[i] for i in range(12))]
    assert preserving == [tuple(range(12))]
    # Exact truth has C2 while this accepted point estimate has only identity.
    # Thus estimates must not be used as a new exact symmetry premise.
    assert len(classification.fixed_reference_stabilizer(modular_clock.parameters())) == 2


def test_scalar_clock_compression_does_not_preserve_the_moving_code():
    core = (0, 2)
    v = transport.ideal_code(core)
    weights = modular_clock.parameters().port_weights
    diagonal = np.ones(512, complex)
    parameter = .7
    for state in local_bias.basis():
        product = np.prod([weights[i] for i, n in enumerate(state) if n])
        diagonal[circuit.encode(state)] = np.exp(1j*parameter*np.log(product))
    evolved = diagonal[:, None]*v
    compressed = v.conj().T@evolved
    assert np.max(abs(compressed-compressed[0, 0]*np.eye(2))) < 1e-14
    leakage = evolved-v@compressed
    row = next(r for r in transport.clock_compatibility()['rows'] if r['core'] == list(core))
    predicted = float(F(row['modular_code_leakage_sine_squared_coefficient']))*np.sin(parameter*np.log(25/8)/2)**2
    assert predicted > 0
    assert np.max(abs(leakage.conj().T@leakage-predicted*np.eye(2))) < 1e-14
