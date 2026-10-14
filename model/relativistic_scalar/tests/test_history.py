from fractions import Fraction as F

import numpy as np
import pytest

from model.relativistic_scalar import bounds, exact, experiment
from model.relativistic_scalar.lattice import ScalarLattice, ScalarParameters


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('connected', [False, True])
def test_actual_source_complete_history_against_independent_sparse_evolution(sign, connected):
    lattice = ScalarLattice(experiment.LEVELS[0])
    actual = experiment.history(lattice, 4, sign, connected)
    reference = experiment.history(lattice, 4, sign, connected, source_execution=False)
    assert [row['bits'] for row in actual] == [[0, 0], [0, 1], [1, 0], [1, 1]]
    assert abs(sum(row['probability'] for row in actual)-1) < 1e-11
    assert experiment.tv(actual, reference) < 2e-8
    assert experiment.tv(actual, reference) < float(2*bounds.product_error(
        lattice.dimensions, experiment.INTERVAL, 4))+1e-11


def test_committed_first_record_changes_later_probabilities_without_reset():
    lattice = ScalarLattice(experiment.LEVELS[0])
    connected = experiment.history(lattice, 4)
    disconnected = experiment.history(lattice, 4, connected=False)
    assert experiment.tv(connected, disconnected) > 1e-3
    # The committed first-bit law is independent of later feedback.
    for first in (0, 1):
        left = sum(row['probability'] for row in connected if row['bits'][0] == first)
        right = sum(row['probability'] for row in disconnected if row['bits'][0] == first)
        assert abs(left-right) < 1e-11


def test_native_attachment_retains_every_coin_and_quantum_branch():
    histories = [experiment.history(ScalarLattice(experiment.LEVELS[0]), 4, sign)
                 for sign in (-1, 1)]
    rows = experiment.native_records(histories)
    assert len(rows) == 1056
    assert abs(sum(row['probability'] for row in rows)-1) < 1e-11
    assert len({(row['initial'], *row['records']) for row in rows}) == 1056
    detached = experiment.native_records(histories, connected=False)
    assert [row['probability'] for row in rows] == [row['probability'] for row in detached]
    assert any(row['final_port_zero'] != control['final_port_zero']
               for row, control in zip(rows, detached))


def test_quartic_pointer_witness_survives_source_product_error():
    result = experiment.witness(experiment.LEVELS[0], 8)
    assert result['strictly_positive_certified_margin']
    assert result['contrast'] > result['ideal_executed_contrast_lower_float'] > 0
    full = exact.interaction_witness_certificate(experiment.WITNESS_TIME)
    assert full['contrast_lower'] > 0


def test_particle_creation_and_local_tensor_cutoff_are_not_number_conservation():
    lattice = ScalarLattice(experiment.LEVELS[0])
    state = lattice.vacuum_state()
    advanced = lattice.evolve(state, float(F(1, 128)))
    occupation_four = lattice.basis_state((4,)+(0,)*7)
    assert abs(np.vdot(occupation_four, advanced)) > 1e-5
    free = ScalarLattice(experiment.LEVELS[0], ScalarParameters(coupling=0))
    # Its vacuum-to-four derivative vanishes; it can appear at higher order.
    assert abs(np.vdot(occupation_four, free.hamiltonian@state)) < 1e-14
    assert abs(np.vdot(occupation_four, lattice.hamiltonian@state)) > 1e-3


def test_resource_ledger_counts_actual_words_and_every_tensor_site():
    lattice = ScalarLattice(experiment.LEVELS[0])
    one = experiment.source_cost(lattice, 1)
    two = experiment.source_cost(lattice, 2)
    assert one['native_pulses'] > 0
    assert two['native_pulses'] > one['native_pulses']
    assert one['proper_code_service_events'] > 0
