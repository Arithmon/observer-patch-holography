"""Independent matrix controls for all existing CKM/PMNS readout entry points."""
import importlib
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for directory in (ROOT, ROOT/'particles/flavor', ROOT/'particles/neutrino'):
    sys.path.insert(0, str(directory))

READERS = [
    ('neutrino.build_pmns_from_shared_flavor_basis', '_standard_pmns_parameters'),
    ('neutrino.derive_neutrino_dimensionless_law_candidate_audit', '_pmns_parameters'),
    ('neutrino.derive_neutrino_two_parameter_exact_adapter', '_pmns_parameters'),
    ('neutrino.derive_neutrino_weighted_cycle_repair', '_pmns_parameters'),
    ('neutrino.derive_neutrino_physical_majorana_phase_theorem', '_standard_pmns_parameters'),
    ('neutrino.derive_neutrino_weighted_cycle_shared_basis_representation', '_standard_pmns_parameters'),
    ('flavor.derive_quark_d12_mass_branch_and_ckm_residual', '_standard_ckm_parameters'),
    ('flavor.sigma_ud_orbit_provider', '_standard_ckm_parameters'),
    ('flavor.enumerate_quark_local_basis_orbit_diagnostic', '_ckm_tuple'),
    ('flavor.derive_quark_transport_frame_diagnostic_orbit', '_ckm_tuple'),
]


@pytest.mark.parametrize('module,name', READERS)
def test_every_reader_rejects_a_nonunitary_matrix(module, name):
    reader = getattr(importlib.import_module('particles.'+module), name)
    with pytest.raises(ValueError):
        reader(np.full((3, 3), .2, dtype=complex))
