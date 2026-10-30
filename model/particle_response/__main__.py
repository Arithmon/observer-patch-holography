"""Execute the current complete-reference particle and probe checks."""
import json

from .dimer import certificate, setup
from .source import probe_readout
from .spectrum import FieldSpectrum
from ..full_fock_field.field import FieldGeometry


def report():
    thermal = FieldSpectrum()
    ground = setup()
    return dict(
        scope='FINITE_COMPLETE_FIELD_PARTICLE_AND_HOLE_RESPONSE',
        half_filled_particle_candidates=certificate(),
        faithful_beta_one_mu_zero_response=thermal.certify(),
        half_filled_pure_ground_response=ground.certify(),
        source_ground_single_probe=probe_readout(ground),
        source_thermal_single_probe=probe_readout(thermal),
        source_atomic_thermal_pair_probe=probe_readout(
            FieldSpectrum(FieldGeometry((1, 1, 1))), 'pair', steps=1),
        absolute_physical_mass_prediction=False,
        topological_defect_assignment=False,
        Lorentz_dispersion_identified=False,
        physical_particle_identification=False)


if __name__ == '__main__':
    print(json.dumps(report(), indent=2, allow_nan=False))
