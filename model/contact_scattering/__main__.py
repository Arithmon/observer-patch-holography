"""Execute the current finite interaction and fixed-broadening controls."""
import json
from fractions import Fraction as F

from .collision import collision_certificate, execute
from .scattering import response, volume_enclosure


def main():
    rows=[]
    for side in (8,32,128):
        value=response([0,0,0],.125+.25j,side=side)
        bound=volume_enclosure(side,F(1,4))
        rows.append({'side':side,
                     'tau_real_imag':[value['tau'].real,value['tau'].imag],
                     'optical_residual':value['optical_residual'],
                     'free_green_error_bound_display':float(bound['free_green_error']),
                     'tau_error_bound_display':float(bound['tau_error'])})
    print(json.dumps({'scope':'SUPPLIED_WEYL_CONTACT_FIELD_INTERACTION',
                      'retarded_frequency_real_imag':[.125,.25],
                      'fixed_broadening_volume_rows':rows,
                      'source_interaction_echo':execute(2)['report'],
                      'zero_interaction_source_control':execute(2,interaction=0)['report'],
                      'rational_signal_certificate':collision_certificate(),
                      'physical_particle_mass_or_cross_section_identified':False},
                     indent=2,allow_nan=False))


if __name__=='__main__':
    main()
