"""A resource certificate must not admit schedules rejected by execution."""
from fractions import Fraction

import pytest

from . import source
from .lattice import WeylLattice


@pytest.mark.parametrize('name,value', [
    ('steps', -1), ('steps', 0), ('steps', 1.5), ('steps', True),
    ('steps', Fraction(3, 2)), ('segments', -1), ('segments', 0),
    ('segments', 1.5), ('segments', True), ('duration', float('nan')),
    ('duration', float('inf')), ('duration', True),
])
def test_invalid_schedule_cannot_emit_a_resource_certificate(name, value):
    arguments = dict(duration=1., steps=1, segments=1)
    arguments[name] = value
    # Admission must fail before accessing a lattice or compiling a source word.
    with pytest.raises(ValueError):
        source.resources(object(), **arguments)


def test_valid_source_counts_are_integral_positive_and_add_under_concatenation():
    lattice = WeylLattice(1)
    forward = source.resources(lattice, .25, 2, segments=1)
    backward = source.resources(lattice, -.25, 2, segments=1)
    twice = source.resources(lattice, .25, 2, segments=2)
    for key in (*source.COUNT_KEYS, 'completed_hopping_blocks',
                'completed_contact_blocks', 'routing_FSWAP_calls'):
        assert type(forward[key]) is int and forward[key] > 0
        assert backward[key] == forward[key]
        assert twice[key] == 2*forward[key]
