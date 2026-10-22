"""Exact seam-symmetry selection and a temporal-scheduling countercontrol.

The native edge orbit fixes a symmetry-invariant *one-step* seam law. It does
not fix correlations between successive seam choices. All probabilities below
are exact and refer to the binary total-two sector and native attempt index.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

from ..common_source import slow_geometry
from ..observer import native


# Concrete source-port automorphisms. Each is checked against every seam
# before its orbit is used; no unverified coordinate transformation enters.
GENERATORS = (
    (7, 0, 3, 4, 9, 1, 2, 10, 11, 5, 6, 8),
    (0, 1, 2, 3, 6, 7, 4, 5, 9, 8, 11, 10),
)


def edge_orbit(seams, generators=GENERATORS):
    """Certify the orbit of the first declared seam under graph symmetries."""
    seams = tuple(tuple(edge) for edge in seams)
    edge_set = set(seams)
    if len(edge_set) != len(seams) or not seams or any(
            len(edge) != 2 or not 0 <= edge[0] < edge[1] < 12 for edge in seams):
        raise ValueError('distinct increasing native seams required')
    for permutation in generators:
        if tuple(sorted(permutation)) != tuple(range(12)):
            raise ValueError('a complete twelve-port permutation required')
        if {tuple(sorted((permutation[a], permutation[b]))) for a, b in seams} != edge_set:
            raise ValueError('a proposed generator does not preserve native seams')
    seen, frontier = {seams[0]}, [seams[0]]
    while frontier:
        a, b = frontier.pop()
        for permutation in generators:
            image = tuple(sorted((permutation[a], permutation[b])))
            if image not in seen:
                seen.add(image)
                frontier.append(image)
    return frozenset(seen)


@lru_cache(maxsize=1)
def _native_two_current_coefficients():
    """Directly sum two actual endpoint/coin attempts over all 66 states."""
    seams = slow_geometry.source().seams()
    states = native.states(binary=True)
    if len(seams) != 30 or len(states) != 66 or len(edge_orbit(seams)) != 30:
        raise ArithmeticError('the source carrier or its symmetry orbit changed')
    # Every fixed-seam kernel preserves the uniform configuration measure.
    # Hence the measure remains uniform even conditional on a seam history;
    # correlated choices of seams do not change the stationary preparation.
    index = {state: i for i, state in enumerate(states)}
    for seam in seams:
        incoming = [0]*len(states)
        for state in states:
            for coin in (0, 1):
                after, _ = native.native(state, seam, coin)
                incoming[index[after]] += 1
        if set(incoming) != {2}:
            raise ArithmeticError('a fixed native seam no longer preserves uniform preparation')
    first, iid_both, repeated_both = F(), F(), F()
    for state in states:
        for first_edge in seams:
            for first_coin in (0, 1):
                after, (_, first_current) = native.native(state, first_edge, first_coin)
                if not first_current:
                    continue
                first += F(1, 66*30*2)
                for second_edge in seams:
                    for second_coin in (0, 1):
                        _, (_, second_current) = native.native(after, second_edge, second_coin)
                        if second_current:
                            iid_both += F(1, 66*30*30*4)
                            if second_edge == first_edge:
                                repeated_both += F(1, 66*30*4)
    if (first, iid_both, repeated_both) != (F(5, 33), F(23, 990), F(5, 66)):
        raise ArithmeticError('the native current law changed')
    return first, iid_both, repeated_both


def persistent_seam_current_law(repeat_probability):
    """Stationary invariant seam Markov law: repeat, otherwise redraw uniformly.

    Seam at the first attempt is uniform. At later attempts the previous seam
    is repeated with probability q; with probability 1-q a fresh uniform seam
    is drawn. For 0 <= q < 1 this scheduler still exposes every seam and is
    irreducible and aperiodic on its edge state. The coin remains fair.
    """
    if (type(repeat_probability) is not F or
            not F() <= repeat_probability < F(1)):
        raise ValueError('a rational repeat probability in [0,1) required')
    rate, iid_both, repeated_both = _native_two_current_coefficients()
    return {'repeat_probability': str(repeat_probability),
            'one_step_seam_probability': '1/30',
            'stationary_moving_current_probability_each_attempt': str(rate),
            'two_successive_moving_currents_probability': str(
                (1-repeat_probability)*iid_both+repeat_probability*repeated_both),
            'all_one_step_seam_marginals_uniform': True,
            'seam_choices_independent': repeat_probability == 0}


def repeat_probability_from_current_pair(two_move_probability):
    """Identify q from exact retained bits, conditional on this scheduler family."""
    if (type(two_move_probability) is not F or
            not F(23, 990) <= two_move_probability < F(5, 66)):
        raise ValueError('an exact two-current probability in this scheduler family required')
    return (990*two_move_probability-23)/52


def run():
    seams = slow_geometry.source().seams()
    orbit = edge_orbit(seams)
    iid = persistent_seam_current_law(F())
    persistent = persistent_seam_current_law(F(1, 2))
    if (orbit != frozenset(seams)
            or iid['two_successive_moving_currents_probability'] != '23/990'
            or persistent['two_successive_moving_currents_probability'] != '49/990'
            or repeat_probability_from_current_pair(F(49, 990)) != F(1, 2)):
        raise ArithmeticError('native source symmetry or scheduling contrast failed')
    return {'edge_orbit_size': len(orbit),
            'native_edges_checked_under_both_generators': 30,
            'normalized_symmetry_invariant_one_step_seam_weight': '1/30',
            'iid_uniform': iid, 'uniform_marginal_persistent': persistent,
            'same_stationary_uniform_configuration_law': True,
            'same_one_attempt_native_current_law': True,
            'two_attempt_native_current_laws_equal': False,
            'symmetry_or_stationarity_selects_iid_schedule': False,
            'repeat_probability_identifiable_from_exact_two_current_law_within_declared_family': True,
            'finite_sample_uncertainty_or_scheduler_family_selected': False,
            'physical_time_or_particle_mass_identified': False}
