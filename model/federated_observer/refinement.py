"""Exact observer restrictions and their measured dynamical refinement defect.

Identical twelve-port carriers remain the elementary objects.  ``parent``
retains the local port label and forgets the four-child cell address.  This
is a declared readout map, tested against the actual edge-dual source graph.
No physical clock or continuum identification is assigned to it.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from .source import federation, binary_step, cell_parents
from .source import SOURCE_PATHS as SOURCE_DEPENDENCIES


SOURCE_PATHS = (
    *SOURCE_DEPENDENCIES,
    'model/federated_observer/refinement.py',
    'model/federated_observer/source.py',
    'Lean/Tower/ConsensusTower.lean',
    'Lean/QFT/CPRestrictionNet.lean',
    'Lean/QFT/SourceRegionalNet.lean',
    'Lean/QFT/LocallyCovariantLimit.lean',
    'Lean/QFT/LimitRestriction.lean',
    'Lean/ObservableNormalForms/ObservableNormalForms/Refinement.lean',
    'Lean/Screen/DiscreteGauss.lean',
    'Lean/Screen/RegionalContinuity.lean',
)


def _levels(coarse, fine):
    if (type(coarse) is not int or type(fine) is not int
            or not 0 <= coarse <= fine <= 2):
        raise ValueError('ordered source levels 0 through 2 required')


def parent(vertex, coarse, fine):
    """Source face ancestry with the local twelve-port address unchanged."""
    _levels(coarse, fine)
    if type(vertex) is not int or not 0 <= vertex < federation(fine).ports:
        raise ValueError('a fine global port is required')
    return 12 * cell_parents(coarse, fine)[vertex // 12] + vertex % 12


@lru_cache(maxsize=6, typed=True)
def parent_map(coarse, fine):
    _levels(coarse, fine)
    cells = cell_parents(coarse, fine)
    assert cells == tuple(c // (4 ** (fine - coarse)) for c in range(len(cells)))
    return tuple(12 * cells[v // 12] + v % 12 for v in range(federation(fine).ports))


def lift(values, coarse, fine):
    """Unital *-homomorphism J on the commutative coordinate-site algebra."""
    _levels(coarse, fine)
    if len(values) != federation(coarse).ports:
        raise ValueError('one observable value per coarse port required')
    return tuple(values[p] for p in parent_map(coarse, fine))


def restrict(values, coarse, fine):
    """Uniform conditional mean R over descendants of each port address."""
    _levels(coarse, fine)
    if len(values) != federation(fine).ports:
        raise ValueError('one observable value per fine port required')
    totals = [0] * federation(coarse).ports
    for p, value in zip(parent_map(coarse, fine), values):
        totals[p] += value
    size = 4 ** (fine - coarse)
    return tuple(F(value, size) if isinstance(value, (int, F)) else value / size
                 for value in totals)


def static_certificate(coarse, fine):
    """Exact sparse receipts for J, R, E=JR and the uniform site state.

This is an algebra of functions on port labels, not the configuration
algebra of the conserved-particle experiment.  Uniform site measures do
not assert compatibility of its complete fixed-k occupation law.
"""
    pmap = parent_map(coarse, fine)
    nc, nf = federation(coarse).ports, federation(fine).ports
    size = 4 ** (fine - coarse)
    counts = Counter(pmap)
    assert tuple(sorted(counts)) == tuple(range(nc))
    assert set(counts.values()) == {size}
    assert nc * size == nf
    # Each column support of J is one disjoint fibre.  Its R sum is one;
    # multiplying columns gives the same indicator or zero.  These sparse
    # identities certify all observables, not a sampled vector test.
    rj_diagonal = tuple(F(counts[p], size) for p in range(nc))
    uniform_pairing = tuple(F(counts[p], nf) for p in range(nc))
    assert set(rj_diagonal) == {F(1)}
    assert set(uniform_pairing) == {F(1, nc)}
    return {
        'coarse_ports': nc, 'fine_ports': nf, 'fibre_size': size,
        'all_fibres_checked': len(counts),
        'restriction_times_lift_identity': True,
        'lift_preserves_products_and_adjoint': True,
        'conditional_expectation_idempotent_unital_positive': True,
        'uniform_site_state_preserved': True,
        'full_fixed_particle_configuration_state_claimed': False,
        'fixed_two_static_projection_control': {
            'fine_probability_two_particles_share_one_parent_port':
                _fraction(F(size - 1, nf - 1)),
            'coarse_binary_probability_two_particles_share_one_port': _fraction(F(0)),
            'block_count_alphabet': list(range(size + 1)),
        },
    }


def graph_parent_certificate(coarse, fine):
    """Check both cell ancestry and the stronger port-edge graph morphism."""
    pmap = parent_map(coarse, fine)
    cg, fg = federation(coarse), federation(fine)
    coarse_edges = set(cg.seams)
    coarse_cell_edges = {tuple(sorted((a // 12, b // 12)))
                         for a, b in cg.inter_seams}
    categories = Counter()
    witnessed = set()
    cell_bad, port_bad = [], []
    for a, b in fg.seams:
        pa, pb = pmap[a], pmap[b]
        pair = tuple(sorted((pa, pb)))
        if pa == pb:
            categories['collapsed'] += 1
        elif pair in coarse_edges:
            categories['coarse_edge'] += 1
            witnessed.add(pair)
        else:
            categories['not_coarse_edge'] += 1
            if not port_bad:
                port_bad.append({'fine_seam': [a, b], 'mapped_ports': [pa, pb]})
        ca, cb = pa // 12, pb // 12
        if ca != cb and tuple(sorted((ca, cb))) not in coarse_cell_edges:
            cell_bad.append((a, b))
    return {
        'coarse_level': coarse, 'fine_level': fine,
        'fine_edges_checked': len(fg.seams),
        'port_edge_categories': dict(sorted(categories.items())),
        'coarse_edges_with_direct_lift': len(witnessed),
        'coarse_edges_without_direct_lift': len(coarse_edges - witnessed),
        'cell_parent_graph_morphism': not cell_bad,
        'cell_parent_counterexample_count': len(cell_bad),
        'port_parent_graph_morphism': not port_bad,
        'first_port_edge_failure': port_bad[0] if port_bad else None,
    }


def aggregate_laplacian_rows(coarse, fine):
    """Sparse integer rows of L_f J, without constructing dense matrices."""
    fg = federation(fine)
    pmap = parent_map(coarse, fine)
    rows = []
    for v, neighbours in enumerate(fg.neighbours):
        row = Counter({pmap[v]: len(neighbours)})
        for w in neighbours:
            row[pmap[w]] -= 1
        rows.append({p: value for p, value in row.items() if value})
    return tuple(rows)


def defect_rows(coarse, fine, rate=F(1)):
    """D=(L_f J-c J L_c)/2 for per-edge attempt intensity, exactly."""
    if not isinstance(rate, (int, F)) or isinstance(rate, bool) or rate < 0:
        raise ValueError('a nonnegative rational coarse intensity rate required')
    cg = federation(coarse)
    rows = []
    for p, fine_row in zip(parent_map(coarse, fine),
                           aggregate_laplacian_rows(coarse, fine)):
        row = {q: F(value, 2) for q, value in fine_row.items()}
        row[p] = row.get(p, F(0)) - rate * F(len(cg.neighbours[p]), 2)
        for q in cg.neighbours[p]:
            row[q] = row.get(q, F(0)) + rate * F(1, 2)
        row = {q: value for q, value in row.items() if value}
        assert sum(row.values(), F(0)) == 0
        rows.append(row)
    return tuple(rows)


def defect_norm(coarse, fine, rate=F(1)):
    return max((sum(map(abs, row.values()), F(0))
                for row in defect_rows(coarse, fine, rate)), default=F(0))


def nonlumpability_witness(coarse, fine):
    """Two fixed-two-particle states share every block count and differ later.

The exact first-moment closure of binary exchange makes a common spectator
cancel from the difference, including edges between occupied sites.
"""
    pmap = parent_map(coarse, fine)
    fg = federation(fine)
    rows = aggregate_laplacian_rows(coarse, fine)
    first = {}
    for v, p in enumerate(pmap):
        if p not in first:
            first[p] = v
            continue
        u = first[p]
        for target in sorted(set(rows[u]) | set(rows[v])):
            if target != p and rows[u].get(target, 0) != rows[v].get(target, 0):
                spectator = next(w for w in range(fg.ports) if w not in (u, v))
                cu = -rows[u].get(target, 0)
                cv = -rows[v].get(target, 0)
                before = tuple(sorted((u, spectator)))
                after = tuple(sorted((v, spectator)))
                assert block_counts(before, coarse, fine) == block_counts(after, coarse, fine)
                return {
                    'fine_ports': [u, v], 'common_parent_port': p,
                    'target_parent_port': target,
                    'outgoing_target_edge_counts': [cu, cv],
                    'initial_fixed_two_states': [list(before), list(after)],
                    'one_source_attempt_block_mean_difference':
                        _fraction(F(cu - cv, 2 * len(fg.seams))),
                    'common_block_counts': True,
                    'no_scalar_time_rescaling_can_repair_exact_lumpability': True,
                }
    return None


def _multiply_lazy_numerator(values, graph, denominator):
    """(denominator I-L) times an integer-valued column observable."""
    return tuple((denominator - len(row)) * values[v]
                 + sum(values[w] for w in row)
                 for v, row in enumerate(graph.neighbours))


def finite_horizon_certificate(coarse, fine, steps=(1, 2, 4, 8)):
    """Exact indicator errors for a common uniformized Markov step.

    Q_l=I-L_l/(2*dmax) has common auxiliary intensity h=1/dmax.
    It is a uniformization kernel for the per-edge Poisson generator,
    rather than one source attempt.  The common comparison is explicit.
"""
    if (not steps or any(type(n) is not int or n < 1 for n in steps)
            or tuple(sorted(set(steps))) != tuple(steps)):
        raise ValueError('increasing positive integer comparison horizons required')
    cg, fg = federation(coarse), federation(fine)
    witness = nonlumpability_witness(coarse, fine)
    target = witness['target_parent_port'] if witness else 0
    dmax = max(max(map(len, cg.neighbours)), max(map(len, fg.neighbours)))
    denominator = 2 * dmax
    coarse_vector = tuple(int(v == target) for v in range(cg.ports))
    fine_vector = lift(coarse_vector, coarse, fine)
    norm = defect_norm(coarse, fine)
    results = []
    for n in range(1, max(steps) + 1):
        coarse_vector = _multiply_lazy_numerator(coarse_vector, cg, denominator)
        fine_vector = _multiply_lazy_numerator(fine_vector, fg, denominator)
        if n in steps:
            expected = lift(coarse_vector, coarse, fine)
            error = F(max(abs(a - b) for a, b in zip(fine_vector, expected)),
                      denominator ** n)
            # Signed residual rows sum to zero; a [0,1] observable costs
            # half the row L1 norm. Both kernels preserve that interval.
            bound = min(F(1), F(n, dmax) * norm / 2)
            assert error <= bound
            results.append({'steps': n, 'auxiliary_intensity': _fraction(F(n, dmax)),
                            'max_indicator_error': _fraction(error),
                            'telescope_indicator_bound': _fraction(bound)})
    return {'common_step_intensity': _fraction(F(1, dmax)),
            'common_kernel_denominator': denominator, 'target_block': target,
            'finite_horizons': results,
            'continuous_semigroup': 'exp(-t L/2)',
            'continuous_tv_bound': 'min(1, t*norm_inf(D)/2)',
            'parameter': 'auxiliary per-edge attempt intensity; no physical duration'}


def block_counts(occupied, coarse, fine):
    """Fine occupation totals visible to the declared parent-port blocks."""
    pmap = parent_map(coarse, fine)
    if (type(occupied) is not tuple or tuple(sorted(set(occupied))) != occupied
            or any(type(v) is not int or not 0 <= v < len(pmap) for v in occupied)):
        raise ValueError('canonical distinct occupied fine ports required')
    counts = [0] * federation(coarse).ports
    for v in occupied:
        counts[pmap[v]] += 1
    return tuple(counts)


def continue_counts(initial, records, coarse, fine, after_event=0):
    """Read only initial block counts and authentic complete signed records.

All incident crossing records are required; no fine configuration or
endpoint state is read. Internal currents cancel. Event IDs give order.
"""
    pmap = parent_map(coarse, fine)
    size = 4 ** (fine - coarse)
    if (type(initial) is not tuple or len(initial) != federation(coarse).ports
            or any(type(x) is not int or not 0 <= x <= size for x in initial)):
        raise ValueError('valid initial block counts required')
    if type(after_event) is not int or after_event < 0:
        raise ValueError('a nonnegative committed checkpoint event required')
    known_edges = set(federation(fine).seams)
    counts = list(initial)
    previous = after_event
    for record in records:
        if type(record) is not tuple or len(record) != 5:
            raise ValueError('record=(event,a,b,coin,current) required')
        event, a, b, coin, current = record
        if (any(type(x) is not int for x in record) or event <= previous
                or (a, b) not in known_edges
                or (coin, current) not in ((0, 0), (0, 1), (1, -1), (1, 0))):
            raise ValueError('ordered authentic-source-format record required')
        previous = event
        counts[pmap[a]] -= current
        counts[pmap[b]] += current
        if any(x < 0 or x > size for x in counts):
            raise ValueError('record contradicts the known block capacity')
    return tuple(counts)


def current_certificate(coarse, fine):
    """Every seam/endpoint/coin case, then continuing fixed-k histories."""
    graph = federation(fine)
    pmap = parent_map(coarse, fine)
    checked = 0
    # The native operation is identity off its two endpoints. Checking all
    # eight local cases is therefore an exhaustive symbolic spectator proof.
    for seam in graph.seams:
        a, b = seam
        for x, y, coin in product((0, 1), repeat=3):
            initial = tuple(v for v, bit in zip(seam, (x, y)) if bit)
            final, current = binary_step(initial, seam, coin)
            delta = Counter(pmap[v] for v in final)
            delta.subtract(pmap[v] for v in initial)
            expected = Counter()
            expected[pmap[a]] -= current
            expected[pmap[b]] += current
            assert all(delta[q] == expected[q] for q in set(delta) | set(expected))
            checked += 1
    crossing = next(edge for edge in graph.inter_seams if pmap[edge[0]] != pmap[edge[1]])
    internal = next(edge for edge in graph.seams if edge != crossing)
    spectator = next(v for v in range(graph.ports) if v not in crossing)
    initial = tuple(sorted((crossing[0], spectator)))
    tape = (crossing, crossing, internal, crossing, internal, crossing)
    initial_counts = block_counts(initial, coarse, fine)
    histories, checkpoints, nonzero = 0, 0, 0
    for coins in product((0, 1), repeat=len(tape)):
        state = initial
        records = []
        for event, (seam, coin) in enumerate(zip(tape, coins), 1):
            state, current = binary_step(state, seam, coin)
            nonzero += int(current != 0)
            records.append((event, *seam, coin, current))
            counts = continue_counts(initial_counts, records, coarse, fine)
            assert counts == block_counts(state, coarse, fine)
            if event == 3:
                checkpoint = counts
        resumed = continue_counts(checkpoint, records[3:], coarse, fine, after_event=3)
        assert resumed == counts
        checkpoints += 1
        histories += 1
    # Omitting an actual crossing current destroys the otherwise exact law.
    moved, current = binary_step(initial, crossing, 0)
    assert current == 1
    assert continue_counts(initial_counts, (), coarse, fine) != block_counts(moved, coarse, fine)
    return {
        'native_local_cases_checked': checked,
        'complete_coin_histories': histories,
        'record_events_per_history': len(tape),
        'probability_of_each_complete_coin_word': _fraction(F(1, 2 ** len(tape))),
        'zero_current_events_retained': histories * len(tape) - nonzero,
        'nonzero_current_events': nonzero,
        'prefix_block_readout_equalities': histories * len(tape),
        'checkpoint_continuations_equal': checkpoints,
        'initial_fixed_two_occupation': list(initial),
        'seam_tape': [list(edge) for edge in tape],
        'omitted_boundary_record_control_detects_mismatch': True,
        'access': 'known block preparation plus all authentic signed crossing currents',
        'markov_closure_from_boundary_records_claimed': False,
    }


def _fraction(value):
    value = F(value)
    return {'numerator': value.numerator, 'denominator': value.denominator}


def run():
    """JSON-safe exact refinement subsection for the combined packet."""
    result = []
    for coarse, fine in ((0, 1), (1, 2), (0, 2)):
        witness = nonlumpability_witness(coarse, fine)
        assert witness is not None
        result.append({
            'coarse_level': coarse, 'fine_level': fine,
            'static_site_algebra': static_certificate(coarse, fine),
            'source_graph_parent': graph_parent_certificate(coarse, fine),
            'nonlumpability': witness,
            'generator_defect_norm_inf': _fraction(defect_norm(coarse, fine)),
            'coarse_fourfold_rate_defect_norm_inf': _fraction(defect_norm(coarse, fine, F(4))),
            'finite_horizon': finite_horizon_certificate(coarse, fine),
            'retained_boundary_currents': current_certificate(coarse, fine),
        })
    p02 = parent_map(0, 2)
    p01, p12 = parent_map(0, 1), parent_map(1, 2)
    assert p02 == tuple(p01[p] for p in p12)
    return {'schema': 'oph.federated_observer.refinement.v1',
            'comparison_parameter': 'auxiliary per-edge attempt intensity',
            'elementary_ports_per_carrier': 12,
            'ancestry_composition_checked_ports': len(p02),
            'state_scope': 'uniform state on the commutative coordinate-site algebra',
            'exact_dynamical_tower_claimed': False,
            'continuum_limit_claimed': False,
            'refinements': result}


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
