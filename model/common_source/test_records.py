"""Operator-level checks of the checkpoint adapter and its creation boundary."""
from fractions import Fraction as Q
from itertools import product
import unittest

from model.observer import native
from model.common_source import records as r
from model.common_source import gauge_matter as g


def unit(size, row, column, value=r.ONE):
    return tuple(tuple(value if (i, j) == (row, column) else r.ZERO
                       for j in range(size)) for i in range(size))


def plus(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = g.matrix_add(result, matrix)
    return result


def trace(matrix):
    return (sum((matrix[i][i][0] for i in range(len(matrix))), Q()),
            sum((matrix[i][i][1] for i in range(len(matrix))), Q()))


def record_projector(record, count, dimension):
    return plus(*(unit(count*dimension, record*dimension+k,
                       record*dimension+k) for k in range(dimension)))


def branches():
    # Actual common-source exact unitaries, with all two-outcome QND branches.
    return tuple(tuple(g.matmul(g.occupation_projector(0, outcome), operation)
                       for outcome in (0, 1))
                 for operation in (g.ideal_electric(), g.ideal_echo()))


class CheckpointRecordsTests(unittest.TestCase):
    def test_product_access_is_exactly_the_first_mode_parity_commutant(self):
        result = r.product_access_certificate()
        self.assertEqual(result['rational_closure_lower_rank'], 128)
        self.assertEqual(result['exact_first_mode_parity_upper_dimension'], 128)
        self.assertEqual(result['disconnected_coupling_access_dimension'], 16)
        self.assertEqual(result['extra_first_mode_pointer_access_dimension'], 256)
        self.assertFalse(result['full_B_C16_access_from_declared_menu'])
        self.assertFalse(result['fixed_two_read_program_algebra_identified'])
        self.assertTrue(result['local_source_service_menu_supplied'])

    def test_product_access_uses_actual_source_terms_and_mode_order(self):
        import numpy as np
        from model.interacting_band_extension.construction import ProductBand, normalized
        from model.interacting_band_extension import algebra
        lattice, data = ProductBand(), r.product_access_generators()
        metric = np.asarray(data['bargmann_metric'], dtype=float)
        def orthonormal(matrix):
            return np.asarray(matrix, dtype=float)*np.sqrt(metric[:, None]/metric[None, :])
        for term, matrix in zip(lattice.terms, data['terms'], strict=True):
            np.testing.assert_allclose(orthonormal(matrix), lattice.embed(term).toarray(),
                                       atol=1e-14, rtol=0)
        coordinate = orthonormal(data['coordinate'])
        np.testing.assert_allclose(coordinate, lattice.apply(np.eye(16), (1,),
            normalized(algebra.polynomial(4, 1))), atol=1e-14, rtol=0)
        onsite = orthonormal(data['terms'][1])
        self.assertGreater(np.linalg.norm(coordinate@onsite-onsite@coordinate), 1)
        # The parity off-diagonal matrix unit is invisible to every allowed
        # algebra element, whereas the excluded first-mode pointer sees it.
        self.assertEqual(data['coordinate'][0][1], 0)
        self.assertNotEqual(data['excluded_first_coordinate'][0][1], 0)

    def test_actual_product_pointer_pulls_dynamic_effect_into_durable_record(self):
        import numpy as np
        from scipy.linalg import block_diag
        from model.interacting_band_extension.construction import ProductBand, normalized
        from model.interacting_band_extension import algebra
        from model.scalar_observations import source
        lattice = ProductBand()
        q = normalized(algebra.polynomial(4, 1))
        values, basis = np.linalg.eigh(q)
        self.assertLess(max(abs(values)), np.pi/2)
        pointer, _ = source.instrument(4, 1, np.pi/2, 2)
        full = tuple(lattice.apply(np.eye(16), (1,), matrix) for matrix in pointer)
        write = np.vstack(full)
        np.testing.assert_allclose(write.conj().T@write, np.eye(16), atol=2e-10, rtol=0)
        sine = lattice.apply(np.eye(16), (1,), (basis*np.sin(values))@basis.T)
        record = np.diag([0.]*16+[1.]*16)
        np.testing.assert_allclose(write.conj().T@record@write,
                                   (np.eye(16)+sine)/2, atol=2e-10, rtol=0)
        kick, _ = source.feedback(4, .25, 2)
        control = block_diag(np.eye(16), lattice.apply(np.eye(16), (1,), kick))
        np.testing.assert_allclose(control.conj().T@record@control, record,
                                   atol=2e-10, rtol=0)
        # The new output event is not a durable input observable: its input
        # effect is nonprojective and changes under the actual onsite service.
        effect = (np.eye(16)+sine)/2
        self.assertGreater(np.linalg.norm(effect@effect-effect), .1)
        onsite = lattice.embed(lattice.terms[1]).toarray()
        self.assertGreater(np.linalg.norm(onsite@effect-effect@onsite), .1)

    def test_native_all_input_gram_and_intertwiner_certificate(self):
        receipt = r.native_checkpoint_certificate()
        self.assertEqual(receipt["all_phase_column_counts"],
                         {"0": 78, "1": 468, "2": 2808})
        self.assertEqual(receipt["source_tape_dimension"], 3354)
        self.assertEqual(receipt["postcommit_dimension"], 3276)
        self.assertTrue(receipt["postcommit_dual_projectors_fixed"])
        self.assertFalse(receipt["full_operational_observer_inhabitant_claimed"])
        observer = receipt['postcommit_operational_observer']
        self.assertTrue(observer['all_OperationalObserver_clauses_realized_on_declared_postcommit_tower'])
        self.assertFalse(observer['Lean_inhabitant_compiled'])

    def test_actual_checkpoint_state_public_readback_and_native_generator(self):
        state = r.checkpoint_preparation()
        self.assertEqual(len(state), 132)
        self.assertTrue(all(left == right and value == (Q(1, 132), Q(0))
                            for (left, right), value in state.items()))
        identity = {}
        for record in native.LETTERS:
            projector = r.checkpoint_record_element(record)
            self.assertEqual(len(projector), 546)
            identity.update(projector)
            self.assertEqual(r.checkpoint_public_readback(projector), projector)
            self.assertEqual(r.checkpoint_generator(projector), {})
        self.assertEqual(r.checkpoint_dual(identity), identity)
        basis = r.checkpoint_basis()
        a, b = basis[0], basis[-1]
        operator = {(a, a): (Q(2), Q(1)), (a, b): (Q(3), Q(-1))}
        public = r.checkpoint_public_readback(operator)
        self.assertEqual(len(public), 546)
        self.assertEqual(set(public.values()), {(Q(1, 273), Q(1, 546))})
        self.assertEqual(r.checkpoint_public_readback(public), public)
        self.assertEqual(r._sum(public.values()), (Q(2), Q(1)))
        self.assertEqual(r.checkpoint_public_readback({(a, b): r.ONE}), {})

    def test_checkpoint_adjoint_matches_original_coherent_native_map(self):
        # Independent forward/dual pairing on matrix units includes coherent
        # coin cross terms and the separate phase Kraus branches.
        def pairing(effect, state):
            return r._sum(r._multiply(value, state.get((right, left), r.ZERO))
                          for (left, right), value in effect.items())
        states = [native.State(native.states()[j], (record,))
                  for j in (0, 17) for record in native.LETTERS]
        for a, b in product(states, repeat=2):
            operator = {(a, b): (Q(2, 7), Q(-3, 11))}
            forward = r.native_tick_operator(operator)
            c, d = next(iter(forward))
            effect = {(d, c): (Q(1, 3), Q(2, 5))}
            self.assertEqual(pairing(effect, forward),
                             pairing(r.checkpoint_dual(effect), operator))
            phase_cross = {(a, c): r.ONE}
            self.assertEqual(r.checkpoint_dual(phase_cross), {})
        a = states[0]
        self.assertEqual(r.checkpoint_generator({(a, a): r.ONE}),
                         {(a, a): (Q(-1), Q(0))})

    def test_evidence_prediction_clauses_do_not_identify_forecast_with_event_mass(self):
        labels = tuple(product((0, 1), native.LETTERS))
        state = r.checkpoint_preparation()
        read = lambda label: r.native_record_weight(state, label[1])
        for x, y in product(labels, repeat=2):
            self.assertFalse(r.checkpoint_precedes(x, x))
            self.assertFalse(r.checkpoint_precedes(r.checkpoint_predict(x), x))
            if read(x) == read(y):
                self.assertEqual(read(r.checkpoint_predict(x)), read(r.checkpoint_predict(y)))
            for z in labels:
                if r.checkpoint_precedes(x, y) and r.checkpoint_precedes(y, z):
                    self.assertTrue(r.checkpoint_precedes(x, z))
        self.assertTrue(r.checkpoint_precedes((0, (0, 0)), (1, (0, 0))))
        self.assertEqual(read((0, (0, 0))), read((0, (1, 0))))
        self.assertNotEqual(native.FORECAST[(0, 0)], native.FORECAST[(1, 0)])
        certificate = r.checkpoint_observer_certificate()
        self.assertTrue(certificate['record_conditioned_controls_are_nonconstant'])
        self.assertFalse(certificate['record_creation_in_persistence_domain'])

    def test_exact_buffer_refinement_state_pairing_composition_and_generator(self):
        state = r.checkpoint_preparation()
        a = next(iter(state))[0]
        effect = {(a, a): (Q(2), Q(3))}
        def pair(left, right):
            return r._sum(r._multiply(value, right.get((b, a), r.ZERO))
                          for (a, b), value in left.items())
        for bits in (0, 1, 2):
            lifted = r.checkpoint_buffer_lift(effect, bits)
            prepared = r.checkpoint_buffer_lift(state, bits, density=True)
            self.assertEqual(pair(prepared, lifted), pair(state, effect))
            # Applying native continuation independently on each unchanged
            # buffer block intertwines the generator, not just expectations.
            actual = {}
            for address in range(2**bits):
                block = {(x, y): value for ((x, u), (y, v)), value in lifted.items()
                         if u == v == address}
                actual.update({((x, address), (y, address)): value
                               for (x, y), value in r.checkpoint_generator(block).items()})
            self.assertEqual(actual, r.checkpoint_buffer_lift(r.checkpoint_generator(effect), bits))
        once = r.checkpoint_buffer_lift(effect, 1)
        twice = {((a, 2*u+w), (b, 2*v+w)): value
                 for ((a, u), (b, v)), value in once.items() for w in (0, 1)}
        self.assertEqual(twice, r.checkpoint_buffer_lift(effect, 2))
        # Normalizing the observable embedding would erase its state pairing.
        wrong = {pair: tuple(part/4 for part in value)
                 for pair, value in r.checkpoint_buffer_lift(effect, 2).items()}
        self.assertNotEqual(pair(r.checkpoint_buffer_lift(state, 2, density=True), wrong),
                            pair(state, effect))

    def test_persistence_observer_rejects_blank_and_invalid_interfaces(self):
        blank = native.State(native.states()[0])
        for function in (r.checkpoint_dual, r.checkpoint_public_readback, r.checkpoint_generator):
            with self.assertRaisesRegex(ValueError, 'blank sector'):
                function({(blank, blank): r.ONE})
        for label in ((False, (0, 0)), (2, (0, 0)), (0, (True, 0)), (0, (1, 2))):
            with self.assertRaises(ValueError):
                r.checkpoint_predict(label)
        for bits in (True, -1, 4):
            with self.assertRaises(ValueError):
                r.checkpoint_buffer_lift({}, bits)

    def test_native_first_write_is_not_a_global_fixed_record(self):
        state = native.State((1, 1) + (0,)*10)
        before = {(state, state): r.ONE}
        after = r.native_tick_operator(before)
        self.assertEqual(r.native_record_weight(before, (0, 0)), r.ZERO)
        self.assertEqual(r.native_record_weight(after, (0, 0)), (Q(1, 2), Q(0)))
        # A coherent coin write retains cross-coin matrix entries.
        self.assertEqual(len(after), 4)
        self.assertEqual(len(r.native_block_readback(after)), 2)

    def test_committed_projector_expectations_on_off_diagonal_matrix_units(self):
        # Include unreachable first-letter combinations and imaginary operators.
        source = native.states()[0]
        other = native.states()[17]
        for a in native.LETTERS:
            for b in native.LETTERS:
                for right_loads in (source, other):
                    left = native.State(source, (a,))
                    right = native.State(right_loads, (b,))
                    operator = {(left, right): (Q(0), Q(3, 7))}
                    after = r.native_tick_operator(operator)
                    for record in native.LETTERS:
                        self.assertEqual(r.native_record_weight(after, record),
                                         r.native_record_weight(operator, record))
                    self.assertEqual(r.native_block_readback(after),
                                     r.native_tick_operator(r.native_block_readback(operator)))

    def test_full_sector_halts_even_on_coherences_and_phase_cross_terms_dephase(self):
        a = native.State(native.states()[2], ((0, 0), (1, -1)))
        b = native.State(native.states()[3], ((1, 1), (0, -1)))
        operator = {(a, b): (Q(-2, 3), Q(5, 7))}
        self.assertEqual(r.native_tick_operator(operator), operator)
        blank = native.State(a.loads)
        self.assertEqual(r.native_tick_operator({(blank, a): r.ONE}), {})

    def test_reachable_quantum_record_coherence_changes_but_classical_values_survive(self):
        def loaded(*ports):
            loads = [0]*12
            for port in ports:
                loads[port] += 1
            return native.State(tuple(loads))
        a, b = loaded(0, 2), loaded(1, 2)
        initial = {(x, y): (Q(1, 2), Q(0)) for x in (a, b) for y in (a, b)}
        first = r.native_tick_operator(initial)
        second = r.native_tick_operator(first)
        first_record = r.native_first_record_matrix(first)
        second_record = r.native_first_record_matrix(second)
        i = native.LETTERS.index((0, 1))
        j = native.LETTERS.index((0, 0))
        self.assertEqual(first_record[i][j], (Q(1, 4), Q(0)))
        self.assertEqual(second_record[i][j], (Q(1, 8), Q(0)))
        self.assertEqual([first_record[k][k] for k in range(6)],
                         [second_record[k][k] for k in range(6)])

    def test_block_and_public_readbacks_have_different_fixed_algebras(self):
        work_coherence = unit(8, 1, 2, (Q(2, 3), Q(-3, 5)))
        record_coherence = unit(8, 1, 6, (Q(1, 3), Q(2, 5)))
        projector = record_projector(1, 2, 4)
        operator = plus(work_coherence, record_coherence, projector)
        self.assertEqual(r.block_readback(work_coherence, 2), work_coherence)
        self.assertEqual(r.public_record_readback(work_coherence, 2), unit(8, 0, 0, r.ZERO))
        self.assertEqual(r.block_readback(record_coherence, 2), unit(8, 0, 0, r.ZERO))
        self.assertEqual(r.public_record_readback(projector, 2), projector)
        block = r.block_readback(operator, 2)
        public = r.public_record_readback(operator, 2)
        self.assertEqual(r.block_readback(block, 2), block)
        self.assertEqual(r.public_record_readback(public, 2), public)
        self.assertEqual(r.public_record_readback(block, 2), public)
        self.assertEqual(trace(operator), trace(block))
        self.assertEqual(trace(operator), trace(public))

    def test_pinching_bimodule_law_with_noncommuting_work_operators(self):
        operator = plus(unit(4, 0, 3), unit(4, 2, 3, (Q(0), Q(1))))
        left = plus(unit(4, 0, 1), unit(4, 2, 3))
        right = plus(unit(4, 1, 0), unit(4, 3, 2))
        self.assertEqual(r.block_readback(g.matmul(g.matmul(left, operator), right), 2),
                         g.matmul(g.matmul(left, r.block_readback(operator, 2)), right))

    def test_common_gauge_instrument_is_complete_and_preserves_old_records(self):
        family = branches()
        certificate = r.controlled_certificate(family)
        self.assertEqual(certificate["work_dimension"], 4)
        self.assertEqual(certificate["retained_outcomes"], 2)
        # Exhaust the operator basis: diagonal tests alone miss coherent errors.
        for row in range(8):
            for column in range(8):
                operator = unit(8, row, column, (Q(2, 7), Q(-3, 11)))
                outputs = r.controlled_instrument(operator, family)
                total = plus(*outputs)
                self.assertEqual(trace(total), trace(operator))
                for record in (0, 1):
                    projector = record_projector(record, 2, 4)
                    self.assertEqual(trace(g.matmul(projector, total)),
                                     trace(g.matmul(projector, operator)))

    def test_controlled_branch_cross_record_coherences_match_independent_blocks(self):
        family = branches()
        # X_(r=0,s=1)=i|2><1|. The full result is not a classical mixture.
        operator = unit(8, 2, 5, (Q(0), Q(1)))
        outputs = r.controlled_instrument(operator, family)
        for outcome, actual in enumerate(outputs):
            block = g.matmul(g.matmul(family[0][outcome], unit(4, 2, 1, (Q(0), Q(1)))),
                             g.dagger(family[1][outcome]))
            expected = tuple(tuple(block[a][b-4] if a < 4 <= b else r.ZERO
                                   for b in range(8)) for a in range(8))
            self.assertEqual(actual, expected)

    def test_record_projectors_survive_but_work_operators_need_not(self):
        unitary_family = ((g.ideal_electric(),), (g.ideal_echo(),))
        original = unit(8, 4, 4)
        after, = r.controlled_instrument(original, unitary_family)
        self.assertNotEqual(after, original)
        self.assertEqual(trace(g.matmul(record_projector(1, 2, 4), after)), r.ONE)
        self.assertEqual(r.block_readback(after, 2), after)

    def test_zero_probability_outcome_is_retained_and_no_postselection_occurs(self):
        identity_family = tuple(tuple(g.occupation_projector(0, value) for value in (0, 1))
                                for _ in (0, 1))
        output = r.controlled_instrument(unit(8, 0, 0), identity_family)
        self.assertEqual(len(output), 2)
        self.assertEqual(trace(output[0]), r.ONE)
        self.assertEqual(trace(output[1]), r.ZERO)

    def test_incomplete_or_mutated_interfaces_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "exactly complete"):
            r.controlled_certificate(((g.occupation_projector(0, 0),),))
        with self.assertRaisesRegex(ValueError, "outcome menu"):
            r.controlled_certificate((branches()[0], (g.identity(4),)))
        with self.assertRaisesRegex(ValueError, "dimensions differ"):
            r.controlled_instrument(g.identity(2), branches())
        for count in (True, 0, 3):
            with self.assertRaises(ValueError):
                r.block_readback(g.identity(4), count)
        with self.assertRaises(ValueError):
            r.native_tick_operator({(native.State((2,)+(0,)*11), native.State((2,)+(0,)*11)):
                                    (0.5, Q(0))})
        with self.assertRaisesRegex(ValueError, "not been committed"):
            state = native.State((2,)+(0,)*11)
            r.native_first_record_matrix({(state, state): r.ONE})


if __name__ == "__main__":
    unittest.main()
