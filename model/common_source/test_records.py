"""Operator-level checks of the checkpoint adapter and its creation boundary."""
from fractions import Fraction as Q
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
    def test_native_all_input_gram_and_intertwiner_certificate(self):
        receipt = r.native_checkpoint_certificate()
        self.assertEqual(receipt["all_phase_column_counts"],
                         {"0": 78, "1": 468, "2": 2808})
        self.assertEqual(receipt["source_tape_dimension"], 3354)
        self.assertEqual(receipt["postcommit_dimension"], 3276)
        self.assertTrue(receipt["postcommit_dual_projectors_fixed"])
        self.assertFalse(receipt["full_operational_observer_inhabitant_claimed"])

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
