"""Exact counterexamples for the consensus paper's probabilistic interfaces.

All probabilities below are rational and are derived from the displayed
transition laws, independently of simulations or serialized receipts. They
separate block-indexed contraction from calendar sampling, predictable kernel
selection from output selection, and record persistence from stationarity.
"""

from fractions import Fraction as F
from itertools import product
from pathlib import Path


def push(law, kernel):
    return tuple(sum(p * row[j] for p, row in zip(law, kernel))
                 for j in range(len(law)))


def test_state_dependent_block_lengths_violate_old_calendar_time_constant():
    # Block-endpoint chain: E[D_next | D] = D/2 + 1/4 exactly.
    block = ((F(3, 4), F(1, 4)), (F(1, 4), F(3, 4)))
    lam, epsilon = F(1, 2), F(1, 4)
    for distance, row in enumerate(block):
        assert row[1] == lam * distance + epsilon
    assert push((F(1, 2), F(1, 2)), block) == (F(1, 2), F(1, 2))

    # Unit-distance blocks last two ticks; zero-distance blocks last one.
    # States are (zero, first high tick, second high tick). The distance
    # stays constant within every block: A=1, beta=0, maximum length L=2.
    calendar = ((F(3, 4), F(1, 4), F(0)),
                (F(0), F(0), F(1)),
                (F(1, 4), F(3, 4), F(0)))
    stationary = (F(1, 3),) * 3
    assert push(stationary, calendar) == stationary
    # Irreducible and aperiodic, since all states communicate and 0 loops.
    assert calendar[0][0] > 0
    assert calendar[0][1] * calendar[1][2] * calendar[2][0] > 0
    assert stationary[1] + stationary[2] == F(2, 3)
    assert F(2, 3) > epsilon / (1 - lam)

    # A deterministic initial state also violates the purported limsup
    # constant; the stationary start is not essential to the example.
    law = (F(1), F(0), F(0))
    for _ in range(100):
        law = push(law, calendar)
    assert abs(law[1] + law[2] - F(2, 3)) < F(1, 10**8)


def test_deterministic_block_lengths_preserve_expected_tube():
    kernel = ((F(3, 4), F(1, 4)), (F(1, 4), F(3, 4)))
    for initial_distance in (0, 1):
        law = (F(1 - initial_distance), F(initial_distance))
        for block_index in range(30):
            bound = (F(1, 2) ** block_index * initial_distance
                     + (1 - F(1, 2) ** block_index) * F(1, 2))
            assert law[1] == bound
            # Every deterministic within-block tick sees the same law.
            for _ in range(2):
                assert law[1] <= bound
            law = push(law, kernel)


def test_future_selected_kernel_does_not_preserve_marginal_certificate():
    # Each fixed map sends a fair bit to a fair bit.
    output = lambda block_type, bit: bit if block_type == 0 else 1 - bit
    for block_type in (0, 1):
        assert sum(F(output(block_type, bit), 2) for bit in (0, 1)) == F(1, 2)
    # Choosing after seeing the bit instead guarantees unit output.
    selected_mean = sum(F(output(1 - bit, bit), 2) for bit in (0, 1))
    assert selected_mean == 1
    # Predictable randomized selection independent of future noise works.
    independent_mean = sum(F(output(choice, bit), 4)
                           for choice, bit in product((0, 1), repeat=2))
    assert independent_mean == F(1, 2)


def test_stationary_marginals_do_not_certify_record_persistence():
    independent = {(a, b): F(1, 4) for a, b in product((0, 1), repeat=2)}
    persistent = {(a, b): F(int(a == b), 2)
                  for a, b in product((0, 1), repeat=2)}
    for joint in (independent, persistent):
        first = tuple(sum(joint[a, b] for b in (0, 1)) for a in (0, 1))
        second = tuple(sum(joint[a, b] for a in (0, 1)) for b in (0, 1))
        assert first == second == (F(1, 2), F(1, 2))
    assert sum(p for (a, b), p in independent.items() if a != b) == F(1, 2)
    assert sum(p for (a, b), p in persistent.items() if a != b) == 0


def test_lueders_conditioning_does_not_average_a_higher_rank_block():
    # Diagonal exact density matrices suffice; P=diag(1,1,0) is central in
    # the block algebra M_2 + C. Its pure input remains pure after conditioning.
    rho, projector = (F(1), F(0), F(0)), (1, 1, 0)
    probability = sum(r * p for r, p in zip(rho, projector))
    conditional = tuple(p * r * p / probability for r, p in zip(rho, projector))
    normalized_projector = tuple(F(p, sum(projector)) for p in projector)
    assert conditional == rho
    assert conditional != normalized_projector
    assert sum(p * p for p in conditional) == 1
    assert sum(p * p for p in normalized_projector) == F(1, 2)


def test_paper_retains_the_counterexample_premise_boundaries():
    paper = (Path(__file__).resolve().parents[2]
             / "paper/reality_as_consensus_protocol.tex").read_text()
    assert r"If the block endpoints \(\tau_m\) are deterministic" in paper
    assert "certified conditional law" in paper
    assert r"\Pr\!\bigl[Y_U(t+h)\ne Y_U(t)\bigr]\le \eta" in paper
    assert r"\emph{partition-averaged} conditional state" in paper
