"""Exact tiny-chain check for the stationary repair autocorrelation boundary."""
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def balanced_transition(n, raised, edges, swap_probability=F(1, 2)):
    states = [tuple(int(i in selected) for i in range(n)) for selected in itertools.combinations(range(n), raised)]
    where = {state: i for i, state in enumerate(states)}
    matrix = [[F(0) for _ in states] for _ in states]
    for index, state in enumerate(states):
        for i, j in edges:
            changed = list(state); changed[i], changed[j] = changed[j], changed[i]
            matrix[index][index] += (1 - swap_probability) / len(edges)
            matrix[index][where[tuple(changed)]] += swap_probability / len(edges)
    return states, matrix


def autocorrelations(matrix, observable, lags):
    size = len(matrix)
    mean = sum(observable) / F(size)
    centered = [F(v) - mean for v in observable]
    current = list(centered)
    answer = []
    for _ in range(lags + 1):
        answer.append(sum(a * b for a, b in zip(centered, current)) / size)
        current = [sum(p * v for p, v in zip(row, current)) for row in matrix]
    return answer


def build():
    states, matrix = balanced_transition(4, 2, [(0, 1), (1, 2), (2, 3)])
    assert all(sum(row) == 1 for row in matrix)
    assert all(matrix[i][j] == matrix[j][i] for i in range(6) for j in range(6))
    assert all(matrix[i][i] >= F(1, 2) for i in range(6))
    # K=2P-I is stochastic and symmetric; consequently spec(P) is in[0,1].
    k = [[2 * p - int(i == j) for j, p in enumerate(row)] for i, row in enumerate(matrix)]
    assert all(min(row) >= 0 and sum(row) == 1 for row in k)
    observable = [state[0] - state[-1] for state in states]
    correlations = autocorrelations(matrix, observable, 20)
    assert all(v >= 0 for v in correlations)
    assert all(a >= b for a, b in zip(correlations, correlations[1:]))
    _, nonlazy = balanced_transition(2, 1, [(0, 1)], F(1))
    adversarial = autocorrelations(nonlazy, [1, -1], 8)
    assert adversarial == [F((-1) ** t) for t in range(9)]
    return {"schema": "oph.codex-observables.stationary-control.v1",
            "spec_sha256": hashlib.sha256((HERE / "stationary_spec.json").read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "states": states, "transition_matrix_rationals": [[str(p) for p in row] for row in matrix],
            "stationary_measure": "uniform1/6", "observable": observable,
            "autocorrelation_rationals": [str(v) for v in correlations],
            "autocorrelation_float": [float(v) for v in correlations],
            "nonlazy_adversarial_autocorrelation": [str(v) for v in adversarial],
            "exact_checks": {"stochastic": True, "symmetric": True, "diagonal_at_least_half": True,
                             "twice_transition_minus_identity_stochastic": True,
                             "autocorrelation_nonnegative_nonincreasing": True},
            "general_argument": "P=(I+K)/2 with K symmetric stochastic gives spec(P) subset[0,1]; spectral expansion C_f(t)=sum_j a_j^2 lambda_j^t is nonnegative and nonincreasing."}


if __name__ == "__main__":
    result = build()
    (HERE / "stationary_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
    print("exact stationary autocorrelation and nonlazy adversarial controls passed")
