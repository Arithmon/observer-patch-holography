"""Source-compiled Palm modular flow on the entire native M66 algebra.

The native reference has sigma=diag(degree/600), hence
sigma^(-is)=exp(is log(60))*exp(is log(5/4) P_degree8). Dropping the
common native scalar fixes the target gauge. All thirty degree-eight
addresses receive the phase; the other 482 encoded addresses receive one.
The 446 nonbinary/padding addresses specify an extension, not a modular
state there. Intermediate source words retain all 512 addresses.

A closed Gray walk on each greatest-bit parity computes the complete
nine-bit Walsh expansion with the existing phase-exact CNOT, Z and I source
requests. It emits 512 phase requests and 510 three-pulse CNOTs. No arbitrary
matrix gate, new scalar-phase service, or projection onto the native shell
is inserted. The underlying proper-code load/export and coherent controls
remain the existing admitted services.

For each exact rational |s|<=16, a rational atanh series encloses log(5/4).
A degree-64 Taylor polynomial and its unitary integral remainder enclose
the target phase. Every original binary component of the executed 512x512
map is then compared with that enclosure. eps>=||A-U|| implies half-diamond
error eps*(2+eps)/2 and norm gain <=(1+eps)^2, with arbitrary references.
Each requested word is checked; an unexecuted numerical word is not certified.

The finite preparation theorem gives t=66 delta_n and modular stability
b=min(1,|s|t/(1-t)) on the successful native work algebra. Triangle inequality
therefore bounds the compiled channel against the actual prepared state's
modular channel by delta_word+b, without constructing that state's logarithm.
This compares channels on native probes, including entangled references;
it does not claim faithfulness of a flag algebra or of the padding extension.
The complete herald experiment also obeys error<=delta_n+p_max*delta_word
against its ideal Palm-output experiment, retaining its failure branch.

The full-word certificate is a terminal logical-channel certificate. Native
proper-code export failures vanish by the existing exact invariant-code
identities; physical service errors and floating state propagation remain
separate. Midword measurements are not part of the declared operation.
The parameter s is dimensionless and is not passive repair or physical time.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import construction as corner, preparation

compiled = corner.compiled
compiler = compiled.compiler
exact = corner.exact
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("full_flow_results.json")
DIMENSION = 512
WORD_ERROR_GATE = F(1, 10**9)
MIXER_ERROR_GATE = F(1, 10**11)
PARAMETERS = (F(0), F(1), F(16))


def parameter(value):
    if type(value) not in (int, F) or abs(value) > 16:
        raise ValueError("an exact rational modular parameter in [-16,16] required")
    return F(value)


@lru_cache(maxsize=1)
def support():
    degree, sigma = preparation.anchor()
    codes = tuple(compiled.circuit.encode(state) for state in corner.quantum.STATES)
    mask = tuple(sorted(code for code, d in zip(codes, degree, strict=True) if d == 8))
    if len(set(codes)) != 66 or len(mask) != 30 or max(codes) >= DIMENSION:
        raise ArithmeticError("full native Palm encoding changed")
    if tuple(F(d, 600) for d in degree) != tuple(sigma):
        raise ArithmeticError("full modular spectrum is not the native Palm spectrum")
    return codes, mask


@lru_cache(maxsize=1)
def walsh_counts():
    mask = support()[1]
    return tuple(sum((-1)**((code & parity).bit_count()) for code in mask)
                 for parity in range(DIMENSION))


def support_certificate():
    codes, mask = support()
    counts = walsh_counts()
    for code in range(DIMENSION):
        value = sum(count*(-1)**((code & parity).bit_count())
                    for parity, count in enumerate(counts))
        if value != DIMENSION*int(code in mask):
            raise ArithmeticError("complete Walsh expansion lost a native or padding address")
    old = set(corner.reference()[3])
    if len(set(mask)-old) != 22:
        raise ArithmeticError("old active-corner omission changed")
    return {"native_dimension": 66, "register_dimension": DIMENSION,
            "degree8_addresses": list(mask), "degree10_addresses": len(codes)-len(mask),
            "previously_omitted_degree8_addresses": sorted(set(mask)-old),
            "identity_extension_addresses": DIMENSION-len(codes),
            "walsh_coordinates_verified": DIMENSION,
            "intermediate_padding_coordinates_retained": True}


@lru_cache(maxsize=1)
def logarithm_interval():
    """log(5/4)=2 atanh(1/9), with a positive rational tail bound."""
    x, terms = F(1, 9), 32
    lower = 2*sum((x**(2*j+1)/(2*j+1) for j in range(terms)), F())
    tail = 2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
    return lower, lower+tail


def phase_enclosure(s):
    """A Gaussian-rational center and an absolute complex phase error."""
    s = parameter(s)
    lo, hi = logarithm_interval()
    theta = s*(lo+hi)/2
    real, imag, power_real, power_imag = F(1), F(0), F(1), F(0)
    for k in range(1, 65):
        power_real, power_imag = -theta*power_imag/k, theta*power_real/k
        real += power_real; imag += power_imag
    radius = abs(theta)**65/math.factorial(65)+abs(s)*(hi-lo)/2
    return real, imag, radius


def word(s):
    s = parameter(s)
    lo, hi = logarithm_interval()
    theta = float(s*(lo+hi)/2)
    counts = walsh_counts()
    tape = [compiler.pulse("I", theta*counts[0]/DIMENSION, [0])]
    for pivot in range(9):
        count = 1 << pivot
        for index in range(count):
            gray = index ^ (index >> 1)
            mask = (1 << pivot) | gray
            tape.append(compiler.pulse("Z", theta*counts[mask]/DIMENSION, [pivot]))
            following = (index+1) % count
            change = gray ^ (following ^ (following >> 1))
            if change:
                tape.extend(compiler.cnot(change.bit_length()-1, pivot))
    if len(tape) != 2042:
        raise ArithmeticError("full native modular source schedule changed")
    return tape


@lru_cache(maxsize=8)
def _executed(s):
    array = compiler.execute(word(s), 9)
    return array.dtype.str, array.tobytes()


def matrix(s):
    dtype, payload = _executed(parameter(s))
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(DIMENSION, DIMENSION)


def certify(s, array=None):
    s = parameter(s)
    array = matrix(s) if array is None else array
    if type(array) is not np.ndarray:
        raise ValueError("an original base ndarray required before dispatch")
    array = exact._array(array, (DIMENSION, DIMENSION))
    real, imag, radius = phase_enclosure(s)
    marked = set(support()[1])
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            ar, ai = exact._component(value)
            wr, wi = ((real, imag) if i in marked else (F(1), F())) if i == j else (F(), F())
            squared += (ar-wr)**2+(ai-wi)**2
    eps = exact._sqrt_interval(squared)[1]+exact._sqrt_interval(F(30))[1]*radius
    return {"operator_error_upper": eps,
            "half_diamond_upper": eps*(2+eps)/2,
            "channel_norm_upper": (1+eps)**2,
            "phase_enclosure_error_upper": radius}


@lru_cache(maxsize=8)
def _check_source_word(s):
    actual = certify(s)
    corner.report._gate(actual["operator_error_upper"], WORD_ERROR_GATE,
                        "full native modular source word")


def checked_bound(s):
    _check_source_word(parameter(s))
    return {"operator_error_upper": WORD_ERROR_GATE,
            "half_diamond_upper": WORD_ERROR_GATE*(2+WORD_ERROR_GATE)/2,
            "channel_norm_upper": (1+WORD_ERROR_GATE)**2}


def composed_bound(s, attempts=369):
    """Every successful work reference; every native input/probe reference."""
    s = parameter(s)
    prepared = preparation.budget(attempts, abs(s))
    stability = prepared["modular_channel_half_diamond_error_upper"]
    if stability is None:
        raise ValueError("finite preparation does not certify a faithful work reference")
    source = checked_bound(s)
    return {"warmup_native_attempts": attempts,
            "herald_native_attempts": 1,
            "raw_source_half_diamond_upper": source["half_diamond_upper"],
            "prepared_reference_modular_half_diamond_upper": stability,
            "source_to_prepared_modular_half_diamond_upper": source["half_diamond_upper"]+stability,
            "complete_flagged_preparation_and_source_half_trace_upper":
                prepared["warmup_half_trace_error_upper"]
                +prepared["herald_success_interval"][1]*source["half_diamond_upper"],
            "herald_success_interval": prepared["herald_success_interval"],
            "expected_native_preparation_attempts_upper":
                prepared["expected_attempts_with_full_warmup_repeated_each_trial_upper"]}


def mixer_word(sign):
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("an integer mixer sign required")
    return [compiler.pulse("Y", sign*math.pi/4, [7])]


def mixer_bound(sign):
    """All-input two-dimensional bound, invariant under spectator/reference lift."""
    mixer_word(sign)  # Validate before compiling the spectator-free certificate.
    actual = compiler.execute([compiler.pulse("Y", sign*math.pi/4, [0])], 1)
    lo, hi = exact._sqrt_interval(F(1, 2))
    middle, radius = (lo+hi)/2, (hi-lo)/2
    target = ((middle, sign*middle), (-sign*middle, middle))
    squared = F()
    for i in range(2):
        for j in range(2):
            real, imag = exact._component(actual[i, j])
            squared += (real-target[i][j])**2+imag*imag
    eps = exact._sqrt_interval(squared)[1]+2*radius
    corner.report._gate(eps, MIXER_ERROR_GATE, "native omitted-pair Ramsey mixer")
    return {"half_diamond_upper": MIXER_ERROR_GATE*(2+MIXER_ERROR_GATE)/2,
            "channel_norm_upper": (1+MIXER_ERROR_GATE)**2}


def witness(s):
    """Actual source Ramsey execution on an omitted native degree-eight state.

    The initial native configuration is code16. The admitted Y service makes
    its coherent superposition with code144; both states are outside the old
    active pair experiment. The final admitted work-register bit measurement
    retains all 512 coordinates in its two outcomes. Nonnative support is a
    numerical diagnostic, not a separately measured membership projector.
    Initialization and projective readout are specified apparatus; their
    physical costs are not included in the source pulse count.
    """
    s = parameter(s)
    codes, mask = support()
    if 16 not in mask or 144 not in codes or 144 in mask:
        raise ArithmeticError("the independent omitted-pair witness changed")
    initial = np.zeros(DIMENSION, complex); initial[16] = 1
    prepared = compiler.execute(mixer_word(-1), 9, initial)
    phased = matrix(s)@prepared
    analyzed = compiler.execute(mixer_word(1), 9, phased)
    groups = (tuple(code for code in range(DIMENSION) if not (code >> 7) & 1),
              tuple(code for code in range(DIMENSION) if (code >> 7) & 1))
    vectors = []
    for group in groups:
        vector = np.zeros(DIMENSION, complex)
        vector[list(group)] = analyzed[list(group)]
        vectors.append(vector)
    probabilities = [float(np.vdot(vector, vector).real) for vector in vectors]
    real, _, radius = phase_enclosure(s)
    ideal_interval = ((1+real-radius)/2, (1+real+radius)/2)
    delta = corner.composite_bound(checked_bound(s), mixer_bound(1))
    mixer_bound(-1)
    target_probability = float((1+real)/2)
    half_l1 = (abs(probabilities[0]-target_probability)
               +abs(probabilities[1]-(1-target_probability)))/2
    # Raw numerical maps need not preserve trace exactly: half trace norm
    # controls the complete record's half L1, not each individual probability.
    if half_l1 > float(delta)+float(radius)/2:
        raise ArithmeticError("source Ramsey record exceeds its complete-word comparison")
    residual = abs(sum(probabilities)-float(np.vdot(analyzed, analyzed).real))
    if residual > 2e-15:
        raise ArithmeticError("register-bit pointer discarded terminal probability")
    outside = sum(abs(analyzed[code])**2 for code in range(DIMENSION) if code not in set(codes))
    return {"initial_native_code": 16, "coherent_partner_native_code": 144,
            "pointer_outcomes": ["register_bit7_zero", "register_bit7_one"],
            "ideal_match_probability_interval": [compiled._outward(ideal_interval[0], upper=False),
                                                  compiled._outward(ideal_interval[1])],
            "executed_complete_probabilities": probabilities,
            "complete_record_half_l1_error_diagnostic": half_l1,
            "raw_total_probability_diagnostic": sum(probabilities),
            "outside_native_support_probability_diagnostic": float(outside),
            "whole_Ramsey_half_diamond_upper": corner.compiled._outward(delta),
            "native_modular_source_pulses": len(word(s)),
            "Ramsey_mixer_source_pulses": 2,
            "declared_configuration_and_coherent_services": True,
            "declared_projective_register_bit_readout": True,
            "initialization_and_readout_physical_costs_included": False,
            "raw_output_renormalized": False}


def run():
    support_row = support_certificate()
    rows = []
    for s in PARAMETERS:
        bound = composed_bound(s)
        rows.append({"modular_parameter": str(s),
                     "source_word_sha256": hashlib.sha256(corner.canonical(word(s)).encode()).hexdigest(),
                     "resources": compiler.resources(word(s)),
                     "finite_error_budget": preparation._compact(bound),
                     "omitted_pair_Ramsey": witness(s)})
    paths = sorted(set((*corner.SOURCE_PATHS, "model/native_modular_reference/preparation.py",
                        "model/native_modular_reference/full_flow.py",
                        "model/record_laws/burn_in.py", "Lean/QFT/StructuralInheritance.lean")))
    hashes = {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}
    return {"schema": "oph.native_modular_reference.full_flow.v1", "all_pass": True,
            "source_dependency_count": len(paths),
            "source_dependency_map_sha256": hashlib.sha256(corner.canonical(hashes).encode()).hexdigest(),
            "support": support_row, "executions": rows,
            "clock": {"reference": "diag(degree/600)", "frequency": "log(5/4)",
                      "Heisenberg": "sigma^(is) A sigma^(-is)",
                      "state_action_up_to_native_scalar": "exp(i s log(5/4) P_degree8)",
                      "parameter_domain": "exact rational |s|<=16; each source word freshly certified",
                      "source_word_gate": str(WORD_ERROR_GATE)},
            "scope": {"full_M66_modular_channel_implemented": True,
                      "all_512_native_and_padding_input_columns_certified": True,
                      "arbitrary_native_input_reference_entanglement_covered": True,
                      "actual_prepared_state_logarithm_controller_needed_for_this_error_bound": False,
                      "coherent_load_export_and_projective_bit_read_services_admitted": True,
                      "physical_proper_code_service_errors_bounded": False,
                      "floating_state_propagation_is_a_directed_certificate": False,
                      "intermediate_pulse_readouts_covered_by_terminal_certificate": False,
                      "native_passive_repair_equals_modular_clock": False,
                      "physical_seconds_or_temperature_identified": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        corner.compare_receipt(strict_json.load(args.check), result)
        print("Full native Palm modular flow and finite preparation bound verified.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(corner.canonical(result)+"\n")
    else:
        print(corner.canonical(result))


if __name__ == "__main__":
    main()
