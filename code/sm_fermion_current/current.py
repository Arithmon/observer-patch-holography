"""Exact CAR orbital expectations and semiclassical hypercharge edge factors.

Each internal component occupies a distinct one-particle block. Stored complex
numbers are coefficients of occupied orbitals, not commuting fermion fields.
Electric feedback preserves expectation Gauss; it is not a quantum gauge field.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "oph.sm_fermion_current.v1"
STEP = F(1, 2)
HOPPING = F(1)
TARGET = 21


@dataclass(frozen=True)
class C:
    real: F = F(0)
    imag: F = F(0)

    def __post_init__(self):
        if type(self.real) not in (int, F) or type(self.imag) not in (int, F):
            raise TypeError("Gaussian rational coefficients required")
        object.__setattr__(self, "real", F(self.real))
        object.__setattr__(self, "imag", F(self.imag))

    def __add__(self, other):
        other = complexq(other)
        return C(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.real, -self.imag)

    def __sub__(self, other):
        return self + -complexq(other)

    def __rsub__(self, other):
        return complexq(other) + -self

    def __mul__(self, other):
        other = complexq(other)
        return C(self.real * other.real - self.imag * other.imag,
                 self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = complexq(other)
        denominator = other.norm()
        if not denominator:
            raise ZeroDivisionError("Gaussian rational division")
        value = self * other.conj()
        return C(value.real / denominator, value.imag / denominator)

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("integer power required")
        if exponent < 0:
            return (C(1) / self) ** (-exponent)
        answer, value = C(1), self
        while exponent:
            if exponent % 2:
                answer = answer * value
            value = value * value
            exponent //= 2
        return answer

    def conj(self):
        return C(self.real, -self.imag)

    def norm(self):
        return self.real * self.real + self.imag * self.imag


def complexq(value):
    return value if isinstance(value, C) else C(value)


I = C(0, 1)
PHASE = C(F(3, 5), F(4, 5))


def encode(value):
    if isinstance(value, C):
        return [str(value.real), str(value.imag)]
    if type(value) in (F, int):
        return str(value)
    raise TypeError("unsupported register value")


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("ascii")


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def load_parent(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def multiplets():
    parent = load_parent("code/sm_local_action/jet_action.py", "sm_current_parent_jets")
    return [{"name": name, "color_dimension_signed": color,
             "weak_dimension": weak, "multiplicity": abs(color) * weak,
             "hypercharge": str(charge), "integer_charge": int(6 * charge)}
            for name, color, weak, charge in parent.FIELDS]


def geometry():
    parent = load_parent("code/sm_abelian_reduction/pilot.py", "sm_current_parent_geometry")
    original = parent.geometry()
    return {"q": original["q"], "orbit_order": original["orbit_order"],
            "sites": [{key: row[key] for key in
                       ("index", "address_Qphi", "source_record", "source_seam_currents")}
                      for row in original["sites"]],
            "edges": [{"ends": row["ends"], "axis": row["axis"]}
                      for row in original["edges"]],
            "parent_geometry_scope": "same golden sites and oriented adjacency; scalar masses, conductances and plaquette law are not reused"}


def inner(left, right):
    return sum((a.conj() * b for a, b in zip(left, right, strict=True)), C())


def spin_action(axis, phase, spinor, adjoint=False):
    """K=i sigma_axis U^q, or its Hermitian adjoint."""
    x, y = spinor
    if axis == 0:
        result = [I * y, I * x]
    elif axis == 1:
        result = [y, -x]
    elif axis == 2:
        result = [I * x, -I * y]
    else:
        raise ValueError("axis must be 0, 1 or 2")
    multiplier = -phase.conj() if adjoint else phase
    return [multiplier * value for value in result]


def cayley(left, right, axis, phase):
    r = STEP * HOPPING / 2
    a, b = (1-r*r)/(1+r*r), -2*I*r/(1+r*r)
    kr = spin_action(axis, phase, right)
    kl = spin_action(axis, phase, left, adjoint=True)
    return ([a*x+b*y for x, y in zip(left, kr, strict=True)],
            [a*x+b*y for x, y in zip(right, kl, strict=True)])


def orbital_port(name, site, spin):
    return f"orbital/{name}/{site}/{spin}"


def spinor(values, name, site):
    return [values[orbital_port(name, site, spin)] for spin in range(2)]


def probability(values, name, site):
    return sum((z.norm() for z in spinor(values, name, site)), F(0))


def observables(values, geo, fields):
    charge, variance = [], []
    for site in range(len(geo["sites"])):
        rho, var = F(0), F(0)
        for row in fields:
            p = probability(values, row["name"], site)
            q, multiplicity = row["integer_charge"], row["multiplicity"]
            rho += multiplicity*q*p
            var += multiplicity*q*q*p*(1-p)
        charge.append(rho)
        variance.append(var)
    residual = [-rho for rho in charge]
    currents, electric = [], []
    for edge, row in enumerate(geo["edges"]):
        u, v = row["ends"]
        momentum, phase = values[f"electric/{edge}"], values[f"link/{edge}"]
        residual[u] += momentum
        residual[v] -= momentum
        current = F(0)
        for field in fields:
            q, name = field["integer_charge"], field["name"]
            ku = spin_action(row["axis"], phase**q, spinor(values, name, v))
            current -= 2*HOPPING*field["multiplicity"]*q*inner(spinor(values, name, u), ku).imag
        electric.append(momentum)
        currents.append(current)
    norms = {row["name"]: sum((probability(values, row["name"], site)
                              for site in range(len(geo["sites"]))), F(0))
             for row in fields}
    return {"charge": [str(x) for x in charge],
            "charge_variance": [str(x) for x in variance],
            "electric": [str(x) for x in electric],
            "gauss_residual": [str(x) for x in residual],
            "instantaneous_current": [str(x) for x in currents],
            "orbital_norms": {name: str(value) for name, value in norms.items()}}


def execute(geo, fields, cohort):
    if cohort not in ("baseline", "phase_intervention", "gauge_copy"):
        raise ValueError("unknown cohort")
    values, versions, writers = {}, {}, {}
    events, checkpoints = [], []
    previous = "0"*64

    def emit(operation, read_names, writes, metadata):
        nonlocal previous
        names = sorted(set(read_names))
        event = {"id": len(events), "operation": operation, "metadata": metadata,
                 "reads": [{"port": name, "version": versions[name],
                            "writer": writers[name], "value": encode(values[name])}
                           for name in names],
                 "parents": sorted({writers[name] for name in names}),
                 "writes": [], "previous_hash": previous}
        for name, value in sorted(writes.items()):
            values[name] = value
            versions[name] = versions.get(name, -1)+1
            writers[name] = event["id"]
            event["writes"].append({"port": name, "version": versions[name],
                                    "value": encode(value)})
        event["event_hash"] = digest(event)
        previous = event["event_hash"]
        events.append(event)

    for field in fields:
        for site in range(len(geo["sites"])):
            emit("prepare_orbital", [],
                 {orbital_port(field["name"], site, 0): C(F(3, 40)),
                  orbital_port(field["name"], site, 1): C(F(1, 10))},
                 {"multiplet": field["name"], "site": site})
    for edge in range(len(geo["edges"])):
        emit("prepare_link", [], {f"link/{edge}": C(1), f"electric/{edge}": F(0)},
             {"edge": edge})
    if cohort != "baseline":
        names = [orbital_port("e_c", TARGET, spin) for spin in range(2)]
        emit("phase_intervention", names, {name: PHASE*values[name] for name in names},
             {"multiplet": "e_c", "site": TARGET})
    if cohort == "gauge_copy":
        names, updates = [], {}
        phases = [PHASE**(site % 3-1) for site in range(len(geo["sites"]))]
        for field in fields:
            for site in range(len(geo["sites"])):
                for spin in range(2):
                    name = orbital_port(field["name"], site, spin)
                    names.append(name)
                    updates[name] = phases[site]**field["integer_charge"]*values[name]
        for edge, row in enumerate(geo["edges"]):
            u, v = row["ends"]
            name = f"link/{edge}"
            names.append(name)
            updates[name] = phases[u]*values[name]*phases[v].conj()
        emit("gauge_transform", names, updates,
             {"site_exponent_rule": "site_mod_3_minus_1"})

    def checkpoint(label):
        # The public record is derived from precisely these decoded registers.
        names = [name for name in values if not name.startswith("factor_current/")]
        result = observables(values, geo, fields)
        if any(value != "0" for value in result["gauss_residual"]):
            raise ArithmeticError("expectation Gauss constraint")
        if any(value != "1" for value in result["orbital_norms"].values()):
            raise ArithmeticError("occupied orbital normalization")
        emit("readout", names, {}, {"checkpoint": label})
        checkpoints.append({"event": len(events)-1, "readout": result})

    checkpoint("initial")
    for edge, row in enumerate(geo["edges"]):
        u, v = row["ends"]
        reads = [f"link/{edge}", f"electric/{edge}"]
        updates = {}
        current = F(0)
        for field in fields:
            name, q = field["name"], field["integer_charge"]
            left, right = spinor(values, name, u), spinor(values, name, v)
            phase = values[f"link/{edge}"]**q
            new_left, new_right = cayley(left, right, row["axis"], phase)
            midpoint_left = [(x+y)/2 for x, y in zip(left, new_left, strict=True)]
            midpoint_right = [(x+y)/2 for x, y in zip(right, new_right, strict=True)]
            j = -2*HOPPING*field["multiplicity"]*q*inner(
                midpoint_left, spin_action(row["axis"], phase, midpoint_right)).imag
            before = sum((z.norm() for z in left), F(0))
            after = sum((z.norm() for z in new_left), F(0))
            if field["multiplicity"]*q*(after-before) != -STEP*j:
                raise ArithmeticError("midpoint charge/current identity")
            current += j
            for site, output in ((u, new_left), (v, new_right)):
                for spin, value in enumerate(output):
                    port = orbital_port(name, site, spin)
                    reads.append(port)
                    updates[port] = value
        updates[f"electric/{edge}"] = values[f"electric/{edge}"]-STEP*current
        updates[f"factor_current/{edge}"] = current
        emit("edge_factor", reads, updates, {"edge": edge})
    checkpoint("final")
    return {"cohort": cohort, "events": events, "checkpoints": checkpoints,
            "final_event_hash": previous}


SCOPE = {
    "finite_CAR_bilinear_and_occupied_Slater_orbitals": True,
    "complete_one_generation_representation_census": True,
    "hypercharge_edge_factors_on_prepared_q5_sites": True,
    "same_midpoint_variational_current_and_electric_feedback": True,
    "exact_discrete_charge_continuity_and_expectation_Gauss": True,
    "authenticated_register_read_write_replay": True,
    "nonabelian_gauge_or_Yukawa_dynamics": False,
    "operator_Gauss_constraint_or_quantized_electric_field": False,
    "continuous_time_trajectory_error_bound": False,
    "quantum_anomaly_cancellation_from_regulator_or_measure": False,
    "source_selected_action_population_or_clock": False,
    "physical_spin_signal_order_or_laboratory_identification": False,
    "inherited_Cartan_scalar_or_Whitney_trajectory_bound": False,
}

PREMISES = [
    "Prepared q5 golden addresses, oriented nearest-neighbor graph and fixed Pauli frame.",
    "Registered all-left-Weyl multiplets with one generation; CAR and occupied Slater state are supplied.",
    "Uniform hopping k=1 and local Cayley duration dt=1/2 are declared independently of scalar conductances.",
    "Compact U1 link acts with integer charge q=6Y; its angle and electric momentum use that declared normalization.",
    "Only number-conserving hypercharge hopping factors evolve; Higgs, Yukawa, gauge-electric and plaquette factors are absent.",
    "Electric momentum is a classical mean-field coordinate driven by the CAR expectation current.",
    "Exact rational arithmetic, protected records and classical read/write access are supplied.",
]

LAW = {
    "step": str(STEP), "hopping": str(HOPPING),
    "spin_link": "K_e=i*sigma_axis*U_e^q, with Pauli axes x,y,z",
    "Hamiltonian_factor": "k*(c_u^dagger*K_e*c_v+c_v^dagger*K_e^dagger*c_u)",
    "one_particle_equation": "i*(psi_new-psi_old)/dt=H_e*(psi_new+psi_old)/2",
    "orbital_state": "one normalized occupied orbital for each distinct internal component; 15 orthogonal blocks",
    "electric_feedback": "E_new=E_old-dt*J_mid, J_mid=d_theta<E_mid> with U->exp(i*theta)*U",
    "midpoint_covariance": "sum of midpoint orbital outer products, a positive contraction; not a normalized 15-particle Slater state or a continuous-time midpoint",
    "current": "-2*k*sum_species(multiplicity*q*Im(mid_u^dagger*K_e*mid_v))",
    "charge": "sum_species(multiplicity*q*norm_squared(orbital_at_site))",
    "gauss": "outward_divergence(E)-expectation_charge",
    "charge_variance": "sum_species(multiplicity*q^2*p_site*(1-p_site)); electric field is a number",
    "schedule": "one ascending sweep over all 144 edges; each factor consumes the current endpoint orbitals",
    "initial_orbital_at_each_site": [["3/40", "0"], ["1/10", "0"]],
    "intervention": {"multiplet": "e_c", "site": TARGET, "phase": encode(PHASE)},
    "gauge_rule": "g_site=phase^(site_mod_3_minus_1); psi->g^q*psi; U_uv->g_u*U_uv*conj(g_v)",
    "operator_boundary": "nonzero local charge variance excludes operator Gauss for the declared classical E",
    "readout_current": "instantaneous current uses checkpoint orbitals; factor_current registers use that edge factor midpoint",
}


def produce():
    geo, fields = geometry(), multiplets()
    runs = [execute(geo, fields, cohort)
            for cohort in ("baseline", "phase_intervention", "gauge_copy")]
    reference, changed, gauge = [run["checkpoints"][-1]["readout"] for run in runs]
    if changed != gauge:
        raise ArithmeticError("public gauge-copy mismatch")
    difference = {key: [str(F(b)-F(a)) for a, b in zip(reference[key], changed[key], strict=True)]
                  for key in ("charge", "electric", "instantaneous_current")}
    if not any(F(value) for value in difference["electric"]):
        raise ArithmeticError("vacuous intervention")
    return {"schema": SCHEMA, "scope": SCOPE, "premises": PREMISES,
            "law": LAW, "geometry": geo, "multiplets": fields, "runs": runs,
            "comparisons": {"gauge_invariant_public_records_equal": True,
                            "intervention_minus_baseline_final": difference,
                            "initial_local_charge_variance": "945/512",
                            "generation_anomaly_sums_integer_charge": {
                                "mixed_gravity_U1": sum(f["multiplicity"]*f["integer_charge"] for f in fields),
                                "U1_cubed": sum(f["multiplicity"]*f["integer_charge"]**3 for f in fields),
                                "SU3_squared_U1_twice_Dynkin": 2*fields[0]["integer_charge"]+fields[1]["integer_charge"]+fields[2]["integer_charge"],
                                "SU2_squared_U1_twice_Dynkin": 3*fields[0]["integer_charge"]+fields[3]["integer_charge"],
                                "SU3_cubed": 2-1-1,
                                "SU2_doublet_count": 3+1,
                                "SU2_global_parity": (3+1) % 2,
                            }}}
