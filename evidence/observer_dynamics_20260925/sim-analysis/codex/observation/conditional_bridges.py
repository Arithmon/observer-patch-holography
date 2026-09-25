"""Small analytic controls for explicitly assumed observational bridges.

No parameter is fitted here. The OU field is a proposed nonconserved-noise
extension, not a result of the existing integer repair simulation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def ou_power(eigenvalue, time, diffusion=1.0, noise=1.0, *, conserved=False):
    """Variance for da=-diffusion*lambda*a dt+sqrt(2*noise)dW, a(0)=0."""
    if eigenvalue < 0 or time < 0 or diffusion <= 0 or noise < 0:
        raise ValueError("invalid OU parameters")
    if eigenvalue == 0:
        return 0.0 if conserved else 2 * noise * time
    driving = noise * eigenvalue if conserved else noise
    return driving * (-math.expm1(-2 * diffusion * eigenvalue * time)) / (diffusion * eigenvalue)


def sphere_ou_shape(ell, tau):
    """Amplitude-free C_l shape; overall D/nu is fitted separately if used."""
    if ell < 1 or int(ell) != ell:
        raise ValueError("exclude the monopole")
    return ou_power(ell * (ell + 1), tau)


def chemical_potential_energy_factor():
    # Number-conserving Bose-Einstein perturbation: Q/rho = (4a/3-b)*mu.
    zeta2, zeta3, zeta4 = math.pi**2 / 6, 1.2020569031595942854, math.pi**4 / 90
    return 1 / (4 * zeta2 / (3 * zeta3) - zeta3 / zeta4)


def independent_history_power(k, dimension, diffusion=1.0):
    """Integral_0^infinity t^(d/2-1) exp(-2Dk^2t) dt; variance weights."""
    if k <= 0 or dimension <= 0 or diffusion <= 0:
        raise ValueError("positive k, dimension and diffusion required")
    return math.gamma(dimension / 2) / (2 * diffusion * k*k)**(dimension / 2)


def finite_history_zero_mode(t_min, t_max, dimension):
    """Finite-history IR limit, showing loss of the infinite-history pole."""
    if t_min < 0 or t_max <= t_min or dimension <= 0:
        raise ValueError("invalid time interval or dimension")
    half = dimension / 2
    return (t_max**half - t_min**half) / half


def result():
    factor = chemical_potential_energy_factor()
    fraction = 16 / 55  # ensemble large-N expectation for original uniform0..5 loads
    return {
        "schema": "oph.observation.conditional-analytic-bridges.v1",
        "proposed_sphere_OU": {
            "equation": "da_lm=-nu*l(l+1)*a_lm dt+sqrt(2D)dW_lm; zero initial state and fixed/removed monopole",
            "power": "C_l(T)=D/[nu*l(l+1)]*(1-exp(-2nu*l(l+1)T))",
            "parameters": "amplitude D/nu and shape time tau=nu*T; not fitted by this producer",
            "not_current_repair_law": True,
            "templates": [{"tau": tau, "rows": [{"ell": ell, "Cl_shape": sphere_ou_shape(ell, tau)}
                                                    for ell in range(2, 41)]}
                          for tau in (0.001, 0.01, 0.1, 1.0)],
            "limits": ["At fixed finite T, lambda->0 gives 2DT, not a scale-invariant IR pole.",
                       "For each fixed positive lambda, T->infinity gives D/(nu*lambda).", 
                       "Conserved noise with mode variance proportional to lambda instead gives white stationary power."],
        },
        "local_conserved_forcing_control": {
            "equation": "dy=-nu*L*y dt+sqrt(2D)*B*dW_edges, with B*B^T=L and conserved total",
            "mode_power": "C_lambda(T)=D/nu*(1-exp(-2nu*lambda*T)) for lambda>0",
            "stationary_power": "D/nu, independent of lambda",
            "comparison_at_tau_one": [{"ell": ell,
                "nonconserved_power": sphere_ou_shape(ell, 1.0),
                "conserved_power": ou_power(ell*(ell+1), 1.0, conserved=True)}
                for ell in (1, 2, 4, 8, 16, 32)],
            "not_current_repair_law": "This linear local Gaussian control explains the effect of noise conservation; it is not an exact replacement for integer exclusion dynamics.",
        },
        "independent_scale_history_mixture": {
            "assumption": "Independent white-source histories with variance measure t^(d/2-1)dt, each diffused for age t",
            "identity": "integral_0^infinity t^(d/2-1)*exp(-2D*k^2*t)dt=Gamma(d/2)/(2D)^(d/2)*k^-d",
            "dimension_three_requires": "variance weight sqrt(t)dt, an extra assumption not supplied by current OPH repair",
            "dimension_two_requires": "constant variance weight dt, as for infinite-duration additive white forcing",
            "finite_window_IR": "bounded, tending (2/d)*(t_max^(d/2)-t_min^(d/2)); no exact IR pole at fixed finite history",
            "rows": [{"dimension": d, "k": k, "power": independent_history_power(k, d),
                      "k_to_d_times_power": k**d * independent_history_power(k, d)}
                     for d in (2, 3) for k in (0.125, 0.25, 0.5, 1, 2)],
            "nonclaims": ["Variance weights are not amplitude weights; coherently reusing the same source changes the formula.",
                         "No exponent or history distribution was inferred from the archived simulations or fitted to Planck.",
                         "Infinite history is a singular limit; physical models require finite age and UV/IR regularization."],
        },
        "conditional_FIRAS_heat_bound": {
            "primary_limit_source": "https://arxiv.org/abs/astro-ph/9605054",
            "thermal_history_reference": "https://arxiv.org/abs/1109.6552",
            "mu_per_fractional_energy_release": factor,
            "mu_era_fractional_energy_limit_approx": 9e-5 / factor,
            "y_era_fractional_energy_limit_approx": 4 * 1.5e-5,
            "uniform_0_to_5_initial_V_per_port": 55 / 6,
            "half_integer_balanced_V_per_port": 13 / 2,
            "large_N_expected_released_fraction_of_initial_V": fraction,
            "mu_era_eta_times_reservoir_to_photon_energy_limit_approx": (9e-5 / factor) / fraction,
            "y_era_eta_times_reservoir_to_photon_energy_limit_approx": (4 * 1.5e-5) / fraction,
            "assumptions": [
                "E_reservoir=epsilon*V_initial is a separately postulated physical energy identification.",
                "eta is the fraction of epsilon*DeltaV actually deposited as photon heat.",
                "The entire injection occurs in the named efficient distortion era, with standard thermalization and small distortions.",
                "The bound constrains eta*E_reservoir/E_gamma, not epsilon, observer density or cosmic time separately.",
                "The initial/terminal V fractions shown are ideal large-N values; a prediction uses the actual ledger and redshift-dependent visibility.",
                "No photon heat, redshift, or energy-unit map is measured by the existing simulation.",
            ],
        },
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    path = ROOT / "conditional_bridges.json"
    data = result()
    if args.verify:
        if json.loads(path.read_text()) != data:
            raise SystemExit("bridge receipt mismatch")
        print("verified analytic bridge controls")
    else:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        print(path)


if __name__ == "__main__":
    main()
