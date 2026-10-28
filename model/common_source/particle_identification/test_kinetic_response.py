"""Independent free-action, finite-difference and generalized-eigenvalue checks."""
from dataclasses import replace
import hashlib
from pathlib import Path, PureWindowsPath

import numpy as np
import pytest
from scipy.linalg import eigvalsh

from ..mass_response import model
from . import kinetic_response as k, response


def test_report_source_keys_are_portable_under_windows_paths(monkeypatch):
    relative_to = Path.relative_to
    with monkeypatch.context() as patch:
        patch.setattr(Path, "relative_to", lambda self, *args, **kwargs:
                      PureWindowsPath(relative_to(self, *args, **kwargs)))
        actual = k.report()["source_sha256"]
    root = Path(k.__file__).resolve().parents[3]
    names = ("model/common_source/particle_identification/kinetic_response.py",
             "model/common_source/particle_identification/response.py",
             "model/common_source/mass_response/model.py")
    assert actual == {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in names}


def test_free_inverse_matches_one_particle_action_and_not_neutral_density():
    packet = response.finite_response(interaction=False)
    z = .23+.31j
    jet = k.inverse_jet(packet, z)
    sites = [1 << i for i in range(4)]
    action = model.ring_hamiltonian(interaction=False)[np.ix_(sites, sites)]
    np.testing.assert_allclose(jet["inverse"], z*np.eye(4)-action, atol=2e-14)
    np.testing.assert_allclose(jet["inverse_frequency_derivative"], np.eye(4), atol=3e-14)


@pytest.mark.parametrize("interaction", [False, True])
def test_inverse_derivative_matches_independent_finite_difference(interaction):
    packet = response.finite_response(interaction=interaction)
    z, step = .31+.23j, 1e-5
    analytic = k.inverse_jet(packet, z)["inverse_frequency_derivative"]
    numeric = (np.linalg.inv(response.green_function(packet, z+step))
               -np.linalg.inv(response.green_function(packet, z-step)))/(2*step)
    np.testing.assert_allclose(analytic, numeric, atol=2e-9, rtol=0)


def test_nonunitary_coordinates_need_generalized_kinetics_to_preserve_free_poles():
    packet = response.finite_response(interaction=False)
    z = .31+.23j
    a = np.array([[2, .2j, 0, 0], [0, .7, .1, 0], [0, 0, 1.3, -.15j], [.05, 0, 0, .9]])
    changed = k.inverse_jet(packet, z, coordinates=a)
    actual_sum = sum(pole.residue for pole in (*packet.addition, *packet.removal))
    np.testing.assert_array_equal(changed["CAR_integrated_residue"], a@actual_sum@a.conj().T)
    np.testing.assert_array_equal(changed["expected_canonical_CAR_Gram"], a@a.conj().T)
    kinetic = changed["inverse_frequency_derivative"]
    mass = z*kinetic-changed["inverse"]
    # The free pencil is z Z-M; changing coordinates alone cannot move its roots.
    np.testing.assert_allclose(kinetic, np.linalg.inv(a@a.conj().T), atol=2e-14)
    sites = [1 << i for i in range(4)]
    expected = np.linalg.eigvalsh(model.ring_hamiltonian(interaction=False)[np.ix_(sites, sites)])
    np.testing.assert_allclose(eigvalsh(mass, kinetic), expected, atol=2e-14)
    assert np.max(abs(np.linalg.eigvalsh(mass)-expected)) > .1


def test_interacting_frequency_dependence_is_not_claimed_as_static_mass():
    packet = response.finite_response()
    one = k.inverse_jet(packet, .31+.23j)
    two = k.inverse_jet(packet, -.19+.37j)
    assert np.max(abs(one["inverse_frequency_derivative"]-two["inverse_frequency_derivative"])) > 1e-3
    report = k.report()
    assert not report["response_inferred_from_observer_records"]
    assert not report["positive_physical_kinetic_metric_established"]
    assert not report["physical_mass_or_species_identification"]


@pytest.mark.parametrize("coordinates", [np.eye(3), np.zeros((4, 4)), np.eye(4)*np.nan,
                                        np.ma.array(np.eye(4), mask=False)])
def test_incomplete_or_noninvertible_coordinate_maps_are_refused(coordinates):
    with pytest.raises(ValueError):
        k.inverse_jet(response.finite_response(), .2+.3j, coordinates=coordinates)


@pytest.mark.parametrize("damage", ["missing", "doubled", "negative", "nonhermitian", "nan_frequency", "wrong_sign"])
def test_tampered_spectral_packets_cannot_claim_canonical_CAR(damage):
    packet = response.finite_response()
    first = packet.addition[0]
    if damage == "missing":
        packet = replace(packet, addition=packet.addition[1:])
    else:
        edits = {
            "doubled": {"residue": 2*first.residue},
            "negative": {"residue": -first.residue},
            "nonhermitian": {"residue": first.residue+1j*np.eye(4)},
            "nan_frequency": {"frequency": float("nan")},
            "wrong_sign": {"frequency": -first.frequency},
        }
        packet = replace(packet, addition=(replace(first, **edits[damage]), *packet.addition[1:]))
    with pytest.raises(ValueError):
        k.inverse_jet(packet, .2+.3j)


@pytest.mark.parametrize("frequency", [0., .2-.1j, complex(float("nan"), .2), .2+complex(0., float("inf")), True])
def test_frequency_keeps_the_retarded_upper_half_plane_contract(frequency):
    with pytest.raises(ValueError):
        k.inverse_jet(response.finite_response(), frequency)
