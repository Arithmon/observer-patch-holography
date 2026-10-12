# Testing the readouts before assigning physical units

The executable checks compare three readouts with the mathematical relations
of the specified construction. They can reject an incorrect implementation.
They do not use independent observations of nature, and they do not identify
the model's clocks, coordinates or energies with laboratory instruments.

Run the checks and retain their full inputs and outputs with:

```sh
python3 -m model.common_source.smooth_bridge --check model/common_source/smooth_results.json
python3 -m pytest -q -W error model/common_source/test_calibration.py model/common_source/test_units.py
python3 -m pytest -q -W error model/common_source/test_sanity.py
```

The combined receipt contains explicitly separate experiments. The clock
comes from the original sharp-readout history; the finite-resolution
instrument and its energy ledger use the new Ramsey history. Sharing a
receipt does not transfer the original clock law to the new detector.

## Geometric distances

The exact native rank-three projector supplies a Gram matrix. Squared
distances inferred from that matrix can be compared with the registered
three-coordinate realization, including the distances associated with native
seam currents. This checks that the algebraic slow band and the coordinate
readout agree, independently of a coordinate rotation.

Both descriptions refer to the same finite carrier. A scale in metres,
spatial translations, a signal propagation law and a growing spatial family
are additional requirements. An icosahedral carrier distance is not yet a
measured distance between physical laboratories.

## Modular durations

The actual source clock has a known binary probability law on its declared
monotone phase window. Repeated, freshly prepared trials at each programmed
setting give simulated binomial counts. Inverting the observed frequency
estimates the modular parameter; simultaneous Hoeffding intervals show the
sampling uncertainty. The seed, sample size, all settings and all outcomes
remain in the receipt.

This checks operational readability in the modular parameter. The samples
are generated from the executed source probabilities, so they test the
readout and inference procedure, not a theory against independent data.
They do not establish the duration of sequential interventions, eliminate
periodic clock aliases, or calibrate a second. The original interacting
Hamiltonian's model-time parameter is also distinct from this modular
parameter.

## Energy transfer

The smooth source history keeps the energy expectation at every preparation,
evolution, recording and feedback checkpoint. The finite Hamiltonian is the
actual compressed interacting Hamiltonian, including its boundary terms.
An independent Kraus calculation checks the emitted detector's energy
transfer. Split-evolution drift and physical apparatus energy are separate.

In the infinite canonical oscillator model the unconditional detector
increase is `1/8` for a normalized input of finite form energy. The separate
vacuum diagnostic compares finite oscillator matrices at dimensions
2, 4, 8, 16 and 32 with that limit. These are analytic matrix diagnostics;
dimensions above 4 are not executions of the common interacting source.
Small errors in this special vacuum diagnostic do not replace the full
interacting-history error certificate. No claim of an energy increment of
`1/8` on each conditioned branch is made.

## What fixing the speed of light would determine

Let one model length, model time and model energy unit correspond to `L`
metres, `T` seconds and `E` joules. **If** a common operational propagation
law establishes a dimensionless speed `c_model`, choosing the SI value of
the light speed imposes only

```text
L / T = c_SI / c_model.
```

There are two remaining freedoms among `L`, `T` and `E`. In particular,
choosing `c` does not fix an energy unit. The SI fixes
`c_SI = 299792458 m/s` exactly and defines the metre using the second;
see the [BIPM metre definition](https://www.bipm.org/en/si-base-units/metre).

If the same model additionally has an operational action constant
`hbar_model`, choosing its physical value adds

```text
E T = hbar_SI / hbar_model.
```

One overall scale remains. For example, with the illustrative choices
`c_model = hbar_model = 1`, both `T = 10^-15 s` and `T = 10^-12 s`
are compatible with the same physical `c` and `hbar`; their length units
differ by a factor of 1000 and their energy units by the reciprocal factor.
The [unit check](units.py) verifies this freedom rather than silently choosing
a physical scale. The SI value `h = 6.62607015 × 10^-34 J s` is exact;
`hbar = h/(2 pi)`. These are external inputs from the
[BIPM defining constants](https://www.bipm.org/en/measurement-units/).

The common model does **not** yet supply `c_model` connecting its native
distance readout with its modular clock. Nor does the supplied Hamiltonian
establish an energy-frequency relation for that modular clock. Consequently
the conditional conversion examples cannot calibrate the existing receipts.

## Relations that would constitute a physical test

The [sanity checks](sanity.py) also diagonalize the existing scalar family's
actual weighted Dirichlet operator on its declared unit cube. For four fixed
standing-wave labels, the largest relative frequency error against its
supplied massive dispersion law decreases as follows:

| Grid parameter | Maximum selected frequency error |
| --- | --- |
| 13 | 2.137% |
| 34 | 0.2875% |
| 89 | 0.04116% |

The target is `omega^2 = 1 + pi^2 (nx^2 + ny^2 + nz^2)`, with the inherited
mass and propagation coefficient set to one. Its high-frequency limit is
approximately linear in wavenumber when the grid resolves the mode. These
are standing-wave labels on a box, not measured momentum eigenvalues. The
frequency and oscillator energy spacing follow from a supplied quadratic
action. This verifies that existing discretization; it does not derive
Lorentz symmetry or extend the common interacting source to three spatial
dimensions.

The same check rejects one attempted clock identification in the common
source. For the actual final sharp-readout state, the local density changes
under the common Hamiltonian: derivative norms are about `0.2086` at scalar
dimension 2 and `0.6603` at 4. The modular flow of that same local reference
state preserves it. A norm-one local observable witnesses the difference,
and direct finite differences verify the Hamiltonian derivative. Therefore
no scalar rescaling makes these two flows identical on all local
expectations in this state. This is a numerical obstruction to that specific
identification, not a rejection of other observer clocks or of OPH.

A useful next physical attachment must measure propagation, frequency,
energy and spatial translations on one family of observation histories.
For a relativistic sector, one could then calibrate one propagation speed
and test the dispersion relation over other momenta, including comparison
of independently reconstructed sectors. Dimensionless speed ratios, mass
ratios and detector probabilities cannot be repaired by changing metres,
seconds or joules. They provide meaningful rejection tests once the model
and observational identification are specified.

Recovering these operational relations would meet the intended requirement;
fundamental spacetime or fields need not be postulated. Assigning SI names
to separate parameters does not establish those relations.
