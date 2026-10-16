# Identifiability and calibration checks

These executable controls establish what the supplied action, detector and
geometry benchmarks can identify. They retain exact counterexamples,
independent finite-matrix checks and a retrospective empirical comparator.
They do not select a physical action or establish OPH-specific empirical
agreement. The [source-selection guide](../SOURCE_SELECTION.md) states the
corresponding physical obligations.

## Action and geometry ambiguities

The compressed four-level calibration is exactly insensitive to

```text
B(N) = N(N-1)(N-2)(N-3)/24.
```

Adding this local positive term leaves every operation inside that calibrated
space unchanged but changes its eight-level extension. Greater precision on
the same calibration cannot remove the ambiguity. The
[identification experiment](../identification_design/README.md) addresses it
with a coherent preparation involving level four and a phase-sensitive read;
its declared two candidate laws differ by total variation `0.7071067812`.
The action family and coherent services remain supplied.

Likewise, the unit-edge cube graph admits a planar embedding with the same
existing graph dynamics and observation records. Inferring rank three needs
the declared orthogonal-cube metric convention or additional metric-sensitive
observations. Connectivity rank and cycle count alone are not a physical
spatial dimension.

The scalar sampling grids retain the same seven modes. Increasing their
sampling density changes a gradient coefficient but does not enlarge the
field-mode family. The exact-gradient endpoint is appropriate for fixed-band
comparisons; omitted modes require a separate enlargement and shared
observable domain.

[architecture_counterexamples.py](architecture_counterexamples.py) executes
these controls; [architecture_results.json](architecture_results.json)
retains their exact statements and numerical contrasts.

## Independent empirical comparator

The [anharmonic transfer](../anharmonic_transfer/README.md) tests the supplied
spin/pulse model against 125 positive-duration position-estimator values from
[McGarry et al.](https://arxiv.org/html/2603.04744v1), using their
[public deposition](https://zenodo.org/records/18752287). An independent,
ancilla-free comparison uses

```text
H = 2*pi*500*N + 4000*cos[(pi/3)*(a+a_dagger)+potential_phase],
L = sqrt(18)*N.
```

It uses the same transferred preparation, potential phases, literal `h=.4`
detector and published simulated-duration coordinate, with no additional
fit. DOP853 integrates the finite Lindblad equation without normalizing or
projecting the state. This is a supplied effective-Hamiltonian comparator;
it does not execute native source words.

Its aggregate estimator RMSE is **0.317681**, compared with **0.287643** for
the detailed pulse model: an improvement of **9.46%** for the latter. The
pulse model improves four of the five cases; the simpler model is better
in case e. The zero-prediction control has RMSE **0.645911**. The comparison
is retrospective, with observed responses available when it was designed.
It is neither a reserved-data test nor a distinction between OPH and
ordinary quantum mechanics.

The 32-versus-48-level target calculation changes an estimator by at most
`1.62e-7`; trace residuals are below `1.7e-15`. These are numerical checks,
not unrestricted-occupation certificates or apparatus-error bounds. The
[occupation certificate](../occupation_control/README.md) supplies a separate
controlled comparison for the specified pulse model.

The empirical inputs and complete output are retained in
[empirical_probe.json](empirical_probe.json). The compressed
[transfer reference](data/transfer_reference.json.gz) fixes the preparation
and comparison payload, independently of a live receipt regeneration. Its
decompressed SHA-256 is
`71403f474c39705d20a6c30dc2f36b7556943dd4f1734a5a1a395184744d2398`.
The manifest binds the compressed bytes; the producer validates the
uncompressed digest before parsing. No empirical values are modified.

## Measurement provenance and interpretation

All 130 deposited position estimates, including the five zero-duration
values, satisfy

```text
|63.6*x - round(63.6*x)| < 3.13e-8.
```

This printed quantization is compatible with hundredth-sized signals divided
by `.636`, whereas the documented conversion uses `sqrt(2)*.4`. The pattern
does not authenticate a detector gain, displacement coordinate, accepted
trial count or extraction procedure. It provides no basis for rescaling the
observations after inspecting residuals. The published `h=.4` comparison is
retained with its uncertainty explicitly unresolved.

Preparation transfer also remains an empirical assumption: the largest
zero-duration estimator residual is **0.340550**. The largest positive-duration
residual is about **1.0823**, exceeding the selected ideal finite-model radius
of **0.001517**. That radius does not include apparatus uncertainty or removal
of the occupation cutoff. Refining the time splitter alone therefore cannot
settle the observed discrepancy.

Selection dwell, detector-pulse durations, phase tracking and preparation
validity need independent laboratory metadata. The current comparison uses
the stated ideal selection and readout interfaces. Aggregate means without
authenticated accepted counts do not support a statistical rejection.

For normalized characteristic functions, Hermiticity requires
`chi(-beta)=conjugate(chi(beta))`. For a real-displaced thermal state,
`chi(u+iv)=exp[-(nbar+1/2)(u²+v²)+2i*alpha*v]`: conjugating beta conjugates
the characteristic value. Exact reflected copies in deposited grids are
not independent observations. The [readout controls](../readout_controls/README.md)
and [two-order Gram constraints](../characteristic_constraints/README.md)
apply these identities with their stated detector and preparation premises.

## Reused mathematical results

The [finite environment derivation](../anharmonic_transfer/DERIVATION_ENVIRONMENT.md)
uses positive Gaussian-phase quadrature and CPTP contraction. The
[splitting and conditioning bounds](../anharmonic_transfer/DERIVATION_BOUNDS.md)
retain the finite oscillator commutator
`[a,a_dagger]=I-d|d-1><d-1|`. Replacing it by the unrestricted identity
invalidates a uniform finite-matrix bound. Complete accepted/rejected outputs
give conditional TV at most `epsilon/q`, where `q` is the reference
acceptance; small acceptance amplifies error.

[test_math_probes.py](test_math_probes.py) checks these algebraic ingredients,
the exact cutoff-boundary residual, finite forced-oscillator counterexamples,
and complete-source number/parity word reductions. The implemented source
words are documented in the apparatus derivations. The
[Etienney–Robin–Rouchon occupation theorem](https://arxiv.org/html/2501.09607v6)
is attached with a controlled trajectory integral in the
[occupation derivation](../occupation_control/DERIVATION.md). Agreement
between two cutoffs alone does not substitute for that theorem.

Finite instrument realization reuses
[Stinespring dilation](https://community.ams.org/journals/proc/1955-006-02/S0002-9939-1955-0069403-4/S0002-9939-1955-0069403-4.pdf)
and [universal finite gates](https://arxiv.org/abs/quant-ph/9503016), conditional
on admitted coherent services. Compilation shows that the chosen law can be
implemented; it does not choose that law over other compilable laws.
The [quantum-network framework](https://arxiv.org/abs/0904.4483) and
[operational Markov criterion](https://arxiv.org/abs/1801.09811) require
comparisons under the same allowed interventions, including surviving memory
and reference outputs. The [joint clock construction](../joint_clock_process/README.md)
provides an explicit native example where a reduced instantaneous field state
is insufficient but its history-dependent description is exact.

Relativistic reconstruction separately requires a compatible correlation
hierarchy with the hypotheses of the corrected
[Osterwalder–Schrader theorem](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-42/issue-3/Axioms-for-Euclidean-Greens-functions-II-with-an-Appendix-by/cmp/1103899050.pdf).
A finite GNS construction or positive transfer matrix does not establish
those hypotheses. The known [four-dimensional Ising/phi4 scaling result](https://arxiv.org/abs/1912.07973)
also constrains the specified conventional scalar families; it is neither
a general exclusion of OPH nor a proof of a nontrivial continuum limit.
See [observer-relative recovery](../OBSERVER_RELATIVE_RECOVERY.md) for the
current finite construction and remaining compatibility obligations.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identifiability_checks.architecture_counterexamples
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identifiability_checks.empirical_probe
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/identifiability_checks
```

The tests retain exact architecture comparisons, source/data bindings,
corruption controls and the complete empirical replay. Numerical diagnostics,
ideal mathematical bounds and physical uncertainty have separate meanings.
