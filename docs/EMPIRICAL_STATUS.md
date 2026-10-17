# Empirical status and repair-rule comparison

**No retained laboratory comparison selects native repair rule C or E.**
The repository contains reproducible comparisons with physical measurements,
conditional postdictions, hardware implementation tests and executable designs
for distinguishing repair rules. These have different evidential roles. None
establishes that the executable observer reproduces the full range of physical
observations or supplies an OPH-specific empirical confirmation.

The target is agreement between observers' measurement statistics and their
effective interpretations. It does not require a fundamental spacetime object.
It does require identifying the same preparation, intervention and readout in
the model and the experiment before comparing their probabilities.

## How many comparisons?

There are **three broad external comparison families in the executable-model
work** if the fine-structure reference is included:

1. Oscillator observations from the [McGarry et al. deposition](https://doi.org/10.5281/zenodo.18752287).
   Gaussian readout, driven dynamics, changed-potential transfer and their
   measurement diagnostics reuse that release. They are separate analyses,
   not independent experimental datasets.
2. The inherited coupling endpoint against the
   [NIST CODATA 2022 constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).
   This is a reference-value comparison, not an observer-trial dataset.
3. Bell observations from the [Hensen et al. deposition](https://doi.org/10.4121/uuid:6e19e9b2-4a2d-40b5-8dd3-a660bf3c0a31).
   This is a second independent experimental deposition. The replay tests a
   local classical response interpretation, not the full quantum observer.

Thus “three families” is a useful inventory label; “three successful,
independent tests of the native observer” is not. The wider OPH repository also
contains the particle, cosmology and quantum-hardware comparisons below.
Reusing records, overlapping cosmological likelihoods, or reproducing a fitted
reference does not multiply independent evidence.

## Executable-model comparisons

The values below are retained analysis outputs. RMSEs in different measurement
coordinates cannot be ranked against one another. A small numerical simulation
error is not a bound on physical apparatus error.

| Comparison and evidence | Observations and result | Supplied inputs and interpretation | Selects C/E? |
| --- | --- | --- | --- |
| [Source-executed Gaussian readout](../model/scalar_observations/experimental_readout/EXECUTABLE_READOUT.md) | 21 complex calibration probes determine two preparation parameters; 100 reserved complex probes give 200 spin means. Primary spin-mean RMSE `0.06647166`; alternative coordinate extent gives `0.07178838`. | Supplied oscillator, displaced-thermal preparation family and detector; retrospective calibration/evaluation split. Coordinate conventions and experimental postselection remain unresolved. | No. |
| [Analytic Gaussian readout](../model/scalar_observations/experimental_readout/README.md) | The same 100 reserved complex probes give RMSE `0.07639048`. No interval violations under the declared sampling model. | Measured-axis factorization predicts interior probes. Its broad sampling intervals assume 250 independent bounded observations per mean; accepted counts are not authenticated. `not_rejected_under_assumed_sampling` is not an empirical validation. | No; same records, another analysis. |
| [Driven anharmonic dynamics](../model/anharmonic_observations/README.md) | 50 complex zero-duration probes calibrate preparation; 707 positive-duration scalar means give primary RMSE `0.11116365`. Harmonic control `0.28563357`, no-dephasing control `0.12061688`, frozen-state control `0.34737797`. | Supplied pulse Hamiltonian, preparation model, time convention and dephasing channel; no fit to positive-duration means. These controls compare those effective models. | No native C/E likelihood is evaluated. |
| [Changed-potential transfer](../model/anharmonic_transfer/README.md) | Five preparations and 125 positive-duration derived position estimates; RMSE `0.28764279` in dimensionless estimator units. The [primitive runtime](../model/primitive_observer/README.md) reproduces this score on the same records. | Declared changed potential, preparation and detector map. Material trajectory discrepancies remain. The runtime replay is not another dataset. | No. |
| [Simpler-action comparator](../model/identifiability_checks/README.md) | The same 125 estimates give cosine-Hamiltonian RMSE `0.31768062`, versus pulse-model `0.28764279`. | Both Hamiltonians are supplied. Their difference does not identify a native repair operation. | No. |
| [Readout identities](../model/readout_controls/README.md) and [quantum characteristic constraints](../model/characteristic_constraints/README.md) | Reused oscillator means have action-independent residuals. Under nominal pairing the probability-RMSE floor is `0.03419381` on the 707-point target. Gram witnesses require correction to at least one complex mean by at least `0.02439` or `0.03129` for coordinate extents 4 or 5. | Deterministic consistency diagnostics on noisy exported means. No accepted trial counts or authenticated detector-error model establish a statistical rejection. Changing the oscillator action cannot remove an identity-level discrepancy. | No. |
| [Inherited electromagnetic endpoint](../model/predictions/README.md#inherited-coupling-comparison) | Native observer attachment uses inverse alpha `136.994835177413`, versus CODATA `137.035999177(21)`: a `300.388218 ppm` deficit. Its numerical enclosure misses the reference five-standard-uncertainty interval. | Excludes treating this particular output as the complete endpoint with zero missing correction. Missing physical transport corrections have no certified error law; this is not an OPH-wide significance or rejection. | No; the endpoint is shared input, not a C/E discriminator. |
| [Bell-record replay](../model/predictions/BELL_RECORDS.md) | All 4,746 candidate rows are retained; the published selection gives 245 eligible trials, 196 CHSH-game wins, `S=2.42249987046...` and calibrated fixed-horizon null tail `p=0.03917464356...`. | Reproduces the published local-response null calculation with its setting allowance, locality/access and stopping premises. It rejects the local-factorized classical extension at the 5% threshold under that contract. | No. The full quantum C/E instruments are outside that classical null. |

The Bell number is a tail bound under a declared null, not a probability that
OPH is true or false. The code does not authenticate locality, settings,
event-ready timing or data custody. The CHSH-game win fraction `196/245` is
not the “plus” frequency in the native interference protocol; substituting
one for the other would compare different observables.

[Host-clock captures](../model/native_clock/HOST_CLOCK.md) are a separate
measured software benchmark. The [joint clock law](../model/native_clock/CLOCK_BRIDGE.md)
accounts for finite incident/current-count conversion uncertainty and retains
baseline acceptance and identical-record workload rejection. These timestamps
do not identify native ticks with physical proper time or modular phase.

[Preparation probes](../model/empirical_diagnostics/README.md),
[calibrated measurement plans](../model/calibrated_measurement/README.md) and
random-setting Gram tests specify future observations. Their synthetic
controls and unfilled count templates are not experimental measurements.

## What the repair-rule comparison establishes

The [candidate guide](../model/CANDIDATES.md) owns the complete model definitions.
The labels C and E refer to how the native quantum instrument records an
event; they are not labels for the alternative fine-structure closure maps.

| Candidate | Executed distinction | Empirical status |
| --- | --- | --- |
| N1 / C: current record | Same load repair as E; retains coherences that its current label does not distinguish. | No matched physical dataset selects it. |
| N1-E / E: contrast record | Same classical transition law; records extra contrast information and gives different coherent readouts. | No matched physical dataset selects it. |
| N2: local six-level response | Preserves the native load/record marginal; a declared internal-response assay distinguishes it from N3. | Load/current data alone cannot select it; physical internal-response attachment missing. |
| N3: flat endpoint response | Preserves the same load/record marginal and has trivial closed-loop response under its global coordinate change. | Same empirical boundary as N2. |
| N4: reciprocal load–response exchange | Private preparations affect bare-occupation and memory assays. This differs from the N1 current-detector protocol. | No identified physical preparation/readout dataset tests those distinctions. |
| S1 and supplied field benchmarks | Executed free/interacting field and detector laws with declared inertia, action, state and time. Oscillator comparisons evaluate particular supplied effective models. | Conditional benchmarks; not a selection of C, E or a unique microscopic repair law. |

For every one of the 78 load basis states, the C and E instruments have the
same joint law of coin, current and post-repair load when E's extra contrast
label is summed out. The equality continues by induction through finite
adaptive classical protocols from diagonal preparations with the same
populations, using common coin/current access and outcome-by-outcome
diagonal-preserving controls, with the same scheduler and stopping rule. Under
that shared observation contract, any admissible classical record dataset
has likelihood ratio `L_C/L_E=1` wherever defined. Feedback using E's extra
contrast label, or coherent interventions, is a different contract.

For the public oscillator and Bell records there is an additional problem:
no identified native C/E apparatus map supplies the two laboratory
likelihoods. “Not specified” is not a computed Bayes factor of one. Better
oscillator fits, classical agreement, and a Bell-null rejection cannot fill
that missing map.

The exact [source discrimination](../model/observer/source_discrimination/README.md)
also identifies E as C followed by retained seam-parity measurement. Forgetting
the extra record does not undo its dephasing. A measured E-like response can
therefore distinguish the C/E instrument hypotheses only with the relevant
ancillary/environment records and errors fixed; unrestricted unmodelled
parity recording can make C look like E.

## A discrimination experiment that is ready mathematically

The [fixed interference protocol](../model/predictions/README.md) uses one
native repair on seam `(0,1)` after a specified coherent preparation of two
load configurations. Its shared binary projector gives:

| Preparation | C: plus probability | E: plus probability |
| --- | ---: | ---: |
| Phase 0 | `5/8` | `3/8` |
| Phase pi | `1/8` | `3/8` |
| Phase pi/2 | `3/8` | `3/8` |
| Equal incoherent mixture | `3/8` | `3/8` |

This is an executable, falsifiable instrument distinction. The fixed design
requires 8,192 fresh independent attempted trials at each setting, every joint record bin and
abort, and an independently established apparatus probability-error budget
of `1/64`. Its simultaneous sampling-failure bound is below `0.000012606`.
Both or neither candidate may survive; survival is not proof of correctness.
The actual retained empirical evaluator returns **`no_data`**, with no
candidate verdict. The C and E demonstration inputs are explicitly synthetic.

The [continuing-carrier protocol](../model/closed_prediction_protocol/README.md)
removes exact stationary resets and actual trial independence. It uses 379
unread preparation attempts before each of 3,927 trials and three unread
attempts in a useful probe. Useful-branch matching probabilities are `341/432`
for C and `1343/2000` for E; herald failure `28/33` remains in the denominator.
The exact-source wrong-choice bound is `3/400`. A 1% physical bound additionally
requires independently authenticated conditional whole-trial service error
at most `1/1570800`; this is an unachieved physical tolerance, not a measured
gate fidelity. The retained two-trial executions are not an experimental
3,927-trial campaign.

An empirical distinction therefore requires identifying and calibrating that
coherent preparation/repair/readout operation, then collecting complete records
under one fixed interpretation for both candidates. It needs no fundamental
spacetime reconstruction or universal clock. Implementing a chosen C or E
channel on a quantum computer would test that engineered instrument; selecting
an unknown natural process additionally requires identifying it independently
of the candidate being tested. **There is no justified empirical winner among
the native rules on the retained data.**

## Other empirical comparisons in OPH

These families are part of the wider research programme. They are not extra
replications of the native observer experiments. The linked maintained
ledgers retain the individual rows, inputs, less successful outcomes and
reference editions; their numerical agreement is not automatically a
likelihood or an independent prediction.

| Family and authoritative evidence | Established comparison and boundary |
| --- | --- |
| [Fine structure and tau](research/INDEPENDENT_POSTDICTION_COMPARISON.md) | Structured inverse-alpha output `136.994835177413`; the gauge-width variant gives `137.035660136947`, a `2.474095 ppm` deficit. Closer agreement does not select the physical continuation. Conditional tau output `[1776.968991,1776.969063] MeV` is about `0.434` reported reference standard uncertainties from PDG 2024, using measured electron/muon masses and balance/order premises. It is exactly the same forward map as the balanced Koide reference, so that observable cannot distinguish the two formulations. |
| [Charged-lepton scale and HVP accounting](POSTDICTION_LEDGER.md) | Target-anchored closure intervals and the conditional response architecture retain their measured ancestry. One KNT19 hadronic accounting replay is compatible; the independent multi-class comparison is not evaluable. Neither the required correction nor its uncertainty is derived by agreement alone. |
| [Electroweak, QCD and nucleon comparisons](POSTDICTION_LEDGER.md) | Conditional Higgs/top comparisons, threshold-dependent QCD scale and an external-lattice-ratio nucleon estimate retain supplied transport/theory inputs. W/Z running chart values are not evaluable as physical pole masses without the missing same-quantity map. No source-only absolute quark masses are emitted. |
| [Flavor candidates](../code/particles/RESULTS_STATUS.md) | The target-informed weighted-cycle neutrino candidate is rejected by the declared correlated NuFIT 6.1 gate (`Delta chi2=20.11955` or `18.43528`, versus threshold `11.82916`). Fitted absolute splittings are not predictions. The [postdiction ledger](POSTDICTION_LEDGER.md) also retains the rejected declared down-quark Clebsch route. These do not identify C/E. |
| [Cosmology postdictions](../code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md) | Primordial tilt/running, tensor, isocurvature, non-Gaussianity, curvature, capacity and dark-energy comparisons retain their premise and reference ancestry. Edge tilt near `0.966021` is compatible with several subsets but has a `-2.43` marginal standardized residual for the retained ACT+Planck+lensing+DESI DR1 combination. Seen-data diagnostics are not frozen scores. |
| [Full CMB readout](../evidence/observer_dynamics_20260925/sim-analysis/codex/observer_cmb/measured_comparison.json) | Planck TT/TE/EE marginal standardized-residual RMS is approximately `1.01363/1.01976/1.01944` on `2471/1967/1967` multipoles. Amplitude, tilt, background and standard transfer are imported; covariance and nuisance treatment are absent. No joint likelihood or independent native prediction follows. |
| [WMAP angular-template exposure](../evidence/custody/falsification/frozen_targets/fz11_2026-07-31/frozen_target_spin_six_primitive_port_prediction_2026-07-31.md#prior-related-exposure) | The retained custody document reports a preregistered A5 search with null family-wide `p=0.64`. This is a documented exposure, not an independently replayed raw-data campaign in the current comparison packages. Its data and CMB likelihood class are excluded from an FZ-11 verdict; the frozen dispersion-coefficient relations have no qualifying physical comparison. |
| [Fixed-capacity dark-energy comparison](../code/cosmology/fixed_capacity_wlaw/README.md) | Retained DESI DR2 Gaussian diagnostics for `(-1,0)` range from `2.60` to `3.93` standardized units across overlapping combinations. These are not frequentist exclusions or C/E tests; chain-based interpretation and eligibility are governed by the declared protocol. |
| [SPARC and dark-sector diagnostics](../code/cosmology/rar_deep_regime/README.md) | Deep-regime acceleration/BTFR comparisons fit an acceleration constant and expose tension under declared weighting and mass-to-light choices. Full-RAR calibration reuses the same data and a published interpolation; `a0` is fitted. Lensing/quadrupole calculations with uncertified integration or supplied field laws are not native physical confirmation. |
| [IBM quantum hardware](../code/ibm_quantum_cloud/README.md) | State/readout/layout checks have mixed outcomes and incomplete raw custody for some legacy summaries. The blinded record-gated controller passes its frozen primary nulls on development and held-out QPUs. This supports the programmed implementation. Standard quantum mechanics predicts the circuits; neither OPH-versus-QM discrimination nor native C/E selection follows. |
| [Frozen physical targets](FROZEN_PREDICTION_LADDER.md) | Cosmology, angular invariants, conditional tau, dispersion, record-count, neutrino and Kerr targets have their own frozen, pending, deferred or void status. A registered prediction, custody proof or synthetic sensitivity test is not an acquired empirical verdict. Their immutable contracts remain authoritative. |

The [falsification programme](OPH_FALSIFICATION_PROGRAM.md) and
[observation ledger](registers/OBSERVATION_LEDGER_V3.md) distinguish these
bounded exclusions and conditional comparisons from programme-wide claims.
This report is a navigation and interpretation summary, not a replacement
prediction register or a statistical combination of dependent results.

## Reproduction and decision status

From the repository root:

```sh
python -m model.predictions --check model/predictions/results.json
python -m model.predictions.bell_records --check model/predictions/bell_results.json
python -m model.closed_prediction_protocol --check model/closed_prediction_protocol/results.json
python -m model.replay_empirical
```

The [reproduction guide](../REPRODUCE.md) gives the pinned environment and
broader theorem, source and empirical replays. Exact source-discrimination,
complete-count validation and synthetic C/E controls are covered by the
existing tests. Passing those checks establishes reproducibility of the
stated analysis; it does not turn synthetic records into acquired data or
authenticate a physical apparatus.
