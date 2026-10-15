# Empirical and apparatus audit, 15 October 2026

The current numerical receipts support approximate reproduction of a supplied
quantum apparatus model. This review found no verified factor-of-two,
dephasing-rate, terminal-selection, or score-arithmetic error in the declared
anharmonic implementation. It did find a stronger missing comparison and
substantial unresolved identification problems between that implementation
and the deposited experiment. Those problems limit physical interpretation
more than the numerical integration errors do.

This audit is retrospective: all existing outcomes were available. It changes
neither frozen predictions nor their scores. New diagnostic calculations are
in [empirical_probe.py](empirical_probe.py), with complete outputs and input
hashes in [empirical_probe.json](empirical_probe.json).

The transfer input is the [compressed original receipt](data/transfer_receipt_93e59bfb.json.gz)
from baseline `93e59bfb`. Its decompressed SHA-256 remains
`71403f474c39705d20a6c30dc2f36b7556943dd4f1734a5a1a395184744d2398`.
The probe validates that digest before parsing; its `transfer_receipt` input
hash refers to these original bytes, and the manifest binds the compressed
artifact. Subsequent live source-custody or numerical replay changes do not
replace the historical preparation, scores or zero-time residuals used here.

## 1. The strongest missing simple comparison

Frozen and harmonic controls remove the effect the apparatus was designed to
produce. They establish that the anharmonic dynamics matters, but give little
information about whether the detailed spin/pulse/environment construction is
needed to reproduce these particular observations.

I independently integrated the considerably simpler, ancilla-free law

\[
 H=2\pi(500)N+4000\cos[(\pi/3)(a+a^\dagger)+\Phi],\qquad
 L=\sqrt{18}N.
\]

This is the target effective Hamiltonian already identified in the frozen
packet. It uses the same transferred preparations, potential phases, literal
`h=.4` detector and published simulated-duration coordinate. There are no
additional fits. The independent DOP853 integration uses directly constructed
matrices, no source compiler, and no state normalization or positivity repair.
This is a simpler scientific comparator, not an alternative native execution.

| Case | Executed physical-pulse model | Simpler target Hamiltonian |
|---|---:|---:|
| a | 0.357529 | 0.394860 |
| b | 0.283991 | 0.335332 |
| c | 0.269777 | 0.300591 |
| d | 0.302882 | 0.330582 |
| e | 0.201735 | 0.191322 |
| All 125 observations | **0.287643** | **0.317681** |

Values are RMSE in the deposited position-estimator units. The detailed model
improves aggregate RMSE by **9.46%** against this comparator, rather than the
approximately 69% against frozen/harmonic controls. It improves four cases;
the simpler model does better in case e. The zero-prediction control gives
0.645911. These retrospective comparisons change neither the original test
protocol nor its valid arithmetic.

The target-model occupation comparison, 32 versus 48, changes an estimator
by at most `1.62e-7`; trace residuals are below `1.7e-15`. These are numerical
diagnostics, not certified continuum error bounds. The direct target model
has no terminal spin rejection because it omits the spin altogether; that is
part of the simplification being tested.

**Consequence:** retain the pulse model as the more detailed apparatus model,
but use this effective-Hamiltonian comparator in future protocol designs.
The data distinguish these two supplied descriptions only modestly at the
present level of agreement. Neither comparison distinguishes OPH from quantum
mechanics.

## 2. Readout normalization is a real unresolved provenance problem

The paper specifies `Im chi(i*h)/(sqrt(2)*h)`, `h=.4`, and
“200 measurements per point” for Fig. 5. Its Appendix C specifies two 50 us mid-circuit selections.
Its Appendix D uses `L=sqrt(18)*N`. Its Fig. 5 mapping puts the centred
preparation in d and the larger asymmetry in e.
([Primary paper, Sections II.2–III.2 and Appendices C–D](https://arxiv.org/html/2603.04744v1))

The local [deposited README](../anharmonic_transfer/data/SOURCE_README.md)
reverses d/e. Following the figure/prose and the zero-time centred row is
reasonable, but remains a documented resolution of conflicting metadata.

The audit extends the previous zero-time observation: **all 130 deposited
numbers**, including the five initial values, satisfy

\[
 |63.6\,x-\operatorname{round}(63.6\,x)|<3.13\times10^{-8}.
\]

Thus their printed quantization is compatible with hundredth-sized signals
divided by `.636`. This does not authenticate a probe coordinate, a gain,
accepted counts, or even the exact extraction procedure. In particular,
`.636` is not `sqrt(2)*.4`. A post hoc rescaling selected from residuals would
not resolve this problem.

This calculation uses the unchanged 1,770-byte deposited CSV, SHA-256
`83cdb519902985b53381e5310894f7f424b2cbbca86b1ed51f2f4015dd7d4624`.
The existing [provenance receipt](../anharmonic_transfer/provenance.json)
also binds its published MD5; the audit output binds the exact local input
and probe source. No spreadsheet parsing or newline conversion intervenes
in the hash check.

The current comparison correctly preserves the published `h=.4`, scores the
derived observations directly, and refrains from inferring raw counts
([data parser](../anharmonic_transfer/data.py#L60),
[score](../anharmonic_transfer/experiment.py#L160)). Nevertheless,
the readout law's identification with these data is unresolved. It should be
resolved before claiming precision agreement, an empirical confidence level,
or a physical discrepancy with the supplied quantum theory.

## 3. Preparation transfer is not an authenticated common initial state

The primary transfers two parameters from a different zero-time scan. That
is an honest, fixed transfer test, but it is not evidence that every run had
the same displacement magnitude and thermal occupation. Its predicted minus
deposited zero-time estimator residuals are:

| a | b | c | d | e |
|---:|---:|---:|---:|---:|
| +0.340550 | +0.073254 | -0.261933 | +0.015723 | +0.104701 |

The largest mismatch occurs before evolution. The
[calibration](../anharmonic_observations/calibration.py#L15) correctly fits
only the earlier zero-time characteristic scan, while the
[secondary inverse](../anharmonic_transfer/calibration.py#L22) correctly uses
the finite-difference formula. Neither identifies non-Gaussian preparation,
detector gain, coordinate drift, or thermal occupation from one new number.
The inverse also has multiple branches before the imposed principal-branch
restriction. Consequently the secondary analysis worsening later RMSE does
not show that the old preparation was physically correct; it shows that this
particular one-number correction did not improve transfer.

**Consequence:** obtain independent preparation scans and their uncertainties
before adding more preparation parameters. Keep calibration and subsequent
evolution tests separate.

## 4. Complete finite model is not the complete timed laboratory circuit

The runner's continuous spin history and terminal-only postselection match
the intended selection semantics. It keeps the rejected mass and both
pointer outcomes ([source](../anharmonic_transfer/source.py#L137)). The finite
carrier construction also correctly separates published simulated duration
from drive elapsed duration. It assumes the detuning phase advances through
carrier gaps ([source](../anharmonic_transfer/source.py#L35)); the current
protocol already marks the laboratory compensation convention unauthenticated.

Additional selection dwell and detector-pulse durations are not evolved in
the current source. Preparation, selection and readout are ideal interfaces.
This is a declared scope limitation rather than a hidden numerical bug.
Once carrier durations are explicitly modeled, a complete timing audit should
also specify the detector's phase origin, phase tracking through selections,
finite readout displacement duration, and noise during those intervals.
Adding an arbitrary drift for these intervals without the phase convention
would be as unjustified as assuming it vanishes. No observed residual is
attributed to this missing information here.

The displayed paper equations contain inconsistent rotation-sign/axis
conventions, already acknowledged in the existing derivation. The source's
target-based convention is internally consistent, but cannot certify the
exact laboratory pulse program without independent pulse metadata. The ideal
finite-model splitting bound does not cover these apparatus choices.

## 5. Validate symmetry before using the archived two-dimensional scans

The paper's printed symmetry statements in Section II.2.3 and Appendix E
need correction before they can be used as reconstruction rules.
([Primary paper](https://arxiv.org/html/2603.04744v1))

They are printed as `chi(-beta)=chi(beta*)` and, for the initial state,
`chi(beta)=chi(beta*)` respectively.

Directly from `D(beta).dagger=D(-beta)` and Hermitian `rho`,

\[
 \chi(-\beta)=\overline{\chi(\beta)}.
\]

For the real-displacement thermal family, explicitly

\[
 \chi(u+iv)=e^{-(\bar n+1/2)(u^2+v^2)+2i\alpha v}.
\]

Thus even momentum gives `chi(u+iv)=chi(-u+iv)`, whereas conjugating beta
conjugates chi; it does not leave a nonzero real displacement unchanged.
These elementary counterexamples disprove literal use of the printed
alternatives. The present one-dimensional packets do not depend on that
reflection rule. Before using the archived 2D arrays, authenticate their axes,
which quadrant was measured, which entries were reflected, and whether any
background was subtracted. Do not count reflected entries as independent
observations.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python model/audit_20261015/empirical_probe.py
```

All existing scores are descriptive. The fixed preparation, apparatus,
environment and readout remain premises. Passing implementation tests and
small ideal finite-model error bounds do not authenticate those premises.
