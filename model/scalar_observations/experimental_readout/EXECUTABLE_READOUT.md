# Physical readouts from the executable observer

This comparison calibrates two preparation parameters from 21 complex axis
measurements, then uses native source operations to predict 100 other complex
measurements. The predicted binary probabilities are the norms of the executed
pointer branches. The analytic Gaussian formula estimates the preparation
parameters and independently checks the result; it does not supply the
predicted readouts.

The observations are the zero-simulated-duration trapped-ion oscillator data
of [McGarry et al.](https://arxiv.org/abs/2603.04744v1), deposited under
[DOI 10.5281/zenodo.18752287](https://zenodo.org/records/18752287).
This is an executable, calibrated prediction for physical measurements
excluded from the parameter fit. The public data had already been examined
in the [axis-factorization comparison](README.md), so it is a
retrospective test rather than a prospective or blinded prediction.

## Retained result

| Coordinate convention | Fitted alpha | Fitted nbar | Spin-mean RMSE | Binary-probability RMSE | Ideal occupation probability bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| Paper extent 4, primary | 0.96341309 | 0.12777123 | 0.06647166 | 0.03323583 | < 3.82e-7 |
| Deposit extent 5, alternate | 0.76974579 | 0, boundary optimum | 0.07178838 | 0.03589419 | < 1.59e-7 |
| Zero-response diagnostic | — | — | 0.12810121 | 0.06405060 | — |

Each score uses all 200 separately reported real/imaginary spin means in the
reserved quadrant. The probability RMSE is the spin-mean RMSE divided by
two. It is an average discrepancy, not a guarantee for each point. The
occupation bound is uniform over these settings at the fitted parameter
point; it does not cover uncertainty in that fit or laboratory errors.

For a concrete unselected example, the first reserved grid point under the
primary convention is `beta=0.4+0.4i`. The source predicts outcome-zero
probabilities **0.79341921** for the real readout and **0.78493547** for the
imaginary readout. The reported experimental frequencies are **0.716** and
**0.828**, respectively. All other predictions and discrepancies are retained
with equal status in the receipt.

In the retained run, the maximum raw full-phase source-operation residual
is below `9.98e-11`; the corrected residual is below `5.42e-13`. Maximum
probability-mass residual is below `3.2e-15`. These are numerical diagnostics.
The independent unrestricted Gaussian calculation differs from the emitted
probabilities by less than `8.2e-14`; it is a reference check after execution.
The explicit implementation uses at most **3,489,068 native pulses per shot**,
including preparation and both coordinate-basis passes, with physical
initialization, storage and final readout costs outside that count.

## What executes

The declared single-mode state is `D(alpha) rho_thermal(nbar) D(alpha)†`, with
real displacement `alpha` and nonnegative thermal occupation `nbar`. The
[calibration](executable_calibration.py) fits those two parameters on the
axes only, with no fitted visibility, offset or origin renormalization.
The noisy origin and all other axis values remain in the fit. A fitted
thermal parameter is not an independent laboratory temperature measurement.
The 33 fixed phase starts provide numerical basin seeds. An 80-decimal-digit
stationary refinement reevaluates the objective from the original axis values
and selects the final binary64 parameter point, including the thermal boundary.
This makes the declared point reproducible without certifying a global optimum
or a parameter-confidence region; optimizer start identities are not the
scientific fit output.

The [source implementation](executable_source.py) uses 64 occupation levels:

1. A supplied classical draw selects a capped geometric thermal occupation.
   Native words prepare that Fock code from the blank register. No draw is
   discarded or conditioned on acceptance.
2. A source-compiled real displacement prepares the finite analogue of the
   calibrated state. A compiled number rotation selects the measurement
   quadrature.
3. The full coordinate-basis words and conditional phases implement the
   smooth two-outcome pointer detector. Phases zero and minus pi/2 read the
   real and imaginary characteristic components.
4. Norms of the resulting pointer branches give both binary probabilities;
   averaging over the supplied classical preparation weights gives the
   predicted readout. The source function accepts settings and preparation
   parameters, but no observed responses.

Every basis, displacement, rotation and detector phase is contracted from
the existing twelve-port source pulses. The Gray compilation is extended
locally to six field qubits using the same gates; a seventh qubit is the
pointer. Dense basis matrices are cached only after executing their words.
The requested words are generated in a private 80-decimal-digit arithmetic
context from the exact ratios of the declared binary64 parameters. The
coordinate eigensystem, signs and full Gray/Walsh schedule are evaluated there;
only angles passed to the existing compiler become binary64. Independent
SciPy targets continue to check the executed operations. Fixed arithmetic
defines replayable words, not an interval accuracy or physical-service claim.
The receipt counts the full word at every measurement and binds its hash.
Independent finite-matrix targets check execution; they do not replace the
executed operations.
The million-pulse operations use fixed blocks of at most 4,096 native pulses
for host contraction. The shared compiler's polar thresholds are unchanged:
they apply to every executed block and to the product of the corrected blocks.
The `1e-10` raw-defect and `5.1e-11` polar-correction guards apply to
individual blocks and the corrected-block product; they do not bound the
complete raw-word defect or accumulated correction.
Both the complete raw product and the final corrected product must pass the
`2e-10` phase-sensitive target gate. The receipt reports accumulated numerical
corrections and raw drift. Blocking preserves the requested source word and
its pulse/service counts; these diagnostics are not outward roundoff or
physical-error bounds.
The [derivation](EXECUTABLE_DERIVATION.md#long-source-word-contraction) gives
the exact matrix-product telescoping identity and its floating-evaluation scope.

This is a computation through the source model, not execution on laboratory
hardware. The blank register, classical thermal draw, coherent code services
and final pointer readout remain supplied apparatus. Their real physical
costs and errors are not established by counting native pulses.

## Approximation and statistical scope

The [occupation derivation](EXECUTABLE_DERIVATION.md) bounds the actual capped
thermal preparation, finite displacement and full pointer dilation against
the unrestricted oscillator. Both outcomes and the disturbed terminal field
are covered. Exact rational moment bounds use the exact binary values of
the fitted parameters and measurement coordinates. Numerical contraction
residuals and floating thermal-weight errors remain separate from that ideal
mathematical certificate.

Two coordinate extents are retained: four from the paper and five from the
dataset description. Neither is selected by its agreement on reserved
measurements. Reflection structure fixes the array axes; a fitted
displacement sign is in array coordinates and does not resolve the physical
direction. The two calibrated models also have different occupation errors.

The RMSE is descriptive. Calibration-parameter uncertainty is not propagated,
and the deposited aggregate means do not authenticate accepted shot counts
or apparatus errors. The axis-product confidence calculation does
not apply to this different two-parameter fit. The laboratory observations
remain conditional on its preparation selection, even though this executable
finite preparation has no discarded branch.

This establishes a preparation/readout prediction within the supplied
oscillator model. The oscillator family and apparatus mapping are premises;
OPH has not selected them independently. Agreement therefore does not
distinguish OPH from ordinary quantum mechanics. There is no intervening
interacting evolution here, no test of the seven-mode scalar dynamics and
no calibration of a physical clock.

## Reproduce

The [protocol](executable_protocol.json) fixes the family, split, both scales,
dimension and scoring rule. The [receipt](executable_results.json) retains
every source operation ledger, both predicted outcomes, every observation,
residual and rational occupation certificate. `all_pass` concerns software
and numerical contracts, not a scientific significance verdict.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.scalar_observations.experimental_readout.executable --check model/scalar_observations/experimental_readout/executable_results.json
python3 -m pytest -q -W error model/scalar_observations/experimental_readout
python3 -m model.verify --hashes-only
```

Replay is offline and can take tens of minutes because the full dense source
words are executed. The original dataset, Gaussian-product protocol and receipt,
and seven-mode receipts remain unchanged. The executable protocol and receipt
were revised for this numerical replay policy; the
[original executable receipt at `df89a7e5`](https://github.com/FloatingPragma/observer-patch-holography/blob/df89a7e5/model/scalar_observations/experimental_readout/executable_results.json)
remains available in repository history. Dataset authorship, CC BY 4.0
license and original-file checksums are recorded in [provenance](provenance.json).
