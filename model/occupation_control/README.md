# Controlled occupation removal for the oscillator apparatus

This packet implements the occupation-error attachment recommended by the
[mathematical audit](../audit_20261015/MATH.md#3-a-directly-applicable-occupation-truncation-theorem).
It bounds a continuous, unrestricted oscillator's records using an
interval-enclosed finite computation. The bound includes thermal preparation,
preparation displacement, the complete pulse history, the pointer detector,
time discretization and evaluator roundoff.

The result is useful early in the evolution and increasingly loose later.
For case `a` of the supplied [anharmonic apparatus](../anharmonic_transfer/README.md),
the retained 64-level run has complete half-trace error below `1.32e-8` at
0.4 ms and `0.002041` at 2 ms. At 10 ms it is about `0.19946`, too loose to
justify a precision claim. The 32-level bound exceeds one at 10 ms and does
not certify that the selected reference branch even has positive mass. Small
differences between two cutoffs had not established these guarantees.

| Published coordinate | 32-level half-trace radius | 64-level half-trace radius |
|---|---:|---:|
| 0 ms | 6.75e-13 | 9.51e-13 |
| 0.4 ms | 0.009815 | 1.319e-8 |
| 2 ms | 0.247765 | 0.002041 |
| 5 ms | 0.688188 | 0.026703 |
| 10 ms | 1.674805 | 0.199459 |

These are errors relative to the declared apparatus law, not experimental
error bars. They are conservative integrated residual bounds. A large bound
does not demonstrate a large actual truncation error. The long-window bound
is dominated by occupation removal; making the time integrator more accurate
alone will not fix it. The retained run covers case `a`; the other four
potential phases are supported by the evaluator but have no retained full
window certificate in this packet.

For the 64-level run at 2 ms, the selected two-point estimator is
`-0.1434118369` with certified numerical/occupation radius `0.016264164`.
At 0.4 ms the corresponding radius is below `9.905e-8`. These radii include
terminal conditioning and final scalar evaluation; they are wider than the
complete-state half-trace radii because selection and estimator scaling
amplify error.

## What runs and what is measured

The evaluator receives the previous, unchanged preparation parameters
`alpha=-1.0081292326480422`, `nbar=0.13314125172459518`, and the supplied pulse,
noise and detector settings. It never reads empirical outcome files. Its
computation is:

1. Prepare a capped thermal mixture, charge its exact tail and rounding, and
   evolve the real-displacement preparation.
2. Retain the spin while evolving the negative/positive/negative sideband
   pulses with number dephasing and the two intervening spin rotations.
   Finite carrier gaps are included for nonzero phase cases.
3. At each requested endpoint, attach a pointer in `|+>`, execute opposite
   half-displacements and read `-Y_pointer`, keeping both final spin branches.
4. Report the spin-selected detector mean and its two-point position
   estimator `Im chi(i*0.4)/(sqrt(2)*0.4)`. This is not exactly `<x>`.
5. Bound every polynomial segment's full differential-equation defect and
   endpoint discontinuity with outward arithmetic, including the single
   occupation edge crossing out of the retained space.

The retained paths use order-12 polynomials, five cells per 50 microsecond
sideband, eight cells for each displacement, and two for each spin rotation.
The integration parameter is supplied drive elapsed duration; displacement
and instantaneous rotation segments use their own normalized control
parameters. Host runtime has no physical-clock interpretation.

**This is an independent semantic evaluator and error certificate.** It does
not execute native observer source words. The apparatus, preparation family,
quantum structure, noise and clock remain supplied. Source-derived action
selection is unchanged. This packet does not replace native execution
receipts. For orientation, its case-`a` final estimator differs from the
corresponding historical finite prediction by about `2.4e-7` at both cutoffs; that numerical
agreement is separate from the error certificate.

## Guarantee and assumptions

The [derivation](DERIVATION.md) maps the existing
[Etienney–Robin–Rouchon theorem](https://arxiv.org/html/2501.09607v6) to this
apparatus and constructs its conservative unrestricted propagator. The
finite polynomial need not be positive. Its full continuous defect, rather
than a collection of sampled boundary populations, controls the comparison.

The interval kernel assumes ordinary IEEE-754 binary64 round-to-nearest
arithmetic, gradual underflow and correctly rounded elementary operations
and square root. It avoids BLAS in the enclosure calculation, widens each
operation, and fails on nonfinite bounds. Target numeric parameters are the
exact binary values in the receipt, not exact symbolic multiples of pi.

The complete error is `half_trace_error_upper`. Final selected records have
separate probability, mean and estimator radii, including scalar evaluation
rounding. Conditional certificates are unavailable when the bound cannot
guarantee positive reference acceptance. The packet additionally reports
`unrestricted_vs_exact_finite_half_trace_error_upper` for composition with
independent finite-source execution bounds.

The certificate excludes experimental parameter uncertainty, detector
normalization provenance, preparation-transfer validity, the laboratory's
unmodeled selection dwell, native compilation error and physical service
error. It establishes no different empirical score or OPH-specific law.

## Next steps and reproduction

Use the early-window bound when testing proposed source laws; there is now
an explicit numerical uncertainty to keep below the discriminator size.
For longer windows, enlarge the retained occupation space adaptively or use
a tighter boundary-row nuclear-norm enclosure. The current boundary is
computed only at the last occupied level, so adaptive refinement can reuse
the same theorem and sparse solver. Native compilation can be reserved for
the candidate protocols that survive these inexpensive semantic tests.
Larger semantic cutoffs do not by themselves extend the current native
compiler's 64-occupation limit.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/occupation_control
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.occupation_control --check model/occupation_control/results.json
```

To generate a fresh result, replace `--check` with `--output <path>`.
The two retained windows take a few minutes on the development host. Replay
recomputes every bound; it does not infer a certificate from agreement with
the stored rounded numbers. Each receipt's certificate belongs to its own
computed binary path. Cross-platform replay tolerance is a reproducibility
check, not permission to attach a stored radius to different computed means.
Tests include exact rational arithmetic checks,
independent high-precision generator enclosures, actual occupation escape,
a deliberately inaccurate integration step, detector phase conventions and
failure to certify terminal selection.
