# Conditional detector prediction

This protocol predicts an experiment with the declared scalar apparatus. It
is not a new observed constant or an empirical validation of OPH. No detector
output, measured coupling or fitted response sets the prediction.

## Fixed experiment

Prepare the native source uniformly on its 66 binary total-two states. Repair
seam `(0,1)` with a fair coin and retain its current `j`. In an independently
prepared scalar vacuum, apply momentum impulse `4*j*f` in the left slab.
Evolve the existing massive Dirichlet scalar action and measure the centered
binary Weyl instrument at the right slab. Retain its outcome `y=+1` or `-1`.
Only y selects the next native seam: `(0,5)` for +1 and `(0,1)` for -1.
Use a fresh fair coin, retain the second native current `j2`, and halt.

Use the existing golden grids q=13,34,89,233, their actual compact smearings
and mass one. The primary model time is `39/40`; `19/20` and `1` are additional
fixed settings. Model time is `c*t/L`, with no laboratory clock calibration.
No exposure, amplitude, mass, grid coefficient or time is fitted to outcomes.

If `r_q(t)=Pr(y=+1 | j=+1)-1/2`, the full protocol predicts

```
E[j2] = 5*r_q(t)/132
E[final port-0 load] = 1/6 - 5*r_q(t)/132
E[y*j] = 10*r_q(t)/33
Pr(y=+1) = 1/2
```

For q=233 the inherited outward certificate predicts a negative mean j2
throughout `[19/20,1]`, with magnitude at least
`(5/132)*(34438518047/250000000000)`, approximately 0.005218.
The continuum mean is within
`(5/132)*(12513633763/250000000000)`, approximately 0.001896,
of its finite-grid counterpart. These bounds are mathematical properties of
the ideal supplied apparatus; floating point diagnostics do not certify
implementation error.

## Controls and decision rule

Reversing the preparation sign must reverse E[j2]. Zero impulse or zero
detector coupling gives E[j2]=0. Independently shuffling donor detector bits
keeps both seam frequencies at one half but gives E[j2]=0. A controller
disconnected from the detector and fixed to `(0,1)` also gives E[j2]=0.
Both bit values and all native branches are retained; no postselection is
permitted. The premeasurement detector marginal alone cannot detect the
feedback effect: its probability remains one half.

The software run uses 200,000 fresh trials per grid at the primary time, seed
20261012. It computes its forecast before drawing outcomes and checks samples
against that forecast; it does not train on them. An independent dense-action
oracle checks the small grid, while coherent-state and truncated-oscillator
tests check the instrument and its backaction. Fixed endpoint times test the
same dynamics without fitting. These are software checks, not held-out
observations of nature.

For an independent apparatus test, fix trial count n, familywise significance
level delta, settings and service-error budget before collecting data. For
each of m bounded current means in [-1,1], use the Hoeffding half-width
`sqrt(2*log(2*m/delta)/n)`. A comparison with a continuum interval must include
the stated finite-grid bound and independently bounded preparation/control/
readout error. Reject this apparatus model if the resulting intervals are
disjoint. Numerical diagnostics require a separate numerical-error allowance.
Do not diagnose OPH from a mismatch before checking those apparatus premises.

## Continued-field protocol

The [two-read experiment](CONTINUATION.md) fixes times `19/20` and `1`, gamma=4
and the same four grids. It keeps the first conditional field and native
loads, adds a second field observation and uses its bit for a third native
repair. Both outcome signs and every source path are retained. Its frozen
software execution uses 50,000 trials per grid and seed 20261012.

At q=233 the ideal mathematical model predicts a zero-current detector-bit
correlation above `0.045619388905`. This distinguishes a fresh-field reset
but also occurs in a correlated commuting Gaussian control. The separate
measurement-disturbance prediction is stronger: on the positive-current
preparation, performing the first read raises the second plus probability
by more than `0.000022281162` compared with omitting it. The commuting control
has zero such disturbance. These are two distinct experimental questions.

The q=233 complete-history TV bound relative to the declared continuum is
`0.025709913643`. This controls probabilities of all complete recorded
events; an expectation in [-1,1] has the generic error bound twice that
number. It is not a bound on conditional probabilities after arbitrarily
rare postselection. The separate positive-current disturbance certificate
uses the declared preparation directly.

The software samples do not resolve the small disturbance difference at a
statistical confidence level. A prospective apparatus comparison must fix
independent performed/omitted-first-read batches, their preparation contract,
sample budget, rejection thresholds and service errors before collecting
data. Apply the stated concentration rule to the difference with both batch
errors included. A model-versus-reset or model-versus-commuting discrepancy
rejects that specified comparison, not all alternative theories or OPH.
