# An executable modular clock for the native record-transfer description

The native current records determine a finite transfer description whose
modular flow can be executed through the existing coherent source compiler.
No oscillator action, mass or external clock calibration enters the phase
generator. The coherent code, probe preparation and selected modular intervals
remain apparatus choices. This is a finite operational connection, without a
claim that native dissipative waiting is unitary evolution or that a continuum
spacetime has been recovered.

## From actual records to the reference and generator

The [native clock](../native_clock/README.md) proves marked strong lumping of
all 66 native configurations into three token-separation classes. Its exact
transition matrix and stationary class weights are

```text
P = [[14/15, 1/15, 0], [1/15, 9/10, 1/30], [0, 1/6, 5/6]],
pi = (5/11, 5/11, 1/11).
```

Reversibility makes `S=D^(1/2) P D^(-1/2)` symmetric, where `D=diag(pi)`.
The existing full transfer has exact bounds `(2/3)I <= K <= I`; its invariant
class-constant restriction is positive as well. The source characteristic
polynomial factors as

```text
(x-1)(x^2-5x/3+311/450).
```

Thus the three eigenvalues are `1`, `lambda_+=(25+sqrt(3))/30` and
`lambda_-=(25-sqrt(3))/30`. On this reconstructed three-dimensional Hilbert
space, `sigma=S/Tr(S)=3S/8` is faithful and

```text
H = -log S,                 log sigma = -H - log(8/3) I,
sigma^(is) A sigma^(-is) = exp(-isH) A exp(isH).
```

The dual state evolution is `exp(isH) rho exp(-isH)`. This fixes the orientation
convention; the source phase is `+s E_j`, with `E_j=-log(lambda_j)`.
The reference is derived from the positive native transfer. This does not
prove that OPH uniquely selects this reference among all admissible ones.

The reconstruction uses the finite reflection-positive transfer argument in
[the native transfer derivation](../primitive_observer/DERIVATION_TRANSFER.md),
the exact marked reduction, and finite inner modular flow as in
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean). The additional
work is an explicit compatible source program and its numerical certificate;
no new abstract modular theorem or Lean proof is claimed.

This `M3` is a representation of the reduced history, not the native `M66`
working-state algebra. Its reference differs from the diagonal Palm reference
used by [the native modular probe](../native_modular_reference/README.md).
An equality of their physical interpretations is not assumed.

## Complete source execution

The spectral basis is encoded in addresses 0, 1 and 2 of two work qubits.
Address 3 is explicit padding. The complete source target is
`diag(1, exp(is E_+), exp(is E_-), 1)`. It implements the modular action on the
whole effective `M3` algebra, with a declared padding extension. On the block
algebra `M3 direct_sum C`, the positive extension `S direct_sum 1`, normalized
by `11/3`, produces the same action. This is not a new spatial dimension.

Each probe starts from the declared blank code and executes:

1. An emitted balanced source word preparing `(ket(0)+ket(j))/sqrt(2)`.
2. The complete Walsh phase word, whose phases come from the source spectrum.
3. An emitted X or Y analysis word and a complete four-outcome computational
   read, followed by the emitted inverse analysis word.
4. For the continuing experiment, another modular segment and another X read
   on that same disturbed state, retaining both records and all padding.

The ideal single-read probabilities at address zero are

```text
p_X(0;s,j) = [1+cos(s E_j)]/2,
p_Y(0;s,j) = [1+sin(s E_j)]/2.
```

At modular parameter `s=8`, the X probabilities are approximately `0.801845656`
and `0.277055154` for the two modes. There is no physical-seconds conversion.
Both quadratures give a phase reading modulo a period; a single random
outcome does not determine elapsed time. The listed phase settings define a
programmed probe, rather than an autonomous clock mechanism.

The continuing `(4,1)` experiment retains 16 two-read records and uses 430
primitive source pulses, or 860 proper-code service events. For
`q_k=[1+cos(s_k E_j)]/2`, its nonzero ideal record probabilities are
`q1*q2`, `q1*(1-q2)`, `(1-q1)*(1-q2)` and `(1-q1)*q2`. Repreparing the original
state between reads gives a different law and is a regression control.

Every operation is carried by the shared validated `runtime.Program`.
A scoped pure-state Kraus evaluator applies its original binary arrays to
vectors and retains every output, including tiny padding leakage. This avoids
subtractive cancellation at nominally zero branch weights. It neither repairs
nor relaxes the generic runtime's relative positivity guards. An independent
generic-density execution checks the complete unconditioned terminal state.
The scoped checkpoint serializes the exact binary vectors and resources;
resumption verifies them against a fresh complete prefix. It is not asserted
to be the generic runtime's checkpoint format.

## Connection to measured native correlations

The same two modes give the connected covariance of two actual native
movement bits separated by `n` unread attempts:

```text
C_E(n) = w_+ lambda_+^n + w_- lambda_-^n,
w_+ = 1/7260 + sqrt(3)/16335,
w_- = 1/7260 - sqrt(3)/16335.
```

The producer checks this against exact full 66-state marked path sums at
`n=0,1,4,8`. The new X/Y source readouts reconstruct
`C_M(s)=sum_j w_j exp(is E_j)`, the corresponding analytic continuation.
The marked transfer matrices are **not POVM effects**; their measurement
semantics remain the native current instrument. The new coherent probes
measure spectral phases in the reconstructed code. Analytic continuation
plus its executable simulation is a supported correspondence, while agreement
of that interpretation with independent physical observations remains open.

## Error and scope

Exact rational square-root and atanh-series intervals enclose the source
energies. A real-argument complex Taylor remainder encloses every target phase.
Each original binary source matrix is compared on all four input columns with
the complete ideal unitary. Its Frobenius bound supplies an operator bound
`epsilon <= 1e-11`, half-diamond bound `epsilon(2+epsilon)/2` and channel-norm
bound `(1+epsilon)^2`. The certificate retains amplification through every
source use and arbitrary untouched reference systems. The full largest
seven-gate experiment has half-diamond allowance at most `1e-10` against the
declared ideal code experiment; the intervening exact computational read
contracts the comparison.

That certificate covers raw compiled gates. Subsequent binary state propagation
and scalar probability residuals are diagnostic, and physical apparatus
service errors are unbounded. Tiny raw normalization defects are retained as
record weights. No empirical counts, interacting continuum limit or Lorentz
metric are supplied here. A further physical identification must preserve
observer-accessible correlations and intervention responses, without requiring
one universal clock or fundamental fields.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.transfer_modular_probe --check model/transfer_modular_probe/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/transfer_modular_probe
```
