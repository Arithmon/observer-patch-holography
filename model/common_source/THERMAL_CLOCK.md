# A fixed thermal reference for the changing common history

Stationarity of a modular reference state does not prevent another state from
changing under its modular flow. The earlier clock comparison tested a more
specific identification: it used the actual post-measurement scalar marginal
as its own modular reference and compared that local modular flow with the
interacting Hamiltonian. That identification fails on the tested state. The
programmed local modular clock remains valid, and a different, fixed reference
can align modular evolution with the same finite Hamiltonian.

Here the reference is a faithful Gibbs state on the **full common quantum
work algebra**. The actual prepared and measured work state remains the one
already executed by the source. The reference is a declared mathematical
choice; it is never substituted for, mixed into or prepared as the actual
work state.

## Keep the reference and the actual state distinct

Let `H` be the existing finite interacting Hamiltonian on the work space of
dimension `4d^2`, with `d=2` or `4`. Choose a fixed finite `beta > 0` and put

\[
 \sigma_\beta=\frac{e^{-\beta H}}{Z_\beta},\qquad
 Z_\beta=\operatorname{Tr}(e^{-\beta H}).
\]

In finite dimension this is strictly positive and normalized. Its logarithm is
`log(sigma_beta)=-beta H-log(Z_beta)I`. In the convention already used by the
observer package,

\[
 \sigma_s^{\sigma_\beta}(A)
 =\sigma_\beta^{is}A\sigma_\beta^{-is}
 =e^{-i\beta sH}Ae^{i\beta sH}
 =\alpha_{-\beta s}(A),
 \quad \alpha_t(A)=e^{itH}Ae^{-itH}.
\]

The correspondence is therefore **`t=-beta s`**, including the orientation
sign. The scalar partition-function phase cancels. This is the familiar
finite Gibbs/modular identity, used here as an explicit matrix construction.
It is consistent with the distinction between a statistical state and the
flow it defines in the
[Connes--Rovelli thermal-time proposal](https://arxiv.org/abs/gr-qc/9406019).
No new general modular theorem is claimed.

The actual density `omega` need not equal `sigma_beta`. Under model-coordinate
time its state is

\[
 \omega(t)=e^{-itH}\omega e^{itH},\qquad
 \dot\omega(t)=-i[H,\omega(t)].
\]

For any detector effect `E`,

\[
 \operatorname{Tr}(\omega(t)E)
 =\operatorname{Tr}(\omega\,\alpha_t(E))
 =\operatorname{Tr}(\omega\,\sigma_{-t/\beta}^{\sigma_\beta}(E)).
\]

These expectations can change when `[H,omega]` is nonzero. In contrast,
`[H,sigma_beta]=0`, and expectations in the reference itself are stationary.
There is no conflict: the reference defines the flow; it need not be the
changing preparation whose observations read that flow.

The current receipt uses `beta=1` as a declared dimensionless convention.
Choosing a different positive `beta` rescales the modular parameter required
for the same model-time evolution. Neither this choice nor the identity
determines a laboratory temperature, second or energy unit.

For a fixed full matrix algebra and a specified `H`, the Gibbs form is also
necessary for this exact identification of automorphism groups. Equality of
their generators requires `[log(sigma)+beta H,A]=0` for every matrix `A`.
The commutant of the full matrix algebra consists of scalars, so normalization
forces `sigma=exp(-beta H)/Z`. This elementary observation does not select
`sigma`; it identifies the reference-selection obligation precisely.

## The observer algebra and old records

The reference now lives on the full matter/link/two-scalar algebra `M_(4d²)`.
That is a larger declared access choice than the earlier site-0 scalar
observer. The interacting `H` generally moves local scalar observables into
the full work algebra. A restriction of the full Gibbs state to one site
does not, in general, have a modular flow equal to the restriction of this
full Hamiltonian flow. The full-work identity must not be silently reassigned
to every regional observer.

The construction can retain the old record labels explicitly. On

\[
 \mathcal A=\bigoplus_r M_{4d^2},\qquad
 \Sigma_\beta=\bigoplus_r q_r\sigma_\beta,
 \quad q_r>0,\quad\sum_rq_r=1,
\]

the scalar factors `q_r` cancel inside each modular conjugation. The modular
flow acts as `alpha_(-beta s)` in every work block and fixes each central
record projector. The chosen reference record weights may differ from the
actual observed record probabilities; the executed history is not reweighted.
In the receipt the reference weights are `(1/10,2/10,3/10,4/10)`, explicitly
declared and positive.

If `v_r` are the four existing unnormalized terminal work vectors, the actual
recorded state remains `direct_sum_r |v_r><v_r|`. Its branches may be pure and
nonfaithful on the full matrix blocks. This causes no logarithm problem for
the fixed faithful reference: we take `log(Sigma_beta)`, not the logarithm of
those branch densities, and add no artificial noise to the actual state.

The earlier repair/readback construction uses a tracial work reference. Its
matrix density is `I/D`, whose modular flow is the identity. It cannot generate
the present nontrivial `H` through this identification. The new Gibbs clock
reference is a different declared reference; the existing repair process has
not been shown to generate or select it. In particular, preservation of the
old tracial reference does not imply preservation of the new Gibbs reference.

## A detector that actually changes

[thermal_clock.py](thermal_clock.py) starts from the unchanged complete
sharp-readout history in [experiment.py](experiment.py). For each of the
declared appended model times `t=0,1/2,1`, it continues a fresh copy of every
existing terminal branch with the common source evolution and the actual
source gauge-number pointer on bit zero. The modular label is `s=-t/beta`.
There are eight retained outcomes: the four old record pairs and both new
pointer bits. No branch is selected or normalized to make the detector tick.

The Gibbs reference supplies a second description of these unitary segments;
it is not a fresh preparation. The source pointer probabilities are computed
from the actual continued work vectors. Their comparison with the modular
and continuous-Hamiltonian expressions tests the alignment, and their
variation tests that this preparation and detector provide a nonconstant
finite readout.

These settings are separate experiments on freshly continued complete
histories. They are not sequential measurements on one copy of the system.
A nonconstant probability response alone supplies neither an exact time from
one random bit nor a globally invertible or autonomous clock. Observation
windows, statistical inference, state preparation, control services and a
physical calibration remain separate operational requirements.

## Continuous identity and executed product words

The Gibbs/modular equation is exact for the finite continuous generator `H`.
The source program executes symmetric product words, which approximate
`exp(-itH)` at a declared number of steps. They are not asserted to equal
that exponential at a finite step count. Their all-input ideal discrepancy
is bounded by the existing finite temporal certificate; appending a pointer
and retaining the classical history contracts that error.

Numerical source matrices and the reference spectral calculation also carry
floating-point defects. The receipt reports the source guard/comparison
diagnostics separately from the ideal product-formula bound. Finite matrix
positivity and logarithm checks use their stated numerical thresholds;
they are not outward spectral enclosures for arbitrary temperatures or
occupation limits. Coherent source services, storage and pointer readout
remain supplied physical services.

No occupation or spatial continuum theorem for this Gibbs reference is being
added here. The old sharp-readout and the separate smooth-readout occupation
results retain their own precisely stated protocols and limits.

## Why the previous result is still correct

The prior [flow comparison](sanity.py) chooses the actual terminal site-0
state `rho0` as its modular reference on `M_d`. Its own modular flow preserves
that state. But continuation of the actual full work state under `H` gives a
nonzero derivative of the site-0 marginal and of a local observable. No
constant rescaling can identify those two flows on all local expectations
for that state. This remains a valid obstruction to that particular choice.

The existing [programmed modular clock](modular_clock.py) also remains valid:
it perturbs the actual site-0 ensemble and implements the modular flow of
its declared local reference through an additional source pulse program.
It was not a claim that the original interacting `H` already generated that
local modular flow.

The present construction changes both the reference and the accessible
algebra. It aligns a fixed Gibbs reference on the full work algebra with the
same `H`, while letting the actual state change. It neither erases the old
counterexample nor requires the actual observation state to be Gibbs.

## Interventions and the completed history

The common experiment has two Hamiltonian segments separated by record
writing and conditional feedback. The **same fixed reference** describes
both continuous `H` segments through `s=-t/beta`. The gauge pointer, conditional
pulse and final detector remain their actual intervention maps between those
segments.

Measurement with committed outcome records is generally not an automorphism,
and a conditional pulse need not belong to this fixed one-parameter group.
Accordingly the entire prepared/measured/controlled history is not claimed
to be one uninterrupted modular unitary evolution. Its event ordering and
instruments remain part of the declared observer protocol, with old records
fixed during subsequent blockwise evolution.

Arbitrary gauge recording and conditional feedback need not preserve
`sigma_beta`. No claim is made that all existing repair, writing or control
channels have become Gibbs-preserving. The exact common-reference identity
applies to the between-event `H` segments; it does not supply a complete Lean
observer or an inhabited common-reference physical tower for all operations.

The [history repair construction](HISTORY_REPAIR.md) solves compatibility of
the complete branching history; it does not force its actual work state to
equal the Gibbs reference or turn its intermediate event states into one
stationary instantaneous state. Given the same action and chosen reference,
the modular labeling is also independent of the auxiliary solver schedule.
The host repair counter, Hamiltonian model coordinate and modular parameter
continue to denote different things, with only the explicit conditional
relation above asserted between the latter two.

## Scope of the repair

For any specified finite Hamiltonian one can choose such a faithful Gibbs
reference. This makes the reference-flow correspondence concrete and useful,
but also explains its limitation: it does not show that OPH selects this
Hamiltonian, this algebra, this `beta` or this reference from self-consistency.
It is not a proof of Gibbs equilibration or autonomous reference preparation,
and it does not produce an observer clock merely by naming an algebraic flow.

Removing reference selection as a supplied premise would require the actual
construction to produce its reference and establish the claimed relation to
the observed dynamics. Alternatively one may choose a faithful reference and
define `H=-(1/beta)log(sigma)` up to an additive scalar as a hardware design.
That is an allowed construction choice, but it repackages the chosen input;
it does not independently derive the dynamics or agreement with nature. The
present bridge does not close this reference-selection gap.

The demonstrated outcome is a fixed-reference modular description of the
same finite changing work process, with a continued source detector and
retained records. BW identification, physical clock calibration, spatial
propagation and agreement with nature remain open attachments.
