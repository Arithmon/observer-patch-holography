# Complete-field modular histories

The same finite interacting field now has a source-operated modular flow on
its **complete Fock algebra**, together with complete observation histories.
The construction uses the existing [Gibbs preparation](../field_thermal/README.md),
[field action](../full_fock_field/README.md), source compiler and
[complete export executor](../reference_recovery/README.md). It extends the
earlier two-level clock attachment to the actual one/two-site prepared field.
No particle-number sector is selected.

The supplied action, inverse temperature, initialized work bits, coherent
controls and source services remain explicit inputs. This is a finite
implementation result for CP-02, not physical clock calibration, passive
repair dynamics or an accepted integrated completion gate.

## The parameter and sign are fixed by the same reference

Let H be the full finite Weyl/contact Hamiltonian, beta>0, and
rho=exp(-beta H)/Z on the full Fock algebra. Finite dimension makes rho
faithful. Functional calculus gives

```
log rho = -beta H - log(Z) I,
sigma_s(A) = rho^(is) A rho^(-is)
           = exp(-i beta s H) A exp(i beta s H).
```

With the Heisenberg convention tau_t(A)=exp(itH) A exp(-itH), this is exactly
tau_(-beta s)(A). Thus the Schrödinger source word must approximate
U_s=exp(+i beta s H), not its adjoint. Its scalar Gibbs partition phase
cancels in every conjugation. The group law, preservation of rho and common
parameter conversion hold on **every** full-field observable, without
choosing a two-level corner or assigning separate species conversions.
All derivatives, increments and source durations use this convention.

The exact implementation claim is the existing ideal source-word theorem
applied at t=-beta s. If its operator error is eta, then

```
||S_s* A S_s - sigma_s(A)|| <= 2 eta ||A||
```

for the ideal unitary split word S_s. The same argument applies with any
untouched reference because the operator norm is unchanged by tensoring
with its identity. The implementation retains an exact rational split
bound from the actual support-overlap count; it does not estimate this
bound from a fitted trajectory. Rational inputs specify s and beta. The
compiled angles and matrix comparisons use floating arithmetic, whose
residuals are separately checked, not outward-certified error bars.

## Complete channels and observed histories

`certificate` executes every input basis column on the same native word.
Write S for its successful output map and F_j for its first-failure Kraus
maps, including all complement basis values in each coarse abort record.
The exact source/export model obeys

```
S* S + sum_j F_j* F_j = I.
```

The executor freezes every full six-dimensional source primitive and the
word, retains each failure weight and its full input Gram matrix, and
reconstructs every failure map on request. This avoids materializing a
large failure tree without dropping its quantum outputs. Failure Kraus
branches sharing an abort record are combined **incoherently**. The sender
resets on export. A failure terminates in its retained record; successful
outputs are never renormalized.

Embed an ideal unitary U in the success output and zeros in the abort
outputs. A Stinespring comparison gives the conservative bound

```
||V_actual - V_ideal||
 <= sqrt(||S-U||^2 + ||sum_j F_j* F_j||).
```

For a complete basis this controls every input vector, including every
untouched reference, rather than just the thermal state. Its numerical
value is reported as a diagnostic; the algebraic completeness identity and
ideal splitting bound have their separate mathematical proofs. Taking a
partial trace of the retained dilation cannot increase trace distance.

`history` composes one to four modular increments with complete occupation
reads. Each read appends a distinct initialized pointer bit to bounded
working memory. The two ideal maps are P_0 U_s and P_1 U_s. For m reads,
every one of the 2^m histories is retained as an unnormalized operator, in
chronological bit order. No equal-outcome history is privileged. Replacing
one unitary at a time in this isometric composition bounds the entire
history error by the **sum** of its segment errors, independently of the
number of outcomes. Noncommuting reads disturb the state, so this is not
a claim that every observed history remains in the Gibbs reference.

`prepared_readout` executes Gibbs preparation, modular continuation and a
number read in one source word. It retains the field purifier and both
pointer branches. The ideal marginal before the read is rho; the ideal
read outputs are P_b rho P_b. Their sum generally differs from rho. The
joint purification is pure and is not a faithful state on the complete
field-plus-purifier algebra. Faithfulness belongs to the field restriction.

The execution domain is one/two two-component sites, at most sixteen
field basis states, up to four field-only history pointers or one pointer
with the complete purifier, and |beta s|<=4. Each increment uses at most
128 subdivisions. The existing executor's pulse limit also applies.
Resource reports count all native pulses and load/export events, including
zero-angle requests. Initialized bits, retained records, the purifier and
abort tape are accounted separately. No host runtime, event count, memory
cost or pulse count is assigned a physical duration or energy.

## Regional clocks are a separate statement

The [reference observation tower](../reference_observation/README.md)
does not in general intertwine the global modular flow with the modular
flow of a correlated reduced state. The tests retain a counterexample on
this same interacting Gibbs state. The exact
[regional retraction obstruction](../regional_modular/README.md) remains
intact. Implementing the global automorphism therefore does not supply
geometric modular diamonds, a common physical clock, arbitrary regional
modular covariance, a continuum limit or the Einstein premises.

The tracial control sets this same action to zero and gives identity
modular evolution. A nonzero interacting action gives nontrivial evolution
of an occupation observable. Independent logarithms of the full density
check the sign and the action on complex observables; complex reference
columns, all read outcomes, deliberate source-sign corruption and a leaking
unitary test prevent a stationary-state-only or postselected false pass.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.modular_flow.flow
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/modular_flow
```
