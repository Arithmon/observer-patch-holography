# Certifying the computed output of the continuing observer

The [prepared-history enclosure](SMOOTH_HISTORY.md) already bounds the exact
finite continuous history and its infinite-occupation counterpart. Its
comparison with the source executor was a floating diagnostic. The output
adapter closes that numerical comparison for a specified execution: it treats
the returned binary amplitudes as exact numbers and encloses their distance
from the independently computed polynomial output.

[`exact_output.py`](exact_output.py) supplies `compare_vectors` and
`compare_densities`; [`output_certificate.py`](output_certificate.py) runs
the actual source histories and their comparison targets. The object being
certified is the complete computed output, with every
committed record and its quantum state. This is a certificate of the accuracy
of those represented numbers for the declared preparation and program. It is
not a hardware error bound, a certificate for every input to the source words,
or a replacement for their execution. The source program still produces the
amplitudes; the polynomial method supplies an independent comparison target.

The [primitive observer runtime](../primitive_observer/README.md) shares these
source maps but executes separate output arrays. The bounds here certify the
arrays freshly replayed by `protocol` and `open_protocol`; the primitive
runtime's floating compatibility diagnostics do not transfer this certificate
to its returned arrays, native transfer reconstruction or laboratory histories.

## Fixed output and comparison target

For a pure-history execution, let the returned vector for record `r` be
`v_r`. Every real and imaginary binary floating component is converted
exactly to a rational number. No decimal formatting, probability floor,
branch selection, or normalization intervenes. The represented operator is

\[
 O=\bigoplus_r |v_r\rangle\langle v_r|,
 \qquad m=\operatorname{tr}O=\sum_r\|v_r\|^2.
\]

It is positive, but its exactly computed mass `m` need not be one. A runtime
normalization check within a floating tolerance does not change this fact.
The certificate retains this mass discrepancy.

The target `rho` is the normalized ideal output of the same prepared
continuous program, either at the declared finite occupation cutoff or in
the infinite Fock space. Operations, durations, coupling, current, feedback
flags, record order and preparation must match the target being certified.
The finite preparation retains the declared total occupation through `d-2`;
it is not silently replaced by a different rectangular vacuum projection.
The infinite preparation has its separate truncation and preparation-tail
bound. The two polynomial constructions have different initial normalizations.

This comparison uses a complete snapshot of the output. A mutable buffer,
altered program or changed record menu cannot redefine that snapshot after
certification. A digest can detect accidental changes, but is not proof that
an externally supplied error claim or state was produced by the executor.
Published execution claims require the actual source run and independently
recomputed certificate, with their source and program provenance.
The retained source bindings include both JSON inputs actually consumed by
the native compiler: its twelve generators and its exact control-coefficient
receipt, as well as the code that validates and executes them.

## Exact embedding of the polynomial output

The existing polynomial enclosure writes a branch as

\[
 u_r=\frac{1}{D_r\sqrt S}
     \sum_{n_0,n_1,g}c_{r,n_0,n_1,g}
          |n_0,n_1,g\rangle_B,
 \qquad
 \|\,|n_0,n_1,g\rangle_B\|^2=n_0!n_1!.
\]

Here the coefficients are Gaussian integers, `D_r` is a positive integer,
and `S` is the exact rational initial squared norm shared by all branches.
Thus its amplitude in the normalized Fock basis is

\[
 (u_r)_{n_0,n_1,g}
   =\frac{c_{r,n_0,n_1,g}}{D_r}
       \sqrt{\frac{n_0!n_1!}{S}}.
\]

The factorials cannot be omitted or absorbed into a floating conversion of
`S`. The finite runtime index is `g+4*(n0+d*n1)`. Its vector is embedded by
zero outside `0 <= n0,n1 < d`; the polynomial remains on its entire retained
support, including occupations outside that rectangle.

For each record define

\[
 p_r=\|v_r\|^2,\qquad
 q_r=\|u_r\|^2,\qquad
 z_r=\langle v_r,u_r\rangle.
\]

Both squared norms are exact rationals. The overlap is a finite sum of
rational complex coefficients times square roots of positive rationals.
Directed rational square-root intervals enclose its real and imaginary
parts, hence `|z_r|^2`. Only the overlap sum is restricted to the runtime
rectangle. In particular, `q_r` includes the entire outside-support mass;
discarding it would certify a different output.

## Complete quantum and classical distance

For arbitrary, possibly unnormalized vectors, the rank-one identity is

\[
 \frac12\bigl\||v\rangle\langle v|
                    -|u\rangle\langle u|\bigr\|_1
 =\frac12\sqrt{(\|v\|^2+\|u\|^2)^2
                         -4|\langle v,u\rangle|^2}.
\]

It follows by restricting the Hermitian difference to the span of the two
vectors. Orthogonal record blocks make the complete norm the sum of these
branch norms. Directed evaluation therefore gives an interval `[L,U]` for

\[
 D(O,P)=\frac12\|O-P\|_1
       =\frac12\sum_r\sqrt{(p_r+q_r)^2-4|z_r|^2},
 \qquad P=\bigoplus_r|u_r\rangle\langle u_r|.
\]

Cauchy--Schwarz proves the true radicand nonnegative. An interval lower
endpoint may be replaced by zero using that proved fact; a negative upper
endpoint instead indicates an invalid enclosure and must be rejected. This
is different from clipping an incorrectly computed exact squared distance.
Zero and arbitrarily small branches remain in the sum.

The earlier enclosure proves

\[
 D(P,\rho)\le E,
 \qquad E=\tau+\epsilon+\epsilon^2/2,
\]

where `epsilon` includes its complete-vector Taylor, support-removal and
coefficient errors, and `tau` is its preparation trace-distance bound. The
triangle and reverse-triangle inequalities then give

\[
 \boxed{\max(0,L-E)\le D(O,\rho)\le U+E.}
\]

There is no cap at one: `O` is a positive represented operator with mass `m`,
not an exactly normalized state. The universal positive-operator upper
bound is `(m+1)/2`, if that separate bound is used. Capping at one is valid
for comparisons of two normalized ideal targets, not for this raw output.

The direct finite-continuous comparison includes the actual discrepancy of
preparation, splitting, source compilation and floating contraction from that
finite target. The direct infinite comparison additionally includes the
occupation/preparation difference. Adding the analytic temporal allowance
again to either direct comparison would count it twice. Those analytic
all-input and asymptotic theorems remain useful separate results; the direct
certificate neither replaces them nor isolates an error for each stage.

A scalar phase on an entire committed branch leaves its rank-one output and
this distance invariant. Relative phases between components of one branch
are observable and affect the distance. Comparing only branch probabilities
would miss that distinction. The certificate concerns committed records; it
does not certify a coherent superposition of alternative record histories.

## Record effects and a matched feedback intervention

Let `F` be an effect, `0 <= F <= I`, on the complete record/work output and
let `eta` bound `D(O,rho)`. The Jordan decomposition of the Hermitian
difference gives

\[
 |\operatorname{tr}(FO)-\operatorname{tr}(F\rho)|
 \le \eta+\frac{|m-1|}{2}.
\]

Its positive and negative trace masses are
`(||O-rho||_1 +/- (m-1))/2`. This proves the displayed bound and explains why
the usual normalized-state probability rule cannot simply use `eta` alone.
For example, `O=(1+delta)*rho` has distance `delta/2` but total-mass error
`delta` when `delta>0`.
The ideal probability interval may be intersected with `[0,1]`, because the
target is normalized. A probability computed from the raw output must itself
be recomputed exactly or enclosed; a separately rounded diagnostic is not
an exact expectation merely because its input amplitudes were certified.

The retained feedback control uses `d=4`, current zero, and sets both
feedback angles to zero. It compares two complete executions with the same
preparation, current, cutoff, source splitting and record menu. Only the
declared feedback enable flag changes. Each execution receives its own
matching infinite polynomial target; the baseline also has a separate finite
continuous comparison. This is an intervention
comparison between two programs, not two outcomes selected from one history.
Checkpoints cannot be transplanted across the changed program identity.
Each zero-angle source request still executes and is charged as one native
pulse with its load/export services; disabling the feedback does not remove
those requests from the resource ledger.

If exact-output comparison gives a lower bound `delta` on `D(O_on,O_off)`
and the respective ideal-target errors are `eta_on,eta_off`, then

\[
 D(\rho_{\rm on},\rho_{\rm off})
 \ge\max(0,\delta-\eta_{\rm on}-\eta_{\rm off}).
\]

A strictly positive result certifies a complete-output feedback effect.
Alternatively, disjoint intervals for a specified readable record effect
certify that effect directly. Neither conclusion follows just from two
different floating diagnostics. The finite target can give a useful
intervention certificate even when a coarse infinite-target error is too
large to resolve that same contrast.

## Mixed represented outputs

A mixed block need not admit a single-vector formula. For a returned finite
matrix `B_r`, retain each original binary component exactly and embed the
matrix by zero outside its finite `D`-dimensional code. The difference from
`|u_r><u_r|` has rank at most `D+1`: its range lies in the finite code plus
the one polynomial vector. Consequently

\[
 \frac12\|B_r-|u_r\rangle\langle u_r|\|_1
 \le\frac{\sqrt{D+1}}2
       \|B_r-|u_r\rangle\langle u_r|\|_F.
\]

If the comparison polynomial lies wholly inside the finite code, `D`
replaces `D+1`. A directed enclosure of the squared Frobenius norm uses

\[
 \|B_r\|_F^2+q_r^2
       -2\operatorname{Re}\langle u_{r,\mathrm{in}},
                                      B_r u_{r,\mathrm{in}}\rangle.
\]

The full `q_r` is retained. The formula and rank bound hold even when the
computed matrix is not exactly Hermitian. Sum the branch bounds and add the
polynomial target error as above. There is no projection onto positive or
Hermitian matrices and no trace renormalization. For such general computed
bytes the safe effect error is `2*eta`; the positive-output mass formula is
not asserted without its hypotheses.

The lower bound uses `||X||_1 >= ||X||_F` and
`||X||_1 >= |Tr X|` for each difference block `X`. The latter retains both
the real and imaginary trace discrepancy. These lower bounds concern the
distance to the polynomial; the complete target error is subtracted before
reporting a lower bound against the normalized ideal history.

Comparing a noisy computed output with a noiseless polynomial certifies its
distance to that noiseless target, including the effect of noise. It does not
certify numerical error relative to a matched noisy infinite process. That
different target requires its own enclosure. The Frobenius estimate can be
loose; a finite result is not automatically a nonvacuous comparison.

The retained mixed executions use the existing canonical finite environment
schedule. The producer first compares their matrices with the exact
**noiseless finite continuous** history, then adds the existing upper bound
between that finite history and the noiseless infinite history. This is a
triangle inequality through two normalized ideal targets; the represented
mixed matrix itself is not normalized or repaired. The finite polynomial
error appears in both comparisons, which is conservative. No Gaussian
infinite-history certificate is inferred.

## Connected results, resources and remaining obligations

This adapter composes the existing factorial polynomial construction,
complete-instrument error induction and rank-one trace-norm identity with
the actual canonical executor. It adds no operation grammar, preparation,
reset, clock or physical service. Hamiltonian time is still the sum of the
declared evolution durations; pointer parameters, native pulses, verification
work and host runtime remain distinct. The source resource ledger describes
the executed history. Rational verification work is additional host work,
not extra pulses in the physical protocol.

All three currents and both executed occupation cutoffs require their own
output comparisons. A maximum of the **pure-output bounds** over currents
can then pass through the existing **ideal** native attachment for independently prepared work and an
arbitrary native/reference input: the first native record chooses a current,
its weights sum to one, and later common CPTP stages contract the Hermitian
trace norm. There is no branch-count multiplier. This statement does not
certify the numerical execution of the surrounding native repair, encoding
or recovery maps. Nor does it cover initially entangled native/work states,
arbitrary work preparations or all-input channel error.

Physical preparation, source controls, proper-code services, record storage,
native logical storage, detector response and clock calibration still need
their own contracts. Their absent bounds remain unknown. In particular, the
new state-specific numerical result must not be entered as zero error for
every service in the physical ledger.

Exact arithmetic certifies the supplied output numbers; it does not make
different platforms produce identical floating numbers. A replay must
recompute a certificate for its own snapshot and validate the published
comparison contract. A stored radius or digest must not be silently attached
to different output bytes. Fixed conservative publication envelopes may be
checked against freshly recomputed exact bounds; their acceptance rule must
remain distinct from exact source/provenance and internal-consistency checks.

Independent controls should exercise factorial weights, outside-code support,
unequal total masses, zero and tiny records, within-branch phase changes and
exact branch-global phase invariance. A mismatched current, feedback flag,
or record allocation must not become an accepted same-program certificate.
The adapter can compare arrays of any origin with its fixed prepared target;
the producer separately establishes that its arrays came from the declared
source run. Nonfinite inputs, mutable aliases and altered certificate
envelopes must fail before a success result is emitted.

## Retained executions and reproduction

The [receipt](output_certificate_results.json) executes the canonical four-read history with 512 product steps
per interval for each current `-1,0,1` and each cutoff `d=2,4`. Every one of
the sixteen records is compared, both for the pure source output and for the
mixed source output with the canonical environment. The pure outputs have
separate finite-continuous and infinite-target certificates. The additional
`d=4`, current-zero execution sets both feedback angles to zero while preserving
all reads and Hamiltonian elapsed time.

The following complete-output upper bounds take the maximum over all three
currents. They concern the represented work output and its named target;
the mixed column includes the effect of the supplied environment.

| Cutoff | Pure output to finite continuous target | Pure output to infinite target | Mixed output to noiseless infinite target |
| --- | ---: | ---: | ---: |
| `d=2` | `0.000001` | `0.263093` | `0.285928` |
| `d=4` | `0.000001` | `0.048204` | `0.136189` |

All entries bound half the complete trace norm. The first column of bounds
includes actual preparation, finite splitting, source compilation and
floating contraction against the exact finite continuous process. The
infinite columns additionally compare different occupation descriptions.
They are not errors of a matched noisy infinite process or physical services.

For the feedback intervention the independently enclosed infinite clock
probability for `r3=1` lies in `[0.688428, 0.688449]` with feedback and
`[0.411894, 0.411920]` without it. The separately enclosed difference lies in
`[0.276509, 0.276555]`. These are outward intervals from the
polynomial calculation, including its preparation and complete-history
errors. The source-output certificates separately bound each executed
clock-event mass against its matching ideal interval. This establishes a
readable feedback effect in the supplied model; it does not calibrate the
last gauge bit as physical elapsed time.

The published numerical guarantees are rounded outward to a rational grid
of `10^-6`. Replay executes the source again, recomputes the exact output
comparison and checks these reported guarantees with strict JSON types.
Raw array hashes and exact binary masses are not asserted to match between
numerical platforms. Source bindings, program specifications, record counts
and resource counts remain exact.

```sh
python -m model.common_source.output_certificate --check model/common_source/output_certificate_results.json
python -m pytest -q -W error model/common_source -k 'exact_output or output_certificate'
```
