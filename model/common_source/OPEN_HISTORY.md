# A complete observer with unobserved environment operations

The canonical specification is still [`protocol.Program`](COMPOSITION.md).
`open_protocol.py` extends its state representation to a density matrix for
each visible record, so that finite unobserved operations can compose with
preparation, disturbing reads, feedback and continuation. It reuses the same
source preparation, current kick and operation maps. An empty environment
reproduces the existing pure histories on their entire quantum outputs.

The completed connection is between the finite source phase environment in
`anharmonic_transfer/environment.py`, the common four-read observer, its
occupation/temporal bounds, and its corrected native memory. It does not
identify the two-site observer with the trapped-ion experiment: the source
operation is reused, while its variance here is a supplied dimensionless
parameter. The published-data package retains its own calibration contract.

## State, operations and accessible records

For dimension d=2 or 4 the work dimension is 4d², with index
`gauge + 4*(n0 + d*n1)`. An optional finite reference is the fastest index.
The state after a visible prefix is the classical/quantum direct sum

\[
  \Omega=\bigoplus_r \rho_r,\qquad \rho_r\succeq0,
  \quad\sum_r\operatorname{tr}\rho_r=1.
\]

Blocks are unnormalized. Zero and tiny outcomes remain in the record menu.
With no supplied state, the declared vacuum is prepared once and receives
the native current kick. A supplied work/reference density replaces that
preparation and receives the same kick; it does not inherit a vacuum bound.
The supplied state must be normalized, finite, Hermitian and positive to
the runtime's stated numerical tolerance. Those checks are diagnostics, not
certificates of exact positivity of arbitrary floating input.

Each original operation acts as `K rho K*`, with the identity on the
reference. A read appends its bit; feedback may use only an already committed
bit. The environment is an ordered schedule of `(after, site, variance,
order)` entries, unique per boundary/site. Boundary zero is after preparation
and the current kick. At each entry the runtime applies

\[
 \mathcal E(\rho)=\sum_{k=1}^M w_k
       (e^{-i\theta_k N_s}\otimes I_R)\rho
       (e^{i\theta_k N_s}\otimes I_R).
\]

Weights and executed phase words come from the existing finite source
channel. Every invocation has a fresh independent hidden alternative.
The Schur contraction used to sum these executed words is the same CP
map; it is not an average of amplitudes. Draws are discarded, never added
to the visible tape and never available to feedback. Reusing one hidden
draw across boundaries defines a different process and is a retained
negative control.

The runtime accepts the original finite grammar, including sharp reads.
The occupation theorem still rejects sharp-read histories. A finite source
execution and a valid occupation certificate have different domains.
The bound of eight reads and 64 operations is unchanged; the density
implementation additionally refuses a complete output exceeding 2,000,000
complex entries. It checks the requested reference/record allocation before
source operations. Extreme exact phase certificates also have an explicit
serialization budget, described in [OPEN_ENVIRONMENT.md](OPEN_ENVIRONMENT.md).

## Composition and persistence

For each prefix, summing its complete next operation gives a trace-preserving
map: the old instruments satisfy their completeness identity, and each
environment map is a convex mixture of unitary conjugations. Induction
therefore proves completeness of every prefix, including an arbitrary
reference. No independence of successive visible outcomes is assumed.
Summing future outcomes recovers the probabilities of every older visible
prefix. It generally does not recover the earlier quantum state: subsequent
operations can disturb it.

The environment applied at cursor k belongs to that cursor's completed
state. Resumption starts at operation k+1, without repeating the environment,
kick or preparation. JSON checkpoints preserve the complete density blocks,
program, environment, backend, chronology and resources. Import reexecutes
the prefix and compares its output; a replacement digest does not make a
forged density a valid history. In-memory objects carry integrity checks
and are detached before continuation. Neither interface measures, clones
or reconstructs a physical unknown quantum state from classical records.

## Connected results and comparison targets

| Existing result | Instantiation and consequence |
| --- | --- |
| Canonical source operations, `protocol.py` | The same current kick, Hamiltonian words and read/feedback maps act on each mixed block. All visible records remain. |
| Finite phase source, `anharmonic_transfer/environment.py` | Executed local phase words supply a multi-Kraus operation. Source tables and Schur contraction are checked before use; the analytic Gaussian is a separate comparison backend. |
| Complete temporal bound, `continuation_bounds.py` | Contractivity of each intervening CPTP operation telescopes the existing temporal errors through noise and later reads, for arbitrary finite work/reference inputs. |
| Uniform occupation bound, `continuation_bounds.py` | Every phase commutes with the fourth-number moment and cutoff projection. Its insertion adds zero moment growth and zero projection leakage, so the same recurrence bounds the matched noisy infinite history. |
| Prepared output enclosure, `smooth_history.py` | Recomputed for the original noiseless preparation/history, then combined with the finite noise-to-identity error. This target is the **noiseless** infinite history. |
| Corrected native memory, `ENCODED_MEMORY.md` | The complete work instrument is CPTP and acts only on work. The existing logical identity theorem therefore applies even though work outputs are mixed. [OPEN_NATIVE.md](OPEN_NATIVE.md) supplies the actual source attachment and passive control. |

For J environment invocations and their half-diamond bounds εj, the finite
source-to-continuous-Gaussian comparison is bounded by
`min(1, temporal + sum εj)`. The finite-to-infinite matched Gaussian
comparison additionally charges the occupation error. The latter remains
vacuous at the executed cutoffs under the uniform theorem; it is not silently
replaced by the prepared noiseless bound. The separate nonvacuous prepared
comparison is

\[
 D(\text{noisy finite},\text{noiseless infinite})
 \le D(\text{noiseless finite},\text{noiseless infinite})
       +\sum_j \|\mathcal E_j-I\|_\diamond/2.
\]

All these are statements about the ideal source maps and complete outputs,
not certified errors of the floating simulation. The inherited service
ledger requires every work stage and every environment invocation before
returning a physical total. Native attachment additionally requires its
input, encoder, storage, repair/control and recovery services. Missing
coverage stays unknown; neither roundoff nor service error is assigned zero.

At growing cutoff d, an exact Gaussian variance v can be approximated with
`M >= max(d, ceil(2*v*(d-1)^2))`. The factorial remainder obeys
`d*M!*[v*(d-1)^2]^M / (2*(2*M)!) <= d*2^(-M-1)`.
For a fixed finite schedule this tends to zero. The reported bound rounds
outward on a grid `1/(10^12*(d+2)^6)`, whose added allowance also tends to
zero. Combined with the existing occupation/temporal family this gives a
matched noisy ideal refinement; it is an analytic family, not a claim to
execute d>4 or orders above 64. Binary64 phase parameters additionally need
the explicitly charged conversion error; fixed-order noise is not asserted
to converge uniformly as d increases.

## Reproduce and audit

```sh
python -m model.common_source.open_observer --check model/common_source/open_observer_results.json
python -m pytest -q -W error model/common_source -k open_
```

The receipt executes d=2,4, all three native currents, 512 product steps per
interval and four site-zero channels of variance 1/64 and order 8, after
operations 1,4,6,9. It retains all 16 visible outputs, verifies checkpoint
continuation and earlier record probabilities, compares a separate dense
continuous/Gaussian history, and distinguishes fresh draws from reused draws
and from no environment. C/E protected and passive native attachments retain
their entire raw outcome menus, with encoded Bell and phase comparisons.
The receipt contains compact summaries; replay reconstructs the full outputs.

Resource ledgers distinguish the realized path from expansion of its entire
visible tree. Each continuing path still has one work preparation and four
pointer blanks; four discarded environment draws add eight native phase
pulses at d=2 and twelve at d=4. A hidden finite menu is a randomized choice of one phase word, not
the sequential physical execution of all its alternatives. The native
ledger adds the actual encoder, repair and selected recovery words and
reports a compatible path attaining its pulse maximum. Host density storage
and unquantified physical services are outside pulse counts.

Independent tests compare full complex blocks with direct matrix formulas,
including a flagged Choi input and a nontrivial reference. Wrong tensor
order, a selected hidden draw, and reversed feedback can all preserve trace;
the oracle rejects each. Tests also cover pure-history agreement, nonzero
tiny branches, original-type validation, allocation refusal, mutable caller
aliases, resealed checkpoint forgeries, complete native endpoints and exact
resource counts. Their floating residuals are not substituted for ideal
proof bounds.

The Hamiltonian model time remains 5/2. Pointer interactions, source pulse
counts and supplied boundary variances are separate quantities. No physical
clock conversion or simultaneous Hamiltonian/Lindblad evolution is inferred.
Physical calibration, noise on native memory/records, repeated repair cycles,
and the broader effective-field/gravity target remain obligations of the
[construction guide](../CONSTRUCTION.md). This PR completes an executable
finite composition and its comparison contracts, not those physical claims.
