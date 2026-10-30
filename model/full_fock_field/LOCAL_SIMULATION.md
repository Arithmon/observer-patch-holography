# Local observer simulation of the interacting field

The existing full-Fock source schedule has a volume-linear error coefficient,
and a local observation needs only a finite padded part of the field. Both
statements apply to the supplied Weyl/contact action at fixed spatial spacing,
with its Hamiltonian model-time parameter. Neither chooses a physical action,
prepares a vacuum/reference, nor takes a chiral spatial continuum limit.

## The actual source schedule and a smaller split bound

List the existing hopping terms in `FieldGeometry.groups` order, followed by
the onsite contact terms. A hopping term has even CAR support on two modes
and norm `c/2`; contact has even two-mode support and norm `U`. Contacts on
different sites commute, so the two forward contact halves in `evolve` equal
the contact part of this symmetric ordered product. No source word changes.

For `H=sum_i H_i` and `B_i=sum_(j>i) H_j`, the finite-dimensional unitary
bound of Childs, Su, Tran, Wiebe and Zhu,
[Proposition 16, equation (152)](https://arxiv.org/pdf/1912.08854v3#page=39), is

```
||S(t)-exp(-itH)|| <= |t|^3 sum_i {
  ||[H_i,[H_i,B_i]]||/24 + ||[B_i,[B_i,H_i]]||/12 }.
```

It holds at finite time, rather than only for a truncated series. Iterating
unitary steps divides the cubic coefficient by `n^2`. Negative durations
follow by taking adjoints.

Disjoint even CAR terms commute. In the first double commutator only pairs
`i<j` whose supports overlap contribute. In the second, also require `k>i`
to meet the union of the supports of `i,j`. Each remaining norm is bounded
by four times the product of its three term norms. The executable ordered
support count uses the actual bounded field catalogue and exact rational
arithmetic. It counts possible contributions, without claiming every
overlapping commutator is nonzero. No Fock matrix is allocated.

There is also a closed bound for arbitrary finite rectangular boxes. Each
mode meets at most six hopping terms and one contact; their norm sum is
`W=3c+U`. A hopping support meets weight at most `2W`; after an overlapping
second term, its union has at most three modes and meets weight at most
`3W`. The outer contact terms contribute zero because their suffix consists
only of commuting contacts. With `E` spatial edges, the resulting coefficient is

```
K <= E*c*W*(c/6 + 2W),
||S(t/n)^n-exp(-itH)|| <= K*|t|^3/n^2.
```

This is linear in volume for fixed `c,U`, and vanishes when hopping is absent.
The norm-sum coefficient `(cE+UV)^3/3` remains valid but is unnecessary here.
The improvement follows from locality on the complete Fock algebra, not from
a fixed particle count, Gaussian approximation or a momentum-sector cutoff.

## Complete local histories by nested buffers

Retain the factorial tail in
[the locality construction](README.md#local-thermodynamic-dynamics-at-fixed-spacing),
rather than weakening it to an exponential in distance. With `d=R-r>=1`,
the same Duhamel comparison gives the sharper boundary estimate

```
delta(R,r,s,t) <= min(2, s*(2R+1)^2/4 * sum_(n>=d+1) (24c|t|)^n/n!).
```

Indeed, integrate the existing commutator tail from zero to `|t|`; the factor
`c` in each of at most `6(2R+1)^2` crossing bonds cancels the factor `24c`
from that integral. This also gives zero at `c=0` by continuity. The code
bounds the remaining tail by an exact finite sum and a geometric remainder
after successive-term ratios are at most one half.

This bounds the operator difference for every norm-one observable on `s` sites
inside the radius-`r` cube, comparing radius `R` with any containing box.
The same bound holds with an arbitrary untouched reference tensor factor:
only field bonds enter the commutator chains. Hence it is a completely bounded
norm estimate, and duality gives half-diamond error at most `delta/2` for
the channel whose output is the declared local CAR algebra. This compares
restricted outputs of common input states; it is not a distance between
global states of different volumes.

For one terminal observation, add the finite-box split error `eta`; the
complete output channel error is at most `delta/2+eta`. Any terminal
instrument and retained reference are covered by contractivity.

A repeated history requires more than reusing a single-observable bound.
Fix a finite chronological list `t_1,...,t_m`, a terminal observer support
`X`, and local instruments on `X` and retained observer memory. Use even
local Kraus operators, so their extensions to the CAR field are local;
number-preserving instruments additionally preserve total particle number.
All raw outcomes and memory remain in the output, with no conditioning on
success. Initial field states may be correlated with the exterior and with
arbitrary references.

Construct nested cubes backwards. The last evolution uses a buffer around
`X`; each earlier evolution uses a new buffer around the **whole cube**
needed by the later approximate observable. Forward execution therefore
evolves successively smaller cubes, restricting to the next cube after each
segment, while retaining observer records and memory. The initial input is
the restriction of the common ambient state to the largest cube. It is a
supplied input, not an automatically prepared state.

To prove the history bound, telescope the Heisenberg composition backwards.
At each replacement the later approximate effect belongs to precisely the
smaller declared cube, tensored with the retained memory/reference. Apply its
boundary estimate there; all preceding dual channels are unital completely
positive contractions. The same argument at every matrix level and channel
duality give

```
half_diamond(complete histories) <= sum_i (delta_i/2 + eta_i).
```

There is no outcome-count multiplier or division by a success probability.
The buffers and source schedules depend on local support, the finite history,
coefficients and requested error, and not on ambient volume. The same bound
passes to the fixed-spacing thermodynamic dynamics. It does not claim an
unrestricted global-field channel approximation with volume-independent cost.

## Executable budget and explicit limits

`local_observer_schedule` allocates equal error budgets to the segments,
chooses strict buffers using outward rational exponential bounds, and chooses
the smallest positive split count satisfying its chosen coefficient bound.
It returns exact errors and integer source counts. The local observer's
initial restricted state, restriction services, interventions and readout must
still be supplied and charged separately. Floating and hardware errors are
outside the ideal bound. Source pulse counts do not define physical time.

The exact planner admits radii through 4096 and `48c|t_i|<=64`; time subdivision
is required outside that exponential-enclosure domain. Larger mathematical
buffers exist but are not claims about this bounded calculator. The ordered
support enumerator uses the existing `FieldGeometry` limit of 4096 sites.
Its closed degree bound and integer counts also describe larger finite boxes
without constructing geometry lists or exponential arrays. Returned flags
separate those mathematical schedules from the geometry API and the explicit
16-mode, `2^20`-entry state/reference evaluator. These flags concern shapes;
floating parameter conversion and numerical admission must still be checked.
Even an admitted geometry does not promise that its full Fock state can be
materialized.

For site-order strides `d_j=product_(k>j)L_k`, a forward traversal has
`sum_j E_j(8d_j-4)` FSWAPs. A symmetric `n`-step schedule has `4nE` hopping
blocks, `2nV` contact blocks, and twice `n` times that FSWAP count. Completed
hopping/FSWAP words each use 22 single-code and four pair-code pulses;
contact uses eight and two. These count formulas are checked against the
actual source templates and existing routed resource counter.

For `c=1/8,U=1`, the ordered coefficient on `2x2x2` is `629/256`, over one
hundred times smaller than the norm-sum cubic coefficient. A one-site local
observation at `t=1/8` with requested half-diamond error `1/100` uses radius
four and 14 steps under this planner, giving an exact bound below `0.005479`.
Its 340,803,288 pulses describe a large supplied coherent schedule, not a
small executed state array: the cube has 729 sites and 1,458 field qubits.
No ambient-field volume enters that schedule or bound.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/full_fock_field/test_local_simulation.py
```
