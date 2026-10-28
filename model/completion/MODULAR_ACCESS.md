# Modular action on the constructed reference: an access obstruction

CP-02 / [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138).
The actual finite observer has a stationary reference, but the frozen native
operations cannot implement its nontrivial modular action on the fixed full
matrix algebra of a nontracial record block. An exact test input gives a
trace-distance implementation error of at least **1/2**. This is a bounded
algebra/access obstruction, not a rejection of every possible observer clock.

The closest existing result is the
[native real-operation boundary](../native_interface/README.md#the-exact-reachability-boundary).
That result already proves the admitted native grammar preserves diagonal
inputs and, even with a specified real mixer, cannot distinguish complex
conjugate states by its real effects. The addition here is its attachment to
the **actual full-record stationary reference**, with a specific modular
parameter and an exact implementation-error witness. No new general theorem
about real quantum operations is claimed.

## Fixed algebra and actual reference

Use the [stationary full-record observer](observer_reference_results.json)
at repair probability one, uniform exploration probability one half and one
retained FIFO slot. The classical record labels define direct-sum blocks.
Within a positive block restrict the source space explicitly to the positive
support of its diagonal reference. Let its conditional density be `sigma`.
It is faithful on that support. This is a mathematical full matrix algebra;
the campaign has not established observer access to every one of its effects.

The [retained witness](modular_access_results.json) locates unequal positive
weights inside the same actual block for each endpoint:

| Law | Complete record | Positive support | Source pair in the fixed H78 ordering |
| --- | --- | ---: | --- |
| C | explore seam `(0,1)`, coin `0`, current `0` | 56 | `(76,61)` |
| E | explore seam `(0,1)`, coin `0`, contrast `-1`, decoded current `0` | 10 | `(59,61)` |

The exact rational weights and their ratio are retained, together with the
block's canonical hash. Choose the ordering so `q=p_i/p_j>1`. Conditioning
on the block divides both weights by the same positive record probability,
which cancels from `q`.

## The real-operation invariant

Every complete C/E branch is a positive real scalar times a native partial
injection. Identity holds and halts are real. A classical memory update,
record export or classical random choice has a real basis-map realization.
The declared diagonal independent ancillas are real joint states.

Consequently, every finite admitted history sends a real joint density
matrix to a real subnormalized density matrix. Classical sums and partial
traces preserve this property. Conditioning on an outcome of positive
probability divides by a positive real number, so a heralded branch remains
real. The same condition is closed under trace-norm limits. The proof uses
the operational grammar and all native branches, rather than a sampled
list of histories. The implementation checks every native primitive support
and separately replays signed real coherent probes through the actual
observer's branches and adaptive continuation.

The ancillary premise concerns the **joint** ancillary state. Individually
real marginals would not exclude hidden complex joint correlations. A new
complex phase service is not part of the frozen campaign.

## An exact modular witness

Use the convention

```text
sigma_s(A) = sigma^(i s) A sigma^(-i s),
sigma_s,*(rho) = sigma^(-i s) rho sigma^(i s).
```

The second expression is the state transformation dual to the first. For
the real pair probe `|+>=(|i>+|j>)/sqrt(2)`, choose the **dimensionless modular
parameter**

```text
s_star = pi / (2 log(q)).
```

Its transformed off-diagonal entry is
`exp(-i s_star log(q))/2 = -i/2`. In the ordered pair `(i,j)` the target is

```text
rho_target = [[1/2, -i/2],
              [i/2, 1/2]].
```

Let `Y=[[0,-i],[i,0]]` on that pair and zero elsewhere, including other
record blocks. It has operator norm one. The target has `Tr(Y rho_target)=1`;
every normalized real density matrix has `Tr(Y rho_real)=0`. Trace-norm
duality gives

```text
dist(rho_target,rho_real)
  = (1/2)||rho_target-rho_real||_1 >= 1/2.
```

The lower bound also applies if a proposed implementation leaks outside the
pair, changes the classical record or conditions on a positive-probability
native history. The real state `diag(1/2,1/2)` attains distance one half, so
the real-state comparison is sharp. This does not assert that every
minimizing real state is natively preparable.

An explicitly excluded control `diag(1,i)` maps the probe to the target
exactly. It demonstrates the load-bearing restriction: a complex phase
resource breaks the real-operation invariant.

## What this decides

The admitted operations cannot implement the stated nontrivial modular map
on this **fixed source-basis full support algebra**. The probe tests an
all-input map on that algebra; it is not claimed to be preparable from the
campaign's diagonal initial data. `Y` is a mathematical separation witness,
not an executed or admitted native detector. Restricting the comparison to
native classical record effects changes the contract and need not expose
this quantum-state error.

On the accessible commutative record algebra, the modular automorphism flow
is trivial and identity implementation remains available. A different
quantum encoding, extra coherent registers, an added complex service,
different reference/algebra or carrier refinement requires another map and
is outside this exclusion. Neither a physical-clock no-go nor rejection of
the broader executable hypothesis follows. The CP-02 gate remains open.

Reproduce from the repository root:

```bash
python -m model.completion.modular_access --write
python -m model.completion.modular_access --check
python -m unittest model.completion.test_modular_access
```
