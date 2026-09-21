# SOURCE-CURRENT-TOMOGRAPHY-0

## Purpose

This packet asks a deliberately narrower question than the conditional
`PORT-CURRENT-INNER` construction:

> Do the registered source-bound artifacts contain enough
> information to reconstruct the current generators, their bracket, and
> same-current overlap holonomy from ordered response histories?

The Stage 0 answer is **no**. This is a fail-closed information audit, not
a no-go theorem for other source data and not a rejection of the existing
conditional port-current algebra.

## Source-bound data

The charged-response semantic artifact fixes, target-blind:

- the twelve-port carrier binding;
- the maximal-distance antipodal response;
- the four relative sector signs on `1 + 5 + 3 + 3'`;
- the orientation convention;
- the finite refinement persistence maps.

These facts do not reconstruct the current.

## Absent reconstruction data

The present artifacts do not contain:

1. ordered two-sided port-response histories;
2. twelve source-reconstructed skew-adjoint infinitesimal generators;
3. the positive pairing reconstructed from those histories;
4. a commutator reconstructed from the same histories;
5. closed overlap words built from the same source packet;
6. same-current projective implementers establishing internal holonomy.

Those are the operative source-tomography obligations inherited from the
scientific boundary of the former issue #566.

## Exact negative control

The existing B14 oriented-face construction is a useful control against an
easy mistake.

From the twenty oriented faces it forms the equal-weight cyclic bracket

`[e_a,e_b]=e_c`, `[e_b,e_c]=e_a`, `[e_c,e_a]=e_b`.

This construction respects the committed incidence structure and proper
order-60 action, but its Jacobi tensor has exactly 240 nonzero coordinates
and squared norm 240.

Therefore incidence naturality and A5 compatibility alone do not force a Lie
current.

B14 also compares this failed bracket to compact bracket families, but that
comparison is explicitly conditional on a Euclidean/Hilbert-Schmidt metric
and a minimum-distance repair rule that are not source-derived. It therefore
cannot be promoted into the missing source selector.

## Stage-0 verdict

```text
INSUFFICIENT_SOURCE_DATA
```

This verdict means only that the registered source-bound packet does not
determine the physical current at the required altitude.

## Verification

```bash
python3 code/a5_closure/source_current_tomography_stage0.py verify

python3 code/a5_closure/tests/test_source_current_tomography_stage0.py
```

The audit intentionally consumes existing committed artifacts rather than
introducing a replacement current fixture.

## Admissible reconstruction input

A positive reconstruction requires a versioned source packet of ordered
two-sided response histories. Without gauge-group or particle labels, the
contract is:

```text
ordered histories
    -> tangent response operators
    -> 12-dimensional source span
    -> exact pairing
    -> commutator closure
    -> overlap words
    -> projective implementers modulo centralizer
```

Two scientifically useful exits satisfy this contract:

1. the reconstruction is unique up to the explicitly allowed basis/gauge
   equivalence and matches the conditional current fixture after the fact; or
2. two inequivalent current realizations survive the same source packet, in
   which case the packet emits an explicit non-identifiability certificate.

No downstream Standard Model label, coupling, mass, or laboratory target may
participate in the selector.
