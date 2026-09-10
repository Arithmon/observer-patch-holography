# Conditional Cartan scalar and current packet

This package joins a local classical reduction of the pinned SM jet action
to a small charged scalar/Maxwell execution on the prepared q5 golden
addresses. The observer-like sites have scalar state, incident connection
ports, current readback, versioned records and recurrent electric feedback.
The addresses and action are supplied; native source selection is not claimed.

Run from the repository root:

```sh
python3 code/sm_abelian_reduction/build_abelian.py --check
python3 code/sm_abelian_reduction/verify_abelian.py
python3 -m pytest -q code/sm_abelian_reduction/test_abelian.py
```

`build_abelian.py --write` regenerates the receipt. The verifier imports no
producer. It independently checks trace normalization, symbolic currents,
the frozen Yukawa monomial census, golden geometry and source incidence,
then replays every local operation with 60-decimal arithmetic. Source pins
identify both this package and immutable parents. `validate_custody` checks
only schema/scope/source custody and explicitly returns
`mathematical_replay=False`; it cannot substitute for `verify`.

The receipt separates three preparations: baseline, a charge-neutral local
phase intervention with connections fixed, and a gauge copy changing both
matter phase and incident real connections. It records two completed
Strang steps over dimensionless model time [0,1/100]. A first-kick checkpoint
is an intermediate split state, not a new exact physical timestamp.
Static geometry is immutable law input; each event lists all dynamic ports
read and written, their versions, actual writer IDs and semantic parents.
The serial hash chain supplies integrity, not a signal-order definition.
`source_seam_currents` encode the algebraic address preparation. They are
distinct from the dynamical Maxwell current and do not record an executed
native seam-preparation sequence.

The exact first-edge current 2/125 and feedback 1/25000 concern ideal
arithmetic. Whole-history fields and readouts are binary64, checked within
an explicitly numerical tolerance; there is no validated roundoff or
continuous-time error enclosure. Energy is reported, not asserted exactly
conserved by the split. The producer's reversed-potential-order difference
is a separately marked diagnostic; independent tests check disjoint subflows.

The reduction is local and classical, includes all discarded equations and
retains arbitrary Yukawa coefficients at zero fermion fields. The discrete
parent uses a logarithmic plaquette chart, not a Wilson-cosine action. Real
unwrapped connections are essential; generic kinetic coefficients do not
select a compact subgroup of the full charge lattice. No global gauge,
quantum truncation, electromagnetic identification, laboratory units,
count-clock attachment or q233 detector bound is inherited.
