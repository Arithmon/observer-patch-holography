# A source realization of coherent codes and an operational clock

The existing twelve-port source has enough native response control to
perform an entangling gate: the sum of its generators is the observable
relative phase `diag(i I3,0)`. Coherent proper-code transport then supplies
the entire massive-clock circuit. The construction includes a full-algebra
instrument, spherical support/refinement, source response completeness,
finite preparation, positive control time and an executed serial readout
instrument, including all no-click branches.

Classical record agreement is insufficient. In the complete code-channel
family, A3 with tracial Choi reference selects dephasing under one-basis
agreement and identity under two complementary-basis agreement. The noisy
optimizer and its refinement semigroup are solved exactly. A classical
source completion obeys the declared A1--A3 grammar and RG but cannot
transmit an unknown qubit: its four complementary-test average is at most
3/4. Coherent code gluing (CCG) is an explicit proposed extra requirement,
not a silently adopted fourth core axiom.

The compiler preserves the operational clock on vacuum and one-particle
states. A finite analytic witness retains probability swing above **.6752**
after preparation, control-time and phase-noise costs. Phase noise has a
population-independent signal bound on that sector, but its rate must
scale as `O(a)` for bounded Floquet accounting energy at refinement.
Population-independent stability does not hold for occupation-flip noise.
Ordinary register flights do not automatically supply fermionic signs on
many-particle states; the executable countercontrol keeps that boundary.

Read the [contract](CONTRACT.md) and
[complete proof and resource ledger](../../extra/COHERENT_SOURCE_CLOCKS.md).
This is an explicit mathematical source completion, not a simulator capture
or selection of empirical masses, physical energy or a unique representation.

## Reproduce

With the repository's pinned `requirements.txt`:

```sh
python code/m1_source_realization/build.py
python code/m1_source_realization/verify.py
python -m pytest -q code/m1_source_realization
```

The checker imports no new producer functions. It reconstructs source
generators from the parent's independent exact formula, verifies the exact
control coefficients without solving the producer's inverse problem,
checks channels against analytic Bell spectra, and checks spatial evolution
against direct path entries. Preparation is compared to every requested
leaf amplitude, including zero subtrees. The entire vacuum/one-particle
detector effect is replayed from its serial destructive instruments, and
integer schedules expose every serialized buffer access. The same readout
with dephasing at every code access loses its phase dependence. Noise is
checked against the independent sector-channel identity. The parent clock's outward interval
certificate is executed, and resource arithmetic is repeated at 65 digits.
The source-support tower is checked as oriented chains and actual retained
factor maps. These finite tests support the analytic all-size proofs;
they are not substituted for them or claimed Lean formalizations.

The receipt is compact and reproducible. It contains all small execution
outputs, not a tape for the enormous analytic witness. Strict JSON types,
catalogues and source/claim commitments reject malformed and substituted
inputs. Matrix residuals have absolute tolerance `2e-10`; other nonzero
values use relative `2e-9`, with `1e-12` for expected zero. Tiny positive
time, rate and probability budgets cannot be replaced with zero.
Hostile controls also exercise the CLI with Python assertions disabled.

## Audit corrections and contract coverage

The audit corrected preparation and readout time that had counted accesses
to one processor as concurrent. Initial blanking, all eight departures at
each tree node and all eight output reads are now charged serially. Exact
integer event counts reject missing operations even when their tiny service
times would fit inside the tolerance on the larger flight time. The signal
lower bound remains above .6752.

The readout is executed as a complete destructive instrument, including
every no-click output. Its removal of the particle required extending the
phase-noise proof to vacuum plus one particle. The A3 constructor now clamps
loose fidelity constraints at uniform guessing and rejects invalid errors.
The source-transfer constructor admits only the declared proper codes. Pure basis preparation is an
explicit irreversible intervention; the selected code service is fixed and
cannot be used as an extra programmable gate. The proof also treats limits
of reversible responses at fixed cutoff and the correct finite-workspace
compression of the parent's accounting observable.

| Contract deliverable | Proof and executable evidence |
| --- | --- |
| Complete source, instruments and refinement | Proof sections 1--3; exact controls, all matrix-unit images, support chains and retained-factor maps |
| Complete classical countermodel | Sections 3--4; all-program entanglement-breaking proof, complementary-read controls and the same detector with dephased code accesses |
| Full finite A3 selection problem | Section 4; Bell-spectrum replay, generic feasible-channel controls, loose-bound and semigroup tests |
| Native gates and charged compilation | Sections 1, 5--6; entangled-code controls, full spatial matrices, nontrivial final coin phases and integer schedules |
| Informative clock and domain boundaries | Sections 6--7; parent interval certificate, full serial detector effect, phase generator, ultraviolet cost and occupation/fermionic countercontrols |
| Independent, hostile verification | `check.py` and `verify.py`; strict schema, source and claim hashes, omitted-event and altered-physics mutations, optimized-Python CLI |

The exit is a proved realization within the named source grammar plus a
precise extra transport principle and countermodel. CCG is not derived from
the three axioms; the enormous finite clock remains an analytic witness,
and agreement on vacuum and one-particle states is not an all-Fock result.
