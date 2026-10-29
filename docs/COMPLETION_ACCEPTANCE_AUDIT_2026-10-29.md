# Completion acceptance and the inherited particle-clock obstruction

The finite observer service accounting is stronger than the earliest CP-01
summary suggests: CP-01 permits services that are **constructed or explicitly
counted**. Missing physical hardware alone therefore does not defeat its
acceptance contract. The substantive question is whether its complete
observer, reference, operations and refinement interface belong to the same
frozen construction as the field and particle claims.

This review supplies an exact negative answer for one natural proposed join.
Adding independent copies of the existing native reference cannot extend the
eta two-state clock to the moving pair band while preserving its modular
flow. A positive, size-independent error separates the two generators. This
removes that proposed shortcut from the active integration route and specifies
which premise a larger supplied field-clock construction must change.

The [checker](../model/completion/acceptance_audit.py) binds the reviewed source
and receipts. The [tests](../model/completion/test_acceptance_audit.py) check
the witness against both momentum-space fibers and the independently assembled
8,128-dimensional position-space CAR action. No canonical physical premise or
integrated gate is discharged by this finite exclusion.

## Exact source and field comparison

Let `q=25/8`. The canonical native binary total-two reference has eigenvalues

```text
sigma_x = (64/7505) q^m(x),
m(x) = x_0+x_1 in {0,1,2},
multiplicities = {0:45, 1:20, 2:1}.
```

For `n` independent copies, the eigenvalues are `(64/7505)^n q^m`, with integer
`0 <= m <= 2n`. The exact multiplicities are the coefficients of
`(45+20z+z^2)^n`. Every modular-generator difference is an integer multiple of
`log q`. This is a statement for every finite `n`; the evaluator checks finite
examples without allocating their tensor products. Restricting to a subspace
invariant under the full modular group preserves this spectral property.

The existing [eta clock map](../model/particle_clock/README.md) has field
energies `0,U` and fixes `beta U=log q`. In these same energy units the native
inherited generator, modulo a scalar, has integer-spaced levels. Neither more
copies nor a common change of time and energy units changes the spacing ratio.

Consider the **same** periodic Weyl/contact action at `c/U=1/8`. At total
momentum zero, write `H=T+U|d><d|`, where `d` is the normalized constant singlet.
For every even periodic side `L>=4`, discrete Fourier orthogonality gives

```text
<d,Td> = 0,
<d,T^2 d> = 6 c^2 = 3 U^2/32,
<d,T^3 d> = 0.
```

Indeed, `h(-k)=-h(k)`, the singlet obeys
`(sigma_a tensor I + I tensor sigma_a)d=0`, and each uniform momentum average
of `sin(k_a)^2` is `1/2`. Odd powers cancel under `k -> -k`. The singlet and
`Td` satisfy the fermion exchange constraint; this is a physical two-fermion
subspace. The checker independently evaluates the `L=4` identities using
Gaussian integer arithmetic, with no rounded trigonometry or fitted pole.

The compression to the orthonormal span of `d,Td` is

```text
H_compressed/U = [[1, sqrt(3/32)], [sqrt(3/32), 0]].
```

Its characteristic polynomial at `13/12` is `-1/288`. The variational
principle therefore gives `E_pair/U > 13/12`. The free norm bound
`||T||/U <= sqrt(3)/4 < 1/2` gives `E_pair/U < 3/2`. The rank-one perturbation
has only one eigenvalue above the free band, so the eigenvalue bounded here
is the existing isolated pair pole, not an unrelated excitation. In particular,

```text
13/12 < E_pair/U < 3/2,
distance(E_pair/U, integers) > 1/12.
```

The numerical retained `L=4` value is approximately `1.08859449`; the exclusion
uses only the strict rational bounds.

## What the exclusion quantifies

For any finite number of independent canonical native copies, there is no
isometry from a field domain containing both the eta dark state and the
zero-momentum pair that intertwines their full field evolution with the
inherited native modular evolution at the eta calibration. The contradiction
is the noninteger energy difference above. It occurs on two actual
field eigenstates, so retaining a larger spatial field algebra cannot cure it.

There is also a robust bound. Let `A` have the inherited integer spectrum in
units `U`, let `V` be any isometry, and allow any scalar energy offset `a`.
If

```text
||(A-aI)V - V H/U|| <= epsilon,
```

then applying this estimate to each of the two normalized field eigenstates
puts its eigenvalue within `epsilon` of the shifted integer spectrum. Their
difference is consequently within `2 epsilon` of an integer. Hence
`epsilon > 1/24`. This bound does not assume that `V` is invariant.
Matching `V* A V` to the field Hamiltonian alone is insufficient: the tests
give an explicit compression with the correct expectation and nonzero
off-subspace leakage.

An operational countercontrol uses the equal superposition of the eta dark
state and the zero-momentum pair. At Hamiltonian model time `|t|=2 pi/U`,
the inherited native modular evolution is scalar for every finite copy count.
The field assigns a different phase to its pair component. A two-outcome
measurement onto the initial superposition and its orthogonal complement has
probability difference

```text
sin^2(pi E_pair/U) > 1/36.
```

The half trace distance exceeds `1/6`; appending an untouched reference
preserves the witness. These conservative inequalities use
`sin(pi delta) >= 2 delta` for `delta` between `0` and `1/2`.
The corresponding native modular parameter is `|s|=2 pi/log q`. This is an
analytical comparison of ideal flows, not an executed long-parameter source
word or a physical duration. The existing record-clock experiment's declared
parameter cap is separate. The generator obstruction applies even when the
requested exact intertwining interval is small.

The exclusion does **not** cover a different reference, correlated copies,
different port weights, a new law, nonisometric simulation with additional
services, or a supplied actuator implementing a separate Hamiltonian. In
particular, a supplied global field reference
`rho_field=exp(-beta H)/Tr exp(-beta H)` gives the desired modular relation by
construction. It changes the reference premise. Its preparation, selection,
global normalization and relation to the native process must then be recorded;
the old product-reference preparation proof does not transfer to it. This is
the distinction between inheriting a clock from a source state and using the
source's inferred rate to program a supplied field actuator.

## Capacity is a second, independent condition

The complete two-fermion space on a two-component periodic cube has dimension
`choose(2L^3,2)`. A faithful isometric encoding into `n` native `H66` factors
requires `66^n >= choose(2L^3,2)`.

| Periodic side | Complete two-fermion dimension | Native H66 factors required |
| --- | ---: | ---: |
| 4 | 8,128 | at least 3 |
| 6 | 93,096 | at least 3 |
| 8 | 523,776 | at least 4 |
| 12 | 5,970,240 | at least 4 |
| 16 | 33,550,336 | at least 5 |

These are necessary capacities for the complete stated algebra and arbitrary
reference-preserving states. They are not encoder constructions or memory
requirements for every restricted observable menu. The compiled 512-address
work register is an explicit coherent service; its unused addresses are not
extra native stationary `H66` states. Classical exported records also do not
provide coherent Hilbert-space capacity. Satisfying this dimension bound does
not solve the independent spectral obstruction.

## Acceptance review and the smallest missing joins

The contracts are [CP-00 #1135](https://github.com/FloatingPragma/observer-patch-holography/issues/1135),
[CP-01 #1136](https://github.com/FloatingPragma/observer-patch-holography/issues/1136)
and [CP-13 #1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137).
Their dependency ordering is part of the contract: CP-01 requires CP-00;
CP-13 requires both. A completed theorem map can keep future physics premises
conditional, but it must map to an actual accepted predecessor construction.

| Criterion | Evidence and remaining concrete work |
| --- | --- |
| CP-00: law and permitted implementation class | Fixed C/E and bounded endpoint-biased finite classes are explicit. Their contracts expressly exclude a full carrier/refinement freeze. A common class must state whether independent source products, correlated references, field actions and coherent actuators are primitives or derived constructions, and enumerate its allowed variation. |
| CP-00: targets, calibration, exposure and stopping | The common hypothesis declares the charged spectrum, neutrino and gravity targets. Existing public-data comparisons retain calibration ancestry and failures. The proposed V2 contract binds a common target/observable menu, uncertainty rule, allowed calibrations, computation budget and decision thresholds for independent review. It cannot relabel earlier comparisons prospective. |
| CP-00: actual source and refinement interface | The finite observer and field services have concrete interfaces. A tensor-product source declaration alone fails the inherited clock test proved here. A replacement must specify its actual state/reference and embedding or simulation maps, including how source and field records compose under growth. |
| CP-01: continuing complete observer | The retained H66 process prepares, records, exports and continues; coherent probes, readout, error bounds and native acquisition are explicit. The record-clock receipt counts 31,458,883 native attempts and the supplied coherent realization separately. These counted services are legitimate evidence under the issue's wording. |
| CP-01: common process, resources and premises | The continuing observer must be an actual member of the adopted source class, with each finite theorem premise proved on its stated interface. CP-01 does not require CP-06's full spatial field clock. The eta M2 map and its obstructed product extension constrain that later join. The theorem ledger must distinguish proved observer instances from the undischarged global source-selection obligations. |
| CP-01: original q=3 recipe | Its original common-code preparation/history stays unexecuted. The active program must either execute that recipe or explicitly retire its adoption through the full CP-00/CP-13 contract. A new benchmark cannot silently perform this retirement. |
| CP-13: statements, maps and cleanup | The finite attachment review assigns all 84 existing premise rows and records the actual domains of 13 reviewed attachments. This review additionally excludes the inherited product-clock extension. Acceptance requires the full declared class and accepted observer interface; it does not require every later field or gravity theorem premise to be solved in advance. |

A CP-00 freeze need not prove that its class meets the physical targets. It
must, however, define the actual admissible class and the test that will decide
that question. A union of existing benchmark names does not supply composition
rules or a common refinement. Conversely, the ability of a supplied universal
actuator to encode an arbitrary desired action does not establish the native
selection required by the common hypothesis. A constructive next contract
should bind a specific service grammar and growth family, preserve this
selection obligation, and classify all earlier fits and exclusions before
future validation.

The [V2 class and observer proposal](EXECUTABLE_PROGRAM_CLASS_V2_PROPOSAL.md)
addresses the predecessor acceptance route explicitly, including a typed
source grammar and independent review. The inherited-product clock join is
excluded on its stated domain regardless of the administrative gate decision.
The finite component campaign itself is not substituted for the complete
program.

Run the source-custody and exact checks from the repository root:

```sh
python3 -m model.completion.acceptance_audit --check
python3 -m pytest -q -W error model/completion/test_acceptance_audit.py
```
