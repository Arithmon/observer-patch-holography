# Yang--Mills collar-gap certificates

## Source-contract verifier

`verify_collar_gap_certificate.py` checks the declared finite row arithmetic used
by the issue-306 theorem witness. Both certificate interfaces use
`collar_gap_contract.py` for exact rational accounting. They require source
descriptions, positive declared rate bounds, refinement targets inside the
declared table, and normalized conditional-probability pairs with exact
total-variation distances below the supplied bounds. Every influence and its
integer multiplicity contributes to the row sum. The reported expression is

```text
gap_lower = c_floor * (1 - eta_upper).
```

Run the bundled contract witness with:

```bash
python3 code/yang_mills/verify_collar_gap_certificate.py \
  code/yang_mills/certificates/issue_306_theorem_contract_witness.json
python3 -m pytest code/yang_mills/test_collar_gap_certificate.py
```

The bundled JSON is marked `theorem_contract_witness`. It tests the arithmetic
and the explicit constant `3/8`. Source descriptions, rate bounds and the
completeness of conditional rows remain supplied premises: comparing one pair
does not establish that it maximizes influence over a source kernel. The
checker does not reconstruct a Gibbs measure or prove the rate, refinement or
continuum hypotheses of the paper's theorem.

Both interfaces reject `scope: physical_source_receipt`, including fixtures with
asserted provenance or continuum fields. Neither implements the source,
transfer, OS/noncollapse or continuum verifications needed for that status.
Accepted outputs have `physical_clay_receipt: false`.

## Finite calibration fixture

`finite_collar_gap_certificate.py` expands and checks a finite Ising
calibration family. Its exact rational table has 244 active types,
`c_floor = 1`, `eta_upper = 1/2`, and `gap_lower = 1/2`. It is deliberately
not a physical compact-simple-gauge Yang--Mills receipt: the physical
placeholder manifest must fail closed for missing source, continuum, and
transfer evidence.

```bash
python3 code/yang_mills/finite_collar_gap_certificate.py certify \
  --manifest code/yang_mills/manifests/atomic_4d_ising_calibration.json \
  --output code/yang_mills/receipts/atomic_4d_ising_calibration.receipt.json
python3 code/yang_mills/finite_collar_gap_certificate.py verify \
  --manifest code/yang_mills/manifests/atomic_4d_ising_calibration.json \
  --receipt code/yang_mills/receipts/atomic_4d_ising_calibration.receipt.json
python3 -m pytest code/yang_mills/tests/test_finite_collar_gap_certificate.py
```

A compact family repeats the complete template with self-targets and independent
copies of its data. Supply either this family or an explicit table. A template
must omit generated IDs and targets; it cannot override them. Boolean or
fractional multiplicities, malformed source fields, and duplicate JSON keys
are rejected. Receipt replay checks every field and its JSON type, including
the full manifest and expanded-table hashes.

The retained audit under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
reproduces a lost-influence defect at `de60560b`: compact bounds `3/4` and `1/4`
were reduced to `1/4`, so a noncontractive table returned a positive gap.
The witness interface also truncated fractional multiplicities and both
interfaces permitted unsupported physical promotion. Independent controls
maximize probability differences over all events, check heterogeneous tables,
and construct a two-spin heat-bath operator with a complete exact eigenbasis.
They verify attained gap bounds at both signs of the correlation, within
`2^-100` of the contractivity boundary, and at rate scales beyond binary64.
Three-spin controls derive every conditional row from the full supplied law
and check the entire weighted operator inequality by exact Schur elimination,
with unequal rates and two influences per site. A singular-support control
constructs a nonconstant centered zero-energy mode and requires rejection at
`eta = 1`; its zero gap is calculated from the conditional expectations.

```bash
python3 -m pytest -q code/yang_mills/tests/test_collar_certificate_accounting.py \
  code/yang_mills/tests/test_finite_collar_gap_certificate.py \
  code/yang_mills/test_collar_gap_certificate.py -W error
```

The 244-type receipt reproduces byte-for-byte; its `1/2` floor and the witness's
`3/8` floor are unchanged. The finite Z2 transfer diagnostic, frozen evidence,
paper and book claims, and registry payloads are unchanged. These corrections
establish no physical mass gap.

## Finite Z2 transfer-receipt diagnostic

`z2_finite_transfer_receipt.py` evaluates the paper's finite
ground-state-transform and cross-fiber receipt on Z2 lattice gauge theory
(L x L periodic spatial torus, Gauss-law sector, one heat-bath collar per
spatial link). The ground-state transform is the Doob transform by the Perron
vector. Two transfer objects are tested: the Wilson transfer matrix with
`H = -log(T / lambda_max)` and the Kogut-Susskind Hamiltonian.

Result on the committed receipt (`receipts/z2_finite_transfer_receipt.json`):

* the receipt is exact at `beta_s = 0` (every rate equals `log coth beta_t`, the dual coupling);
* for every tested interacting Wilson point
  `beta_s = beta_t in {0.1,0.3,0.5,0.7,1}` it fails: the best constant-rate
  fit leaves a 4% to 24% relative residual at `L = 3`;
* for the Kogut-Susskind Hamiltonian the single-flip form is exact with
  fiber-dependent rates `c_l(o) = lam (r + 1/r)`, `r = Omega(o)/Omega(X_l o)`,
  so the scalar cross-fiber equality fails (rate spread 1.04 to 2.05 at
  `L = 3`) while `c_l(o) >= 2 lam` remains an analytic positive floor;
* among the interacting tested points, the Dobrushin sum `eta_*` is below one
  at the weakest Wilson and Kogut-Susskind couplings.

This is a finite-grid diagnostic on a toy gauge system (`physical_clay_receipt:
false`), not a universal no-go or a compact-simple-gauge receipt. Anisotropic
couplings, the anisotropic Wilson-to-Hamiltonian limit, alternate transfer
objects, and fiber-dependent positive rates remain viable. The next proof
target is the quotient-space variable-rate approximate-tensorization lemma.

```bash
python3 code/yang_mills/z2_finite_transfer_receipt.py --L 2 3 \
  --output code/yang_mills/receipts/z2_finite_transfer_receipt.json
python3 -m pytest code/yang_mills/tests/test_z2_finite_transfer_receipt.py
```
