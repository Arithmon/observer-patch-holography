# Exact McKay certificate for the certified SL(2,F5) spin doublet

## Scope

This note documents `sl2f5_mckay_e8_certificate.py`.

The executable starts from the exact `SL(2,F5) -> PORT-SPIN-LIFT` bridge already
certified by `sl2f5_port_spin_bridge_certificate.py`.  It does not assume a
character table and does not identify a graph by group order or an element-order
profile.

## Exact derivation

The certified faithful two-dimensional spin representation is reconstructed
first.  Its exact character over `Q(sqrt(5), i)` is then used to generate
symmetric-power characters through

`S_(n+1) = chi_2 S_n - S_(n-1)`.

Exact finite-group character inner products recover nine pairwise orthonormal
irreducible characters with dimensions

`1, 2, 2, 3, 3, 4, 4, 5, 6`,

whose squared dimensions sum to `120`.  Conjugacy classes are derived directly
from the supplied `SL(2,F5)` multiplication law.  The central element `-I`
splits the recovered irreducibles into five center-trivial and four spinorial
representations.

For every recovered irreducible `V_i`, the certificate computes

`N_ij = <chi_2 chi_i, chi_j>`

exactly.  It checks nonnegative integral multiplicities, dimension conservation,
symmetry, absence of loops, simple-lacedness, connectedness, the tree condition,
and the dimension-vector identity `A d = 2 d`.

Only after this fusion graph has been derived is it compared with a separately
encoded affine `E8` graph.  A dimension-preserving graph isomorphism is required.

## Galois control

The certificate applies the field automorphism

`sqrt(5) -> -sqrt(5)`

while fixing `i`.  The conjugate doublet is checked to remain faithful and is
identified with the second recovered two-dimensional irreducible.  Its McKay
fusion graph is independently recomputed and is again isomorphic to affine
`E8`.

The certified spin character itself is not fixed by this Galois action and has
explicit `sqrt(5)`-sensitive trace values.  Therefore affine-`E8` graph type
alone does not select one of the two real embeddings of `Q(sqrt(5))`.

A reducible two-dimensional control, `1 direct-sum 1`, does not reproduce the
affine-`E8` graph.

## Claim boundary

This certificate proves an exact executable McKay construction for the faithful
`SL(2,F5)` spin doublet already supplied by the conditional `PORT-SPIN-LIFT`
packet.

It does **not** prove or select `phi` rather than its Galois conjugate, state
`27^phi` or any mass relation, physically identify the finite spin action,
source-select the charged-double-triplet current fixture, or formalize McKay or
affine `E8` in Lean.

## Verification

```bash
python3 code/a5_closure/sl2f5_mckay_e8_certificate.py verify
python3 -m pytest -q code/a5_closure/tests/test_sl2f5_mckay_e8_certificate.py
```

The mandatory A5 certificate suite also exercises a focused smoke reconstruction
through `test_coset_carrier_certificate.py`.
