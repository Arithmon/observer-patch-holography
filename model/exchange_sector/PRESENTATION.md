# A finite presentation of native exchange

The original native two-complex admits a direct noncommutative certificate:

```
235 generators, 495 native face/square relators
  -- 234 verified elementary eliminations -->  <g | g^2>.
```

The trace is particularly small: 167 eliminations set a generator to the
identity, and 67 replace it by one other generator or its inverse. No longer
replacement word is needed for this native presentation.

Thus its fundamental group is `Z/2`. The surviving generator has odd
ordered-particle cover parity, and the native triangular exchange reduces
to its nontrivial class. This proves the exhaustive positive/sign character
classification on the declared deformation domain using only that finite
two-complex. The [deleted-product proof](STATISTICS.md) remains an independent
geometric explanation; no auxiliary higher cells or continuous surface
classification is needed by this certificate.

## What the verifier proves

Choose the existing 65-edge spanning tree of the connected configuration
graph. Collapsing the tree gives one free generator per remaining edge.
An oriented native polygon gives the ordered word of its non-tree edges,
with inverse letters on reversed edges. This is the usual CW fundamental
group presentation, before abelianization. Every one of the 495 admitted
local polygons is retained initially.

A relator containing generator `g` exactly once can be cyclically rotated
to `g^s w=1`, where `s=+1` or `-1` and `w` contains no `g`. Solve it as
`g=w^-1` or `g=w`, substitute in every other relator, and remove that
generator and its defining relator. This elementary Tietze transformation
preserves the presented group. Free cancellation, inversion, cyclic
conjugation, and deletion of duplicate or trivial relators also preserve
the normal closure. Letters are never commuted.

The producer chooses a deterministic elimination order. The verifier does
not repeat that choice: it checks that each selected relator is currently
present, that the generator occurs exactly once, that substitution solves
the equation by free cancellation, and that the replacement uses only the
other surviving generators. It reconstructs every subsequent presentation.
The final generator inventory includes any unused free generator; such a
generator cannot disappear merely because it occurs in no relator.

The final presentation contains exactly one generator and exactly its square
relator. This establishes the group, rather than inferring it from `H_1`.
The independent ordered-pair lift sends that generator to the nontrivial
element of `Z/2`; the actual exchange witness maps to it as well. Removing
the triangle relations leaves a different presentation and does not pass
the native certificate. A commutator word whose exponent sums vanish is
also retained unless a valid group transformation removes it.

The trace is an exact finite proof artifact generated from current sources.
It is not a Lean kernel proof, a bounded native rewrite algorithm for every
input path, a physical selection of the sign character, or a spin-statistics
theorem. Native operations, retained records and resource costs are unchanged.

```
python -m model.exchange_sector.presentation
python -m pytest -q -W error model/exchange_sector/test_presentation.py
```
