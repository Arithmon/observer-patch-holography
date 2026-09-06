---
status: work-plan
claim_level: planning
physical_claim: false
branch: arithmon/golden-lift-2-port-group
base: c711eef134ea4290519759bc1524bad9f9004f75
---

# GOLDEN-LIFT-2A — type the twelve-port rotation group

## Goal

Close the explicit group-theoretic gap left by GOLDEN-LIFT-0/1.

Current exact chain:

- `A5PortSixAxesBridge` gives a pointwise `60 × 6` identification between the antipodal quotient of the sixty committed twelve-port rotations and `A5SixAxes.L60`, under the fixed axis relabeling and `rowEquiv`.
- `PSL2F5SixAxesBridge` gives a typed group isomorphism
  `PSL2F5 ≃* SixAxisGroup`, where `SixAxisGroup` is exactly the subgroup carried by `A5SixAxes.L60`.

Missing edge:

```text
PortGroup ≃* SixAxisGroup
```

where `PortGroup` is the actual subgroup of `Equiv.Perm (Fin 12)` carried by the sixty committed port rotations.

Headline target:

```text
PSL2F5 ≃* PortGroup
```

This is a purely finite algebraic bridge. It does not identify the group with abstract `A5`, construct `2I`, invoke McKay, transport the golden sectors as typed representations, select `φ`, or make a physical claim.

## Fresh base

This branch was created from the fork's synchronized `main` at

```text
c711eef134ea4290519759bc1524bad9f9004f75
```

which is byte-identical at the ref level to `FloatingPragma/main` at branch creation time.

## Recommended module

```text
Lean/Screen/A5PortGroupBridge.lean
```

Suggested namespace:

```lean
namespace OPH.A5PortGroupBridge
```

Imports:

```lean
import A5PortSixAxesBridge
import PSL2F5SixAxesBridge
```

Do not modify the older bridges merely to make the new theorem easier unless a tiny reusable lemma is genuinely cleaner there.

## 1. Reuse the existing row data rather than reclassifying the group

Use:

```lean
OPH.A5PortSixAxesBridge.portEl   : Fin 60 → Equiv.Perm (Fin 12)
OPH.A5PortSixAxesBridge.rowEquiv : Fin 60 ≃ Fin 60
OPH.A5SixAxes.mulT               : Fin 60 → Fin 60 → Fin 60
OPH.A5SixAxes.invT               : Fin 60 → Fin 60
```

Define multiplication and inverse indices for the port rows by transport through the already certified row equivalence:

```lean
def portMulIndex (i j : Fin 60) : Fin 60 :=
  rowEquiv.symm (OPH.A5SixAxes.mulT (rowEquiv i) (rowEquiv j))

def portInvIndex (i : Fin 60) : Fin 60 :=
  rowEquiv.symm (OPH.A5SixAxes.invT (rowEquiv i))
```

The load-bearing new check is that these transported indices reproduce multiplication and inverse of the *actual twelve-port permutations*.

Targets:

```lean
theorem portEl_zero : portEl 0 = 1

theorem portEl_mul (i j : Fin 60) :
  portEl i * portEl j = portEl (portMulIndex i j)

theorem portEl_inv (i : Fin 60) :
  (portEl i)⁻¹ = portEl (portInvIndex i)
```

Preferred proof style: pointwise equality of permutations followed by kernel-checked finite tables. Avoid asking `decide` to compare `Equiv.Perm` proof structures directly at quadratic scale.

For example, prove `portEl_mul` with `Equiv.ext` and a closed finite check over `(i,j,k) : Fin 60 × Fin 60 × Fin 12`. If one monolithic `decide` is too expensive, copy `A5SixAxes.rowF_mul`'s banded strategy (5 rows of `i` per private lemma, then dispatch by bounds). This is acceptable and preferable to a cardinality argument.

## 2. Define the actual port subgroup as the range of the committed rows

Do not make membership depend on an arbitrary generated closure if it can be avoided. Once the three row-law lemmas above are available, define:

```lean
def PortGroup : Subgroup (Equiv.Perm (Fin 12)) where
  carrier := Set.range portEl
  one_mem' := ⟨0, portEl_zero.symm⟩
  mul_mem' := by
    rintro _ _ ⟨i, rfl⟩ ⟨j, rfl⟩
    exact ⟨portMulIndex i j, portEl_mul i j⟩
  inv_mem' := by
    rintro _ ⟨i, rfl⟩
    exact ⟨portInvIndex i, portEl_inv i⟩
```

Adjust equality directions as Lean requires.

This definition makes the scientific object exact: `PortGroup` contains exactly the sixty committed permutations, not a larger generated subgroup inferred afterward.

## 3. Prove row indexing is faithful

Reuse the existing antipodal quotient faithfulness rather than a fresh brute-force distinctness theorem if convenient:

```lean
theorem portEl_injective : Function.Injective portEl
```

Route:

- if `portEl i = portEl j`, then `quotientAxis i = quotientAxis j` pointwise;
- apply `A5PortSixAxesBridge.quotient_action_faithful`.

Then define a canonical row index for a subgroup element using its range witness:

```lean
noncomputable def portIndex (g : PortGroup) : Fin 60 :=
  Classical.choose g.property

theorem portIndex_spec (g : PortGroup) :
  portEl (portIndex g) = g.1 := Classical.choose_spec g.property
```

and derive the expected uniqueness lemma for row elements.

No scientific claim depends on the choice: injectivity proves the witness unique.

## 4. Build the homomorphism to the typed six-axis subgroup

`PSL2F5SixAxesBridge.SixAxisGroup` is already the exact subgroup carried by `A5SixAxes.L60`.

Define:

```lean
noncomputable def portToSix :
    PortGroup →* OPH.PSL2F5SixAxesBridge.SixAxisGroup := ...
```

Underlying map:

```lean
g ↦ ⟨OPH.A5SixAxes.el (rowEquiv (portIndex g)), OPH.A5SixAxes.el_mem _⟩
```

For the homomorphism law, first expose the row-level six-axis multiplication lemma if needed:

```lean
theorem sixEl_mul (i j : Fin 60) :
  OPH.A5SixAxes.el i * OPH.A5SixAxes.el j =
    OPH.A5SixAxes.el (OPH.A5SixAxes.mulT i j)
```

This should follow structurally from `A5SixAxes.el_apply` and `A5SixAxes.rowF_mul`; do not recompute the 60×60 table.

The port index of a product should be forced by `portEl_mul` and `portEl_injective`.

## 5. Prove bijectivity by explicit row realization

Targets:

```lean
theorem portToSix_injective : Function.Injective portToSix

theorem portToSix_surjective : Function.Surjective portToSix
```

Injectivity:

- equality in `SixAxisGroup` gives equality of the corresponding `A5SixAxes.el` rows;
- use `A5PortSixAxesBridge.six_axis_rows_injective` / `A5SixAxes.el_apply` as appropriate;
- `rowEquiv` and `portEl` are injective.

Surjectivity:

- for `g : SixAxisGroup`, use `A5SixAxes.mem_iff_el g.property` to obtain its unique row `j`;
- take the port row `i = rowEquiv.symm j`;
- package `portEl i` into `PortGroup` by its range witness;
- show it maps to `g`.

Then define:

```lean
noncomputable def portGroupEquivSixAxisGroup :
    PortGroup ≃* OPH.PSL2F5SixAxesBridge.SixAxisGroup :=
  MulEquiv.ofBijective portToSix
    ⟨portToSix_injective, portToSix_surjective⟩
```

## 6. Compose with GOLDEN-LIFT-1

Final exact interface:

```lean
noncomputable def pslEquivPortGroup :
    OPH.PSL2F5SixAxesBridge.PSL2F5 ≃* PortGroup :=
  OPH.PSL2F5SixAxesBridge.psl_equiv_six_axis_group.trans
    portGroupEquivSixAxisGroup.symm
```

Check the orientation of `.trans` / `.symm` in Lean and use the shortest typed expression that compiles.

The final theorem should mean exactly:

```text
Mathlib PSL(2,F5)
    ≃* committed six-axis subgroup
    ≃* committed twelve-port rotation subgroup.
```

It is not an abstract classification theorem.

## 7. Useful corollaries, only after the equivalence compiles

Good optional facts:

```lean
theorem portGroup_card : Fintype.card PortGroup = 60

theorem every_port_group_element_is_committed_row (g : PortGroup) :
  ∃ i : Fin 60, g.1 = portEl i
```

The latter is definitionally true from `Set.range`; keep it concise.

Cardinality is a corollary, never the proof of the isomorphism.

Do not use `|PortGroup| = |SixAxisGroup| = 60` as the load-bearing identification.

## 8. Deliberate scope boundary

This revision MUST NOT claim or implement:

- abstract `PSL(2,5) ≅ A5` classification;
- `SL(2,5) ≅ 2I` / binary-icosahedral identification;
- an `SU(2)` or quaternion realization;
- McKay correspondence or the affine `E8` graph;
- typed golden `3` / `3'` representations of `PSL2F5` (that is the next revision);
- derivation or selection of `φ`;
- `27^φ`, Koide, or any mass law;
- physical rotations, physical Spin, or an observable.

The only promotion is from a pointwise row/action bridge to a typed isomorphism of the two already committed finite groups.

## 9. Suggested two-revision workflow

### Revision 1 — proof surface

Implement only `Lean/Screen/A5PortGroupBridge.lean` and the minimum umbrella/lake registration required to build it.

Required checks before any scientific documentation promotion:

```bash
cd Lean
lake build A5PortGroupBridge
lake build OPHScreen
```

Run `#print axioms` on at least:

- `portEl_mul`
- `portEl_inv`
- `portToSix_injective`
- `portToSix_surjective`
- the final group equivalence or its supporting theorem if the definition itself is not printable in the same way.

No `sorry`, `admit`, `axiom`, or `native_decide`.

Suggested commit:

```text
Add typed twelve-port rotation group bridge
```

### Revision 2 — promotion / bookkeeping

Only after Revision 1 builds:

- import the module from `OPHScreen.lean`;
- register it in `Lean/lakefile.lean` if not already done in Rev.1;
- update `Lean/docs/PROOF_INDEX.md`;
- update the golden-sector/projective-cover claim wording only to the exact new level;
- update novelty/falsification matrices consistently;
- update the owning paper paragraph, replacing the old caveat that the port bridge is only pointwise;
- regenerate the active-surface inventory;
- refresh theorem-count floors if the threshold crosses;
- run structural gates and the mandatory suite.

Suggested commit:

```text
Document typed port-group identification
```

## 10. Falsifiers / no-cheating checks

GL2A fails if any of the following occurs:

- `portEl_mul` disagrees for one committed row pair and port;
- `portEl_inv` disagrees for one committed row and port;
- the transported row map fails the monoid law;
- `portToSix` has a kernel or misses a committed six-axis row;
- the proof uses only equal cardinalities to infer identification;
- the new claim is promoted to abstract `A5`, `2I`, McKay, golden-representation transport, or physics.

If `rowEquiv` does not respect multiplication, that is a scientifically useful negative result: keep the existing pointwise bridge and stop rather than altering `rowEquiv` post hoc to force the theorem without documenting the change.

## 11. Next step, explicitly out of this PR

If GL2A closes, GOLDEN-LIFT-2B can transport the already certified golden `3` and `3'` sectors to typed `PSL2F5` representations and prove their pullbacks to `SL2F5` are center-trivial (`-I` acts as identity). That separation should be established before any binary-spinor or McKay step.
