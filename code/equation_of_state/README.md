# Finite field pressure and native work identification

These packages inspect observer-like bounded software patches with local
state, ports, readback, records and repair or feedback operations. They retain
the source histories and independently check the readouts derived from them.
The public evidence distinguishes three interfaces:

| Package | Readout | Required supplies |
| --- | --- | --- |
| [native](native/README.md) | Exact scalar pair-mean traces and witnesses that a fixed-volume history does not select a work response | Canonical scalar mean branch and declared schedule; candidate work laws are separate inputs |
| [scalar](scalar/README.md) | Instantaneous stress, continuous-action pressure enclosures, and a finite-mode Bose thermal EoS | Scalar action, geometry, canonical variables, dilation, preparation, model clock and thermal state |
| [maxwell](maxwell/README.md) | Maxwell stress, displacement-current execution and finite-mode classical Gibbs pressure | Maxwell action, spatial geometry, constitutive coefficients, canonical phase space and ensemble |

`w = p/rho` is undefined when energy density vanishes. A field's mean stress
does not establish isotropy, equilibration, a bath, or the energy and pressure
of the entire observer apparatus. No package identifies native scalar repair
with a momentum-conserving collision gas, selects a physical Hamiltonian or
volume from that repair, or calibrates its model units against nature.

The exact work-extension and field-scaling lemmas are in
[`Lean/Thermodynamics/EquationOfState.lean`](../../Lean/Thermodynamics/EquationOfState.lean).
The detailed scientific account is in
[`FINITE_EQUATION_OF_STATE.tex`](../../paper/tex_fragments/FINITE_EQUATION_OF_STATE.tex).
Package READMEs give the executable interfaces and interpretation boundaries.
JSON receipts retain complete declared case sets, source hashes, raw values
and controls. CSV tables are inspection aids; exact arithmetic and enclosure
claims should be checked against their JSON parents and independent verifiers.

From the scientific repository root, run all adversarial checks with:

```bash
python3 -m pytest -q code/equation_of_state
cd Lean
lake build EquationOfState
```

Python dependencies are NumPy, SciPy, SymPy, mpmath and pytest. Individual
receipts or the export manifest record the versions used. Exact-native
verification needs only Python's standard library. Reproducing the native
producer additionally needs the sibling `oph-physics-sim` checkout; its
source bytes are pinned in the native receipt.
