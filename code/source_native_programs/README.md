# Stored computation from pair means

Computed integers can be stored in one bank while another bank stays readable.
Switching the banks lets later layers use earlier results. Transfers, addition,
retirement and reuse use supported scalar pair means with finite precision.
The source has to select the prepared records and temporal controller.

The [contract](CONTRACT.md) defines acceptance and scope; the
[derivation](DERIVATION.md) gives the construction and bounds; the
[family receipt](families_receipt.json) distinguishes complete small native
executions from large-population route certificates.

From the repository root, set `PYTHONPATH=code` (PowerShell:
`$env:PYTHONPATH='code'`), then run:

```text
python -m source_native_programs.build
python -m source_native_programs.verify
python -m pytest -q code/source_native_programs
python -m source_native_programs.family_build --workdir temp/native-program-witnesses --output temp/native-program-families.json
python -m source_native_programs.family_verify --workdir temp/native-program-witnesses --families temp/native-program-families.json
```

Family regeneration takes several minutes and retains its large witnesses in
the work directory. `--reuse-routes` reuses existing witnesses after complete
validation. Without that flag every route is regenerated. Verification compares
the compact committed receipt; it does not replace it. Intentional receipt
regeneration uses `--write-receipt` on the corresponding verifier.

```text
cd Lean
lake build Geometry.SourceNativeProgramsAxiomAudit
```
