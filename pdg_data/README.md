# PDG legacy snapshot

`particle_masses.json` and `particle_masses.csv` are a legacy export of particle
masses made through the PDG Python package (`tools/fetch_pdg_data.py`). They are
kept byte for byte and hash-pinned in
`code/audit/external_data_provenance_registry.json`, which records them as legacy
duplicates. No claim surface reads them; each comparison value in the papers and
ledgers carries its own provenance row.

Two properties follow from how the snapshot was exported. Several values (tau,
Z, top, Higgs) carry the PDG API's floating-point serialization, with more
digits than the published numerals. The electron row is the CODATA 2022 value
rather than the PDG average, and the files do not label it; the PDG 2026
listing carries the same numerals, 0.51099895069(16) MeV. The upstream PDG
edition was not recorded at export time.
