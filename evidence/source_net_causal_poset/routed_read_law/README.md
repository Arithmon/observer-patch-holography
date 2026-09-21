# Complete routed record-read histories

The four q13/q21 baseline and centre-source +1 intervention receipts retain
their input, resource and complete decoded event-stream commitments. Each
phase manifest retains its decoded checksum and the byte sizes/checksums of
the original compressed parts. The production `.input` and `.tape` payloads
are regenerated on demand outside Git. The q3 controls remain local fixtures.

`regenerate.py` reconstructs the exact inputs, executes every native producer
event and feeds the complete stream into independent native replay. It checks
the original phase and whole-stream hashes, full semantics, all metric
decisions and every completed logical value. It does not rewrite the receipts
or substitute summary checks for full replay. Its checks authenticate decoded
events; they do not attest durable custody of original compressed files.

The [routing package](../../../code/source_read_routing/README.md) gives the
compiler, conditional theorem, resource model and reproduction commands.
The parent `archive_manifest.json` inventories only retained repository files.
The segment part lists describe optional materializations, not local files
or an external archive location. The supplied read/feedback law, physical
clock and quantum realization remain separate scientific obligations.
