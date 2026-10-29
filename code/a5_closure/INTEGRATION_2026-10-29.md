# Reviewed source custody after integration

The operational attachment probes were reviewed together with the current
native-model changes. Their original work baselines remain in each receipt.
The original four receipts, producers and independent verifiers are preserved
at commit `78c4eb137cf3bf20e7c1d907561f2f2abe22cbf4`; for example
`git show 78c4eb137cf3bf20e7c1d907561f2f2abe22cbf4:code/a5_closure/receipts/proper_code_attachment_probe.json`
recovers the exact earlier receipt. Current source and prerequisite pins are
advanced separately from those historical work baselines.

The complete imported-source comparison found only these changes relative to
the probes' scientific source pins:

- `model/native_clock/configuration.py`: `success_probabilities` now returns
  fresh `Fraction` objects with the same exact numerator and denominator.
  This prevents callers changing the cache through mutable protected slots.
  Endpoint maps, probability values, source controls and preparation laws
  are unchanged. Four binary/full78 and scan/retry cases were checked for
  exact value equality and independence from mutated returned objects.
- `model/common_source/QUANTUM_NATIVE.md`, used only by the coherent-access
  probe: a paragraph clarifies that the old expanded pulse totals are
  historical, identifies the existing saved 9,738 pulses per flip and links
  the already current receipt. No source code, pulse schedule or number was
  changed by that paragraph.

All other pinned scientific inputs were byte identical. In particular the
native branch, compiler, tomography, pointer and coherent routing code were
unchanged. Current scientific inputs are pinned to
`3f650db7a625f0df66495bdaea4332861507f5d4`. Dependent prerequisite pins identify
the subsequent immutable commits containing their reviewed current receipts.
The original Galois drift diagnostic retains its original source label and
diff hash; integration does not erase that older unresolved provenance check.

Each refreshed receipt is compared with its original after removing only
source/prerequisite/implementation custody fields. Every mathematical result,
matrix-unit map, resource count, output distribution, obstruction scope and
verdict must remain exactly equal. Producers and independent verifiers retain
their strict byte comparisons against the new pinned commits. An arbitrary
future source modification is still rejected; this is no hash bypass.

The retained conclusions remain bounded: proper-code export alone does not
identify spatial loads with quantum states; coherent native access does not
provide a source-derived spatial decoder; and the clean-register decoder
obstruction applies to the fixed finite template grammar. No spectral-gap
operator map, physical particle or calibrated mass is added by integration.
