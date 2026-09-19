# DRAFT frozen target: scalar tilt (unanchored)

Status: DRAFT. Not registered, not externally timestamped, carries no verdict
weight. Anchoring (OpenTimestamps or equivalent over the adopting commit) and
adoption into `reverse-engineering-reality/falsification/frozen_targets/` are
owner actions. Until anchored, this file is a proposal, and any comparison
against data published before the anchor date is retrospective.

Registry row: FZ01-R01 in `fz01_freeze_registry_2026-07-17.json` (same
directory). This draft freezes together with that row.

## Target statement

Point prediction, evaluated at the pixel value P = 1.630968209403959:

```
n_s = 1 - P/48 = 0.96602149563741752083...
```

(P/48 = 0.03397850436258247916...; the printed digits are exact for the
stated P; further digits require the certified P enclosure.)

## Declared candidate menu (trials counted)

The corpus circulates two tilt formula candidates. Both are printed here;
the selection between them is open and uncounted, and any future significance
statement must charge the menu size M = 2.

- Candidate A: `n_s = 1 - P_star/48 = 0.9660214956` (this target's point
  value; Planck 2018 pull +0.267 sigma).
- Candidate B: `n_s = 1 - e * alpha * sqrt(pi)`, cited exactly from
  `reverse-engineering-reality/cosmology/oph_cosmology_finite_source_cmb_program.tex`
  lines 1078-1081: "Two candidate tilt formulas are on record,
  \(n_s=1-P_\star/48\) and \(n_s=1-\mathrm e\,\alpha\sqrt{\pi}\); the
  selection between them is not counted in any significance statement."
  Evaluated at alpha = 1/137.035999177 this gives
  `n_s = 0.96484114303` (Planck 2018 pull -0.014 sigma). The same source
  (`oph_inflation_without_inflaton_observer_screen_synchronization.tex`
  lines 169-172) records the equivalent form `theta = e (P_star - phi)` as a
  distinct diagnostic hypothesis for the tilt exponent.
- Recorded non-member: the finite-repair-clock diagnostic emits the
  uncertified value `n_s = 0.97893075698` with diagonal
  Delta chi2 ~ +39.26 against the same 83-bin PR3 TT baseline
  (`oph_cosmology_finite_source_cmb_program.tex` lines 1099-1103). That is a
  negative result for the clock/selector diagnostic and it must never be
  silently replaced by the favorable analytic branch after inspecting data.
  It is pinned here so the swap is impossible.

Candidate separation: A - B = 0.0011804.

Canonical status of candidate A: analytic candidate, kappa_rep contraction
certificate pending; the +0.267 sigma Planck comparison locates a basin and
scores nothing.

## Comparator (current, informational only)

Planck 2018 scalar-tilt summary: `n_s = 0.9649 +/- 0.0042`. This value is
recorded for reference; it predates any anchoring and carries no verdict.

## Decision policy (freezes with the target)

- Decision data: the next combined-likelihood tilt release from ACT DR6 /
  SPT-3G / CMB-S4.
- Combination freeze: the scoring object is the marginalized baseline-LCDM
  `n_s` posterior from a single named likelihood combination. Frozen order of
  scorers: (1) ACT DR6 primary CMB combined with Planck (the released joint
  primary-CMB chain of that collaboration pairing), (2) SPT-3G main survey
  combined with Planck as an independent second pass, (3) CMB-S4 baseline as
  the terminal pass. Adding BAO, SNe, or any external dataset to the scorer,
  or changing the combination, requires a new registration before the data
  opening.
- Verdicts per candidate, z = |n_s_post - n_s_candidate| / sigma_post:
  - z <= 2: COMPATIBLE for that candidate.
  - z > 3: FAIL for that candidate.
  - 2 < z <= 3: INCONCLUSIVE, carried to the next scorer in the frozen order.
- Lane verdict: the tilt lane FAILs when both candidates FAIL under the same
  scorer. A COMPATIBLE for exactly one candidate resolves the menu selection
  empirically and retires the other candidate in place.
- Internal discrimination: the two candidates separate at 2 sigma when
  sigma_post <= 0.00059 and at 3 sigma when sigma_post <= 0.00039. Above
  that precision the menu stays open even under a COMPATIBLE verdict.
- No verdict is claimed at current precision; the unanchored draft carries no
  verdict weight in any direction.

## Integrity

sha256 of the tracker section-8 text this target descends from:
`54b216b52afe25895552a04b450c3cb33e303802e99618f24361b2bf2b9e0737`
(bytes 73283-75997 of
`oph-physics-sim/docs/OPH_SIGNATURE_EXPERIMENT_TRACKER.md`, pinned copy
`tracker_section8_pinned_text.md`). sha256 of the fenced target block above is
computed at anchoring time over the adopted file.
