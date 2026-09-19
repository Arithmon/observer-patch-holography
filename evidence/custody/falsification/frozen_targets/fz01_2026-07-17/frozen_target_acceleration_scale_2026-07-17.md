# DRAFT frozen target: galaxy acceleration scale (unanchored)

Status: DRAFT. Not registered, not externally timestamped, carries no verdict
weight. Anchoring (OpenTimestamps or equivalent over the adopting commit) and
adoption into `reverse-engineering-reality/falsification/frozen_targets/` are
owner actions. Until anchored, this file is a proposal, and any comparison
against data published before the anchor date is retrospective.

This target is absent from tracker section 8; it is registered as a standalone
draft alongside the FZ-01 registry
(`fz01_freeze_registry_2026-07-17.json`, same directory).

## Target statement

Two readings of the Z6-reserve effective acceleration scale are on record in
`reverse-engineering-reality/cosmology/oph_dark_matter_paper.tex`. Both are
printed here; the paper designates the presence reading as current and labels
the Poisson reading superseded for its stored diagnostics. Menu size M = 2 is
declared and charged.

```
presence reading:  lambda_c = 1 - P/24 = 0.932042991274835
                   a_eff = 1.184737388e-10 m/s^2
Poisson reading:   lambda_c = e^(-P/24) = 0.934300639489386
                   a_eff = 1.179018696e-10 m/s^2
```

Source lines, verified against the paper text:

- presence: `\lambdac=1-\frac{P}{24}=0.932042991274835\ldots` (line 1218,
  eq. z6lambda) and `\aeff=1.184737388\times10^{-10}` (lines 1223 and 2214).
- Poisson: "the superseded mean-count/Poisson scale
  \(\aeff=1.179018696\times10^{-10}\,\mathrm{m\,s^{-2}}\)
  (\(\lambdac=e^{-P/24}\))" (lines 2243-2244).
- separation: "The two scales differ by \(0.485\%\)" (line 2246).

Pixel value: P = 1.630968209403959. Context value, no verdict attached: the
capacity benchmark `a_OPH = 1.029186271e-10 m/s^2` (line 2197) carries the
measured H0 and Lambda through SL-4/SL-5, and the paper classifies every
comparison of it against the empirical scale as a consistency display
(lines 2199-2205); it is recorded here so it cannot later be promoted to a
prediction.

## Comparator (current, informational and retrospective)

Common empirical reference: `a_0 ~ 1.2e-10 m/s^2` (McGaugh-Lelli-Schombert
RAR fit to SPARC; paper lines 2206-2209). SPARC is not held out: the
interpolating function is fitted to SPARC itself, so any current agreement is
retrospective and scores nothing. Current offsets: presence -1.272 percent,
Poisson -1.75 percent, against the 1.2e-10 reference.

## Decision policy (freezes with the target)

- Decision data: a sub-percent determination of a_0 from a held-out dataset.
  Held-out means: a galaxy sample disjoint from SPARC (or a probe class
  independent of rotation curves), with selection function, error model, and
  interpolating-function convention frozen before scoring. SPARC and any
  SPARC-derived recalibration are ineligible as scorers.
- Verdicts per reading, z = |a_0_post - a_eff_reading| / sigma_post:
  - z <= 2: COMPATIBLE for that reading.
  - z > 3: FAIL for that reading.
  - 2 < z <= 3: INCONCLUSIVE, carried to the next qualifying determination.
- Branch verdict: the Z6-reserve acceleration branch FAILs as formulated when
  both readings FAIL under the same scorer. A COMPATIBLE for exactly one
  reading resolves the presence/Poisson selection empirically.
- Internal discriminator: the readings differ by 0.485 percent. A held-out
  determination separates them at 2 sigma when sigma(a_0) <= 0.24 percent and
  at 3 sigma when sigma(a_0) <= 0.16 percent. A determination at ~0.5 percent
  precision begins to weigh the presence/Poisson selection; coarser
  determinations test only the common scale.
- Declared exposure: a held-out result centered at 1.2e-10 with
  sigma = 0.5 percent places the presence reading at z ~ 2.5 (INCONCLUSIVE,
  adverse) and the Poisson reading at z ~ 3.5 (FAIL). The paper records that
  an exact hit on 1.2e-10 would require about 9 percent more protected
  inactive reserve than the P/24 behind lambda_c (lines 2237-2240); that
  repair path is charged as description length and requires a new
  registration, it cannot rescue this target in place.
- No verdict is claimed at current precision; the unanchored draft carries no
  verdict weight in any direction.

## Integrity

sha256 of the tracker section-8 text this registry cycle pins:
`54b216b52afe25895552a04b450c3cb33e303802e99618f24361b2bf2b9e0737`
(bytes 73283-75997 of
`oph-physics-sim/docs/OPH_SIGNATURE_EXPERIMENT_TRACKER.md`, pinned copy
`tracker_section8_pinned_text.md`). sha256 of the fenced target block above is
computed at anchoring time over the adopted file.
