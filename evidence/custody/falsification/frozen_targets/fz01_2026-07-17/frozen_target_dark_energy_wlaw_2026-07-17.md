# DRAFT frozen target: DK-01 fixed-capacity w-law vs DESI DR3 (unanchored)

Status: DRAFT. Not registered, not externally timestamped, carries no
verdict weight. Anchoring (OpenTimestamps or equivalent) and adoption into
`reverse-engineering-reality/falsification/frozen_targets/` are owner
actions. Until anchored, this file is a proposal.

Theorem on file: `proof/epic_wins/dk01_wlaw/DK01_FIXED_N_WLAW.md`
(fixed-N w-law, premises P1 through P4, drift map
`w(a) = -1 + (1/3) d ln N / d ln a`).

## Target statement

Conditional on premises P1 through P4 of the theorem file (fixed total
capacity N over the scored redshift range, the de Sitter capacity relation
`Lambda l_P^2 = 3 pi / N`, constant `l_P`, no repair-rate drift and a closed
dark sector), the OPH fixed-capacity branch predicts, with zero continuous
freedom,

```
(w0, wa) = (-1, 0)    exactly, in the CPL parameterization
```

over the full range any DESI-class BAO + CMB + SNe combination scores. The
prediction is parameter-free: the capacity value `N_star`, the open closure
equation `N = F(N)`, and the 6.6 percent capacity-vs-display gap all cancel
out of `w(a)`.

Companion negative on file (executed 2026-07-14, re-executed 2026-07-17):
the capacity gap projected through a biased-H0 `w0waCDM` fit lands in the
freezing direction, `(w0, wa) = (-1.07, +0.33)` with the H0 display prior,
opposite to the DESI thawing direction in both coordinates and weak
(`Delta chi2 = 1.5`). No display bias can move a constant-Lambda truth into
the thawing quadrant. Artifact:
`epic-wins/analysis/dk01_capacity_gap_w0wa_mock.py`.

## Comparator record (state of the anomaly at draft time)

Values as the tracker records them
(`oph-physics-sim/docs/OPH_SIGNATURE_EXPERIMENT_TRACKER.md`, DK-01):

- DESI DR2 BAO + CMB + SNe prefer a thawing track over LambdaCDM at 2.8 to
  4.2 sigma depending on the SNe compilation (DESY5, Union3, Pantheon+).
- Recorded preference direction: `(w0, wa) ~ (-0.72, -1.0)`.
- Data named by the tracker: DESI DR2 BAO likelihoods, Planck PR4,
  DESY5/Pantheon+ SNe.

Per-compilation numeric cells (owner action; the tracker does not record
these to publication precision):

| Combination | published (w0, wa) | significance vs LCDM |
|---|---|---|
| DESI DR2 BAO + CMB + DESY5 | OWNER SLOT | OWNER SLOT |
| DESI DR2 BAO + CMB + Union3 | OWNER SLOT | OWNER SLOT |
| DESI DR2 BAO + CMB + Pantheon+ | OWNER SLOT | OWNER SLOT |

## Decision data

DESI DR3 BAO likelihoods, combined with a CMB likelihood and one SNe
compilation named at anchoring time.

## Combination freeze requirement

The scoring combination (BAO release, CMB likelihood, SNe compilation, and
the collaboration chains or the re-analysis pipeline used to read the
contour) must be named in the anchored version of this file, before the DR3
cosmology release. Verdicts are scored against exactly that frozen
combination. A change of combination after anchoring voids the verdict and
requires a new registration. If the frozen combination is unavailable at
decision time, the fallback is the DESI collaboration's own headline
DR3 + CMB + SNe `w0waCDM` result, and that fallback must also be named at
anchoring.

## Decision policy

All distances below are measured on the frozen combination's `(w0, wa)`
posterior. "Excluded at k sigma" means the point `(-1, 0)` lies outside the
two-dimensional credible region whose probability mass corresponds to k
Gaussian sigma.

- **Kill (lane).** `(-1, 0)` excluded at more than 3 sigma AND the posterior
  mean lies in the thawing quadrant (`w0 > -1`, `wa < 0`): the fixed-N lane
  fails this registration. Per the corpus kill-condition protocol, the
  audit/re-audit/repair ladder runs before the failure is banked as a
  framework kill; the only in-corpus repairs are the declared premise exits
  (a horizon-growing N law counts as a prediction only if stated and
  anchored before the DR3 release; drift in `l_P` or a second dark
  component are separate lanes with independent constraints).
- **Kill (branch, tracker bar).** Evolving w confirmed at more than 5 sigma
  on the frozen combination: the fixed-capacity branch is falsified
  outright, matching the kill line the tracker carries for DK-01.
- **Forward pass.** `(-1, 0)` inside the 1-sigma region of the frozen
  combination: the fixed-N lane scores a pass, and the DR2 thawing
  preference is banked as having faded, as the fixed-N branch requires.
- **No verdict.** `(-1, 0)` between 1 and 3 sigma: exposure is carried
  forward to the next release; no pass, no kill.
- **Sign check rider.** If the DR3 posterior mean lies in the freezing
  direction (`w0 < -1`, `wa > 0`) at any significance, the result is scored
  against the display-bias mechanism of the mock before any physical
  reading; the mock's predicted bias direction under a wrong H0 prior is
  freezing, and this rider prevents a bias artifact from being banked as
  physics in either direction.

## What this draft is

This draft states the prediction, the comparator, the decision data, and
the thresholds. It is unanchored: it carries no verdict weight, it
timestamps nothing, and it binds no one. OTS anchoring, the combination
freeze, the per-compilation comparator cells, and adoption into the
falsification register are owner actions.

## Integrity

sha256 of the target statement block (exact bytes between the first and
second code fence above) is to be computed at anchoring time over the
adopted file. The prediction is reproducible with no numeric input: it is
the constant `(-1, 0)` forced by the theorem file from premises P1 through
P4.
