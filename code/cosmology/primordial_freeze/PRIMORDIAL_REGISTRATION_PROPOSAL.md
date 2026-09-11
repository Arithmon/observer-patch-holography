# FZ-15 to FZ-17 registration proposal: primordial branch, pre-SO/CMB-S4/LiteBIRD

Status: PROPOSAL. Nothing in this document is a freeze, a frozen prediction,
or a score. Adoption of the cells below as register rows, anchoring, the
numeric kill-band cells, the `content_sha256`, the `frozen_utc` stamp, and the
custody and attestation records are owner actions under the custody
discipline of the DK-01 draft target. Adding rows to
`claims/frozen_prediction_register.json` triggers the register pin cascade
(forecast contract, discriminator strata, angular certificate, closure
preflight) and is likewise an owner action. Until the owner freezes and
anchors a row, this document binds no one.

## Sources of record

- Edge-center tilt theorem (conditional): `paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex`,
  theorem "edge-center tilt and repair-clock reconciliation":
  `theta = P_star/48`, `n_s = 1 - P_star/48`, under the reserve-generator
  receipt (full-collar density `P_star/24`) and the orientation-reversal
  coarea identity. The finite source generator receipt is not supplied.
- Rank-one single-clock branch (conditional): `paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex`,
  definition "rank-one single-clock branch" and theorem "exact and
  approximate growing-mode coherence": with zero initial orthogonal component
  and no orthogonal forcing, the source-side isocurvature fractions vanish
  and the primordial covariance has rank one.
- Homogeneous radial source family and thin-shell lift:
  `theorem homogeneous radial source family` in the screen-spectrum package:
  the source dilation cocycle forces an exact power law
  `Delta_zeta^2(k) = A_zeta (k/k_star)^(-theta)`, hence zero running.
- Screen field: the primordial source is the scalar collar-volume field
  `q_r = Pi_{>=2}[(1/3) log(J_X/J_bar_X)]`; no transverse-traceless source is
  emitted by the construction, so the primordial tensor amplitude of the
  branch is zero.
- Comparison pixel: `P_C = 1.6309682094039593` (CODATA-located, target-informed
  through alpha) and `P_fwd = 1.630972095858897` (forward map), from
  `extra/fine_structure_constant_derivation.tex`.
- Retrospective ledger: `code/cosmology/postdiction_ledger/` (seen data only).

## Proposed rows

### FZ-15: scalar spectral index

- id: FZ-15
- owning_issue: 742
- milestone: freeze and anchor before the first Simons Observatory or CMB-S4
  primary-CMB parameter release that is not public at anchoring time.
- content (proposed, verbatim-ready): Conditional edge-center tilt, not an
  independent derivation of the reserve-generator receipt: under the
  full-collar density `P_star/24`, the orientation-half identity, the
  dilation cocycle and the rank-one single-clock branch, the scalar spectral
  index is `n_s = 1 - P_star/48`. At the comparison pixel `P_C` this is
  `n_s = 0.9660215`; at `P_fwd` it is `0.9660214`; the branch spread is below
  `1e-7` and the row carries no continuous freedom. The `e(P_star - phi)`
  clock branch of the inflation paper (`n_s = 0.96484` at `P_C`) is a
  separately named diagnostic alternative; this row freezes the edge-center
  branch only, and the alternative may not be substituted after exposure.
- comparison_protocol (proposed): post-freeze published primary-CMB posteriors
  of the frozen experiment combination (named at anchoring: Simons
  Observatory large-aperture telescope with Planck, or CMB-S4 with Planck,
  with the named fallback), base LambdaCDM `n_s` at the pivot
  `k = 0.05 Mpc^-1`, marginalized 68 percent interval; no reanalysis and no
  target-dependent nuisance change. Planck 2018, ACT DR6, SPT-3G and every
  combination published before the freeze are seen data, postdiction
  bookkeeping only. A shift of the pivot or of the running prior is a
  declared no-verdict condition unless the frozen protocol names the
  conversion.
- kill_band (proposed structure; numeric cells OWNER SLOTS): FZ15-R01 (FAIL):
  the frozen value lies outside the frozen combination's two-sided credible
  interval at the declared level. Direction-neutral: a measurement above or
  below the value fails the branch equally. Compatibility supplies no
  confirmation credit for the reserve-generator receipt, which is not
  supplied; it records only that the receipt-conditional value survives.
  Seen-data context for the owner: the value sits at `+0.27 sigma` from
  Planck 2018, `+0.69 sigma` from SPT-3G with Planck, `-0.79 sigma` from
  the SPT+Planck+ACT combination, and `-1.28 sigma` from ACT DR6 with Planck,
  so the seen data do not select between the edge-center value and the ACT
  central value and a future two-sigma band of width `0.002` will.

### FZ-16: primordial tensor amplitude

- id: FZ-16
- owning_issue: 742
- milestone: freeze and anchor before the first LiteBIRD or CMB-S4 B-mode
  tensor-to-scalar release that is not public at anchoring time.
- content (proposed): Conditional exact zero: the primordial source of the
  branch is the scalar collar-volume screen field and no transverse-traceless
  source is emitted, so the primordial tensor-to-scalar ratio is `r = 0`
  exactly at the source; the second-order tensor background sourced by scalar
  modes (of order `1e-6` or below in the standard transfer) is not an
  exception to the row. This differs from generic single-field slow-roll
  inflation, which gives `r > 0` (Starobinsky: `r` near `0.003` at fifty-five
  e-folds), so a null at the LiteBIRD design sensitivity (`delta r` below
  `0.001`) discriminates, and a detection kills the scalar-only branch.
- comparison_protocol (proposed): post-freeze published tensor-to-scalar
  posteriors or upper limits of the frozen experiment at the pivot
  `k = 0.05 Mpc^-1`; seen data: BK18 (`r < 0.036`), the fourth Planck release with BK18 and
  BAO (`r < 0.032`), Planck 2018 with BK15 (`r_0.002 < 0.056`).
- kill_band (proposed structure): FZ16-R01 (FAIL): a detection of `r > 0` at
  the declared significance (OWNER SLOT, at least five sigma proposed) with
  the frozen foreground and lensing-delensing protocol. A tightened upper
  limit is COMPATIBLE and earns no confirmation credit beyond the
  discrimination against models with `r` above the limit.

### FZ-17: running and isocurvature exact zeros

- id: FZ-17
- owning_issue: 742
- milestone: freeze and anchor before the first CMB-S4 primary-CMB parameter
  release that is not public at anchoring time.
- content (proposed): Conditional exact zeros of the same branch: the running
  `dn_s/dln k = 0` exactly (dilation-cocycle power law) and the source-side
  isocurvature fraction `beta_iso = 0` exactly (rank-one single-clock
  branch). Slow-roll inflation gives a running of second order in
  `n_s - 1`, of order `-5e-4` for the measured tilt, below the projected
  CMB-S4 sensitivity, so the running row discriminates only against models
  with running of order `1e-2` and the row states that explicitly; the
  isocurvature row is shared with adiabatic single-field models.
- comparison_protocol (proposed): post-freeze published running posterior
  and isocurvature bounds of the frozen experiment; seen data: Planck 2018
  running `-0.0045 +- 0.0067`, ACT DR6 P-ACT-LB running `0.0062 +- 0.0052`,
  Planck 2018 X uncorrelated CDI `100 beta_iso < 3.5, 3.8, 3.9` at the three
  reference scales.
- kill_band (proposed structure): FZ17-R01 (FAIL): running excluded from
  zero at the declared level in either direction; FZ17-R02 (FAIL): a detected
  nonzero isocurvature fraction at the declared level.

## Eligible freeze boundary

An owner may freeze the conditional values, their named premise exits, the
post-freeze data release, the pivot and parameter conventions, a
direction-neutral nuisance protocol and complete decision thresholds before
accessing the comparison. Freezing cannot turn the reserve-generator receipt,
the rank-one single-clock branch, the dilation cocycle or the scalar-only
source into derived facts. Per the corpus kill-condition protocol, the audit,
re-audit and repair ladder runs before any failure is banked as a framework
kill.

## Custody, attestation, content_sha256

OWNER SLOTS.
