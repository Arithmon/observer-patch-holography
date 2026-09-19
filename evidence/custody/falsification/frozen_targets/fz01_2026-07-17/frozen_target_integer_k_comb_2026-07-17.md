# FROZEN TARGET DRAFT: Integer-k Kerr Transition Comb (Ringdown)

Date: 2026-07-17.
Status: UNANCHORED DRAFT. This document carries no verdict weight. Anchoring
(hash commit of this text plus the completed likelihood artifacts, before any
decision data is opened) is the owner action. Until anchored, nothing here is
a registered prediction.

Companion statement: `INTEGER_K_COMB_STATEMENT.md` (same directory).

## 1. Supersession record

This draft SUPERSEDES the void MM-01 alpha = 4 line and the 2026-07-16
alpha = 4 displacement audit.

Reason: the canonical OPH relation is alpha = 4 ln k with integer k >= 2,
because horizon entropy is a record count and admissible dominant transitions
divide an integer edge dimension by an integer factor. Reading alpha = 4
requires k = e, which is no integer dimension ratio; alpha = 4 is therefore
inadmissible as an OPH forecast. The historical MM-01 arithmetic (nominal
-5.9 sigma on GW150914 under frequency-error-only propagation) and the
2026-07-16 stacked ln BF = -0.01 displacement construction are preserved as
audit artifacts only, per
`reverse-engineering-reality/docs/OPH_FALSIFICATION_PROGRAM.md` (Live
Tensions rows: VOID/MISATTRIBUTED and VOID as an OPH forecast) and tracker
row MM-01 (`oph-physics-sim/docs/OPH_SIGNATURE_EXPERIMENT_TRACKER.md`).
The displacement audit additionally substituted Schwarzschild spacing and a
uniform nearest-tooth offset for the derived Kerr transition spectrum; no
part of it transfers to this registration. Post-hoc bundles, secondary
lines, or fitted normalizations cannot rescue or promote the alpha = 4 line.

## 2. Observable

The offset-subtracted transition ratio, pinned to the integer-k ladder. For
a Kerr remnant with inferred (M, chi) and azimuthal family m, candidate
spectral features f_a, f_b above the rotation line satisfy, under the comb
hypothesis,

    (f_a - m Omega_H/(2 pi)) / (f_b - m Omega_H/(2 pi)) = ln k_a / ln k_b

for integers k_a, k_b >= 2, equivalently fixed universal-coordinate
positions x_k = ln k/(8 pi) with
x = (GM/(c^3 g(chi))) (omega - m Omega_H). Reference ladder against the
k = 2 tooth: ln 3/ln 2 = 1.585, ln 4/ln 2 = 2.000, ln 5/ln 2 = 2.322. The
ratio is identical for every remnant mass and spin and is redshift-free.
Secondary structure available to the likelihood: the (k-1)/k KMS weight
hierarchy and the mass-independent fractional linewidth
64 pi^2 p_0/(a ln k) with declared a in [1, 10].

## 3. Decision data

- Loudest available ringdown events of the GW250114 class (per-event
  frequency precision near 1.2% of f220; recorded in
  `verdicts/2026-07-16/V2_ringdown_comb.md`), individually and stacked in
  the universal x coordinate.
- Next-generation detectors (Einstein Telescope / Cosmic Explorer class),
  where single-event ringdown frequency precision near 0.3% makes one event
  decisive under the recorded Monte Carlo power estimates (expected
  |ln BF| approx 1.8 per event at 0.3%).
- Event selection, catalog versions, and data cuts are part of the freeze
  and must be fixed before any decision data is opened.

## 4. Decision policy

Preconditions (all frozen before data, per the falsification program's
ringdown-template-completion requirements):

1. Frozen strain likelihood. The registration requires a derived
   strain/asymptotic-readout likelihood: the map from the transition comb
   to observable strain, with exterior-background bridge assumptions
   declared. A frequency-list comparison without a strain likelihood does
   not qualify and returns no verdict.
2. Frozen prior/weighting across the allowed k lines (the corpus candidate
   is the KMS (k-1)/k hierarchy times GR greybody factors; the frozen
   artifact must normalize it into a proper pre-data prior).
3. Frozen remnant covariance treatment (full (M, chi, f) covariance;
   frequency-error-only propagation is inadmissible), nuisance model,
   selection effects, and trials accounting.
4. Frozen confidence level: 99% (approx 3 sigma equivalent) for both the
   kill and the detection directions, stated here as the draft value and
   fixed at anchoring.

SNR / precision gate: a per-event verdict contribution requires ringdown
SNR >= 30 and at least two resolved post-rotation-line spectral features,
each with frequency uncertainty below one third of the local k = 2 tooth
spacing. Events below this gate contribute to the stack only through the
frozen likelihood and can never singly trigger a verdict. If no event or
stack in the decision dataset clears the gate at the frozen confidence
level, the verdict is INCONCLUSIVE. The gate values are draft numbers and
are fixed at anchoring.

Verdict conditions:

- KILL: a resolved transition pair (both features clearing the gate, in one
  event or in the coherent stack) whose offset-subtracted ratio excludes
  every integer-k ratio ln k_a/ln k_b with k_a, k_b in {2, ..., 12} at the
  frozen 99% confidence, under the frozen likelihood and covariance. This
  kills the discrete-horizon integer-division continuation template. It
  does not by itself kill the derived area-spectrum statement; that scope
  boundary is recorded in PAPER.tex.
- DETECTION CANDIDATE: coherent stacking of independent events at the
  predicted x_k positions, with the offset-subtracted ratio of the two
  strongest features consistent with a single integer pair, at the frozen
  confidence with trials correction. A single-event pair passing the gate
  with ratio on the ladder is reportable as a candidate, never as a
  confirmation, until an independent event repeats it.
- INCONCLUSIVE: everything else, including all datasets in which only one
  feature clears the gate. The recorded power estimate is that a decision
  at |ln BF| approx 2 needs approximately 4 to 5 GW250114-class events or
  one ringdown at 0.3% frequency precision.

## 5. Conditionality

The comb is a continuation-level template: it adds the integer-division
selection rule and, for strain contact, standard semiclassical exterior
inputs. A negative verdict kills the template branch. The derived content
of the axioms is the discrete area spectrum Delta A = 4 l_P^2 ln(d'/d),
which a comb non-detection leaves standing; a resolved pair off the integer
ladder (the KILL condition above) contradicts the integer-transition
reading directly.

## 6. Owner actions to reach anchored status

- Produce the derived strain/asymptotic-readout likelihood artifact.
- Normalize the (k-1)/k greybody weighting into a pre-data prior.
- Fix event list, catalog versions, SNR gate, confidence level, and trials
  accounting.
- Hash-commit this document plus the artifacts; record the hash in the
  tracker; open the data once.
