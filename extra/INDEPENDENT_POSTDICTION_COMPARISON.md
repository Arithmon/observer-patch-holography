# Independent alpha/tau reproduction and the first bounded natural comparison

Issue #1011 is resolved at its explicitly permitted **positive reproduction /
precise NOT_READY** exit. Both certified alpha roots and the conditional tau
window are independently reproduced. The selected prospective route is the
existing FZ-10 target, owned by #546. No new natural outcome has been evaluated,
no new prediction registry is created, and no source-selection parent is closed.

The executable evidence is [the independent packet](../code/independent_postdictions/README.md).
Its [contract](../code/independent_postdictions/CONTRACT.md) fixes the bounded
three-candidate menu and records that the first replay preceded ranking.
Its receipt retains the approximate and less successful rows, not just the
closest endpoint. These are independent implementations by the same assistant,
not a claim of institutional independence or blindness to historical results.

## 1. What is reproduced, and from which inputs?

Primary mathematical sources are `code/P_derivation/paper_math.py`, its retained
July 2026 interval certificate and approximate trunk, and the positive-circulant
identity and tau artifact in `code/particles/leptons/` and `runs/leptons/`.
The fresh producer imports none of their code or output numbers. The verifier
uses different numerical formulations and reads the old certified enclosures
only to compare the newly reconstructed roots. Source hashes and narrow
projections of the canonical claim and frozen-target rows are in `receipt.json`.

| Input | Value and standard uncertainty | Role |
|---|---|---|
| CODATA 2022 inverse alpha | 137.035999177 +/- 0.000000021 | Historical comparison; also the explicitly marked mixed diagnostic |
| CODATA 2022 electron mass | 0.51099895069 +/- 0.00000000016 MeV | Conditional tau calibration |
| CODATA 2022 muon mass | 105.6583755 +/- 0.0000023 MeV | Conditional tau calibration |
| PDG 2024 tau mass | 1776.93 +/- 0.09 MeV | Exposed historical comparison only |

These are attributed transcriptions of the [NIST CODATA table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)
and [PDG 2024 tau listing, page 1](https://pdg.lbl.gov/2024/listings/rpp2024-list-tau.pdf).
They reproduce the historical artifact's edition, rather than choosing a newer
reference after seeing its residual. No covariance matrix is supplied here.
Input marginal standard uncertainties do not define a joint probability law.

The forward alpha calculation has no measured alpha or electroweak-scale input.
It does have declared running coefficients, closure, continuation masses,
transport, and discrete selectors. Solving these equations is not fitting their
output to CODATA; selecting this map historically is a separate ancestry issue.
The tau calculation consumes two measured masses and the balance and ordering
premises. Its code does not consume measured tau; the **balance premise itself
has measured-triple ancestry**, as retained in the old artifact.

## 2. Dimensionless alpha reconstruction

Write `g=alpha_U`, `z=m_Z/v`, `phi=(1+sqrt(5))/2`. Eliminating the energy unit
from `v=P^(-1/2) exp(-pi/(2g))` and `M_U=exp(-2pi) P^(1/6)` gives

```
L = log(M_U/m_Z) = -2pi + (2/3)log(P) + pi/(2g) - log(z)
alpha_i = [1/g + b_i L/(2pi)]^-1,   b=(33/5,1,-3).
```

The three simultaneous equations are

```
z^2 = pi (alpha_2 + 3 alpha_1/5)
ell_SU2(4pi^2 alpha_2) + ell_SU3(4pi^2 alpha_3) = P/4
(P-phi) A(P,g,z) = sqrt(pi).
```

Here `ell(t)=sum d_R log(d_R) exp(-t C_R)/sum d_R exp(-t C_R)`.
SU(2) uses `d=n+1, C=n(n+2)/4, n=0..120`; SU(3) uses
`d=(p+1)(q+1)(p+q+2)/2, C=(p^2+q^2+pq+3p+3q)/3`, `p,q=0..90`.
The fresh solver collects degenerate Casimirs; the verifier separately sums
every representation. The retained interval certificate supplies the existing
rigorous tail/enclosure control for the first two maps; this floating-point
replay is not a new interval proof or a proof of global uniqueness.

The declared stage-five square-root profile is
`r_k=1+sqrt(2)cos(2/9+2pi k/3)`, sorted increasingly. Holding muon and tau
exponents at 4 and 3, the complete local electron menu is 5 through 11. The
maximum pairwise logarithmic mismatch selects 7. With `(n_e,n_mu,n_tau)=(7,4,3)`,

```
m_k/v = r_k^2 * 2^(1/6) / [product_j(r_j^2 sqrt(2) 6^n_j)]^(1/3).
m_q/v = 6^(-n_q)/sqrt(2),  n_(u,c,t,d,s,b)=(6,3,0,6,4,2).
```

Thus all internal masses in this replay are **dimensionless ratios to v**.
Old variable names ending in `_gev` in the Planck-normalized calculation do not
make their numerical values physical GeV. No measured Z mass is inserted.

The structured readout is

```
A = 1/alpha_2 + 5/(3alpha_1)
    + sum_leptons K(m/v) + (1-3alpha_3/pi) sum_(u,c,d,s,b) w_q K(m_q/v)
K(r) = (2/pi) integral_0^1 x(1-x) log[1+(z/r)^2 x(1-x)] dx
w_(u,c)=4/3, w_(d,s,b)=1/3.
```

The gauge-width map adds `g` and is solved again self-consistently. The
asymptotic map replaces each kernel by `[log(z^2/r^2)-5/3]/(3pi)`.
The producer uses a closed integral formula; the verifier performs direct
quadrature. It reconstructs the running through `M_U/(zv)` rather than copying
the producer's eliminated equations. Printed roots satisfy all equations to
`1e-38`, with fresh roots inside both old certified enclosures.

| Declared map | Fresh inverse-alpha output | Absolute residual, ppm of reference |
|---|---:|---:|
| Structured | 136.9948351774129372952894294644 | 300.3882179448 |
| Structured plus gauge width | 137.0356601369465765688727304625 | 2.47409480326 |
| Converged asymptotic | 136.9940258961185127104409023921 | 306.2938288739 |

The historical approximate trunk printed `136.994020662724205139718642793`.
It differs from the converged asymptotic solve by `-0.00000523339430757`;
its printed `(P,A)` pair has signed relative outer-equation defect
`-0.000005245187842595`. This preserves its already-declared approximate status
instead of inventing an exact root at the old display.

The separate mixed diagnostic gives `137.0359595136085677774670579934` by
adding `g(P_C)` to the structured fixed-point output, with
`P_C=phi+sqrt(pi)/alpha_inverse_measured`. **This consumes the measured endpoint**
and is neither the forward gauge-width fixed point nor an independent alpha
prediction. It is retained because confusing these quantities would overstate
the agreement.

The nearest forward result is still about **16,145 reference standard
uncertainties** from CODATA. That number is not a theory significance: no
physical truncation/hadronic/readout error law is supplied. If a symmetric
additive theory remainder `[-R,R]` were proposed, mere overlap with the reference
two-standard-uncertainty interval would require
`R >= max(|A-A_ref|-2 sigma_ref,0) = 0.000338998053423...`.
This is a necessary overlap bound, not a derived uncertainty or permission to
invent that remainder. The certified arithmetic error is a different quantity.

Known trials include the seven electron exponents, the two certified map
variants, and the approximate/asymptotic continuation. The existing
`selection_accounting.json` and `tracking/null_model_scorecard.md` retain an
audit grid of 48 constant-pair labels, 42 distinct pairs, and 41 alternatives;
zero alternatives meet its after-the-fact `2.5e-6` threshold. That grid is
neither exhaustive nor a precomparison probability measure. This task does
not rerun its certificates or reproduce a purported exhaustive 112-model
search. The full historical search is unknown; **no global chance probability
or positive uniqueness weight follows**.

## 3. Tau: one constraint, two inputs, and an exact reference equivalence

On the positive chamber, the declared circulant identity is
`Q=1/3+(2/3)(rho/a)^2`. The supplied balance `rho/a=1/sqrt(2)` gives
`Q=(e+m+t)/(sqrt(e)+sqrt(m)+sqrt(t))^2=2/3`.
Put `x=sqrt(t)`, `p=sqrt(e)+sqrt(m)`, `s=e+m`. Clearing the positive denominator
gives `x^2-4px+3s-2p^2=0`, hence

```
t = [2p + sqrt(6p^2-3s)]^2                         (ordered branch).
```

The plus branch is selected by `t>m` for these inputs; the other positive root
is `3.31735654442055 MeV`, below the muon mass. The fresh producer instead
solves the original ratio equation, at the center and all four input corners.
The verifier independently evaluates the quadratic and analytic derivatives.

The result is `1776.9690272931573454044081080674 MeV`. All first derivatives of
the plus-root formula with respect to positive `e,m` are positive: writing
`a=sqrt(e), b=sqrt(m), D=sqrt(3a^2+12ab+3b^2)`,

```
J_e = sqrt(t)/a * [2+(3a+6b)/D] = 304.5561230227...
J_m = sqrt(t)/b * [2+(6a+3b)/D] = 15.3451267855...
```

Therefore the rectangular marginal-error input box has its exact extrema at
the both-low and both-high corners. Rounding outwards reproduces
`[1776.968991,1776.969063] MeV`, width **72 eV**. This is a deterministic box
enclosure, **not a joint 68% confidence interval**. The measured triple gives
`Q=0.6666644634026367...`, consistent with the old displayed comparison.

There is a covariance-independent *first-order* sensitivity bound. For every
positive semidefinite input covariance with the reported marginal errors,
Cauchy-Schwarz gives `sigma_prediction <= |J_e|sigma_e+|J_m|sigma_m`, here
`0.0000353425205864 MeV`. Representing covariance as Gram vectors and applying
the triangle and reverse-triangle inequalities also gives

```
sigma_tau - sigma_prediction_max <= sigma_(tau-prediction)
                                <= sigma_tau + sigma_prediction_max.
```

Even permitting correlations between the historical tau measurement and the
inputs, the standardized residual lies between `0.43346637` and `0.43380695`
in this linearization. This is a sensitivity calculation, not a Gaussian tail
area, exact nonlinear error distribution, or fresh evidence for the premise.

**Reference theorem.** Conditional OPH balance plus ordering and the balanced
Koide reference map send the same `(e,m)` to exactly the same `t`. Given the
same observation/error/nuisance model, their likelihoods are equal for every
possible observation; their likelihood ratio is identically one. Integrating
identical input/nuisance laws preserves that equality. No precision improvement
in this observable can distinguish these two formulations. The formula also
contains no P, so this test independently constrains no pixel constant.

Against a free tau-mass parameter (the unfixed charged-lepton Yukawa baseline),
there is one output constraint from two inputs. Under a **declared illustrative
normal error model**, `2 log(L_free/L_balanced)=(y-t)^2/sigma^2`, which is
`0.188040692739...` for the historical PDG point with fixed calibration center.
This uses the same data and sigma for both reference models. No model prior,
selection penalty, null tail area or Bayes factor is inferred. The constraint
is useful for rejecting the shared balanced relation; historical agreement is
a standard-relation postdiction, not distinct OPH confirmation.

## 4. The bounded candidate decision

Apply the declared lexicographic rule: defined observable/physical dictionary,
independently constrained relation, uncertainty treatment, eligibility, then
precision and cost/time. Exactly these three routes were ranked after replay:

| Rank | Candidate and discrimination | Attachment / uncertainty / readiness |
|---|---|---|
| 1 | FZ-10 tau: one mass constraint against a free tau parameter; exactly equivalent to Koide | Physical mass dictionary explicitly conditional; a frozen numerical rule exists. Needs a new eligible dedicated release and uncertainty/exposure dossier. **Selected, NOT_READY.** |
| 2 | Declared alpha/P maps: fixed numerical outputs against a free electromagnetic coupling | Missing map selection, same-quantity/Thomson attachment and physical transport error budget. Tiny experimental error cannot fill these gaps. NOT_READY. |
| 3 | FZ-11 photon dispersion: linked coefficient manifold against the zero-coefficient Lorentz-invariant baseline | Canonical register says unarmed: same physical sector, polarizations, frame/boost, carrier isolation, nuisance likelihood and finite scale remain required. No justified sensitivity/date forecast. NOT_READY. |

For FZ-10 the physical chain is: declared finite positive-chamber tracial-GNS
balance contract -> balanced square-root-mass profile -> identification with
charged-lepton pole masses -> calibrated electron/muon inputs -> ordered tau
output -> a dedicated mass estimator with its experimental calibration and
uncertainty. This states the physical consumer and its missing source premise;
it does not identify an abstract GNS coordinate as a particle pole by fiat.
No apparatus or self-reading patch simulation is consumed. Source/outcome
admission remains with #730, calibration with #736, common model integration
with #740, and the selected target with #546.

The [Belle II 2023 measurement](https://arxiv.org/abs/2305.19116), already exposed
in the PDG table, quotes `0.08 MeV` statistical and `0.11 MeV` systematic errors.
Its systematic error alone exceeds FZ-10's `0.045 MeV` compatibility gate.
Under unchanged systematics, more luminosity alone cannot meet that gate.
This is a concrete precision obstacle, not a forecast of a future experiment's
uncertainty or a claim that the target is permanently unreachable. Local
comparison arithmetic is cheap; no supported date for a decisive release is
available. No hardware campaign or funding is proposed.

## 5. Dataset discovery, protocol, and what is actually blocked

The checked-in INSPIRE snapshot requests only titles, identifiers, dates,
document type and collaboration metadata. The post-freeze literal `tau mass`
query returns zero; importantly, the same spelling also misses the known 2023
Belle II title. The broader `tau / Greek tau / lepton` query returns four rows:

| INSPIRE / arXiv | Metadata disposition |
|---|---|
| 3200778 / 2609.05689 | Vector-like top-partner recoil-mass search; other observable |
| 3193587 / 2608.19277 | Charged-lepton hierarchy / Koide theory; not a dedicated mass measurement |
| 3096034 / 2512.22756 | First-public date precedes freeze; not a new tau measurement |
| 3193022 / 2408.05903 | First-public date precedes freeze; not a new tau measurement |

The historical `lepton mass` control has 33 results, of which the requested
first 25 include Belle II `2663717 / 2305.19116`. That historical control is
truncated, explicitly; the four-row future-date query is not. The index's date
can reflect a later publication, so it cannot substitute for first-public date.
These controls justify **no eligible measurement identified by this bounded
discovery**, not completeness over all public, private or unindexed data.
The complete queries, retrieval time and retained metadata are versioned.

`protocol.json` is a concrete **future data-release contract**, not a claim that
a presently unidentified dataset is frozen. Before outcome access, publicly
register its exact commit/content hashes and timestamp evidence. Nominate the
first dedicated experimental charged-tau mass release first made public after
that registration, breaking ties by canonical identifier. Nomination is made
from metadata, never from value, uncertainty or agreement. Exclude world
averages, theory, reprints and updates of previously released measurements.
If that first candidate was exposed or fails admission, retain it and stop;
do not search for another favorable candidate.

Admission must lock public dataset/source bytes, first-public evidence,
event/run overlap against exposed samples, every participant's exposure
declaration, pole-mass/estimator convention, calibration and systematic
covariances, and the official combined standard uncertainty. A collaboration
aggregate is used unchanged; this protocol does not invent an error combination,
pool experiments, or silently turn asymmetric errors into a symmetric sigma.
If those conditions cannot be met, the outcome is NOT_READY. Metadata
authenticity and undisclosed exposure cannot be established by an offline JSON
checker; they require externally inspectable evidence before admission.

Use the immutable FZ-10 center, window and historical input calibration.
Do not reanchor to a later CODATA edition. The ordered numbered-bullet rule is:
distance above three reported standard uncertainties -> FAIL; otherwise distance
at most two uncertainties and sigma at most `0.045 MeV` -> COMPATIBLE; otherwise
INCONCLUSIVE. The first FAIL bullet is unconditional; the precision condition
belongs to COMPATIBLE. This interpretation is explicit before any new outcome,
including the synthetic coarse-error failure control. Historical PDG data yield
**INCONCLUSIVE**, despite the small normalized residual.

A separate full-window check reports FAIL only if every point of the frozen
window is more than three uncertainties away, and COMPATIBLE only if every point
is within two and the precision gate holds. At very small future errors a
rounded-center failure need not exclude the whole conditional window. Both
labels must be published; this diagnostic never rewrites the historical rule.
There is one release, one primary observable, one primary verdict, and a stop
after any outcome, including failure, inadequate precision or blocked admission.
Corrections retain the original result and receive a separate version.

Current blockers are specific: no eligible untouched dataset identified; no
dataset-specific overlap, exposure and uncertainty dossier; and no external
registration evidence yet for this new dataset-selection contract. The old
prediction freeze alone does not supply the new registration. No physical
comparison has occurred. Publishing this PR does not silently assert that
those admission obligations have been met.

## 6. Verification and exit coverage

The independent verifier checks equations, branch ordering, quadrature, original
enclosure agreement, analytic sensitivities, historical inputs, all retained
outcomes, the bounded discovery controls and normative protocol. Rebuilding
receipt hashes does not repair a false number or relax the reviewed contract.
Tests attack units, uncertainties, exact threshold edges, selectors, readout
components, mode changes, target leakage, omitted rows, policies, JSON duplicates,
nonfinite/oversized inputs, outcome-dependent selection and optimized-Python
execution. Synthetic controls are tests of code, never observations of nature.

| #1011 obligation | Delivered evidence / exit |
|---|---|
| Independent initial alpha and tau replay | Fresh simultaneous roots / direct Q solve, independent direct-sum and integral / quadratic checker; declared unit and selector ledger |
| Same-data reference comparison | Exact OPH/Koide equivalence; free-mass likelihood illustration; alpha residual and required-remainder diagnostic, with unknown search scope retained |
| At most three candidates after replay | Ranked table and prior contract; FZ-10 selected with concrete precision obstacle |
| Pre-outcome contract and exposure | Frozen target reused; protocol fixes future selection, errors, exclusions, multiplicity, stopping and outcomes; no claimed dataset freeze |
| Independent replay and hostile controls | Offline verifier, test suite, Windows/Linux CI; no synthetic physical evidence |
| Honest natural comparison or blocker | NOT_READY with a specific future-release acquisition and admission contract; no hardware or cloud expenditure |

This closes the finite reproduction/readiness task, not the experimental test,
physical mass-selection theorem, alpha endpoint attachment, or M1 source law.
