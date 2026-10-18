# OPH predictions and public-data comparisons

The target is a prediction for our observed Universe that adds information
beyond the Standard Model and General Relativity, followed by an executable
comparison with existing observations. Building a new laboratory realization
of the finite source is not a prerequisite. The observable dictionary, physical
premises and uncertainty model must nevertheless be explicit.

Two kinds of claim need separate interpretation. Fixing a parameter that a
standard model leaves free is a testable restriction of that model. Predicting
an intrinsic vacuum dispersion or another incompatible response changes its
dynamics. A successful parameter comparison does not establish incompatible
dynamics, and agreement alone does not uniquely identify OPH.

This audit is a retrospective assessment. It does not arm or change any
[frozen prediction](../FROZEN_PREDICTION_LADDER.md), replace an eligible future
dataset, or assign a prospective verdict to already exposed observations.

| Candidate | Additional prediction and physical premise | Public-data disposition |
| --- | --- | --- |
| Primordial scalar tilt | `n_s = 1 - P/48 = 0.9660214956374175...`, with the inherited alpha-calibrated `P = 1.6309682094039593`. The physical primordial generator, mode identification and standard CMB transfer remain conditional. | Executed the official ACT combined posterior comparison below: the fixed value lies at its lower 0.612% position. This is retrospective tension, not a new prediction score. |
| Charged-tau mass | Balanced-circulant mass interpretation and measured electron/muon masses give approximately `1776.969027 MeV`. This relation is algebraically the Koide relation; the balance premise has observed-mass ancestry. | A conditional mass restriction, not a change to SM dynamics or a unique OPH signature. The historical Belle II 2023 comparison below is compatible at its reported precision. |
| Intrinsic photon dispersion | Equal-weight primitive-port or seam-edge propagation gives linked fourth- and sixth-order coefficients below, if these directions are physical translations of the photon field. | A conditional departure from minimal locally Lorentz-invariant vacuum propagation. Published timing limits are far too weak at the proposed scale; the Auger transfer is source- and interaction-dependent. |
| Native clock ratios and noise | Same-history local movement counts have definite waiting distributions, correlations and covariance-derived spatial relations. | No identified atomic-clock, pulsar or photon observable has this record law. These laws admit ordinary stochastic realizations; matching them alone would not distinguish OPH. |
| C/E coherent record law | One fixed coherent preparation gives plus probabilities `5/8` and `3/8` under the two instruments. | These are ordinary quantum instruments. The distinction selects a record law, not OPH over standard quantum mechanics. |
| Joint timed coherent response | The charged total-two source history fixes a joint cost/clock/detector law. Under one explicitly synthetic calibration, its two settings have identical `0.52` deadline mass and distinct `0.2451`/`0.2715` joint detector masses. | An executable conditional forecast with shared timing coefficients and all timeouts retained. No measured system or physical duration is identified; these finite quantum instruments are not unique to OPH. |
| Fine-structure constant | Several declared endpoint maps exist, with incomplete physical transport and map-selection obligations. | Withhold a new prediction claim. The [existing endpoint audit](../../model/predictions/README.md#inherited-coupling-comparison) already finds a 300.388218 ppm gap for its named map, without a justified theory-error distribution. |
| Higgs mass | The historical near-match used a target-conditioned calibration. The [replacement electroweak producer](../../code/particles/calibration/source_ew_vev_matching.py) uses additional criticality, matching and energy-scale premises; it is not derived from the finite observer. | Withhold a distinctive prediction claim: the physical matching and theory-error budget remain open. The replacement's [retrospective comparison](../../code/particles/runs/calibration/source_ew_vev_matching_comparison.json) should retain those qualifications. |

## Executed CMB comparison

The [executable reanalysis](../../code/cosmology/tilt_posterior/README.md)
uses the collaboration's actual ordered posterior samples, not a Gaussian
distribution manufactured from its published mean and standard deviation.
We downloaded the official `p-actbase-l-b_lcdm_camb` archive: full
multifrequency ACT DR6 plus Planck primary CMB, ACT/Planck lensing and DESI
2024 BAO. The model is flat six-parameter power-law LambdaCDM with the
collaboration's nuisance treatment and priors, including `n_s` uniform on
`[0.85,1.1]`, at pivot `0.05 Mpc^-1`. These are conditional cosmological
inferences, not model-independent direct measurements of a repair event.
[ACT release](https://act.princeton.edu/act-dr6-data-products),
[likelihood and parameter paper](https://arxiv.org/html/2503.14452v1).

The [official burn-in instructions](https://portal.nersc.gov/project/act/dr6.02/chains/README_chains.html)
specify 50% of stored rows for this CAMB product. Applying that fraction
separately to all four chains leaves 379,175 rows with total multiplicity
508,985. The result is:

| Quantity | Actual weighted-chain result |
| --- | ---: |
| Posterior mean `n_s` | `0.974258242` |
| Posterior standard deviation | `0.003377503` |
| Fixed edge-center branch | `0.9660214956374175...` |
| Posterior mass at or below that branch | `3114/508985 = 0.00611805849...` |
| Historical clock alternative `1-e(P-phi)` | `0.96484114303...` |
| Posterior mass at or below that alternative | `1020/508985 = 0.00200398833...` |

The primary branch lies in the lower 0.612% of the empirical posterior.
This indicates tension with this particular dataset combination and its
cosmological assumptions. It is not a frequentist p-value, the probability
that OPH is true, a likelihood ratio or a Bayes factor. The alternative is
retained as an exposed diagnostic, not substituted for the fixed branch.
No parameters or likelihoods were fitted to improve either agreement.

The mean and standard deviation reproduce the published rounded
`0.9743 +/- 0.0034`, while the measured tail differs from a Gaussian shortcut.
The four individual edge-branch tail positions range from about 0.00373 to
0.01031. Only 85 contiguous tail excursions occur in the retained ordering;
the many repeated `n_s` values arise during nuisance-parameter updates.
Rows and multiplicities are not independent observations. Per-chain and
leave-one-chain-out sensitivity is retained; no rare-tail effective-sample-size
or confidence certificate is claimed. The release's ordinary convergence
checkpoint is not such a certificate.

The 163.5 MB source archive remains external. A 567 kB lossless projection
retains every original weight/`n_s` token, all chain ordering and burn-in rows,
with archive and original-member hashes. The independent archive parser and
exact-rational replay agree on both integer tail totals. Reproduce locally:

```sh
python code/cosmology/tilt_posterior/audit.py --check
python -m pytest -q code/cosmology/tilt_posterior
```

These observations were released in 2025 and were already exposed before
FZ-15. Its registered comparison concerns a specified future Simons
Observatory product, with additional convergence and eligibility conditions.
This analysis changes none of those terms and issues no FZ-15 verdict.
It supplies an actual-data constraint on the declared primordial branch;
deriving that branch from the runnable observer's joint records remains open.

## Conditional tau comparison

For `p = sqrt(m_e) + sqrt(m_mu)` and `s = m_e + m_mu`, the ordered Koide
branch is `m_tau = [2p + sqrt(6p^2 - 3s)]^2`. The
[independent executable reproduction](../../code/independent_postdictions/README.md)
gives `1776.9690272931573 MeV` from the declared electron/muon calibrations;
its first-order propagated uncertainty is at most `0.000035343 MeV` for any
input covariance with the stated marginal errors.

[Belle II's 2023 measurement](https://arxiv.org/abs/2305.19116) reports
`1777.09 +/- 0.08(stat) +/- 0.11(syst) MeV`. Combining the two reported
uncertainties in quadrature gives `0.1360 MeV`; the prediction is lower by
`0.12097 MeV`, about `0.89` of that experimental uncertainty. This is a
retrospective residual diagnostic, not a new likelihood analysis or an
eligible post-freeze measurement. The balanced OPH and Koide formulations
predict exactly the same mass for every allowed input; with the same
measurement model their likelihood ratio is one. The
[existing comparison](INDEPENDENT_POSTDICTION_COMPARISON.md) details this
equivalence and the premise's observed-mass ancestry.

## Why the native clocks and C/E test are not yet distinguishing predictions

The [single-carrier clock comparison](../../model/native_clock/COMPARISON.md)
already gives a nontrivial internal spacetime correspondence. If `n` counts
global movements and `K_v` counts movements incident to port `v`, then
`sum_v K_v = 2n` on every history and `K_v/n -> 1/6` in probability. The
asymptotic count covariance distinguishes three graph-distance classes and
recovers the finite seam geometry. These are dimensionless event counts;
neither external physical timestamps nor observer federations are needed to
derive them. They do not yet identify a measured physical clock or establish
Lorentzian spacetime.

The classical histories are finite Markov processes. More generally, for any
finite instrument with Kraus maps `K_r` and `sum_r K_r^* K_r = I`, the map

```text
V |psi> = sum_r K_r |psi> tensor |r>
```

is an isometry, extends to a unitary with an ancillary system, and reproduces
the instrument when that system is measured. This applies to both executed
[C/E laws](../../model/predictions/README.md#physical-scope). Their differing
interference probabilities therefore admit ordinary quantum realizations.
This is a statement about these finite laws, not a theorem excluding every
possible OPH distinction or every ultraviolet completion.

The [charged clock-history attachment](../../model/laboratory_clock/README.md#compiled-operations-require-separate-timing-coefficients)
strengthens the internal test: one shared two-resource calibration predicts
the joint timed response without choosing separate clocks for the two
coherent settings. It preserves the complete stopping distribution and
proves that calibration workloads with fixed pulse-to-repair ratio cannot
identify the two operation durations separately. Its additional attempt
counter is a declared controller readout, and its timestamps are synthetic.
The [protocol](NATIVE_CLOCK_EXPERIMENT_PROTOCOL.md#a-joint-timed-coherent-prediction-with-charged-operations)
states the exact conditional forecast and the physical correspondence needed
to use it. It changes no frozen prediction or empirical verdict.

Usable public clock observations exist. King et al.'s
[argon optical-clock paper](https://arxiv.org/abs/2205.13053) has
[public corrected argon/ytterbium ratio traces](https://zenodo.org/records/6901524),
including UTC timestamp columns verified by directly reading a file header.
Those are ratios of identified atomic transitions, not native incident
counts. No source-derived atomic-transition map or ratio prediction is
available here, so fitting an arbitrary time scale to those traces would
not test this source. By contrast, the
[2021 BACON publication](https://doi.org/10.1038/s41586-021-03253-4) supplies
its raw data on request rather than as an immediately public trace deposit.
The executed ACT comparison above is retained as the actual public-data
constraint; clock-data availability does not turn a missing observable map
into an additional physical comparison.

## Conditional dispersion: exact branch relation and actual sensitivity

In natural units, write the candidate physical dispersion as

```text
omega^2 = k^2 + C4 k^4 + [B0 + B6 I6(direction)] k^6 + O(a^6 k^8).
```

`I6` is the fixed normalized icosahedral angular invariant, not particle spin.
The [primitive-port producer](../../code/a5_fingerprint/spin_six_primitive_port_prediction.py)
and [seam-edge producer](../../code/a5_fingerprint/seam_current_edge_prediction.py)
give distinct alternatives:

| Coefficient | FZ-11 primitive ports | FZ-12 seam edges |
| --- | --- | --- |
| `C4` | `-a^2/20` | `-a^2/20` |
| `B0` | `a^4/840` | `a^4/840` |
| `B6` | `2 a^4/7875` | `-a^4/12600` |
| `B0/C4^2` | `10/21` | `10/21` |
| `B6/C4^2` | `32/315` | `-2/63` |

The finite geometry does not select the physical translation action, sole
propagation support, photon identification, kinetic term, physical clock or
preferred-frame transport. Those are premises of these branches. In
particular, taking the square root of a positive spatial operator as a wave
frequency is additional dynamics; it is not a consequence of the finite
occupation process relaxing by that operator.

Under the existing additional scale assumption `a = sqrt(P) l_P`, with
`l_P = 1.616255e-35 m` and `hbar c = 1.973269804e-7 eV m`, direct substitution
gives

```text
d = a^2/20 = 5.470944672629916e-58 eV^-2,
E_QG,2 = 1/sqrt(d) = 4.275322088547057e19 GeV,
v_group = 1 - (3/2) d E^2 + O(a^4 E^4).
```

The LHAASO collaboration's GRB 221009A analysis reports the subluminal
quadratic bound `E_QG,2 > 6.9e11 GeV` at 95% confidence using its ML/MINOS
procedure; its calibrated variant gives `7.2e11 GeV`. These are published
results under the paper's light-curve, spectral, propagation and detector
assumptions, not a new fit to the events here. [LHAASO, Table 1](https://arxiv.org/html/2402.06009v2)

For the quoted `6.9e11 GeV` bound, the proposed coefficient divided by the
bound's coefficient is `(6.9e11 / 4.275322088547057e19)^2 = 2.60472e-16`.
Thus this timing comparison is compatible but has negligible power to test
the assumed scale. It supplies no empirical preference over zero dispersion
and no test of the linked sixth-order angular coefficients. No claim that
this is the strongest currently available timing analysis is needed.

The [Auger electromagnetic analysis](https://arxiv.org/abs/2112.06773)
instead reports `delta_gamma,2 > -1e-58 eV^-2` under a high-energy
proton-containing source scenario, and no electromagnetic bound for its
reference composition without protons beyond `1e19 eV`. In the photon-only
branch `delta_gamma,2 = -d` to leading order, so the candidate magnitude is
5.47094 times that quoted limit. This is not an OPH exclusion: source
composition, lepton dispersion, interactions and shower response determine
whether that bound applies. The [existing feasibility calculation](DISPERSION_FEASIBILITY.md)
shows that shared photon/lepton dispersion substantially changes the pair
threshold; free shells alone do not supply the interaction rates.

A decisive test needs an isolated intrinsic propagation response and the
joint coefficient uncertainties. A resolved positive `C4`, or incompatible
linked sixth-order coefficients, can reject the specified cosine branch.
A null can reject a fixed-scale branch only with sufficient sensitivity.
If arbitrarily small positive `a` remains allowed, finite-precision null
measurements cannot exclude the entire family, whose closure contains zero.
Even a matching nonzero pattern would require comparison with other
Lorentz-violating theories and systematic explanations before supporting a
distinctive OPH interpretation.
