# ACT DR6 retrospective fixed-tilt comparison

The conditional edge-center branch gives
`n_s = 1 - P_C/48 = 463690317905960407/480000000000000000`, or
approximately **0.966021496**. In the official ACT DR6 P-ACT-LB baseline
posterior, approximately **0.612%** of the retained multiplicity is at or below
this value. The historical clock alternative gives approximately **0.964841143**,
with **0.200%** at or below it. These are posterior positions conditional on the
published cosmological model and priors; they are not p-values, Bayes factors,
probabilities that OPH is true, or evidence that the native observer model has
derived our Universe.

This is a reproducible reanalysis of measured-Universe data already used in the
[cosmology postdiction ledger](../postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md).
It adds exact weighted marginal calculations from the released ordered chains,
including per-chain and leave-one-out sensitivity, rather than a Gaussian
approximation to rounded published summaries. It does not add independent data
to that ledger or refit a likelihood. Both exposed formulas remain visible;
the clock alternative does not replace the edge-center commitment.

## Data and provenance

The [ACT collaboration release page](https://act.princeton.edu/act-dr6-data-products)
links the official [NASA LAMBDA archive listing](https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_chains_lcdm_get.html)
and [NERSC chain documentation](https://portal.nersc.gov/project/act/dr6.02/chains/README_chains.html).
The selected product is exactly `p-actbase-l-b_lcdm_camb`, not the smaller CLASS
or foreground-marginalized ACT-lite alternatives. It combines Planck and the
full multifrequency ACT DR6 primary CMB likelihood, ACT/Planck CMB lensing,
and DESI 2024 BAO. No Pantheon+ supernova likelihood is included.

The [ACT DR6 parameter paper](https://arxiv.org/html/2503.14452v1), section VIII.3
and Appendix H, supplies the physical interpretation. The dimensionless scalar
spectral index describes the primordial scalar power law with amplitude pivot
`k_* = 0.05 Mpc^-1`. The chain assumes baseline six-parameter flat Lambda-CDM
and the collaboration's standard radiation transfer, recombination, foreground,
calibration, and nuisance models. Its configuration gives a uniform `n_s` prior
on `[0.85, 1.1]`. The normal reference distribution in that configuration only
initializes the sampler; it is not the prior. Marginalizing the released joint
posterior retains the effect of its fitted parameter correlations, without
inventing a diagonal observational covariance.

The archive contains four ordered Cobaya chains, configuration, and convergence
records. This package retains only the original ASCII `weight` and `ns` tokens
from every numbered-chain row, together with member sizes/hashes and the
checkpoint. There is no thinning, rounding, or removal of initial rows in the
stored projection. The full 163.5 MB archive remains outside the repository.

| Artifact | Size in bytes | SHA256 |
| --- | ---: | --- |
| Official archive | 163494761 | `8f2801a69663ef1142bb05660e0118ef7b1e0bdf3d05afba434a0569023b98ab` |
| Retained `act_dr6_ns_weight.json.gz` | 567477 | `13dab540b652300bd000cf1ac3b197af55c7c411de215682482ea55993f40485` |

The archive is available from
[NASA LAMBDA](https://lambda.gsfc.nasa.gov/data/act/chains/lcdm/p-actbase-l-b_lcdm_camb.tar.gz)
and its [NERSC mirror](https://portal.nersc.gov/project/act/dr6.02/chains/compressed/lcdm/p-actbase-l-b_lcdm_camb.tar.gz).
The listing dates delivery to 18 March 2025. The FZ-15 stamp is
4 October 2026, and this reanalysis is dated 18 October 2026. Consequently the
data are retrospective relative to that commitment. The projection is a
lossless projection onto two columns, not a copy of the complete posterior
state or a replacement for the full likelihood.

## Calculation and its limits

The official README prescribes **50% burn-in for this exact CAMB product**.
The calculation discards `floor(0.5 * stored_row_count)` initial rows separately
from each chain, then applies the remaining original integer multiplicity
weights. It does not remove half the total weight and does not apply the 10%
burn-in recommended for the separate CLASS product. Of 758350 source rows,
379175 remain, with total multiplicity 508985. Comparisons with the fixed
rational edge prediction use the exact decimal tokens, without a float
comparison at the decision boundary. Quantiles invert the weighted empirical
CDF without interpolation.

| Retrospective result | Value |
| --- | ---: |
| Pooled posterior mean of `n_s` | 0.9742582419612366 |
| Pooled posterior standard deviation | 0.0033775031447940 |
| Edge-center mass at or below prediction | `3114/508985 = 0.00611805849` |
| Historical clock-alternative mass at or below prediction | `1020/508985 = 0.00200398833` |
| Edge-center per-chain lower-tail mass range | 0.00372945 to 0.01031207 |
| Clock-alternative per-chain lower-tail mass range | 0.00077226 to 0.00377020 |

The mean and standard deviation reproduce the published rounded
`0.9743 +/- 0.0034`. The fixed branch lies in the lower posterior tail, so this
comparison supplies tension rather than a numerical confirmation. A tail mass
alone does not quantify a model's relative evidence.

The released checkpoint reports `converged: true` with
`Rminus1_last = 0.009799812902692016` and four MPI chains. This ordinary
convergence record does not certify rare-tail accuracy. Many consecutive rows
have identical `n_s` during fast nuisance-parameter updates. The 2319 retained
rows below the edge prediction occur in only 85 contiguous excursions, so
neither stored rows nor summed multiplicities count independent observations.
Per-chain and leave-one-out sensitivity is reported without treating it as a
confidence interval. This package supplies no rank/folded split-Rhat or tail-ESS
certificate and cannot satisfy the complete
[FZ-15 adjudication requirements](../../../docs/FROZEN_PREDICTION_LADDER.md).
No prospective FZ-15 verdict is issued or changed.

The model-side coordinate `P_C = 1.6309682094039593` is inherited from a
CODATA-alpha calibration. A physical primordial generator, mode embedding,
and standard CMB transfer remain additional hypotheses. No joint uncertainty
propagation over that calibration and these missing theoretical premises is
claimed. Standard Lambda-CDM treats `n_s` as a fitted primordial input; the
conditional OPH branch fixes that input, but this comparison does not establish
its derivation from the runnable native observer construction.

## Reproduce

From the repository root, obtain the pinned archive outside the repository:

```sh
curl --fail --location \
  https://portal.nersc.gov/project/act/dr6.02/chains/compressed/lcdm/p-actbase-l-b_lcdm_camb.tar.gz \
  --output /tmp/p-actbase-l-b_lcdm_camb.tar.gz
python3 code/cosmology/tilt_posterior/project_archive.py \
  /tmp/p-actbase-l-b_lcdm_camb.tar.gz \
  /tmp/act_dr6_ns_weight.json.gz
```

The generator checks the archive size and pinned SHA256 before processing. It
streams regular tar members without extracting paths or executing configuration
code, hashes the original member bytes, preserves original token strings,
and writes sorted compact JSON with gzip timestamp zero and no filename header.
The rebuilt projection reproduces the retained bytes and digest. Only Python's
standard library is required for this projection step.

Replay the retained comparison using the repository's Python dependencies:

```sh
python3 code/cosmology/tilt_posterior/audit.py --check
```

To analyze the independently rebuilt projection without changing the retained
receipt:

```sh
python3 code/cosmology/tilt_posterior/audit.py \
  --projection /tmp/act_dr6_ns_weight.json.gz \
  --output /tmp/act_dr6_tilt_results.json
```

`results.json` records the input and analysis hashes, exact comparisons,
weighted posterior summaries, and sensitivity checks. CI can replay it from
the small retained projection without a network download.
