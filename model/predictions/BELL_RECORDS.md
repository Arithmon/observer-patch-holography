# Classical record responses and published Bell measurements

The diagonal native field branch is a classical stochastic process. It cannot
be a complete model of the Bell observations under a **local factorized
response** interpretation: after independently chosen settings, each side
uses only its setting, local state and records, with no communication before
the outputs. This statement concerns that interpretation and access contract.
It does not impose a fundamental spacetime ontology, exclude the full quantum
C/E instruments, or identify native repair distance with laboratory separation.

## Bound from the recorded branch

The existing binary endpoint map is identity or a local exchange. Diagonal
preparations remain diagonal under its recorded operations and local classical
feedback. Shared initial states, randomness and past records can be collected
in a variable lambda. Under the stated local-access and setting-independence
conditions, a trial has response

```text
P(x,y | a,b,lambda) = P_A(x | a,lambda) P_B(y | b,lambda).
```

Local random choices can be included in lambda. The resulting deterministic
response tables assign two binary outputs to each side; there are sixteen
tables. Each wins either one or three of the four settings in the game
`(-1)^(ab) x y = +1`. Convex mixtures therefore satisfy `P(win)<=3/4` for
uniform independent settings, equivalently `S<=2`. This is the established
CHSH argument applied to this source branch, not a new physical principle.
An executable control allowing one side to read both settings wins every
trial and demonstrates why that access condition matters.

The native control executes all 1,056 setting/coin/input histories on the
66-state carrier. Settings select `(0,1)` or `(0,5)` at the first region and
`(3,8)` or `(3,9)` at the other. Disclosed initial occupations and local signed
currents supply both output counters. Reversing the two disjoint repairs
preserves every result. With the uniform fixed-total-two preparation the
protocol has `S=26/33` and winning probability `79/132`. These values describe
this example; the more general upper bound covers all admitted local response
strategies. No laboratory realization of the example is claimed.

## Public experimental comparison

[Hensen et al., Nature 526, 682–686 (2015)](https://doi.org/10.1038/nature15759)
publish a Bell experiment with independent setting generation and separated
readouts. Their [deposited data](https://doi.org/10.4121/uuid:6e19e9b2-4a2d-40b5-8dd3-a660bf3c0a31)
include 4,746 candidate event-ready rows and an analysis example. The retained
[numerical file](data/hensen2015.txt) preserves every original byte from
`bell_open_data.txt`; its archive and file hashes appear in the receipt.
The original data remain subject to the repository's
[4TU General Terms of Use](https://doi.org/10.4121/resource:terms_of_use),
separately from this repository's code license. Credit belongs to the dataset
authors. No authors' analysis code or manuscript text is vendored.

The [independent replay](bell_records.py) uses their published event windows,
invalid-marker rule and excitation-click eligibility, then their readout
windows. It does not scan alternative thresholds. Eligibility does not inspect
the settings or final spin outcomes; both detected and undetected readouts
count as outcomes. This reproduces:

| Quantity | Result |
| --- | ---: |
| Eligible trials | 245 |
| CHSH-game wins | 196 |
| CHSH estimate | `16036094/6619647 = 2.42249987046...` |
| Published setting-calibration allowance `tau` | `27/2500000` |
| Conditional local win bound | `q=3/4+3(tau+tau^2)` |
| Fixed-horizon tail upper bound | `0.0391746435588...` |

The tail is calculated exactly as
`sum_(j=196)^245 binom(245,j) q^j (1-q)^(245-j)` and exported with outward
rational rounding. The statistical argument needs a conditional win bound
on every trial, given its past. Independent device trials are unnecessary:
induction on the remaining trials gives binomial-tail domination even with
device memory. The horizon must be fixed or have a separately justified
stopping rule. This package reproduces the authors' fixed-count analysis;
it supplies no optional-stopping certificate or new trial-registration claim.

The allowance for setting predictability is an external calibration premise,
not estimated from the winning records. The code also cannot authenticate
the physical separation, setting independence, event-ready timing or file
custody. Under that declared experimental/null contract the local-factorized
extension is rejected at the 5% threshold. This retrospective reproduction
is not a new prospective experiment or a posterior probability that a theory
is false.

The result prevents promoting the diagonal clock/field demonstrations to a
complete description of these measurements while retaining the local
factorization contract. The full OPH construction allows private quantum
states and incompatible instruments, and is not covered by that null. This
comparison neither selects C versus E nor establishes their physical
implementation, physical spacetime or a laboratory-clock calibration.

```sh
python -m model.predictions.bell_records --check model/predictions/bell_results.json
python -m pytest -q -W error model/predictions/test_bell_records.py
```
