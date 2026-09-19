# DRAFT frozen target: corrected capacity band (unanchored)

Status: DRAFT. Not registered, not externally timestamped, carries no verdict
weight. Anchoring (OpenTimestamps or equivalent) and adoption into
`reverse-engineering-reality/falsification/frozen_targets/` are owner
actions. Until anchored, this file is a proposal.

## Target statement

Conditional on the readback map F, premise CP-2, and the corrected balance
being one of the three enumerated one-term candidates (declared menu M = 40,
selection recorded in
`reverse-engineering-reality/code/capacity_readback/CP1_CORRECTED_BALANCE_CANDIDATES_2026-07-17.md`),
the OPH capacity closure places the record capacity and the cosmological
constant in the band

```
N           in [3.30836, 3.31137] x 10^122
Lambda*l_P^2 in [2.84618, 2.84878] x 10^-122
```

evaluated at the certified enclosure midpoints
P_fwd = 1.63097209585889737696451390350695562985390,
alpha_U = 0.041124247441816685140889933889659717292128290516, with the three
candidate values

```
seed 15pi/16:        N = 3.311373322e122,  Lambda*l_P^2 = 2.846184059e-122
(pi/24)a_U per port: N = 3.311176978e122,  Lambda*l_P^2 = 2.846352830e-122
pi/48:               N = 3.308356957e122,  Lambda*l_P^2 = 2.848779042e-122
```

The uncorrected bridge value N = 3.532131543e122 is not in the band and is
carried as the null alternative.

## Decision policy (to freeze with the target)

- Reference observable: Omega_Lambda*h^2 posterior from the frozen likelihood
  combination Planck 2018 TT,TE,EE+lowE+lensing (the combination the corpus
  consumes; the +BAO sensitivity is recorded in
  `code/capacity_readback/planck_posterior/`). Any change of combination
  requires a new registration.
- Verdicts at future posteriors, z computed in ln N against the band edges:
  - band excluded at z > 3 with posterior width sigma_lnN <= 0.01: the
    corrected-balance family dies; the residual burden returns to the balance
    derivation and the Lambda side.
  - uncorrected bridge excluded at z > 3 while the band survives: the
    correction term is established as physical, pending its derivation.
  - both excluded: the capacity lane dies as formulated.
- No verdict is claimed at today's precision (band sits at 0.03 sigma; the
  test begins at sigma_lnN <= 0.01).

## Integrity

sha256 of the target statement block (the fenced band values above, exact
bytes between the first and second fence): to be computed at anchoring time
over the adopted file; this draft records the values to full stated digits so
the adopted file can be reproduced independently from
`proof/nclosure/scan.py` and the certified enclosures.
