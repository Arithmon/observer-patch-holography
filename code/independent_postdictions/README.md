# Independent postdictions (#1011)

Fresh numerical implementations reproduce the two certified alpha/P roots and
the balanced, ordered tau interval. The selected prospective target is canonical
FZ-10; readiness is **NOT_READY**. No held-out physical measurement is evaluated.

See [the derivations and exit audit](../../extra/INDEPENDENT_POSTDICTION_COMPARISON.md),
[execution contract](CONTRACT.md), and [future-release protocol](protocol.json).

From the repository root, with Python 3.12 and `requirements.txt` installed:

```sh
PYTHONPATH=code python -m independent_postdictions.verify
PYTHONPATH=code python -m pytest -q code/independent_postdictions
PYTHONPATH=code python -m independent_postdictions.build
```

In PowerShell set `$env:PYTHONPATH='code'` first and omit the shell prefix.
`build` explicitly regenerates only the new packet's `receipt.json`; verification
never rewrites evidence. There are no network calls in build, verify or tests.
The fresh solver imports neither historical producer nor its numbers. The
checker independently sums representations and integrates the one-loop kernel;
tau uses the original implicit Q equation versus the checker's quadratic.

Expected inverse-alpha roots are `136.99483517741294` and `137.03566013694658`;
the separate converged asymptotic value is `136.99402589611851`. Conditional tau
is `1776.96902729315735 MeV`, with the historical outward box
`[1776.968991,1776.969063] MeV`. The checked receipt also retains the historical
approximate trunk and the explicitly measured-alpha-consuming mixed diagnostic.

The alpha mathematical enclosures are not physical theory errors. Tau's input
box is not a joint confidence interval, the balance premise has measured-triple
ancestry, and the result is identical to the Koide reference prediction. No
chance probability or unique OPH confirmation follows. Historical PDG data
produce FZ-10's INCONCLUSIVE verdict. A window-robust sidecar preserves the
difference between a rounded-center verdict and exclusion of the whole window.

`admission.nominate` consumes metadata only and never skips an exposed first
candidate for a fresh second one. It cannot authenticate attestations or admit
a physical dataset. `comparison.compare` implements arithmetic only; callers
must first meet the external admission contract. It does not mark input JSON
as a measurement or bypass NOT_READY. All controls in tests are synthetic.

The decimal interface accepts magnitudes up to `1e12`, at most 80 fractional
places, and standard uncertainties at least `1e-30 MeV`. Within that arithmetic
domain the 110-digit context makes threshold arithmetic exact; rounded division
is display-only. Precision, rounding, exponent limits and traps are isolated
from caller settings. An independent rational interval oracle checks 297
scheduled center/window boundary cases. Coarse uncertainties remain valid inputs, including above
`100 MeV`; they are not a selection filter. An unsupported numerical format or
range needs a reviewed implementation update before evaluation, never silent
rounding or substitution of another dataset.

Optional metadata discovery, **only when separately needed**, writes a new file:

```sh
PYTHONPATH=code python -m independent_postdictions.fetch_catalogue new_snapshot.json
```

The checked snapshot includes unsuccessful literal-tau queries, the broader
query, and a historical positive control. Its conclusion is bounded discovery,
not exhaustive absence of measurements. The contract, historical inputs and
exact classified catalogue snapshot have separate review pins; blindly
rebuilding custody cannot loosen them. The inherited interval certificate also
has an immutable digest: its bounds cannot be widened during a packet rebuild.
These controls bind the audited material, not the truth of an external index or
an experiment's exposure declaration. Metadata titles can reveal outcomes;
unexpected exposure must be recorded and blocks the first candidate's admission.
Historical freezes, canonical prediction rows and OTS proofs are unchanged.
