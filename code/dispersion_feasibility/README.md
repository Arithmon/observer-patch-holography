# Fixed-scale dispersion feasibility

The photon-only and shared photon/electron/positron hypotheses have different
pair-production thresholds. A direction-uniform analytic bound controls their
full three-dimensional thresholds; high-precision collinear witnesses check the
arithmetic independently. The calculation supplies no photon flux or empirical
rejection. Its physical scale, clock, field attachment and massive dispersion
are explicit hypotheses.

The [decision contribution and proof](../../docs/research/DISPERSION_FEASIBILITY.md)
state the inputs, exact domain and implication for issue #1025.

From the repository root with the pinned requirements and `PYTHONPATH=code`:

```text
python -m dispersion_feasibility.bounds
python -m dispersion_feasibility.kinematics
python -m dispersion_feasibility.check
python -m pytest -q code/dispersion_feasibility
```

`kinematics` prints the complete report. `check` reads the retained report or an
explicit path. It imports no producer, reconstructs the thirty directions by a
different formula, and verifies energy conservation, stationarity, positive
collinear curvature and the analytic bracket. Numerical witnesses are not
global optimization certificates; the proof supplies the global bound.

This contribution occupies one bounded feasibility session. Numerical work is
limited to 15 cumulative local CPU minutes, with no cloud run, flux transport,
future-data access or Lean expansion. It does not assert completion of #1025,
activate #1026 or amend the immutable FZ-12 registration. The P value comes
from the exposed alpha calibration; the length conversion and electron mass
are declared nominal inputs without a propagated calibration uncertainty.
