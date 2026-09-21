# Finite Maxwell pressure and classical Gibbs exports

This package measures canonical metric work and integrated field stress for
the supplied Whitney Maxwell action on the existing twelve-port tetrahedral
cone. It also evaluates a separately supplied classical Gibbs ensemble of
the same finite action. It retains all 30 physical cone modes; these differ
from the 19 curl modes of the boundary-only carrier.

The geometry, Maxwell action, constitutive coefficients, canonical state
interpretation, initial data and model time are inputs. The Gibbs branch
additionally supplies the ensemble and temperature (`k_B=1`). Neither branch
derives these inputs from pair-average repair. A finite classical cutoff is
not a quantum blackbody or continuum thermodynamic limit.

The original source records belong to bounded self-reading computational
patches with local state, ports, readback, retained records and feedback.
The new pressure readout authenticates and consumes their decoded potentials.
It reports **field energy only**, excluding source, clock, observer memory
and material boundary energies. The fresh source-free execution retains full
states but is not itself a new observer-instrument implementation.

From the RER repository root, with Python, NumPy, SciPy, SymPy and pytest:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 code/equation_of_state/maxwell/maxwell_eos.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 code/equation_of_state/maxwell/verify_maxwell_eos.py --output code/equation_of_state/maxwell/runs/verification.json
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m pytest -q code/equation_of_state/maxwell/test_maxwell_eos.py
```

`runs/maxwell_eos_receipt.json` contains the source and implementation hashes,
oriented mesh, complete field forms, 81 canonical states, all 80 Maxwell and
work checks, metric-deformation stencils, original decoded-field readouts,
controls, and the Gibbs covariances and partition functions. The neighboring
CSV files expose trajectory diagnostics and the nine Gibbs parameter cases.
The JSON is authoritative; CSV files are convenience exports.

The independent verifier uses positive simplex quadrature instead of the
producer's analytic moment assembly. For the Gibbs calculation it uses a
different canonical spanning-tree/cotree gauge, retaining the Liouville
measure, and checks energy, all stress components and partition-function
volume work. The implementation is binary64, checked to numerical
tolerances, not a rigorous interval enclosure.

`mean_pressure` is one third of the stress trace. Directional stresses may
include tension. The Gibbs branch separately checks isotropy of its prepared
ensemble; it does not show thermalization of the dynamical trace. The ratio
is undefined at zero field energy. Adding other energy sectors can change
the total pressure/energy ratio, as the explicit comparison controls show.

The metric-scaling derivation and its held-fixed variables are stated in
the producer's module and Gibbs-function docstrings. Original sourced slab
readouts use instantaneous field energy at the declared interpolation
midpoint; they must not be confused with the original integrator's modified
conserved energy.
