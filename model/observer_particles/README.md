# Particles and observer fields

Exact response operators recover one generation's internal charge
multiplets under the declared spin lift. One Higgs doublet admits all three
charged-fermion interaction channels. Conditional scalar, fermion and vector
constructions share explicit spatial operators or metric dictionaries.
These finite witnesses connect the representation and field descriptions;
physical species, masses and a responsive gravitational geometry require
their separate selection and observation tests.

The entangled-pair exchange construction and its fine-structure candidate
are due to [A. D. Wissner-Gross](https://alexwg.org/publications/ChainTheoryAlpha.pdf),
*A Parameter-Free Candidate for the Fine-Structure Constant from Entangled
Qudits* (2026). This package independently implements that finite construction
and adds OPH comparison, complete-outcome and normalization controls. Its
physical photon identification is a separate conjecture.

Read the [mathematical contracts](DERIVATION.md), the
[particle paper](../../paper/deriving_the_particle_zoo_from_observer_consistency.tex)
and the [verification receipt](results/verification.json).

From the repository root:

```sh
python -B model/observer_particles/run.py --verify-only
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B model/observer_particles/run.py
python -m pytest -q model/observer_particles
```

Python requires NumPy, SciPy, SymPy and mpmath. The second command regenerates
this directory's receipts and operator artifacts. Each receipt binds its
producer and computational inputs; exact algebra and floating diagnostics
have distinct scopes. The accompanying Lean module is
`Screen/ObserverFieldInterfaces.lean`.
Run `lake build ObserverFieldInterfaces` from `Lean/` for its exact interface
proofs. The larger representation and interaction producers retain their
separate executable-algebra receipts.
