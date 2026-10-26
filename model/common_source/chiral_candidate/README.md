# A source-executed chiral matter candidate

This candidate joins the existing OPH internal matter representation to an
all-left Weyl kinetic symbol, its verified one-Higgs interaction tensors,
and the native nonscalar family feature. A continuing source experiment
reads three charged rest gaps against quantum clocks. It retains the full
90-mode Gaussian field and independently executes the probed Fock factor
with its clock purification references.

The representation and chiral continuum symbol have the minimal Standard
Model fermion types. **The trial masses fail the observed hierarchy.** The
finite spectral menu also supplies no local nonperturbative chiral gauge
regulator. Neither result is a physical identification of our Universe.

## The action and the choices that define it

[action.py](action.py) reads the retained public `M15` representation and
Higgs tensors. In all-left notation the internal types are
`Q,u^c,d^c,L,e^c`; the three singlets are charge conjugated. The construction
uses `2 Spin components ×15 internal components ×3 native family axes`.
Every component remains, including the neutrinos.

The chosen trial assigns

```text
Y_up = a F,  Y_down = b F,  Y_charged_lepton = J F,
a=log(25/8), b=log(8/5), J=log(125/64), Omega=log5,
H=(0,v), v=Omega, c=Omega/2.
```

`F` is the exact native family feature computed in
[particle_identification](../particle_identification/SM_MATCH.md).
Its eigenvalues are `61/450` and `41/225±7 sqrt(5)/750`.
The choices of Yukawa map, fixed neutral scalar background, coordinate speed
and supplied spatial attachment are explicit architecture prescriptions.
They are legitimate candidates to test and contain no empirical mass fit.
They are not consequences of the representation theorem. The scalar is
not dynamically generated or quantized in this trial.

The Higgs vertices give a **symmetric internal all-left mass kernel** `K45`.
The creation pairing is `Delta90=epsilon_spin tensor K45` for this real
background. Its transpose is antisymmetric, as CAR requires. This uses
the gauge-compatible Higgs lines rather than an arbitrary bare mass between
unrelated chiral representations. The producer rechecks the current/Higgs
identities and anomaly traces, and preserves model electric charge in the
neutral background. [ACTION.md](ACTION.md) states the conventions and exact
hierarchy bounds.

## Chirality, dispersion and the finite cutoff

[kinetic.py](kinetic.py) uses `h(p)=-c sigma.p tensor I45` and the
opposite-momentum Nambu block

```text
B(p) = [[h(p), Delta], [Delta†, -h(-p)^T]].
```

Spin conjugation of the holes gives the standard massive chiral symbol;
its squared energies are `c²|p|²+s_j(K)²`. Color and electric charge,
Spin/helicity, the current action and Nambu redundancy are checked on the
same operators. Counting colors and families gives 21 charged Dirac fields
and three massless left-Weyl neutrinos. The neutrinos have no right-handed
partner in `M15`. Nambu zero eigenvalues do not add sterile species.
This minimal neutrino branch does not reproduce measured neutrino masses.

The spatial symbol is a declared continuum/spectral construction. It is
not the previous eight-component vectorlike regulator relabeled as chiral
matter. The naive lattice control retains eight Weyl nodes with balanced
handedness. Removing them from a plot would not remove them from the law.
See [KINETIC.md](KINETIC.md) and the
[regulator routes](REGULATOR_ROUTES.md) for the exact domain and relevant
primary literature. A local interacting chiral regulator, dynamical gauge
fields and a physical Lorentzian reconstruction remain separate work.
The [frozen-metric attachment](spacetime_attachment.py) now verifies an
exact native Clifford/chiral basis change, tetrad pullback, common scalar
characteristic cone and lapse/clock conversion. It also records the
retained incidence zero modes and radial boundary obstruction. This is a
local principal-symbol bridge, not a complete common field history.
The [shareable research note](../../../docs/research/OPH_CHIRAL_WEYL_HIGGS_CANDIDATE_2026-10-26.md)
develops the derivation and these limits in one place.

## One continuing clocked field history

[gaussian_source.py](gaussian_source.py) executes actual native two-mode
Fock words on **all four local input columns**, recovers their completed
Majorana action, and composes the finite CAR transformation. It includes
fermionic SWAP routing with its double-occupation sign and charges every
routing pulse. No `2^90` state vector is allocated. The finite Gaussian
identity supplies the many-mode extension. References and all measurement
weights remain explicit; there is no polar, covariance or probability repair.
[SOURCE.md](SOURCE.md) gives the representation proof and numerical limits.

[history.py](history.py) keeps all 42 charged pair blocks and six neutral
zero modes. An initial fermionic wire layout places three probed lepton
blocks first; no modes are discarded. At rest these blocks factor exactly.
The experiment:

1. Prepares blank field modes and two calibrated modular-reference clocks
   once, retaining both clock purifiers.
2. Advances the full field and clocks for modular parameter `tau=1/4`.
3. Commits three lepton-mode occupation reads. The first record can trigger
   a neutral phase on another probed block.
4. Continues the same field and clocks for another `tau`, reads that block
   again, and reads the two clock quadratures.

There are 64 six-bit terminal records. The full field has 16 conditional
Gaussian outputs before joining the clock records. The source experiment
also executes the six selected modes and two clocks as an explicit
reference-bearing Fock calculation. Every earlier record marginal is
preserved, and the two implementations agree with an independent continuous
target. The clock calibration remains a declared programmed modular service.

Starting from an empty pair, its occupation oscillation has frequency twice
the charged quasiparticle gap. Its measured cosine and the final clock
phase determine `m/Omega`; the decoder does not consult the Hamiltonian or
simulator interval. The three recovered ratios, in native family-axis order,
are approximately

```text
0.1080141452, 0.1359561376, 0.0907450442.
```

These are tree-level trial rest energies in the specified background and
clock, not physical mass predictions in MeV. Blank field preparation is not
a proof that repair prepares the interacting physical vacuum. The full
continuing law and its measured rest gaps are defined; uncomputed gauge and
scalar dynamics are not charged to “simulation cost.”

## What the mass test says

All three trial Yukawa matrices are proportional to the same `F`. They have
the same family ratios, no tree-level quark mixing or CP phase, and every
charged tree mass lies within an exact factor `15/4` of every other one.
The charged-lepton family spread is below `3/2`. A common scale or a family
permutation cannot turn those ratios into observed hierarchies.

[hierarchy_controls.py](hierarchy_controls.py) also tests three fixed
nonlinear texture prescriptions using the certified Hubbard threshold
function, with `t=Omega/2` and `U=a F_j,b F_j,J F_j`. This imports a
**numerical texture function**, not a three-dimensional dynamical mass-gap
theorem. No coefficient is fitted or scanned against measured masses.
Ascending ratios are approximately:

| Prescription | Lightest | Middle | Heaviest |
| --- | ---: | ---: | ---: |
| Linear `F` | 1 | 1.1903 | 1.4982 |
| Hubbard function, `U=a F_j` | 1 | 204.97 | 65,743.94 |
| Hubbard function, `U=b F_j` | 1 | 353,640.05 | 354,937,934,800 |
| Hubbard function, `U=J F_j` | 1 | 8,075.31 | 136,826,308.27 |

After computing predictions, a separate comparison checks both independent
charged-lepton ratios using broad windows derived from the
[PDG summary](https://pdg.lbl.gov/2025/download/db2024.pdf). All four fixed
prescriptions fail. The nearby middle ratio in one nonlinear control cannot
be accepted while its heaviest ratio fails. Other actions are not excluded.

The specific remaining physics is the mass-generating interaction and its
vacuum/kinetic readout. The exact internal chiral types and an available
native symmetry breaker are now supplied. These alone do not select the
successful Yukawa texture, a Higgs scale or a dressed observed particle pole.

## Replay and retained evidence

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/chiral_candidate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -u -m model.common_source.chiral_candidate.run --output /tmp/chiral-candidate.json
.venv/bin/python -m model.common_source.chiral_candidate.run --verify model/common_source/chiral_candidate/results.json
```

The [receipt](results.json) binds the source and comparison inputs. Its
[NPZ](results.npz) retains the literal probe's unnormalized work/reference
outputs and the full normalized Gaussian descriptors with their separate
unnormalized probability weights. A zero-probability conditional state has
no normalized descriptor; its outcome still belongs to the instrument.
Verification checks custody and retained state/weight consistency; a fresh
output path performs the actual replay.

Floating source-block and covariance residuals are diagnostics, not outward
physical-service bounds. The hierarchy function uses outward rational
intervals. Quantum Higgs and gauge fields, nonperturbative chiral-regulator
construction, physical clock units, observed masses and responsive geometry
are not established by this trial. All earlier retained packages and r2053
artifacts remain unchanged; this candidate is unreleased.
