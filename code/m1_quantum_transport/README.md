# Charged quantum transport

The main result completely characterizes the **untilted first-order Weyl
velocity matrices** of finite native register flights with onsite scattering:
the sum of their singular values is at most the native record-flight speed c.
Every such matrix is attainable. In three isotropic dimensions the sharp
field speed is c/3. Saturation forces a weighted spherical two-design;
the unique minimal solution has four equally weighted tetrahedral directions.
An explicit four-channel unitary realizes it as one charged flight.
Four counts directions across the full word, and internal channels only
for a single flight. Multi-stage words can reuse two internal channels.
Its full four-channel dynamics has a uniform native-time error bound, retaining
the alternating high-band phase at odd ticks. The very same walk sends exact
occupation records at speed c from explicit finite preparations, so the
factor-three difference is an observable property of the complete process.

A broader theorem needs no flight decomposition: the velocity ellipsoid of
any translation-invariant finite-range unitary lies inside its displacement
polytope. With N nonzero displacements an isotropic field obeys
N(1-v/c)>=2 for v>0. The tilted extension requires a strictly positive
centered radius; stationary and purely translated bands are not excluded.
Fixed finite stencils cannot even approach the full round-cone
speed under refinement. A growing direction set is necessary; it is not a
sufficiency claim. The stronger c/3 bound remains for every timed flight
lowering, regardless of how its directions or internal dimension grow.

The explicit three-flight walk attains c/3. Its entire spectrum has eight
Weyl cones, four of each chirality. The derivation proves an explicit dynamic
error bound, global spectral coercivity and the full fermionic thermal limit;
it never discards extra species. A separate finite-site theorem excludes exact
positive propagation delays for interacting bounded Hamiltonians. Covariant
noncommuting Pauli velocities give a standard continuum escape with speed c,
and identify the operational structure a microscopic realization must supply.

Read [the contract](CONTRACT.md) and [the proof](DERIVATION.md). This is a
sharp obstruction and optimal construction for a declared operation class.
It does not certify a complete twelve-port A1--A3 source, adopt a new axiom,
derive a physical energy calibration, or identify gauge SU(2) with spin.
The finite reference is a declared quasienergy vacuum, fixed throughout.

## Reproduce

From the repository root, with `requirements.txt` installed:

```sh
python code/m1_quantum_transport/build.py
python code/m1_quantum_transport/verify.py
python -m pytest -q code/m1_quantum_transport
cd Lean
lake build Geometry.M1QuantumTransportAxiomAudit
```

The compact receipt is regenerated only after independent reconstruction.
The producer uses symbolic matrix products and explicit occupied-bit gates.
The checker imports neither the producer nor SymPy: it uses its own rational
quadratic field, path enumeration, closed spectral formula and exterior-power
minors. It reconstructs every catalog before comparing exact JSON types and
values; hashes bind regenerated full matrices and histories, not trusted inputs.

Coverage includes all six flight orders with 84 intermediate branch events,
576 exact momentum matrices, all eight zero and eight pi points, 6,400 momenta
at two inverse temperatures, and 688 basis-stage executions over complete
16- and 64-dimensional fermionic Fock spaces. The latter are finite gate
demonstrations of the lift, not a simulation of a macroscopic three-torus.
Each final state is also run through the actual reversed instructions,
adding 688 independently reconstructed inverse stages. Reversed trajectories,
conjugated onsite gates and the return of every basis state are checked,
in addition to matrix unitarity and all number sectors. Four
velocity resolutions, all 144 two-qubit Pauli probe pairs, three actual
remote-preparation responses and 80 coherent spinors are also retained.
The minimal tetrahedral process adds 64 exact four-channel momentum matrices,
6,400 complete thermal momenta, both rank-two bands, all sixteen onsite Fock
inputs, four lossless native record experiments and twenty trajectory samples.
The general-stencil evidence checks
1,664 exact directional Loewner inequalities by differentiating the actual
Fourier kernel independently of the producer's gate derivatives. Four support
polytopes are reconstructed by independent linear-system and supporting-plane
methods, including the tetrahedron that attains both speed bounds.
Three exact scope controls differentiate stationary, purely translated and
positively centered tilted unitaries, rejecting a false zero-radius cap bound.

Hostile controls cover retiming, lost histories, hidden species, changed
phases despite identical probabilities, regenerated signless many-body
programs, invalid positive effects, superluminal velocities, stale source
custody and missing catalogs. Selected forgeries run through the optimized
Python CLI, where assertions cannot enforce acceptance. Dynamics and the
coercive bound are tested on all cones and full shifted momentum grids.

Thirteen Lean lemmas certify finite algebra reductions, with a transitive axiom
audit that rejects `sorryAx` and compiler-trust declarations. The compressed
operator resolution, SVD characterization, analytic rigidity, thermal limit
and continuum field arguments are analytic proofs, not fully formalized
operator theorems. Finite numerical checks support those proofs and do not
replace their all-level quantifiers.
