# Physical CMB Theorem Boundary

This unpublished research interface specifies the scientific conditions for a physical CMB claim: a bounded finite observer-patch source on an inhabited physical source-causal refinement family, a covariant stress tensor, a physical scale map, gauge-invariant initial data, an Einstein--Boltzmann transfer system, and a complete observable likelihood. OPH-FPE CMB curves remain diagnostics.

A declared conservative-record population with complete local neighbour reads and a layer-duration law has a controlled flat causal-order and normalized count-volume limit. For interior timelike diamonds and a nondegenerate reference interval in the same history, complete authenticated ancestry and anchors give a fourth-root count ratio converging to the proper-time ratio without coordinates, timestamps or density in the readout. Its reference interval fixes relative units. The separate Lorentz-cone selection theorem requires a nonzero closed convex pointed cone, operational covariance under the proper icosahedral rotations and every rapidity of one boost axis, and a time orientation. These hypotheses concern histories, interventions and clock readouts; an algebraic source frame does not supply them.

On the same family the spatially flat expanding background has an exact finite form. The layered read law in comoving coordinates is conformally invariant, so the flat order is the causal order of every spatially flat expanding history; the scale factor enters only as the record-density weight `a^4` per comoving cell and conformal tick that a constant physical event density assigns; the count clock reads proper time up to the ratio of the largest and smallest scale factor on the diamond; and the redshift between two epochs is the fourth root of the ratio of record counts of congruent comoving diamonds. These statements are proved in `paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex` and machine-checked in `Lean/Geometry/SourceNetConformalRecordDensity.lean`. The record-production profile is a supplied datum on this family; the Friedmann equation is, in this language, a law for that profile, and no statement derives it.

On a supplied fixed spatial cone, the charged scalar action supports controlled real-sector continuum trajectories, a neutral interacting Hilbert space with explicit initial states, and a Jacobi clock on regular nonturning solutions. On the golden source population itself, a positive source-metric scalar action has a finite-time continuum limit with an explicit detector bound, a canonical quantization with a unique ground state, and a certified joint spatial and temporal refinement of one compact detector. For a massless scalar the limiting characteristics are the null cone, and the high-frequency limit is the null-geodesic Liouville flow, which supplies the free-flow generator of the transport paper on the flat family for a scalar record field. These constructions use declared geometry, action parameters and energy; cosmological source selection, expanding geometry and physical calibration require the conditions below.

## Supplied constructions

- **Flat causal order and count law** on the declared golden-record family, with the retrospective count clock and its finite enclosure.
- **Spatially flat expanding background as record density** on the same family: conformal invariance of the order, the `a^4` count weight, the proper-time reading of the count clock, and the redshift as a count ratio. The profile is supplied.
- **Scalar continuum and free-flow generator** on the flat family, with the null cone identified with the read-law cone under the common-model condition that the action's wave speed is the read-law ratio.
- **Seam-count area law**: on the layered family the number of unordered read pairs crossing a surface per layer converges to `(pi/4) n^2 a^4` per unit area, an exact area law whose Planck-area dictionary is `n^2 a^4 = 1/(pi l_P^2)`; see `paper/tex_fragments/CROSSING_READ_AREA_LAW.tex` and `Lean/Geometry/CrossingReadAreaLaw.lean`.
- **Conditional primordial values**: `n_s = 1 - P_star/48` under the reserve-generator receipt, `r = 0`, zero running and zero source-side isocurvature under the rank-one single-clock branch and the dilation cocycle, `Omega_K = 0` under the small-loop holonomy identification, and `(w_0, w_a) = (-1, 0)` on the fixed-capacity branch.

## Required constructions

- **Finite source:** source-derived scalar amplitude, shape, clock, release surface, and background capacity, including the reserve-generator receipt that the tilt theorem consumes.
- **Record-production profile:** the law that selects how many records a comoving cell produces per conformal tick, which is the expansion history on the sourced family; the flat family supplies the form of the count measure and no profile.
- **Physical causal carrier:** physically identified events and links on a compatible refinement, with exact order reflection or a uniformly vanishing causal discrepancy, calibrated count-to-volume convergence, and compatible source directions. The conservative-record flat limit uses a specified population, local read law and layer duration. Its ancestry-count clock reconstructs relative geometric duration; selecting the law and identifying a common physical clock remain additional.
- **Continuum promotion:** independently calibrated count-to-volume density, independent dimension and manifoldlikeness tests, stable topology, and a unique Lorentzian limit. Smooth Einstein promotion additionally needs either same-family tensor-curvature convergence or the continuum small-ball/null-balance identification; scalar-curvature convergence alone is diagnostic.
- **Physical scale:** maps to physical wavenumber, angular multipole, scale factor or redshift, and a common mode basis.
- **Covariant stress:** total stress-energy closure, recipient stress for nonzero exchange, causal response, and refinement convergence.
- **Dark-sector kernels:** abundance, pressure, sound speed, anisotropic stress, exchange current, and modular collar response from one covariant source/action attachment. Conserved comoving charge alone supplies only homogeneous inverse-volume dilution, which on the sourced family is the layer-independence of the charge per comoving cell.
- **Transfer:** regular gauge-invariant initial modes, recombination and collision terms, the two-polarization photon sector and the record-process clauses of the refinement Liouville theorem, numerical convergence, and recovery of declared limiting cases.
- **Likelihood:** datasets, masks, covariance, nuisance parameters, priors, and combination rules.

The radial mathematics proves the one-shell obstruction and two sufficient
conditional uniqueness routes: physical dilation and cross-covariance tomography. No finite source map
that passes either radial route is supplied, nor is an inhabited physical
source-causal continuum family with the required physical order and volume
approximation, manifoldlikeness, topology, and tensor-curvature or
small-ball/null-balance certificate. The
relativistic repair stress, dark abundance, Boltzmann bridge, and joint
likelihood remain unsupplied. CMB outputs are diagnostics and do not establish or falsify a
cosmological OPH claim.

## Comparison with public data

The conditional primordial, background and dark-sector values are compared with public measurements in the retrospective ledger `code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md`, generated by its producer and checked by an independent verifier, with the exact rational brackets of its rows machine-checked in `Lean/ObserverPatchHolography/EinsteinBranch/CosmologyLedgerBrackets.lean` and the cocycle-to-tilt step behind `n_s = 1 - P_star/48` in `Lean/ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltCocycle.lean`. Every comparison there is on seen data and carries a class label: conditional-theorem postdiction, shared baseline, closure-candidate display, fitted comparison value, or not evaluable. Rows that distinguish the branch from the standard baseline are the exact tensor zero, the exact zero running, the specific tilt value and the fixed-capacity equation of state. The proposal `code/cosmology/primordial_freeze/PRIMORDIAL_REGISTRATION_PROPOSAL.md` states verbatim-ready register rows for the tilt, the tensor zero, and the running and isocurvature zeros, with every numeric kill cell an owner slot; the fixed-capacity row FZ-13 and the ringdown row FZ-14 keep their pending status.

The detailed interfaces live in:

- `cosmology/oph_cosmology_finite_source_cmb_program.tex`
- `cosmology/oph_cosmology_data_likelihood_contracts.tex`
- `cosmology/oph_boltzmann_transport_derivation.tex`
- `paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex`
- `paper/tex_fragments/CROSSING_READ_AREA_LAW.tex`
- `paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex`
- `paper/tex_fragments/PHYSICAL_SCALE_BRIDGE_THEOREMS.tex`
- `paper/tex_fragments/FINITE_COVARIANT_PARENT_THEOREMS.tex`
- `paper/tex_fragments/WHITNEY_REAL_CONTINUUM.tex`
- `paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex`
- `paper/tex_fragments/WHITNEY_EPHEMERIS_CLOCK.tex`
- `paper/tex_fragments/SOURCE_NET_CAUSAL_LIMIT.tex`
- `paper/tex_fragments/SOURCE_COUNT_CLOCK.tex`
- `paper/tex_fragments/SOURCE_METRIC_SCALAR_CONTINUUM.tex`
- `paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex`
