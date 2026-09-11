# OPH Cosmology And Dark Gravity

This directory contains unpublished research on the cosmological continuation of Observer Patch Holography: modular-charge dark gravity, finite-source primordial structure, observer-screen synchronization, Boltzmann transport, and comparison with public cosmological data.

The main papers provide finite observer and matter structures, controlled real-field continuum trajectories on a supplied domain, an interacting Hilbert space, and an internal action clock. These companions develop the source laws, expanding geometry, transfer maps, calibration, and likelihood contracts needed to turn such objects into cosmological observables.

## Paper Map

| Paper | Main contribution |
| --- | --- |
| [Dark Gravity](oph_dark_matter_paper.pdf) ([source](oph_dark_matter_paper.tex)) | Modular collar charge, quantitative recovery saturation, conditional deep-galaxy law, and reproducible galaxy comparisons |
| [Finite-Source CMB Program](oph_cosmology_finite_source_cmb_program.pdf) ([source](oph_cosmology_finite_source_cmb_program.tex)) | Conditional source-screen spectrum, finite receipt contract, transfer requirements, and CMB promotion path |
| [Inflation Without an Inflaton](oph_inflation_without_inflaton_observer_screen_synchronization.pdf) ([source](oph_inflation_without_inflaton_observer_screen_synchronization.tex)) | Observer-screen synchronization, conditional edge-center tilt, one-shell radial obstruction, physical dilation and tomography routes, flatness, and horizon coherence |
| [Cosmological Vacuum And Structure Formation](oph_cosmological_vacuum_and_structure_formation.pdf) ([source](oph_cosmological_vacuum_and_structure_formation.tex)) | Vacuum boundary, fluctuation ensembles, proto-objects, worldlines, and structure seeds |
| [Cosmology Data And Likelihood Contracts](oph_cosmology_data_likelihood_contracts.pdf) ([source](oph_cosmology_data_likelihood_contracts.tex)) | Boundary between the conditional source theorem, finite source instantiation, transfer, nuisance treatment, and official likelihood comparison |
| [Boltzmann Transport Derivation](oph_boltzmann_transport_derivation.pdf) ([source](oph_boltzmann_transport_derivation.tex)) | Finite transport interface between OPH sources and observable distribution functions |
| [Black-Hole Information Ledger](oph_black_hole_information_ledger.pdf) ([source](oph_black_hole_information_ledger.tex)) | Proved finite unitary dynamics and sign-unitary record channels, with a conditional entropy and thermality formulation |

The [formal radial-lift theorem fragment](../paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex)
contains the one-shell non-identifiability proof, physical source-dilation
theorem, radial cross-covariance tomography, exact amplitude conversion, and
finite-window bound. The
[simulator contract](../code/cosmology/SIMULATOR_RADIAL_CONTRACT_330.md)
specifies the fail-closed evidence split used by finite runs.

The [physical CMB theorem program](physical_cmb_theorem_program.md) specifies the supplied and unsupplied physical-event, order-faithful-placement, count--volume, continuum, source, lift, stress, abundance, transfer, and likelihood inputs in one place.

## Expansion, Area Law, and Public Data

On the golden source-record family of the spacetime construction, the expansion history is a record density: the comoving read law generates the same finite causal order for every scale-factor profile, the scale factor enters only as the count weight `a^4` per comoving cell and conformal tick, the count clock reads proper time, and the redshift between two epochs is the fourth root of a record-count ratio ([fragment](../paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex), [Lean](../Lean/Geometry/SourceNetConformalRecordDensity.lean)). The number of unordered read pairs crossing a surface per layer is an exact area law, `(pi/4) n^2 a^4` per unit area, with the Planck-area dictionary `n^2 a^4 = 1/(pi l_P^2)` ([fragment](../paper/tex_fragments/CROSSING_READ_AREA_LAW.tex), [Lean](../Lean/Geometry/CrossingReadAreaLaw.lean), [receipt](../code/causal_refinement/crossing_read_area_law_receipt.json)).

The conditional primordial, background, capacity and dark-sector values are compared with public measurements in the [cosmology postdiction ledger](../code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md), produced deterministically and checked by an independent verifier, with its brackets machine-checked in [Lean](../Lean/ObserverPatchHolography/EinsteinBranch/CosmologyLedgerBrackets.lean); every row is on seen data and carries its class. Verbatim-ready freeze proposals for the tilt, the tensor zero, and the running and isocurvature zeros are in the [primordial registration proposal](../code/cosmology/primordial_freeze/PRIMORDIAL_REGISTRATION_PROPOSAL.md); every numeric kill cell is an owner slot.

## Dark-Gravity Structure

Under the collar-stress attachment, recovery defects bound a candidate dark source. Conserved comoving charge gives the familiar inverse-volume dilution. Two explicit scaling and composition premises select a spherical deep-acceleration law with baryonic Tully–Fisher scaling.

The physical interpretation links this source to relativistic stress, abundance, lensing, cosmic expansion, and calibrated readouts on the same branch. The papers state the hypotheses and compare the resulting formulas with galaxy data.

## Reproducibility

Companion radial-lift and clock-certificate code lives in
[`../code/cosmology/`](../code/cosmology/). The larger simulator and
visualization surfaces are:

- [OPH physics simulator](https://github.com/muellerberndt/oph-physics-sim)
- [Interactive simulation](https://simulation.floatingpragma.io)

The cosmology papers are focused research companions and are not part of the core release bundle.
