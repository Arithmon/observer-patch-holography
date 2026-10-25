# OPH Paper Index

This directory is the canonical publication surface for the main OPH papers. Each paper is kept as a TeX source beside its built PDF.

The papers distinguish a proposed complete-history interpretation from finite
repair theorems and conditional physical reconstruction. Internal computation
and timeless closure describe the same proposed structure; they do not equate
a solver counter, causal precedence, modular parameter or local clock reading.
The [working terms](../docs/INTERPRETATION_AND_TERMS.md) summarize these
definitions. Mathematical background supplied by a benchmark remains an
explicit premise even when the programme aims to reconstruct its physical role.

The [paper set](../docs/PAPER_SET.md) records the model milestone around which
these papers are organized and the role of each. A new paper release follows
the owner's approval under the
[release policy](../model/CONSTRUCTION.md#research-programme-and-paper-release-hold);
audits, correctness fixes, dated errata and review builds continue between
releases, and preview CI does not authorize publication. The
[basis](../docs/BASIS.md) states the closure hypothesis, the no-outside rule,
the agreement condition and the executable specification that the papers share.

For a first reading, begin with the flagship synthesis,
[*Observers Who Read Their Universe into Existence*](../flagship/from_observer_consensus_to_standard_physics.pdf).
It gives the account, develops selected results and states the closure
hypothesis. The construction paper,
[*An Executable Observer for Observer Patch Holography*](executable_observer_construction.pdf),
states the specification, its acceptance criteria, its certified components
and the premise map. The remaining papers form the theorem library: *Observers
Are All You Need* for observers, records and the closure equations, the two
derivation papers for spacetime and gravity on one side and Standard Model
gauge structure on the other, and the consensus, particle and screen papers
for the finite repair theorems, the matter sector and the carrier
architecture.

| Paper | Role |
| --- | --- |
| [Observers Who Read Their Universe into Existence](../flagship/from_observer_consensus_to_standard_physics.pdf) ([source](../flagship/from_observer_consensus_to_standard_physics.tex)) | Flagship synthesis: the account of the theory, selected results, the construction scoreboard and the closure hypothesis; released on [PhilPapers](https://philpapers.org/rec/MUEFOC) |
| [An Executable Observer for Observer Patch Holography](executable_observer_construction.pdf) ([source](executable_observer_construction.tex)) | Construction paper: specification, acceptance criteria, certified components, common field operators, premise map and open connections |
| [Reality as a Consensus Protocol](reality_as_consensus_protocol.pdf) ([source](reality_as_consensus_protocol.tex)) | Theorem library: finite repair, protected records, quotient normal forms and provenance order |
| [Observers Are All You Need](observers_are_all_you_need.pdf) ([source](observers_are_all_you_need.tex)) | Theorem library: observers, records, quantum identities, thermodynamics, matter continuum and the closure equations |
| [Recovering Observer Spacetime and Einstein Dynamics](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) ([source](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)) | Theorem library: provenance compiler, source-derived causal order, ambient 1+3 carrier and the conditional continuum and Einstein chain, common metric dictionary and the static-product stress obstruction |
| [Deriving Standard Model Gauge Structure](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) ([source](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex)) | Theorem library: categorical and finite-carrier gauge routes, conditional response-derived matter image, spin-lift premises and global-form boundary |
| [Deriving the Particle Zoo](deriving_the_particle_zoo_from_observer_consistency.pdf) ([source](deriving_the_particle_zoo_from_observer_consistency.tex)) | Theorem library: internal matter decoding, one-Higgs interactions, unselected family coefficients, hierarchy coordinates and the executable pipeline |
| [Federated Echosahedral Screen Microphysics](screen_microphysics_and_observer_synchronization.pdf) ([source](screen_microphysics_and_observer_synchronization.tex)) | Theorem library: twelve-port source selector, central records and observer synchronization |

The flavor theorem is developed in [The Positive-Chamber Koide Identity for
Icosahedral Face Circulants](../extra/koide_identity_from_positive_c3_face_circulants.pdf),
and the finite de Sitter theorem in [The de Sitter Time-Advance Sign from a
Finite Screen with Fixed Capacity](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf).
Focused papers are listed in the [supplement index](../extra/). Papers that
left the active set keep their sources and last built PDFs in
[`archive/`](archive/); the [paper set](../docs/PAPER_SET.md) lists them.

Shared TeX lives in [`tex_fragments/`](tex_fragments/). The [BFT/QECC appendix](appendix_B_bft_qecc_extensions.tex) is included by *Reality as a Consensus Protocol*.

The [fine-structure supplement](../extra/fine_structure_constant_derivation.pdf)
retains the certified trial-map roots and separately reconstructs the
[entangled-qudit candidate of A. D. Wissner-Gross](https://alexwg.org/publications/ChainTheoryAlpha.pdf).
His exchange construction and exact scalar receive direct attribution.
The OPH complete-outcome, support-alignment and normalized-readout controls
retain the independent physical coupling and same-cell premises.

## Reproducibility

The shared release identifier lives in [`release_info.tex`](release_info.tex).
[`paper_release_manifest.json`](paper_release_manifest.json) records the paper
artifacts and the builder-stamped canonical book receipt.

From the repository root:

```bash
python3 tools/refresh_paper_release.py --preview
```

This rebuilds review PDFs and their manifest without requiring a version bump,
including when the visible identifier belongs to an existing release. After
corrections, main-branch PDFs can therefore be current-source previews with
the previous visible identifier; tagged or deposited release artifacts retain
their original bytes. After
an explicit maintainer decision under the release policy and manual review,
run `python3 tools/bump_paper_release.py` followed by
`python3 tools/refresh_paper_release.py --publication` to prepare a candidate
under a new release identifier. Publication mode rebuilds the canonical book
and binds its bytes to that identifier in the same manifest.
