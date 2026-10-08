# OPH proof and paper audit October 2026

This audit reviews the OPH proof and paper stack at base commit `1ba8a011`,
dated 2026-10-08. It corrects mathematical implications that could produce
wrong observables or predictions, along with the numerical implementations
used to support those implications. The audit covers all 23 registered paper
roots, their relevant shared fragments, the Lean inventory, and selected
proof-to-code interfaces. Coverage depth is recorded below: it is not an
independent reproof of every theorem or an empirical validation of OPH.

The principal errors concern missing hypotheses, invalid sampling arguments,
loss of small spectral values, and promotion of unverified inputs. Correcting
them narrows several conclusions. No new physical prediction is established.
Historical frozen targets and custody artifacts retain their original bytes.

## Findings that change mathematical or physical conclusions

| Finding | Counterexample or consequence | Repair and executable evidence |
| --- | --- | --- |
| Random repair blocks were substituted into a fixed-index expectation bound | A chain with block-endpoint mean distance `1/2` and block lengths one or two has calendar-time mean `2/3`, outside the claimed `1/2` tube | Keep the stochastic theorem at fixed block indices; require deterministic endpoints for its stated all-time constant. `code/consensus/test_sampling_and_record_boundaries.py` |
| Adaptive block selection need not preserve the certified kernel | Select between `U` and `1-U` after seeing a fair bit: each marginal mean is `1/2`, the selected output is always one | Require block-start measurability and the conditional-kernel identity; the same exact test module checks the obstruction |
| Marginal stationarity was mistaken for record persistence | Persistent and independently resampled fair bits have the same marginals; their two-time disagreement probabilities are zero and `1/2` | Define persistence by joint two-time disagreement on identified record labels |
| Conditioning was mistaken for subsequent partition averaging | A pure state inside a rank-two event stays pure under Lüders conditioning; normalized block averaging has purity `1/2` | Correct the consensus consumer to match the existing Lean partition-average theorem and require nonzero event probability |
| A bounded covariance was assigned a nonzero dilation exponent globally | Unitary conjugation preserves operator norm. For nonzero bounded `C`, `D_s* C D_s = exp(-theta*s) C` forces `theta=0`; `A k^-theta` is unbounded on all positive wavenumbers otherwise | Separate bounded-global convergence from a common-core/local-band criterion; require a nonzero positive covariance. `paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex` |
| Radial theorem statements omitted normalization, feasibility and inversion domains | `W=2 delta_R` defeats the shell error bound; at `r0=0`, all `j_l(k r0)` vanish for `l>=1`; a target outside `Ran(A)` has no exact constrained continuation | State positive probability-window and moment assumptions, `r0>0` and an inversion class, target feasibility, and the separate positivity requirement for a physical spectrum |
| Radial inversion erased resolvable constraints and understated residuals | Normal-equation conditioning and a unit absolute tolerance can accept lost modes at small physical spectrum scales | Integrate independently reviewed PR [#1063](https://github.com/FloatingPragma/observer-patch-holography/pull/1063), preserving Mario Poneder's authorship. Prior-whitened SVD, independent KKT controls and guarded residual evaluation remain numerical policies, not interval proofs |
| Finite-window evaluation returned false zero bounds | At `theta=1e-14`, a radius change from `1` to `1.01` was erased; weights `[1e308,1e308]` normalized to zeros after overflow | Use scaled weight normalization, `log1p`/`expm1`, and the gamma recurrence for the derivative norm. Twelve independent controls reproduce eight failures on the base source: `code/cosmology/test_radial_window_regressions.py` |
| No tensor source was promoted to zero tensor power | `h=sin(k tau)/(k tau)` is a nonzero homogeneous solution of the radiation-era tensor equation | Require zero initial tensor amplitude and velocity, no forcing, positive scalar pivot power, and linear primordial scope. Three live ledger rows lose generic slow-roll discrimination credit; their numerical comparison verdicts do not change |
| A finite tensor upper limit was said to discriminate generic slow roll | The generic alternative admits arbitrarily small positive tensor power | Restrict exclusion to specified alternatives above the achieved sensitivity; semantic mutations of the live ledger are rejected even after rehashing |
| Small nonzero Majorana masses became zero | `diag(1e-200,2e-200,3e-200)` returned zero masses; within the actual intrinsic family, `a=1+1e-8,rho=1` lost both masses `a-rho`. The downstream gap producer also underflowed squared gaps and used a dimensionful cutoff for their ratio | Replace squared normal-matrix/cubic readout with scaled direct SVD, phase and congruence checks. Preserve the independent cubic diagnostic. Compute gaps exactly on the supplied binary64 masses before one rounding, refuse nonrepresentable gaps, and read their ratio without an absolute unit cutoff. `code/particles/neutrino/test_takagi_numerics.py` |
| Changing physical units erased Majorana phases | The absolute imaginary-part cutoff selected a real shortcut at `m_star=1e-20` although the relative complex phase remained nonzero | Use the real shortcut only for exactly real inputs; retain scale and subprocess/JSON congruence controls |
| Alpha endpoint metadata falsely authorized theorem promotion | Inconsistent component/total intervals, a `not_positive` string, malformed spectra and an arbitrary certificate label passed | Distinguish `contract_satisfied` from a replayed certificate and physical promotion. Check exact necessary interval relations, supports, residues and budgets. The current backend has no certificate replay and cannot promote a contract. Tests live in `code/P_derivation/test_thomson_spectral_transport.py` and its two endpoint consumers |
| Fine-to-coarse cost control did not ensure a coarse optimum was reachable | Coarse costs `A=0,B=1`, with only a fine lift of `B`, defeat the claimed implication even at zero distortion | Require a feasible contraction and a fine lift of the coarse optimum; give the three-inequality proof. `code/particles/test_paper_proof_boundaries.py` |
| Fine-structure existence and spectral support premises were missing | `x/2` maps `(0,1)` contractively into itself without a fixed point; repeated support nodes collapse a 24-site Jacobi realization; `m=0` annihilates the resolvent difference | Require a nonempty closed interval, distinct support nodes and `m>0`. Remove the unsupported claim that register count fixes the residual scale; fixed total weight does not bound the inverse-support moment |
| A Gibbs law and a settled prior were called equivalent without a fiber-mass condition | A uniform three-state prior with normal-form fibers of sizes one and two induces `(1/3,2/3)`, rather than uniform two-state counting | State support, positive normalization and matching fiber masses explicitly. `code/audit/test_quotient_ensemble_math.py` |
| Global centering was used where all conserved records must be removed | A nonconstant protected record is globally centered and has zero Dirichlet energy | Define the Poincaré constant on the complement of the full common fixed space. Constants-only centering requires ergodicity. The same tests also distinguish compressing `log(sigma)` from taking the logarithm after compression in boundary quantum MaxEnt |
| Approximate Lyapunov descent was promoted to convergence of all limit points and readouts | A smooth periodic orbit satisfies `Phi'=-||grad Phi||^2+2` while `||grad Phi||^2=2-cos(t)/2` repeatedly exceeds two and its readout oscillates | Prove only the time-average and lower-limit bounds. Add the exact annulus counterexample, strict tie-preserving Hopfield descent plus fairness, and the local smooth-branch premises and sign for equilibrium propagation. `code/audit/test_recurrent_proof_boundaries.py` |
| A global odd-square identity was promoted to local superconformal field theory | `Q=sigma_x`, parity `sigma_z`, `H=I` satisfy `Q^2=H` without local field content or an operator-product expansion | Require an actual local continuum supercurrent/OPE receipt; central charge and current rank do not identify the full graded heterotic theory or internal lattice. `tools/test_string_continuation_boundaries.py` |
| A signed supersymmetric torus amplitude was used to read a vacuum exponent | Its factor `theta3^4-theta4^4-theta2^4` vanishes | Use ordinary unprojected vacuum characters; keep the modular spin-sector amplitude separate. Exact finite-series controls accompany the cited all-order Jacobi identity |
| Finite gauge spectra and Doob readouts could silently lose positive values | Ill-conditioned transfer operators and unresolved ground populations could produce false numerical gap/transition claims | Integrate independently reviewed PR [#1062](https://github.com/FloatingPragma/observer-patch-holography/pull/1062), preserving Mario Poneder's authorship. Full-space/high-precision controls and explicit resolution refusals support the numerical repair |

Smaller theorem corrections specify `0<epsilon<=1/2` where the two-site
heat-bath gap is written as `2 epsilon`, handle zero-noise concentration
without dividing by zero, and restore the Gaussian mismatch factor `1/2`
and fixed-precision/entropy qualifications. Live cosmology summaries now
distinguish finite-window clock estimates from proper-time limits and give
the actual adopted versus unadopted freeze status. The matching FLRW claim,
novelty summary and observation ledger use the physical weighted counts, the
two-diamond scale enclosure and the constant-profile or shrinking-window
conditions; raw flat cardinality alone does not reconstruct an expanding
proper-time ratio.

## Paper coverage

The following paths are relative to the repository root. “Targeted” means
the listed inference chains were reviewed; it does not assert that every
paragraph or imported continuum theorem was independently reproved.

| Paper root or group | Review depth |
| --- | --- |
| `flagship/from_observer_consensus_to_standard_physics.tex` | Targeted synthesis: quantum probabilities/instruments, collar entropy, cap-normal versus frame geometry, conditional Einstein chain, closure status, sixth/eighth-order dispersion normalizations and ratios |
| `paper/observers_are_all_you_need.tex` | Targeted synthesis and shared proof consumers; corrected shared ensemble and noisy-consensus qualifications |
| `paper/reality_as_consensus_protocol.tex` | Deep targeted theorem review: confluence/fairness, stochastic blocks, record persistence, conditioning, concentration and appendix consumers |
| `paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex` | Targeted shared-fragment and Lean interfaces: timelike/null tensor reconstruction, `4pi/15` and `8pi^2/15` coefficients, common physical-domain assumptions, covariance and ensemble claims |
| `paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex` | Targeted compact Lie-type classification, center and representation assumptions, distinction from physical matter/current/global group realization |
| `paper/deriving_the_particle_zoo_from_observer_consistency.tex` | Deep targeted mass/phase readout and refinement-natural optimization proofs; EW and physical-label boundaries |
| `paper/screen_microphysics_and_observer_synchronization.tex` | Targeted finite quotient laws, common fixed-space gap, native occupation/filter covariance and normalization; shared ensemble corrections |
| `paper/paradise_as_fixed_point_consensus.tex` | Scope review of fixed-point/normal-form assertions; its theological and identity interpretations remain speculative and are not mathematical or empirical predictions |
| `extra/fine_structure_constant_derivation.tex` | Deep targeted fixed-point existence, spectral/Jacobi realization, residual-scale argument and endpoint verification boundary |
| `extra/koide_identity_from_positive_c3_face_circulants.tex` | Detailed algebra and positive-chamber/balance review against `KoideCirculant.lean`; no new defect established |
| `extra/machine_checked_finite_event_algebras.tex` | Targeted partition pinching/averaging, nonzero conditioning, instrument and Born assumption boundaries |
| `extra/observable_normal_forms.tex` | Targeted exact fibers, finite Markov expectations, refinement telescoping and pseudoinverse assumptions |
| `extra/observer_patch_holography_as_string_vacuum_selector.tex` | Deep targeted criticality, local OPE, torus character and heterotic identification chain |
| `extra/thinking_as_patch_net_fixed_point_search.tex` | Deep targeted finite/continuous descent, Hopfield fairness, Gaussian objective, equilibrium propagation and fixed-point sensitivity |
| `extra/yang_mills_gap_clay_problem.tex` | Targeted finite-gap/continuum-premise separation, counterexample domains and integrated numerical transfer evidence |
| `extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex` | Detailed finite KL maximum, gradient/Hessian, capacity-transfer and graph-spectrum algebra; continuum shock and physical dictionary remain separate inputs |
| `cosmology/oph_cosmology_finite_source_cmb_program.tex` | Deep targeted radial covariance, tilt, window/tomography/inverse domains, primordial tensor inference and comparison classification |
| `cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.tex` | Targeted power-law normalization, tensor initial data, observational discrimination and clock/expansion scope |
| `cosmology/oph_cosmology_data_likelihood_contracts.tex` | Targeted tensor premises, retrospective-versus-frozen classification and model-specific exclusion |
| `cosmology/oph_cosmological_vacuum_and_structure_formation.tex` | Targeted finite vacuum versus cosmological release, ensembles, supplied action and clock/continuum qualifications |
| `cosmology/oph_boltzmann_transport_derivation.tex` | Targeted Thomson, hierarchy, acoustic, diffusion `16/15` and line-of-sight normalizations; no new defect established; no independent Boltzmann solver constructed |
| `cosmology/oph_dark_matter_paper.tex` | Targeted positive-kernel, conserved-charge, mass/BTFR exponents and missing stress/lensing attachment; no new defect established; no new SPARC fit |
| `cosmology/oph_black_hole_information_ledger.tex` | Targeted finite Stone, support-aware entropy, conditional thermal/Page and emission signs; no new defect established; no new event likelihood |

## Formal proof scope and trust

The base inventory contains 636 Lean modules and 12,324 public theorem/lemma
declarations. The audit changes no Lean theorem implementation or premise to
make a paper claim pass. The paper errors identified here are mostly outside
the formalized statements or at the point where a consumer strengthens those
statements. Existing mathematical and physical premises stay explicit.

Thirteen finite table proofs use `native_decide`; their trust includes the
native compiler/runtime and their inventory remains unchanged. A successful
Lean build checks the formal statements under their actual hypotheses. It
does not prove a physical realization, an unformalized analytic limit, or a
paper's interpretation of a theorem. The default Lake target and the optional
`OphGap` carrier target are reported separately in the validation record.

## Validation and reproduction

The mandatory runner includes the new exact counterexamples and numerical
regressions. The tests use rational probability kernels, finite matrix
identities, independent KKT solutions, high-precision original-input controls,
and mutations that recompute hashes before testing scientific validity.

```sh
python3 tools/run_mandatory_suite.py
python3 tools/check_claim_registry.py
python3 tools/check_axiom_consistency.py --check-inventory
python3 tools/check_lean_native_decide_inventory.py
python3 tools/check_lean_theorem_count.py
python3 tools/check_lean_docstring_style.py
python3 tools/build_fz_registry.py --check
python3 tools/refresh_paper_release.py --preview
cd Lean
lake build
```

Recorded local results (the groups overlap and are not summed):

- 206 focused proof-domain and numerical regression tests passed.
- 440 cosmology, radial, capacity and lensing tests passed with warnings as errors.
- 137 independently reviewed finite-gauge precision/verifier controls passed.
- 81 final particle/endpoint regressions, 110 neutrino tests, 69 Koide/A5
  controls, and 62 hadronic-consumer tests passed. The P suite's scratch-copy
  fixture failures were resolved by 36 tests on the complete tree.
- Downstream review added eight mass-gap consumer controls; the expanded
  Majorana/gap module passed all 48 tests with warnings as errors. Three
  ratio/cancellation controls fail against the original gap producer.
  Current splitting and nested export receipts were regenerated together.
- 34 consensus/string and existing paper-interface checks passed.
- 60 mandatory-runner sharding/workflow tests passed.
- The final scientific collection imports 18,090 tests successfully.
- All 23 registered TeX roots compiled; all warning budgets and the preview
  release manifest passed. The two existing neural-paper underfull-warning
  anchors were moved with their paragraphs; their counts and badness limits
  were not increased.
- The claim registry, public quantitative surfaces, axiom inventory and
  theorem-count/native-trust gates passed. The default Lake build and the
  164-step standard mandatory suite are running at draft creation.

The full-suite completion and any execution limits are recorded on the PR. No cloud simulation, new observation, physical calibration,
refitting of experimental data, or change to a frozen outcome rule is part
of this repair. Numerical tolerance checks are not outward-rounded error
certificates; Takagi congruences failing the declared normwise residual
checks are refused, and unreplayed endpoint proofs cannot authorize promotion.

The equilibrium-propagation premise/sign correction agrees with the original
[Scellier–Bengio derivation](https://arxiv.org/html/1602.05179v4), particularly
its locally continued strict minimum. The signed string-amplitude obstruction
is explicit in [Alvarez-Gaumé and Vázquez-Mozo](https://arxiv.org/abs/hep-th/9212006v1),
equations (598) and (601). The new finite counterexamples are derived directly
in the linked repository tests.
