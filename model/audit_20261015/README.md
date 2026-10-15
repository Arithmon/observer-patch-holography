# Audit of the executable OPH approach — 15 October 2026

**The checked finite computation is stronger than the physical identification
it currently supports.** The audit found one replay-check defect, exact
counterexamples to broader recovery claims, a substantially stronger empirical
comparator, and useful existing theorems. It found no counterexample to the
checked ideal finite environment, splitting or conditioning bounds. This is
an audit of the current executable approach, not a certification of every OPH
paper or proof.

Baseline: public model commit `93e59bfb8de1501174c08c800912635159ae8615`.
Three independent reviews examined mathematics, empirical provenance and
architecture, with additional inspection of source execution and replay.
New probes are retrospective diagnostics. The original protocols, source
words, data and result receipts remain unchanged.

## Findings that affect what we can conclude

| Finding | Evidence | Consequence |
|---|---|---|
| The empirical advantage was compared with weak controls. | A simpler supplied cosine Hamiltonian gives RMSE **0.317681**, versus **0.287643** for the detailed pulse model. | The gain is **9.46%** against this comparator; the roughly 69% gain over frozen/harmonic controls remains arithmetically correct. This does not distinguish OPH from conventional quantum mechanics. |
| The deposited measurement conversion is not authenticated. | All **130** position estimates satisfy `63.6*x ≈ integer` within printed rounding. The documented conversion uses `sqrt(2)*.4`, not `.636`. | Resolve readout provenance before precision or significance claims. No post hoc rescaling is justified by this pattern alone. |
| Calibration does not uniquely identify its extension to higher occupations. | Adding `N(N-1)(N-2)(N-3)/24` locally is exactly invisible to the complete compressed four-level calibration, but changes the existing eight-level continuation. | The cross-cutoff prediction tests the supplied action family. Repeating the same calibration more precisely cannot exclude this alternative. |
| The graph does not determine physical spatial dimension. | The same unit-edge cube graph has a planar embedding and produces identical existing graph dynamics and records. | Rank three requires the declared orthogonal-cube metric convention. Additional metric-sensitive observations are needed. |
| Larger sampling grids do not add field modes. | Every retained scalar grid contains the same seven modes; only one gradient coefficient changes. | Use the already available exact-gradient endpoint for fixed-band work; enlarge the mode family to investigate omitted field degrees of freedom. |
| A portable replay comparison rejected valid numerical variation. | Linux and Windows CI failed at the derived occupation-refinement position difference. A test reproduces failure of the old criterion while each underlying prediction remains within tolerance. | Fixed by validating summaries against their own predictions, then comparing the predictions with the unchanged tolerance. No physical prediction was altered. |

Details and reproducible controls:

- [Empirical review](EMPIRICAL.md): effective-Hamiltonian comparison,
  preparation mismatch, readout normalization, timing and published symmetry
  errors. [Probe](empirical_probe.py), [complete result](empirical_probe.json).
- [Architecture review](ARCHITECTURE.md): exact calibration and geometry
  ambiguities, fixed-band simplification, identification and reconstruction
  theorem mappings. [Controls](architecture_counterexamples.py),
  [complete result](architecture_results.json).
- [Mathematical review](MATH.md): checked bounds, exact source reductions,
  occupation-truncation theorem and finite-cutoff pitfalls.
  [Independent tests](test_math_probes.py).

## What survives, and what still requires an attachment

The runner compiles a supplied apparatus law into native source words,
contracts those words, evolves the spin and oscillator with a supplied finite
noise instrument, retains selection and detector branches, and compares a
derived readout with published observations. The ideal finite-model
certificates control specified approximations to that computation. The audit
found no confirmed runtime sign, dephasing-rate, selection or scoring defect
in that declared model.

That achievement does not select the supplied Hamiltonian, preparation,
noise, detector or calibrated clock from the observer architecture. The same
admitted universal controls can implement other supplied finite laws. The
source-selection tests already expose this distinction; it should govern
which next experiments and proofs receive effort. The time-independent
complete-history description and its intrinsic computation remain two views
of the proposed structure. Neither description removes these supplied
benchmark premises.

The empirical limitation is concrete. The transfer starts with zero-time
estimator residuals as large as **0.340550**, and selection/readout timing and
phase conventions remain incomplete. The largest later residual is about
**1.0823**, whereas the primary ideal finite-model selected-estimator radius
is at most **0.001517**. That radius excludes occupation removal, floating
arithmetic and apparatus/preparation uncertainty. Refining the splitter alone
will not resolve the empirical disagreement.

The host evaluates the finite instrument law using cached contracted
matrices. Resource ledgers count the formal source applications of a
realization; they are not a claim that laboratory hardware or the host
literally applies billions of pulses individually. The source contraction
uses explicitly reported floating stabilization. Its numerical diagnostics
are not certified arithmetic error bounds. These distinctions are already
present in the implementation and remain necessary in public descriptions.

## Existing results that can save work

**Occupation error:** the a posteriori Lindblad truncation theorem of
[Etienney, Robin and Rouchon](https://arxiv.org/html/2501.09607v6) applies once
the unrestricted propagator/domain and trajectory integral are controlled.
For each sideband, our independently checked half-trace residual reduces to

```
Omega*sqrt(d)/2 * ||(I_spin ⊗ <d-1|) rho_d(t)||_1.
```

Number dephasing and carrier gaps contribute no cutoff-boundary residual.
This offers a direct route beyond comparing two cutoffs. Preparation,
readout, trajectory integration and roundoff still need their own bounds.

**Native implementation:** universal-gate and Stinespring results can carry
the general finite-instrument existence argument once the actual native
services are established. Quantum-comb results provide the continuing
observer's operational object: a sequence of memory channels tested through
the allowed intervention interface. The [architecture review](ARCHITECTURE.md)
maps these hypotheses explicitly. Repeated compilation of supplied laws is
useful implementation verification, but does not repeatedly solve selection.

**Spacetime reconstruction:** corrected Osterwalder–Schrader reconstruction
is the appropriate attachment after obtaining a common Euclidean correlation
hierarchy with the requisite positivity, covariance and regularity. A finite
Hamiltonian or GNS representation does not supply those hypotheses. The
known four-dimensional scalar triviality result also constrains the ordinary
scalar scaling route; it does not rule out all finite effective models or
all OPH constructions. See the scoped references in the architecture review.

## Simplifications verified without rewriting evidence

At occupation dimension 64, exact native number drift needs **7 pulses
instead of 838**, and controlled occupation parity needs **2 instead of
2,054**, preserving the full scalar phase. These reduce two particular
operations; they are not a comparable reduction of the total dense compiler
cost. Adopt them with new word/resource receipts when changing the runtime.

The current 16-node phase mixture also greatly exceeds the ideal precision
needed relative to splitting. An outcome-independent calculation permits
4 nodes for sideband half channels and 6 for carrier gaps with environment
radius below `6.973e-6` in the specified d64 refinement. This is a proposed
budget choice, not a replacement run. An exact finite phase decomposition is
available but requires more branches here.

The unrestricted forced oscillator has an analytic displacement/rotation
factorization. **Do not substitute truncated ladder matrices and call it the
same finite propagator.** A retained counterexample agrees on the vacuum to
under `1e-8` while differing in full operator norm by more than `1e-2`.

## Recommended next work

1. Authenticate the empirical readout conversion and complete pulse/phase
   schedule, or choose a dataset with raw counts and sufficient metadata.
   Include the effective-Hamiltonian comparator in the next locked protocol.
2. Implement the boundary-residual occupation calculation and then validate
   its integral, preparation and detector comparisons. Report diagnostics
   separately until the complete bound is established.
3. Expand the identifiable action family and choose interventions that see
   its null directions, starting with the exact invisible polynomial above.
   For spatial recovery, add observables that constrain a metric beyond the
   existing graph.
4. Use the exact seven-mode endpoint; spend spatial refinement effort on
   additional modes and matched observation maps. Keep occupation, volume,
   spacing and momentum bandwidth distinct.
5. Establish a source-selected restriction that rules out an otherwise
   admissible competing law, then target an observation that distinguishes
   them. This is the route from successful supplied-law execution to an
   OPH-specific empirical test.

## Validation and preservation

```sh
python -m model.verify --hashes-only
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/audit_20261015 model/anharmonic_transfer
```

The audit manifest binds its reports, probes and results. Replay tests
recompute both new diagnostic receipts; mathematical probes exercise complete
matrices and independent generator constructions. The existing transfer
checker additionally tests amplified summary differences and rejects altered
summaries inconsistent with their underlying predictions. Its exact protocol,
word hashes and rational certificates remain enforced. Existing retained
prediction/result bytes have not been regenerated.
