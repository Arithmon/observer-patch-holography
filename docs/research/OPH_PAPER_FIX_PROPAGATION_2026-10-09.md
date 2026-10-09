# Paper correction propagation, 2026-10-09

This follow-up checks the active paper consumers of the October 8 proof audit,
the model-era corrections `a20d4736` and `7521d5c3`, and PRs #1065–#1069,
against main `fd847bafdae30b0b6c30f437fe874fb1f036f914`. It includes the
standalone flagship, which does not import the shared proof fragments.
The earlier audit reviewed that paper but did not change its source or PDF.
The earlier correction-disposition report was bounded; it did not establish
that every duplicate paper statement had been corrected.

## Missed consumers corrected

| Consumer | Missing condition or incorrect provenance | Correction |
| --- | --- | --- |
| Flagship thermodynamics section and *Observers Are All You Need* thermodynamics section | A faithful reference and arbitrary conserved moments were said to yield a finite-multiplier Gibbs optimizer | Require a faithful feasible state on the moment surface. Boundary optima can have smaller support and need a restricted or limiting treatment. |
| Observer-paper Higgs/top table and two synthesis-fragment passages | The displayed boundary-scale candidate was called frozen or fixed before comparison | Use comparison-exposed, consistent with the existing candidate audit and the observer paper's comparison paragraph. |
| Gravity paper's shared Kerr-comb section and black-hole information paper | The normalized coordinate and offset-subtracted tooth ratios omitted the nonzero thermal-scale domain | State positive mass and strictly subextremal spin. At extremality the template teeth collapse to the rotation line and those ratios are undefined. |

The Gibbs obstruction is elementary: with reference `I/2`, conserved observable
`|1><1|`, and required expectation zero, the only feasible density matrix is
`|0><0|`. Every finite Gibbs multiplier instead has positive excited-state
probability. The existing exact tests in
`code/audit/test_quotient_ensemble_math.py` retain this boundary case, and the
shared `PAPER.tex` maximum-entropy statement already had the needed premise.
The Lean conditional identity assumes a positive exponential state; its
theorem body needs no repair.

Higgs/top provenance is recorded in `code/particles/MASS_CANDIDATE_STATUS.md`
and `code/particles/runs/calibration/d11_boundary_scale_selection_audit.json`:
the candidate was identified after the implied-scale calculation. A narrower
freeze before a future higher-order calculation does not make the displayed
two-loop comparison prospective. No mass coordinate or comparison verdict
changes.

The ringdown domain already appears in the repaired numerical contract from
#1066. Its explicit paper qualification is mirrored in the mutable pending
FZ-14 description and the live OL-B4 observation note, with both Markdown
surfaces regenerated. Historical frozen rows, custody bytes, likelihood
policies, and scientific statuses are unchanged.

## Consumers requiring no further correction

The flagship already requires one common refinement tower with the named
vanishing remainders for its Einstein implication, and supplies its matter
representation explicitly. Those qualifications remain consistent with the
collar, modular-transfer, and declared-charge corrections.

The reviewed consensus/compiler consumers require the specified functional
update and compiler correspondence. Checkpoint consumers retain the finite
D=24 scope and the missing full source-class/physical attachment. Native
moment consumers use unchanged short-lag rational matrices and an independent
enumeration. Whitney consumers distinguish analytic action and convergence
results from numerical evaluation checks. No additional false paper claim
was established from #1065, #1067, #1068, or #1069.

## Validation and limits

- All 23 registered TeX roots rebuilt, including the flagship; all warning
  budgets passed with no overfull boxes or reference, citation, glyph, or font
  errors. The release manifest was regenerated and validated.
- Four PDFs changed: flagship, observers, gravity, and black-hole information.
  Other rebuilt paper artifacts were unchanged. The manifest tracks the eight
  core publication PDFs; the black-hole paper remains a supplemental artifact.
- Visual checks cover flagship page 52, observers pages 63, 78, 80, 148 and 149,
  gravity page 98, and black-hole information page 7.
- The final focused run passed 348 tests: quotient/MaxEnt mathematics, ringdown
  arithmetic, public surfaces, release manifests, and observation-ledger gates.
  The first invocation used a nonexistent test filename and was corrected.
  An intermediate manifest check ran during rebuilding and correctly detected
  the temporary mismatch; the completed build and subsequent full focused run
  passed with synchronized artifacts.
- Reader style, claim registry, active axiom inventory, source-current inventory,
  frozen-register custody/history, observation and premise registers, and
  whitespace checks passed.
- Initial remote mandatory checks on `ffcbe37f` found three live generated
  records still bound to the old complete prediction-register hash: forecast
  contract state, discriminator stratum 643, and the dependent closure
  preflight. Their canonical producers regenerated only hashes and one byte
  count; all scientific payloads and verdicts are unchanged. Producer and
  independent replay checks pass, as do 108 focused tests and eight subtests.
  The paper-preview CI on `ffcbe37f` passed. The follow-up commit's mandatory
  CI checks the synchronized records.

This is a focused correction-propagation review, not an independent reproof
of every paper. It establishes no new physical prediction or joint model of
the Universe. The book and served website were not re-audited or rebuilt.
No local model campaign or its frozen evidence was published or modified.
Remote CI is reported separately for the actual pushed commit.
