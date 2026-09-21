# Native checkpoint selection

This package derives a finite record-publication schedule from native
observation constraints and an explicit information-projection problem.
It extends [temporal recovery](../source_temporal_acceptance/README.md).

The main new path result is analytic: with two independent unknown records
at ports 0 and 1 of a calibrated path 0--...--d, publishing both at d needs
exactly **2d-1 means**. The shortest histories form a two-wave precedence
order with **Catalan(d-1)** linear extensions. Every positive stationary
edge-product reference selects the uniform law on those histories. The
constraint grammar, feasibility, general finite selection/policy theorem,
native commutation and decoder/error algebra have Lean proofs; the full
path classification is not Lean formalized.

On captured W12 seams 0--1--14--23--45, there are exactly five shortest
seven-mean histories. The derived online kernel publishes both inputs on
every selected history; it does not assume a successful word or an IID
proposal stream. A complete census checks 194,377 words across twelve
uniform/tilted/aggregate/deadline controls. Thirty policy trajectories each
replay all nine signed preparations. Separately, 1,969,761 words exhaust the
shortest-horizon language through path depth five, and all 6,918 positive
histories through depth ten are independently replayed.

The reference, complete *chosen* move grammar, public record domain,
deadline and full-history entropy cover are explicit inputs. Source
selection of those inputs, M1's metric radius, physical precision/control,
immutable-version semantic order and geometric identification remain open.
Optimizers for independently chosen deadlines need not be compatible.

See [CONTRACT.md](CONTRACT.md) for Objective, Deliverables and Exit, and
[DERIVATION.md](DERIVATION.md) for the proof and exact source boundary.

```powershell
$env:PYTHONPATH='code'
python -m source_checkpoint_selection.build
python -m source_checkpoint_selection.verify --write-receipt
python -m source_checkpoint_selection.verify
python -m pytest -q code/source_checkpoint_selection
```

The producer uses response-space dynamic programming and Dyck enumeration.
The verifier uses exhaustive scalar basis runs, determinant tests and
independent prerequisite-poset enumeration. It never imports the producer.
Only compact controls and census commitments are retained; expanded words
are recomputed. Source pins include the complete local Lean import closure,
the canonical axiom specification, captured support and executable code.
