# Exact federation under the canonical law: tower levels five to ten

Settlement receipts of the production gluing of the geodesic icosahedral tower at levels
5 to 10 (20,480 to 20,971,520 twelve-port carriers; 645,120 to 660,602,880 seams), produced
on 2026-09-24 by the lean engine `oph_exact/federation_huge.py` of the simulation repository:
levels 8, 9 and 10 on an r7i.16xlarge host, levels 5, 6 and 7 on a laptop (level 5 with four
schedules, the others with sixteen). The level-six archive
`evidence/exact_federation_L6_canonical_20260909` remains the canonical public run of the law;
the level-six receipt here reproduces its quotient hash and its sweep range, and this package
extends the run in scale under exactly the same conventions.

## What was run

For each level: the loads `default_rng(20260909 + level).integers(0, 6)` on every port; sixteen
schedules of the integer nearest-agreement law with seeds `909000 + 1000 L + 100 + k`, each
schedule `|S|` i.i.d. seam draws per sweep followed by `|S|` tie coins from one PCG64 stream,
run to the exact attempt at which `V = sum x^2` reaches the balanced-class minimum; one
budgeted schedule of the real seam-mean law (64 sweeps); the per-carrier normalized response
kernels `K_n = 12 C_n / tr C_n` on 64 sample cells (one per pentagonal vertex, the rest seeded)
at `n = 1, 5, 30, 100`, and at `n = 300` on the first four cells.

The conventions (seam order, load and schedule seeds, draw law, move laws, termination,
canonicalizers) are those of the level-six archive and are restated in every `receipt.json`.
The engine differs from the archive's producer only in representation (an `int8` state, seams
addressed by index arithmetic, memory-mapped geometry, draws taken in chunks that reproduce the
single-call stream bit for bit, response kernels propagated on the ball of graph radius
`2n + 1`, which is exact because the operator moves one hop per step). At level six the
engine reproduces the archive's schedule 915100 on every field, including the terminal vector
and all 62 per-sweep draw digests, and the archive's kernels to the last bit.

## Contents

- `L{5,...,10}/receipt.json`: the level receipt (geometry digests and gluing statistics,
  expectations from the loads, all integer schedules with their `V` ledgers, per-sweep draw
  digests and per-sweep state digests, the mean schedule, the kernels, timings, platform and
  module pins).
- `L{5,...,10}/integer_{seed}.json`, `mean_{seed}.json`: the per-schedule records with the
  per-sweep move counts.
- `L{5,...,10}/verify_{k}.json`: the record of the full independent verification run on the
  producing host for schedule `k` (see below).
- `L{5,...,10}/texture.json`: the settling readouts of `oph_exact/federation_texture.py`
  (descent curve, transport cost, retention of the coarse initial imbalance by scale, angular
  power spectrum of the settled field).
- `L{6,8,9,10}/readback.json`: the per-carrier slow-band readback of the settled states
  (`oph_exact/federation_readback.py`): band energy shares over the eleven-dimensional
  zero-sum port space (bands of dimension 3, 5, 3 against the isotropic 3/11, 5/11, 3/11),
  and the retention of the slow and fast fields by the settlement.
- `reserve/reserve_L{5,...,9}.json` with `reserve/verify_federation_reserve_independent.py`:
  the reserve-generator readouts of `oph_exact/federation_reserve.py` (research-repository
  issue 985): per-sweep unit transfers, swaps, descents and waits by seam class (intra-carrier,
  intra-face, collar) and orientation, the tick survivals and generators of the face-load
  reserve, the excess contraction, and the covariance survival under one refinement step read
  from the textures. Every reserve schedule replays the engine's own draws and reproduces the
  settlement receipt's terminal digest, sweep count, ledger and move totals.
- `verify_tower.py`: a standalone verifier that imports no simulator code.
- `test_verify_tower.py`: the verifier passes on the package and fails on mutations.
- `manifest.json`: SHA-256 of every file and a summary per level.

The geometry caches (up to 1.3 GB) and the terminal arrays (252 MB per schedule at level ten)
stay on the producing host; their SHA-256 digests are in the records.

## Verification

Two layers. On the producing host, `oph_exact/verify_federation_huge_independent.py`
(shares nothing with the engine but the carrier seam template) recomputed the cell count and
the components from the cached inter-seam table, checked every port glued at most once and the
port-pairs digest, recomputed the loads and from them the balanced minimum and the expected
component-multiset hash, and for every integer schedule checked the terminal array's digest,
dtype and range, its membership in the balanced class, `V`, conservation, the recomputed
multiset hash, the monotone ledger, every per-sweep draw digest re-drawn from
`default_rng(seed)` with single-call draws, and a layered numpy replay from the last retained
checkpoint through every intermediate state digest to the terminal state; for the mean
schedule the digest, `Phi`, `V`, the deviation and the withheld-hash rule; and for the kernels
symmetry, trace, the slow-band shares and a dense recomputation on the local ball. The
`verify_{k}.json` records hold what was checked (draw sweeps checked, sweeps replayed).

In this package, `python3 verify_tower.py` recomputes what the receipts can be held to without
the host data: the tower counts, the loads and the expectations derived from them, the
balanced-class terminal values, the expected hash against every schedule's quotient hash, the
ledgers, the attempt indices, the digest chains, the per-schedule records against the receipt,
the slow-band shares from the kernel matrices with the carrier's projector rebuilt from its
thirty seams, the host verification records and their coverage, and the manifest; it then runs
the standalone reserve verifier on every reserve receipt against the settlement receipts and
textures of the package.

## Reserve-generator readings

The collar at depth `m` is identified with the `30 * 2^m` inter-face seams of the production
gluing, the certificate's oriented slot count. The depth-scaled collar hazard per sweep is one
function of the sweep index at every level (0.025 to 0.026 at sweep 1, 0.012 once settled);
the forward share of collar crossings is 0.498 to 0.504 (the loads carry no orientation); the
descent excess contracts by 0.81 to 0.83 per sweep. Under the certificate's dyadic tick of
`2^m` sweeps the tick outlives the settlement from depth 6 upward; where a tick fits, the
face-load generator reads 0.012 (level 5) or 0.006 (half-length ticks, levels 5 to 7), against
the targets `P*/24 = 0.068` and `P*/48 = 0.034`; the settled field's covariance survives one
refinement step at one quarter in the angular basis (`theta` about 2, as for the i.i.d. loads),
against the required 0.977 (`theta = 0.034`). The declared loads and repair law do not emit the
reserve generator; the obstruction is the record ensemble, not the repair law.

## Nonclaims

The gluing is a declared convention, not source-derived; no refinement or continuum limit,
no field attachment, no physical identification and no physical clock are asserted. The real
seam-mean law is budgeted and does not terminate at these levels within the budget; its
lattice-snap hash is withheld. The spectral gap is not computed at these levels. The
response-kernel readings at depth are properties of the gluing convention; cofinal gluing is
not constructed.
