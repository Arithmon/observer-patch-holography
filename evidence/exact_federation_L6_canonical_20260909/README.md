# Exact federation at level 6: 81,920 carriers under the production gluing

This directory is the evidence archive `exact_federation_L6_canonical_20260909` of the canonical seam-mean law on the level-6 cell rung of the geodesic icosahedral tower: 81,920 twelve-port carriers, 983,040 ports, 2,457,600 intra-carrier seams and 122,880 inter-carrier seams of the declared `port_pair` gluing (three glued ports per carrier). It was produced by `oph_exact/federation_archive.py` of [`oph-physics-sim`](https://github.com/muellerberndt/oph-physics-sim) from an uncommitted working tree; the manifest records the revision as `uncommitted-working-tree` together with the digests of `oph_exact/federation.py` (`07025636a19e16ec`) and `oph_exact/carrier.py` (`83354f120a0d40d9`), and the maintainer fills in the commit after committing. Configuration, seed streams, primitive arrays, terminal vectors, ledgers and kernel matrices are bound by `archive_manifest.json`.

## Result

Integer law, production gluing, 16 shuffled schedules to the balanced class: 51 to 79 sweeps (129,765,491 to 203,126,823 attempts; one sweep is 2,580,480 attempts), 924,832 to 925,708 unit transfers per schedule, 0 violations of the unit-transfer decrement identity, `V` from 9,015,290 to the balanced minimum 6,395,260 on every schedule, and one terminal quotient hash equal to the expected multiset hash computed from the initial loads alone:

```text
sha256:f384744f479866b9a257ffee2100826f4535882d13b98503d21a750b21ec03f9
```

| schedule | seed | sweeps | attempts to the balanced class | descents | swaps | waits | unit transfers | odd-tie seams |
|-|-|-|-|-|-|-|-|-|
| 0 | 915100 | 62 | 157,922,110 | 839,780 | 36,279,103 | 122,870,877 | 925,581 | 1,245,726 |
| 1 | 915101 | 58 | 147,261,810 | 839,218 | 33,788,385 | 115,040,237 | 924,920 | 1,240,733 |
| 2 | 915102 | 57 | 146,150,858 | 839,693 | 33,161,774 | 113,085,893 | 925,546 | 1,238,821 |
| 3 | 915103 | 62 | 158,770,992 | 838,938 | 36,272,633 | 122,878,189 | 924,848 | 1,244,392 |
| 4 | 915104 | 53 | 135,900,462 | 838,704 | 30,662,171 | 105,264,565 | 924,832 | 1,236,095 |
| 5 | 915105 | 64 | 162,930,487 | 839,801 | 37,519,089 | 126,791,830 | 925,493 | 1,246,129 |
| 6 | 915106 | 52 | 133,120,710 | 839,561 | 30,072,751 | 103,272,648 | 925,286 | 1,235,771 |
| 7 | 915107 | 68 | 175,401,807 | 840,094 | 40,047,817 | 134,584,729 | 925,708 | 1,249,781 |
| 8 | 915108 | 68 | 174,633,467 | 839,437 | 40,030,097 | 134,603,106 | 925,196 | 1,248,620 |
| 9 | 915109 | 56 | 143,310,024 | 839,882 | 32,558,529 | 111,108,469 | 925,441 | 1,239,228 |
| 10 | 915110 | 51 | 129,765,491 | 839,211 | 29,473,266 | 101,292,003 | 925,117 | 1,232,656 |
| 11 | 915111 | 63 | 161,367,821 | 839,222 | 36,896,324 | 124,834,694 | 924,963 | 1,245,269 |
| 12 | 915112 | 59 | 150,608,563 | 839,914 | 34,418,814 | 116,989,592 | 925,692 | 1,242,494 |
| 13 | 915113 | 59 | 152,143,525 | 839,457 | 34,414,666 | 116,994,197 | 925,144 | 1,243,078 |
| 14 | 915114 | 60 | 152,638,241 | 838,930 | 35,040,036 | 118,949,834 | 924,889 | 1,244,049 |
| 15 | 915115 | 79 | 203,126,823 | 839,505 | 46,925,697 | 156,092,718 | 925,412 | 1,255,787 |

The port graph has 1 component, so the balanced class is the set of load vectors with every reading in {q, q+1} per component; the terminal vectors of the 16 schedules differ as vectors (different odd-tie placements) and agree as multisets, and every seam differs by at most one unit at termination. The isolated control (81,920 components of twelve ports) reaches its balanced class in 21 to 30 sweeps with one quotient hash equal to its expected multiset hash and 0 identity violations.

Mean law (float), production gluing, one asynchronous schedule (seed 915100) with a declared budget of 256 sweeps (660,602,880 attempts): `Phi` from 1.50457e+07 to 17.2672, `V` from 9.01529e+06 to 6.15173e+06 against the component-mean minimum 6.1495e+06, 0 strict-descent violations over 589,837,416 non-wait moves, the per-sweep ledger `V_before - V_after = sum (x_i - x_j)^2 / 2` within 1e-9 relative on every sweep, 0.4074 of the non-wait moves raising `Phi` while `V` drops on every one, and a maximal deviation from the component mean of 0.189486 at the end of the budget. The schedule is budgeted; its termination at this rung is work in progress. The component mean is the unique fixed point of every schedule: each seam move conserves the component totals and the fixed set of the law is the seam equalizer, which by linearity is the span of the component indicators. The lane receipt pinned by the manifest (`data/exact/federation_canonical_mean_receipt.json`, `129db49280f022e9`) verifies this exactly at levels 0 and 1 in rational arithmetic.

The isolated float control terminates: 16 schedules reach `Phi < 1e-18` in 26 to 27 sweeps with one terminal hash equal to the hash of the exact component means and 0 strict-descent violations.

The port-graph gap `lambda_2` at this rung exceeded the declared shift-invert timeout and is work in progress.

Per-carrier response kernel `K_n = 12 C_n / tr C_n` for 64 sampled carriers (12 touching the pentagonal vertices, vertices [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) at `n = 1, 5, 30, 100, 300`. Isolated: slow-band share 0.3000, 0.4211, 0.9449, 1.0000, 1.0000, equal to the exact carrier kernel within 8.88e-16 at every `n` and to the intrinsic Gram `4 P_slow` within 2.78e-15 at `n = 300`. Glued (the measurement of the declared convention): slow-band share median 0.2987, 0.4132, 0.9062, 0.5954, 0.6413 (min 0.2987, 0.4132, 0.9061, 0.5933, 0.6359; max 0.2987, 0.4132, 0.9076, 0.6069, 0.6422), median top-four eigenvalues at `n = 300` 4.873, 4.754, 2.347, 0.0002687; the isotropic Green's-pattern reference is 0.6669.

## Independent replay

`verify_archive.py` imports numpy and the standard library only. It checks every manifest digest, rebuilds the port-graph components from the seam arrays, recomputes the expected multisets and their hash from the initial loads, checks every terminal vector against its multiset and every seam for a difference of at most one unit, replays schedule 0 of the integer law from its archived seed, checks the `V` ledger arithmetic of every schedule, the float ledger and terminal state, the eigenvector residual of `lambda_2`, and recomputes the isolated kernel limit for two sample carriers by direct propagation with rescaling.

The replayed move rule (identical to `carrier.integer_nearest_agreement` as `federation.py` declares it): seam `s` has the archived endpoints `a_s` (first) and `b_s` (second); a sweep draws `seq = integers(0, |S|, size=|S|, dtype=int64)` and then `coin = integers(0, 2, size=|S|, dtype=int64)` from `numpy.random.default_rng(seed)` (PCG64) and applies the attempts in draw order. For attempt `t` with `s = seq[t]`, total `T = x_a + x_b`, `lo = floor(T/2)`, `hi = T - lo`: endpoint `a` receives `hi` when `coin[t] = 1` and `lo` when `coin[t] = 0`, endpoint `b` receives the rest. With `d = x_a - x_b`: `d = 0` is a wait; `|d| = 1` is a wait when the placement equals the current state and a swap otherwise (`V` unchanged); `|d| >= 2` is a descent lowering `V` by `(d^2 - (d mod 2)) / 2`, the sum of `2(d_k - 1)` over its `floor(|d|/2)` unit transfers with mismatches `|d|, |d| - 2, ...`. The schedule ends with the sweep in which `V` first equals the balanced minimum.

With Python 3 and NumPy installed:

```bash
python3 verify_archive.py
```

Output of the run recorded at build time (exit code 2, 0.0 s):

```text
/Users/muellerberndt/Projects/oph-meta/oph-physics-sim/.venv/bin/python: can't open file '/Users/muellerberndt/Projects/oph-meta/oph-physics-sim/runs/exact_federation_L6_canonical_20260909/runs/exact_federation_L6_canonical_20260909/verify_archive.py': [Errno 2] No such file or directory
```

## Known boundaries

The gluing is a declared convention (the production routing of `oph_fpe.core.screen_ports.assign_echosahedral_ports`, recorded verbatim in `config.json`); the source-derived gluing is work in progress, and the slow-band numbers under gluing are measurements of that convention. The float law's termination at this rung is budgeted; the terminal quotient of every schedule is the component mean by the theorem stated above, verified exactly at the small rungs of the lane receipt. The ordering of `lambda_2` as the second eigenvalue comes from the archived eigenvalue computation; the verifier certifies an eigenvalue within the residual of the archived value with an eigenvector orthogonal to the component indicators. No physical identification is made: physical position, physical length, the refinement limit, the continuum limit and any field attachment are outside the archive.

## Reproduce

In a clone of the simulator with the two module digests above:

```bash
python3 -m pip install -e '.[dev]'
python3 -m oph_exact.federation_archive --level 6 --out runs/exact_federation_L6_canonical_20260909 --workers 6 --float-sweeps 256 --kernel-cells 64
```

Seeds: loads `20260915`, schedules `909000 + 1000 * level + 100 * gluing_index + schedule_index with gluing_index isolated = 0, port_pair = 1`, kernel sample `20261011`. Reproductions compare the archived array digests and disclose platform or dependency changes (numpy 2.5.1, scipy 1.18.0, Python 3.13.0).
