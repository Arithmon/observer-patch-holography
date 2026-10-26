# Source-clocked charged response: an unreleased candidate experiment

This experiment computes a charged excitation frequency, an interaction
response and propagation to a neighboring site in **one continuing quantum
history**. A native repair record controls a later local intervention; the
resulting occupation record controls the next native repair event. The work
is prepared once. Every detector outcome and scratch abort remains in the
law, and two clock purifiers remain untouched reference columns.

The result is a finite interacting charge gap in a declared architecture.
It is not yet an identification of a particle in our Universe. In particular,
the executed four-site system, the separate free three-axis extension and
the alternative Hubbard charge-gap theorem have different domains. Joining
their conclusions would require an additional construction.

## What fixes the numbers

[model.py](model.py) calls the public native
[reference producer](../../native_modular_reference/local_preparation.py),
rather than importing measured particle data. On the 66 binary load
configurations, core occupation `q=0,1,2` has multiplicities `45,20,1` and
per-configuration weights

```text
(w0,w1,w2) = (2/225, 1/36, 2/45).
a = log(w1/w0) = log(25/8)
b = log(w2/w1) = log(8/5)
g = a-b       = log(125/64)
Omega = a+b   = log(5)
```

The architecture assigns `a,b` to alternating hopping, `g` to a centered
quartic density interaction and `Omega` to the clock generator. This
assignment is an explicit construction choice. Once it is fixed, no mass
fit or empirical coupling enters the calculation. “Computed from source”
does not mean that native endpoint repair uniquely selects this architecture.

These weights are the existing **exact ideal preparation target**. This
experiment does not execute a finite native warmup or propagate its error.
Known-state preparation, coherent pulse control, fresh pointers and the
modular-clock calibration are declared services. The committed native load
register and the encoded fermion occupation register are different systems.

## One history and its observer readout

The work consists of two identical four-mode fermion rings. Each ring has
alternating hopping and interaction

```text
H = -sum_edges t_ij (c_i^dagger c_j + c_j^dagger c_i)
    + g sum_edges (n_i-1/2)(n_j-1/2).
```

The implementation reuses the public occupation-sign CAR representation
in [gauge_matter.py](../gauge_matter.py). Its exchange signs are fermionic
by construction. The conserved occupation is a global model charge;
electromagnetic charge, physical spin and a Standard Model species are not
identified by this fact.

Two arms of a neutral interferometer have occupations `(2,2)` and `(3,1)`.
Both carry total charge four, and both cells evolve under the same law.
The measured pair frequency is `E3+E1-2E2`. A common chemical potential
cancels from it. Identical cells and particle-hole symmetry make this
frequency twice the single charged addition/removal gap `Delta4`.

Each realized history follows this order:

1. Sample the declared uniform binary native prior and execute the native
   `E` event with its fair placement coin. Only its committed current reaches
   the controller. Nonzero current enables a later kick; its exact
   probability is `5/33`.
2. Prepare the neutral charge-transfer state and two faithful-reference
   clocks once. Each clock has populations `(1/6,5/6)` and modular gap
   `Omega`. Their purification registers remain inaccessible to the reader.
3. Advance the work and clocks for modular parameter `tau/2`, commit a weak
   charge read with contrast `1/9`, then advance for another `tau/2`.
4. Commit the charge-transfer return read. Apply the native-controlled
   occupation phase at site 0, then continue the same work and clocks for
   `tau`. There is no work reset or new vacuum preparation.
5. Read occupation at neighboring site 1, then the two clock quadratures.
   The occupied outcome chooses seam `F` for the next actual native event;
   the other chooses `E`. The native integer loads continue from step 1.

[experiment.py](experiment.py) executes the quantum part;
[native_attachment.py](native_attachment.py) enumerates its complete native
attachment. There are 32 normal quantum records and 12 retained source
scratch-abort records. Including the native configurations and coins gives
8,448 normal joint paths and 1,584 terminal abort paths. Aborts receive no
invented final native event. Exact classical path weights and branch labels
retain the full quantum output without normalization or copying it into a
different preparation. This attachment starts from a diagonal native prior;
it makes no coherent `C/E` discrimination claim.

The observer decoder uses return and clock record probabilities. It does
not consult Hamiltonian eigenvalues or the simulator interval. If `v` is
the weak-read visibility and `r=-2/3` the calibrated clock contrast, then

```text
cos(phi) = (2 P(return)-1)/v
theta = atan2((2 P(clock_Y+)-1)/r, (2 P(clock_X+)-1)/r)
Delta4/Omega = phi/theta.
```

Here `phi=2 Delta4 tau` and `theta=2 Omega tau`, since the clocks are read
after the second interval. The declared phase window removes aliases. A
noncircular clock response, appreciable abort probability or an inadmissible
phase makes the point inverse inconclusive. Finite-shot confidence and
calibration uncertainty are not certified by this infinite-ensemble inverse.
The clock parameter is programmed modular flow, not laboratory seconds.

## Results and controls

The retained replay is [results.json](results.json). For the quantum program
at each intervention setting, before averaging over the native current,
representative values are:

| Quantity | Result |
| --- | ---: |
| Free finite charged gap / clock frequency | 0.4159406516 |
| Interacting finite charged gap / clock frequency | 0.6391024998 |
| Independent analytic interacting ratio | 0.6391025996 |
| Receiver occupation after the site-0 kick | 0.8281102 |
| Receiver occupation with no kick | 0.6250000 |
| Receiver occupation with work pinned after the mass read | 0.6250000 |
| Receiver change when only the later quartic interaction is disabled | 0.0005221 |

The kick, no-kick, pinned and late-interaction-off runs share the same
earlier mass-read prefix. Pinning the work leaves both clocks advancing.
The remote response therefore includes propagation and a measurable
interaction contribution after the committed gap measurement. The native
feedback also changes the final native distribution compared with fixed
seam `E`, by total variation about `0.00769`. A transcript using only native
coin/current records and the receiver bit distinguishes those feedback
policies by total variation about `0.05504`. This exceeds its ideal
splitting comparison budget. These are effects of the declared controller;
they do not derive the work Hamiltonian from unread native repair.

The gap is computed from the continuing coherent law. It is not an
eigenvalue of a native repair relaxation operator. Its free part is a
structural dimer gap, and the quartic term dresses it. This primary
experiment does not demonstrate spontaneous mass generation.

## What the finite theorems add

[DERIVATION.md](DERIVATION.md) proves the complete finite ring spectrum
independently of the numerical occupation matrices. It also treats a
**separate free eight-component three-axis extension** with

```text
E(k)^2 = 3(a-b)^2 + 4ab sum_j sin(k_j/2)^2.
```

Its gap is volume independent. A specified refinement has an explicit
`O(1/n^2)` bound toward `E^2=3g^2+Omega^2 |p|^2/4`. All eight components
remain present. This result supplies no interacting three-dimensional pole,
species reduction or physical Lorentzian identification for the executed ring.

[HUBBARD_GAP.md](HUBBARD_GAP.md) gives an alternative uniform spinful
one-dimensional law with `t=Omega/2,U=g` and no bare dimer mass. The exact
[Lieb–Wu result](https://arxiv.org/html/cond-mat/0207529) and its
[positive-integral threshold formula](https://arxiv.org/html/cond-mat/9911281)
yield a genuinely interaction-generated **thermodynamic charge threshold**.
Our quadrature-free rational certificate encloses

```text
0.0006382623246826815 < Delta_mu/Omega < 0.0006382625321756352.
```

The individual holon threshold is half this pair threshold. The spin sector
is gapless, and an electron pole is not obtained by renaming the charge gap.
The cited low-energy relativistic description is 1+1 dimensional. Small
periodic chains have large shell effects; the receipt retains them as
controls and does not extrapolate them to the certified thermodynamic gap.

## Execution and error accounting

From the repository root, using the pinned environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/mass_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -u -m model.common_source.mass_response.run --output /tmp/mass-response-replay.json
.venv/bin/python -m model.common_source.mass_response.run --verify model/common_source/mass_response/results.json
```

The runner refuses to overwrite a receipt and rejects source changes during
execution. Verification checks source custody and retained structural
checks; it does **not** rerun the experiment. Use a fresh `--output` for a
full replay. Source hashes bind the implementation, public compiler inputs
and mathematical derivations. The r2053 commit is recorded as a baseline;
this experiment is unreleased and changes no frozen r2053 artifact.

The source backend contracts actual compiled pulse words on all local
input columns and executes all instrument branches. An independent backend
uses continuous occupation matrices. Comparisons include work, all records,
untouched references and aborts in the complete classical/quantum output.
The maximum is 13 active source qubits plus two untouched reference qubits.
Pulse costs count a realized history; caching an operator contraction does
not reduce those costs. Total occupation four is preserved at completed
logical operations after preparation, not at every intermediate processor
pulse.

An exact rational ideal splitting bound covers both cells. For each
interval `T` split into `n` palindromic steps it is
`min(2, 2*4^3*T^3/(3*n^2))`; instrument contractivity and telescoping combine
the interval bounds. All three measured intervention differences exceed
the sum of their ideal splitting budgets. Floating contraction residuals
are diagnostics, not outward numerical certificates. Physical coherent
services, finite native preparation and sampling errors require separate
budgets. Only the Hubbard threshold's rational interval is an outward
numerical enclosure here.

The remaining physical bridge is concrete: an interacting spatial
continuation must retain the charged response and its measured dispersion
under volume/refinement control, with its observer clock and services
accounted. Its measured charge, statistics and response must then justify
particle identification. The finite gap, free three-axis theorem and Hubbard
threshold cannot substitute for that shared-history construction.
