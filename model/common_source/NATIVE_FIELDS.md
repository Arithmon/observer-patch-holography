# Field relations of the existing native repair

The twelve-port rule already gives an exact local continuity law and a
conditional current law. The existing uniform-seam geometry calibration
then gives graph diffusion of the mean load. These statements follow from
the executed integer repair; no Hamiltonian, new clock convention or
continuum field equation is added here.

The calculation reuses the rule in [native.py](../observer/native.py), the
source chart in [slow_geometry.py](slow_geometry.py), and the existing
[scalar-seam identities](../../Lean/ObserverPatchHolography/ScalarSeamRepair.lean)
and [repair-spectrum identities](../../Lean/Screen/PortGramRepairCovariance.lean).
The contribution is their explicit attachment to actual committed currents
and the full finite transition kernel. This is not a new general diffusion
theorem. The [finite graph transport interface](../../Lean/Thermodynamics/GraphDiffusion.lean)
already separates conservation from a physical interpretation of the
transported quantity.

## Domain and the local compatibility condition

There are twelve nonnegative integer loads with conserved total two, hence
78 states (66 have binary loads). The source chart has thirty seams. For a
selected seam `(a,b)`, compatibility means that the new endpoint loads
differ by at most one and have the same sum as before; all other loads and
earlier committed records are preserved. This condition is specified
independently of the update.

Let `d = x_a - x_b`, let `r = |d| mod 2`, and let `c` be the existing placement
coin. The actual recorded current and endpoint update are

```text
j = (d - (2c-1)r)/2,
x'_a = x_a - j,    x'_b = x_b + j.
```

Thus for the incidence column `b_e = unit_a - unit_b`, every individual
event satisfies

```text
x' - x = -b_e j.
```

This is a discrete continuity equation on the existing graph, including
zero-current events. It holds for the actual feedback-selected seams as
well as for the geometry calibration. Summing it over any region gives
the change in that region's load as its net boundary current. It identifies
conserved native load; calling that load electric charge, energy or particle
number would require a further operational identification.
The currents are committed observer records. Full load configurations and
their means are working-state audit data; this calculation does not give
the controller access to an otherwise hidden density field.

Conditional on the complete past and the selected seam, the existing fair
coin gives `E[j] = d/2` and `Var(j) = r/4`. The expected transfer therefore
opposes the endpoint load difference. Individual odd-sum events are rounded
integer repairs, not exact equal-endpoint scalar averaging.

## The existing uniform calibration closes the mean equation

In the uniform all-seam calibration, each seam has probability `1/30`
independently of the working loads and earlier outcomes. With incidence
matrix `B` and `L = B B^T`, the mean `mu_k` obeys

```text
mu_(k+1) - mu_k = -L mu_k / 60.
```

Here `k` counts native attempts in that calibration protocol. It is neither
the auxiliary history-repair solver counter nor an identified modular clock.
The actual two-event observer uses its committed record to select `E` or
`F` and then halts. It retains the local continuity and conditional current
laws, but does not have this autonomous uniform mean evolution.

An exact check builds the 78-state transition kernel `K` directly from all
4,680 state/seam/coin branches. If `X` has the load states as its rows and
`T = I-L/60`, the identity `K X = X T` proves closure for every initial
distribution. Induction gives `K^k X = X T^k`. This is stronger than fitting
a mean decay curve from sampled histories.

For `z = mu - (1/6)1`, the mean deviation `V = ||z||^2/2` satisfies

```text
V_next - V = -z^T L z / 60 + z^T L^2 z / 7200.
```

The graph's nonzero eigenvalues are `5-sqrt(5)`, `6`, and `5+sqrt(5)`.
Consequently this deviation decreases strictly unless the mean is constant.
It is not a measured thermodynamic energy. In particular, the uniform law
on the 66 binary states is stationary, individual repairs continue, and
every such trajectory has `||x-(1/6)1||^2 = 5/3`. Relaxation of the mean
does not turn a continuing history into one instantaneous fixed load state.

## What this says about wave laws

The three-dimensional slow band has mean multiplier

```text
q = (55+sqrt(5))/60 = 0.953934466291663...
```

All nonconstant multipliers are real, positive and less than one. A slow
component `f_k = q^k f_0` has

```text
f_(k+1) - 2f_k + f_(k-1) = (q + 1/q - 2) f_k,
q + 1/q - 2 > 0.
```

An undamped oscillator sampled at uniform intervals has instead coefficient
`-4 sin(omega*dt/2)^2 <= 0`. A positive scalar conversion of attempt count
to clock time cannot turn the decaying native mean into such a wave mode.
The restriction is to this closed mean-load sector with fixed linear
readouts. It excludes no general OPH construction, coherent record
observable, feedback history or quantum-memory sector.
In particular, repair-attempt count has not been identified with subjective
modular time. These decay statements do not rule out oscillatory modular
correlations of an observer's restricted state on a completed history.

Twelve ports make the finite verification exhaustive. The rank-three slow
band supplies the existing geometric quotient, not a spatial refinement
family. To obtain a wave law from more of the existing process, one must
identify its actual phase or memory observables, their closed evolution,
and the observer's clock reading on the same histories. Defining a wave
Hamiltonian from `L` would install that law rather than derive it from these
transitions. The supplied Hamiltonian in the separate interacting source
is therefore not relabeled as a result of this calculation.

## Existing quantum memory does not yet supply wave motion

The [N4 construction](../CANDIDATES.md#n4-a-reciprocal-coherent-exchange)
already retains private quantum memory. Its positive memory statistic
excludes a single closed heat mode, but does not demonstrate oscillation.
On its binary load carrier, the existing forgotten-record channel `C_e`
is self-adjoint for the Hilbert--Schmidt inner product: its current Kraus
family is closed under adjoint, with diagonal zero-current maps and paired
opposite-current maps. The retained single-seam identities

```text
E_e^2 = E_e,    E_e C_e = C_e E_e = E_e,
C_e^n = E_e + 2^(1-n)(C_e-E_e),    n >= 1
```

then give eigenvalues in `{0,1/2,1}`. In particular `C_e` is a positive
contraction as an operator on matrix space, not necessarily a projection.
The N4 seam channel is
`Ad_W (C_e tensor id) Ad_W^*`, a unitary conjugate in that same inner product.
Its uniform average remains self-adjoint with spectrum in `[0,1]` for the
already declared exchange parameter. No new choice of that parameter is
needed for this restriction.

For the existing stationary tracial reference and a Hermitian probe `f`,
the two-endpoint correlation `tau[f Phi^n(f)]` is a nonnegative weighted
sum of real decays. The [certified memory statistic](../../code/native_response/README.md#exact-memory-reduction-and-the-decisive-counterexample)
can therefore measure a mixture of decay rates without being a wave mode.
This argument concerns the forgotten-record uniform channel and its stated
correlations, not every coherent experiment on retained records. The
existing M6 source can execute reversible control words, but choosing their
sequence as a field dynamics is still part of its supplied apparatus.

Run the exact checks and retained receipt:

```sh
python3 -m pytest -q model/common_source/test_native_fields.py
python3 -m model.common_source.native_fields --check model/common_source/native_fields_results.json
```
