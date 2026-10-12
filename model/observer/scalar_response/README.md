# A scalar detector inside the native observer protocol

This runnable construction connects the native record-and-feedback process
to the existing massive free-scalar continuum family. A first native current
sets a local field preparation, a separated binary detector reads the evolved
field, and its retained bit selects the second native repair. The field action,
vacuum, geometry and model time are explicit apparatus inputs.

The new result is the executed connection and its error transfer. The existing
scalar Hamiltonian is not derived from integer load averaging. This is a
conditional free-field construction with a testable apparatus prediction,
not an observed-universe simulation or empirical confirmation of OPH.

## Run

Install the [model requirements](../../requirements.txt), then run from the
repository root:

```sh
# Four grids, 200,000 independent complete trials per grid.
python3 -m model.observer.scalar_response --check model/observer/scalar_response/results.json

# A new seed and a fresh receipt; existing files are never overwritten.
python3 -m model.observer.scalar_response --trials 10000 --seed 42 --write /tmp/scalar-observer.json

# Local field updates, temporal refinement and native feedback.
python3 -m model.observer.scalar_response.local_evolution --check model/observer/scalar_response/local_results.json

# Independent field, quantum-instrument and native-history checks.
python3 -m pytest -q model/observer/scalar_response

# Replay the owning directed-interval continuum certificate.
python3 code/source_scalar_packet/verify_source_common_scalar.py

# Independent temporal certificate, with its spatial parent replayed.
python3 code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py
```

Replay compares source bindings, exact rational bounds, record counts and
structure exactly. Only floating diagnostics have a `2e-10` absolute replay
tolerance for platform variation. This tolerance is not a certified numerical
error bound. The parent interval verifier independently certifies its ideal
mathematical source; our numerical diagonalization does not replace it.

## The actual process

1. Prepare one of the 66 native binary total-two states uniformly. Repair
   `(0,1)` with a fair coin and retain `(coin,j)`.
2. Prepare a fresh scalar vacuum on the same golden tensor grid used by the
   parent theorem. Apply the compact left-slab momentum kick `4*j*f` and
   evolve its supplied local Hamiltonian to model time `t`.
3. Couple a centered pointer to the compact right-slab field and read a
   binary result y. Both outcomes have explicit Kraus operators. The program
   retains their conditional states as two Weyl-displaced coherent components,
   including their relative phase and Gram matrix.
4. Commit y. The actuator receives only that bit, selecting seam `(0,5)` for
   +1 and `(0,1)` for -1. Execute the second native repair, append its record,
   and halt. The postmeasurement field is terminal.

This is two native repair events and one field observation per trial. Each
trial starts with fresh source and field preparations; a batch is not a
continuing observer with unlimited memory. Later feedback preserves both
earlier classical records. No outcome is discarded.

`field.py` independently assembles the local mass and stiffness forms from
the actual golden addresses and evolves every oscillator mode. Tensor
separation reduces memory to two-dimensional slices; no nine-mode projection
or Fock occupation cutoff is substituted for the full field. At q=233,
12,487,168 mathematical modes contribute. Small-system coherent-state and
occupation-basis tests independently check the quantum instrument and its
backaction. The software does not emit a native quantum-control pulse word
or allocate a physical oscillator network of that size.

## Local evolution and temporal refinement

`local_evolution.py` executes neighboring-site canonical kick–drift–kick
updates at q=13 and q=34. It evolves the prepared mean and two detector
pullback columns by the local tensor stencil. Contracting those columns with
the original vacuum covariance includes the squeezing caused by split
evolution. Keeping the continuous vacuum variance unchanged would describe
a different detector experiment.
Terminal branch records retain the two Weyl-shifted squeezed components,
their coefficients and their Gram overlaps.

The model time is the nearest certified update to `39/40`, retained exactly
as the step count times the rational step size. At q=13 an additional run
halves the step at the same actual model time. A spectral split propagator
independently checks the local recurrence. The computed detector probability
then drives the same native record-and-feedback protocol and controls.

The [local receipt](local_results.json) retains checkpoint summaries and
sampled native records, not an authenticated log of every field read and
write. This is numerical execution of the supplied scalar dynamics; no
native quantum-control pulse word or laboratory oscillator array is executed.
q=89 and q=233 are covered by the inherited analytic temporal certificate,
but their local update loops are not executed in this receipt. The existing
four-grid spectral receipt and mathematical certificates retain their scope.

## Prediction and continuum error

Let `r_q(t)=Pr(y=+1 | j=+1)-1/2`. Exact enumeration of all 528 histories gives

```
E[j2] = 5*r_q(t)/132
E[final port-0 load] = 1/6 - 5*r_q(t)/132
```

The marginal detector bit is fair. Its correlation with the first current
causes the feedback effect; an independent donor bit has the same action
frequencies and no mean current. Reversing the field preparation flips the
signal. Zero preparation or zero detector coupling removes it. A separate
fixed-seam control disconnects the detector and also has zero mean current;
that control does not preserve the action frequencies.

At model time `39/40`, the computed detector-plus probability for j=+1 is
approximately 0.334515, 0.319421, 0.317093 and 0.316771 on q=13,34,89,233.
The finest-grid prediction for mean j2 is approximately **-0.00694050**.
The retained execution tests fresh Born outcomes against predictions computed
before sampling, and retains paired reversal, independent-donor and fixed-seam
controls. Endpoint times `19/20` and `1` are checked without fitting.

The inherited bound `|r_q-r_continuum| <= epsilon_q` transfers to

```
TV(complete recorded histories) <= 5*epsilon_q/33
|E[j2]_q - E[j2]_continuum| <= 5*epsilon_q/132
```

At q=233 the second bound is below **0.001897**, uniformly over
`19/20 <= t <= 1`. The ideal finite-grid mean is at most **-0.005217**
throughout that window. The first two grids have only the parent's trivial
error bound; q=89 and q=233 tighten it. The parent theorem proves the
continuum limit on its stated refinement family; these four certificates
give concrete finite bounds, not a convergence rate inferred from samples.

For local split evolution, replace `epsilon_q` in the transfer formulas by
the parent's combined spatial and temporal probability bound. Comparing to
a target between update times also includes its timing allowance. At q=233
the certified temporal term is at most `0.000008071512` and the combined
grid-time probability error is at most `0.050062606564`; these are analytic
bounds for the ideal split dynamics. At the locally executed q=13 and q=34,
the inherited spatial bound remains trivial. Floating agreement and step
refinement check the numerical implementation and do not certify its roundoff
error or remove the supplied field law.

The [derivation](DERIVATION.md) supplies the composed instrument and exact
transfer proof. The [prediction protocol](PREDICTION.md) fixes settings,
controls and a prospective statistical rejection rule. It is independently
testable given the declared apparatus and its service-error budget, but
there are no laboratory data here and no OPH-specific new physical constant.
P/N closure is not used to fit this unitless free-scalar benchmark. Physical
units, calibration, interactions and hardware errors remain further work.

## Files

| File | Purpose |
| --- | --- |
| [field.py](field.py) | Full-mode dynamics, Weyl branches and coherent-state oracle |
| [experiment.py](experiment.py) | Native composition, exact error transfer and sampled trials |
| [results.json](results.json) | Four-level execution and inherited rational bounds |
| [test_field.py](test_field.py) | Independent action and quantum-instrument oracles |
| [test_experiment.py](test_experiment.py) | Independent history, control, replay and bound checks |
| [local_evolution.py](local_evolution.py) | Local stencil evolution, squeezed detector probability and native continuation |
| [local_results.json](local_results.json) | Executed q=13,34 comparisons, step refinement and inherited temporal bounds |
| [test_local_evolution.py](test_local_evolution.py) | Independent local-action, split-map, Gaussian-state and refinement controls |
