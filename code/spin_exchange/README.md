# Two-port exchange discriminator

This interface compares declared CAR, CCR and distinguishable-particle models
using the same two-port Cayley hopping coefficients. It imports the `e_c`
channel coefficients from `sm_fermion_current/quantum_link_receipt.json`:
`U = (3 I - 4 i X)/5`, supplied hopping strength 2 and step 1/2.
Only those coefficients are reused. The parent's electric-flux dressing,
operator Gauss constraint, initial state and execution history do not transfer.

From the research repository root:

```sh
python3 code/spin_exchange/build.py --write
python3 code/spin_exchange/verify.py
python3 -m pytest -q code/spin_exchange/test_spin_exchange.py
```

The producer uses creation polynomials with occupation-factorial inner products
for CCR and anticommuting creation for CAR. The independent verifier uses the
four-dimensional labeled tensor space and symmetric/antisymmetric states. Both
use rational arithmetic. The receipt pins the parent, its two source modules and
scientific fragment, and this package's producer, verifier and tests.
The tests mutate statistics, normalization, source binding, error bounds and
scope, including coherent normalized wrong-statistics rewrites.

The new preparation puts **two excitations in the same internal/spin channel**,
one in each spatial input port. The channel uses declared local spin frames;
the raw spin vectors at the ports need not coincide. Its theoretical coincidence probabilities are
1 (CAR), 49/625 (CCR) and 337/625 (orthogonal hidden labels). Every model has
one-body density equal to the identity. On conserved one-excitation-per-channel
sectors, both algebras restrict to the same matrix units: one-body transport and
mean-current feedback cannot discriminate them. A supplied single-particle
central spin sign also agrees between the models. Neither check selects CAR.

For partially overlapping hidden internal states, the squared overlap `eta`
reduces the coincidence gap to `576 eta/625`. Exact controls at `eta=0,1/4,1`
are replayed by propagating a reduced labeled-tensor density matrix. Orthogonal
hidden states erase the exchange signal; identical-channel preparation is an
essential premise, not something inferred from these probabilities.

The declared preparation trace-distance budget and readout effect norm budget
are each 1/100. Their sum bounds coincidence probability error for the **exact
supplied transport**. These are assumed budgets, not calibrated uncertainties;
transport implementation error and finite-sample uncertainty are not included.
The receipt contains ideal probabilities, never sampled detector outcomes.

For an OPH operational attachment, a bounded self-reading patch must supply
local preparation state, two ports, retained coincidence readback, records and
feedback or repair, together with a public evidence bundle for those error
bounds. This calculation provides the mathematical discriminator and its
interpretation boundary. It supplies neither that physical operation nor a
source-derived statistics law, continuum spin--statistics theorem or an
impossibility theorem for the observer axioms.
