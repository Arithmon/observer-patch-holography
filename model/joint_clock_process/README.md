# Joint local clocks and quantum fields with explicit memory

This packet supplies an exact finite correspondence between native repair
histories and a declared observer's effective field, two local clocks and
interventions. It preserves complete outcome probabilities and conditional
quantum outputs, including an untouched reference. It extends the earlier
[first-moment memory identity](../federated_observer/MEMORY.md) to quantum
instruments. It also finds a simpler, exact description at this experiment's
completed clock events: a four-mark classical predictor, together with the
retained conditional reference operators.

The construction is a finite operational reduction. It does not establish a
continuum field theory, a relativistic spacetime, or an identification of native
event order with modular flow. The field below is an observer-chosen accessible
algebra; its coordinates are not fundamental fields or a derived spatial
dimension.

## One source, two local clocks

The existing [federated quantum packet](../federated_quantum_sector/README.md)
supplies the source path with active ports `(0,1,72,82)` and seams
`(0,1), (1,82), (72,82)`. Ports 0 and 1 belong to one elementary twelve-port
carrier; ports 72 and 82 belong to another. The controlled experiment chooses
these three seams uniformly and uses the original fair native placement coin.
It retains one particle on the four active ports, with all other 236 source
ports empty. Every selected operation preserves that support, including all
superpositions and reference correlations. This schedule is a declared service;
it is not the uniform schedule on the entire federation.

Clocks A and B read the signed native occupation changes at ports 0 and 72.
Their joint mark is one of

```text
no event, A exit, A entry, B exit, B entry.
```

Only actual nonzero repair currents tick a clock. The seam and coin remain
hidden CP alternatives. The no-event map Q includes genuine remote motion
between ports 1 and 82. The declared joint record interface makes both clocks'
records available for comparison; message delivery times and physical
simultaneity are not inferred. Native attempt number is a resource counter.
The two tick counts are observer-local records. Neither parameter is asserted
to be a modular parameter or calibrated seconds.

Let R_m be a nonzero-mark CP map. Direct source evaluation on all sixteen
matrix units gives `Q*(I)=(5/6)I`. Thus the chance of no union event within k
attempts is exactly `(5/6)^k` for every initial state and untouched reference.
The eventual instrument is

```text
S_m = R_m (I-Q)^(-1) = sum_(n>=0) R_m Q^n.
```

The uniform survival bound proves convergence and completeness. The expected
cost of the next union event is six native attempts. Dividing R_m by the
one-step event probability is still wrong: Q contains remote motion and is not
`(5/6) Id`. Independently inserted identity holds cancel from the completed
instrument, by the same resolvent identity as the
[single local-clock construction](../relational_event_clock/README.md).
They change finite-cap probabilities and resource costs.

## An exact event-boundary simplification

In the ordered active basis, the four completed maps have the exact form

```text
S_m(rho_WR) = |d_m><d_m| tensor Tr_W[(E_m tensor I_R) rho_WR],

d = (1, 0, 82, 72),
E_Aexit  = diag(1, 0,   0, 0),
E_Aentry = diag(0, 2/3, 0, 1/3),
E_Bexit  = diag(0, 0,   1, 0),
E_Bentry = diag(0, 1/3, 0, 2/3).
```

The addresses in d are native port addresses, not matrix indices. Each final
nonzero-current Kraus branch has rank one on this one-particle support. The
source calculation checks the full measure-and-prepare identity, not only
traces. It therefore holds with every untouched reference.

The next-mark transition matrix, with columns denoting the preceding mark in
the order above, is

```text
       A exit  A entry  B exit  B entry
A exit     0       1       0       0
A entry  2/3       0     1/3       0
B exit     0       0       0       1
B entry  1/3       0     2/3       0
```

The stationary mark probabilities are each 1/4. Consecutive union events use
different clocks with probability 1/6, rather than the 1/2 from independent
uniform clock identities. Starting after any event, the probability of no A
tick during the next 2k union events is at most `(2/3)^k`; the same bound holds
for B. Thus each local clock continues almost surely. These counts and costs
do not supply physical clock rates. This chain is periodic; a stationary law
does not mean convergence of every initial mark distribution.

Four marks suffice for future work-record probabilities under fixed or
record-controlled work-only instruments. Such an intervention updates the next
transition using its action on the known state `|d_m><d_m|`. A controller that
consults more history must retain that classical history. Complete reference
outputs additionally require the conditional reference operator in the formula;
four scalar probabilities cannot replace the reference assemblage. Operations
on a correlated reference or additional hidden controller memory are outside
this four-mark closure.

C and E have exactly the same completed event instrument here. Each completed
event destroys work coherence in this chosen one-particle experiment. This is
a property of this instrument and sector, not a claim about all observer
clocks. Coherence survives in unresolved waits and can be read before an event;
discarding those waits would remove the discriminator below.

## Full quantum field histories require memory

Identify the ordered four-dimensional support with `field qubit tensor hidden
qubit`: the field bit names the carrier and the hidden bit names its within-
carrier address. The effective field algebra is `M2 tensor I2`. Coherent field
operations are declared joint apparatus, not inferred spatially local gates.
Define the state restriction R by tracing the hidden factor and choose the
trace-preserving lift `J(x)=x tensor I2/2`. Then `RJ=I`. For a retained record
word, let L_j be its actual outcome map, including its declared intervention.
With `F=I-JR`, define

```text
A_j = R L_j J,    B_j = R L_j F,
C_j = F L_j J,    D_j = F L_j F,
x_0 = R rho_0,   y_0 = F rho_0.
```

Block substitution gives, for every finite word,

```text
x_(n+1) = A_(n+1) x_n + B_(n+1) D_n ... D_1 y_0
          + sum_(j=0)^(n-1) B_(n+1) D_n ... D_(j+2) C_(j+1) x_j.
```

An empty product is the identity. This follows by inserting
`y_(j+1)=C_(j+1)x_j+D_(j+1)y_j` repeatedly into
`x_(n+1)=A_(n+1)x_n+B_(n+1)y_n`. No Markov, Gaussian, equilibrium or unitarity
assumption is introduced. Outcome maps are unnormalized; their traces are the
full word probabilities. The identity includes arbitrary CP interventions
whose selection depends on retained records. Applying it to each reference
matrix block proves the same full quantum output identity with any untouched
reference. The finite executable guards the recurrence at eight maps; the
algebraic identity itself holds at every finite length.

The memory terms and hidden forcing are generally not positive maps. Their
complete sum is the restriction of the original CP history. This predictor
still uses source operator blocks and the declared initial hidden component;
it does not prove an efficient finite memory window or reconstruct an unknown
preparation from field records alone.

Two checks make that boundary concrete. Preparations at native ports 0 and 1
have the same reduced field state. Their completed first-event laws are
`(1,0,0,0)` and `(0,2/3,0,1/3)`, with total variation one. Even if the initial
hidden factor is maximally mixed, removing the memory kernel changes the
word `(no event, A exit)` from exact probability 5/72 to 1/16. The full memory
recurrence restores every reduced matrix entry, not just these probabilities.

## What executes and what it predicts

The shared primitive runtime executes three actual native attempts on continuing
four-by-four matrices, then a complete field-X measurement. The chosen coherent
preparation is `(port 0 + port 72)/sqrt(2)`. A controlled variant applies a field
Hadamard after each recorded A exit. This supplied apparatus can change
occupation; the clock is the specified native-current detector and does not
silently count the apparatus as another native tick. Its visible trigger and
resource charge remain in the experiment. Identity feedback makes no change.

Every five-symbol attempt word and both detector outcomes are retained:
`5^3 times 2 =250` terminal records per candidate and feedback choice. Records
with fewer than one or two events explicitly retain their timeout status and
full native state. No state is reset at an event. A checkpoint after the first
attempt and its feedback is serialized, imported by prefix replay and
continued. Two-attempt entangled-reference runs retain 50 complete records
each. The exact memory recurrence is checked against every terminal reduced
field/reference matrix, and every full work/reference matrix is compared with
the actual Kraus runtime.

Without feedback, the field-plus probabilities after three attempts are
**35/54 for C** and **14/27 for E**. The common probability of no clock event
is 125/216. Within that no-event branch the corresponding conditional plus
probabilities are 189/250 and 133/250. The complete law keeps that branch and
all successful-event branches; these are not postselected success scores.
Both candidates have the same completed clock law, so the field measurement
adds information absent from the clocks alone.

These are conditional model predictions. The declared support, preparation,
schedule, coherent apparatus and detector need an experimental correspondence
before comparing them with laboratory data. No public dataset is identified
with this protocol, and no empirical confirmation is asserted. Floating runtime
residuals are numerical diagnostics; this packet supplies no compiled pulse
words or hardware-error certificate.

## Theorem attachments and remaining obligations

The [native signed-difference reduction](../observer/QUANTUM_MIXING.md) and
[federated quantum packet](../federated_quantum_sector/README.md) supply the
source Kraus maps, invariant support and arbitrary-reference extension. The
event resolvent and irrelevant-hold identity reuse the
[local-event construction](../relational_event_clock/README.md).
[Regional continuity](../../Lean/Screen/RegionalContinuity.lean) identifies each
native signed current with the recorded port change; it does not describe the
separately supplied coherent apparatus. The
[time-order ledger](../../Lean/Time/TimeOrderLedger.lean) and
[clock readout algebra](../../Lean/Time/ClockReadout.lean) preserve the
distinction between event order, modular parameter and calibrated comparison.

The extension of the existing block-elimination identity is from first moments
to record-conditioned quantum maps and their complete output states. Its proof
is the finite substitution above, not a new Lean theorem or an invocation of a
Gaussian theorem on this non-Gaussian process. The four-mark event closure is
an additional source-specific simplification checked on the complete operator
space. This discharges the finite joint-law correspondence for this declared
experiment. Compatible refinement, controlled continuum limits, physical
calibration and relativistic spacetime recovery remain separate obligations.

```sh
python -m model.joint_clock_process --check model/joint_clock_process/results.json
python -m pytest -q -W error model/joint_clock_process
```

For the [complete local-history truncation bound](../native_action/README.md),
the observed region must contain all four active ports: the declared field
Hadamard and phase readout span both carriers. The locality theorem does not
treat those services as operations confined to one elementary observer.
