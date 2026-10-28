# Continuing exported memory with a transported reference

The actual [continuing C/E observer](../completion/OBSERVER.md) can retain the
reference metadata from the
[recorded-history transport construction](../native_modular_reference/history_transport.py)
through repeated FIFO exports. The additional metadata is an unordered
two-port core: 66 possible locations, encoded in 7 logical bits. It updates
from each committed native current before the older record is exported.
This adds a finite reference controller to the existing ideal instrument;
it does not execute another preparation or a new coherent gate.

Exact retained runs cover 24 model attempts for C, E and the delayed tagged
C/E mixture, with 23 exports and one retained record. Every outcome enters
the complete instrument. Both endpoint laws reach all 66 core locations by
attempt 4. The four clock words previously compiled for a two-attempt policy
therefore do not cover this continuing process. An uncovered request receives
`unavailable_service`; its inventory query leaves the source unchanged and
the repair history continues. It is not counted as a successful clock action.

## Contract and actual maps

The carrier is the invariant binary total-two shell, dimension 66, in the
existing native basis order. The total-two dimension-78 instrument is reused
by restriction; double occupations are outside this transport contract.
Allowed histories contain the existing repair outcomes, tagged C/E mixtures,
identity holds, classical record feedback and positive uniform exploration.
An arbitrary untouched reference system may be correlated with the source.
There are no coherent interventions in this continuation policy.

The independently specified source obligations are conserved total two,
binary occupancy, unchanged spectators and native floor/ceiling agreement at
the attempted seam. The record obligation is the full common committed tape,
including selector, seam, law, coin, record label and decoded current. The
exported prefix concatenated with the retained FIFO must equal that tape.
The metadata obligation is the reference extension determined by the complete
record history, whether or not older records remain in the working FIFO.

For an unordered core `A`, define the faithful full-shell reference

```text
sigma_A(x) = a[q_A(x)],       q_A(x) = sum_{p in A} x_p,
(a[0],a[1],a[2]) = (2/225, 1/36, 2/45).
```

At `A={0,1}` this is exactly the predecessor's locally prepared reference,
with no newly fitted weights. These references have unit trace on binary66.
The original preparation and its success/failure contract remain prerequisites;
the exact long runner starts from its ideal reference once.

The new kernel calls `completion.observer.branches` and preserves its branch
weights, memory writes, exports and semantic zero outcomes. On a nonzero
native current at seam `(a,b)`, binary occupancy means the source branch
exchanges those endpoint bits. The controller therefore updates `A` by the
same port transposition. A zero current, hold or timeout leaves `A` unchanged.
The update reads the native record, not the source configuration. It does not
physically route the ports or assert a graph automorphism.

For each supported basis transition `i -> j`, the executable certificate checks

```text
sigma_A(i) = sigma_A_after(j).
```

Consequently, for every real modular parameter `s`, every admitted branch
`K_b=sqrt(w_b) M_b`, and every untouched reference `R`,

```text
(sigma_A_after^(is) tensor I_R) (K_b tensor I_R)
    = (K_b tensor I_R) (sigma_A^(is) tensor I_R).
```

This identity holds on all matrix units and composes through any finite
admitted history. The scalar branch weights may depend on prior classical
records. Holds, selector probabilities and tagged law mixtures introduce no
new phase ratio. The actual reference corner for history `h` is
`tau_h = w_h P_range sigma_A_h P_range`. Its normalized form differs only by
a scalar on its support. The faithful full-shell extension `sigma_A_h` is
neither the full-shell posterior nor a replacement for the separately solved
stationary source/FIFO reference. No stationarity of these transported
extensions is asserted.

The all-core certificate checks 522,720 supported columns, implying the
identities on 21,985,920 matrix units, and independently checks 3,960 native
input-column completeness sums. Tests reconstruct every native outcome using
floor/ceiling arithmetic, exercise complex entangled source/reference inputs,
and compare a complete two-step tagged tree against the compressed marginal.
The noncommuting-reference countercontrol remains excluded by the same domain
premise as in the predecessor.

## What survives export, and what still costs storage

The source and last record alone cannot determine the transported reference.
For C, compare the complete native histories

```text
((0,5,coin=0,current=0), (6,7,coin=0,current=0))
((0,5,coin=0,current=1), (6,7,coin=0,current=0)).
```

After the first record exports, both histories have the same last record and
share the same output matrix unit at binary coordinates `(15,51)`. Their
transported cores are `{0,1}` and `{1,5}`, with reference ratios `8/25` and
`5/8` on that unit. Both histories have positive probability under the actual
exploration floor. E has the same ambiguity with first contrast labels `-1`
and `1`, and final contrast zero. Thus keeping only the last FIFO record
forgets a required reference distinction. The new core metadata retains it.

The core and the feedback sign have at most `66*2=132` states, within the
declared 256-state finite-controller boundary. Eight logical bits suffice
for their joint encoding. The FIFO is accounted separately: ten bits encode
an upper bound of 513 possible complete semantic records across the two
feedback proposals and exploration, including zero-probability labels and
hold. The sign can also be decoded from the retained record. These are logical
storage bounds; no native hardware implementation of the memory/controller
or its workspace is supplied here.

At attempt `n>0`, the one-slot export policy has one retained record and
`n-1` environmental history records. The full history therefore has growing
total storage. The clock inventory status is determined by the core and is
not an additional independent controller state. Retaining copies of queried
status histories would also require storage. Under the alternative halt
policy, a full FIFO produces the existing timeout and preserves source,
memory and core; no hidden erase or recovery restarts it.

The quantity `q_A(x)` is preserved when both source and metadata core are
transported. This is a joint source/frame quantity. It supplies no new
source-only conserved charge, localized topological defect or matter energy.

## The four-word menu and the longer runs

The predecessor's compiled clock inventory at modular parameter `s=4` is
derived from its actual two-attempt branches. Its four cores are `{0,1}`,
`{0,8}`, `{1,5}` and `{5,8}`.
A single exploration branch at seam `(0,7)`, C coin zero and current one,
has probability `53/21600` in the endpoint campaign and moves the core to
`{1,7}`. That outcome is retained and its clock service is explicitly
unavailable. The diagnostic does not silently compile, substitute or certify
another word. Even an available inventory entry only identifies a predecessor
word; the present continuation runner does not execute clock gates.

The long runner evolves exact integer counts for the **diagonal marginal**
on source, core and last-current sign. It groups full records only after
their native actions, because future control uses that sign. It traces
earlier exports in this marginal; `branches` and `replay` still specify the
complete instrument and reconstruct individual exported histories.

| Member | Attempt 24 repairs / holds (expectation) | Positive joint states | Probability of a core outside the four-word menu |
| --- | --- | --- | --- |
| C | 24 / 0 | 8514 | approximately 0.34244946 |
| E | 24 / 0 | 8514 | approximately 0.34244946 |
| Delayed tagged C/E | 18 / 6 | 8514 | approximately 0.27840256 |

The receipt retains exact rational probabilities, source hashes, distribution
hashes and intermediate attempts 0, 1, 2, 4, 8, 12 and 24. Each marginal has
exact mass one. C and E coincide only in this diagonal compressed marginal;
their complete instruments retain different record and coherence behavior.
Separate 24-attempt recorded-word replays check every compatible initial
source state and reconstruct all 23 exports, without re-preparation.

This joins continuing memory to the transported-reference metadata. Native
storage/export hardware, further clock compilation, autonomous clock waiting,
observer clock calibration and physical energy remain separate obligations.
The parameter here is a discrete model-attempt count, not observer duration.
No integrated CP-01 gate or physical comparison is accepted by this packet.

Run from the repository root:

```sh
python3 -m model.clock_defect_energy.continuation --check
python3 -m pytest -q -W error model/clock_defect_energy/test_continuation.py
```
