# Independent reconstruction of the observer instrument

The [written observer construction](PROOF.md) and its
[native-event ancestry](../native_memory_clock/PROOF.md) define the audited
law. The [contract](observer_independent_contract.json) pins the upstream
revision, source/input bytes, admitted inputs, domains and stopping rules.
The source85 witness, proof, protocol and thresholds were read before this
implementation. Independence means a separate derivation and implementation.
It supplies no blind test, external acceptance or physical measurement.

## Native rule and total gates

The edge order is `(01,02,03,12,13,23)`. Each edge contributes `+1` at its
lower vertex and `-1` at its upper vertex. Its flux is the odd-address bit
minus the even-address bit. Addresses in each pair point forward, then
backward. The loads are that incidence product plus the declared background.
For example source1 at background0 has loads `(-1,1,0,0)`; address0 changes
source1 to source0, and reversing the addressed convention changes it to
source3. The strict directed-current tests retain no-operations. The uphill
control negates the load difference before applying those tests.

For fixed unaddressed bits `s`, let `g(s,b)` be the addressed successor bit.
The total map sends `(s,b,r)` to `(s,g(s,b) XOR r,b)`. Given `(s,c,d)`, restore
the input source bit to `d`, evaluate `g(s,d)`, then recover memory as
`c XOR g(s,d)`. Substitution proves both inverse identities for every bit
pair. Blank memory is required to interpret the result as a native event.

The independent program interprets the four instruction lists in the proof.
It retains five source buffers, four first-address registers, three couriers,
two two-bit pointers, two old-bit memories and the action address. It
executes each XOR, modular translation and swap, including blank complements.
The inverse reverses the instruction list, subtracts translations and uses
the independently derived event inverse. Each instruction is a total
permutation, which proves full-carrier unitarity algebraically. The finite
checks cover every prepared prefix and a bounded deterministic dirty domain;
they are not an enumeration of the enormous full carrier.

On the native blank preparation, gate1 leaves `y=f_a(x)` and the old bit
`x_a`. Gate2 writes `r=y_a` and validity into A and the courier. Gate3 moves
that courier. Gate4 gives B the same pair, adds `6r` to the initially zero
action address and performs `z=f_(6r)(y)`, retaining `y_(6r)`. The implementation
checks this independently evaluated endpoint against the literal instruction
execution. The first address is a coherent uniform twelve-state register;
the second address is adaptive and has no independent uniform preparation.

## Sparse operators and arbitrary references

For prefix `j` let `d_j(x,a)` denote its complete output register record.
Its isometry column is

```text
E_j |x> = (1/sqrt(12)) sum_a |d_j(x,a)>.
E_(j,p) |x> = (1/sqrt(12)) sum_(a:pointer(d_j(x,a))=p) |d_j(x,a)>.
```

All 4096 input columns, twelve addresses, five prefixes, two backgrounds and
five declared variants are enumerated with integer records. An injective
mixed-radix encoding covers every register. Distinctness of all 49,152
records per prefix proves disjoint column support; each column has twelve
entries of squared weight exactly `1/12`. Thus `E_j* E_j=I` exactly.
The public pointer partitions are exhaustive and disjoint, giving
`sum_p E_(j,p)* E_(j,p)=I`.

The adapter compares every complete output record and public pointer
partition with the producer, using the same symbolic positive amplitude
`1/sqrt(12)`. This is equality of the full public isometry columns, including
every private register. It implies equality of every unnormalized block
`(E_(j,p) tensor I_R) rho (E_(j,p)* tensor I_R)` for arbitrary input density
and arbitrary untouched reference R. No finite selection of input vectors
or dense Choi matrix is used to establish that statement. Each column
comparison concerns the executable basis-permutation law; it does not test
an alternative phase-decorated architecture.

The public block sums private histories before taking the outer product.
The complete computational readout measures those histories separately.
These are different instruments. The declared exploratory diagnostic uses
source0 and source3 with reference amplitudes `(1,i)/2` and `(1,1)/2`. Its
coefficient determinant is `(1-i)/4`, so it is entangled. The diagnostic
verifies normalization and the untouched reference and finds a nonzero
coherence between distinct archives within one public outcome. Archive
dephasing preserves public probabilities and changes that retained density.

## Clock and unfinished outputs

The scalar clock has hopping `1/2`. If the edge instructions are total
unitaries, the block-diagonal history gauge conjugates the joint Hamiltonian
to the scalar path matrix tensor identity. This proves program-independent
clock weights and the full sine spectrum `cos(k*pi/(N+1))`, each with data
multiplicity `2^72 * 12^5`. The proof applies to each declared control.
Clock arrival alone supplies no receiver validity under disabled
communication or publication. The independent enumeration derives those
pointer flags from the executed instructions for each control.

The packet calculation sums binomial integers to recover normalization,
position moments, adjacent overlaps and zero symmetrized covariance. It
then forms the one-sided variance bound and verifies the signs giving
monotonic decrease over the whole readiness window. The finite-boundary
tail uses the independently computed distance and rational geometric
enclosure. No claimed constant or result array supplies these values.

The ideal failure bound is `226536448/278318778205` plus the finite-boundary
tail. A computed dyadic enclosure bounds that tail by `2^-731`. This is
tighter than the elementary bound in the construction proof. An independent
integer-square-root interval bounds the binary packet vector error; twice
that error bounds normalized preparation error. The combined bound is below
`1/1000`. This representation bound does not certify numerical propagation
roundoff. The adapter checks both signed binary preparations against their
declared phases and support.

For native blank preparation, B.valid is true exactly at prefix4. The
clock-independent isometries and scalar probabilities therefore give

```text
Phi_t = p_t Phi_complete + (1-p_t) Phi_unfinished,
p_t = sum_(n>=p+4) |alpha_n(t)|^2.
```

The unfinished map includes all four unfinished prefixes with their
relative weights when `p_t<1`. If `p_t=1`, that map can be any channel because
its coefficient vanishes. Completed and unfinished outputs have orthogonal
public validity flags. Convexity bounds half-diamond distance to the ideal
completed instrument by `1-p_t`; measuring the flag attains that value.
The argument holds with any untouched reference and any input. It proves
the uniform single-terminal-snapshot claim, not a joint trajectory or
permanent record probability.

The source85 feedback and uphill thresholds concern the conditional
completed S4 source law. Their exact TVs are `1/2` and `1`. For the
unconditional S4 marginal, every unfinished prefix has S4 blank under all
five programs, and the clock weights agree. Consequently the TVs are
`p_t/2` and `p_t`. This conversion specifies the buffer and readout: it is
not a claim about the complete retained quantum-channel distance. The
conditional completed law is reported together with its completion weight.

The separate N48 diagnostic builds literal forward and adjoint columns
using the independent instructions, verifies closure and Hermiticity,
compares producer columns and propagates the complex entangled input.
Scalar 48-by-48 diagonalization provides a separate clock calculation.
Unnormalized pointer blocks are compared with the five-prefix mixture,
retaining every unfinished weight. Its numerical residuals and completion
weight apply only to N48. No numerical result is transferred to N4096.

## Coverage and reproduction

The producer oracle imports the production native rule; the existing quantum
composition tests use eight input sources. This reconstruction adds an
import-isolated native rule and full 4096-source public isometry comparison.
It leaves production, protocol, provenance, frozen evidence and thresholds
unchanged. Directory discovery in the observer CI lane collects the focused
test file without a lane change.

From the repository root, using the installed numerical dependencies:

```sh
python3 -B model/native_observer_archive/independent.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B -m model.native_observer_archive.compare_independent
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B -O -m unittest model.native_observer_archive.test_independent
```

The standalone command blocks production imports and denies every repository
file read except the contract and protocol, including attempted receipt
reads. It moves no tracked file. It requires only the Python standard
library. The adapter validates the audited code/input byte bindings before
importing production, using the shared pin-scope policy. Proof digests in the
contract document the audited snapshot; editorial files are exempt from live
custody comparison. The exact checker retains integer records and rational weights;
the composition diagnostic uses the declared `2e-10` numerical gate.

Mandatory mutations reach register-column equality, inverse identity,
retained coherent-block equality and complete-law accounting predicates.
Hash validation is separate custody coverage. Communication and publication
controls exclude invalid receiver records, and the feedback and uphill
controls discriminate an actual later source law on source85. No difference
on every input is claimed; source0 is a permitted nondiscriminating input.

This validates the finite supplied construction. It establishes no physical
time or energy normalization, particle dictionary, continuum limit,
permanent autonomous archive, CCF-LT-2 equivalence or OpenWave physical row.
The signed-graph energy-clock diagnostic has a separate source and supplies
no premise here.
