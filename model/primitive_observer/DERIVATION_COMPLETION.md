# Native fixed-seam repair and entropy completion

The native binary C/E channels have the same entropy-completion endpoint and
different finite-event observation laws. Their difference has an exact,
reference-uniform norm. This attaches the existing finite completion theorem
to the actual native transition rule and gives a bounded computational
shortcut. It does not select that rule as physical dynamics.

The parameter `n` below is the number of native events on the **fixed seam
(0,1)**. It is neither host execution time nor calibrated physical duration.
Every event's coin and repair record is forgotten. Retained-record histories,
feedback and changing seams are different channels and cannot use this
compression. The executable API rejects requests for those uses.

## Domain and independent retained algebra

Use the existing 66-dimensional invariant binary shell: twelve binary loads
with total two. On the fixed seam, its Hilbert space decomposes as

```
H = H_equal direct-sum (C^2_orientation tensor C^10_spectator).
dim(H_equal) = choose(10,2)+1 = 46.
B = M46 direct-sum (I2 tensor M10).
```

The equal block consists of states with both seam endpoints empty, or both
occupied. The unequal block has one particle on a seam endpoint and its
other particle at one of ten spectator sites. B retains every equal-block
observable and every spectator observable independent of endpoint
orientation. Choosing this retained algebra, the binary domain and the
reference `I66/66` are premises. B is specified by these observables before
identifying a native repair channel with its completion.

For any matrix A, define E to keep the entire equal block, erase cross-block
entries and replace the unequal block by

```
(I2/2) tensor Tr_orientation(A_unequal).
```

This is the tracial conditional expectation onto B, both on observables and
states. Its Kraus representation consists of the equal-block projection and
the four operators `(1/sqrt(2)) |a><b| tensor I10` on the unequal block.
Their adjoint products sum to the identity. The map is completely positive,
unital, trace preserving, Hilbert-Schmidt self-adjoint and idempotent; it
fixes B pointwise. The tracial reference has trivial modular flow, so B is
modular invariant. The existing
[finite entropy-completion theorem](../../code/quantum_information/ALGEBRAS_AND_REPAIR.md#the-complete-finite-completion-theorem)
therefore identifies E as the unique relative-entropy completion of **all**
B data relative to that reference, including boundary states.

## Attachment to actual native C/E channels

For each input configuration x and fair coin c, the native primitive places
the endpoint loads at the floor and ceiling of their mean, in the order
specified by c. Its current is `j=x0-y0`. C groups amplitudes by the complete
record `(c,j)`; E groups by `(c,x0-x1)`. The record labels are then traced
out, producing the two channels studied here.

On the binary shell, the native E channel is exactly the independent block
expectation above. The native C channel additionally retains one half of
each equal/unequal cross-block matrix unit. Write

```
Q(A) = P_equal A P_unequal + P_unequal A P_equal.
C = E + Q/2.
E^2 = E, Q^2 = Q, EQ = QE = 0.
```

`completion.audit()` constructs the native record sums using exact rational
coin weights and compares them with the block formulas on **all 4,356 matrix
units**. It independently executes a second native C event on every unit.
Complex linearity extends these identities to the entire operator algebra;
tensoring with the reference identity gives the same identities for any
finite reference. This is more than a check on diagonal preparation states.

The projection ranges have dimensions

```
rank(E) = 46^2 + 10^2 = 2216,
rank(Q) = 2*46*20 = 1840,
dim(ker(E) intersect ker(Q)) = 66^2-2216-1840 = 300.
```

Thus C has eigenvalues `1`, `1/2`, `0` with those multiplicities, and

```
C^n = E + 2^(-n) Q,                 integer n >= 1.
```

Both laws fix B and have the same completion endpoint. Their finite-event
laws differ. Completion of a specified retained algebra does not select
the speed or complete recorded instrument of a repair process.

## Exact complete-input distance and its observed witness

Let `Z=P_equal-P_unequal`. It is unitary, and
`Q=(id-Ad_Z)/2`. Every channel has diamond norm one, so `||Q||_diamond<=1`.
Choose unit vectors a and u in the two different blocks and prepare
`(|a>+|u>)/sqrt(2)`. Q maps its density to a two-dimensional off-diagonal
matrix with eigenvalues `+1/2,-1/2`, whose trace norm is one. This attains
the bound, without any reference or postselection. Therefore

```
(1/2) ||C^n-E||_diamond = 2^(-n-1).
```

The diamond upper bound includes arbitrary reference correlations. After
either channel, any common subsequent CPTP processing and complete readout
can only reduce the resulting trace distance and classical total variation.
This statement compares the terminal forgotten-record channels; it does
not supply the eliminated event records to subsequent processing.

For an explicit native witness, a occupies ports `(2,3)`, u occupies
`(0,2)`, and v occupies `(1,2)`. Prepare the same a/u superposition and read
the complete binary instrument with effects `|a+u><a+u|/2` and its full
complement. The plus probabilities are

```
E:   3/8,
C^n: 3/8 + 2^(-n-1).
```

The minus probability is always retained. At one event the complete laws
are `(3/8,5/8)` and `(5/8,3/8)`, with total variation `1/4`. This is an
internal prediction under the specified preparation and measurement; no
experimental observations or calibrated apparatus enter the attachment.

## Executable shortcut and boundaries

`apply_compressed` applies the exact formula to sparse Gaussian-rational
operators on system times reference, using the index `system*r+reference`.
It retains the complete terminal operator and supports unnormalized,
non-Hermitian inputs as a linear-map interface. It does not normalize a
state or discard any terminal reference component.

For the generic floating channel executor, the same map has the positive
mixture realization

```
C^n = (1-2^(1-n))*E + 2^(1-n)*C,    n >= 1.
```

`binary_kraus` builds matrices from the native transition records;
`compressed_kraus` returns their square-root-weighted mixture. It rejects
binary64 underflow of a positive mixture weight. These floating adapters
are numerical executions, not outward roundoff certificates, emitted native
pulse words or physical service-error bounds. The computational shortcut
does not assert that n physical repairs can occur in one physical event.

The full 78-dimensional total-two carrier is a necessary negative control.
Both seam-concentrated states `(2,0,...)` and `(0,2,...)` map to `(1,1,...)`.
Consequently the full E channel maps the identity's diagonal values at
those three addresses from `(1,1,1)` to `(0,0,3)`. It remains trace preserving
but is not unital. Its action cannot be identified with the binary tracial
conditional expectation on that larger domain.

The additional contribution is the exact native attachment, finite-event
decomposition, saturated complete-input distance and executable compression.
The conditional-expectation theorem and channel-norm arguments are reused
with their domains stated. None of these results selects the retained
algebra, reference, interface or physical repair law from unrestricted OPH
candidates, and none establishes a physical clock or a field Hamiltonian.
