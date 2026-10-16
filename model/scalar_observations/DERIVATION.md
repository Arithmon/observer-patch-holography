# One physical detector across an interacting spatial family

The spatial and occupation comparisons have distinct scopes. The spatial
[receipt](results.json) bounds the fixed-band spatial discretization.
The [occupation derivation](OCCUPATION_DERIVATION.md) and its
[receipt](occupation_results.json) certify the principal smooth history and
an unrestricted interaction-dependent record event. A spatial bound alone
does not control occupation error.

This construction fixes the physical region, field band, initial preparation,
detector smearing, feedback and Hamiltonian model-time intervals before
comparing spatial grids. It gives an interacting comparison of the
**complete observation history** within that common domain. It reuses the
twelve-port source compiler and the complete smooth-instrument method in
[the common-source construction](../common_source/SMOOTH_REFINEMENT.md).

There are two different mathematical comparisons. Finite occupation executions
have an all-input spatial bound and a temporal splitting bound. Separately,
seven unrestricted oscillators have a finite-energy spatial bound through the
same physical detector and its feedback. The latter does **not** bound the
occupation error of the former. Removing the fixed field band, identifying a
physical clock or selecting this action from native repair is not established.

## 1. Fix the physical observation before changing the grid

The supplied region is the torus `[0,2*pi)^3`, of volume `V=(2*pi)^3`.
Keep its seven real orthonormal field modes,

    phi(x) = V^(-1/2) [Q0 + sqrt(2) sum_i
                      (Qcos_i cos(x_i) + Qsin_i sin(x_i))].

Thus the fixed physical smearing `f(x)=sqrt(2/V) cos(x_1)` gives
`integral f(x) phi(x) dx = Qcos_x`. It denotes the same observable on every
grid. This differs from the earlier occupation-four witness, which was
designed to detect a quartic matrix element in a selected oscillator basis.
The new readout is still a supplied apparatus choice, not a measured detector
of nature or a derivation of an observer's accessible algebra.

With `r_i^2=Qcos_i^2+Qsin_i^2`, direct spatial integration gives

    W(Q) = V integral phi(x)^4 dx
         = Q0^4 + 6 Q0^2 sum_i r_i^2
           + (3/2) sum_i r_i^4 + 6 sum_(i<j) r_i^2 r_j^2.

The spatial sampling grid has `n>=5` points in each direction and spacing
`a=2*pi/n`. Every quartic trigonometric monomial has frequency at most four
in each coordinate. Its discrete mean therefore equals its torus integral;
none of its nonconstant Fourier modes aliases to a constant. The mass and
kinetic inner products are exact by the same argument. Only the supplied
nearest-neighbor gradient changes:

    kappa_n = [sin(pi/n)/(pi/n)]^2,
    H_n = (1/2) sum_j P_j^2 + (m^2/2) sum_j Q_j^2
          + (kappa_n/2) sum_(j=1..6) Q_j^2 + (g/24) W(Q),
    [Q_i,P_j]=i delta_ij,     g=lambda/V.

`H_infinity` sets `kappa=1` on these **same seven modes**. The executed
parameters are `m=g=1`, unit action constant and gradient speed. Consequently
the physical quartic coefficient in this convention is `lambda=V`; it is
not a weak-coupling assertion. The coefficient `m` is a bare positive mass
parameter, not an inferred interacting pole mass. There is no normal ordering
or hidden mass counterterm.

The band is fixed while the quadrature and finite-difference grid refine.
There is no addition of high-momentum modes, ultraviolet limit or changing
physical box. Projection onto this band does not preserve exact spacelike
locality. Agreement here therefore establishes neither a full local
relativistic QFT nor its vacuum or scattering theory. Its significance is
that an interacting field description and an actual instrument can now be
compared on the same fixed observational domain.

## 2. Preparation, measurements and continuing feedback

The preparation is the product blank in the oscillator basis `Omega=2`.
It is a fixed state on seven modes, not the interacting ground state.
A signed momentum kick `exp(i k0 Qcos_x)` prepares the intervention. Its
sign can depend on the retained first native repair current. After the first
Hamiltonian interval, record both outcomes of

    K0 = cos((eta Qcos_x + theta)/2),
    K1 = sin((eta Qcos_x + theta)/2).

Their squared effects sum to identity. The first bit may trigger
`exp(i k Qcos_x)`. The disturbed state continues through a second Hamiltonian
interval and a second read. All outcomes, surviving quantum states and native
records remain in the comparison; no branch is normalized or postselected.
The final bit controls the common native continuation.
The executor keeps the disturbed vectors while continuing each history. The
receipt stores its classical history probabilities, not every terminal quantum
vector; the complete-output statement is the stronger mathematical bound.

The blank characteristic function supplies the declared zero-duration
calibration comparison,

    <exp(i eta Q)> = exp(-eta^2/8),
    Pr(bit=0) = [1+cos(theta) exp(-eta^2/8)]/2.

The initial momentum kick commutes with this detector, so it does not change
that zero-duration calibration. These are analytic properties of the supplied
blank and apparatus. Simulated calibration counts do not turn them into
independent physical observations or determine an interacting mass.

All times below are the supplied Hamiltonian **model coordinate**. Source
pulse counts measure requested control resources; host runtime measures the
cost of simulation. Neither calibrates a physical duration.

## 3. Correct finite occupation source and its spatial bound

Let `Pi_d` retain a rectangular set of occupations in the fixed `Omega=2`
basis, and put `H_(n,d)=Pi_d H_n Pi_d`. Full `Q^2`, `P^2` and `Q^4` are
formed before compression. In particular `Pi Q^4 Pi != (Pi Q Pi)^4`.

The finite detector instead is deliberately defined by `cos((eta Q_d+theta)/2)`
and `sin((eta Q_d+theta)/2)`, where `Q_d=Pi Q Pi`. Likewise the finite kicks
are functions of `Q_d`. These are valid complete finite instruments; they
are not asserted to equal compressions of their unrestricted counterparts.
At each fixed occupation code they are exactly the same on all spatial grids.

The source words execute the seven onsite terms and the twenty-one pair
terms. An onsite quartic coefficient is `g/24` for the zero mode and `g/16`
for each nonzero mode. The coefficient of `Q_i^2 Q_j^2` is `g/8` for the
sine/cosine pair on one axis and `g/4` for all other pairs. Every pair
commutes with every other pair, including pairs sharing a mode, because
every appearance uses the same compressed full `Q_j^2`.

Write `X_j=Pi_j Q_j^2 Pi_j`. The only difference between finite grid and
exact-gradient Hamiltonians is

    H_(n,d)-H_(infinity,d) = (kappa_n-1)/2 sum_(j=1..6) X_j.

The unitary Duhamel identity gives

    ||exp(-it H_(n,d))-exp(-it H_(infinity,d))||
      <= |t| |1-kappa_n|/2 sum_(j=1..6) ||X_j||.

Common instrument dilations are isometries; common record-conditioned
feedback and native continuation are channels. Telescoping the two evolution
intervals therefore bounds the complete recorded-output trace distance and
classical history TV distance by

    min(1, (t1+t2) |1-kappa_n|/2 sum_(j=1..6) ||X_j||).

This finite-code bound is uniform over all initial states and arbitrary
untouched references. It does not depend on the number of record branches.
Trace distance and classical TV both use the factor one half.

`bounds.py` encloses each Hermitian or anti-Hermitian norm by absolute row sums after exact
conversion from factorial-weighted Bargmann monomials to normalized Fock
coordinates. Square roots are rounded upward using integer arithmetic.
For `n>=5`, the elementary sine bound and `pi<22/7` give

    0 <= 1-kappa_n <= pi^2/(3n^2) < 484/(147n^2).

Thus no floating evaluation of sine or pi enters the certificate.

## 4. Temporal splitting with exact local commutator bounds

Put `A=sum_i a_i` for the onsite terms and
`B=sum_(i<j) c_ij X_i X_j` for the pair terms. One symmetric step is
`exp(-ihA/2) exp(-ihB) exp(-ihA/2)`. Centering each `a_i` by its blank
expectation removes the same scalar phase from exact and split evolutions.
The code independently retains all scalar phases in the actual source word.

A conservative fallback uses `S>=sum_i||a_i-centered_i||+sum c_ij||X_i||||X_j||`.
Exact and product evolutions agree through second order, and each third
derivative has norm at most `S^3`. For `N` steps the vector/operator error is
at most `S^3 |t|^3/(3N^2)`.

A substantially smaller bound uses local commutators. Define

    C_i=[a_i,X_i],  D_i=[a_i,C_i],  E_i=[X_i,C_i],
    R_i=sum_(j!=i) c_ij X_j.

Commutation between different modes gives the exact identities

    [A,[A,B]] = sum_i D_i R_i + sum_(i,j) c_ij C_i C_j,
    [B,[A,B]] = sum_i E_i R_i^2.

The corresponding triangle and row-sum bounds only require local matrices.
`C_i` and `E_i` do not depend on `kappa`, since the varying onsite term is a
multiple of `X_i`. `D_i` is affine in `kappa`, so convexity bounds its norm
on the entire interval `[0,1]` by its two endpoint bounds. The same endpoint
argument applies to the centered onsite norm.

For bounded self-adjoint `A,B`, differentiating the symmetric product and
subtracting its exact evolution equation gives a defect bounded at duration
`s` by

    s^2 ||[A,[A,B]]||/8 + s^2 ||[B,[A,B]]||/4.

To obtain this bound, express each commutator with an exponential as an
integral of conjugated `[A,B]`; the first-order terms cancel. The remaining
differences of conjugations are bounded by the appropriate nested
commutator times their integration interval. Integrating the defect and
telescoping the `N` unitary steps gives

    error <= |t|^3/(24 N^2)
             (||[A,[A,B]]|| + 2||[B,[A,B]]||).

`product_error` takes the smaller of this bound and the cubic fallback.
The errors for the two evolution intervals add through the full instrument,
then are clipped at one for recorded-output trace distance. At occupation
two, all compressed onsite and squared-coordinate matrices are diagonal;
the temporal split is exact. This degeneracy is not a validation of the
occupation-two code as a physical scalar approximation.

The spatial and splitting bounds can be added to compare ideal source
execution with the exact-gradient **finite-code** target. Floating compiler
roundoff, matrix contraction error, stabilization, analogue precision and
physical preparation/storage/readout errors are separate. Agreement with a
floating sparse exponential is a numerical diagnostic, not an outward
certificate for those errors.

## 5. An unrestricted interacting spatial theorem through both reads

For `g>0`, the nonnegative confining polynomial potential defines the
Friedrichs Hamiltonian on `L2(R^7)`. Schwartz functions give a common core
for the displayed calculations; the energy estimates extend to its closed
quadratic-form domain. The crucial polynomial inequality is

    W(Q) >= (sum_(j=0..6) Q_j^2)^2.

Indeed, subtracting the right side leaves `4 Q0^2 sum r_i^2`,
`(1/2) sum r_i^4` and `4 sum_(i<j) r_i^2 r_j^2`, all nonnegative.
Since the kinetic and quadratic potentials are also nonnegative,

    ||sum_(j=1..6) Q_j^2 psi||^2
       <= (24/g) <psi,H_infinity psi>.

For a normalized finite-energy state with energy `E`, Duhamel comparison
against the exact-gradient evolution therefore gives

    ||(exp(-itH_n)-exp(-itH_infinity))psi||
       <= |1-kappa_n| |t| sqrt(24E/g)/2.

The potential difference is bounded from the energy form domain into `L2`
and is infinitesimally `H_infinity`-bounded. The Kato--Rellich perturbation
argument therefore gives a common Hamiltonian operator domain. Duhamel's
identity first holds there and extends by finite-energy approximation;
`H_infinity` evolution conserves that energy.
This argument needs a first energy moment, not a Gaussian interacting state
or a graph-norm estimate. It also applies to vector-valued wavefunctions
carrying records and untouched references.

For the actual smooth detector, `sum_b K_b^2=1`,
`sum_b K_b K_b'=0` and `sum_b (K_b')^2=eta^2/4`. Applying the product rule
to the kinetic quadratic form gives the exact unconditional identity

    sum_b <K_b psi,H_infinity K_b psi>
       = <psi,H_infinity psi> + eta^2 ||psi||^2/8.

The multiplication operators preserve the potential form domain and their
bounded derivatives preserve the kinetic form domain. Each branch has
finite energy. This is an identity after summing **both** unnormalized
branches, not an energy increment asserted on every conditioned branch.

A kick `exp(i k Q)` changes momentum by `k`. For a direct sum of retained
branches where any subset receives that kick, Cauchy--Schwarz and
`<P^2>/2<=E` give

    E_after_feedback <= E + |k| sqrt(2E) + k^2/2
                      = (sqrt(E)+|k|/sqrt(2))^2.

The bound does not discard the untriggered branches or require a lower
bound on a trigger probability. For the product `Omega=2` blank followed by
one signed initial kick, the exact-gradient initial energy is exactly

    E0 = 7/2 + (7m^2+6)/8 + 49g/128 + k0^2/2.

Its momentum mean is zero before the kick. The quartic expectation used
here is `<W>=147/16`; no interacting-vacuum approximation is made. Set

    E1 = (sqrt(E0+eta^2/8)+|k|/sqrt(2))^2.

Telescoping the two unitary intervals around the common detector and
conditional feedback now proves the complete unrestricted history bound

    TV <= trace_distance <= min(1,
       |1-kappa_n|/2 [t1 sqrt(24E0/g)+t2 sqrt(24E1/g)]).

The second read and subsequent native feedback are common channels, so they
add no comparison error. The theorem retains every disturbed quantum output
and record; the first detector's backreaction is explicitly included in
`E1`. The initial kick sign and detector phase do not alter these bounds.
Directed rational roots and the rational `kappa` envelope implement the
theorem in `unrestricted_spatial_history`.

This supplies an `O(n^-2)` interacting spatial observation comparison for
the seven **unrestricted** modes, at fixed physical smearing and model times.
Its positivity argument requires `g>0`; the interaction-deletion control
`g=0` is covered by the separate finite-code bound, not by dividing this
theorem by zero.

## 6. The remaining occupation and physical obligations

The unrestricted theorem and finite source execution have matching physical
definitions, but no quantified occupation comparison between them has been
proved in this packet. Their initial kicks, smooth detectors and conditional
kicks create infinite occupation support in the unrestricted model. Agreement
between occupation-two and occupation-four numerical outputs, and the
zero-duration blank calibration, do not supply that missing bound.

The finite-code grid and temporal certificates must consequently retain a
finite-code target. The unrestricted spatial certificate must remain a
separate theorem. The earlier eight-site occupation certificate cannot be
transferred: it used another action, preparation and detector.

An additional occupation theorem, wider physical momentum bands, compatible
state preparation and calibrated operational time are subsequent tasks.
Calibration observations and later comparison observations must have explicit
roles. The current numerical choices are exploratory; their response to the
supplied interaction or feedback is not a preregistered prediction or an
independent reproduction of laboratory data. Selecting the source law and
matching actual observed quantum-field, classical and gravitational regimes
remain the larger OPH construction obligations.
