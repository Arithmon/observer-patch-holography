# A finite source implementation of the two-site scalar instrument

This construction emits and executes gates from the existing twelve-port
source response. It connects prepared oscillator registers, field evolution,
a binary pointer and the existing native record controller. Its Hamiltonian,
site interpretation, coherent code service, resets and model time are declared
design inputs. It does not require that this design be uniquely selected.

The independent compatibility condition is that a decoded source-code
operation agrees with its stated finite oscillator or pointer operation,
including relative phase and every readout outcome. The schedule is a fixed
finite sequence of proper-code services and pulses. Its legitimate result is
a complete observation-and-feedback history, not an instantaneous fixed
field state. Model time in the equations below is the supplied Hamiltonian
parameter; source pulse counts and host execution time are different quantities.

## 1. Physical target and finite-code target

The two-site infinite-oscillator target is

```
kappa = 9/32,
H = (P0^2 + P1^2 + Q0^2 + Q1^2)/2 + kappa*(Q0-Q1)^2/2.
```

The symmetric and antisymmetric normal coordinates have frequencies 1 and
5/4. This is a chosen finite action, not a spatial refinement of the golden
scalar family. Its independent physical identification remains a further task.

Let `Pi_d` project each oscillator onto number states `0,...,d-1`, and let
`Q_d` and `P_d` be the truncated coordinate matrices. The executed Hamiltonian
uses the **compressed physical squares**:

```
Q2_projected = Q_d @ Q_d + (d/2)*|d-1><d-1|,
P2_projected = P_d @ P_d + (d/2)*|d-1><d-1|,
H_on = P2_projected/2 + (1+kappa)*Q2_projected/2,
H_d = H_on tensor I + I tensor H_on - kappa*Q_d tensor Q_d.
```

These identities follow by retaining the single excursion from level `d-1`
to level `d` and back. Therefore `H_d = Pi_d H Pi_d` on the finite code.
Deleting the top-level term changes the Hamiltonian. In particular, at
`d=2`, both `Q_d^2` and `P_d^2` equal `I/2`; that uncorrected onsite
Hamiltonian has no oscillator gap. Also

```
[Q_d,P_d] = i*(I-d*|d-1><d-1|).
```

Finite-code evolution must consequently be checked as a quantum matrix
operation. Treating it as an exact canonical Gaussian shear would be wrong.
The boundary projector does not commute with `Q_d`, so the corrected
potential cannot be grouped into one exactly commuting quadrature layer.

## 2. Preparation from blanks and its state bound

Only the antisymmetric normal mode requires squeezing. With `tau=1/9`, its
infinite squeezed vacuum has amplitudes

```
a_(2n) = sqrt(4*sqrt(5)/9) * (-tau)^n * sqrt(binomial(2n,n))/2^n,
a_(2n+1) = 0.
```

For `K=d-2=2J`, keep occupations through `K`. The retained probability is

```
c_J = (4*sqrt(5)/9) * sum_(n=0)^J binomial(2n,n)/(4^n * 81^n).
```

The site-basis vector follows by expanding
`(a0_dagger-a1_dagger)^(2n)/sqrt(2)^(2n)`. Its coefficient on `|k,2n-k>` is

```
a_(2n) * (-1)^(2n-k) * sqrt(binomial(2n,k))/2^n / sqrt(c_J).
```

Every occupied site level is at most `K<d`. A number-preserving normal-mode
rotation therefore embeds this entire prepared sector exactly in the chosen
onsite code. The program compiles the resulting real normalized site vector
into a deterministic source state-preparation word from reset blanks; an
eigenvector supplied as the executor's initial state would be a different test.

The full infinite-oscillator vacuum overlap is exactly `sqrt(c_J)`, giving
the pure-state trace distance

```
D(prepared, original_vacuum) = sqrt(1-c_J).
```

The recurrence
`p_(2n+2)=p_(2n)*(2n+1)/(81*(2n+2))` supports directed rational bounds on this
tail. A simple independent bound follows from
`binomial(2n,n)/4^n <= 1`:

```
1-c_J <= (1/81)^(J+1)/(1-1/81).
```

This is a preparation bound on a common infinite state space. It survives
any common subsequent channel by trace-distance contraction. It is not, by
itself, a bound comparing finite-code evolution with infinite field evolution.
Resetting the work registers preserves an untouched reference marginal, but
does not preserve entanglement previously held by the reset registers.

### Compiling larger signed states from the source controls

The eight-qubit two-site `d=16` register uses a binary amplitude tree. At
level `k`, the higher bits select a prefix `s`. Let `A_s0,A_s1` be the norms
of the target amplitudes beneath its two child branches. A rotation with
`theta_s=atan2(A_s1,A_s0)` divides the parent norm between them. At the last
level, use the two signed amplitudes themselves in `atan2`; this preserves
negative leaves and global phase. An empty prefix has angle zero.

The conditional rotation is compiled using

```
P_s = 2^(-m) * sum_z (-1)^(s dot z) Z^z,
exp(-i*Y_k*sum_s theta_s*P_s)
 = product_z exp(-i*c_z*Y_k*Z^z),
c_z = 2^(-m)*sum_s theta_s*(-1)^(s dot z).
```

Here `m` is the number of higher control bits. Every factor at this one
level commutes with the others. The target bit carries Y and the higher
bits carry the indicated Z factors, so these are the existing native Pauli
words, with their phases and proper-code services charged. Levels are
applied from the highest bit to the lowest. Induction on the prefix length
proves the resulting blank-input amplitudes equal the specified real target;
the emitted word also defines a unitary on all other inputs.

The full schedule includes every Walsh coefficient, including zero. A mask
with weight `w` uses `3+6*w` native pulses: two Y-basis changes, the central
Z rotation and forward/inverse native CNOT parity words. Summing over masks
and levels gives `3*((n-1)*2^n+1)` pulses for `n` qubits: 5,379 at `n=8`.
The existing Gray preparation words through six qubits are preserved.
The orthogonal-basis compiler's same Gray elimination extends to four
qubits for the local `d=16` basis changes.

## 3. Evolution, pointer and committed feedback

Write `A=H_on tensor I+I tensor H_on` and
`B=-kappa*Q_d tensor Q_d`. For a rational step `h`, execute

```
S(h) = exp(-i*h*A/2) exp(-i*h*B) exp(-i*h*A/2).
```

The two onsite factors commute with each other. Each is compiled as a
finite matrix operation; diagonalization provides control coefficients,
not an uncharged basis-change gate. The bond factor is diagonal in the
product quadrature basis, so its diagonal phases commute exactly. Repeating
64 or 128 steps reaches the same declared model time 1. Independent direct
finite exponentials test the finite evolution and its step refinement.
The separate rational certificate below bounds ideal finite split error;
measured agreement is not a certified floating-error or continuum bound.

For native current `j`, the preparation impulse is `exp(i*j*Q0_d)`.
Prepare the pointer in `|+>`, apply `exp(-i*Z_pointer*Q1_d/2)`, and measure
Pauli Y. The field Kraus operators for outcomes `y=+1,-1` are

```
K_y = (exp(-i*Q1_d/2) - i*y*exp(i*Q1_d/2))/2,
K_y^dagger K_y = (I+y*sin(Q1_d))/2,
sum_y K_y^dagger K_y = I.
```

Both terminal field branches and their Born weights belong to the result.
The committed pointer bit selects the second native seam. Subsequent
feedback preserves the earlier classical records; it need not preserve
coherence between distinct record values. The field is terminal after this
read, and a later trial uses fresh preparation.

Parity maps `j` to `-j` and changes `sin(Q1_d)` to its negative. The prepared
state, Hamiltonian and split factors preserve that symmetry. Thus
`Pr(y=+1|j)=1/2+j*r_d` for `j in {-1,0,+1}`. The existing 528-path native
enumeration then gives

```
E[j2] = 5*r_d/132,
E[final port0 load] = 1/6 - 5*r_d/132.
```

An independent donor bit preserves the fair marginal detector distribution
but removes its correlation with the first current, giving zero mean second
current. For two parity-symmetric responses `r` and `r'`, complete-history
total variation is exactly `5*abs(r-r')/33`. These are finite classical
postprocessing identities, not a physical continuum identification.

The ideal infinite-oscillator reference has detector variance `9/20` and
mean for unit impulse at model time `t`

```
mu(t) = (sin(t) - (4/5)*sin(5*t/4))/2,
Pr(y=+1|j) = (1+exp(-9/40)*sin(j*mu(t)))/2.
```

This independently computable reference is an informative comparison, not
a replacement for the emitted finite gate execution or its error budget.

## 4. Finite splitting and source-service error contract

For the two bounded finite matrices, the symmetric split has the uniform
operator-norm bound

```
||S(t/n)^n-exp(-i*t*(A+B))||
 <= (t^3/n^2)*(||[A,[A,B]]||/24 + ||[B,[A,B]]||/12).
```

Here is a direct defect proof of the constants. Set `X=-i*A`, `Y=-i*B`,
`C=[X,Y]`, and `T(s)=exp(s*X/2) exp(s*Y) exp(s*X/2)`. Writing
`Ad_U(C)=U C U^-1`, differentiation and cancellation of the first-order
commutator terms give

```
(T'(s)-(X+Y)*T(s))*T(s)^-1
 = Ad_exp(s*X/2)(
     integral_0^(s/2) (Ad_exp(-u*X)(C)-C) du
     - (1/2)*integral_0^s (Ad_exp(u*Y)(C)-C) du).
```

Unitary conjugation preserves norm, and integrating each inner derivative
therefore bounds this defect by
`s^2*||[X,C]||/8+s^2*||[Y,C]||/4`. Variation of constants integrates it
to the one-step constants `1/24` and `1/12`. Telescoping `n` unitary steps
gives the displayed bound. Tensoring a reference identity changes none of
these operator norms. The resulting bound also controls state trace
distance and any detector probability after a common preparation/readout.

`bounds.py` bounds each nested-commutator operator norm by its Frobenius
norm using exact rational arithmetic. In the unnormalized number basis
`(a0_dagger)^n0 (a1_dagger)^n1 |00>`, all entries of `A` and `B` are
rational: annihilation has entry `n`, creation has entry `1`, and their
quadratic coefficients are rational. With metric weights
`w_(n0,n1)=n0!*n1!`, a raw matrix `M` has

```
||M_normalized||_F^2 = sum_(row,column) M_(row,column)^2*w_row/w_column.
```

Both nested commutators are real in these coordinates. Integer square-root
rounding gives rational upper bounds on their norms. For `d=8`, the
squared Frobenius norms are respectively `123616882431/134217728` and
`23943536253/67108864`. At model time 1 the resulting finite split bound
is below `0.000694` for 64 steps and `0.000174` for 128 steps. These bounds
compare the ideal split and continuous evolution of the **same finite
compressed action**. Under their shared parity, the detector error further
transfers to classical native-history total variation at most `5/33` times
the bound and second-current mean error at most `5/132` times the bound.
Those smaller coefficients concern the classical native records, not the
complete terminal quantum field state.

The ideal source decomposition is an algebraic identity in finite
dimensions, including phase. For actual pulses with uniform operator
errors `eta_l`, plus proper-code service, reset and readout channel errors,
telescoping gives the sum of those errors as an upper bound on the final
state/instrument error, provided the service bounds hold with arbitrary
reference systems. The demonstration's floating residual tolerance tests
implementation consistency; it does not certify those physical errors.

The numerical executor first contracts each complete subroutine from the
full source exponentials, then takes a guarded polar factor of that executed
matrix before repetition. This is numerical stabilization of an exactly
unitary word, not an extra physical gate or an inserted target matrix.
Both raw and corrected residuals and the correction size are retained.
Large raw defects are rejected, and independent target/branch comparisons
remain unchanged. No physical failure probability or measured branch is
renormalized by this procedure. Its floating residuals do not replace a
certified roundoff budget.

## 5. Full occupation bound for the declared two-oscillator instrument

The following comparison covers the ideal continuous finite instrument:
normalized prepared vacuum, current impulse, evolution under the compressed
physical Hamiltonian, both pointer arms, readout and terminal field states.
It does not add an interacting field, a spatial continuum, laboratory
identification or physical source-service bounds. The Hamiltonian parameter
`t` is model time, restricted here to `0 <= t <= 1`; `j` is one of `-1,0,1`.

Let `Pi` project both site occupations onto `0,...,d-1`, and let
`N = n0+n1` refer to the original unit-frequency number basis. The comparison
uses the same two-site action and pointer specified above. Denote its infinite
vacuum by `Omega`. The common normal mode is unsqueezed; the difference mode
has squeezing parameter `tau=1/9`. Number-preserving normal-mode rotation
preserves `N`.

### Uniform number tails along the complete infinite history

During the current impulse and subsequent Hamiltonian evolution, the infinite
state is, up to a phase, `D(alpha) Omega`. The covariance is stationary under
this Hamiltonian. Conservation of the displacement's classical oscillator
energy and the lower frequency bound 1 imply

```
|alpha|^2 = (|q|^2+|p|^2)/2 <= j^2/2 = R0.
```

During the impulse, replace `j` by its intermediate value; the same bound
holds. At the end of model time `t`, the site's mean momentum is

```
p1(t) = j*(cos(t)-cos(5*t/4))/2,
|p1(t)| <= |j|*t/8.
```

This uses the global Lipschitz constant 1 of cosine. A pointer arm
`exp(+-i*s*Q1/2)`, `0 <= s <= 1`, adds momentum `+-s/2` on site 1, so every
state in either pointer arm obeys

```
|alpha_pointer|^2 <= R0 + |j|*t/16 + 1/8 = R1.
```

For a real number `1 < z < 9`, the exact Gaussian generating function is

```
<z^N> = sqrt(80/(81-z^2)) * exp((z-1) * (
    |alpha_plus|^2
    + 10/(9+z) * (Re alpha_minus)^2
    +  8/(9-z) * (Im alpha_minus)^2)).
```

One derivation uses the difference-mode Bargmann wavefunction

```
C * exp(-tau*w^2/2 + (alpha_minus+tau*conj(alpha_minus))*w),
C = (1-tau^2)^(1/4)
    * exp(-|alpha_minus|^2/2-tau*conj(alpha_minus)^2/2).
```

Replacing `w` by `sqrt(z)*w` computes the number-generating norm. The complex
Gaussian integral separates into real and imaginary directions and gives the
displayed coefficients. The common-mode coherent-state contribution is
`exp((z-1)*|alpha_plus|^2)`. Because `8/(9-z)` bounds all three coefficients,

```
<z^N> <= M(z,R)
       = sqrt(80/(81-z^2)) * exp((z-1)*8*R/(9-z)).
```

For `d >= 2`, `z >= 5/2`, and `k=0,1,2`, the sequence `n^k/z^n` decreases
for all integers `n >= d`: its successive ratio is at most
`(1+1/d)^2/z <= (9/4)/(5/2) < 1`. Therefore

```
E[N^k * 1_(N>=d)] <= d^k * M(z,R)/z^d.
```

`cutoff.py` evaluates finitely many rational choices `z=5/2,3,...,8`, each a
valid bound, and takes their minimum. Every root and exponential is rounded
upward by exact rational operations. The exponential uses a Taylor sum with
a positive geometric bound on its remaining tail; the ratio bound is
`x/(m+2) < 1` after term `m`. No measured endpoint population or numerical
quadrature enters these uniform bounds.

### Boundary residuals of the actual projected generators

For a vector `chi=(I-Pi) psi`, an outside number configuration has `N >= d`.
A coordinate can re-enter `Pi` only through its annihilation term, giving

```
||Pi Q0 (I-Pi) psi|| <= sqrt(E[N*1_(N>=d)]/2).
```

The pointer generator is `Q1/2`, so its bound has one additional factor
`1/2`. The Hamiltonian expands as

```
H = (1+kappa/2)*(N+1)
    + (kappa/4)*sum_i(a_i^2+a_i_dagger^2)
    - (kappa/2)*(a0+a0_dagger)*(a1+a1_dagger).
```

The diagonal term cannot cross the projection. Pure creation terms cannot
bring an outside configuration back inside. The remaining terms obey

```
||(a0^2+a1^2) chi|| <= ||N chi||,
||a0*a1 chi|| <= ||N chi||/2,
||(a0_dagger*a1+a1_dagger*a0) chi|| <= ||N chi||.
```

For the first bound use the number-preserving complex mode rotation
`b_plus=(a0+i*a1)/sqrt(2)`, `b_minus=(a0-i*a1)/sqrt(2)`: then
`a0^2+a1^2=2*b_plus*b_minus` and `2*sqrt(n_plus*n_minus) <= N` on each total
number sector. The second bound follows from `sqrt(n0*n1) <= N/2`. The
hopping operator is the second quantization of a two-by-two Hermitian matrix
with eigenvalues `+-1`, so its norm on the total-`N` sector is at most `N`.
Orthogonality of different input total-number sectors extends these bounds
to arbitrary vectors in the domain of `N`. The triangle inequality gives

```
||Pi H (I-Pi) psi|| <= kappa * sqrt(E[N^2*1_(N>=d)]),
kappa = 9/32.
```

Thus the constants refer to the physical compressed action, including the
occupation-boundary diagonal terms, and to the actual impulse and pointer.

### Projected Duhamel comparison and terminal branches

Let `psi` follow the infinite history and `psi_d` the finite history, embedded
in the infinite number space. On every interval with generator `G`, the
projected difference `e = Pi psi - psi_d` obeys

```
i de/ds = (Pi G Pi) e + Pi G (I-Pi) psi.
```

Variation of constants bounds its growth by the integral of the last term.
The finite generator is bounded and self-adjoint; the Gaussian infinite
states are in the required generator domains. The same `Pi` applies across
all intervals, so the comparison proceeds through the whole history without
introducing intermediate projection, resetting or renormalization.

The initial finite state retains total occupations through `d-2` and
normalizes that vector. If its infinite retained probability is `c`, then

```
||psi_prepared - Pi Omega|| <= ||psi_prepared-Omega||
                            = sqrt(2-2*sqrt(c))
                            <= sqrt(2*(1-c)) = e0.
```

This is a vector-norm bound; the pure-state trace-distance preparation bound
is not silently substituted for it. Let `T0` and `T1` be the rational bounds
`min_z M(z,R0)/z^d` and `min_z M(z,R1)/z^d`. Integrating the uniform residuals
over the impulse, model evolution and unit-length pointer parameter gives

```
e_projected <= e0
    + |j|*sqrt(d*T0/2)
    + t*kappa*d*sqrt(T0)
    + sqrt(d*T1/8).
```

The pointer's two Z-basis arms are orthogonal components of one Stinespring
state. Each arm satisfies these same bounds, so the normalized coherent
superposition does too. At the endpoint, the outside component is orthogonal
to the projected difference, and its squared norm is at most `T1`. Hence

```
e_full <= sqrt(e_projected^2 + T1),
D(instrument_d, instrument_infinite) <= min(1,e_full).
```

Here `D` is trace distance of the actual output states for the declared
prepared history, retaining the classical outcome and each unnormalized
terminal field branch. It is not an all-input diamond-norm bound. Readout,
record-preserving feedback and discarding registers are common channels and
cannot increase this distance. In particular the detector-probability error
is bounded by the same number. The bound never divides by a Born probability
and does not assert that every normalized conditional branch is close when
its probability is arbitrarily small.

At model time 1 and unit current, the rational certificate gives upper bounds
approximately 0.799447, 0.0771708 and 0.000195217 for `d=4,8,16`, respectively.
The zero and negative currents obey the same uniform bound. The finite
Strang bound can be added by triangle inequality to compare the ideal split
instrument with the infinite two-oscillator instrument. Physical service and
numerical execution bounds require separate terms.

For any fixed rational `z` in the admissible interval, the bound tends to zero
like a polynomial in `d` times `z^(-d/2)`, together with the squeezed-vacuum
preparation tail. Taking a minimum preserves this convergence. Directed
root and preparation arithmetic increase their precision with `d`, so the
reported bound has no fixed rounding floor. This proves occupation
refinement for this specified apparatus. It supplies no spatial refinement
or identification with a field measured in nature.

## 6. Rational detector prediction and the fair-response alternative

The supplied infinite two-oscillator model has

```
mu = (sin(1)-(4/5)*sin(5/4))/2,
p_infinite = (1+exp(-9/40)*sin(mu))/2.
```

For arguments in `[0,5/4]`, the alternating sine series has decreasing term
magnitudes. Consecutive partial sums therefore bracket the exact sine.
Directed rational rounding encloses `sin(1)` and `sin(5/4)` and hence `mu`.
The resulting rational mean interval lies within `[0,1]`, where sine is
increasing, so evaluating the two mean endpoints gives a valid enclosure
of `sin(mu)`.

The existing positive-exponential Taylor bound gives an upper bound `E`
on `exp(9/40)`. Its tail is at most `10^-30` and its final upward rounding
at most `10^-30`, so `E-2*10^-30` is a lower bound. Reciprocation gives
positive upper and lower bounds on `exp(-9/40)`. Combining those intervals
in the displayed probability formula yields the exact rational reference
interval retained by `reference_probability_interval()`. No floating
trigonometric value or measured detector rate enters this certificate.

Let `epsilon` be the sum of the occupation-instrument and finite Strang
bounds, capped at one. If `[L,U]` is the infinite reference interval, the
ideal finite source-word instrument has

```
max(0,L-epsilon) <= Pr(Y=+1|j=1) <= min(1,U+epsilon).
```

For the executed `d=16,n=128` design, the certified probability interval is
contained in `[0.514957066204, 0.517885937293]`, which excludes the
fair response `1/2`. The receipt retains its exact rational endpoints. Both
`j=0` fairness and negative-current reversal follow by the declared parity.
This separates the specified ideal coupled apparatus from a fair-detector
alternative. It is a construction prediction and supplies no identification
with observations of nature. The source executor's floating result is
checked against independent finite matrices; certified numerical execution
and physical-service errors require separate terms.

The nearest existing results are the [source proper-code
construction](../../../code/m1_source_realization/README.md), the
[source parity circuits](../../../code/m1_fermionic_source/README.md), and
the [native scalar feedback identities](../scalar_response/DERIVATION.md).
The added result is their explicit composition through prepared finite
scalar and pointer registers. This two-site action has no inherited
golden-grid spatial error, laboratory calibration, interactions or observed
physical data.
