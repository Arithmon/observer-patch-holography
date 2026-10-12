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

## 4. Error contract and remaining cutoff obligation

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

There is a precise route to a full oscillator comparison. For the finite
embedding `iota`, a normalized finite state `psi_d`, and its continuous
projected evolution,

```
||exp(-i*t*H)*iota*psi_d - iota*exp(-i*t*H_d)*psi_d||
  <= integral_0^t ||(I-Pi_d)*H*iota*exp(-i*s*H_d)*psi_d|| ds.
```

This is Duhamel's identity and unitarity. Quadratic couplings make the
residual depend only on the highest number layers. A complete bound must
also handle the impulse, split evolution and pointer factors on the same
history. End-point occupation, agreement of two cutoffs and unbounded
numerical quadrature do not certify this integral. The present preparation
bound therefore does not close a full infinite-oscillator or spatial
continuum error budget.

The nearest existing results are the [source proper-code
construction](../../../code/m1_source_realization/README.md), the
[source parity circuits](../../../code/m1_fermionic_source/README.md), and
the [native scalar feedback identities](../scalar_response/DERIVATION.md).
The added result is their explicit composition through prepared finite
scalar and pointer registers. This two-site action has no inherited
golden-grid spatial error, laboratory calibration, interactions or observed
physical data.
