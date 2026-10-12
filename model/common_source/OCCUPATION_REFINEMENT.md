# Occupation limit of the same interacting observer history

The finite common construction has an occupation limit with its actual final
field instrument: the complete unnormalized classical/quantum history
converges in trace norm to a specified interacting two-oscillator history.
There is also an explicit finite-cutoff error bound **before** the final sharp
field readout. A finite-cutoff error rate for that readout is not supplied.
These are different conclusions and the executable certificate keeps them
separate.

This is a fixed two-site occupation limit. The two sites, their spring,
coupling, model-time evolution, preparation and readout are declared inputs.
It establishes neither a spatial continuum nor agreement with nature. Actual
source execution in this package remains at onsite dimensions 2 and 4; the
proof defines the larger occupation family without claiming those source
words have been executed.

## The family and the common observable interface

Use the same basis as [scalar_coupling.py](scalar_coupling.py), with gauge bits
0 and 1 denoting right occupations. On
`C^4 tensor Fock(C^2)`, let

\[
 N=n_0+n_1,\quad Q_j=(a_j+a_j^\dagger)/\sqrt2,\qquad
 H_g=X_0+X_1+\frac{I-Z_0Z_1}{2},
\]
\[
 H=H_g+N+1+\frac{\kappa}{2}(Q_0-Q_1)^2
                  +g(a_0+a_0^\dagger)Z_0,
 \qquad \kappa=\frac9{32},\quad g=\frac14.
\]

Thus `H_g` is the existing gauge-invariant matter/link action and electric
energy, and the interaction is exactly `g sqrt(2) Q0 (nL0-nR0)`.
Here `N` is the total oscillator occupation operator, not OPH's cosmological
capacity parameter.
For even `d`, let `Pi_d` retain `0 <= n0,n1 < d` and define
`H_d = Pi_d H Pi_d`. In particular the squares in `H` are compressed after
multiplication. Squaring a truncated quadrature would give another family.

The initial scalar vector is the ground state `Omega` of the free coupled
scalar Hamiltonian. Its common normal mode has frequency 1 and its difference
mode has frequency 5/4. The finite preparation `Omega_d` is its normalized
projection to **total** occupation `N <= d-2`, exactly the preparation already
used by the source compiler. It belongs to the rectangular `Pi_d` code.
The gauge preparation is `|00>` with the declared optional native-current
gauge kick; that kick commutes with every function of `N`.

The history evolves for model time 1, writes the gauge-bit-0 outcome `a`,
applies the specified conditional gauge rotation, evolves for model time 1/2,
then writes `b` using the positive/negative spectral projectors of `Q0`.
Its finite version uses the spectral projectors of `Q0,d = Pi_d Q0 Pi_d`.
All branches remain unnormalized. The gauge-bit projectors and conditional
gauge operations are the same operators at every cutoff and commute with
`Pi_d` and `N`. The final field projector is treated separately below.

The compatibility being tested is convergence of these same branch operators
under the number-basis embedding, not agreement between separately fitted
outcome tables. No extra scheduling limit is taken: the causal event order
and the two model-time intervals remain fixed.

## Domains and the limiting unitary

Finite linear combinations of number states, with arbitrary gauge vectors,
form a dense invariant domain `D`. On states supported in `N <= m`,

\[
 \|H\psi\|\le4(m+2)\|\psi\|.
\]

For example `||H_g|| <= 3`; putting
`b=(a0-a1)/sqrt(2)`, the spring is
`kappa (b+b^dagger)^2/2`. Its norm on this subspace is at most
`2 kappa (m+2)`, the linear term at most `2 |g| (m+2)`, and `N+1` at most
`m+2`. Since `3 <= (3/2)(m+2)`, the coefficient is at most
`1+2 kappa+2|g|+3/2 = 57/16 < 4`.

Each application of `H` increases total occupation by at most two. Consequently

\[
 \|H^k\psi\|\le4^k\prod_{j=0}^{k-1}(m+2j+2)\|\psi\|.
\]

The analytic-vector series has a positive radius, so the analytic-vector
criterion makes `H` essentially self-adjoint on `D`. This uses the standard
criterion stated as Theorem 1.7 in
[Simon, Topics in Functional Analysis](https://math.caltech.edu/SimonPapers/R4.pdf).
Hereafter `H` denotes that self-adjoint closure.

Extend `H_d` by zero outside its finite code. For each `psi in D`,
`H_d psi = H psi` once `d` is large enough. Core convergence therefore gives
strong resolvent convergence, and hence strong convergence of the unitary
groups at each fixed model time. These are the core and functional-calculus
results in Lemma 6.36 and Corollary 6.33 of
[Teschl, Mathematical Methods in Quantum Mechanics](https://diana.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf).
Their hypotheses have thus been checked on this particular Hamiltonian.

## A fourth-moment bound through the recorded feedback

Let `w(N)=(N+2)^4` and `M(t)=<w(N)>`. For a recorded history, the expectation
means the sum over its unnormalized branches, or equivalently the expectation
on its coherent record dilation. The same estimates hold with an untouched
reference system.

The number-conserving terms in `H` commute with `w`. In the normal-mode basis
only `kappa(b^2+b^dagger^2)/2` and `g(a0+a0^dagger)Z0` contribute to the
commutator. Resolve a vector into total-number sectors and put
`x_n=(n+2)^2 ||psi_n||`. Cauchy--Schwarz gives
`sum_n x_n x_(n+k) <= sum_n x_n^2 = M`.

For squeezing, set `x=n+2 >= 2`. The normalized coefficient multiplying
`x_n x_(n+2)` is bounded by `8 kappa`, since

\[
 \frac{[(x+2)^4-x^4]\sqrt{(n+1)(n+2)}}{x^2(x+2)^2}
 \le\frac{(x+2)^4-x^4}{x(x+2)^2}\le8,
\]

where the last numerator difference is
`8x(x+2)^2-((x+2)^4-x^4)=8x^2-16 >= 0`.
For the linear term the analogous coefficient is at most `4|g|`, because

\[
 \frac{[(x+1)^4-x^4]\sqrt{x-1}}{x^2(x+1)^2}\le2.
\]

An exact polynomial certificate for the squared inequality, with `y=x-2`, is

\[
 4x^4(x+1)^4-(x-1)[(x+1)^4-x^4]^2
 =4y^8+64y^7+440y^6+1692y^5+3976y^4
   +5884y^3+5500y^2+3175y+959\ge0.
\]

Hence, in every finite compression and in the infinite history,

\[
 |\dot M|\le C M,\qquad C=8\kappa+4|g|=\frac{13}{4},
 \qquad M(t)\le M(0)e^{Ct}.
\]

Compression cannot increase the number-sector operator norms used above.
One can first prove the differential estimate in finite dimension, then pass
to the strong unitary limit using lower semicontinuity of the positive moment.
The gauge measurement has `sum_a P_a^* w P_a=w`; each conditional gauge
rotation also preserves `w`. Thus their **combined unnormalized** branches
preserve this moment at the checkpoint. No division by an outcome probability
is needed, and the growth time in the full history is `T=1+1/2`.

In the initial free vacuum, `N=2R`, with probability generating function

\[
 \mathbb E[z^R]=\sqrt{\frac{80}{81-z}},\qquad
 f_k=\mathbb E[(R)_k]=\frac{(1/2)_k}{80^k}.
\]

Here `(R)_k` is the falling factorial and `(1/2)_k` is rising. Expansion gives

\[
 M_0=16(1+15f_1+25f_2+10f_3+f_4)
    =\frac{143748821}{8192000}.
\]

Conditioning this distribution on `N <= d-2` decreases the expectation of
the increasing weight, so the same initial bound holds for `Omega_d`.

## Quantitative comparison before the field readout

Only number-lowering or occupation-exchange terms can map the complement of
the rectangular cutoff back into it. The pair-annihilation term is
`kappa(a0-a1)^2/4`, bounded by `kappa N/2`. The exchange term
`-kappa(a0^dagger a1+a1^dagger a0)/2` is bounded by `kappa N/2`.
The remaining term is `g a0 Z0`, bounded by `|g| sqrt(N)`. Consequently

\[
 \|\Pi_d H(1-\Pi_d)\psi\|
 \le\kappa\|N1_{N\ge d}\psi\|
       +|g|\|\sqrt N1_{N\ge d}\psi\|
 \le\sqrt{M(t)}\left(\frac{\kappa}{d+2}
                  +\frac{|g|}{(d+2)^{3/2}}\right).
\]

The use of `N >= d` only enlarges the complement. Likewise
`||(1-Pi_d) psi|| <= sqrt(M(t))/(d+2)^2`.

Write `A=sqrt(M0)`,

\[
 I_T=A\frac{e^{CT/2}-1}{C/2},\qquad
 L_d=\frac{\kappa}{d+2}+\frac{|g|}{(d+2)^{3/2}}.
\]

For `r=d/2`, the actual preparation's omitted probability obeys

\[
 p_{\rm tail}\le\frac{81}{80}\frac{\binom{2r}{r}}{4^r81^r},
 \qquad e_0=\|\Omega_d-\Omega\|\le\sqrt{2p_{\rm tail}}.
\]

The tail estimate follows from the ratio of successive pair probabilities
being at most `1/81`. Projected Duhamel comparison on each time interval,
joined through the commuting gauge-record isometry, gives

\[
 e_{\rm proj}\le e_0+I_T L_d.
\]

The terminal outside component is orthogonal to that projected difference.
The full record-dilation vectors immediately before the field readout satisfy

\[
 \|\Psi_d-\Psi\|
 \le\left[(e_0+I_T L_d)^2+
              \frac{M_0e^{CT}}{(d+2)^4}\right]^{1/2}.
\]

This also bounds their pure-state trace distance, and after record dephasing
it bounds the pre-readout classical/quantum trace distance. The unclipped
bound tends to zero as `O(1/d)`. [refinement.py](refinement.py) implements it
using exact rational arithmetic and directed square-root/exponential bounds.
Those arithmetic utilities are reused from the scalar package; its free
Gaussian evolution bound is **not** imported into this interacting history.

## The actual sharp readout also converges

`sign(Pi_d Q0 Pi_d)` is generally different from `Pi_d sign(Q0) Pi_d`.
For example at `d=2`, the first is `X` on site 0, while the latter is
`sqrt(2/pi) X`. Applying a common-channel contraction bound would therefore
be invalid.

Instead consider the one-site truncated Jacobi matrix `Q_d`, extended by zero
on the outside subspace. Finite Hermite vectors form a core for the
self-adjoint multiplication operator `Q`, and `Q_d f=Q f` for every such
vector once `d` is large enough. Thus `Q_d -> Q` in strong resolvent sense.
The spectral measure of multiplication by `x` has no atom at zero for any
`L2(R)` vector. Its positive/negative spectral projectors therefore converge
strongly. This is precisely the non-atomic-endpoint condition of Lemma 6.39
in [Teschl's text](https://diana.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf).
For even `d` the finite matrix has no zero eigenvalue; the artificial zero
eigenspace outside the finite code disappears strongly as `Pi_d -> I`.

Tensoring with site 1, the gauge and the retained record preserves this
convergence. If `B_(d,b)` and `B_b` are the finite and infinite projectors,

\[
 \|B_{d,b}\Psi_d-B_b\Psi\|
 \le\|\Psi_d-\Psi\|+\|(B_{d,b}-B_b)\Psi\|\longrightarrow0.
\]

There are finitely many retained histories, so their full unnormalized
classical/quantum instrument outputs converge in trace norm. Subsequent
native-only record-conditioned updates preserve this conclusion by channel
contractivity. This covers both feedback enabled and feedback disabled, all
three declared native currents, and their finite classical mixtures.
In particular, the 1056-outcome history of [native_join.py](native_join.py)
is obtained from those conditional quantum histories by adjoining the same
finite native transitions and keeping their labels. Convex mixing and that
common classical/quantum channel preserve trace-norm convergence; no extra
cutoff factor of 1056 is required.

The spectral convergence argument supplies no finite-`d` numerical rate here.
Accordingly the pre-readout error number must not be reported as the final
detector error or as an error on the complete recorded instrument. A certified
finite-section approximation of the sharp spectral projectors, or an explicitly
changed finite-resolution detector, would be needed for that stronger bound.

Also, a sharp sign projector generally turns a smooth wavefunction nonzero
at the origin into one with a jump. The post-readout kinetic energy need not
be finite. The present history ends its quantum field operations at this
readout; the moment estimate cannot automatically be reused through another
sharp field measurement and another evolution segment.
The separately appended finite modular-clock experiment is also outside this
occupation theorem: no cutoff limit of its modular logarithm or clock action
is asserted here.

## What this adds and how to check it

The existing finite split estimate controls source product evolution against
`exp(-it H_d)` at a fixed occupation cutoff. The result here controls the
interacting occupation family itself and proves convergence of the same
recorded gauge-feedback/field-readout history. The underlying operator
convergence theorems are standard; the Hamiltonian mapping, moment constants,
committed-record composition and cutoff residual above are candidate-specific.

For the ideal palindrome at any fixed cutoff the same generator bound is
`S_d <= 3+(1+2 kappa)(2d-1)+|g| sqrt(4d-2) = O(d)`. Its operator error over
a fixed model-time interval is at most `S_d^3 t^3/(3 n^2)`. Taking, for example,
`n(d)=d^2` split steps on each of the two intervals makes their summed error
`O(1/d)`. At each `d` those two implementations use the **same** finite final
readout, so its channel is contractive for this split-error comparison.
Combining this with the proved occupation limit gives a simultaneous ideal
time/occupation limit for the same history. It still supplies no final sharp
readout rate and does not assert source execution above the retained `d=2,4`.

Run from the repository root:

```bash
python -m model.common_source.refinement
python -m pytest -q model/common_source/test_refinement.py
```

The tests independently check the vacuum moment, coefficient inequalities,
the compressed Hamiltonian's moment commutator, checkpoint preservation,
directed bounds, and the distinction between the two sign operators.
Numerical checks corroborate those calculations; the all-cutoff and strong
convergence claims rest on the proof above. Finite splitting, floating source
execution, physical source services and spatial refinement remain separate.
