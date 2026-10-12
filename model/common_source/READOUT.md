# Infinite coordinate-sign readout

The finite source pointer measures `B_d = 1_(Q_d > 0)`. An occupation-limit
comparison must also account for its difference from the infinite oscillator
projector `P = 1_(Q > 0)`: functional calculus does not commute with compression.
Already at `d=2`, `B_2=(I+X)/2`, whereas
`V_2* P V_2=I/2+X/sqrt(2*pi)`. Their operator-norm difference is
`1/2-1/sqrt(2*pi)`, approximately `0.10105772`.

`halfline.py` evaluates the original infinite projector directly on finite
polynomial approximants. For normalized Hermite functions `h_n`, integration
of the oscillator differential equations gives, for unequal indices,

```text
integral_0^infinity h_m(x) h_n(x) dx
  = [h_m(0) h_n'(0)-h_m'(0) h_n(0)]/[2(n-m)].
```

The diagonal is `1/2`; other equal-parity entries vanish. For `m=2j,n=2k+1`
the entry is

```text
(-1)^(j+k) sqrt((2j)! (2k+1)!)/
[sqrt(2*pi) 2^(j+k) j! k! (2k+1-2j)].
```

The standard Hermite recurrence, derivative and polynomial formulas are
tabulated in [NIST DLMF §18.9](https://dlmf.nist.gov/18.9) and
[§18.5](https://dlmf.nist.gov/18.5). The displayed half-line formula follows
by the boundary integration above.

In the unnormalized Bargmann basis `|n>_B=(a†)^n|0>`, multiplying by
`sqrt(m! n!)` removes every factorial radical. Thus a finite vector with
Gaussian-rational coefficients has positive weight `S/2+T/sqrt(2*pi)`:
`S` is its factorial-weighted norm and `T` is an exactly computed rational.
The untouched second oscillator contributes its factorial norm; different
spectator and gauge labels do not interfere. A separately supplied rational
initial normalization divides both terms. Current branch mass is never used
as a normalization denominator.

The implementation encloses pi by Machin's identity
`pi=16*atan(1/5)-4*atan(1/239)` and the signed remainder of each alternating
Taylor series. Integer square roots then enclose `1/sqrt(2*pi)`. No floating
quadrature, target detector probability, or observed convergence supplies
these bounds. Both signs are retained, and their exact weights sum to `S`.
Unnormalized inputs of mass greater than one retain that mass.

The complete postread vector is `P_b u`, represented by multiplying the
polynomial Gaussian wavefunction by the appropriate half-line indicator.
Its density kernel is `1_b(x0) 1_b(y0) u(x) u(y)*`. This retains all
coherences within each sign subspace and all earlier classical records.
The resulting vector generally has infinite occupation support. The routine
does not project it back into the polynomial support. Consequently a proved
trace-distance bound for the complete pre-read history passes through this
same infinite projective measurement by contractivity.

For a separate comparison with the executed finite pointer, let `A=V_d*P V_d`
and `B=B_d`. For any finite unnormalized pure input `u` (including arbitrary
untouched reference factors), put `p_b=<u|A_b|u>`, `q_b=<u|B_b|u>` and
`z_b=<u|A_b B_b|u>`. The exact half trace norm between the two complete
postread outputs is

```text
(1/2) sum_b sqrt((p_b+q_b)^2-4*abs(z_b)^2).
```

Here `A_1=A`, `A_0=I-A`, and likewise for `B`. Sum this expression over
earlier committed records. This identity includes the infinite output tail;
checking only the positive probability does not certify that quantum output.
It follows from the trace norm of a difference of two rank-one operators.
