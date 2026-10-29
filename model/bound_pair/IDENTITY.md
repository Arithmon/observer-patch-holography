# Identity and currents of the same bound pair

This CP-05 component serves [#1141](https://github.com/FloatingPragma/observer-patch-holography/issues/1141)
and the common excitation/action joins in [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
and [#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140).
It uses the same supplied Weyl-contact CAR action as the
[bound-band calculation](README.md). The new attachment is a conserved
current and a spin/statistics map for that excitation, rather than a
representation label copied from a separate matter model.

## Exact CAR charge and composite statistics

On the full Fock space let `b_x^*=c_(x,0)^* c_(x,1)^*` and
`N=sum_(x,a) c_(x,a)^* c_(x,a)`. CAR gives

```text
[N,b_x^*]=2 b_x^*,       [(-1)^N,b_x^*]=0,
(b_x^*)^2=0,            [b_x,b_x^*]=1-n_(x,0)-n_(x,1),
[b_x^*,b_y^*]=[b_x,b_y^*]=0  for x != y.
```

Thus the on-site pair carries **two units of the declared particle number**,
has even fermion parity, and disjoint pair operators commute. The on-site
algebra is hard-core; it is not the canonical boson oscillator algebra.
`identity.py` verifies these identities and all underlying CAR relations
on every input of an independent 16-dimensional four-mode Fock construction.
The charge equation is not established by evaluating `N=2 I` after sector
restriction: it is checked before restricting particle number.

The bound eigenstate belongs to that exact two-fermion sector and therefore
has number two at every momentum. This does not identify its charge with
electric charge. CAR statistics, the local spin components, the lattice and
the action remain supplied. No Standard Model representation, generation,
physical Lorentz spin, creation/decay process or real-space defect invariant
is assigned by these calculations. Extending a commuting on-site pair algebra
to an asymptotic multiparticle statistics theorem would require the separate
many-pair/scattering construction.

The countercontrols distinguish parity from number: a pair drive
`b_x^*+b_x` preserves fermion parity but violates number conservation; a
single-fermion drive also violates parity. Neither drive belongs to the
number-preserving field evolution tested here.

## Number transport on the actual action

Let `A_xy=sum_(a,b) c_(x,a)^* h_(xa,yb) c_(y,b)` for ordered neighboring
sites `x<y`. Define the outward number current

```text
J_(x->y)=i(A_xy-A_xy^*),       J_(y->x)=-J_(x->y).
```

The complete action, including its contact term, satisfies

```text
dN_x/dt=i[H,N_x]=-sum_y J_(x->y),
[U sum_x n_(x,0)n_(x,1),N_y]=0.
```

These are operator identities on the complete finite two-fermion space,
valid for every state and untouched reference. Time means Hamiltonian model
time with the supplied convention `hbar=1`. On the retained open cube the
independent construction checks all 351 unordered site pairs and all 27
local equations on 1,431 dimensions using Gaussian-integer sparse matrices.
Currents have denominator 16. Summing over any region cancels internal
currents and gives its complete boundary flux. A half-cube check retains all
nine crossing number currents. Removing the contact interaction leaves the
number-current **operators** unchanged, but changes evolved states and hence
their current expectations. This does not contradict the interaction-dependent
[energy currents](../field_stress/README.md), which are different operators.

An external Peierls link phase has Hamiltonian bond
`exp(i A) A_xy+exp(-i A) A_xy^*`; its derivative at zero is this same current.
Site phases act by `exp(i sum_x alpha_x N_x)`, with the corresponding
endpoint phase difference on each hopping matrix element. Tests verify both
the response sign and the full action's covariance, including unchanged
contact energy. This is a coupling to an externally specified probe of the
number symmetry, not a derived dynamical electromagnetic field or a gauge
sector selected by native repair.

## The singlet label has a restricted meaning

All three on-site internal spin generators commute with `b_x^*`; the
two-component pair is their singlet. At the eta momentum
`K=(pi,pi,pi)` the full bound eigenstate is exactly that singlet. Away from
eta, Weyl hopping mixes its spin content. The implementation computes the
positive triplet weight of the same isolated eigenstate at finite-torus
momentum `(0,pi,pi)` while checking its physical fermion exchange relation.

The global standalone internal spin generators
`S_j=(1/2) sum_x c_x^* sigma_j c_x` do not commute with the actual action.
Their nonzero full-sector commutators are checked with exact integer
numerators. Consequently an on-site singlet label does not establish a
conserved physical spin-zero particle throughout this band. A continuous
spacetime rotation representation and its physical spin assignment remain
separate obligations; no identification by a shared representation dimension
is made.

The finite identities, the numerical fiber diagnostic and their physical
interpretation have separate scopes. The tests also compare number flux
against an independent time derivative and reject an action mutation before
it can inherit the conservation certificate.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/bound_pair/test_identity.py
```
