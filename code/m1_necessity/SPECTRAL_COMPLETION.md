# Completing the free-field replacement beyond selected probes

Smooth fixed-mode agreement alone is too weak for a field-theory replacement.
The explicit T graph admits a uniform-weight action with extra low-energy
modes. The final replacement uses the following specified positive action
on the very same read graph to close that problem.

## A real failure of the uniform-weight control

Work first on the periodic q^3 box, with the T parameters of
BALANCED_REPLACEMENT.md. Write A_T for its uniformly weighted, exact-second-
moment normalized operator, and alpha=6/(h^2 M2). For every Fourier index
k=(q/m)*ell, all coarse displacements m*u have phase one. There are m^3
such distinct periodic modes. Only the 24t dyadic axis displacements act
on them, so every one has

    0 <= lambda_T(k) <= 48t*alpha.

The coarse ball contains the box K/4<=u_x<=K/2, |u_y|,|u_z|<=K/2.
Consequently M2>=m^2 K^5/64, and

    lambda_T(k) <= (18432/L^2)*t*2^(-7t) -> 0.       (S1)

Thus m^3 modes approach the mass frequency mu for positive mass. At any
fixed inverse temperature beta>0 their normal-ordered thermal energies are
bounded below by a positive constant per mode for all sufficiently large t.
The total diverges. The ordinary positive-mass continuum field in a fixed
periodic box instead has finite thermal energy. This is an explicit failure
of full thermal equivalence despite the uniform action's smooth-mode limit.
The alias-mode multiplicity and eigenvalue upper bound are exact, not a
small-grid inference or a numerical fit.

## The completed action, with no change to reads or entropy

Every comparison stencil already contains all six unit axis steps. Let

    (A_nn f)(x) = h^(-2) sum_(i=1..3) [2f(x)-f(x+h e_i)-f(x-h e_i)].

Use either zero exterior values or periodic boundaries consistently, and
define the final action for each V=S,A,B0,T by

    Ahat_V = (A_V + A_nn)/2.                          (S2)

This fixed one-half combination is part of the explicit construction, not
a fitted parameter or an assumption that a coercive action exists. Every
nonzero V read has positive weight alpha/2; a nearest neighbor has the
additional weight 1/(2h^2). No read is deleted, padded or assigned a different
entropy. The original and nearest-neighbor second moments each equal twice
the identity, so their average has exactly the required normalization.
Positivity, symmetry and the fixed-zero boundary terms are retained. Its
spatial remainder is bounded by

    (a^2+h^2)*M4/8 <= a^2*M4/4,                     (S3)

so the previous field and detector error bounds continue to hold. The
crossing, degree, routing and count theorems are unchanged, since (S2) uses
exactly the existing read graph. The centered-leapfrog obstruction at a/c
is not removed by changing positive weights; exact evolution is the field
comparison here, not an assertion of a native one-tick solver.

## Uniform spectrum control on the periodic box

For every principal Fourier wavevector p, |h*p_i|<=pi. Concavity of sine
on [0,pi/2] gives

    (4/pi^2)*|p|^2 <= lambda_nn(p) <= |p|^2.

Positivity and the exact second-moment inequality 1-cos(x)<=x^2/2 give
0<=lambda_V(p)<=|p|^2. Therefore the final action satisfies the all-mode bound

    (2/pi^2)*|p|^2 <= lambdahat_V(p) <= |p|^2.       (S4)

In particular the nonzero alias modes in (S1) have
lambdahat_T>=8/(h^2*m^2), which tends to infinity. There are no extra
bounded-energy alias species. Modewise convergence follows from (S3), and
(S4) controls the entire remaining spectrum rather than selected probes.

For fixed positive mu,c,hbar and beta, put omega=sqrt(c^2*lambda+mu^2).
Both the thermal logarithmic partition function and normal-ordered energy
are sums of

    -log(1-exp(-beta*hbar*omega)),
    hbar*omega/(exp(beta*hbar*omega)-1).

With delta=1-exp(-beta*hbar*mu)>0, these are bounded above by, respectively,
delta^(-1)*exp(-beta*hbar*omega) and
delta^(-1)*hbar*omega*exp(-beta*hbar*omega). Equation (S4) bounds them by a
constant times (1+|p|)*exp(-gamma*|p|), gamma>0, uniformly in the lattice.
This is summable on the fixed reciprocal lattice: Euclidean norm dominates
the l1 norm divided by sqrt(3), and one-dimensional exponential sums converge.
Dominated convergence proves the full finite-volume thermal partition,
normal-ordered energy and entropy limits, with entropy beta*E+log Z.
Vacuum zero-point energy is not claimed finite or identified by subtraction
of a fitted value; normal ordering is explicitly part of the comparison.

## The same result on the clipped box

The completed action also controls the zero-exterior geometry of the actual
read graph. Its form dominates one half of the nearest-neighbor form. View
grid vectors as equal-cell, piecewise-constant L2 functions extended by zero.
A bound on their nearest-neighbor energy bounds the L2 norms of all three
forward difference quotients. Integer translations telescope; fractional
cell translations give the same estimate up to an O(h) error. Consequently
the L2 translation modulus tends uniformly to zero as the translation and h
tend to zero. Fixed bounded support and the L2 compactness criterion give
strongly convergent subsequences. Summation by parts identifies weak gradient
limits; the zero extension makes each limit an element of H0^1 of the cube,
and gives the lower bound by its continuum Dirichlet energy.

The coarse part has the same full, not merely half, continuum lower bound.
For any compactly supported smooth test f, positivity gives

    <u_n,A_V u_n>_h >= 2<u_n,A_V f_n>_h-<f_n,A_V f_n>_h.

Consistency sends A_V f_n to -Delta f in L2, and the right side tends to
2<u,-Delta f>-||grad f||^2. Taking the supremum over C_c^infinity, dense in
H0^1, yields ||grad u||^2. Applying this to the two positive parts of (S2)
proves the exact lower form limit. Smooth compact test functions give the
upper form limit; their H0^1 density supplies recovery for general finite-
energy references. Thus the limit is the full Dirichlet Laplacian.

For completeness, eigenvalue convergence follows directly by min-max. Smooth
compact approximations to the first j continuum eigenfunctions give the
upper bound on the j-th discrete eigenvalue. The corresponding first j
discrete eigenvectors have bounded energy; the compactness above supplies
strong limits preserving their orthonormality. The lower form bound on
every linear combination gives the reverse min-max inequality. Hence each
ordered eigenvalue converges with its multiplicity. Degenerate eigenspaces
are treated by their whole spectral projections, not arbitrarily chosen
individual eigenvectors. Passing the discrete eigenvector equation against
sampled compact smooth tests identifies each subsequential limit with the
corresponding Dirichlet eigenspace. This proves convergence of the full
finite-rank projections at every spectral gap, as used below.

There is also a uniform bound on the entire ordered spectrum. For q sites
per axis the nearest-neighbor zero-exterior eigenvalues are

    (4/h^2) sum_i sin^2(pi*k_i/(2(q+1))),  1<=k_i<=q.

They are at least |k|^2/L^2 because h=L/q and sin(x)>=2x/pi on this range.
At most r^3 positive integer triples have |k|<=r. Min-max and (S2) therefore
give, for every j<=q^3,

    lambdahat_(n,j) >= j^(2/3)/(2L^2).              (S5)

No upper bound on the j-th frequency is needed for thermal domination:
omega*exp(-beta*hbar*omega/2)<=2/(beta*hbar*e). With a uniform positive
frequency gap, the energy term is therefore bounded by a constant times
exp(-beta*hbar*omega/2), and (S5) gives exp(-gamma*j^(1/3)). That sequence
is summable, as is the logarithmic partition bound. Dominated convergence
proves the full Dirichlet thermal partition, energy and entropy limits.
These Dirichlet statements also allow mu=0: (S5) supplies the uniform
frequency gap sqrt(mu^2+c^2/(2L^2)), so the same estimates use that gap
in place of mu. The periodic statements require mu>0 because their constant
mode has no massless oscillator vacuum.
Finite differences with q sites place the exterior zero nodes a distance h
outside one face; that vanishing shift is accounted for by the zero-extension
compactness argument and does not change the limiting cube.

## Smooth quantum instruments and the causal response

Canonically quantize these defined positive-frequency actions on the clipped
box, allowing mu>=0, in their Gaussian vacuum states and smooth coherent
displacements of those vacua. Ordered
spectral convergence with compactness gives convergence of bounded continuous
spectral functions on sampled smooth smearings. For discontinuous spectral
projections, the cutoff must lie outside the limiting spectrum; eigenvalue
convergence alone does not control a jump at a limiting eigenvalue. The
harmonic functions and inverse-frequency covariance used here are continuous
on the uniformly gapped spectrum. Position covariance uses B^(-1/2),
where B=c^2*Ahat+mu^2, and is uniformly bounded by the inverse frequency
gap just proved. For momentum
covariance B^(1/2), smooth compact f has uniformly bounded ||Ahat f_n||_h
by (S3); its contribution above spatial eigenvalue Lambda is O(Lambda^(-3/2)).
This supplies the missing uniform tail bound for the unbounded covariance.
The same estimates apply to bounded-time harmonic evolution and coherently
displaced smooth initial data. Therefore joint Weyl characteristic functions
converge for every fixed finite collection of smooth compactly supported
position and momentum smearings on the clipped box, not only Fourier probes.

The field commutator has the state-independent coefficient

    i*hbar <g_n, sin((s-t)*sqrt(B_n))/sqrt(B_n) f_n>_h.

Its limit is the Dirichlet Klein-Gordon commutator. Finite propagation of that
wave equation makes the limit zero when |s-t| is smaller than the distance
between the supports divided by c. Replacing the observation times by their
nearest declared layer ticks changes this coefficient by at most
hbar*delta_n*||g_n||_h*||f_n||_h per rounded time, because the derivative is
cos((s-t)*sqrt(B_n)), of norm at most one. Since delta_n->0, the same causal
response limit is obtained at the declared event ticks. This joins the
limiting free quantum response cone to the independently derived read/count
cone; it does not claim exact finite support or a native execution of the
matrix exponential.

All constructions, coefficients and boundary conditions are specified.
The spectral compactness, min-max, domination and quantum arguments here
are analytic proofs, not additional Lean declarations. Finite tests check
the complete spectrum of a small zero-exterior operator, the periodic alias
failure and its correction, and reject deletion of the stabilizing weights.
Thermal Gibbs entropy and the one-nat planar crossing statistic are different
observables; preserving both does not identify them or derive a horizon
entropy bridge. No interacting, gravitational or empirically calibrated
theory is inferred from these completed free-field comparisons.
