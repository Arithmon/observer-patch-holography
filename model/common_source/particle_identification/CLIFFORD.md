# What the free eight-component continuum field represents

This is a classification of the **separate free three-axis continuation** in
[`mass_response/dispersion.py`](../mass_response/dispersion.py). It does not
identify the executed four-mode ring with a spatial Dirac field. It supplies
no observed electron, neutrino or quark identity.

## Continuum operator and spin

The source coefficients are `a=log(25/8)`, `b=log(8/5)`, `J=a-b` and
`Omega=a+b`. For the six Hermitian Clifford matrices already constructed by
the free extension, set

```text
alpha_j=Gamma_(2j),
beta=(Gamma_1+Gamma_3+Gamma_5)/sqrt(3),
M=sqrt(3)J, c=Omega/2,
H(p)=c sum_j alpha_j p_j+M beta.
```

`M` is a candidate rest-energy coefficient, not a physical mass in laboratory
units. The source has `J>0`. The four matrices `(alpha_1,alpha_2,alpha_3,beta)`
are Hermitian, square to one and anticommute. Thus `H^2=(c^2|p|^2+M^2)I8`.
The Lorentz gamma matrices `gamma^0=beta`, `gamma^j=beta alpha_j` obey
`{gamma^mu,gamma^nu}=2 diag(1,-1,-1,-1) I8`.

For cyclic `(i,j,k)`, the rotation-spin generators

```text
S_i=-(i/2) alpha_j alpha_k
```

satisfy `[S_i,S_j]=i epsilon_ijk S_k` and `sum_i S_i^2=3I8/4`.
They commute with `beta` and rotate the three `alpha` matrices as a vector.
This is a spin-one-half Lorentz-spinor representation. It does not derive
physical statistics from a matrix Casimir: a CAR quantum-field realization
and its physical observer dictionary are separate obligations.

The physical-convention chirality matrix is
`gamma5=-i alpha_1 alpha_2 alpha_3`. It commutes with the kinetic alpha
matrices and anticommutes with beta. Since `M>0`, the Hamiltonian mixes the
two chirality subspaces. The result is a massive Dirac-type representation,
not an isolated Weyl field.

## Exactly two equivalent Dirac copies

Define the two transverse mass directions and the Clifford-four grading

```text
eta1=(Gamma_1-Gamma_3)/sqrt(2),
eta2=(Gamma_1+Gamma_3-2Gamma_5)/sqrt(6),
K=alpha_1 alpha_2 alpha_3 beta.
```

`K` and each eta anticommute with every continuum generator, while K commutes
with the eta matrices. The three matrices

```text
T1=K eta1, T2=K eta2, T3=-i eta1 eta2
```

are Hermitian Pauli generators commuting with all four continuum matrices.
Their span with identity is the **complete complex commutant**, of dimension
four and isomorphic to `M2(C)`.

One direct completeness proof avoids a numerical rank assumption. The sixteen
ordered Clifford-four monomials are trace-orthogonal. Multiplying each by
`I,T1,T2,T3` gives sixty-four trace-orthogonal, nonzero matrices, spanning the
entire `M8(C)`. Every nonidentity Clifford-four monomial anticommutes with at
least one generator. Expanding an arbitrary commuting matrix in this full
basis therefore leaves only the four identity-monomial coefficients.
The module checks this complete basis and independently computes the linear
commutant nullity; floating thresholds are diagnostics of the exact argument.

The projectors `(I+T3)/2`, `(I-T3)/2` have rank four, and T1 intertwines
their actions. Each is one equivalent four-component Dirac representation.
The eight-component continuum candidate therefore contains **two degenerate
Dirac copies**, not one minimal Dirac species and not three generations.
A flavor basis change is arbitrary within this multiplicity; it is not an
observed flavor assignment.

For every momentum, `P_+=(I+H/E)/2` and `P_-=(I-H/E)/2` each have rank four.
Intersecting either with either flavor projector gives rank two: two spin
states for each energy sign and each Dirac copy. All these components are
retained.

## Charge and the particle/antiparticle interpretation

Complex-field phase rotation `psi -> exp(i q theta) psi` is a global U(1)
symmetry. Its normalization q is arbitrary here. The matrix Hamiltonian
does not identify it with electromagnetism or establish a local gauge field.

With the present gamma representation, `C=eta1` followed by complex
conjugation satisfies `C^2=I` and maps `H(p)` to `-H(-p)`. This is a
kinematic charge-conjugation structure of the free field. If a CAR field and
filled negative-energy reference are supplied, the usual normal-ordered
interpretation gives positive-energy particle excitations of charge `+q`
and negative-band holes of charge `-q`, with two spin states in each copy.
An empty Fock state is not automatically this reference. Neither the vacuum
nor the CAR/interacting spatial-field attachment follows from the eight
single-particle matrices alone.

Since `M>0`, the free continuation has no massless single-particle branch and
does not supply a photon pole. Above the filled negative-band vacuum, a
nonzero connected neutral fermion-bilinear response creates a particle-hole
pair, costing at least `2M`. This is a lower threshold bound; matrix elements
can suppress its attainment, and the total density at zero momentum is
conserved and has no such dynamical response. Thermal intraband responses
have a different state domain. A two-particle continuum threshold supplies
neither an isolated scalar bound-state pole nor a vector gauge-boson pole.

## Finite lattice symmetry differs

The finite Bloch operator retains all six independent matrices:

```text
H_lat(k)=sum_j[(a-b cos k_j) Gamma_(2j-1)+b sin k_j Gamma_(2j)].
```

Its universal all-momentum algebra is the full Clifford-six algebra `M8(C)`;
its momentum-independent internal commutant is only the scalars. To see this, differences between
the two momenta `+pi e_j/2` and `-pi e_j/2` recover `Gamma_(2j)`; their sum
minus twice the zero-momentum operator recovers `Gamma_(2j-1)`. All sixty-four
Clifford-six monomials are trace-orthogonal and span `M8(C)`. Special very
small periodic momentum sets can have larger accidental commutants; this is
not a claim that every individual eight-by-eight matrix has scalar commutant.
It also does not classify the full finite-space spectral commutant or forbid
momentum-dependent, spatially nonlocal symmetry maps.

Thus the continuum U(2) multiplicity symmetry is not an exact internal
symmetry of the complete finite regulator. Likewise the continuum rotation
and Lorentz interpretation is not an assertion of exact continuous Lorentz
invariance on a finite lattice. The existing controlled refinement connects
the free lattice operator to the classified continuum operator. It does not
complete an interacting continuum field theory.

## Why this does not identify Standard Model species

Under the ordinary internal-symmetry contract, an unbroken color action would
act in the Dirac multiplicity and commute with the continuum propagator.
There is no nontrivial Lie-algebra action of su(3) in this two-dimensional
multiplicity: a nonzero map from simple su(3) to u(2) would be injective,
but their real dimensions are eight and four. This excludes this specific
color assignment; it is not a no-go theorem for OPH with larger carriers,
additional fields or another physical dictionary.

Electrons or neutrinos also require their observed gauge currents, charge
normalization, chiral interactions, mass spectrum, statistics and physical
clock identification. Two equal-mass Dirac copies establish none of those
assignments. Scalar or vector composite poles would need a separate bound-state
and propagation calculation. A neutral finite-ring interference frequency
is not a Higgs or gauge-boson identification.
