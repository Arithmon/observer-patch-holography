# Chiral continuum symbol and its counting

This construction attaches a three-dimensional Weyl kinetic symbol to the exact
all-left internal action in [action.py](action.py). The internal representation,
family feature and Yukawa tensors are checked there. Here the declared spatial
attachment is

```math
h(p)=-c\,\boldsymbol\sigma\cdot p\otimes I_{45},\qquad
c=\Omega/2,\qquad \Omega=\log 5.
```

There are 15 internal all-left components and three existing family-coordinate
axes, with spin outermost. The coordinates, canonical kinetic metric and the map
from the reference frequency to coordinate speed are architecture choices.
The evolution parameter is the programmed modular parameter. Neither physical
seconds nor a speed in metres per second has been selected.

The checked objects are finite matrices evaluated on a finite menu of continuum
momenta. They supply a useful chiral candidate and an independent target for a
CAR source execution. They do not constitute a local nonperturbative regulator
of a chiral gauge theory. In particular, constructing or diagonalizing the target
symbol alone does not execute its local source implementation.

## Mass sign, opposite momentum and dispersion

Let $K=K^T$ be the internal all-left mass kernel, and let
$\epsilon=i\sigma_2$. In the action convention used here,

```math
\mathcal L_m=+\frac12\chi^T(\epsilon\otimes K)\chi+\mathrm{h.c.}
```

The creation coefficient of the CAR Hamiltonian
$H_{\rm pair}=\frac12c^\dagger\Delta c^{\dagger T}+\mathrm{h.c.}$
is $\Delta=-(\epsilon\otimes K)^\dagger=\epsilon\otimes K^*$.
The conjugation matters for a complex symmetric kernel. The default neutral
Higgs and native trial Yukawa matrices are real, so the default equals the
pairing kernel exported by the action.

In the Nambu order $\Xi(p)=(\chi(p),\chi^\dagger(-p)^T)$,

```math
\mathcal H(p)=
\begin{pmatrix}
h(p)&\epsilon\otimes K^*\\
(\epsilon\otimes K^*)^\dagger&-h(-p)^T
\end{pmatrix}.
```

The lower block contains **opposite momentum**. Replacing it by $-h(p)^T$
gives a different Hamiltonian and fails the massive dispersion check. Using
$\sigma_i\epsilon+\epsilon\sigma_i^T=0$, direct multiplication gives

```math
\mathcal H(p)^2=
\begin{pmatrix}
c^2|p|^2I+I_2\otimes K^*K&0\\
0&c^2|p|^2I+I_2\otimes KK^*
\end{pmatrix}.
```

Since $K$ is symmetric, these internal factors are positive semidefinite
and have the squared singular values $s_i(K)^2$. Thus the symbol has

```math
E_i(p)^2=c^2|p|^2+s_i(K)^2.
```

These are genuine gaps of this quadratic propagation law, rather than repair
relaxation rates or chemical potentials. Calling them measured physical masses
would additionally require the physical field attachment, vacuum, interactions,
clock and length calibration, and controlled comparison to observations. The
trial Higgs norm $v=\Omega$ is supplied; spontaneous selection of that norm
has not been demonstrated.

Changing the hole spin basis by
$V=\operatorname{diag}(I,\epsilon\otimes I_{45})$ displays

```math
V\mathcal H(p)V^\dagger=
\begin{pmatrix}
-c\sigma\cdot p\otimes I&I_2\otimes K^*\\
I_2\otimes K&+c\sigma\cdot p\otimes I
\end{pmatrix}.
```

The two blocks therefore have opposite spatial chirality and mix through the
Higgs mass. The right-handed particle component is the conjugate of an
all-left $u^c,d^c,e^c$ component. It is not a newly introduced mirror field.

## Spin, chiral currents and charge

The raw Nambu rotation generators are
$S_i=\operatorname{diag}(\sigma_i/2,-\sigma_i^T/2)\otimes I_{45}$.
They satisfy $[S_i,S_j]=i\epsilon_{ijk}S_k$ and
$\sum_iS_i^2=3I/4$. In the spin-conjugated basis both blocks transform by
$\sigma_i/2$. For $p\ne0$, $S\cdot p/|p|$ commutes with the free
massive symbol. A massless all-left positive-energy particle has helicity
$-1/2$; its positive-energy antiparticle has helicity $+1/2$.
Helicity at zero momentum is undefined. A nonzero mass mixes left and right,
so definite chirality is not a massive energy eigenstate.

Each of the 12 anti-Hermitian internal currents $T$ acts on the kinetic
operator as $I_2\otimes T$ and commutes with it. On the raw Nambu symbol it
acts as $\operatorname{diag}(I_2\otimes T,I_2\otimes T^*)$. A fixed Higgs
background preserves a current exactly when

```math
T^TK+KT=0.
```

The neutral background preserves the nine-dimensional color plus electric
charge algebra. The remaining three current directions also transform the
Higgs; their failure to commute with a fixed neutral kernel is the expected
broken symmetry, not a failure of the full current/Higgs identities checked
in the action. The floating rank in `report()` is a diagnostic of the supplied
current span, not a separate exact rank certificate.

Electric charge is
$Q_N=\operatorname{diag}(I_2\otimes Q,-I_2\otimes Q^T)$, where
$Q=T_3+Y$ in the all-left convention. It commutes with $\mathcal H$.
Charge is classified with its projectors: arbitrary eigenvectors in degenerate
energy subspaces are not reliable charge labels. This internal normalization
does not determine the electromagnetic coupling. Gauge currents alone do not
provide dynamical gauge fields or establish massless gauge-boson poles.

### Visible spectral weights without discarding modes

The kinetic report uses the
[finite spectral diagnostic](../../observer_particles/spectral_response.py)
on the same 180-dimensional $\mathcal H$ and $Q_N$. It fixes four canonical
coordinate interpolators before diagonalization: the first internal
coordinate in each spin and Nambu block, at indices $0,45,90,135$. These
are a small coordinate probe, not a selected set of energy eigenvectors or
an implemented charged detector apparatus.

For exact joint energy/charge projectors the matrix residue is
$B^\dagger P_{E,q}B$. Summing over the entire retained domain gives
$B^\dagger B=I_4$ for this probe. Every signed frequency, zero mode and
invisible subspace remains in the output. No projected weight is
renormalized. A mode dark to these four columns can be visible to another
interpolator; darkness is not a decoupling argument. The original
positive-energy charge census still counts the complete symbol, including
all modes invisible to this probe.

Numerically unresolved frequencies are reported together with their full
frequency lists and intervals. A cluster dimension is not an assertion of
exact degeneracy or a physical particle multiplicity. Residue ranks and
darkness are numerical diagnostics at the reported threshold, with input
Hermitian-part, commutator and block-projection errors exposed. No exact
rank or outward spectral-error certificate is claimed.

This calculation extends the supplied-action check by its interpolator
weights. It does not obtain $\mathcal H$ or $B$ from observer records.
Charged preparation, readout and any required phase-reference service,
as well as the physical state and clock, must be accounted for by a
record reconstruction. The
[field-reconstruction guide](../../../docs/research/OPH_FIELD_RECONSTRUCTION.md)
keeps that inverse separate from finite action spectroscopy.

### A direct weak-chirality discriminator

`weak_chirality_report()` reconstructs the normalized fundamental weak raising
matrix $E_{34}$ as an exact linear combination of the **actual twelve source
currents**. It applies the same coefficients to their exterior lifts. The
charged-lepton doublet, conjugate singlet and neutrino components are decoded
from the actual electric-charge and weak-Casimir projectors. The label of a
component alone is not used as evidence for its current action.

Let $J_+$ denote that internal matrix. Its Nambu action is
$\operatorname{diag}(I_2\otimes J_+,-I_2\otimes J_+^T)$. On the
physical charge $-1$ sector, $J_+^\dagger J_+$ is the projector onto
the charged left-doublet block. It vanishes on the charge-conjugated-singlet
hole block. At zero Higgs norm and nonzero momentum, positive-energy
negative-helicity charge $-1$ states therefore have normalized total
weak-current weight one; the positive-helicity states have weight zero.

For a native family mass parameter $m_i$, the corresponding massive
positive-energy helicity branches have total weights

```math
w_{\lambda,i}=\frac12\left(1-2\lambda\,
\frac{c|p|}{\sqrt{c^2|p|^2+m_i^2}}\right),\qquad \lambda=\pm\frac12.
```

This follows by multiplying the positive-energy and helicity projectors of
the four-dimensional charged particle block, then tracing against
$J_+^\dagger J_+$. The positive-helicity weight approaches
$m_i^2/(4c^2|p|^2)$ in the relativistic limit. The report computes these
weights from the matrices and compares them with the formula. At nonzero mass
both helicities can couple, as required by left/right mixing.

The complete final-state space is summed in $J_+^\dagger J_+$ **before**
projecting onto an incoming branch. Restricting the final state to one
same-momentum neutrino spin can give a kinematic zero that does not measure
total left-current strength. The reported quantity is a normalized vertex
weight, not an on-shell decay rate or scattering cross section. It neither
constructs a propagating weak gauge field nor identifies an observed electron.

## Independent modes, antiparticles and the reference state

The internal kernel has rank 42. Its only zero rows are the three left-neutrino
family components. Its charged singular values appear twice because each
Dirac mass joins two all-left fields of opposite electric charge. These repeated
singular values are not additional families.

Counting color and family, the declared content is nine up-type, nine down-type
and three charged-lepton Dirac fields, together with three left-Weyl neutrinos.
These are representation types of the candidate, not identifications of observed
species. In particular the native family axes do not identify electron, muon and
tau, and the current trial hierarchy fails their observed hierarchy. Massless
neutrinos also fail to account for observed neutrino mass splittings.

At nonzero momentum the 180-dimensional Nambu symbol has 90 positive and 90
negative eigenvalues. The positives contain 84 massive modes: 21 Dirac fields
times two particle spins and two antiparticle spins. The other six are the
three left-Weyl particle helicities and their three antiparticle helicities.
At zero momentum there are 84 positive, 84 negative and 12 Nambu-zero entries.
Those 12 entries encode six complex CAR zero modes, not twelve neutrinos or an
independent sterile right-handed neutrino.

Nambu particle-hole conjugation gives
$\tau_x\mathcal H(p)^*\tau_x=-\mathcal H(-p)$. Negative eigenvalues are
redundant conjugate partners at opposite momentum. They are not an additional
species list. A particle/antiparticle interpretation uses the corresponding
quasifree ground-state reference, with a choice of zero-mode occupations where
needed. It does not use the empty encoded Fock vacuum as the physical vacuum.

For distinct $p$ and $-p$, an actual source has CAR annihilators
$(\chi(p),\chi(-p))$, with one-body matrix
$\operatorname{diag}(h(p),h(-p))$ and antisymmetric pairing

```math
\Delta_{p,-p}=\begin{pmatrix}0&\Delta\\-\Delta^T&0\end{pmatrix}.
```

Its doubled BdG matrix splits into the displayed $p$ and $-p$ symbols
after a permutation. At $p=0$ the original 90 CAR annihilators and their
$\Delta$ suffice; duplicating the zero-momentum source would double the
CAR degrees of freedom by mistake. `momentum_pair_source` rejects this case.

## Explicit regulator failure control

The naive local choice
$h_{\rm naive}(k)=-c\sum_i\sigma_i\sin k_i$ has a zero at every
$k_i\in\{0,\pi\}$: eight nodes. The chirality of its linearization is

```math
\operatorname{sign}\det\!\left[-c\,\operatorname{diag}(\cos k_i)\right]
=-(-1)^{\#\pi}.
```

There are four nodes of each chirality, with zero net chirality. An even torus
contains the corners exactly. An odd torus merely omits their exact momenta;
near-mirror excitations approach zero under refinement. This is checked as a
negative control, and is not advertised as a mirror-free implementation.

A finite spectral menu of the continuum symbol avoids adding those particular
nodes by construction. Its position representation can be nonlocal, and local
gauge transformations can mix retained and discarded momenta. Exact global
current identities and anomaly traces do not establish a local gauge-covariant
cutoff, a positive interacting theory or a controlled continuum limit. Those
regulator and interaction obligations remain before this candidate could serve
as a complete Standard Model simulation.

## The same clock requires a detector-history map

The frozen-chart Clifford and mass-shell identities in
[spacetime_attachment.py](spacetime_attachment.py) align principal symbols.
They do not equate the native instrument's measured response with the
unmeasured quadratic Weyl–Higgs propagator. The actual incident-clock
response gives a quantitative test of that proposed identification.

Use the native collective insertions $X_j=|0,j\rangle\langle5,j|$ on the
binary total-two source. The current/contrast law and mixed-cap selector
remain the ones executed by
[collective_response.py](../hierarchy_diagnosis/collective_response.py).
Let $I$ and $D$ be its unnormalized incident and remote maps on the ten
collective coherences. These are response maps, not work-state probability
matrices. Their actual seam probabilities and coherence losses are

```math
a=\frac5{32},\qquad r=\frac{11}{16},\qquad
k_C=\frac3{64},\qquad k_E=\frac3{32}.
```

Every killing seam is incident on the clock cap. Exact native branch
enumeration gives nonnegative entries and

```math
I\mathbf1=(a-k)\mathbf1,\quad \mathbf1^T I=(a-k)\mathbf1^T,
\qquad D\mathbf1=r\mathbf1,\quad \mathbf1^T D=r\mathbf1^T.
```

For proposal deadline $h$, success means the first cap-incident work event,
including zero-current incidents. Its probability is
$\beta_h=a\sum_{n=0}^{h-1}r^n$, independent of the work state. The
**conditional successful response** is

```math
R_h=\frac{I\sum_{n=0}^{h-1}D^n}{\beta_h},\qquad
\lambda=1-\frac{k}{a}
=\begin{cases}7/10,&C,\\2/5,&E.\end{cases}
```

Both row and column sums of $R_h$ equal $\lambda$. Positivity gives
$\|R_h\|_2\le\sqrt{\|R_h\|_1\|R_h\|_\infty}=\lambda$;
the uniform vector attains equality. Thus its operator norm and spectral
radius are exactly $\lambda$, for every finite positive deadline. The
nine outer coordinates form an invariant block with the same uniform
eigenvector and bound. No floating rank or spectral assignment is used.

The complete Ramsey experiment measures this loss. At $h=2$, success and
selector-timeout probabilities are each $135/512$, and deadline timeout
has probability $121/256$. Summing the actual signed Ramsey outputs on
successful histories gives $\lambda(135/512)$, hence conditional contrast
$7/10$ or $2/5$ from input contrast one. All 21 detector outcomes and all
terminal flags are retained; conditioning is an explicitly stated
diagnostic. Preparation of the coherent probe and the Ramsey apparatus
remain declared services. This audit executes the native response but
does not compile another apparatus pulse sequence.

For any unitary $U$ and any linear response map $J$, the reverse triangle
inequality gives

```math
\|UJ-JR_h\|_{\rm HS}\ge(1-\lambda)\|J\|_{\rm HS}.
```

Consequently a nonzero exact intertwiner cannot map this conditioned
native response to the **unmeasured closed quadratic** Nambu evolution
$U=e^{-i\mathcal H\Delta t}$. The relative defect is at least $3/10$ for C
and $3/5$ for E. Changing the tetrad basis or regraduating a constant clock
step preserves unitarity and does not remove this bound. Tests also apply
it to the actual 180-dimensional candidate symbol. This rejects a direct
closed quadratic response identification; it does not reject an effective
field coupled to detectors, an environment, interactions or a thermal state.

The environment distinction has an explicit control. The complete channel
with $K_0=\sqrt{1-p}\,1$, $K_1=\sqrt p\,Z$ and
$p=(1-\lambda)/2$ preserves charge populations and multiplies its charged
matrix unit by $\lambda$. Its two outcomes are complete. It supplies a
mathematical dilation of the uniform response, with $p=3/20$ for C or
$3/10$ for E; it is not an identification of the entire native spatial
channel or a physical detector. A viable field-history bridge must specify
and execute the corresponding environment and readout, not compare only
the unmeasured Hamiltonian's eigenvalues.

Even the finite tick count is not a complete history. With exactly one
incident tick, the conditional source-site $2$ to detector-site $7$
response is zero at deadline one and $1/540$ at deadline two. Remote
repairs and their waiting records affect the propagation before the tick.
A continuum clock concentration theorem may control this distinction in
its stated limit; replacing these finite records by the count alone is
not an exact equality. None of the contraction factors above is identified
as a physical rest mass, and no Standard Model species assignment follows
from this test.

Run the independent checks from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest \
  model/common_source/chiral_candidate/test_kinetic.py -q
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest \
  model/observer_particles/test_spectral_response.py -q
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest \
  model/common_source/chiral_candidate/test_spacetime_attachment.py -q
```
