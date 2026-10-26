# OPH particle mass response and hierarchy audit

This document states the current chiral Weyl–Higgs candidate, hierarchy
diagnosis and full repair-response result. Its calculations and retained
outputs are in the linked model packages.
The native recorded repair supports a propagating **collective core-relative
charge response**, an exact correlation current, and compatible local detector
records. Exact signed-current histories also recover nonzero finite-frequency
occupation readouts. The collective operators are even matrix-unit transfers
with a diffusive response. A specified 15-dimensional exterior carrier fits
inside the same native 66-state space, but the complete recorded repair does
not preserve any noiseless 15-dimensional work subspace.
Active record correction on multiple native carriers supplies a different
finite construction: the same declared Weyl–Higgs field propagates at zero
and three opposite momentum pairs, with complete continuing detector and
clock outputs. Its representation, trial Yukawa law, spectral cutoff and
coherent services remain supplied.
An explicitly [calibrated action control](#56-calibrated-action-and-matching-control)
also realizes the charged-family matrices and measured mixing inputs on this
representation. Its external Standard Model mass matching is a calibration
check, not a mass measurement from native observer records or a derivation.
The measured Standard Model mass hierarchy has **not** been recovered, and
these collective operators have not been identified with Standard Model fields
in the observer's physical Lorentzian description.

**Operational particle mass.** A mass is inferred from jointly accessible
detector records and interventions in a stated observer clock and spacetime
description. Observers must agree on the same response law and on the
energy–momentum quantity extracted from it. A sharp stable rest mass requires
the corresponding physical pole and dispersion; thresholds and unstable or
confined species require their stated alternatives. Agreement on a detector
outcome after its state update is an essential observation contract. It does
not by itself identify a repair decay rate, modular frequency or graph
eigenvalue with that mass observable.

- [Candidate construction](#candidate-construction)
- [Hierarchy diagnosis](#hierarchy-diagnosis)
- [Full repair response](#full-repair-response)
- [Current conclusions and open tests](#8-current-conclusions)

## Candidate construction

This section describes a concrete OPH-compatible chiral matter architecture and
what its executable tests establish. It combines the existing internal
matter representation, verified Higgs interaction tensors and a native
family-dependent operator. Its fixed-Higgs-background fermion sector is
executed as a 90-mode CAR system, with a continuing source-clock experiment
that reads three charged rest gaps.

The candidate has the minimal Standard Model's chiral fermion representation
types and the expected Weyl–Higgs dispersion. **Its chosen mass prescription
does not reproduce measured particle hierarchies.** Three fixed nonlinear
comparison prescriptions also fail the joint charged-lepton hierarchy test.
No observed particle spectrum, complete interacting Standard Model, or
theory of gravity is claimed as derived by this construction.

The implementation and retained evidence are in
[the candidate package](../../model/common_source/chiral_candidate/README.md),
[its receipt](../../model/common_source/chiral_candidate/results.json) and
[its quantum-state artifact](../../model/common_source/chiral_candidate/results.npz).
The equations below explain their construction independently of the code.

### 1. The question being tested

The [continuing response experiment](../../model/common_source/mass_response/README.md)
joined a finite charged gap, interactions, neighboring-site response and
quantum-clock measurements in one history. Its
[operator-response audit](../../model/common_source/particle_identification/README.md)
showed that this gap has nonzero CAR spectral residue and can be read by an
even, charge-conserving transfer instrument. It also classified the separate
eight-component free continuum extension as two degenerate Dirac copies.

That establishes a charged fermionic response in the tested model. The
chiral candidate combines the existing OPH matter and Higgs representations
with an explicit propagating field, whose gaps are measured against the same
clock convention. The representation attachment and its family mass law are
declared construction choices.

The candidate supplies a finite execution of that representation and a direct
mass test that rejects its simplest trial coupling map.

### 2. Source inputs and construction choices

The calculation has three independent inputs: an internal representation,
a prepared native reference and a native spatial response band. Their
combination into a coherent action is a declared architecture choice.

```mermaid
flowchart TD
    A[OPH internal matter and Higgs tensors] --> D[Declared chiral Weyl-Higgs action]
    B[Native reference: a, b, J, Omega] --> D
    C[Native rank-three family feature F] --> D
    D --> E[CAR source evolution and quantum clocks]
    E --> F[Measured trial gaps and retained records]
    F --> G[Separate comparison with observed hierarchies]
```

#### 2.1 Internal matter and Higgs channels

The existing
[representation decoder](../../model/observer_particles/representation_identification.py)
recovers the nonvacuum even exterior module

```math
\mathcal M_{15}=\Lambda^2(\mathbb C^5)\oplus\Lambda^4(\mathbb C^5)
               =\mathbf{10}\oplus\bar{\mathbf5},
```

under its stated source-linked lift. Its five joint charge/Casimir blocks
are the following, in **all-left Weyl notation**:

| Internal field type | $(SU(3),SU(2),Y)$ | Internal dimension |
| --- | --- | ---: |
| $Q$ | $(3,2,1/6)$ | 6 |
| $u^c$ | $(\bar3,1,-2/3)$ | 3 |
| $d^c$ | $(\bar3,1,1/3)$ | 3 |
| $L$ | $(1,2,-1/2)$ | 2 |
| $e^c$ | $(1,1,1)$ | 1 |

Here $Q_{\mathrm{electric}}=T_3+Y$. The superscript $c$ denotes the
left-handed charge conjugate of a right-handed singlet. Those signs matter
when forming both masses and currents. The minimal module contains a left
neutrino and no independent right-neutrino singlet.

The existing
[one-Higgs construction](../../model/observer_particles/unified_yukawa_action.py)
supplies a single $H=(1,2,1/2)$, its conjugate, and exact invariant tensors
for the channels

```math
Q u^c H,\qquad Q d^c H^*,\qquad L e^c H^*.
```

The new action uses these tensors directly and rechecks their current
covariance and anomaly traces. The internal representation is conditional
on the stated lift; its existence does not force a unique field action or
choose Yukawa coefficients. The
[Standard Model comparison map](../../model/common_source/particle_identification/SM_MATCH.md)
records the precise existing theorem and ledger dependencies. For physical
conventions, see the [PDG electroweak review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf);
its hypercharge variable is twice the $Y$ used here.

#### 2.2 Exact native reference coefficients

The [public reference producer](../../model/native_modular_reference/local_preparation.py)
defines an ideal faithful distribution on 66 binary configurations with
total native load two. A distinguished two-port core has occupation
$q=0,1,2$, multiplicities $45,20,1$, and respective weights **per
configuration**

```math
(w_0,w_1,w_2)=\left(\frac{2}{225},\frac1{36},\frac2{45}\right).
```

The logarithmic source statistics used in the candidate are

```math
\begin{aligned}
a&=\log(w_1/w_0)=\log(25/8),\\
b&=\log(w_2/w_1)=\log(8/5),\\
J&=\log\frac{w_1^2}{w_0w_2}=\log(125/64)=a-b,\\
\Omega&=\log(w_2/w_0)=\log5=a+b.
\end{aligned}
```

For the two core occupations, the modular generator $K=-\log\sigma$ can
be written as a constant plus $-a(n_0+n_1)+J n_0n_1$. These coefficients do not
themselves create a propagating Weyl Hamiltonian. Assigning them to the
coherent channels below is part of the architecture being tested.

This calculation consumes the **exact ideal reference target**. It does
not execute a finite native warmup or propagate its preparation error.
The protected core, coherent controls, blank states and clock calibration
are declared operational inputs.

#### 2.3 A native family symmetry breaker

The prepared reference is not invariant under the full carrier rotation
group. Its one-body occupation means are

```math
\mathbb E[n_i]=
\begin{cases}
29/90,&i\text{ in the two-port core},\\
61/450,&i\text{ outside it}.
\end{cases}
```

Let $B$ span the actual selected rank-three native slow band and
$G=B^T B$. Compressing the occupation matrix gives

```math
F=G^{-1}B^T\mathrm{diag}(\mathbb E[n_i])B.
```

The [exact producer](../../model/common_source/particle_identification/family_attachment.py)
obtains three distinct positive eigenvalues:

```math
\lambda_1=\frac{61}{450},\qquad
\lambda_{2,3}=\frac{41}{225}\mp\frac{7\sqrt5}{750}.
```

They are approximately $0.135556,0.161352,0.203092$. The calculation
checks all 60 rotation-covariance identities, explicitly transports the
different native and retained-family port charts, and obtains the exact
group average $I_3/6$.

This supplies a concrete nonscalar native family operator. The earlier
unbroken-family theorem, whose invariant commutant is scalar, remains
correct: its full-invariance hypothesis is absent for this prepared
reference. Physical family interpretation and a coupling to masses remain
additional attachments. The protected-core preparation has not thereby
selected the electron, muon and tau basis.

### 3. The declared Weyl–Higgs action

The candidate uses two left spin components, 15 internal components and
three native family axes: **90 complex CAR modes** in one momentum block.
The trial family matrices and neutral background are

```math
Y_u=aF,\qquad Y_d=bF,\qquad Y_e=JF,\qquad
H=(0,v),\quad v=\Omega.
```

The kinetic coordinate speed is chosen as $c=\Omega/2$. The value $v$
is the neutral component norm in this convention, not an independently
derived electroweak scale in GeV.

Using the verified Higgs tensors produces a symmetric 45-dimensional
all-left internal mass kernel $K$. With
$\epsilon=i\sigma_2$, the real-background creation pairing is

```math
\Delta=\epsilon\otimes K,\qquad \Delta^T=-\Delta.
```

The transpose antisymmetry is required by CAR. In particular, a Dirac mass
is formed through the appropriate Higgs line and charge-conjugated singlet;
it is not inserted directly between unrelated gauge representations.
The neutral background preserves the assigned electric charge.

For left fields $\chi$, the kinetic symbol is

```math
h(p)=-c\,\sigma\cdot p\otimes I_{45}.
```

The opposite-momentum Nambu operator is

```math
\mathcal B(p)=
\begin{pmatrix}
h(p)&\Delta\\
\Delta^\dagger&-h(-p)^T
\end{pmatrix}.
```

For a general complex action kernel the creation coefficient uses the
corresponding conjugation; the implemented default is real. Spin conjugation
of the hole block gives the standard massive chiral form, with

```math
E_j(p)^2=c^2|p|^2+s_j(K)^2.
```

Negative Nambu energies are conjugate partners, not extra mirror families.
Counting color and family, the massive channels form **21 charged Dirac
fields**, and the minimal module leaves **three massless left-Weyl neutrinos**.
The latter does not reproduce observed neutrino masses. At zero momentum
the 12 Nambu zero eigenvalues represent six complex CAR zero modes, not
twelve neutrino species.

The [action proof](../../model/common_source/chiral_candidate/ACTION.md) and
[kinetic proof](../../model/common_source/chiral_candidate/KINETIC.md) give
the matrix ordering, signs, count and dispersion checks.

### 4. What establishes chirality in this candidate

All-left notation alone is insufficient to establish a weak-current
discriminator. The implementation reconstructs a weak raising operator
from the actual source current generators and checks its action on the
decoded lepton projectors.

In the massless limit, its total charged-current coupling strength on the
physical charge-minus-one energy branches is one for negative helicity and
zero for positive helicity. With a nonzero trial mass, the weights follow

```math
w_{\mp}=\frac12\left(1\pm\frac{c|p|}{E}\right).
```

The check uses the complete $J^\dagger J$ final-state sum; it does not
mistake a same-momentum final-neutrino spin selection for chirality. These
are matrix vertex weights of the declared action, not measured scattering
cross sections. No dynamical $W$ field is executed. A charged current
drive would need its compensating charged apparatus or gauge field.

### 5. Source execution and the observer clock

#### 5.1 Finite CAR computation without a huge Fock vector

The source backend executes native two-mode hopping, pairing and onsite
words on **all four local Fock input columns**. Their completed action on
Majorana operators determines the full finite Gaussian transformation.
Nonadjacent modes use fermionic SWAP routing, including the minus sign on
double occupation. Every routing pulse is counted.

At rest the chosen background separates into 42 charged pair blocks and
six neutral zero modes. All are retained. An initial fermionic wire layout
places the three probed lepton blocks first. Their six-mode Fock factor and
two clocks are also executed literally, with both clock purification
references retained. This independent calculation checks the Gaussian
reduction rather than replacing a missing large-state simulation with an
unsupported assumption.

The finite Gaussian identity concerns completed logical blocks. Individual
processor pulses or intermediate routes need not preserve the physical
charge interpretation. No polar projection, covariance repair, clipping of
probabilities, or selection of successful detector outcomes is applied.
Floating residuals remain diagnostics rather than physical precision bounds.

#### 5.2 Continuing measurement protocol

The two clock references have populations $(1/6,5/6)$, contrast
$r=-2/3$, and modular gap $\Omega=\log5$. Their parameter is programmed
modular flow; it is not calibrated laboratory seconds.

One history prepares the field and clocks once, advances them for
$\tau=1/4$, and records three lepton-mode occupations. A committed result
can trigger a neutral phase on another block. The same field and clocks
continue for another $\tau$, after which a further occupation and two
clock quadratures are read. No field reset occurs. All 64 six-bit terminal
records and the 16 full-field Gaussian conditionals are retained.

For the initial empty-pair apparatus state,

```math
P_i(\text{empty at }\tau)-P_i(\text{occupied at }\tau)
=\cos(2m_i\tau).
```

The final clock quadratures recover

```math
\theta=\mathrm{atan2}\!\left(
\frac{\langle Y\rangle}{r},\frac{\langle X\rangle}{r}\right)
=2\Omega\tau.
```

Consequently the observer decoder uses only records, calibration and the
admitted phase branch:

```math
\frac{m_i}{\Omega}
=\frac{\arccos\langle(-1)^{n_i(\tau)}\rangle}{\theta}.
```

The prior branch condition matters: cosine data cannot exclude arbitrary
frequency aliases. For this trial, exact elementary bounds give
$J<1$, $\Omega<2$, $\max\lambda(F)<1/4$, hence
$2m_i\tau<1/4<\pi$. The decoder does not consult the hidden Hamiltonian
or simulator duration to obtain the displayed ratios.

The recovered ratios, in **native family-axis order**, are

| Native axis | Measured trial $m_i/\Omega$ |
| --- | ---: |
| 0 | 0.1080141452 |
| 1 | 0.1359561376 |
| 2 | 0.0907450442 |

The source and independent continuous probe outputs agree at numerical
residuals of order $10^{-14}$; full-field and probe probability marginals
agree at order $10^{-15}$. An independent state-level reduction also
compares all conditional six-mode covariances and clock references.
These are numerical checks, not outward roundoff or physical-service
certificates. Finite-shot confidence and finite native preparation errors
are not included in the point inverse.

The preparation is an apparatus state. It does not show that native repair
prepares the physical interacting vacuum. The measured frequencies agree
with the trial's charged quasiparticle gaps under its stated quadratic and
particle–hole interpretation.

#### 5.3 Protected nonzero-momentum continuation

The [protected momentum implementation](../../model/common_source/chiral_candidate/protected_momentum.py)
and its [retained output](../../model/common_source/chiral_candidate/protected_momentum_results.json)
retain one zero momentum and the three pairs `±(1/3,0,0)`,
`±(0,1/3,0)`, `±(0,0,1/3)` in one declared periodic coordinate box.
Each momentum carries all 90 Spin, internal and family modes: **630 CAR
modes in total**, with the zero mode included once. The fixed quadratic
background factors into 96 independent Fock blocks of at most eight modes.
Actual source matrices execute every input column of each block, and an
independent full 90/180-mode Gaussian calculation checks the factorization.
The retained complex factors specify every successful field/clock/reference
output without materializing a vector of dimension `2^630`.

Each occupation is the two-level factor of the exact six-state native
record-correctable code; its three-level spectator remains untouched.
Recorded C/E correction acts as identity on the complete encoded field and
arbitrary references. Source words address the decoded occupation wires,
with independently checked routing. This uses many carriers and active
coherent recovery; it does not contradict the single-carrier capacity bound
or turn the uncorrected dissipative response into a unitary field.

One protocol prepares once, evolves, reads an occupation, applies recorded
feedback, evolves again and reads again. All selector timeouts and disturbed
endpoints remain part of its output. One fixed charged-lepton component is
probed; all other components remain in the continuing field. A second
protocol repeats the same fixed elementary evolution word twice per interval. Its first-read marginal
cancels the nonzero-momentum visibility in the ideal response:

```math
P(t)=A\sin^2(E t),\qquad
\frac{P(2t)}{4P(t)}=\cos^2(E t).
```

The decoder uses those two recorded marginals and the final quantum-clock
quadratures under a declared alias window. It does not read Hamiltonian
eigenvalues. The three orientations test the same relation
`(E(p)^2-E(0)^2)/Omega^2=|p|^2/4` of the supplied Weyl–Higgs symbol.
Pointwise source comparisons are numerical diagnostics; a certified
frequency-inversion error or a physical particle pole is not established.
The budget charges all 630 native carriers and source/service calls; the
six-mode construction's smaller error admission is not transferred.
Export-failure blocks have their own flagged error bound, rather than
being discarded or represented as explicitly retained success vectors.
Floating composition diagnostics are not outward or hardware certificates.

This supplies protected propagation and continuing observation for the
chosen finite matter field. The coordinate box, cutoff, coherent services,
apparatus state, programmed modular clock and trial coupling law remain
construction choices. The measured hierarchy still fails; interacting
quantum gauge/Higgs dynamics, a complete chiral gauge theory, physical spacetime
and the joint gravity attachment remain open.

#### 5.4 Finite overlap kinetic regulator

The [overlap operator](../../model/common_source/overlap_regulator.py)
and [replay receipt](../../model/common_source/overlap_regulator_results.json)
use the same retained $M_{15}=\Lambda^2\mathbb C^5\oplus
\Lambda^4\mathbb C^5$ representation on a periodic four-dimensional
Euclidean lattice. All four Wilson transports and a nonconstant
$SU(3)\times SU(2)\times U(1)$ link background are implemented. The
modified left projector obeys the finite Ginsparg--Wilson, Weyl-map and
gauge-covariance checks. All 16 free corners leave exactly one massless Dirac
kernel species, from which the modified left projector supplies one Weyl
kinetic sector. A small-link domain gives a volume-independent ideal local
approximation tail $16(15/17)^{m+1}$ at range $2m+1$; its free symbol has
error at most $a|p|^2$ on $a|p|\le1/4$. Two complete native queries retain
both output branches and all four spinor input columns.

This constructs a finite chiral **kinetic** regulator. The gauge background,
spacing, three-family replication and source service accuracy are supplied.
Its global chiral fermion measure, quantum gauge/Higgs dynamics, interacting
continuum, Lorentzian continuation and link to the protected 630-mode
many-fermion history are not established. The failed Yukawa hierarchy is
unchanged.

#### 5.5 Local completion and a second-shell discriminator

The [identifiability check](../../model/common_source/chiral_candidate/dispersion_identifiability.py)
constructs a local two-hop six-cell-per-axis lattice symbol that equals the
continuum symbol at every momentum in the protected seven-mode field and has
the same derivative at zero. With the same internal matrices, trial mass,
current and clock, every restricted seven-mode quantum history agrees.
The larger local lattice carries additional modes, states and resource
requirements, and still has eight free kinetic nodes with balanced
handedness. The restricted data therefore do not select a unique spatial
completion or a chiral gauge regulator.

At coordinate momentum $2/3$, the two-hop kinetic coefficient divided by
the continuum speed is $1-\sqrt3/\pi$, versus $2/3$ for the continuum law.
The [second-shell receipt](../../model/common_source/chiral_candidate/dispersion_identifiability_results.json)
retains two number reads, feedback, a quantum clock, a neutral pair
quadrature and all terminal outcomes. Exact joint-event intervals for the
two fixed laws separate. The source replay checks that contrast; its source
arithmetic admission and the stipulated event implementation error are not
physical apparatus guarantees. The full 810-mode native correction/export
budget and the larger lattice's other modes are not executed. The additional
momentum read is a concrete discriminator for the proposed field attachment.

#### 5.6 Calibrated action and matching control

The separate [calibrated matter action](../../model/common_source/calibrated_matter.py)
uses the same internal representation, Higgs tensors and current matrices
with one explicit empirical parameter packet. All nine charged Yukawa
coordinates, the Higgs expectation value and the three gauge couplings share
a specified $\overline{\mathrm{MS}}$ boundary at $91.1876$ GeV. The actual
complex down-family matrix and its all-left conjugation convention produce
the supplied unitary CKM matrix in the charged current. Exact tensor, mass
and Weyl-dispersion checks retain the charged color/spin multiplicities.
These identities establish that this representation can carry the calibrated
mass pattern; they do not select its coefficients from the repair rule.

The [matching receipt](../../model/common_source/calibrated_matching_results.json)
records a separate execution of SMDR 1.3, including its source revision,
input conventions, inverse calibration and independent forward call.
The nine charged-fermion masses and the Higgs and Z calibration inputs
return within their quoted marginal standard uncertainties in that numerical
run. Those values are **external Standard Model outputs, not native observer
records**. Running quark masses retain their scale and scheme; unstable
particle masses retain their pole or Breit--Wigner convention. The unfitted
W output is $80.3491825273$ GeV, compared with the selected reference
$80.3692\pm0.0133$ GeV: a $-20.0175$ MeV difference, or $-1.5051$ quoted
standard uncertainties. It misses a one-standard-uncertainty criterion.
No combined input/theory covariance or OPH-specific empirical success follows.
The external matching neglects CKM effects at its declared approximation.

A rounded published parameter table is retained as an input-precision
control: it misses the electron and muon measurement precision. The fixed
full-precision packet resolves that calibration issue; it supplies no new
parameter derivation. The bosonic tree check explicitly chooses
$m^2_{\mathrm{tree}}=-\lambda_{\mathrm{SM}}v^2$ for its comparison potential;
this is not SMDR's loop-matched quadratic coefficient or a computation of
the full quantum effective-potential Hessian. The optional Weinberg sector
uses NuFIT mixing and mass-squared splittings, normal ordering, zero lightest
mass and zero Majorana phases. Its unmeasured choices and unexecuted
low-energy-to-common-scale matching remain explicit. The original native
candidate imports none of these calibration inputs.

#### 5.7 Calibrated masses read from simulated observer records

The separate [observer producer](../../model/common_source/calibrated_spectrum_observer.py)
executes the rest-field source factors and retains complete simulated detector
counts in its [receipt](../../model/common_source/calibrated_spectrum_observer_results.json)
and compressed array artifact. Its blind inverse reads the counts, clock
convention and error contract; it cannot read the parameter packet or call the
source generator. All nine supplied running-action gaps lie in its intervals:

| Field coordinate | Approximate record interval, GeV at the common MZ boundary |
| --- | --- |
| up | [0.00123406519, 0.00123679893] |
| charm | [0.624245749, 0.624253591] |
| top | [169.888809, 169.890205] |
| down | [0.00270061637, 0.00270336216] |
| strange | [0.0537463744, 0.0537495387] |
| bottom | [2.86239022, 2.86241641] |
| electron | [0.000486343406, 0.000489071040] |
| muon | [0.102758037, 0.102761603] |
| tau | [1.74550466, 1.74552170] |

The receipt contains exact rational endpoints; the displayed decimals are
rounded summaries. These are running-action coordinates, not the experimental
pole or low-scale quark masses returned by external matching. Their widths
range from about 2.73e-6 to 1.40e-3 GeV. In particular, this sampling campaign
does not recover electron or muon measurement precision.

There are 21 dyadic phase levels, two quadratures and two million complete
trials per setting: **84 million simulated trials**. Each trial retains the
selector address, nine field bits and two clock bits, including timeout
addresses. Adaptive analyzer offsets depend only on earlier records. A shared
clock gap `log(5)` and the single supplied conversion `1 GeV/log(5)` set the
units; the programmed modular parameter is not laboratory seconds. Exact
interval inversion and a simultaneous conditional sampling bound of about
1.17e-8 cover all 420 real statistics, under independent simulated trials and
the declared numerical, frame and physical-service bias allowances.

The source execution uses all four Fock input columns of 42 commuting charged
pair factors and retains six neutral zero modes. Actual CKM source frames and
a complete 64-column weak down-sector conjugation test the change of basis.
No dense 2^90-state source output is claimed; the complete weak-frame relation
uses a separate uniform error admission, not an outward analytic certificate.
Repeating the declared programs for all sampled trials would charge
341,796,000,000 native pulses and 683,592,000,000 load/export events. Those are
required resources, not an executed physical protection history.

This establishes a **calibrated source/readout sanity check**. The measured
inputs still determine the action. Interacting self-energies, confinement,
scattering, physical clock calibration and a source-derived mass law remain
unproved. Replay the retained counts and their source channels with:

```sh
python -W error -m model.common_source.calibrated_spectrum_observer --check model/common_source/calibrated_spectrum_observer_results.json
```

### 6. Relation to observers' effective 3+1D spacetime

There are two distinct levels of attachment.

**At the declared field-description level**, the action is a 3+1D
Lorentzian Weyl–Higgs candidate with the dispersion and chiral representation
above. The prior OPH metric construction also supplies an exact pointwise
comparison: for spatial principal tensor $A=L L^T$, use coordinate
principal momentum $q=L^T p$. For the retained metric with lapse $N$,
proper orthonormal spatial momentum is $\hat p=q/N$. The Weyl principal
symbol then shares the scalar characteristic cone. The native four-component
Clifford matrices can be put into chiral form by

```math
W=\frac1{\sqrt2}\begin{pmatrix}I&-I\\I&I\end{pmatrix},\qquad
W\alpha_jW^\dagger=\mathrm{diag}(-\sigma_j,+\sigma_j),\qquad
W\beta W^\dagger=\begin{pmatrix}0&I\\I&0\end{pmatrix}.
```

A proper mass carries the coordinate lapse factor
$K_{\mathrm{coordinate}}=N K_{\mathrm{proper}}$. Clock regraduation must transform
the kinetic, mass and clock rates together.

**At the same reconstructed-history level**, the global attachment is not
complete. The previous source spatial CAR construction uses an incidence
operator with 720 node and 2,850 edge modes, including 2,130 persistent
edge-cycle zero modes. It cannot be identified by a unitary basis change
with a minimal two-component Weyl field on those 720 sites. Discarding the
cycles would discard retained outputs.

There is also a terminology distinction: the earlier scalar effective
description's **Weyl operators are bosonic CCR displacement operators**.
That theorem is not already a theorem about left-Weyl fermions. Its scalar
smoothing controls do not supply the missing chiral spin-connection or
field-state intertwiner.

The old radial domain check gives another concrete requirement: imposing
scalar Dirichlet conditions on every component does not make the first-order
Weyl operator self-adjoint. A spin-compatible domain, boundary conditions,
connection and reference must be supplied and checked. On curved space,
the Dirac square also contains spin-connection and curvature contributions;
it is not simply the scalar spatial operator acting on each component.
The frozen comparison coefficient is not asserted to equal the limiting
native tensor. A common modular reference also needs its state/clock map:
under a time regraduation, $\beta H$ must remain fixed. One cannot hold
both reference conventions fixed and freely rescale only the field rates.

Thus a shared characteristic cone is available as an exact local
compatibility result. A complete chiral field history on the observers'
recovered spacetime has not been established. The
[existing geometric construction](../../model/observer_particles/geometric_matter_join.py)
and [new attachment check](../../model/common_source/chiral_candidate/spacetime_attachment.py)
state that boundary explicitly.

The attachment check also uses the **same native incident clock** as the
collective detector. Conditioning on one successful incident event while
retaining the separate selector and deadline timeout probabilities, the
ten-mode charged response has exact operator norm and spectral radius
$7/10$ under C and $2/5$ under E. The complete Ramsey histories reproduce
those contrast factors. For any linear map $J$ to an unmeasured closed
quadratic Weyl–Higgs propagator $U$, the measured response $R$ obeys
$\lVert UJ-JR\rVert_{\mathrm{HS}}\ge(1-\lVert R\rVert_2)\lVert J\rVert_{\mathrm{HS}}$.
Thus no nonzero exact intertwiner attaches that **record-conditioned response**
to the declared closed unitary evolution. Even at one incident tick, the
site-$2$ to site-$7$ conditional response is zero with a one-proposal
deadline and $1/540$ with a two-proposal deadline: the waiting record matters.
A field coupled to a measured environment can contract; a complete dephasing
control verifies that logical possibility without identifying a physical
detector interaction or a Standard Model field.

A finite spectral menu is also not a completed local chiral gauge
regulator. The naive lattice control has eight Weyl nodes with balanced
handedness. Relevant published regulator routes and their hypotheses are
collected in [REGULATOR_ROUTES.md](../../model/common_source/chiral_candidate/REGULATOR_ROUTES.md).

### 7. The hierarchy tests and their outcome

#### 7.1 Linear native family map

The three trial Yukawa matrices are proportional to the same $F$.
Therefore all three charged channels have the same family ratios. The
family spread is strictly below $3/2$; the spread across **all** charged
tree masses is strictly below $15/4$. The canonical tree quark mixing is
identity, up to basis ordering and rephasings, with no generated CP phase.

Those are exact exclusions of this trial prescription as an observed
spectrum. A common Higgs scale or a permutation of the modes cannot change
the ratios. This conclusion does not depend on a costly QCD calculation.

#### 7.2 Fixed nonlinear controls

To test whether a nonperturbative numerical dependence could amplify the
small native anisotropy, three further prescriptions were fixed without
fitting observed masses. They evaluate the half-filled one-dimensional
Hubbard charge-threshold function with $t=\Omega/2$ and
$U=a\lambda_i(F),b\lambda_i(F),J\lambda_i(F)$.

The underlying [Lieb–Wu theorem](https://arxiv.org/html/cond-mat/0207529)
and [positive threshold representation](https://arxiv.org/html/cond-mat/9911281)
support a quadrature-free rational enclosure. Here that function is only a
**declared nonlinear texture prescription**. Its use does not turn a
one-dimensional gap theorem into a three-dimensional dynamical mass theorem.
The very small thresholds are not claimed as resolved measurements in the
finite clock history.

| Prescription | Ascending ratios, lightest normalized to one |
| --- | --- |
| Linear $F$ | $1:1.1903:1.4982$ |
| Hubbard function, $U=a\lambda_i$ | $1:204.97:65{,}743.94$ |
| Hubbard function, $U=b\lambda_i$ | $1:353{,}640.05:3.5494\times10^{11}$ |
| Hubbard function, $U=J\lambda_i$ | $1:8{,}075.31:1.3683\times10^8$ |

The source predictions are computed before a separate comparison with the
[PDG lepton summary](https://pdg.lbl.gov/2025/download/db2024.pdf). The charged
leptons require approximately $1:206.77:3477.23$. The benchmark uses broad
windows $[206,207]$ and $[3470,3485]$, so the rejection is not an argument
about the last measured digit.

**All four fixed prescriptions fail the joint two-ratio test.** The nearby
middle ratio in the first nonlinear control cannot compensate for its much
larger heaviest ratio. No continuous parameter search or empirical mass fit
was performed. Other nonlinear couplings and interacting repair actions are
not excluded by these tests.

### 8. What is established, and what remains

| Item | Result in this candidate |
| --- | --- |
| Internal chiral matter types and Higgs channels | Explicit operators inherited and rechecked under their stated lift |
| Native family-dependent operator | Exact, covariant, nonscalar, with three distinct eigenvalues |
| Weyl/Higgs Spin, charge and dispersion | Checked on the same declared action and Nambu conventions |
| Weak-chirality discriminator | Correct operator vertex weights; no dynamical weak scattering calculation |
| Continuing source readout | Full 90-mode rest field and protected 630-mode seven-momentum field; complete factorized complex outputs and flagged failures |
| Charged rest-gap ratios and propagation | Rest ratios and three-axis moving response read from source records and quantum-clock quadratures; physical identification remains open |
| Spatial completion from seven momentum modes | Nonunique: a distinct two-hop local symbol agrees on the restricted histories; a second-shell read separates the fixed alternatives |
| Observed mass hierarchy | Fails for every fixed prescription tested here |
| Separate calibrated action | Carries supplied charged masses/mixing; external SMDR matching checks calibration, with an unfitted W miss. No source derivation or full physical observer attachment follows. |
| Quantum gauge/Higgs dynamics | Not executed by this fixed-background quadratic trial |
| Local chiral kinetic regulator | Finite Euclidean overlap operator, modified projector and uniform ideal locality tail on an admissible link domain; full quantum chiral theory open |
| Same recovered-spacetime field history | Local principal-symbol compatibility; full domain/history attachment open |
| Physical units and observed particle identification | Not established |
| GR entropy stationarity and responsive geometry | Physical entropy-area/stress dictionary and responsive geometry remain open |

The specific remaining mass problem is to determine and test the
mass-generating coupling of available native family data to the Higgs
channels, their kinetic forms and the physical vacuum. The current example
shows that getting the chiral representation right and computing every
number from a chosen source prescription does not guarantee the measured
hierarchy. A new prescription must earn that agreement through its actual
charged response, including its extra modes and error accounting.

### 9. Reproduction and evidence

From the public repository root, using Python 3.12 and the root
[pinned requirements](../../requirements.txt):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/chiral_candidate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -u -m model.common_source.chiral_candidate.run --output /tmp/chiral-candidate.json
.venv/bin/python -m model.common_source.chiral_candidate.run --verify model/common_source/chiral_candidate/results.json
.venv/bin/python -m model.common_source.chiral_candidate.protected_momentum --check model/common_source/chiral_candidate/protected_momentum_results.json
.venv/bin/python -m model.common_source.chiral_candidate.dispersion_identifiability --check model/common_source/chiral_candidate/dispersion_identifiability_results.json
.venv/bin/python -m model.common_source.overlap_regulator --check model/common_source/overlap_regulator_results.json
```

The replay requires fresh JSON and companion NPZ paths and refuses to
overwrite them. The NPZ retains all literal source work/reference outputs
and the complete Gaussian conditional descriptors with their separate branch
weights. Verification checks source custody and retained state/weight
consistency; it does not rerun the dynamics.

The retained calculation binds the source code, native reference, exact
representation and interaction artifacts, and prior response evidence.
Exact algebra, rational interval certificates and floating execution checks
remain distinct forms of support. The protected propagation receipt has its
own source bindings and complete complex factor artifact; it preserves the
distinct rest experiment and passive repair-response contracts.


## Hierarchy diagnosis

The immediate failure is in the candidate's **mass-generating prescription**.
The observer reads the implemented quadratic action correctly, but that
action makes masses proportional to a nearly degenerate occupation feature.
Its fixed Higgs background supplies no further interaction that could turn
those coefficients into the measured hierarchy. This diagnoses a rejected
architecture, not a contradiction of every OPH construction.

This diagnosis extends the
[candidate construction](#candidate-construction).
The underlying source dynamics and retained quantum-state artifact are the
same ones used in that construction.
New exact calculations, tests and the reproducible diagnosis are in
[hierarchy_diagnosis](../../model/common_source/hierarchy_diagnosis/README.md).

### 1. Locate the failure before changing the readout

The declared action uses

```math
Y_u=aF,\qquad Y_d=bF,\qquad Y_e=JF,\qquad H=(0,v),\quad v=\Omega.
```

Here the native source computes

```math
a=\log(25/8),\quad b=\log(8/5),\quad J=a-b,\quad\Omega=a+b,
```

and the eigenvalues of the actual family feature are

```math
f_1=\frac{61}{450},\qquad
f_2=\frac{41}{225}-\frac{7\sqrt5}{750},\qquad
f_3=\frac{41}{225}+\frac{7\sqrt5}{750}.
```

The step that turns this occupation feature into a left–right mass coupling
is a declared choice. The source derivation of $F$, internal chiral
representation and invariant Higgs tensors does not prove that
$Y_e=JF$ is the physical coupling.

With canonical kinetic terms the charged-lepton mass parameters are
$m_i=vJf_i$. Their ascending ratios are approximately
$1:1.1903:1.4982$. The retained clock records independently recover those
same ratios. The observed target is approximately $1:206.77:3477.23$;
the audit uses the much wider windows $[206,207]$ and $[3470,3485]$,
so experimental uncertainty is immaterial to this exclusion. These are
comparison inputs only. See the PDG
[muon](https://pdgweb.lbl.gov/2025/listings/rpp2025-list-muon.pdf) and
[tau](https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf) listings.

The current 90-mode field law is quadratic in CAR operators once the Higgs
background is fixed. The earlier four-mode response construction had an
explicit density–density interaction. The current chiral extension uses
the number $J$ as a coefficient of a bilinear pairing; it does not retain
that four-mode interacting dynamics as a field interaction. Nor does it
execute quantum Higgs or gauge fluctuations. A logarithmic statistic from
the native reference is not, by itself, a derivation of a fermion self-energy.

Consequently there is no uncomputed many-body mass enhancement inside this
particular quadratic Hamiltonian. Its dispersion is already exactly
$E_i(p)^2=c^2|p|^2+m_i^2$. A genuinely interacting extension could change
that conclusion, but would be a new dynamics to derive and test.

### 2. The clock, pair factor and vacuum cannot repair this action

For a selected rest pair the implemented Hamiltonian is

```math
H_d=d\,a^\dagger b^\dagger+d^*ba.
```

Its four Fock energies are $-|d|,0,0,+|d|$. A charged odd-CAR excitation
above its ground state costs $|d|$, whereas an initially empty pair has
occupation oscillation frequency $2|d|$. The observer inverse already
accounts for this factor of two. The normal CAR response has nonzero
residues at both signed poles; the three family gaps are not hidden by
vanishing detector overlap.

More generally a quadratic Nambu field satisfies

```math
\Xi(t)=e^{-i\mathcal B t}\Xi(0),\qquad
\{\Xi(t),\Xi^\dagger(0)\}=e^{-i\mathcal B t}.
```

Thus the retarded linear-field response has state-independent pole
locations for this fixed law. Changing the reference state can change
other correlation functions and occupations; it cannot create a different
mass spectrum in this response. This statement does not extend to a
state-dependent interacting self-energy. The time parameter here is the
declared coherent model coordinate, compared with the calibrated modular
clock in the retained experiment, not laboratory seconds.

A common lapse, clock regraduation, Higgs normalization or conversion to
physical units multiplies all three masses by the same factor. It leaves
their ratios unchanged. A change of field basis also preserves poles when
both kinetic and mass forms are transformed consistently.
The native band Gram matrix is exactly $(10+2\sqrt5)I_3$, so the
implemented family frame conceals no anisotropic Gram normalization.
Measuring all families at the same nonzero momentum would compress their
energy ratios toward one, not enlarge the rest-mass hierarchy.

### 3. A stronger exclusion: every affine law built from this feature fails

The result is broader than failure of the coefficients $a,b,J$. Let

```math
Y=uI+vF,\qquad u,v\in\mathbb C,
```

with canonical kinetic normalization and one common Higgs factor. Since
$F$ is Hermitian, its complex eigenvalues under this map are
$z_i=u+vf_i$, and the masses are proportional to $|z_i|$.
Writing $\phi=(1+\sqrt5)/2$, the source eigenvalues obey

```math
\frac{f_3-f_1}{f_2-f_1}=\phi^2,
\qquad z_3-\phi^2z_2+\phi z_1=0.
```

The triangle inequality, allowing every assignment of the largest modulus
to the three coefficients, gives the necessary bound

```math
\boxed{m_{\max}\leq\phi^2m_{\mathrm{middle}}+\phi m_{\min}.}
```

To see that it covers every assignment, the absolute coefficients in the
linear relation are $1,\phi,\phi^2$. Placing the largest mass on the
smallest coefficient gives the weakest constraint; assigning the larger
remaining coefficient to the middle mass maximizes its right-hand side.
Any other placement gives a smaller bound. Complex phases do not evade
the triangle inequality.

For the retained broad comparison windows this requires

```math
\frac{m_{\max}}{m_{\min}}
\leq207\phi^2+\phi=543.55106966\ldots,
```

while the observed lower bound is 3470. Therefore **the entire complex
affine class fails**, including arbitrary identity shifts, overall scales,
sign changes and phases. This is an exact exclusion of that class, not a
numerical search.

It also explains why merely changing the same protected-core reference
weights is insufficient. If $w_q$ is the weight per configuration with
core occupation $q$, direct counting gives

```math
\mathbb E[n_{\mathrm{core}}]=10w_1+w_2,\qquad
\mathbb E[n_{\mathrm{outer}}]=9w_0+2w_1.
```

For the same core and band, its compressed mean is always

```math
(9w_0+2w_1)I+(8w_1+w_2-9w_0)S_{\mathrm{core}},
```

where $S_{\mathrm{core}}$ is the compressed core projector. The present
$F$ is a nonscalar affine function of that same projector, so this whole
family of mean-based linear prescriptions remains in the excluded class.
This argument does not cover arbitrary changes to the carrier, core type,
observables or interaction law.

### 4. Nonlinearity can amplify masses without fixing their pattern

The three previously tested Hubbard-function textures are fixed
prescriptions, not fits. All amplify the source eigenvalue differences,
but all give the wrong relative spacing:

| Prescription | Middle/light | Heavy/light | Heavy/middle |
| --- | ---: | ---: | ---: |
| Linear $F$ | 1.1903 | 1.4982 | 1.2587 |
| Hubbard function with $U=af_i$ | 204.975 | 65,743.94 | 320.74 |
| Hubbard function with $U=bf_i$ | 353,640.05 | $3.5494\times10^{11}$ | $1.0037\times10^6$ |
| Hubbard function with $U=Jf_i$ | 8,075.31 | $1.3683\times10^8$ | 16,943.78 |
| Observed charged leptons, approximately | 206.77 | 3,477.23 | 16.82 |

The near-middle ratio in the $a$ control is accompanied by a heavy mass
about nineteen times too large after light-mass normalization. A common
scale cannot correct it. A useful scale-independent shape statistic is

```math
Q=\frac{\log(m_3/m_2)}{\log(m_2/m_1)}.
```

The controls give $Q$ near 1.082–1.084, while the broad empirical windows
give approximately 0.529–0.531. Certified rational logarithm bounds make
these intervals disjoint. The wrong *shape*, not inadequate arithmetic
precision, is the issue.

The mechanism is understandable: the weak-coupling one-dimensional
Hubbard threshold contains an exponential in $-1/U$, so using
$U=sf_i$ amplifies the reciprocal eigenvalue spacings in a specific
pattern. The source theorem is for the half-filled one-dimensional
Hubbard chain, as discussed by
[Jeckelmann, Gebhard and Essler](https://arxiv.org/html/cond-mat/9911281).
Its use as a numerical texture does not derive an interacting 3+1D chiral
mass. These fixed tests do not exclude all nonlinear laws or all values
of a continuously varied coupling.

An arbitrary function $f(F)$ is also not a prediction by itself. Since
$F$ has three distinct eigenvalues, quadratic polynomial interpolation
can assign any three target values to them. Those polynomial coefficients
would encode the desired answer unless an independently specified source
law determines them.

### 5. Small normalization corrections cannot hide the discrepancy

Let $r=\kappa_2(M)<3/2$ be the current mass-matrix condition number. For
new positive constant kinetic forms $Z_L,Z_R$, canonical normalization
gives, in the usual Dirac convention,

```math
M_c=Z_L^{-1/2}MZ_R^{-1/2},\qquad
\kappa_2(M_c)\leq r\sqrt{\kappa_2(Z_L)\kappa_2(Z_R)}.
```

The transpose in an all-left convention does not change this norm bound.
Reaching a heavy/light ratio of at least 3470 requires the product of
kinetic condition numbers to exceed $5.35\times10^6$. If both
$\|Z-I\|_2\leq\epsilon<1$, it requires
$\epsilon>6937/6943$, approximately 0.999136.

Similarly, a constant additive mass correction with
$\|\delta M\|_2=\delta\,s_{\min}(M)$, $\delta<1$, obeys

```math
\kappa_2(M+\delta M)\leq\frac{r+\delta}{1-\delta}.
```

The same target requires $\delta>6937/6942$, approximately 0.999280.
These are necessary bounds, not sufficient constructions or claims about
the size of every possible quantum correction. They rule out explanations
that leave the current action almost unchanged in the specified norms.
Strong dynamics, additional fields and momentum-dependent self-energies
remain outside this argument.

### 6. Family structure must be recovered from more than one mean matrix

Every separate function of $F$ shares its eigenvectors. Consequently
$Y_uY_u^\dagger$ and $Y_dY_d^\dagger$ commute when the Yukawa
matrices are functions of this one Hermitian feature. For nondegenerate
masses and canonical common kinetic terms their CKM matrix is trivial
up to phases and permutations. Choosing more elaborate eigenvalue
functions cannot supply physical quark mixing.

The new audit also computes the full connected occupation covariance

```math
C_{ij}=\mathbb E[n_in_j]-\mathbb E[n_i]\mathbb E[n_j]
```

from all 66 configurations. Its full rank is 11, with the single constant
null direction required by fixed total load. Its compression to the actual
three-dimensional band is positive definite, has spread below $5/4$,
and **commutes exactly with $F$**. Thus replacing the mean by this
equal-reference covariance does not produce either a large linear
hierarchy or an independent family basis. This is a computed negative
control, not a claim that every dynamical response is determined by these
two matrices.

There is a further symmetry restriction on using only this stationary
reference. Its native stabilizer contains a half-turn represented on the
family band by $R=\mathrm{diag}(-1,+1,-1)$. Any covariant
endomorphism built solely from that reference and rotation-invariant
carrier data must commute with $R$. For the same family action on both
fermion factors and a family-neutral Higgs, the allowed matrices have a
common two-dimensional block and an unmixed one-dimensional block.
They can support at most two-family mixing, with no three-family CKM
CP violation. This restriction is conditional on those inputs; it is not
a symmetry theorem for every full repair context.

#### The existing update contains information that the stationary feature loses

The actual
[reference preparation](../../model/native_modular_reference/local_preparation.py)
already uses two distinct completion caps, $(0,1,5,7)$ and
$(0,1,2,8)$. The half-turn fixes the first and moves the second. Its full
update context therefore has less symmetry than its stationary reference.
The new
[native response audit](../../model/common_source/hierarchy_diagnosis/native_family_dynamics.py)
uses the existing finite selector and endpoint repair, without adding a
mass-dependent parameter.

It enumerates the reduced outcomes on every one of the 66 binary inputs,
including the selector's terminal timeout. Each accepted outcome retains
its seam, coin, current, output state and exact probability. The timeout
retains unchanged work with a separate terminal flag. The complete
one-proposal mean response is

```math
T=\sum_e w_e E_e+\tau I,
\qquad \sum_e w_e+\tau=1,
```

where $E_e$ averages the two endpoint occupations. Its diagonal native
reference is exactly stationary. With orthonormal band frame
$U=B/\sqrt{10+2\sqrt5}$, define the compressed response
$K=U^TTU$. The calculation finds

```math
[F,K]\ne0,\qquad [R,K]\ne0.
```

All three off-diagonal links of $K$ are nonzero. Since $F$ has simple
eigenvalues, its spectral projectors $P_i$ are polynomials in $F$.
Each $P_iKP_j/K_{ij}$ is a matrix unit. Therefore $F,K$ generate
the full three-by-three matrix algebra, an exact statement rather than a
numerical rank guess.

This is a concrete positive result: **the existing update context supplies
independent family structure omitted by the mean-based mass ansatz**.
It does not select Yukawa matrices or predict mixing angles. These
retained matrices are real and do not by themselves derive a physical
CP phase. The caps and preparation protocol also remain declared
operational inputs, not a selected autonomous matter dynamics.

There is a second important limitation: the band is not closed under
this update. For $V=(I-UU^T)TU$, the exact calculation gives

```math
V\ne0,\qquad U^TT^2U-K^2=V^TV\ne0.
```

Thus iterating only the compressed matrix would lose the complement's
response. The displayed identity is an algebraic marginal comparison;
it is **not** a protocol that restarts after a terminal timeout. A continuing
reduction must retain the complementary modes or their derived memory,
along with the actual event and timeout rules. The parameter of this
calculation is one completion proposal, not physical time, modular time
or a calibrated count of elapsed native attempts. No eigenvalue of this
relaxing response has been identified with a Lorentzian particle mass.

### 7. What the next mass calculation must actually derive

The next target is the charged response of a specified continuing,
interacting repair construction. It must determine the effective mass
coupling and kinetic response together. Schematically its charged inverse
response has a kinetic part, a Higgs-dependent part and any interaction
self-energy; physical rest energies follow its singularities in the
observer's calibrated clock. An occupation statistic alone does not
supply that complete response.

Three checks make this a concrete research target:

1. Start from the existing update's independent family response and retain
   its complementary modes or derived memory. Derive the actual coupling to charged
   left and conjugated-right fields, rather than choose an eigenvalue
   function to match masses.
2. Compute the same law's kinetic normalization, background and charged
   response, with its state, regulator and clock assumptions explicit.
   Compare action parameters with observable poles or thresholds in the
   appropriate domain; native relaxation rates need their own bridge.
3. Test both independent charged-lepton ratios and nontrivial family
   mixing before moving to costly QCD or making physical-species claims.

This diagnosis does not require proving that one microscopic architecture
is uniquely possible. It requires a specified architecture to produce its
mass prediction without using the target masses to choose its couplings.
It also leaves intact the separate obligation to attach a complete chiral
field history to the observers' reconstructed spacetime; the current
frozen principal-symbol comparison is not that full attachment.


## Full repair response

This response audit follows the
[hierarchy diagnosis](#hierarchy-diagnosis)
by retaining the actual repair complement and investigating the charged
quantum response. The implementation is
[full_repair_response](../../model/common_source/full_repair_response/README.md).

**Outcome:** the measured Standard Model mass hierarchy is not recovered.
The charged rest-energy observable is calibrated in the declared Weyl–Higgs
control, including its momentum dependence and its distinction from neutral
decay. Identifying that physical field with a native continuing repair
history remains unresolved. Full complementary memory and alternative clock
readings do not repair the tested mass assignments.

The objective is to derive masses from the observer's charged response of
one specified continuing repair construction, then compare its ratios with
measurements. Replacing one trial Yukawa matrix by another source statistic
is a conditional model test; it is not automatically a physical derivation.

### Upstream audit

The concrete concern is upstream of the arithmetic: **we have computed
neutral occupation response and then assigned its coefficients to a charged
fermion mass interaction. The physical identification between those two
objects has not been established.** The earlier assignment $Y_e=JF$
and the new harmonic-kernel assignment below are explicit choices. Their
failures reject those assignments in their stated domains; they do not
invalidate unrelated OPH theorems.

| Connection | What is established in these constructions | What is still required |
| --- | --- | --- |
| Native repair to occupation response | Exact finite maps and complete population outcomes | Identification of the observer's charged particle field |
| Rank-three band to family coordinates | A three-dimensional native response band | Proof that its coordinates are three physical generations, with the appropriate poles and currents |
| Internal matter module and Higgs lines | Conditional chiral representation, current identities and invariant tensors | A source-derived physical coupling and state for those fields |
| Lorentzian description | Declared Weyl–Higgs kinetic symbol and an exact frozen local principal-symbol attachment | A complete chiral field history on the same reconstructed spacetime, with its domains and reference |
| Clock | Operational incident-clock maps and separate modular-flow constructions | The charged observable's response in that clock; a clock alone does not select its Hamiltonian |
| GR theorem stack | Conditional implications under stated entropy and geometric inputs | The same matter/history witness satisfying the outstanding inputs |

This makes the audit falsifiable. Check the microscopic law, observable
algebra, projection, state, clock and effective-field identification at each
link. Do not classify an unsuccessful conditional field attachment as a
failed prediction of the full framework. Conversely, do not count a
conditional representation or geometric compatibility result as confirmation
of all observed particle physics.

### 1. Continuing repair with every termination outcome retained

The starting source is the existing protected-core mixed-cap completion
protocol. Its finite selector has accepted work map $A$, state-independent
success weight $q=1-\tau$, and a separate terminal-timeout output
$\tau I$. The declared continuation experiment repeats accepted updates
until the first timeout. It does not restart a failed history.

After $n$ proposals the live work output is $A^n\rho$. The timeout
after $j$ successful updates retains $\tau A^j\rho$, with its own event
label and no later evolution. Their total weight is exactly one. The
normalized map $C=A/q$ is useful for reporting conditional response;
all survival weights remain part of the experiment.

This distinction prevents a false hierarchy. The eigenvalue $q$ on a
conserved work direction measures the selector's termination hazard.
Interpreting $-\log q$ as a small mass would make that mass depend on the
controller's rejection-round setting. It is not a charged rest gap.

The population model has 66 binary total-two configurations. Its accepted
first-moment response is computed on all twelve ports, rather than only
the selected three coordinates. Completion protects each of the two core
occupations separately. The port response consequently has three conserved
directions; the population response has four connected invariant sectors
of sizes $45,10,10,1$.

### 2. Exact complementary memory and the actual observable dimension

Let $B$ be the existing native band frame, $L=(B^TB)^{-1}B^T$,
$P=BL$, and $Q=I-P$. For any specified transfer $M$, its visible
coordinates satisfy

```math
y_{n+1}=R y_n+\sum_{j=0}^{n-1}M_{n-1-j}y_j
              +LMQ(QMQ)^n q_0,
```

where

```math
R=LMB,\qquad M_k=LMQ(QMQ)^kQMB,
```

and $q_0$ is the initial complementary component. No modes are removed
from the underlying source. Equivalently,

```math
L(zI-M)^{-1}B
=\left[zI-R-LMQ(zI-QMQ)^{-1}QMB\right]^{-1},
```

where the inverses exist. This is the finite block-resolvent identity, not
a new general reduction theorem; see
[Dusson, Sigal and Stamm](https://arxiv.org/html/2105.02058) for the
Feshbach–Schur framework.

For the specified population input/output, the exact reachable and
observable realization has dimension **11**, comprising two visible
conserved directions and nine decaying directions. The remaining global
total-load direction remains in the twelve-port source but has zero
residue for this band input/output. The result is a minimal representation
of an observable response, not eleven particle species.

In particular, three readout coordinates do not imply three spectral poles.
In the uniform original source the rank-three band is a degeneracy at one
relaxation eigenvalue. It was not already a theorem establishing three
distinct physical fermion masses.

### 3. Check the observer's incident clock rather than rescale a counter

The source already defines a local incident clock at the first completion
cap. Split each accepted proposal into incident map $J$ and nonincident
map $R_c$. The next-tick success and failure maps are

```math
S=J(I-R_c)^{-1},\qquad F_{\mathrm{timeout}}=\tau(I-R_c)^{-1}.
```

Finite waiting cutoffs retain the unfinished $R_c^h$ output as well as
all incident successes and terminal failures. The matrix formulas preserve
their distinct meanings.
The resolvent poles describe the infinite-wait stopping limit; the earlier
finite-deadline preparation protocol is not reexecuted in this calculation.

The incident-clock transformation need not be a common rescaling of rates:
incident and nonincident maps need not commute. We therefore compute it
explicitly. The full population response has a conditional logarithmic
decay-rate spread of about **14.81** per accepted proposal. The twelve-port
response has spread about **7.46**, and the twelve-port incident-clock
response about **5.07**. These are numerical spectral diagnostics with
their clock parameter named; they are not charged-particle masses or
certified experimental predictions.
The full 66-state incident-clock logarithmic spread is approximately 9.79.

### 4. A source-computed harmonic kernel and its kinetic derivative

The full one-proposal occupation Dirichlet form is

```math
D=I-T=\sum_e\frac{w_e}{2}(e_a-e_b)(e_a-e_b)^T\geq0,
\qquad T=A+\tau I.
```

It is an algebraic form of that single-proposal work marginal. Iterating
$T$ would discard the terminal status and is not our continuing protocol.
The complement-retaining inverse response is

```math
\mathcal S(z)=\left[(B^TB)^{-1}B^T(D+zI)^{-1}B\right]^{-1}.
```

Its zero-frequency limit and first derivative are computed exactly. With
$G=10+2\sqrt5$, the harmonic extension $U_h$ satisfies
$B^TU_h=GI$, and

```math
\mathcal S(0)=\mathrm{diag}(0,0,\lambda),\quad
\lambda\simeq0.03717449835,
```

```math
Z=\mathcal S'(0)=U_h^TU_h/G
  =I+(U_h-B)^T(U_h-B)/G\geq I.
```

The two zero directions follow from the protected core, not numerical
precision. Positive canonical normalization preserves them.

As a conditional field test, insert this kernel in the existing exact
Higgs tensors through $Y_{u,d,e}=(a,b,J)\mathcal S(0)$, with the computed
common kinetic form $Z$. The actual tensor insertion passes 72 current
identities and charge/antisymmetry checks. The internal mass kernel has
rank 14 and the spin-pairing kernel rank 28. Each charged sector has
family ratios **$0:0:1$**: seven massive Dirac fields counting color,
fourteen massless charged Dirac fields, and three massless left-Weyl
neutrinos in the minimal module.

This is a valid gauge-consistent trial attachment. It does **not** recover
the observed nonzero light-family masses. The exact memory or kinetic
derivative cannot lift the protected static zeros within this law.

A second existing source law provides a useful control: uniform repair on
all thirty seams removes the protected preparation boundary and closes
the original band. It gives

```math
\mathcal S(z)=\left(z+\frac{5-\sqrt5}{60}\right)I,\qquad Z=I,
```

hence three degenerate coefficients under the same kind of attachment.
Neither source was adjusted to observed masses.

The central identification remains open in both tests. $\mathcal S$ is a
Laplace inverse response of **neutral occupation means**. It is not yet
the observer's charged CAR inverse propagator. Its conserved directions
are not established massless fermions until that additional field
attachment is made. A positive transfer can define a separate quantum
Hamiltonian under an explicit reconstruction, as in the
[stochastic-matrix-form construction](https://arxiv.org/html/cond-mat/0502068).
That mathematical possibility does not prove that its unitary history and
charged observables are the native observer's physical ones.

### 5. Charged quantum response: test the observable algebra itself

The diagonal occupation algebra is fixed by the reference's modular flow.
The nontrivial modular clock frequencies occur in off-diagonal operators,
precisely the information absent from a population-only reduction.
Consequently a better population projection, even an exact one with memory,
does not supply the missing charged quantum field.

For the native C/E current instruments, the signed configuration difference
$\delta=x-y$ of a matrix unit $|x\rangle\langle y|$ is preserved by
each current branch. Common occupations can move, but the positive and
negative support of the relative-charge coherence does not. A classical
record-dependent choice of later seams preserves this property; retaining
memory does not by itself create charged transport. Diagonal modular
evolution adds phases and also preserves that support.
Here “support” means support on the native ports. This is not an exclusion
of propagation in every possible collective or transformed spacetime field
basis. Such an alternative needs its actual observable map and response.

There is a finite-domain issue too. On the single-carrier total-two space,
$P_2c_iP_2=0$ for a native occupation-changing annihilator. The 66-state
space cannot itself carry the full canonical algebra of twelve CAR modes.
This does not forbid encoded fields on multiple carriers, approximate
continuum fields, neutral composite operators, or a different coherent
construction. It does show why the separately constructed 90-mode chiral
CAR field must have an explicit physical realization map. It is not
already identified by the native occupation band.

#### Adjacent charge sectors and the actual repair rule

The [adjacent-sector calculation](../../model/common_source/hierarchy_diagnosis/charged_sector.py)
extends the existing endpoint primitive unchanged to binary loads zero
through three, with dimensions $1,12,66,220$. It checks all 7,920 retained
total-two C/E input addresses against the existing instrument. The adjacent
sectors suffice to verify
$P_2\{c_i,c_j^\dagger\}P_2=\delta_{ij}P_2$ and to evaluate the two-insertion
charged response about the actual total-two reference. The chosen ordering
in this CAR completion is explicit; this is a test of a possible native
field representation.

This enlargement does **not** restore charged propagation on native ports.
For the same placement coin, equal C records have the same current $j$.
Equal E records have the same initial endpoint contrast and therefore the
same current. Hence every surviving matrix unit obeys

```math
x'=x-j(e_a-e_b),\qquad y'=y-j(e_a-e_b),\qquad x'-y'=x-y.
```

The identity holds for arbitrary nonnegative integer loads, independently
of the total-two restriction. It composes along every finite recorded seam
history, including record-dependent seam choices. An operator carrying one
native-port charge remains in its original signed-difference block. Against
the diagonal source reference, the two-point charged response between
distinct ports is exactly zero. The executable one-proposal calculation
retains both accepted responses and the terminal timeout separately and
verifies every off-port entry as a rational zero. No finite timeout is
restarted or divided away.

Thus the next field construction must change a substantive ingredient:
coherently recombine different current records, or derive a collective or
encoded charged field with its own propagation map. The result does not
exclude those alternatives or propagation of neutral occupation densities.

#### A collective charged response with local observation records

The [collective response calculation](../../model/common_source/hierarchy_diagnosis/collective_response.py)
constructs an observable directly on the existing 66-state carrier. Write
$|a,j\rangle$ for the binary configuration occupied at ports $a,j$. With
protected port $a=0$ and reference port $b=5$, set

```math
X_j=|0,j\rangle\langle5,j|,\qquad j\notin\{0,5\},
\qquad Q_l=n_0n_l\quad(l\ne0).
```

The same source gives $\sum_{l\ne0}Q_l=n_0$ and
$[Q_l,X_j]=\delta_{lj}X_j$. Every completion seam $(u,v)$ preserves $n_0$
and satisfies

```math
Q'_u-Q_u=-n_0j_{\mathrm{native}},\qquad
Q'_v-Q_v=+n_0j_{\mathrm{native}}.
```

Thus the collective field has a local charge and conserved correlation
current. All 2,640 source/coin/seam input addresses are checked exactly.
This is a core-relative charge whose electric identification remains open.
Reading its current additionally requires access to the protected occupation;
that occupation is not inferred from the native current record alone.
Density and Ramsey quadrature readouts are alternative interventions.
The [native current experiment](../../model/common_source/full_repair_response/README.md#a-measured-native-current-response)
also derives the collective current covariance scale $1/36$, alongside the
raw current scale $19/150$, with the source instrument's actual first and
second moments. Its repair-count response and modular response use distinct
parameters; neither covariance scale is an electromagnetic coupling.
The complete signed records determine every nonzero temporal Fourier
coordinate of the occupation field over a finite native-event window,
including open paths; the zero-frequency occupation remains ambiguous.
The exact frequency Ward identity relates this reconstructed density to the
recorded current. In the fixed selector's stationary linear response, all
nonwhite current poles are positive real density-relaxation poles. The eleven
cycle-current directions instead have white two-point spectra and zero
conditional mean; one bounded cycle read has variance $19/1620$. This
operational current field has no oscillatory photon pole or derived Maxwell
normalization. These frequencies are in native event index, not calibrated
physical time.

The signed endpoints $0,5$ remain fixed, while their common occupation $j$
propagates. For the existing one-round completion selector, its accepted
ten-mode response is

```math
T=(q-k)I-\tfrac12L_{\mathrm{common}},\qquad
q=\frac{27}{32},\qquad k_C=\frac3{64},\quad k_E=\frac3{32}.
```

$L_{\mathrm{common}}$ is the graph Laplacian formed from the actual seam weights.
The generating response is $(I-zT)^{-1}$ on its convergence domain, with
all ten modes retained. Averaging signed local return measurements over the
ten separately prepared coordinates gives
$\mathrm{Tr}(I-zT)^{-1}/10=\sum_r[10(1-z\lambda_r)]^{-1}$.
Every transfer eigenvalue therefore has a nonzero residue in this declared
local measurement family, with multiplicities included. These positive real transfer modes describe decay
and diffusion in proposal count. The ideal reference's modular frequency is
$\log(25/8)$ for the nine outer common-occupation coordinates and
$\log(8/5)$ for the isolated $j=1$ coordinate. A modular advance remains a
separately implemented control with its own parameter.

The actual first incident event on cap $(0,1,5,7)$ supplies a stopping clock.
With two work proposals allowed, the complete recorded histories transport
the prepared relative phase from $j=2$ to $j=7$ along $2\to6\to7$ and
$2\to10\to7$. A local Ramsey detector has exact **unconditional signed
success response $1/2048$**, for both C and E. Success has probability
$135/512$, selector timeout $135/512$, and deadline timeout $121/256$.
No timeout is continued or success probability divided away. All 4,881 C
paths and 10,921 E paths retain their native work records, zero probabilities,
twenty Ramsey/location outputs and complementary detector output.
The selector's auxiliary transcripts are marginalized to the accepted seam
or timeout at this interface; these counts are complete work-path counts,
not an enumeration of every auxiliary selector transcript.
Reversing the input phase reverses the transported signal; the incoherent
mixture has zero signal. One allowed work proposal cannot reach that detector.
Preparation and coherent readout are declared apparatus, executed here as
exact CP semantics without a new pulse-level implementation claim.

The $j=7$ detector needs only ports $\{0,5,7\}$, with effects

```math
E_\pm=\tfrac12(|0,7\rangle\pm|5,7\rangle)
                    (\langle0,7|\pm\langle5,7|),\qquad
E_{\mathrm{rest}}=I-E_+-E_-.
```

The actual incident cap and a second declared access region $(0,5,7,10)$
overlap on those ports. Their local effects embed identically on the source,
preserve total occupation, and satisfy $E_vE_u=\delta_{uv}E_u$. Conditional
on a nonzero-probability outcome, the work state becomes
$E_u\rho E_u/\mathrm{Tr}(E_u\rho)$. Two same-event reads therefore have
joint record law $p(u,v)=\delta_{uv}p(u)$, including the rest and timeout
outputs. These are compatible observations after a concrete local state
update. The work state can continue repairing while the classical record
remains protected; a fresh measurement after intervening repair need not
return the same outcome.

The field statistics are also explicit: $X_j^2=0$ and
$\{X_j,X_j^\dagger\}=P_j$, a rank-two projector. Different $j$ act on
disjoint supports and commute. Their algebra is a direct sum of ten
$M_2(\mathbb C)$ blocks, without ten independent CAR modes or a derived
bosonic field algebra. The sum of the $X_j$ realizes one CAR mode on their
20-dimensional support with coordinate multiplicity. These are even
occupation-transfer operators.

Their accepted ten-mode transfer has one double eigenvalue and eight simple
eigenvalues under each C/E law, certified by its exact characteristic
polynomial. Any unbroken continuous $SU(3)$ action confined to these linear
response modes and commuting with that transfer must act trivially: every
invariant eigenspace has dimension at most two, while a nonzero map from the
simple eight-dimensional $\mathfrak{su}(3)$ into $\mathfrak u(2)$ would have
to be injective. This excludes a direct color-triplet reading of these ten
linear modes under the stated symmetry condition, rather than color in an
enlarged or differently encoded field.

A separate [same-carrier exterior test](../../model/common_source/hierarchy_diagnosis/native_exterior_attachment.py)
maps binary two-port configurations to a chosen exterior basis. The native
constant and rank-five bands give a rank-15 subspace split as $10+5$; the
protected core projects to a rank-two plane in the five-band, with exact
Gram determinant $1/6$. Under the additive one-body charge $q_{6Y}=-2$
on its three-dimensional complement and $q_{6Y}=3$ on that plane, the
ten-component charges match the decoded matter module. The five-component
charges are three copies of $-2$ and two of $3$, whereas the decoded
conjugate-five charges are three copies of $2$ and two of $-3$. A single
charge shift of the constant line would require
both $4$ and $-6$ and cannot repair the mismatch. Exterior orientation
signs are an additional representation choice; native unordered-pair repair
does not supply them.

The actual one-round mixed-cap population transfer fails to preserve this
particular rank-15 projector: its commutator has rank 30 and squared
Frobenius norm $961/8192$. That is a population-transfer diagnostic;
the actual quantum Kraus branches also leak from input $Q_{15}/15$, with
accepted-branch outside probabilities $379/2560$ (C) and $703/2560$ (E).
The complete C/E channel yields a stronger
work-state result. Its positive 20-seam support has the existing
protected-core fixed algebra
$\mathbb C I_{45}\oplus(M_2\otimes I_{10})\oplus\mathbb C I_1$;
the finite selector minorizes the ideal connected schedule by $3/4$ and
retains its $5/32$ timeout. Exact bistochastic Kraus counts give possible
fixed-projector ranks $0,1,10,11,20,21,45,46,55,56,65,66$.
Because an invariant work-subspace projector of a bistochastic channel is
fixed, **no noiseless rank-15 work-state code exists in this fixed carrier
under the current repair rule**. This does not exclude an operator-response
module, a record-corrected code, an ancilla, or another coherent action.
A declared signed exterior lift of the same port Laplacian preserves the
$10+5$ carrier exactly; the recorded native rule does not select that lift.

A conditional fluctuation construction gives a precise bosonic route from
these same operators. For $M$ independently prepared ideal reference copies,
set $F_j=(M\delta)^{-1/2}\sum_r X_j^{(r)}$ on the nine outer coordinates,
where $\delta=1/36-2/225=17/900$. Direct tensor calculations give

```math
\langle[F_j,F_k^\dagger]\rangle=\delta_{jk},\qquad
n_B=\langle F_j^\dagger F_j\rangle=\frac8{17},\qquad
\left\langle([F_j,F_j^\dagger]-I)^2\right\rangle
=\frac{29411}{289M},\qquad
\langle(F_j^\dagger)^2F_j^2\rangle
=2\left(1-\frac1M\right)\left(\frac8{17}\right)^2.
```

The [bounded IID ordered-Weyl argument](../../model/observer_time_fields/REFINEMENT.md#source-coded-quantum-fluctuations-and-an-ordered-detector-law)
applies to the actual nonuniform source weights: the ordered characteristic
functions converge to the thermal CCR state with occupation $8/17$ and
$e^\Omega=25/8$. Direct finite tensor exponentials and reversed ordering
check the noncommuting phase. This is convergence of the specified state
and observables, without a finite operator-norm CCR claim. Independent
reference copies, joint control and their preparation/error costs are
additional services; a growing-copy detector is not executed here.
With identical separate-copy evolution in the same declared parameter,
this normalization leaves every two-point transfer pole unchanged.

The hierarchy conclusion is consequently explicit. This linear collective
response contains only two positive modular frequencies,
$\log(8/5)$ and $\log(25/8)$, whose ratio is less than three. It cannot
supply three distinct charged-lepton masses under one common modular-clock
calibration. Even the additional assignment of attenuation rates to masses
fails the hierarchy test for each of these named conventions:

| Transfer whose positive eigenvalues are logged | C: certified rate-spread upper bound | E: certified rate-spread upper bound |
| --- | ---: | ---: |
| Accepted proposal, including survival loss | $1856/455$ | $4$ |
| Proposal conditioned on selector survival | $342/35$ | $99/16$ |
| Two-proposal first-incident success map | $3072/167$ | $16128/229$ |
| Same success map conditioned on its flag | $515/48$ | $635/24$ |

These exact bounds follow from $0<l\le\lambda\le u<1$ and
$1-u\le-\log\lambda\le(1-l)/l$. For the incident map
$J(I+R)$, positivity follows by similarity to
$(I+R)^{1/2}J(I+R)^{1/2}$. Every bound is below 71, whereas the measured
muon/electron ratio exceeds 206. A common scale, field normalization or
choice among these response modes cannot close that difference. The
conditional success map is a diagnostic of the same finite experiment;
it does not continue failed clocks. These are bounded exclusions of the
listed linear response and clock assignments, not of composite fields or
other continuing constructions.

An exact attachment test locates the spacetime gap. The connected outer
graph has one zero Laplacian mode and eight positive modes, while all nine
have the same modular frequency. If that frequency is identified with one
field's energy and a positive multiple of the Laplacian eigenvalue with its
momentum squared, the sharp free relation $E^2=m^2+c^2p^2$ fails for every
nonzero $c$. The actual native transfer also has no oscillation phases.
This rejects that specific free-pole dictionary. It does not identify graph
eigenvalues with measured continuum momenta, turn proposal count into proper
or Euclidean time, or exclude interacting, thermal or other continuum field
descriptions. The positive result is native collective charge transport and
an agreed local readout; the physical energy–momentum map, field statistics
and Standard Model mass spectrum remain unresolved.

#### A modular clock within one shell does not select charged thresholds

There is an independent ambiguity even if coherent dynamics is admitted.
On the full 4,096-state binary carrier, consider the faithful references

```math
\sigma_\eta(x)=\frac{w_{q(x)}\eta^{(N(x)-2)^2}}{Z_\eta},
\qquad \eta>0.
```

Here $q$ is the protected-core occupation and $N$ the total load. Every
$\sigma_\eta$, conditioned on $N=2$, is exactly the existing source
reference. Each completion seam preserves $q$ and $N$ and connects equal
diagonal weights, so the same completion repair preserves all these
references. Their modular generators differ, up to a scalar, by

```math
K_\eta-K_1=-\log\eta\,(N-2)^2.
```

All total-two histories and their within-shell modular clock frequencies
are identical. Yet adjacent-sector addition and removal frequencies each
shift by $-\log\eta$. A charge-conserving two-copy transition
$(2,2)\leftrightarrow(3,1)$ shifts by $-2\log\eta$; a common chemical
potential cancels from this comparison. Exact binomial sums verify two
distinct examples, $\eta=1$ and $\eta=1/2$, without empirical inputs.

These are explicit reference extensions, not selected physical vacua or a
proof that either extension has a local Lorentzian realization. They show
precisely what the existing clock data leave undetermined: the relative
coherent dynamics of sectors visited by a charged field. Selecting the
full reference and continuing action, or deriving an internal
charge-transfer experiment that fixes their relative response, is an
additional obligation. The remaining uncertainty cannot be resolved by
relabeling the already measured within-shell clock frequencies as masses.

The quantum calculation resolves all 4,356 matrix units into exact
signed-difference blocks. Six operator directions are fixed by completion.
The unmodified C/E maps have real positive spectra, without intrinsic
oscillation phases in their proposal parameter. The same reference's
modular frequency menu is exactly
$0,\pm a,\pm b,\pm\Omega$, whose positive-frequency spread is about
3.4243. Composing one modular advance with a repair tick is an explicit
additional prescription; it gives multipliers

```math
z=\lambda e^{i\omega s},\qquad
\text{attenuation per proposal}=-\log|z|,\quad
\text{phase per proposal}=\omega s\pmod{2\pi}.
```

The frequency per modular parameter is $\omega$ after using the specified
advance $s$ and fixing the phase-alias branch. A proposal count and a
modular parameter are different time labels.

Changing C to E or changing to the incident clock changes the attenuation
in the tested probe while preserving its modular frequency. The two
quantities are therefore operationally different. Logging the accepted
map also adds the common selector-survival term $-\log q$; that term
must not be mistaken for a new light species.

Complete two-proposal quantum histories retain 12,962 C paths and 29,042 E
paths, including zero outcomes, selector timeouts and both terminal pointer
results. No physical charged-particle interpretation follows merely from
resolving these quantum records.

There is also an exact finite exclusion for the direct operator-log
assignment. In the six protected operator sectors, rational positive-pivot
LDL factorizations prove a nonfixed conditional-transfer gap of at least
$1/100$. Exact row-sum loss bounds cover the other 4,270 matrix units.
Together with the exact positive transfer floor, every nonzero rate obeys

```math
\frac1{100}<-\log\lambda<2,
\qquad
\frac{\max(-\log\lambda)}{\min_{\lambda<1}(-\log\lambda)}<200.
```

This is a finite rational certificate for both C/E laws, not a fitted
eigenspectrum or an imported exclusion-process theorem. The measured
muon/electron ratio already exceeds 200. Thus no common-scale assignment of
these nonzero proposal-log rates to the charged leptons works. The conclusion
is restricted to that assignment; the rates have not been identified with
physical masses. It does not alter the earlier $1/600$ source gap bound.

For the geometric-wait incident clock, all 4,356 quantum modes are also
retained using blocks of dimension at most 66. The stopped map is similar
to a positive operator. Its resolved finite positive logarithmic spreads
are approximately **13.65 for C** and **22.29 for E**, with six fixed modes
in each. E also has 264 zero transfer modes at the numerical tolerance;
260 are certified by identically zero incident blocks. They are retained
as erased modes with no finite logarithm, not discarded or relabeled as
particles. These floating clock diagnostics have no outward spectral
certificate and do not inherit the proposal-log bound above.

### 6. An independently specified coherent continuation

We also test the existing coherent-record continuation with the ideal
mixed-cap weights. Its isometry $A_c$ coherently prepares all seam/coin
events, and $S_c$ reverses the native branch while retaining its record.
The declared unitary is

```math
W=S_c(2A_cA_c^\dagger-I).
```

This is a weighted extension of the existing executable continuation,
using the standard spectral relation underlying
[Szegedy's quantum walks](https://arxiv.org/abs/quant-ph/0401053). It is a
different coherent operational contract from the finite selector with
committed timeouts; no failed history is restarted. Event records stay
coherent, and the reflection phase is a declared choice.

All 2,640 native arcs and the 1,584 coin-code padding amplitudes are
retained. The coin reflection is executed with the existing pulse compiler
on all 64 input columns. Native branch reversal is executed as the exact
source-defined permutation; a complete pulse compilation of that entire
permutation is not claimed. Two untouched reference columns and both
terminal Ramsey outputs are retained for the continuing walk. Intermediate
rows describe alternative terminal horizons, not a sequence of committed
measurements followed by undisturbed evolution. The controlled-walk Ramsey
instrument is declared mathematically; its entire pulse word is not compiled.

On its cyclic subspace, the principal phases are
$\pm\arccos\lambda(T_{66})$; the complementary zero and $\pi$ modes
are explicitly counted. The existing ideal completion gap bound
$1/600$ gives a nonzero cyclic phase above $1/18$. Since every
principal absolute phase is at most $\pi<22/7$, the spread of all
nonzero elementary principal phases is less than **$396/7<57$**.
Thus identifying those frequencies directly with the three charged-lepton
masses cannot reproduce their hierarchy.

This bound is about principal elementary frequencies relative to zero.
It does not classify arbitrary transition differences, composite operators
or other logarithm branches. Nor does it supply a chiral species assignment,
a 3+1D dispersion or a calibrated physical clock for the walk. It tests a
concrete coherent alternative; it does not exclude all coherent OPH laws.

### 7. Calibrate the mass observable on an actual charged field

An additional test uses the declared chiral candidate as a control whose
charged field and dynamics are fully specified. Charge/weak-Casimir
projectors select its charged-lepton field responses; nonzero spectral
residues and three unequal momentum probes along different axes test
whether the same response has a common rest-energy parameter.

Using its source clock and $c=\Omega/2$, the control recovers

```math
\left(\frac{E_i(p)}{\Omega}\right)^2-\frac{|p|^2}{4}
=\left(\frac{m_i}{\Omega}\right)^2.
```

The independent finite-matrix residual is below $7\times10^{-16}$.
A separate eight-CAR-mode source calculation at paired momenta uses the
existing compiled local words. Refining from four to eight split steps
reduces the rest-energy-square discrepancy from roughly
$3.22\times10^{-6}$ to $8.06\times10^{-7}$. These are numerical
convergence diagnostics. This operator-response calculation is not a newly
compiled, committed charged-reference detector experiment on the native
66-state repair history.

The most useful control compares **different observables of the same
Hamiltonian and ground state**. For the existing rest pair,

```math
C_{\mathrm{charged}}^{E}(\tau)=\tfrac12 e^{-m\tau},\qquad
C_{\mathrm{neutral,connected}}^{E}(\tau)=\tfrac14 e^{-2m\tau}.
```

Its Lorentzian charged response has signed poles at $\pm m$, while the
empty-pair occupation quench oscillates at $2m$. Thus an exponential
decay rate can equal the mass parameter when it is the appropriate charged
Euclidean correlator of that same Hamiltonian. Even in this elementary
control, a neutral decay rate gives a different excitation energy. Generic
repair relaxation has no automatic identification with either correlator.

The operational target is therefore a spectral pole or appropriate
threshold in the identified physical channel, with its charge/statistics,
dispersion, state and clock established. The current finite quadratic
control has isolated poles. General interacting physics requires more
care: charged infrared sectors can lack a sharp mass eigenstate under
Gauss-law assumptions, as shown by
[Buchholz](https://www.sciencedirect.com/science/article/pii/037026938691110X).
Unstable particles require the appropriate pole convention; see the
[PDG electroweak review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf).
Quark mass parameters require a stated renormalization scheme and scale;
see the [PDG quark-mass review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-masses.pdf).
Those qualifications are not mechanisms for changing the failed trial's
coefficients without deriving additional dynamics.

### 8. Current conclusions

The calculation retains the complete complementary response. The incident-clock
transformation is explicit. Both the static mass-seed test and the charged
operator test expose structural limits that were hidden by assigning a
single mean matrix to Yukawa couplings. The current tested routes do not
recover the measured mass ratios.

The [stationary schedule classification](../../model/common_source/full_repair_response/stationary_selection.py)
tests every nonnegative, state-independent seam weight while preserving the
same faithful core-biased native reference. Stationarity forces the weights
of all eight core-to-outer seams to zero; the maximal allowed schedule has
22 seams, including the core-internal seam. Every
admitted recorded history conserves core occupation. Its static Schur family
seed has rank two, and the same declared $M_{15}$ Higgs route then has
14 massive and seven massless charged Dirac fields, counting color. Positive
kinetic normalization cannot remove the zero. Giving one crossing seam
weight $1/64$ restores rank three but changes the reference by the exact
Dirichlet defect $157/5{,}760{,}000$; uniform positive use of all 30 seams
instead selects the uniform diagonal reference with scalar family feature
$I_3/6$. Thus passive schedule changes within this fixed-reference route
cannot supply three nonzero charged-family masses. The result does not
exclude active recovery, other references or other mass laws.

The collective correlation field supplies an explicit positive result: local
charge continuity, transported phase-sensitive response, and compatible
overlapping observer records through the same native history. Independent
reference copies admit a thermal bosonic fluctuation limit under their
declared preparation contract. Their two-point poles retain the native
frequency and attenuation bounds, which do not recover the lepton hierarchy.
The signed records determine the nonzero finite-frequency density response,
while the measured linear current has no photon pole. The same 66-state
carrier admits a specified $10+5$ exterior subspace, but its natural
five-component charge is not the decoded conjugate-five and the complete
native channel has no invariant rank-15 work-state code. The recorded
incident response cannot be intertwined with the unmeasured closed
Weyl–Higgs evolution by a nonzero linear map. A detector/environment
channel or a different encoded history requires its own source derivation.
The physical field and energy–momentum attachment remains open. The positive
detector and continuity results resolve specific finite observation questions;
they do not establish the full physical detector map. The hierarchy failures
exclude the stated action prescriptions and clock/readout assignments. They
do not prove that every possible OPH construction, interacting effective
description or collective field fails.

The protected seven-momentum construction implements a finite encoded
operator map with active record correction. It removes the rest-only
execution restriction for the supplied trial action. It does not select that
action from the native records, alter the failed hierarchy, or establish a
physical spatial/clock dictionary. Its use of many corrected carriers is
outside the single uncorrected rank-15 obstruction. The coherent native
entropy theorem separately supplies a constrained endpoint Hessian and
recovery algebra; no physical area, stress or joint matter/gravity variation
is inferred by placing these results together.

The next physical construction must supply an operator map from a
continuing repair history to the charged left/right field response. It
must preserve charge, CAR or its justified effective limit, propagation,
current vertices, the state and clock, and the complementary information
needed for subsequent measurements. The same map must produce the mass
and kinetic terms; their numerical resemblance to a desired matrix is
insufficient.

This is the point at which to examine an unproved attachment, incomplete
repair context, larger carrier or genuinely different coherent realization.
It is not a reason to tune the successful-decoder branch, count a timeout
hazard as a light particle, or fit an arbitrary function of the family
matrices and call its coefficients derived. A successful construction
would have to pass the joint measured ratios after those choices are fixed.

#### Prioritized open tests

The immediate research target is an effective particle description of the
**same jointly recorded observations**. Settled records need not imply a
static source. Different observers must recover compatible response functions
after the measurement update, including their future conditional predictions.
The following routes remain hypotheses to test, rather than assumed solutions.

| Priority and route | Concrete calculation | Criterion for progress |
| --- | --- | --- |
| 1. Attach the response to spacetime and a clock | Specify a charged field, spatial translations, observer-clock readings and detector/environment channel from one continuing history. Measure several spatial modes and orientations with every terminal flag and waiting record. Test a sharp pole where appropriate, or state an interacting/thermal spectral alternative. | One fixed map preserves the joint detector/intervention law and gives compatible energy–momentum parameters across observers and probes. It must address the exact incident contraction and failed direct free-dispersion map. |
| 2. Determine the continuing coherent implementation | Compare explicit admissible source implementations, including coherent versus committed records and the actual nonzero-frequency current/density readouts. Fix each before examining particle targets, and measure multi-time phase-sensitive responses with its own clock. | The complete source program determines the observable evolution and errors. One incident count cannot replace the waiting record, and committed C/E probabilities alone do not select an unread coherent continuation. |
| 3. Derive interactions and family structure together | Use the continuing repair context and complementary memory to calculate charged response, kinetic normalization and vertices from one action. Show how the decoded conjugate-five and color currents act on the actual field carrier, and test the three-family symmetry breaker. | Distinct physical gaps and residues survive the same charge/statistics and spacetime tests; family coefficients and scale come from that specified construction. The ten-mode direct color action and natural $10+5$ charge action tested here are insufficient. |
| 4. Test controlled enlargement or refinement | For an operator-response code, active record correction, larger carrier, collective sector or copy limit, state preparation, locality, interactions and observer access. Bound complete-record errors as the construction grows. | A common limiting response is obtained with controlled clock and state errors. It must bypass the fixed-work-carrier rank-15 obstruction by an explicit mechanism. Independent copies alone preserve the existing two-point poles. |
| 5. Join the matter and entropy readouts | Extend the [native record entropy result](../../model/native_equilibration/README.md) to the actual coherent matter protocol, then derive geometric area/volume and stress readouts independently on that same history. | The joint construction meets the entropy and geometry contracts without defining curvature or an area term solely to enforce Einstein's equation. The coherent native entropy theorem has a fixed-carrier state-variation domain; it does not close this geometric attachment. |

For any route that passes these checks, freeze its architecture, state,
controls and physical readout before comparing mass ratios with observations.
Use a common clock and justified uncertainties; specify the renormalization
scheme and scale where the physical mass definition requires them. The
comparison must include the simultaneous hierarchy and relevant interaction
tests. Numerical replay tolerances are not physical theory-error budgets.
An inconclusive attachment or a failed specified construction should remain
explicit, with the other admissible routes open.

### 9. Reproduction

The retained machine-readable report and complete coherent-work arrays are
[results.json](../../model/common_source/full_repair_response/results.json) and
[results.npz](../../model/common_source/full_repair_response/results.npz).
The receipt binds the implementation, tests, these notes, pinned dependency
file and supporting evidence by SHA-256. Its verification command recomputes the
report and arrays; an artifact digest alone is not the replay check.
Exact identities and discrete claims are checked exactly. Floating diagnostics
use the stated replay tolerances, which are not physical error certificates.

From the repository root, using Python 3.12 and the root `requirements.txt`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/mass_response model/common_source/particle_identification model/common_source/chiral_candidate model/common_source/hierarchy_diagnosis model/common_source/full_repair_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.full_repair_response.run --verify model/common_source/full_repair_response/results.json
.venv/bin/python -m model.common_source.full_repair_response.stationary_selection --check model/common_source/full_repair_response/stationary_selection_results.json
```

For a new independent output, use `--output` with a fresh JSON pathname;
the producer also writes its adjacent NPZ and refuses to replace either.
The [package README](../../model/common_source/full_repair_response/README.md)
links the supporting packages' separate replay commands. Their source-binding
receipts point to this unified report; the calculated trial predictions are
unchanged. The present research does not promote conditional
field constructions into a released proof of the Standard Model spectrum.

A continuation should begin with an explicit proposed field operator,
carrier/domain, state and continuing law, then specify how its charged
two-point response is read with the observer's clock. Test that proposed
map against the native support and fixed-sector obstructions above before
introducing a new family mass matrix. If an enlarged carrier or coherent
record protocol is needed, declare that change and reproduce both its
observable response and the existing successful construction checks.
Only then compare rest-energy ratios and interaction currents with the
physical species. The unresolved task is **the native-to-physical charged-field
map and the coherent dynamics of the sectors it connects.**
