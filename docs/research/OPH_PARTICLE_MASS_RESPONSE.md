# Particle construction: results, assumptions and research tasks

The executable construction realizes calibrated Standard Model charged-matter
structure and permits observer-record recovery of its supplied mass parameters.
**No measured particle spectrum is derived from the repair law.** The uncalibrated
mass prescriptions fail their stated comparisons. A physical chiral interacting
continuum, its vacuum and its mass/clock calibration are unproved.

The fixed twelve-port carrier is used throughout. The research question is which
repair dynamics and observer implementation produce a physical particle
response on that carrier. Internal representation labels, persistent levels,
repair decay rates and physical particles are different objects until an
explicit observation map connects them.

| Result | What is established | Additional attachment |
| --- | --- | --- |
| Gauge and matter structure | A specified closed twelve-current law, a spin/exterior lift, the fifteen internal states of one Standard Model generation, and three invariant charged-fermion Higgs lines | Complete physical current selection, physical spin/statistics, and an interacting field |
| Family structure | A source-executed protected rank-three band and several explicit symmetry-breaking memories | Identification with three physical generations and selection of their couplings |
| Calibrated realization | Supplied running Yukawas and mixing parameters enter the finite action and are recovered from simulated records | Derivation of those inputs and a physical pole/threshold comparison |
| Uncalibrated coefficients | Complete record laws determine the specified feedback coefficients; blind inverses recover them | The tested coefficients do not reproduce the charged spectrum |
| Propagation and measurement | Finite CAR fields, charge/current identities, complete detector controls, and restricted continuum comparisons | One same-process physical chiral field, controlled interactions and physical units |

Read the [particle paper](../../paper/deriving_the_particle_zoo_from_observer_consistency.tex)
for the mathematical argument. The [model construction criteria](../../model/CONSTRUCTION.md),
[checked premise map](../../model/PREMISES.md) and
[particle provenance ledger](../../code/particles/PARTICLE_PROVENANCE_LEDGER.md)
record the broader acceptance conditions. This page gives the working equations,
decisive exclusions, code entry points and tests needed to continue the research.

## Candidate construction

### Carrier, currents and the finite Standard Model interface

Native port permutations and record-correctable codes supply the geometric
$3\oplus3'$ object. The relevant response constructions have different domains:

| Construction | Exact result | Boundary to preserve |
| --- | --- | --- |
| [Native exchange response](../../model/common_source/native_gauge_object.py) | The tested preserving response algebra has dimension 17 | Twelve ports do not imply twelve independent, closed currents |
| [Specified current law](../../model/common_source/native_sm_current.py) | A closed $\mathfrak u(3)\oplus\mathfrak{so}(3)$ tangent with native rotation covariance and observable relative center | The control menu is specified; complete physical-intervention selection is unproved |
| [Coherent operator filter](../../model/common_source/native_filtered_current_control.py) | An exact invertible twelve-control map into that same algebra, with conditional spin/exterior/Higgs joins | This is an added Hamiltonian control, with distinct port drives; its literal native pulse history is not compiled |
| [Full-current controller](../../model/common_source/native_response_controller.py) | All eighteen exchange-plus-class-sum directions enter the disturbance calculation; a bounded-drive joint limit approaches twelve observable currents | Finite filters are invertible; fixed-bias long-time response has six observable currents. Bias, baths and feedback are supplied |

These results instantiate restricted finite contracts. They do not prove that
unchanged C/E repair autonomously selects the complete physical gauge current.
The [finite theorem comparison map](../../model/common_source/particle_identification/SM_MATCH.md)
and [same-current matter implementation](../../model/common_source/native_matter_attachment.py)
identify the actual maps rather than relying on an abstract Lie-algebra match.

The declared spin lift sends a six-state generator
$D=\operatorname{diag}(A,\widehat w)$ to

```math
L(D)=\operatorname{diag}\left(A,-\frac{i}{2}w\cdot\sigma-
                              \frac{\operatorname{tr}A}{2}I_2\right).
```

The nonvacuum even exterior module
$M=\Lambda^2\mathbb C^5\oplus\Lambda^4\mathbb C^5$ has dimension fifteen.
Source-computed charges and Casimirs identify
$Q,u^c,d^c,L,e^c$ and the compatible central quotient
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$.
The spin cover, exterior-sector choice and physical CAR field are separate
attachments. Even auxiliary exterior parity does not establish physical
fermion statistics. No topological-defect identification is proved.

Both components of each invariant Higgs tensor are checked against all twelve
currents. The allowed lines are $Qu^cH$, $Qd^cH^*$ and $Le^cH^*$.
With three family copies, they admit three arbitrary complex $3\times3$ Yukawa
matrices. The exact gauge/Higgs identities do not determine those 27 complex
entries. The finite family-band theorem selects a rank-three band under its
stated selection rules; it does not select the observed masses. An unbroken
irreducible family action with a family-neutral Higgs instead gives scalar
intertwiners and degenerate families.

### The mass observable and the field convention

For a neutral Higgs component of norm $v$, canonical kinetics give tree-level
action masses $m_{s,f}=v\,\sigma_f(Y_s)$. Here $v=v_{\rm SM}/\sqrt2$ in the
conventional electroweak normalization. For positive kinetic forms, the Dirac
convention gives

```math
M_s=v\,Z_{L,s}^{-1/2}Y_s Z_{R,s}^{-1/2},\qquad
Y_s^\dagger Z_{L,s}^{-1}Y_s u
   =\frac{m^2}{v^2}Z_{R,s}u.
```

Up and down share the left quark metric $Z_Q$. The physical kinetic forms are
not emitted by the representation theorem or by a GNS inner product alone.
Transforming both the mass and kinetic forms consistently preserves the gaps;
a change of coordinates cannot repair a failed spectrum. The
[full-doublet response specification](../../code/particles/flavor/exterior_higgs_response_spec.json)
states the Weyl/Dirac convention and kinetic-whitening contract.

In the all-left convention, orient each internal Higgs tensor from a weak
doublet to a singlet. Its directed family coefficient is $\overline Y$;
the reverse block is its **transpose**, not its adjoint. The resulting symmetric
internal kernel $K_{45}$ has creation pairing
$\epsilon_{\rm spin}\otimes\overline K_{45}$. The
[field handoff](../../model/common_source/native_record_yukawa_field.py)
checks all 180 Nambu branches of the supplied Weyl symbol, including the three
minimal-model neutrino zero directions. Its spatial principal symbol and unit
Higgs background are declared. The
[kinetic guide](../../model/common_source/chiral_candidate/KINETIC.md)
separates a finite momentum symbol from a local chiral regulator.

An operational mass requires a charged response read with the observer's clock:
a stable pole and dispersion where applicable, or a specified unstable,
confined, thermal or threshold alternative. Masses, kinetic residues and
vertices must come from the same process. Dimensionless coherent steps,
modular parameters, repair-attempt counts and laboratory seconds are distinct.

<a id="56-calibrated-action-and-matching-control"></a>
### Calibrated control

The [calibration packet](../../model/common_source/calibrated_matter_parameters.json)
provides nine running charged Yukawas, mixing parameters and the Higgs scale in
stated schemes. The [action](../../model/common_source/calibrated_matter.py)
checks the internal mass/current identities. External
[SMDR matching](../../model/common_source/calibrated_matching.py) checks the
calibration; its outputs are external Standard Model calculations.

Separate [source-record experiments](../../model/common_source/calibrated_spectrum_observer.py)
and [moving-field readout](../../model/common_source/calibrated_transport.py)
recover the supplied running parameters with fixed probe design and a blind
inverse. The Hamiltonian retains the empirical inputs. Decoder blindness does
not remove that ancestry. A running action parameter is not automatically a
laboratory pole mass, and the selected marginal input uncertainties do not
constitute a joint theory-error model.

### Uncalibrated record-memory protocol

The [complete family protocol](../../model/common_source/native_family_protocol.py)
is the compact executable entry point for the filtered current, family memory,
current gate and field handoff. Its [receipt](../../model/common_source/native_family_protocol_results.json)
and adjacent NPZ retain all geometric and duty cases, including zero directions.
Its producer loads no empirical mass packet. Runtime separation does not certify
a prospective prediction or target-independent design history.

The [twelve-state memory](../../model/common_source/native_family_record_memory.py)
retains all 128 native C/E labels and all twelve old-memory inputs: 1,536 final
instrument branches. The advance event is $C=(0,1)$; the second event is
$E=(0,4)$ or $(0,5)$, with probabilities $p=1/12$, $q=1/6$ and joint $1/72$.
Specified seam half-turns update the memory, other outcomes and proposed core
destinations reset it to port 2. All labels, including zero-probability labels,
remain in the calculation.

For the reset-2 law, put $a=p(1-q)$ and $b=pq$. Stationary axis weights are
proportional to $1,a,ab/(1-a)$. Averaging the native rank-one port projectors
produces a positive trace-one coefficient

```math
Y=\sum_k\pi_kP_k,\qquad
\operatorname{eig}(Y)\ \propto\ 1:195.11703757:3579.82624571.
```

These are dimensionless coefficients. The light-axis history C then EC has two
memory updates and three rare instrument events; its weight scales as $p^3$
only when $q$ scales with $p$. An event probability, an amplitude and a tick
count cannot be interchanged.

Fixed-Pauli records at coherent steps $1/8$ and $1/4$ cancel the rank-one
finite-pulse bias exactly. Eight burn-in updates give a conditional operator
allowance of about $4.06\times10^{-10}$, with finite mixing and a stated
binary64 contrast allowance accounted for separately. This is a complete
logical feedback experiment, not a literal compiled native pulse history or
a finite-shot/physical error certificate.

The [current gate](../../model/common_source/native_gauge_family_records.py)
uses the trace-orthonormal six-state current frame lifted to matter. For
$C=\sum_aF_a^\dagger F_a$, its complete instrument is

```math
M_0=\cos(t\sqrt C),\qquad
M_a=iF_a\frac{\sin(t\sqrt C)}{\sqrt C},\quad a=1,\ldots,12.
```

The coarse jump/no-jump instrument is gauge covariant. Individual adjoint-basis
jump labels are not gauge-invariant species labels. Both Higgs components
identify invariant left/right supports, giving the same gate law
$p_s=\sin^2(t\sqrt{C_{L,s}})/12$,
$q_s=\sin^2(t\sqrt{C_{R,s}})/6$ in each sector.

| Channel matched after construction | $(C_L,C_R)$ | Coefficient ratios at the declared $t=\pi/5$ |
| --- | --- | --- |
| Up | $(25/8,4)$ | $1:271.446:6071.036$ |
| Down | $(25/8,3)$ | $1:311.762:6822.880$ |
| Charged lepton | $(9/8,3)$ | $1:675.782:30775.333$ |

All $t\in\{\pi/5,\pi/3,\pi/2\}$ cases use the same reset and family triple.
None reproduces the charged hierarchy; the exact rank-two case is retained.
The two-leg interaction, duration and feedback are declared architecture.

The separate static mean-feedback service $\exp(-ihX\otimes Y)$ has a complete
fixed-coordinate inverse using X/Y contrasts and survival outcomes. It uses
$h=1/4$, with numerical admission restricted to $1/16\le h\le1/2$; complex
probability payloads are rejected. This static experiment is **not identified
with the correlated finite-memory history**. Such an identification needs its
own averaging and source-implementation bounds.

## Hierarchy diagnosis

The following exclusions constrain specified mechanisms. They do not exclude
other repair laws or establish a programme-wide no-go.

### The normalized-coefficient obstruction

For any positive $3\times3$ projector average with trace one,

```math
\frac13\le\lambda_{\max}(Y_s)\le1.
```

The native neutral Higgs tensors have unit nonzero directed singular values.
Color supplies copies, not an extra sector strength. With common canonical
kinetics and one Higgs norm, the three family action masses therefore sum to
$v$ in every sector, counting each family once; sector maxima differ by at
most a factor three. Changing the event rates alone cannot remove this
restriction. The bound is for this mass map, not for arbitrary interacting
poles or different kinetic laws.

Unconditional pulse duties retain the identity on pulse-off branches. For
$D\in\{1,C,E,CE\}$ their first-order coefficients are

```math
Z_D=\sum_{m,r}\pi_mw_r\chi_D(r)P_{\mathrm{new}(m,r)},\qquad
\operatorname{tr}Z_D\in\{1,p_s,q_s,p_sq_s\}.
```

Event and destination are correlated: $Z_C=Y-(1-p)P_{\rm reset}$, generally
not $pY$. The four default-protocol sector-max spans are approximately
$1.0233,2.0572,1.1528,2.3695$. E and CE duties have exact rank two. No
pulsed-branch normalization or zero-mode removal is permitted.

### The static-family-feature obstruction

The [quadratic Weyl–Higgs candidate](../../model/common_source/chiral_candidate/README.md)
uses $Y_u=aF$, $Y_d=bF$, $Y_e=JF$, with
$a=\log(25/8)$, $b=\log(8/5)$, $J=a-b$ and $v=a+b$.
The source feature has eigenvalues

```math
f_1=\frac{61}{450},\qquad
f_{2,3}=\frac{41}{225}\mp\frac{7\sqrt5}{750}.
```

The common family ratios are $1:1.1903:1.4982$. The record inverse recovers
those same values. For every complex affine alternative $Y=uI+wF$,
$f_3-f_1=\phi^2(f_2-f_1)$ implies
$z_3-\phi^2z_2+\phi z_1=0$. The triangle inequality, allowing every assignment
of the largest modulus, gives

```math
m_{\max}\le\phi^2m_{\rm middle}+\phi m_{\min}.
```

The broad charged-lepton comparison windows require the left ratio to be at
least 3470, while this bound gives at most $207\phi^2+\phi<544$.
Thus the entire complex affine class is excluded under that identification.
The [hierarchy diagnosis](../../model/common_source/hierarchy_diagnosis/README.md)
retains the proof, nonlinear controls and normalization perturbation bounds.

For the fixed quadratic Nambu law,
$\{\Xi(t),\Xi^\dagger(0)\}=e^{-i\mathcal Bt}$ fixes the linear-field poles
independently of the reference state. A common clock/Higgs/unit rescaling
preserves mass ratios. Pair occupation oscillates at twice the odd-CAR gap;
that factor is included in the inverse. There is no missing interacting
self-energy inside this specified quadratic Hamiltonian.

### Other tested mechanisms

| Mechanism and primary implementation | Decisive result | Exclusion domain |
| --- | --- | --- |
| [Fixed corrected cycle](../../model/common_source/parameter_selection.py) and [clock control](../../model/common_source/repair_continuation.py) | Logical transfer is identity; local loads/costs compress to scalars | That cycle does not select an unwrapped Hamiltonian or family masses |
| [Five-memory feedback](../../model/common_source/native_repair_mass.py) and [native implementation](../../model/common_source/native_repair_mass_feedback.py) | Unequal coefficients, about $1:30.84:452.62$; the stationary trigger family fails an exact invariant test | Specified projector feedback and physical mass attachment |
| [Two-rotation memory](../../model/common_source/native_mass_record_cycles.py) | Ratios about $1:382.23:5909.43$ fail both broad lepton windows | Specified event and reset law |
| [Loop phases](../../model/common_source/native_mass_holonomy.py) | The tested positive finite-relation candidates fail the hierarchy | Supplied balance and squared-eigenvalue mass relation |
| [Auxiliary dressing](../../model/common_source/native_mass_dressing.py) | The two smallest subthreshold lepton poles have ratio below 32 for every admitted nonnegative dressing strength | The larger declared operator, its retained extra poles and this pole identification |
| [Selector bath](../../model/common_source/native_mass_bath.py) and [family potential](../../model/common_source/native_family_potential.py) | Five low charged responses in the complete bath; the simplest real curvature coupling retains a doublet | Neither dropping modes nor changing their labels yields three observed families |
| [Uniform-reference full entropy](../../code/particles/runs/flavor/full_entropy_w5_selector.json) | At fixed pair norm, every minimizing quadrupole has a repeated eigenvalue, including the closed boundary | Changing the amplitude cannot give three distinct masses through an eigenvalue-only readout; other references, constraints and responses are outside this exclusion |
| [Exchange field](../../model/common_source/native_fields.py) and [staggered extension](../../model/common_source/native_staggered_field.py) | Gap ratios $1:3:12$; free cubic field has two Dirac tastes per sector and both chiralities | Specified actions; no interacting chiral Standard Model limit |
| [Interacting charged response](../../model/common_source/native_charged_spectrum.py) | Dimer poles are action-derived and volume controlled, but fail the lepton ratios | The finite/interacting model does not supply a physical bulk mass |

## Full repair response

The [complete-response package](../../model/common_source/full_repair_response/README.md)
retains waiting, complementary modes and termination flags. Its independent
[stationary-selection calculation](../../model/common_source/full_repair_response/stationary_selection.py)
tests all nonnegative state-independent seam weights at the same faithful
core-biased reference. Stationarity forces all eight crossing seam weights to
zero. The maximal admitted schedule has 22 seams and conserves core occupation;
its Schur mass seed has rank two. Positive kinetic whitening preserves that
zero. A crossing seam restores rank only by changing the reference; uniform
use of all thirty seams selects a scalar family feature.

The same finite carrier supports collective charged continuity and compatible
local detector records. These are even matrix-unit transfers with a diffusive
response, not an identified physical chiral fermion field. The complete
uncorrected recorded channel has no invariant rank-fifteen work code, and its
incident contraction has no nonzero direct intertwiner to the unmeasured
closed Weyl evolution in the tested class. Active correction, collective
operator sectors and other states have different domains; the obstruction
must not be extended to them without proof.

The [protected field source](../../model/common_source/chiral_candidate/SOURCE.md)
executes the declared CAR action at rest and at three opposite momentum pairs.
A different two-hop local completion agrees on that restricted set and differs
on a second shell. Seven-momentum agreement does not determine a unique local
field law. The [same-current transport](../../model/common_source/native_gauge_matter_transport.py)
and [continuum comparison](../../model/common_source/native_gauge_continuum.py)
retain both tastes on a chosen background; they do not derive dynamical chiral
gauge fields.

<a id="515-zero-temperature-modular-response-of-the-same-free-field"></a>
The [free half-space modular calculation](../../model/native_equilibration/wedge_modular.py)
connects the specified staggered vacuum and boost-weighted energy, with a
controlled low-momentum normalization. Its bare planar entropy depends on the
regulator, mass and tastes. The [reciprocal Gibbs calculation](../../model/native_equilibration/reciprocal_vacuum.py)
has a separate finite thermal domain. Neither supplies a renormalized massive
small-ball entropy, independent metric stationarity, Newton coefficient or
interacting Standard Model vacuum.

## Potential sources of error

These checks distinguish a genuine calculation error from failure of a chosen
law. A proposed fix must identify the violated equation or hypothesis.

| Risk | Required check |
| --- | --- |
| Reusing an abstract theorem on another object | Exhibit the actual state sector, generator/control map, representation and both Higgs components; retain the control-menu and spin-cover premises |
| Cooling confused with a Hamiltonian filter | Compute the channel and dual. For rank-three $P$, $\Phi(\rho)=P\rho P+\operatorname{tr}((I-P)\rho)P/3$ has $\Phi^*(B)=PBP+\operatorname{tr}(PB)(I-P)/3$, not $B\mapsto PBP$. Its classical Kraus branches lose relative triplet phase |
| Missing records or conditioning | Retain raw, reset, complement, timeout and zero outcomes. Compare full laws to source replay; normalization or a matching digest alone is insufficient |
| Averaging substituted for execution | Separate correlated finite memory, stationary first-order coefficients and the static mean-feedback service; charge finite mixing, averaging and per-step implementation errors |
| Wrong kinetic or pairing convention | Check all-left conjugation/transposition, physical CAR ordering, both chiralities/tastes and full Nambu multiplicities. Obtain mass and kinetic terms from the same inverse response |
| A decoder artifact mistaken for a mass failure | Test factor-of-two conventions, pole residue, detector overlap, phase range and the full operator error. The static inverse rejects unresolved angles and complex probabilities |
| Phase assigned from an unrelated reversal | $C_3$ orientation reversal exchanges $\delta$ and $-\delta$ without identifying physical antiunitary time reversal. Classical cycle affinity is not a coherent $U(1)$ phase; detailed balance of one kernel does not classify every controlled memory extension |
| Calibration hidden by a blind inverse | Audit the producer and its design ancestry. Empirical Yukawas, Higgs scales, units, target-selected resets and fitted transport maps remain inputs even when the decoder sees records only |
| Incompatible physical comparisons | State mass definition, common scheme/scale, matching, thresholds, confinement/width treatment and joint uncertainty. Replay tolerances are not theory errors; a running chart cannot be assigned a pole-mass pull |
| Independent pieces promoted to one theory | Require one preparation, dynamics, clock and observation regime to satisfy the combined contract. A field symbol, an entropy identity and a gauge table do not prove their physical composition |

## Open challenges and research tasks

The next decisive object is one continuing source process with a **charged
chiral inverse propagator**. It must determine the kinetic derivatives/residues,
Higgs response and mass terms together. Selecting a different Gram matrix or
fitting a function of a family feature does not supply that process.

| Priority | Concrete calculation | Acceptance condition |
| --- | --- | --- |
| 1. Select the physical response domain | Specify carrier/code, actual C/E/coherent implementation, state, scheduler, observer access and complete current family. Test the 17/18-direction and rank-fifteen obstructions in that exact domain | A same-object finite theorem map with every additional control and source service accounted for |
| 2. Derive mass and kinetic terms jointly | Extract the charged two-point response and vertices from that continuing law, with one clock and several spatial momenta. Include complementary modes and family symmetry breaking | Nonzero residues and distinct stable gaps or a stated interacting alternative; shared quark kinetic form and consistent Higgs normalization |
| 3. Resolve sector strength and family hierarchy | Test an explicitly derived mechanism outside the excluded affine/trace-one classes, or supply the physical reason their hypotheses fail | One fixed source program determines the full matrices without empirical mass inputs; all extra modes and zero directions retained |
| 4. Establish field and vacuum limits | Supply chiral measure/anomaly control, locality, state/vacuum and controlled spatial/volume/occupation limits on the same regime | A common interacting chiral field, rather than only a finite symbol, a taste-doubled free limit or a selected background |
| 5. Perform a physical comparison | Freeze the construction and measurement map; propagate physical calibration, matching and model uncertainties before comparing | Joint charged masses/mixing, relevant neutrino and interaction tests in compatible conventions; no per-channel target adjustment |

A useful negative result states its exact excluded class. A useful positive
result preserves the full detector/intervention law, including the state after
measurement and subsequent predictions. These criteria apply to any change of
repair law, coherent memory or effective particle interpretation on the fixed
carrier.

## Reproduction and primary evidence

Use the pinned [model requirements](../../model/requirements.txt). From the
repository root, single-threaded numerical libraries make the small-matrix
replays predictable:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q \
  model/common_source/test_native_filtered_current_control.py \
  model/common_source/test_native_family_record_memory.py \
  model/common_source/test_native_gauge_family_records.py \
  model/common_source/test_native_record_yukawa_field.py \
  model/common_source/test_native_family_protocol.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.common_source.native_family_protocol \
  --check model/common_source/native_family_protocol_results.json
python -m model.verify --hashes-only
```

The complete protocol stores 31 array families and source bindings. The replay
recomputes operators and complete probability laws, rejecting alternative valid
normalized experiments. Hash verification checks custody only. Exact identities
and rank claims are checked separately from numerical diagnostics; displayed
ill-conditioned ratios have a portability allowance without relaxing the
primary-array checks.

| Question | Owning evidence and replay entry |
| --- | --- |
| Native objects and actual current | [Gauge object](../../model/common_source/native_gauge_object.py), [current law](../../model/common_source/native_sm_current.py), [matter joins](../../model/common_source/native_matter_attachment.py) |
| Full-current controlled limit | [Controller and source bindings](../../model/common_source/native_response_controller.py) |
| Complete memory/gate/field candidate | [Protocol JSON](../../model/common_source/native_family_protocol_results.json) and adjacent NPZ; `python -m model.common_source.native_family_protocol --check ...` |
| Calibrated action and records | [Matching receipt](../../model/common_source/calibrated_matching_results.json), [observer receipt](../../model/common_source/calibrated_spectrum_observer_results.json), [moving-field receipt](../../model/common_source/calibrated_transport_results.json) |
| Quadratic trial and exact hierarchy bounds | [Candidate package](../../model/common_source/chiral_candidate/README.md), [diagnosis](../../model/common_source/hierarchy_diagnosis/README.md) |
| Complementary response and passive schedules | [Complete-response package](../../model/common_source/full_repair_response/README.md); `python -m model.common_source.full_repair_response.run --verify model/common_source/full_repair_response/results.json` |
| Compiled failing feedback, loops, dressing and bath | [Mass receipt](../../model/common_source/native_repair_mass_results.json); `python -m model.common_source.native_repair_mass --check ...` |
| Local transport and restricted continuum | [Gauge/matter receipt](../../model/common_source/native_gauge_matter_transport_results.json); `python -m model.common_source.native_gauge_matter_transport --check ...` |

The package READMEs contain implementation-level derivations, individual replay
commands, complete retained arrays and their declared domains. Numerical
failures, successful finite constructions and physical identification are
recorded separately; none is erased by a change of notation or interpretation.
