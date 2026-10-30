# Complete-field particle and hole response

This package identifies actual excitations of the supplied interacting field:
two charge-one doublets and two charge-two singlets about **one specified
half-filled ground state**, with nonzero probe residues, charge, parity and an
explicit particle–hole map. It also measures the corresponding two-point
response through the existing source operations. It advances CP-05/#1141 and
CP-07/#1143; it does not close their physical particle, mass or defect mapping
requirements.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_response
python -m pytest -q -W error model/particle_response
```

## The same field, with two explicitly different references

The Hamiltonian is the complete-Fock Weyl hopping plus onsite contact term from
[`full_fock_field`](../full_fock_field/README.md). Every number sector and every
transition remains in the spectral calculation; there is no one-particle or
low-energy truncation. Numerical evaluation admits up to eight field modes.
The source readout admits one or two complete sites, retains the purifier, and
adds an initialized Ramsey pointer.

`FieldSpectrum` declares its generator `K = H - mu*N` and reference together:

* The faithful `beta=1, mu=0` Gibbs state continues the reference used by
  [`field_thermal`](../field_thermal/README.md). Its signed Lehmann frequencies
  are thermal transition energies, not vacuum-relative masses.
* `dimer.setup()` separately prepares the unique pure ground state at
  `mu=U/2`. The chemical potential and ground reference are supplied benchmark
  inputs. A pure ground state is not faithful on the complete field algebra;
  it cannot replace the Gibbs state in a theorem requiring faithfulness.

The executable runs source probes on both references. It never uses a thermal
pole list as if it were the half-filled ground spectrum.

## Two charged candidates and their pair channels

For a two-site open field along axis `a`, set `t=c/2>0` and
`Delta=sqrt(U^2+16*t^2)`. The onsite change of basis
`d_0=c_0, d_1=i*sigma_a*c_1` turns the hopping into the spin-independent Hubbard
dimer. The internal SU(2) generators in this frame commute with the complete
Hamiltonian. This is a finite internal spin algebra; identifying it with
Lorentz spin in three spatial dimensions requires an additional map.

The half-filled singlet block in the separated-singlet/symmetric-doublon basis
is

```text
H_singlet = [[0, 2t], [2t, U]],
ground(K) = -(U + Delta)/2.
```

If its normalized ground coefficients are `a,b`, then
`a²=(1+U/Delta)/2`, `b²=(1-U/Delta)/2`, and `a*b=-2t/Delta`.
These give the following projectors and nonzero local probe overlaps directly.
The code checks their complete matrices rather than assigning labels to
eigenvalue indices.

| Excitation above the same ground | Excess CAR number | Parity | Internal S² | Multiplicity | Model-frequency gap | Local probe residue |
|---|---:|---:|---:|---:|---|---|
| Hole doublet, N=1 | -1 | odd | 3/4 | 2 | `Delta/2-t` | `(1+4t/Delta)/4` |
| Particle doublet, N=3 | +1 | odd | 3/4 | 2 | `Delta/2-t` | `(1+4t/Delta)/4` |
| Pair hole, N=0 | -2 | even | 0 | 1 | `(U+Delta)/2` | `(1-U/Delta)/4` |
| Pair particle, N=4 | +2 | even | 0 | 1 | `(U+Delta)/2` | `(1-U/Delta)/4` |

At `c=1/8, U=1`, the two charge-one gaps equal **0.4528882032022076**,
with per-mode residue **0.310633906259083**. The pair gap is
**1.0153882032022077**, with local pair residue **0.007464374963667**.
An explicit unitary on all 16 Fock states sends `N` to `4-N`, preserves `K`,
and interchanges the two charged projectors. Their equal gaps are therefore a
model symmetry consequence, not two independently fitted values.

A second charge-one doublet occurs at `Delta/2+t`, with local residue
`(1-4t/Delta)/4`. Its separation `2t=c` from the lower doublet resolves the
kinetic frequency in the same observable. The higher residue vanishes in the
free `U=0` control and is positive for `U>0`. This is a two-site frequency
splitting, not a continuum dispersion relation. The neutral triplet gap
`(Delta-U)/2` is smaller than the charged gap, so the latter is not the full
theory's lowest gap either.

The **conditional electron-like hole / positron-like particle benchmark**
assigns physical charge `e*(N-2)` and requires supplied electromagnetic and
internal-to-Lorentz spin maps plus one common calibrated energy unit. These
are declared candidate assignments, not derived identifications of the
electron and positron. The packet predicts their dimensionless model gap
ratio `1`; it predicts no absolute physical mass. The field coordinate value
`gap/c²` is reported as such. There is no topological-defect attachment,
protection theorem, recovered physical clock or continuum rest-mass claim.
The charge-two channels have no assigned Standard Model species.

## Complete spectral measure and independent checks

For creation probe `O` with `[N,O]=q*O`, the complete arrays are

```text
omega_nm = E_n - E_m,
forward_nm = p_m * |O_nm|²,
reverse_nm = p_n * |O_nm|².
```

They give `G+(s)=Tr(rho O†(s) O)` and
`G-(s)=Tr(rho O O†(s))`, both with phase `exp(-i*omega*s)`.
Here `s` is the Hamiltonian model coordinate. The code retains zero-weight and
forbidden-charge entries, verifies the number-sector support, and checks:

* Zeroth and first moments against direct operator traces and commutators.
* The single-mode CAR sum rule `G+(0)+G-(0)=1`.
* Gibbs detailed balance `reverse=exp(-beta*omega)*forward`.
* Independent Jordan–Wigner matrices and direct time-domain exponentials,
  on all three dimer orientations and a four-site complete Fock space.
* Closed atomic Hubbard weights, the free/interacting extra-pole control,
  the explicit particle–hole unitary and all three internal SU(2) generators.

`Measure.lines()` only groups numerical frequencies for display. Its declared
visibility threshold neither removes transitions from the measure nor proves
exact zeros, exact degeneracy or instrumental detectability. The analytic
dimer formulas provide the nonzero-residue and equal-gap statements used for
the two candidates.

## Executed source readout

`source.prepare` compiles the canonical purification of the declared reference
from an initialized blank, retaining every first export-failure output. Public
`correlator` accepts no replacement state or claimed preparation counters.
The internal reuse bundle is immutable; independent Ramsey settings each count
their own preparation.

For Pauli probes `P,Q`, the source executes controlled `Q`, the existing
Hamiltonian word with the counted commuting `mu*N` phases, controlled `P`,
the inverse word, and pointer readout. The two branches implement
`I` and `U† P U Q`. Two pointer phases measure the real and imaginary parts
of `Tr(rho P(s) Q)`. Scalar phases in the controlled words are retained.
The implementation keeps both unnormalized pointer outputs, the purifier and
all preparation and continuation failures. No success branch or particle
outcome is conditioned on.

Writing a lowering probe as `A=(X+iY)/2`, U(1) invariance gives
`C_XX=C_YY, C_YX=-C_XY`. Consequently

```text
G+ = (C_XX - i*C_XY)/2,
G- = conjugate((C_XX + i*C_XY)/2).
```

Four complete Ramsey settings reconstruct a single-mode response; sixteen
settings reconstruct an onsite pair. The single-mode controlled Majoranas are
odd operators on the encoded field. Their supplied logical access does not
derive a parity-preserving physical detector or reservoir; the net two-point
observable is charge-neutral. Pair probes are even.

For the default half-filled dimer at model coordinate `s=1/4` and four splitting
steps, the run uses nine initialized bits, retains four purifier bits and
counts **52,484 source pulses** including four independent preparations. Its
addition/removal response differs from the complete Lehmann value by less than
`6.2e-7`; the complete output-vector error is below `7.3e-6`. Every setting
retains total output mass one within the declared host tolerance.

The exact rational splitting estimate bounds the ideal source word's operator
error. Finite precision source execution and preparation are separately checked
numerically against complete output vectors; the numerical tolerance is not a
new exact compilation theorem, laboratory error bar or physical clock
calibration. Tests inject genuine export leakage and reject it, along with
external preparation substitution and malformed inputs.

## Periodic momentum diagnostic

`python -m model.particle_response.dispersion` explicitly extends the same
z-axis Weyl/contact action to a periodic **four-site** ring. All 256 Fock states
remain. A signed fermionic translation matrix supplies momentum, and the
eigensystem resolves every number, component and translation sector. The
Fourier convention `c†_k=sum_x exp(-ikx)c†_x/2` gives translation eigenphase
`exp(+ik)` and free symbol `+c*sigma_z*sin(k)`; the opposite Fourier convention
changes this sign without changing the action. No momentum is assigned to an
open-dimer eigenvalue.

The ring's half-filled ground is unique at `U=1`. Charge-one probes give a
lowest doublet-sector gap **0.4175457854** at `k=0,pi`. At `k=pi/2`, the two
components give **0.4723570583** and **0.5944851107**, interchanged at `3pi/2`.
Leading residues are approximately `0.334530`, `0.674414` and `0.325544`.
Addition at `k` and removal at `-k` have equal complete pole residues and gaps;
each mode satisfies the full addition/removal CAR sum rule. Lower dark
spin-3/2 levels remain in the catalogue. The exact odd-number projector
`(15/4-S²)/3` and the singlet-ground probe map select spin 1/2 without deleting
those levels or numerical amplitudes.

The declared lattice benchmark `E(k)²=m²+c²*sin(k)²` fails both with the fixed
dimer mass anchor and with the ring's own zero-momentum gap while holding the
supplied speed `c=1/8`. These are separate comparisons: the ring/dimer rest-gap
ratio is **0.9219621585**, a volume/reference change, whereas the lowest
envelope's same-volume squared-frequency residual exceeds **0.03315**.
A fitted envelope speed is **0.2208544942**, with fitted slope **3.1217092863**
times `c²`. There are only two distinct `sin²(k)` values, so fitting mass and
speed succeeds by construction and supplies no dispersion evidence. It does
not independently recover an effective velocity, a photon or a physical clock.

The component-resolved poles do not share one fitted shell. However, the
untwisted internal transverse spin generators shift crystal momentum by `pi`;
they are not a Lorentz-spin dictionary commuting with translations. The result
rejects this **fixed finite benchmark**. It does not exclude an emergent
Lorentz regime with renormalized velocity, a different finite-volume limit,
the continuum theory or OPH. No source preparation/probe run is claimed for
this larger ring. Independent full Jordan–Wigner matrices, Gibbs references
on all sectors, time-domain response and the free Bloch symbol check the
calculation. The free control retains both zero-energy nodes and explicitly
reports its 16-fold grand-canonical (six-fold half-filled) ground degeneracy.
