# Complete contact spectra and an energy readout

This CP-07 / [#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143)
component supports the common-excitation work under
[#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
and [#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140).
It uses the same supplied Weyl/contact Hamiltonian, CAR dictionary and
coordinate parameters as [the bound-pair construction](README.md). It does
not identify an elementary species or predict a physical mass spectrum.

## The physical finite fiber

On a periodic cube, the first-quantized relative-momentum representation has
four spin amplitudes for every relative momentum k. A physical two-fermion
state satisfies

```text
psi(K-k) = -SWAP psi(k).
```

[spectrum.py](spectrum.py) constructs this subspace before diagonalization.
For a two-element orbit `{k,K-k}`, its four basis columns place `e_a/sqrt(2)`
at k and `-SWAP e_a/sqrt(2)` at its partner. A fixed point permits only the
spin singlet. These columns are orthonormal and exhaust the plus-one
eigenspace of the signed exchange operator. At the three retained L=4
momenta, there are eight fixed points and 28 two-element orbits, giving
120 physical states. The redundant 256-dimensional first-quantized space
is used only as an independent evolution reference.

The tests assemble the free operator independently from
`interacting_field.lattice.symbol`, form its full exchange projector and
verify invariance of the physical space. Momentum must belong to the
actual finite reciprocal grid; a continuous Brillouin-zone parameter is
rejected by this finite-spectrum interface.

## Every spectral weight stays in the observation law

Let d_K be the normalized on-site contact vector. For the complete physical
eigensystem `(E_j,phi_j)`, define

```text
w_j = |<phi_j,d_K>|^2,
mu_K = sum_j w_j delta_(E_j),
A_K(t) = <d_K,exp(-it H_K)d_K> = sum_j w_j exp(-it E_j).
```

The measure is positive and has total weight one. Its first moment is
`<d_K,H_K d_K>=U`; its second moment is `||H_K d_K||^2`. These two moments
are checked against direct matrix action. All physical eigenvalues are
retained, including dark eigenspaces of zero contact weight. For portable
display, numerically degenerate eigenvalues are grouped with a declared
`1e-11` tolerance and their multiplicities and weights are summed. This
grouping is a numerical presentation choice, not an exact degeneracy proof.
Evolution always uses the complete ungrouped eigensystem.

The isolated upper eigenvalue and its weight agree independently with the
rank-one resolvent calculation in `fiber.py`. At L=4 and K=0, its energy
is approximately `1.08859448778` and contact weight `0.920510889037`; the
remaining contact spectral weight is approximately `0.0794891109633`.
At `K=(0,pi,pi)`, the weight is approximately `0.962780719101`.
At `K=(pi,pi,pi)`, the analytic identity gives energy U and weight one.
These are lattice model energies and spectral weights, not measured masses
or detector efficiencies.

The complete contact/complement instrument after evolution has outputs

```text
psi_yes = A_K(t) d_K,
psi_no  = exp(-it H_K)d_K - psi_yes.
```

Both vectors are unnormalized and retained. Their squared norms sum to one,
and their sum is independently checked against a full-matrix exponential.
The two orthogonal projectors also preserve any untouched reference when
both outcomes are kept; an entangled-input test checks this directly.
The retained durations 0, 2, 7 and 20 are Hamiltonian model coordinates.
At K=0 and t=2 the contact and complement probabilities are approximately
`0.769052305842` and `0.230947694158`.

Dropping the scattering amplitude is already wrong at t=0: retaining only
the isolated pole would predict contact probability `Z^2` instead of one.
Renormalizing the projected bound state is a different preparation, whose
contact probability is Z and whose preparation success/failure probabilities
are Z and `1-Z`. Both those filter outputs are included in the report.
They cannot replace the original contact-state experiment. Deleting only
the interaction `U|d_K><d_K|` on the same physical fiber removes the upper
isolated pole; the complete free spectrum lies within its original norm
bound `2 sqrt(3)c`.

These field-level projectors and preparations are supplied mathematical
observables. This component computes their actual finite-Hamiltonian laws
but does not compile a many-site spatial detector into the native source.
The physically distinct executed native clock experiment is joined below.

## What the native clock actually reads

At the eta momentum, a contact-return detector alone has probability one
for every t: it cannot observe the global phase of a single eigenstate.
The [native clock construction](../particle_clock/README.md) instead uses
the invariant pair `(r,d)` with energies `(0,U)`. Both levels have particle
number two and the same total momentum. An ideal coherent preparation
`(r+d)/sqrt(2)` has return probability `(1+cos(Ut))/2` under the actual fiber
Hamiltonian; the test evaluates this with the full matrix exponential.

The supplied source-to-field isometry maps native source states `{1,2}`
and `{7,2}` to r and d. The native reference has corner mass `264/7505`
and conditional density `diag(25/33,8/33)`. Thus on this M2 algebra,
`beta U=log(25/8)` and Hamiltonian model time is `t=-beta s` for modular
parameter s. The report replays the complete existing ideal native
instrument at s=1, including both heralds, all native continuation records,
energy readouts and aborts. It does not divide by herald success.

For C, the unconditional energy-r probability is approximately
`0.0180073391952`; for E it is approximately `0.0174371772182`.
The initial rejected mass is approximately `0.964823451033` and later
outside mass approximately `0.00497446147013`, both retained. Together with
the energy-d and abort outcomes, each complete probability law has mass one.
The separate `particle_clock/results.json` binds the expensive source-word,
record-training and preparation certification and its complete error budget.
This spectral join replays the exact ideal instrument, without pretending
to repeat that source certification or claiming a native clock for every K.

## Volume and continuum boundaries

The report evaluates pole and residue on even torus sides 4, 6, 8, 12 and 16
at the same three reciprocal momenta. At K=0, the L=16 values are
approximately `1.08745032119` and `0.923443831096`. These numerical volume
comparisons disclose the boundary/volume dependence. No error theorem is
inferred from apparent convergence, and the larger-volume rows are pole
diagnostics rather than a replay of every eigenstate or detector.

All these volumes keep the lattice spacing and action fixed. A physical
vacuum, chiral spatial-continuum limit, common physical units and a matched
Standard Model spectrum remain obligations of the owning completion gates.
The result supplies one internally consistent lattice spectral measure,
its complete observation law, and a bounded finite clock attachment.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/bound_pair/test_spectrum.py
```
