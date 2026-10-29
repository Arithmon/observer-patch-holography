# Model foundation audit: histories, references and energy

This second audit starts from merged main `991c23a5d`, following the
[integration audit](MODEL_INTEGRATION_AUDIT_2026-10-29.md). It reviews the
finite model interfaces required before particle masses can be attached.
Independent reviewers challenged the native-history representation, energy
identification, full-field clock budget and regional proof arguments.
The results below have their own explicit domains. No physical particle,
calibrated mass or continuum construction is inferred from finite code tests.

## Corrections

The field resource counter accepted negative or fractional split counts and
Boolean counts, producing negative or fractional operation totals. It checks
the same positive-integer domain as the executed evolution, before counting.
Finite signed Hamiltonian model durations remain valid.

The native occupation restriction accepted negative source/reference indices,
Boolean indices and nonfinite coefficients. It validates the entire source
operator before restriction. Arbitrary exact, unnormalized and nonpositive
operators remain admitted, with complex operators represented by their exact
real and imaginary parts. Validation does not discard an invalid entry merely
because its partial trace would vanish.

The earlier branch-compatibility note is marked historical, with navigation
to the merged audits. Frozen source evidence and scientific outputs are
preserved; current receipts are renewed only after producer replay.

## Exact histories with smaller host storage

[Compiled native histories](../model/compiled_native_history/README.md)
reuse the original interpreter's controller, FIFO, timeout and resource-stop
logic. Each complete native branch is a scalar times a partial injection.
Tagged diagonal input columns determine that map and scalar; the proved
branch form then determines every off-diagonal and reference-valued output.
This inference is invalid for a general quantum channel and is not used for
additional coherent source services.

The retained two-attempt feedback control has 625 raw histories, including
zero outcomes. It needs 1,056 map entries instead of 38,788 explicit Choi
coefficients. Complete Choi comparisons, complex references, product and
federation carriers, holds and terminal outcomes are checked independently.
This reduces host evaluation storage; it changes no physical resource count
and does not remove exponential growth in the number of distinct histories.

## Fewer full-field source pulses

The [reduced clock schedule](../model/band_clock_budget/README.md) uses the
proved centered spectral bound `||H-(35/96)I|| <= 77/96`. At modular parameter
`|s|<=1`, its complete average half-diamond bound is below `0.037991` with
21 symmetric split steps. Twenty steps exceed the former 32-step ceiling.
This is minimal within this bound and ceiling, not a global optimality claim.

The complete schedule drops from 26,238,982 to 17,219,334 supplied source
pulses: 9,019,648 fewer pulses. All 8,128 two-particle states, 64 band momenta,
scattering components, tastes, reference factors and four pointer outcomes
remain. The original schedule remains a reproducible baseline. Acquisition
uses 31,457,280 native attempts separately; global Gibbs preparation, packet
preparation and physical calibration remain premises. The ideal bound does
not include floating contraction or hardware error.

## Energy comparisons needed before mass attachment

[Energy identifiability](../model/energy_identifiability/README.md) separates
three missing comparisons.

1. With the same initial block-diagonal state and sector-preserving evolution
   and instruments, adding an arbitrary scalar in each conserved sector
   leaves every complete history unchanged. A separate clock does not change
   this result. Coherent cross-sector access or an independently established
   Gibbs preparation law with temperature changes the experimental premises.
2. If the generator within each of 64 field fibers is known up to a scalar,
   63 coherent comparisons arranged as a tree determine the relative offsets.
   Extra edges test consistency. These are mathematical matrix-unit access
   conditions; their native source implementation is not supplied.
3. Even connected momentum comparisons cannot see `H -> H+cN` on the fixed
   `N=2` carrier. It adds `2cI` and changes neither normalized Gibbs states nor
   dynamics. The internal eta-pair/dark difference `U` remains observable,
   but is not a vacuum comparison. Dispersion curvature has a separate
   kinematic interpretation and does not establish that missing energy gap.

The native localized coherence and its reference have equal diagonals after
every complete native history. Their same-instrument pushforwards therefore
have zero modular-energy difference on the image support at every finite
history length. This needs no stationary reference on a growing tape.
Positive relative entropy and a contrast lifetime in native attempt count
do not turn that zero into a physical mass.

## Regional reference compatibility

The [regional criterion](../model/regional_modular/README.md) works on the
actual fixed-number direct sum, retaining every allowed regional number.
A faithful state admits a preserving UCP retraction onto the regional
number-preserving algebra exactly when each sector is a product
`r_n tensor tau_n`. For nonzero Gibbs parameter, this is equivalent to a
block-additive Hamiltonian. Boundary hopping violates the necessary center
condition: the exact field witness is `||[N_A,H]||_HS^2 = 117/16`.
It also witnesses failure of exact generator intertwining for any such
retraction. No numerical Gibbs logarithm is needed for this conclusion.

Tracial generator covariance and Gibbs preservation are distinct demands.
The existing tracial expectation remains valid; it cannot silently become
a weighted Gibbs-preserving expectation. Partial traces, actual reduced
densities, regional modular clocks, approximate maps and other algebras are
outside the exclusion. Where nested regions do admit expectations preserving
one faithful reference, uniqueness supplies the tower law. Independently
chosen regional references need not compose.

A faithful-state-preserving UCP idempotent automatically has a fixed range
closed under ordinary products. Calling the same exact range an operator
system supplies no escape. A nonfaithful positive control shows why the
faithfulness hypothesis matters. These finite matrix arguments are explicit
analytic proofs with exact/numerical controls, not new Lean certificates.

## Acceptance boundary

The audit narrows several attachment choices and removes redundant evaluation
work. It does not select the physical field action, prepare its global Gibbs
reference, identify a physical clock, construct a protected moving native
particle or calibrate mass. The canonical premise register and physical
completion gates retain their acceptance requirements. Exact finite
corollaries and source-bound diagnostics are recorded as component evidence.
A finite audit provides evidence against particular defects; it cannot
certify the absence of every future defect.
