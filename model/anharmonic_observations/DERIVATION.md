# Executed pulse dynamics and the empirical comparison

This is a conditional reconstruction of the apparatus in McGarry et al.,
[*Programmable quantum simulation of anharmonic dynamics*](https://arxiv.org/html/2603.04744v1).
The Hamiltonian, preparation family, measurement interface and laboratory
calibrations are supplied. The new result is execution of a complete
anharmonic observation protocol through the existing native source unitaries,
followed by comparison with excluded positive-duration observations. It does
not select this action from OPH or establish a spatial continuum theory.

## Preparation and the quantities measured

Let `a` be the oscillator lowering operator, `N=a†a`, and
`x=(a+a†)/sqrt(2)`. At zero simulated duration the preparation family is

\[
 \rho_0=D(\alpha)\rho_{\rm thermal}(\bar n)D(\alpha)^\dagger,
 \qquad \chi_0(ib)=e^{-(\bar n+1/2)b^2+2i\alpha b}.
\]

Only the explicitly labelled zero-duration column pair in the deposited
one-dimensional characteristic scan enters the two-parameter fit. All fifty
complex measurements participate equally. The fitted displacement is real,
restricted to `[-2,2]`; occupation is nonnegative. The procedure fits no
coordinate scale, gain, offset, potential parameter or duration. The analytic
Gaussian expression is used for preparation calibration, not to supply the
later predictions. The older two-dimensional comparison remains unchanged.

The measurements are binary spin averages, `Re chi(ib)` or `Im chi(ib)`.
Their outcome-zero probabilities are `(1+mean)/2`. The later state is generally
non-Gaussian, and its characteristic function must be calculated from the
executed state. For the position time series, an ordinary least-squares slope
of the thirteen imaginary means against `b`, divided by `sqrt(2)`, is a
derived position diagnostic. It is not an additional independent observation.

## Continuous pulses and conventions

The paper supplies `delta=2*pi*500 s^-1`, `Omega=pi/(150 us)`, pulse duration
`tau=50 us`, virtual spin angle `theta=0.8`, and symmetric-well phase zero.
One complete step has simulated duration `4*tau=200 us`. In the simulated
Schrödinger frame a pulse with sign `s` uses

\[
 H_s=\delta N-s\frac{\Omega}{2}X_{\rm spin}(a+a^\dagger).
\]

This includes the continuously changing drive phase within each pulse. An
instantaneous displacement held at the beginning of an entire step would
omit one of the paper's stated Trotter errors. The chronological sequence is

\[
 H_-\ (\tau),\quad R_z(\theta),\quad H_+\ (2\tau),\quad
 R_z(\theta),\quad H_-\ (\tau).
\]

Here `down=|0>`, `Z|0>=|0>`, and `Rz(theta)=exp(-i theta Z/2)`. At zero
harmonic frequency and small rotation angle the selected-spin generator is
the **positive position cosine**, `cos[2*alpha0*(a+a†)]`, where
`alpha0=Omega*tau/2=pi/6`. The executable tests this sign and quadrature.
The intended approximate scalar Hamiltonian is therefore

\[
 H_{\rm target}=\delta N+B\cos[2\alpha_0(a+a^\dagger)],
 \qquad B=\theta/(4\tau)=4000\ {\rm s}^{-1}.
\]

The prediction executes the pulse sequence rather than substituting this
target Hamiltonian. Equations 12 and 41 in the paper have inconsistent spin
rotation signs under their displayed definitions; its displacement phase
also needs a quadrature offset to realize the stated position potential.
Our convention is fixed by that target and by initialization/readout, before
scoring later data. It is an explicit apparatus interpretation, not a
parameter fitted to a trajectory. Spin conjugation gives an equivalent
description with `down=|1>` and the opposite rotation sign.

For the non-stroboscopic position series, detector quadratures refer to the
simulated Schrödinger frame. At 2 and 4 ms the harmonic frame rotation is
itself the identity. The deposited coordinates are literal `Im(beta)` in
both files; they are not doubled or rescaled by `sqrt(2)`.

## What the native source executes

The occupation cutoff is a power of two. Native Fock words and a compiled
displacement prepare every branch of a supplied capped geometric thermal
mixture, without dropping a preparation branch. The experimental first
selection is represented by this calibrated initial state.

In the spin-X basis, `H_s` has two forced-oscillator blocks related by
oscillator parity. A canonical eighty-digit eigensystem supplies a real
basis for `N-(Omega/delta)*(a+a†)/2`. That basis and its inverse are compiled
into the existing native source pulses; the eigenphases, spin rotations and
controlled parity are likewise executed. Their contracted matrices implement
each finite `exp(-i H_s dt)`. Independent SciPy exponentials check the result;
they do not replace the source in the prediction path. No infinite-dimensional
canonical-commutator identity is assumed at the finite cutoff.

Reusing a contracted matrix is a host optimization. Each application still
has its complete word cost counted. There is no claim that the native pulse
count measures elapsed laboratory or observer time. The model duration comes
from the externally calibrated experimental schedule.

The source retains the full spin–oscillator density matrix across all steps.
At a requested endpoint, let `rho_s=<s|rho|s>`. The selection probability is
`p_s=Tr rho_s`. Both blocks and both masses are retained. For each source
pointer outcome `j`, executed oscillator Kraus operators give

\[
 p(s,j)=\operatorname{Tr}(M_j\rho_sM_j^\dagger),\qquad
 p(j\mid s=0)=p(0,j)/p_0.
\]

Only the conditional mean is compared with the experimentally selected
data. The selection is terminal; an endpoint readout does not collapse the
state used for a later fresh-trial prediction. Intermediate spin projection
or reset would implement a different experiment. Tests explicitly distinguish
these procedures and check joint/conditional normalization.

## Supplied noise and numerical checks

The measured dephasing rate is `gamma=18 s^-1`, with Lindblad operator
`sqrt(gamma)*N`. The exact dephasing-only channel is the Schur multiplier

\[
 (\mathcal D_h\rho)_{sn,rm}
 =e^{-\gamma h(n-m)^2/2}\rho_{sn,rm}.
\]

This channel acts on spin off-diagonal blocks as well as diagonal blocks.
It is supplied environmental physics, not a native-derived noise law.
Symmetric half-dephasing/unitary/half-dephasing splitting approximates the
simultaneous pulse master equation. The primary calculation uses eight
subdivisions per 50 us, with sixteen as a numerical refinement. Occupation
32 is checked against occupation 64. All refinements retain the same fitted
parameters and all detector settings.

An independent reference integrates the unsplit finite Lindblad equation
using directly constructed matrices. Small-dimension tests check zero-noise
agreement and second-order splitting convergence. Cross-cutoff and subdivision
differences are numerical diagnostics, not certified unrestricted error bounds.
Source contraction residuals are also numerical checks, not physical error
certificates. No calibration confidence region, accepted-shot likelihood or
full apparatus-error budget is established here.

## Excluded observations and controls

The protocol fixes all 200 real/imaginary means in the full 2/4 ms scans and
all 507 imaginary means at the 39 positive-duration position-scan times.
Every one is scored. Prediction functions accept settings and zero-duration
calibration only; the complete prediction payload is hashed before validation
loading. The public experiment and its qualitative outcome were already
known, so this is a retrospective postdiction, not a blinded discovery.

Fixed controls are the unevolved prepared state, harmonic evolution with the
same dephasing, and the complete pulse sequence without dephasing. Neither
the controls nor refinements are selected by agreement. Reported scores are
raw spin-mean RMSE, corresponding binary-probability RMSE, and all individual
residuals. Software `all_pass` certifies execution/comparison contracts only.
