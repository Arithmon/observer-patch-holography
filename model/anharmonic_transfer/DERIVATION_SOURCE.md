# Native execution of the asymmetric apparatus

The supplied apparatus acts on a spin tensor a finite occupation space.
Preparations, sideband unitaries and detector instruments reuse the unchanged
native primitives of `anharmonic_observations`. This packet adds finite carrier
rotations and replaces the inserted analytic dephasing multiplier with the
finite source environment described in `DERIVATION_ENVIRONMENT.md`.

## Phase and pulse chronology

Write `q_raw=a+a†`, `N=a†a` and

\[
 H_s=\delta N-s\frac{\Omega}{2}Xq_{\rm raw},\qquad
 R_j(\theta)=e^{-i\theta\sigma_j/2}.
\]

One step executes, in chronological order,

\[
 H_-(\tau),\ R_x(-\varphi),\ R_z(\vartheta),\ R_x(\varphi),\
 H_+(2\tau),\ R_x(\varphi),\ R_z(\vartheta),\ R_x(-\varphi),\
 H_-(\tau).
\]

At zero harmonic drift, to first order in `vartheta`, its spin-zero block is

\[
 I-i\vartheta\cos(2\alpha_0 q_{\rm raw}-\varphi)+O(\vartheta^2),
 \qquad \alpha_0=\Omega\tau/2.
\]

Thus the executor parameter `phi` represents the negative of the potential's
phase `Phi`. This sign follows from the native matrices before comparison
with observations. A test differentiates the executed small-angle sequence
and rejects the opposite sign. The published pulse conventions contain sign
inconsistencies already disclosed in the original symmetric packet; this
extension retains that packet's positive-cosine convention.

## Finite rotations and two duration parameters

For a spin rotation by angle `v`, the declared laboratory carrier rate is
`Omega0=pi/(35 us)`. Its duration is `g=abs(v)/Omega0`. The continuously
advancing drive frame executes

\[
 U_g=\exp[-i\delta Ng]R_x(v),
\]

with both commuting factors contracted from native words. The noise operator
`N` commutes with this gap Hamiltonian. Consequently a full finite environment
channel composed with `U_g` needs no noise/unitary splitting in the gap.
The virtual `Rz` has zero duration.

The published coordinate duration is `4*K*tau`; the chosen continuously
advancing drive duration is `K*(4*tau+4*g)`. These are separately retained.
For `abs(phi)=pi/20` or `pi/10`, each nominal 200 us step lasts 207 us or
214 us on this drive clock. Neither duration is host runtime or native pulse
count. The paper's static detuning motivates continuous phase accumulation,
but compensation during carrier gaps is not authenticated by the deposited
metadata. The instantaneous-rotation sensitivity fixes every `g` to zero
while preserving the spin angles; it is not fitted to later measurements.

## Environment, selection and measurement

Each sideband subpulse executes half a finite environment channel, the native
unitary, then another half environment. There are `8*K*s` such environment
draws for `s` subdivisions per 50 us pulse. Four further full-gap environment
draws occur per step, including the declared zero-duration branch schedule
when a gap vanishes. Every application refreshes its independent classical
draw. The environment is discarded without postselection. The full spin and
oscillator density operator continues through every step.

Endpoint records retain both spin-selection masses and all four joint
selection/detector probabilities. Conditional readout divides only the spin
zero branch by its mass. Recording an earlier endpoint never resets the
continuing state: each requested endpoint represents a fresh-trial experiment.

The detector coordinate is literal `beta=u+i*v`, not a pulse displacement
`beta/2`. For the published two-point finite difference use `u=0`, `v=h=.4`
and `theta=-pi/2`. Its spin mean is `Im chi(i*h)`, and the reported position
estimator is this mean divided by `sqrt(2)*h`. The exact occupation-space
position expectation is retained separately; these quantities need not agree.

## Error scope and execution ledger

Each record composes the exact-rational ideal Gaussian-quadrature radii for
all environment applications with the finite-occupation Strang certificate.
The certificate includes arbitrary retained reference systems. Because both
selection branches remain present, a complete trace-distance radius `epsilon`
and selected mass `q>epsilon` give conditional detector TV at most
`epsilon/q`. The returned acceptance is a floating numerical mass, not a
certified lower endpoint. Its harmless above-one roundoff is capped at one
for the ideal-bound calculation; no floating arithmetic certificate is
claimed. Occupation removal, physical calibration, preparation uncertainty,
and authenticity of the supplied noise model are outside these bounds.

The resource ledger counts every sideband word, virtual rotation, finite spin
gap word and realized environment branch for each application. It separately
reports the host work of contracting all environment branches. Cached native
matrices save repeated host computation without removing word applications
from the declared executable process.

Independent checks integrate complete sidebands with unsplit DOP853 Lindblad
dynamics, exponentiate the spin-gap Hamiltonian directly, compare full joint
and selected readouts, demonstrate second-order splitting convergence, test
the zero-phase reduction to the frozen symmetric runner, and distinguish
continuing spin evolution from per-step selection/reset.
