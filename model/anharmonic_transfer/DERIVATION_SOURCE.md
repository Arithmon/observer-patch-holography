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

## Evaluation of short source pulses

Let `B` denote the ideal oscillator basis word, `D` its eigenphase word,
`A` the spin-axis word and `P` the twice-applied controlled parity word.
The executor contracts each word and its requested inverse independently.
For the ideal words, the inverse identities give

\[
 A P(I_2\otimes BDB^{-1})P A^{-1}
 = I+A\left[P^2-I+P(I_2\otimes B(D-I)B^{-1})P\right]A^{-1}.
\]

The negative-sign pulse applies the further identity
`F U F^-1 = I + F (U-I) F^-1`. Evaluating these increments avoids repeating
the floating mismatch of separately contracted inverse pairs at every short
step. The literal `P^2-I` term remains, including binary-angle phase error.
Both numerical forms use the executed source matrices and retain every
physical stage. However, their rounded inverse matrices need not multiply
exactly to identity. If the executed inverse is `B^dagger+E`, the increment
contains the extra term `B(D-I)E`; it need not be unitary even when the
straight product's factors are numerically unitary. Its error can change
norms as well as phases and accumulate under repeated use. A small
single-pulse target residual alone does not control the complete process's
probability mass. Different subpulse durations must be checked: the factor
`D-I` is much smaller in a finely subdivided history.

The composed increment therefore passes through the existing compiler's
guarded polar representation of an executed unitary. Its fixed limits are
`1e-10` for input unitarity defect, `5.1e-11` for polar change and `2e-13`
for restored defect. The original straight product, the increment before
polar correction and the final matrix must each pass the unchanged `3e-10`
full-phase finite-Hamiltonian check. Both the increment-minus-raw and final
matrix-minus-raw norms must remain below `2e-10`. The receipt records these
separate residuals, the polar correction and the final unitarity defect.
An incorrect unitary or global phase cannot be repaired by polar projection;
the target and complete-change gates remain mandatory.

No target matrix replaces an executed word, and no state, selected branch
or probability is renormalized. Word hashes, physical stages and resource
counts are unchanged. These ideal algebraic identities do not cancel
physical service errors. All reported norms are floating diagnostics, not
outward intervals or physical precision certificates.

The independent [pulse controls](../anharmonic_observations/test_pulse_stability.py)
compare complete finite Hamiltonian matrices and repeated evolution of an
independently constructed energy eigenstate. A nearby unitary inverse with
a small phase offset exposes accumulated mass loss without the final polar
step, despite passing the single-pulse target check. Zero-duration and tiny
nonzero-duration controls also check the retained phase and generator.
The complete detector-law checks remain unchanged through the actual
observation and transfer histories; no improved empirical agreement follows
from repairing this arithmetic.

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
