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

Let `B` denote the executed oscillator basis word, `D` its eigenphase word,
`A` the spin-axis word and `P` the twice-applied controlled parity word.
The ideal inverse-word identities give

\[
 A P(I_2\otimes BDB^{-1})P A^{-1}
 = I+A\left[P^2-I+P(I_2\otimes B(D-I)B^{-1})P\right]A^{-1}.
\]

The negative-sign pulse applies the further identity
`F U F^-1 = I + F (U-I) F^-1`. Evaluating these increments avoids repeating
the floating mismatch of separately contracted inverse pairs at every short
step. The literal `P^2-I` term remains, including binary-angle phase error.
Both forms use the executed source matrices and retain every physical stage.
The original straight product and the increment form must both pass the
finite-Hamiltonian comparison, with a separate fixed representation-change
guard. No target matrix replaces an executed word. These ideal algebraic
identities do not cancel physical service errors, and the floating residuals
are not interval certificates.

The finite increment representation is unitary in exact arithmetic, but its
computed matrix can retain a small Hermitian defect. Production therefore
applies the existing fixed-guard polar representation to that matrix, after
the mandatory raw and increment checks. Its input-defect, correction and
restored-defect limits are unchanged; a further guard limits the total
corrected-versus-raw change to `2e-10`. The corrected full-phase Hamiltonian
comparison is mandatory too. No evolved state, output probability or target
matrix enters the polar calculation. All three representations' diagnostics
remain numerical quantities rather than outward error certificates.

Two source words have exact specialized forms. Writing
`N=sum_b 2^b*(I-Z_b)/2` yields one scalar and one phase pulse per bit.
Occupation parity depends only on the least significant bit, so its spin
control is the existing native `J(-pi)` followed by the scalar `I(pi)`.
These identities preserve full phases and act on the entire finite code;
they do not invoke unrestricted oscillator commutation relations.

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
