# Charged realization, binary instrument and feedback

## Conserved charge on the native branches

Let N_v be the native occupation and choose `Q_v=g(P)N_v` in
Heaviside–Lorentz units. The original integer repair gives exactly
`y-x=j(e_b-e_a)`. Thus the physical transfer on any branch is
`ΔQ=g(P)j(e_b-e_a)` and total charge is conserved. This holds on the full
78-state nonnegative total-two space; the scored preparation uses its
66-state binary shell.

For arbitrary port gauge parameters `lambda_v`, define
`U_lambda|x>=exp(i g sum_v lambda_v x_v)|x>`. Each current-grouped Kraus
operator satisfies

`U_lambda L_(c,j) U_lambda†
 = exp(i g j(lambda_b-lambda_a)) L_(c,j)`.

This follows coefficientwise from the exact displacement equation, not just
the particular rational gauge regression in the code. Therefore each
outcome map and the record-discarded channel are gauge covariant. A coherent
retained Stinespring isometry additionally needs the record/link factor to
transform with the opposite phase; merely declaring an inert environment
would not make it gauge invariant. One may assign that compensating record
phase or include the corresponding Wilson transporter. The scored protocol
uses classical current outcomes and needs branch covariance only.

## Designed electromagnetic probe

The probe is a unit-charge path qubit. Its native-record controlled router
sets the signed closed-loop integral to `j chi`, where the apparatus prepares
chi independently of g. In natural units the gauge-invariant relative phase
is `phi_j=g(P)j chi`. Under `A -> A + d lambda`, the added closed-loop
integral vanishes. The design assumes a coherent two-path source, routing,
a reference phase and Y-basis readout. It does not claim that native charge
itself creates the flux; such an assertion would require a field-source
map and could introduce another factor of g.

Prepare `|+X>`, apply `exp(-i phi_j Z/2)`, and measure Y with outcomes
`y=+1,-1`. The amplitude and effect on the old current register are

```
m_y(j) = [exp(-i phi_j/2) - i*y*exp(i phi_j/2)]/2
M_y = sum_r m_y(j(r)) |r><r|
E_y = M_y† M_y
    = sum_r [1+y sin(phi_j(r))]/2 |r><r|.
```

For every real g and chi, these eigenvalues lie in [0,1] and sum to one.
Consequently `sum_y M_y†M_y=I`: this is a valid binary instrument on every
input density matrix, including coherent records. Appending y to a fresh
orthogonal two-state register makes the detector result a retained classical
fact. All M_y commute with old record projectors, so the old classical
record probabilities are protected. The instrument can change off-diagonal
record coherences; no unknown-quantum-state nondisturbance is claimed.

The further native repair is controlled by y alone. On the enlarged space
its fixed channel acts as the identity on the committed bit within each
bit block and applies the appropriate complete native instrument to the
source. Thus its Kraus completeness follows blockwise, it preserves the bit
projectors, and rereading y is repeatable. Composing the first native CPTP
instrument, the detector and the bit-controlled second instrument is CPTP.
No selected branch is normalized in place of the complete experiment.

At g=0 or chi=0 the bit is fair. Reversing j swaps outcome probabilities.
`d Pr(+|j=1)/d chi` at zero is g/2. For the separate calibrated setting
`chi=pi/(2g)`, the signed j=+1 and j=-1 hypotheses are distinguished with
certainty; j=0 remains fair. That setting is a hardware calibration example,
not the main fixed-exposure comparison and not evidence that P dependence
disappears for fixed apparatus.

## Exact feedback coefficient

Set `s=sin(g(P)chi)`. Since native j is -1,0,1 on this source,
`Pr(y|j)=(1+y j s)/2`. After the first repair let the source be z. The
second seam is `(0,5)` when y=+1 and `(0,1)` otherwise. Averaging its fair
coin gives

`E[j2|z,j]=(2z0-z1-z5)/4 + (s*j/4)(z1-z5)`.

The first event preserves the unconditional uniform binary-shell state,
whose port means are 1/6. Thus the first term averages to zero. First
currents +1 and -1 each have mass 5/66. On either such branch, the other
particle is uniform over ten ports, so its port-5 occupation has mean 1/10.
For positive current, `(z0,z1)=(0,1)` and
`E[j(z1-z5)|j=1]=9/10`; for negative current `(z0,z1)=(1,0)` and that
expectation is 1/10. Hence

`E[j(z1-z5)]=5/66`, and `E[j2]=5s/264`.

Also `E[y]=sE[j]=0` and `E[yj]=sE[j²]=(5/33)s`.
An independent bit with the same fair marginal removes the s-dependent
correlation and restores zero mean second current. The unchanged bit
marginal rules out an explanation in terms of changed seam frequencies.
The weak-exposure expansion gives
`g E[j2]=(5g²/264)chi+O(chi³)=(5pi/66)alpha(P)chi+O(chi³)`.

These are analytic exact identities in s. Decimal evaluation of the
executed source map and trigonometric function is numerical, with an
independent existing certificate comparison for alpha and separate complex
amplitude checks for the probe. No new interval enclosure is asserted.

## What is observation and what is a designed actuator

The probe is a noisy physical transcription of its already available native
input j. By itself it adds no information about j. Its new outcome is a
binary electromagnetic record that actually determines the next repair.
The independent feedback quantity j2 is then read from the next native
event. The detector response is distinct from the native observer's prediction
score; this instrument is not a learned predictor.

The one-bit action controller and the full recorded trajectory have
different accessible algebras. At nonzero unsaturated exposure the complete
audit supports 24 distinct terminal triples `(r1,y,r2)`, but the actuator
receives only y. A one-cell loop aperture is likewise distinct from the area
of the cut supporting either record algebra. The finite protocol halts
after its second native repair. Its charge convention, loop hardware,
measurement and control are part of the accepted design, rather than a
claim that those components were fabricated by the simulation.
