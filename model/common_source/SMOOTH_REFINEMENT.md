# A smooth binary observer with a complete occupation error bound

This is a separate detector variant of the same interacting two-site history.
It replaces the final sharp sign measurement by a Ramsey binary instrument
with effects `(I +/- sin(Q0))/2`. The complete unnormalized classical/quantum
output now has an explicit occupation error bound tending to zero. The new
instrument also preserves finite energy and has a simple unconditional energy
increment in the infinite-oscillator model.

The old sharp-sign protocol and its retained evidence remain unchanged.
Changing the detector here does not prove a sharp-sign error bound. The
separate [complete-output comparison](READOUT.md#complete-comparisons-with-the-actual-finite-pointer)
supplies finite-cutoff sharp-sign enclosures for the zero-current preparation;
it does not supply an interacting all-cutoff asymptotic rate.
This result is a fixed two-site occupation theorem, not a spatial continuum
or an empirical test of OPH.

## The actual smooth instrument

Keep the same Hamiltonian, initial preparation, native-current gauge kick,
first gauge record and conditional feedback as in
[OCCUPATION_REFINEMENT.md](OCCUPATION_REFINEMENT.md). The total Hamiltonian
model time before the detector remains `T=1+1/2`. In particular,

\[
 H=H_g+N+1+\frac{\kappa}{2}(Q_0-Q_1)^2
                  +g(a_0+a_0^\dagger)Z_0,
 \quad \kappa=\frac9{32},\quad g=\frac14.
\]

Append a blank pointer. Apply Hadamard, `exp(-i Z_anc Q0/2)`, `S^dagger`,
then Hadamard and computational pointer readout. Chronological operations are
listed left to right; the total operator is
`H_anc S_anc^dagger exp(-i Z_anc Q0/2) H_anc`. Its two work-space Kraus maps are

\[
 M_0=\frac{e^{-iQ_0/2}-i e^{iQ_0/2}}2
     =e^{-i\pi/4}\cos(Q_0/2-\pi/4),
\]
\[
 M_1=\frac{e^{-iQ_0/2}+i e^{iQ_0/2}}2
     =e^{i\pi/4}\cos(Q_0/2+\pi/4).
\]

Thus `M0* M0=(I+sin(Q0))/2`, `M1* M1=(I-sin(Q0))/2`, and their sum is `I`.
Pointer bit zero denotes the plus effect. This preserves the earlier source
Ramsey convention; it differs from the old sharp detector's meaning for its
bit label. All downstream actions must use the declared bit dictionary of
the new protocol.

The finite instrument replaces `Q0` by its physical number-basis compression
`Q0,d`. Its pointer interaction is exactly the compression of the infinite
pointer generator. Ancilla preparation, basis changes and readout do not act
on scalar occupation. Both outcomes and their quantum outputs are retained;
there is no conditioning or division by their probabilities in the bound.

The source Ramsey construction already exists in
[the scalar source experiment](../observer/source_scalar/experiment.py).
The new [smooth_protocol.py](smooth_protocol.py) attaches it to this common
interacting carrier and its first recorded feedback, rather than importing a
detector probability from a separate model.

## The prior interacting estimate

Let `Pi_d` retain `0 <= n0,n1 < d` and let `N=n0+n1`. The preceding occupation
proof establishes on the same recorded history

\[
 M(t)=\langle(N+2)^4\rangle\le M_0e^{Ct},\qquad
 M_0=\frac{143748821}{8192000},\quad C=\frac{13}{4}.
\]

The expectation includes the sum of all unnormalized gauge-record branches.
It remains uniform over all three native-current gauge kicks and both
settings of the conditional feedback switch. Their operators commute with
`N` and `Pi_d`.

For `r=d/2`, define

\[
 p_d=\frac{81}{80}\frac{\binom{2r}{r}}{4^r81^r},\quad
 e_0=\sqrt{2p_d},\quad A=\sqrt{M_0},\quad A_T=Ae^{CT/2},
\]
\[
 I_T=A\frac{e^{CT/2}-1}{C/2},\qquad
 L_d=\frac{\kappa}{d+2}+\frac{|g|}{(d+2)^{3/2}}.
\]

Here `e0` bounds the vector error of the actual normalized total-occupation
vacuum truncation. The projected pre-pointer error is at most

\[
 e_{\rm prior}=e_0+I_TL_d.
\]

These are the interacting moment and Duhamel estimates from the earlier
proof, not a transfer of its unrelated free Gaussian detector formula.
Self-adjointness, the common finite-number core and the rectangular Galerkin
embedding are established there and unchanged here.

## Include the entire pointer interaction

Use an auxiliary parameter `s in [0,1]` for the pointer coupling. This is a
parameterization of one declared detector unitary, not an additional unit of
Hamiltonian model time or a calibrated physical duration. Its generator is

\[
 G=\frac12 Z_{\rm anc}Q_0
   =\frac{1}{2\sqrt2}(a_0+a_0^\dagger)Z_{\rm anc}.
\]

The number-sector linear-term inequality already proved for `(N+2)^4` applies
to an arbitrary bounded involution on an untouched ancillary system. Since
`1/(2sqrt(2)) <= 1/2`, it gives

\[
 \frac{dM}{ds}\le2M,
 \qquad \sqrt{M(s)}\le A_T e^s.
\]

The number operator acts trivially on that pointer and on the earlier record
register. Accordingly ancilla preparation/basis changes and the final
recording isometry preserve the summed moment. The complete post-readout
history has the explicit finite moment bound

\[
 M_{\rm final}\le M_0e^{CT+2}.
\]

Only the lowering part of `G` can return an outside occupation to the finite
code. As in the earlier cutoff argument,

\[
 \|\Pi_dG(1-\Pi_d)\psi_s\|
 \le\frac12\|\sqrt N1_{N\ge d}\psi_s\|
 \le\frac{A_Te^s}{2(d+2)^{3/2}}.
\]

Integrating over the auxiliary coupling parameter adds the projected
residual

\[
 e_{\rm pointer}\le\frac{A_T(e-1)}{2(d+2)^{3/2}}.
\]

The final infinite outside component satisfies
`||(1-Pi_d) Psi_final|| <= A_T e/(d+2)^2`. It is orthogonal to the projected
difference. Therefore the full pointer-and-record dilation vectors obey

\[
 \boxed{\quad
 \|\Psi_d-\Psi\|\le
 E_d:=\sqrt{\left(e_0+I_TL_d+
                 \frac{A_T(e-1)}{2(d+2)^{3/2}}\right)^2
                      +\frac{M_0e^{CT+2}}{(d+2)^4}}.\quad}
\]

This controls the **actual finite smooth detector**, not the compression of
an infinite outcome projector. Pointer basis changes are common unitaries;
reading both pointer values and dephasing their record is a common channel.
Consequently `min(1,E_d)` bounds the complete recorded-instrument trace
distance, the total variation of its classical record law, and each detector
probability error. The statement is for the declared initial scalar vacuum
family, uniformly over its specified gauge controls; it is not an arbitrary
input diamond-norm estimate.

The convention is `D(rho,sigma)=||rho-sigma||_1/2`. For normalized dilation
vectors, `D(|Psi><Psi|,|Phi><Phi|) <= ||Psi-Phi||`; taking channels and partial
traces can only decrease this distance. Classical total variation uses the
same factor of one half.

The same bound survives the finite classical mixture and record-conditioned
native continuation by convexity and channel contractivity. In particular,
keeping all the native branches does not require multiplying the bound by
their count. The native controller must be the same bit-conditioned channel
on both sides of this comparison; its meaning is the new smooth instrument's
declared record, not an identification with the old sharp-sign outcome.

At fixed `T`, `E_d=O(1/d)`. All constants, including the pointer cost, are
implemented with exact rational arithmetic and directed root/exponential
bounds in [smooth_refinement.py](smooth_refinement.py). There is no positive
rounding floor in the occupation limit.

## Energy and form domains

Remove the two irrelevant constant branch phases and write
`f0(q)=cos(q/2-pi/4)`, `f1(q)=cos(q/2+pi/4)`. They satisfy

\[
 \sum_b f_b^2=1,\qquad
 \sum_b f_b f_b'=0,\qquad
 \sum_b(f_b')^2=\frac14.
\]

For any vector-valued wavefunction in the momentum quadratic-form domain,
the product rule yields

\[
 \sum_b\|\partial_{q_0}(f_b\psi)\|^2
 =\|\partial_{q_0}\psi\|^2+\frac14\|\psi\|^2.
\]

All functions and their first derivatives are bounded. Each branch therefore
remains in that form domain; bounded multiplication also preserves the
quadratic-potential domain. On normalized input, the unconditional kinetic
energy `P0^2/2` increases by exactly `1/8`. In quadratic-form language,

\[
 \mathcal I^*(P_0^2/2)-P_0^2/2=\frac18I.
\]

Every other term of the declared `H` commutes with the `f_b(Q0)`: `P1`, the
coordinate potentials, scalar/gauge coupling and gauge Hamiltonian all have
this property. Thus the unconditional **full** Hamiltonian expectation also
increases by exactly `1/8`, for normalized inputs in its form domain. The
prior history has a finite fourth number moment, which supplies the required
finite kinetic and potential forms. This is a statement about the complete
outcome ensemble. It is not the conditional energy increment of each
individually normalized branch.

The exact `1/8` identity belongs to the infinite canonical oscillator. At a
finite occupation cutoff the canonical commutator has a boundary defect;
the finite instrument is well-defined and has finite energy, but the same
operator identity need not hold. The tests retain this distinction. Energy
injection is stated in this model's dimensionless Hamiltonian convention,
not as a measured laboratory energy or a physical service-cost certificate.

## Temporal refinement, execution and limits

At fixed `d`, the same symmetric source product formula has generator norm
sum bounded by

\[
 S_d\le3+(1+2\kappa)(2d-1)+|g|\sqrt{4d-2}=O(d).
\]

Its error on either model-time interval is at most
`S_d^3 t^3/(3 n^2)`. That comparison uses the **same finite smooth detector**
on both evolutions, so appending the complete instrument contracts the error.
Add the two temporal errors and `min(1,E_d)` by the triangle inequality.
`complete_error_budget(d,steps)` implements this addition for the executed
dimensions `d=2,4`, with floating and physical-service errors excluded.
Its coarse small-register occupation bounds are currently vacuous; having a
controlled limit does not make those retained small registers precise
infinite-detector predictions.

For the ideal family, choosing `n(d)=d^2` on each interval makes the temporal
and occupation errors both `O(1/d)`. Only `d=2,4` source execution is claimed
here. The larger even cutoffs in the rational certificate are mathematical
members of the specified Galerkin family, not unreported execution results.

Run from the repository root:

```bash
python -m model.common_source.smooth_refinement
python -m pytest -q model/common_source/test_smooth_refinement.py
```

The tests check the full ancillary Kraus phases and effect signs, exact
derivative identities, energy injection on independently integrated complex
wavefunctions, the finite-commutator counterexample, directed cutoff
arithmetic and the combined temporal budget. The proof supplies the
all-cutoff statement; finite numerical tests corroborate the formulas.
Neither this variant nor its bound includes the old sharp-sign instrument,
the separately appended modular clock, spatial refinement, physical source
services, or agreement with experimental data.
