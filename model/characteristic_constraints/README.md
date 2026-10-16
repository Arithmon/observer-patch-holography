# Two-order characteristic constraints from executable probes

This packet makes an axis-to-interior measurement constraint without a
Gaussian preparation or an oscillator evolution law. It executes both
orders of two noncommuting source-compiled unitary probes, then attaches the
standard quantum Gram-positivity constraint to the published oscillator
characteristic measurements. The laboratory attachment assumes the reported
Weyl displacement instruments and a common normalized preparation.

The new result is a **state-independent constraint and a descriptive
consistency diagnostic**, not a fitted trajectory or an OPH-specific
empirical success. Both unresolved coordinate extents and both operator
orders are retained. The previous Gaussian and anharmonic scores are
unchanged. The input means have sampling and apparatus errors; they are not
asserted to be exact characteristic values.

## The reused theorem and the new attachment

For any normalized positive state functional and unitaries `U,V`, write
`a=<U>`, `b=<V>`. Positivity of the operator Gram matrix gives

```text
|<UV>-ab| <= sqrt[(1-|a|²)(1-|b|²)],
|<VU>-ab| <= sqrt[(1-|a|²)(1-|b|²)].
```

For the first inequality take the three operators `I,U*,V`. Their Gram
matrix is

```text
G = [[1, conjugate(a), b],
     [a, 1, <UV>],
     [conjugate(b), conjugate(<UV>), 1]].
```

Its determinant is `(1-|a|²)(1-|b|²)-|<UV>-ab|²`. Exchange `U,V`
for the second inequality. This is a small quantum Bochner constraint;
it reuses positivity, not a new axiom about observers. A primary reference
for Weyl relations and quantum Bochner positivity is Lami, Sabapathy and
Winter, [*All phase-space linear bosonic channels are approximately Gaussian
dilatable*, Lemma 2 and Appendix A](https://arxiv.org/pdf/1806.11042).
The finite algebra prerequisite is already represented by
[StateExpectation](../../Lean/EventAlgebra/StateExpectation.lean). This
packet supplies a mathematical attachment and executable checks, not a new
Lean formalization of the continuous Weyl theorem.

For the laboratory convention `D(alpha)=exp(alpha a*-conjugate(alpha)a)`,
put `U=D(u)`, `V=D(iv)` for real `u,v`. Then

```text
UV = exp(-iuv) D(u+iv),   VU = exp(+iuv) D(u+iv).
```

Consequently `chi(u+iv)` lies in the **intersection** of disks centered at
`exp(+iuv) chi(u)chi(iv)` and `exp(-iuv) chi(u)chi(iv)`, each with the radius
above. Changing the unresolved coordinate direction swaps the two disks;
it cannot justify selecting the order giving a smaller residual. This
extends the earlier [collinear moment disk](../measurement_bounds/README.md)
to two noncommuting displacement directions. The disks are generally broad;
they do not determine a unique interior mean.

## What the executable observer actually runs

The finite source demonstration chooses `U=exp(i pi Z/7)` and
`V=exp(i pi X/5)` on a declared proper two-dimensional code. The existing
source compiler constructs these words and their controlled versions from
native pulses. A pointer starts at zero, receives a Hadamard, a controlled
word, a phase shift, and another Hadamard. Its two complete measurement
maps are `(I+exp(i phase)W)/2` and `(I-exp(i phase)W)/2` in exact arithmetic.
Zero phase reads the real mean; phase `-pi/2` reads the imaginary mean.

The runtime executes all eight word/quadrature settings, retains both
outcomes, and checks entangled-reference execution and serialized checkpoint
continuation. The matrix products are read from executed source words;
target Pauli matrices are not substituted for those words. Chronological
`V` then `U` implements `UV`, and both orders are run. Raw numerical maps
are neither polarized nor renormalized. Their mass/unitarity defects and
Ramsey residuals are reported separately from the exact theorem.

The chosen code and coherent loading/control/export services are interface
inputs. This two-dimensional demonstration verifies the general unitary
instrument construction. It does **not** implement the exact continuous
Weyl relations: finite matrices cannot represent the oscillator canonical
commutation relation exactly. The laboratory Weyl identification therefore
remains an explicit apparatus premise; it is not derived from this qubit
experiment. No physical service-fidelity certificate is available here.

## Public data and results analysis

The data are the authenticated, unchanged zero-duration means from
[McGarry et al.'s deposition](https://zenodo.org/records/18752287), described
in their [measurement protocol](https://arxiv.org/html/2603.04744v1).
The existing loader verifies the deposited hashes. Exact decimal strings
are then read as rational numbers. The procedure uses the 21 axis records
and 100 off-axis complex records of one quadrant; exact mirrored copies
are not additional observations. The measured origin is preserved, not
used to rescale the data. Unit Gram diagonals express normalized-state
mathematics; they do not replace the noisy origin measurement.

For each extent, predictions consume only the axes and coordinate indices;
their output is sealed before the scorer consumes interior responses.
This is a retrospective split by a declared rule. The public values were
already available and inspected during development, so this is not a
preregistered or prospective validation. No outcome is fitted and no
coordinate extent is selected by its score.

At extent **4**, one of the 200 order-specific plug-in Gram matrices is
certifiably nonpositive. At extent **5**, three are. In both cases the
strongest retained witness is the `(1,1)` interior point in the positive
product-order convention. A rational witness proves that making its three
complex means jointly physical requires changing at least one by at least
**0.02439** for extent 4, or **0.03129** for extent 5, in complex absolute
value. These are conservative lower bounds on deterministic correction,
not measured detector biases or statistical significance levels.

The certificates evaluate `v*Gv` with rational arithmetic and outward
Taylor intervals for the Weyl phase. If all three complex entries changed
by at most `epsilon`, then `|v*DeltaG v| <= 2 epsilon ||v||²`. A negative
certified quadratic form therefore gives the reported correction floor.
Floating eigenvalues only provide diagnostics. No shot counts are inferred,
no p-value is assigned, and ordinary measurement noise could account for
the incompatibilities. The deposited aggregate means do not provide the
authenticated accepted counts/raw trials needed to decide that question.

## What it means and next steps

The preparation-specific Gaussian relation can now be checked alongside a
weaker constraint applying to every state measured by the declared Weyl
instruments. This helps separate a failure of Gaussian modeling from a
failure to interpret noisy means as exact quantum expectations. It does
not establish that the dataset is unreliable, select an OPH repair law,
or settle the unresolved coordinate extent.

The next empirical step is to obtain accepted counts or raw trials,
coordinate calibration and detector-error bounds. With those inputs the
same Gram constraints can be applied to confidence regions, allowing a
statistical test without fitting an evolution law. The independent
[closed native-law protocol](../closed_prediction_protocol/README.md)
provides a more discriminating future experiment once a physical mapping
of its services is available.

```sh
python -m model.characteristic_constraints --check model/characteristic_constraints/results.json
python -m pytest -q -W error model/characteristic_constraints
```
