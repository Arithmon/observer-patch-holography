# A joint quantum response and observer-clock test

This experiment puts the two observer clocks and a coherent detector response
in one complete native history. An independently enumerated replay agrees with
both representations of that history, including every unnormalized output
matrix, both detector outcomes and every timeout. A control that reproduces
the clock and detector marginals separately fails their joint distribution.
This is a finite operational constraint on a correspondence dictionary. Physical
clock, field and spacetime identification requires further evidence.

The nearest earlier result is the classical stopped correspondence in
[`observer_correspondence`](../observer_correspondence/README.md). The present
package retains all source coherences and an untouched quantum reference. Its
producer establishes the transported-cone identity on every matrix unit;
[`clock_response.py`](clock_response.py) provides an independent endpoint
enumeration and a joint clock/detector test of that identity.

## Declared source, clock and apparatus

The source is native C on the 66 binary states with total occupation two. A
scheduled attempt chooses one of the thirty seams uniformly and a fair coin.
At fixed seam and coin, matrix entries survive together precisely when their
native current records agree. The observer regions are `A=(0,1)` and `B=(0,5)`.
Each retained repair records its incidence on both regions and its occupation
changes on their union `(0,1,5)`. The two clocks count these incident repairs,
including the native feedback repair. These counts are dimensionless.

The initial source state is

\[
 |\psi\rangle=(|0,1\rangle+|0,5\rangle)/\sqrt2,
\]

where each ket names its two occupied ports. The protocol requests a first
union incident repair with a cap of two scheduled attempts. On success it
executes native C feedback at seam `(0,1)` if the first recorded change at port
zero was positive, or at `(0,5)` otherwise. It then applies the declared union
apparatus

\[
 U_\gamma=R_{1,5}(-1)^{\gamma n_0n_1},\qquad\gamma\in\{0,1\},
 \qquad
 R\big|_{|10\rangle,|01\rangle}
 =\begin{pmatrix}3/5&-4/5\\4/5&3/5\end{pmatrix}.
\]

The rotation is the identity on `00` and `11`. The apparatus is supplied and
charged as an intervention; it contributes no native incident count. No claim
about its physical duration follows from that convention. A second union
incident repair is requested with the same cap. Finally, both outcomes of the
`n1` occupation projection are retained, including after a timeout. A first
timeout ends the repair protocol before feedback or the apparatus.

The coherent preparation and this apparatus are explicit premises. Native C
preserves the diagonal subalgebra and does not generate this preparation from
an initially diagonal mixture. The producer also checks that the occupation
phase cannot be represented as a product of one-port occupation phases over
the full sector. Its interaction is an apparatus intervention, not an action
or Hamiltonian selected by native repair.

## Independent enumeration and retained reference

The replay calls the native endpoint rule separately for each basis state,
seam and coin. It pairs bra and ket endpoints only for equal current records,
then sums hidden seam and coin alternatives. It enumerates the actual remote
waits up to each cap, applies the rational apparatus directly, and retains
every terminal density matrix. It imports no producer Kraus groups, resolvent,
operator blocks or transformed instruments. The underlying native endpoint
function and seam geometry remain shared source dependencies.

There are 1,442 terminal output labels, including zero-probability labels.
For each setting `gamma=0,1`, the replay checks the complete matrices against
both the native and effective execution, with and without a two-dimensional
reference. These are eight exact comparisons using rational arithmetic.
The package's all-input matrix identity provides the arbitrary-reference
statement; these explicit preparations are independent execution witnesses.

For the reference preparation
`(|0,1>|0> + |0,5>|1>)/sqrt(2)`, discarding the untouched reference gives the
dephased initial source process exactly. Its source-only clock/detector law
therefore cannot see `gamma`, while its complete source/reference output
changes. The retained entrywise matrix difference is `176/625`; this diagnostic
is not a trace norm or an operational distinguishability claim.

A second reference preparation,
`(|0,2>|0> + |0,3>|1>)/sqrt(2)`, verifies surviving entanglement rather than merely
nonzero off-diagonal entries. Its first-timeout, detector-zero branch has mass
`9/25`. The expectation of its partially transposed, unnormalized density in
`(|0,2>|1> - |0,3>|0>)/sqrt(2)` is `-17/200`, or `-17/72` after conditioning on
that branch. The independent endpoint replay reproduces the whole branch.
The timeout remains in the full experiment.

## A constraint that separate marginal fits miss

Write the complete joint law as
`P_gamma(status, ticks_A, ticks_B, detector_bit)`. Its status masses are
`9/25` for the first timeout, `144/625` for the second timeout and `256/625`
for completion. They sum to one; no successful-run renormalization enters the
joint-law comparisons. The complete status/clock marginal is identical at
both apparatus settings, but

\[
 \operatorname{TV}(P_0,P_1)=\frac{1848}{15625}.
\]

For example, completion with clock pair `(1,3)` has probability `16/625` at
both settings. Its conditional detector-one probability changes from `23/50`
to `7/10`. Clock statistics and coherent response are thus checked on the same
recorded histories, with no independent time calibration for the detector.
The identical clock marginal holds by construction: both clocks tick on
scheduler-selected seam incidences under a state-independent scheduler, so
the status/clock marginal is a function of the seam sequence alone. The joint
law therefore tests the correlation between seam records and the detector,
not an independent clock discriminator.

Construct a control separately at each setting by imposing conditional
independence within each terminal status:

\[
 Q_\gamma(s,c,b)
 =\frac{P_\gamma(s,c)P_\gamma(s,b)}{P_\gamma(s)}.
\]

This control preserves all status-conditioned clock and detector marginals,
including the timeout mass. Nevertheless it differs from the native law:

| Setting | `TV(P_gamma,Q_gamma)` | Conditional clock/detector covariance |
| --- | --- | --- |
| `gamma=0` | `168097/6000000` | `3269/983040` |
| `gamma=1` | `231073/6000000` | `7877/983040` |

The covariance uses the indicators of clock pair `(1,3)` and detector one,
conditional on completion. It vanishes for each control. The full-law total
variation distances use every status. Replacing the initial coherent state by
its diagonal mixture removes the entire phase response, including all
conditional matrices, exactly.

A bijective relabeling of the recorded clock pairs preserves these total
variation distances and the corresponding event covariance. A change of
clock units cannot repair the control's failure. A many-to-one clock
coarsening can erase the distinguishing information: forgetting both clocks
makes this particular control agree. The executable test demonstrates both
facts. For any proposed observational law `D`, the triangle inequality also
gives `max(TV(D,P_gamma), TV(D,Q_gamma)) >= TV(P_gamma,Q_gamma)/2`.
This bound states what an eventual common record dictionary must distinguish;
no observed data have been inserted into it.

## What the result supplies

The experiment supplies one joint clock/field-response contract on a fixed
native source and verifies that its effective description preserves the
whole quantum output. Overlapping observers suffice; a federation is not a
premise. It does not derive the apparatus, the coherent input, an SI clock
scale, Lorentzian spacetime, a continuum field theory or the observed Universe.
Connecting those requires a source-selected observable dictionary that
preserves such joint experiments, together with its stated approximation
domain. Matching isolated clock rates and field frequencies would leave the
correlation tested here unexamined.

Run the independent tests and the package receipt replay from the repository
root:

```sh
python -m pytest model/quantum_observer_correspondence/test_clock_response.py -q
python -m model.quantum_observer_correspondence --check model/quantum_observer_correspondence/results.json
```
