# Validation of the native-record electromagnetic adapter

Mathematical validation, 12 October 2026. The P-closure output
and availability of designed hardware are explicit inputs. The question
here is whether their composition gives a real, testable response of the
declared finite model, with its records and controls kept explicit.

The proposed route does: the first native signed current selects the
orientation of a calibrated probe loop; the accepted coupling determines a
binary detector probability; that detector outcome selects the next native
seam. The detector is a physical transcription of an already retained native
record. It is not claimed to add information beyond an ideal reading of that
record or to learn a predictor.

## Exact feedback coefficient from native arithmetic

Let `j` be the first current on `(0,1)`, `z` the load configuration after that
event, and `s=sin(g(P) chi)`. The detector outcome `y=+1` selects `(0,5)`;
`y=-1` selects `(0,1)`. A fair native repair on `(a,b)` has conditional mean
current `(z_a-z_b)/2`, independently of the odd/even rounding case. Therefore

```text
E[j_2 | z,j]
 = (2z_0-z_1-z_5)/4 + (s j/4)(z_1-z_5).
```

Under the uniform binary total-two preparation, the first native repair
preserves the distribution and every mean load remains `1/6`. The first
term consequently averages to zero. The events `j=+1` and `j=-1` each have
probability `5/66`. Conditional on either, the remaining particle is uniform
over ten other ports, so its probability of occupying port 5 is `1/10`.
On `j=+1`, `z_1=1`; on `j=-1`, `z_1=0`. Hence

```text
E[j(z_1-z_5)] = (5/66)(9/10 + 1/10) = 5/66,
E[j_2] = 5 sin(g(P) chi)/264,
E[N_0,final] = 1/6 - 5 sin(g(P) chi)/264.
```

The first of the last two quantities is directly readable from the existing
second native current record. The second is an equivalent source-state audit;
the model needs no additional occupation meter to read the feedback signal.
Also `E[j]=0` and `E[j^2]=5/33`, giving

```text
E[y] = 0,
E[y j] = (5/33) sin(g(P) chi).
```

An independently shuffled detector bit with the same fair marginal gives
`E[j_2]=0` and `E[N_0,final]=1/6`. Thus the response uses the actual
record/detector correlation, beyond a change in the frequency of selecting
the two seams. A separate floor/ceiling enumerator, without importing the
producer, reproduced all these coefficients as exact rational coefficients
of `s` over all `66*2*2*2=528` source/coin/detector paths.

## Binary instrument and preservation of committed records

With `R_z(phi)=exp(-i phi Z/2)`, preparing `|+X>` and measuring Y gives

```text
phi_j = g(P) j chi,
m_y(j) = [exp(-i phi_j/2) - i y exp(i phi_j/2)]/2,
|m_y(j)|^2 = [1+y sin(phi_j)]/2.
```

This fixes the sign convention as well as normalization. The corresponding
old-record Kraus operator is diagonal: `M_y=sum_r m_y(j(r)) Pi_r`.
Completeness gives `sum_y M_y^* M_y=I`. Because each old classical record
projector commutes with every `M_y`, its Heisenberg image under the
unconditioned detector channel is unchanged. Conditioning on the detector
outcome can update a posterior about the earlier current without changing
the retained letter.

The subsequent native operation acts on the source while retaining both the
old letter and the newly committed detector label. It likewise preserves the
classical record projectors. This supports committed classical records, not
arbitrary quantum-prefix coherence. A separate complex-matrix calculation
constructed `R_z`, the Y eigenvectors and `M_y` directly; it checked the signed
probabilities and completeness, preserved every record diagonal for a coherent
input, and exhibited a changed off-diagonal entry.

For a balanced pair of hypotheses `j=+1` versus `j=-1`, detector success is
`(1+sin(g chi))/2` in the declared positive-visibility regime. This does not
postselect away the `j=0` events in the actual experiment: those events have
fair detector outcomes and remain among the 528 weighted paths.

## Charge convention, gauge transformation and exposure

The declared charge assignment is `Q_v=g(P)N_v`. Native conservation gives
`x'-x=j(e_b-e_a)` on a seam. Therefore for arbitrary local gauge parameters
`lambda_v`, coefficientwise,

```text
sum_v lambda_v (x'_v-x_v) = j(lambda_b-lambda_a),
U_lambda K_(c,j) U_lambda^* = exp[i g j(lambda_b-lambda_a)] K_(c,j).
```

This is branch charge-transfer covariance. Checking one rational vector of
gauge parameters is a regression control; the vector identity supplies the
arbitrary-parameter proof. Kraus phases cancel in the record-discarded branch
law. A coherent enlarged isometry with neutral record basis vectors is not
automatically gauge invariant: its record/environment transformation or link
transport must also be specified if that stronger property is wanted.

The loop exposure `chi` belongs to the calibrated, externally prepared
apparatus field. The native record selects its orientation or zero aperture.
The phase consequently contains one factor of probe charge, `g j chi`.
Calling this a field generated by the source charge itself would require a
separate field equation and calibration and can introduce another charge
factor. No such identification is used here.

For an aperture containing a designed number K of canonical cells, the area
dictionary reads `A_loop=K P ell_star^2`. Holding the loop exposure fixed while
varying the trial P requires field calibration proportional to `1/(K P)`.
That is a fixed-flux comparison, not a fixed-field-density comparison; K is a
design choice and the apparatus need not have one microscopic-cell aperture.
The imported low-energy electromagnetic branch and its units must match the
probe, rather than substituting a unification-width coupling or an mZ anchor
without the corresponding transport.

## Binary interpretation and record capacity

The [native observer](../README.md) uses `b=1[j>0]` as its chosen-seam
label. Retaining only that binary label cannot supply predictive gain beyond
a baseline that already knows the action; the fuller native coin/current
record does supply additional predictive information. For this adapter, the actual detector is
the new binary outcome, and its P-dependent correlation with j is measured.

For generic `|sin(g chi)|<1`, every detector sign has positive probability
for every first current. Independent support enumeration gives **24** reached
terminal triples `(r_1,y,r_2)` and **14** reached two-native-letter words if
the detector label is discarded. Reading only the detector gives two labels.
These are support counts for the changed protocol; the native observer's 11-word
capacity and four-bit prefix encoding do not describe this protocol. A capacity
claim must retain this protocol's own continuation and readout contract.

The response calculation supplies a conditional designed instrument and a
native feedback consequence. It does not fabricate hardware, re-establish the
accepted closure inputs, construct a self-generated Maxwell field, or by
itself identify a unique microscopic physical law.

## Executable validation

The six tests use independent native arithmetic and detector amplitudes. They cover the low-energy
coupling conventions, complex detector amplitudes, exact rational response
coefficients, every retained history row, the bit-only actuator on all 78
source inputs, and the sample/sensitivity controls. The separate calculations
above do not import the adapter's history producer for their native event
arithmetic or detector matrix formulas.

For the declared source branch and unit exposure, the point calculation gives
`g` approximately `0.3028676131`, `P(y=+1|j=+1)` approximately `0.6491292616`,
and mean second native current approximately `0.005648835668`. These are
numerical evaluations of the proved formulas. The finite-cutoff coupling
calculation is compared with the imported certificate's midpoint using the
stated tolerance; it is not a new enclosure or a proof that every displayed
digit lies inside the narrow imported interval.
