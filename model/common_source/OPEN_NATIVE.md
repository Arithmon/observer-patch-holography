# Native memory around a mixed continuing work history

`open_native.py` attaches the existing complete native instruments to the
mixed-state execution in [`open_protocol.py`](open_protocol.py). The same
prepared two-site work register receives the first native current, all four
canonical observer reads, their committed-record feedback, and the declared
unobserved phase channels. The second smooth record selects the final native
seam. There is no work reset or replacement between these operations.

This instantiates the work-only CPTP memory result already proved in
[`ENCODED_MEMORY.md`](ENCODED_MEMORY.md). It supplies executed mixed outputs
for that result; it does not prove a new error-correcting code or protect a
physical native memory against unspecified faults.

## Executed input, complete output and ordering

The input is a normalized native/reference amplitude matrix, with either
78 native input rows (`full_carrier`) or two logical input rows
(`protected_qubit`, the default, and `passive_qubit`). Reference columns
may be complex and entangled with the native input. The work state is the
declared independent preparation. This attachment's public API does not
accept initial native/work correlations or replace that preparation by a
supplied density. The mixed work runtime separately accepts arbitrary
work/reference densities; that broader API is not an extra input contract
for the factorized native attachment.

Each continuing raw output has label

\[
 h=(c_1,a_1,r,c_2,a_2),\qquad r=(r_0,r_1,r_2,r_3).
\]

Both native coin values and every raw record code are retained. A first
invalid raw code halts, with label `(c1, a1, (), None, None)`. It retains
the already prepared, un-kicked work state and pays its preparation cost;
no work operation or environment event, including one at cursor zero,
runs on that halted path. A final raw abort retains the work history
already executed. Native outputs keep all 512 code rows, including
padding. Ideal supported inputs give zero mass to the invalid outcomes;
the numerical output does not discard their residual amplitudes.

For a nonaborted first event, its decoded current \(j_1\) determines the
work kick. The four-read program in [`COMPOSITION.md`](COMPOSITION.md)
then runs with the declared environment schedule. The native register is
retained during this work. The final seam is \((2,3)\) for \(r_2=0\) and
\((0,1)\) otherwise in both encoded modes; the full-carrier alternative
instead uses \((0,5)\) for \(r_2=0\). Protected mode applies the existing
source recovery after each valid native event. Passive mode applies neither
recovery. The attachment stops before terminal native decoding or a
destructive logical read.

Let \(\Xi^{j}_{r}\) be the completely positive work operation for visible
record \(r\), including the current kick, subsequent work operations,
and all unobserved phase mixtures. For a work preparation \(\sigma\),
the retained density is

\[
 \rho^j_r=\Xi^j_r(\sigma),\qquad
 \sum_r\operatorname{tr}\rho^j_r=1.
\]

The environment alternatives are summed as a CP map, not as amplitudes.
They are not additional observer records. The density is generally mixed
and remains unnormalized separately in each visible branch. The source
finite phase mixture and the analytic Gaussian comparison are distinct
backends; their parameters and numerical/ideal bounds are described in
[`OPEN_ENVIRONMENT.md`](OPEN_ENVIRONMENT.md).

Write \(A_h\) for the complete native branch map, including both fair-coin
amplitudes, any encoder and the requested recoveries. `output` streams

\[
       (h,\; A_h\psi_{MR},\;\rho^{j_1}_r).
\]

This represents
\( |A_h\psi\rangle\langle A_h\psi|\otimes\rho^{j_1}_r\), with the
native/reference factor listed first. Neither factor is normalized again,
and neither is replaced by a selected eigenvector. The full tensor density
need not be materialized. Immutable native matrices and work densities
are shared across labels where their actual executed factors coincide.

`audit` compares the complex native columns against independent endpoint
operators and checks

\[
       \sum_h\operatorname{tr}(\rho^{j_1}_r)A_h^\dagger A_h=I
\]

with the halted-path terms included. This is trace preservation on all
native/reference inputs with the specified independent work preparation.
It does not establish an all-input work channel from one prepared density.
The digest binds the complete menu, arrays and work history seals; it
detects changed envelopes, not an adversary who also replaces the digest.

## Existing memory theorem, now with executed CP mixtures

The exact native theorem is stated for the ideal source operators. With
the code isometry \(V\), recorded native operation \(L_a\), recovery
\(R_a\), and recorded probability \(p_a\), it gives

\[
        R_aL_aV=\sqrt{p_a}\,V,\qquad\sum_a p_a=1.
\]

The code, exact \(p_a=k/4\), endpoint support checks, source encoders and
recovery construction are those in `ENCODED_MEMORY.md`. The present
execution reuses their full source maps; endpoint targets remain audit
oracles. No desired recovery matrix is inserted to replace source execution.

Choose any Kraus representation \(B^{j}_{r,\nu}\) of the complete work
instrument, where \(\nu\) includes discarded environment alternatives.
For the protected two-event history the complete ideal joint Kraus map is

\[
 M_{a_1,r,a_2,\nu}
 =\sqrt{p_{a_1}^{01}p_{a_2}^{s(r_2)}}\,
       V\otimes B^{j(a_1)}_{r,\nu}.
\]

The native coin is already part of \(a_i\) and \(p_{a_i}\); there is no
additional coin probability to multiply into the output. Summing the
second native outcome gives one for each visible work record. Summing
\((r,\nu)\) is trace preserving for every first current, and summing
\(a_1\) gives one. Consequently the logical action, before encoding the
output, is \(\mathrm{id}_L\otimes\Psi_W\) for a complete work-to-work/record
channel \(\Psi\). For an arbitrary initial density \(\omega_{LWR}\),

\[
 \operatorname{tr}_{W,\mathrm{records}}
  \mathcal C(\omega_{LWR})
       =(V\otimes I_R)\omega_{LR}(V^\dagger\otimes I_R).
\]

This algebraic statement permits initial logical/work/reference correlations.
The factorized executable API above tests its independent-work specialization.
Conditioning on a retained record may steer a correlated logical/reference
input; the statement concerns the complete, unconditioned marginal.
It requires a CPTP work apparatus acting only on work, correct classical
control and the specified native operations. It does not require the work
Hamiltonian, quadrature, time evolution or preparation to approximate a
particular physical model accurately. A work approximation can therefore
change the full output without changing this exact ideal memory marginal.

`logical_entanglement_fidelity` evaluates the actual complete encoded Bell
overlap after tracing work and records. `phase_distance` computes the
complete CQ trace distance of logical \(|+\rangle\) and \(|-\rangle\)
outputs without terminal decoding. Native branch rank at most two allows
the trace norm to be evaluated by a thin QR decomposition, while the work
factor contributes its trace. No branch probabilities are renormalized.

For passive mode, the earlier endpoint calculation gives exact phase
erasure on the \(01\to23\) path and perfect phase distinguishability on
\(01\to01\). For independently prepared work let

\[
 q_j=\sum_{r:r_2=0}\operatorname{tr}\rho^j_r,
 \qquad \bar q=(q_{-1}+2q_0+q_{+1})/4.
\]

The first-event current probabilities are \((1/4,1/2,1/4)\) for both
native laws. Thus the existing formula remains
\(D_{\rm passive}=1-\bar q\), whereas \(D_{\rm protected}=1\).
Contractivity and the triangle inequality imply that any decoder has
worst-case trace-distance error at least \(\bar q/2\) on these two inputs.
`memory_witness` uses the actual mixed work probabilities in this formula.
Its floating diagnostic is not a directed rational lower certificate.

## Resource and error boundaries

`resources` enumerates compatible complete raw paths, including zero-weight
scheduled paths and first aborts. It reports the componentwise maximum of
each resource and a path attaining the native-pulse maximum. It charges
one actual native encoder in encoded modes, first and selected final native
repair words, valid selected recoveries, the corresponding work path, and
every scheduled independent environment invocation. Source contraction
caches do not reduce the cost of repeated physical invocations. Environment
draws are discarded, so they do not enlarge the retained observer record.
The four work pointers, their four committed bits, and one work preparation
are retained in the work ledger. Native coin/record initializations and
retained record width are reported separately.

An unknown logical input is supplied. No three-pulse diagnostic preparation,
terminal source decoder, logical rotation or terminal read is charged or
executed by this attachment. Host density-tree storage and the physical
cost of unquantified services are not included in the source pulse count.

Native logical storage through the work interval remains an explicit
service requirement. Work-only noise does not model noise on that memory,
record erasure, incorrect controls, faulty recovery, or coupling of the
environment to logical information. Unknown service errors remain unknown.
A physical telescope needs errors for the appropriate whole corrected native
stages and storage/input services on their stated domains, plus the actual
work/noise service errors. Numerical branch checks do not supply these errors.

The local phase variance is a supplied dimensionless finite-register
parameter at a declared operation boundary, not a calibrated physical
diffusion rate or continuous Lindblad evolution during a gate. The exact
noiseless smooth-history certificate is not silently reused as a certificate
of the modified history. Any comparison with the earlier prepared infinite
history must additionally charge the declared environment perturbation and
retain the original preparation, time and domain hypotheses. A finite
phase-to-Gaussian bound compares those finite channels, not an occupation
limit or a laboratory environment.

Native source guards, work completeness checks, endpoint comparisons,
Bell overlaps and phase distances are floating diagnostics. Their tolerances
are not mathematical intervals or physical service bounds. In particular,
the exact memory theorem applies to the ideal source maps; its numerical
instantiation agrees only to the recorded execution diagnostics. No full
operational observer, physical calibration, arbitrary native fault
tolerance or repeated-cycle implementation is asserted here.
