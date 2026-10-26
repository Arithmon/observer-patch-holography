# Why the declared Weyl–Higgs candidate misses the particle hierarchy

**Research diagnosis — 26 October 2026. Unreleased calculations.**

The immediate failure is in the candidate's **mass-generating prescription**.
The observer reads the implemented quadratic action correctly, but that
action makes masses proportional to a nearly degenerate occupation feature.
Its fixed Higgs background supplies no further interaction that could turn
those coefficients into the measured hierarchy. This diagnoses a rejected
architecture, not a contradiction of every OPH construction.

This note extends the
[candidate report](OPH_CHIRAL_WEYL_HIGGS_CANDIDATE_2026-10-26.md).
The original source, note and quantum-state receipts remain unchanged.
New exact calculations, tests and the reproducible diagnosis are in
[hierarchy_diagnosis](../../model/common_source/hierarchy_diagnosis/README.md).

## 1. Locate the failure before changing the readout

The declared action uses

\[
Y_u=aF,\qquad Y_d=bF,\qquad Y_e=JF,\qquad H=(0,v),\quad v=\Omega.
\]

Here the native source computes

\[
a=\log(25/8),\quad b=\log(8/5),\quad J=a-b,\quad\Omega=a+b,
\]

and the eigenvalues of the actual family feature are

\[
f_1=\frac{61}{450},\qquad
f_2=\frac{41}{225}-\frac{7\sqrt5}{750},\qquad
f_3=\frac{41}{225}+\frac{7\sqrt5}{750}.
\]

The step that turns this occupation feature into a left–right mass coupling
is a declared choice. The source derivation of \(F\), internal chiral
representation and invariant Higgs tensors does not prove that
\(Y_e=JF\) is the physical coupling.

With canonical kinetic terms the charged-lepton mass parameters are
\(m_i=vJf_i\). Their ascending ratios are approximately
\(1:1.1903:1.4982\). The retained clock records independently recover those
same ratios. The observed target is approximately \(1:206.77:3477.23\);
the audit uses the much wider windows \([206,207]\) and \([3470,3485]\),
so experimental uncertainty is immaterial to this exclusion. These are
comparison inputs only. See the PDG
[muon](https://pdgweb.lbl.gov/2025/listings/rpp2025-list-muon.pdf) and
[tau](https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf) listings.

The current 90-mode field law is quadratic in CAR operators once the Higgs
background is fixed. The earlier four-mode response construction had an
explicit density–density interaction. The current chiral extension uses
the number \(J\) as a coefficient of a bilinear pairing; it does not retain
that four-mode interacting dynamics as a field interaction. Nor does it
execute quantum Higgs or gauge fluctuations. A logarithmic statistic from
the native reference is not, by itself, a derivation of a fermion self-energy.

Consequently there is no uncomputed many-body mass enhancement inside this
particular quadratic Hamiltonian. Its dispersion is already exactly
\(E_i(p)^2=c^2|p|^2+m_i^2\). A genuinely interacting extension could change
that conclusion, but would be a new dynamics to derive and test.

## 2. The clock, pair factor and vacuum cannot repair this action

For a selected rest pair the implemented Hamiltonian is

\[
H_d=d\,a^\dagger b^\dagger+d^*ba.
\]

Its four Fock energies are \(-|d|,0,0,+|d|\). A charged odd-CAR excitation
above its ground state costs \(|d|\), whereas an initially empty pair has
occupation oscillation frequency \(2|d|\). The observer inverse already
accounts for this factor of two. The normal CAR response has nonzero
residues at both signed poles; the three family gaps are not hidden by
vanishing detector overlap.

More generally a quadratic Nambu field satisfies

\[
\Xi(t)=e^{-i\mathcal B t}\Xi(0),\qquad
\{\Xi(t),\Xi^\dagger(0)\}=e^{-i\mathcal B t}.
\]

Thus the retarded linear-field response has state-independent pole
locations for this fixed law. Changing the reference state can change
other correlation functions and occupations; it cannot create a different
mass spectrum in this response. This statement does not extend to a
state-dependent interacting self-energy. The time parameter here is the
declared coherent model coordinate, compared with the calibrated modular
clock in the retained experiment, not laboratory seconds.

A common lapse, clock regraduation, Higgs normalization or conversion to
physical units multiplies all three masses by the same factor. It leaves
their ratios unchanged. A change of field basis also preserves poles when
both kinetic and mass forms are transformed consistently.
The native band Gram matrix is exactly \((10+2\sqrt5)I_3\), so the
implemented family frame conceals no anisotropic Gram normalization.
Measuring all families at the same nonzero momentum would compress their
energy ratios toward one, not enlarge the rest-mass hierarchy.

## 3. A stronger exclusion: every affine law built from this feature fails

The result is broader than failure of the coefficients \(a,b,J\). Let

\[
Y=uI+vF,\qquad u,v\in\mathbb C,
\]

with canonical kinetic normalization and one common Higgs factor. Since
\(F\) is Hermitian, its complex eigenvalues under this map are
\(z_i=u+vf_i\), and the masses are proportional to \(|z_i|\).
Writing \(\phi=(1+\sqrt5)/2\), the source eigenvalues obey

\[
\frac{f_3-f_1}{f_2-f_1}=\phi^2,
\qquad z_3-\phi^2z_2+\phi z_1=0.
\]

The triangle inequality, allowing every assignment of the largest modulus
to the three coefficients, gives the necessary bound

\[
\boxed{m_{\max}\leq\phi^2m_{\mathrm{middle}}+\phi m_{\min}.}
\]

To see that it covers every assignment, the absolute coefficients in the
linear relation are \(1,\phi,\phi^2\). Placing the largest mass on the
smallest coefficient gives the weakest constraint; assigning the larger
remaining coefficient to the middle mass maximizes its right-hand side.
Any other placement gives a smaller bound. Complex phases do not evade
the triangle inequality.

For the retained broad comparison windows this requires

\[
\frac{m_{\max}}{m_{\min}}
\leq207\phi^2+\phi=543.55106966\ldots,
\]

while the observed lower bound is 3470. Therefore **the entire complex
affine class fails**, including arbitrary identity shifts, overall scales,
sign changes and phases. This is an exact exclusion of that class, not a
numerical search.

It also explains why merely changing the same protected-core reference
weights is insufficient. If \(w_q\) is the weight per configuration with
core occupation \(q\), direct counting gives

\[
\mathbb E[n_{\rm core}]=10w_1+w_2,\qquad
\mathbb E[n_{\rm outer}]=9w_0+2w_1.
\]

For the same core and band, its compressed mean is always

\[
(9w_0+2w_1)I+(8w_1+w_2-9w_0)S_{\rm core},
\]

where \(S_{\rm core}\) is the compressed core projector. The present
\(F\) is a nonscalar affine function of that same projector, so this whole
family of mean-based linear prescriptions remains in the excluded class.
This argument does not cover arbitrary changes to the carrier, core type,
observables or interaction law.

## 4. Nonlinearity can amplify masses without fixing their pattern

The three previously tested Hubbard-function textures are fixed
prescriptions, not fits. All amplify the source eigenvalue differences,
but all give the wrong relative spacing:

| Prescription | Middle/light | Heavy/light | Heavy/middle |
| --- | ---: | ---: | ---: |
| Linear \(F\) | 1.1903 | 1.4982 | 1.2587 |
| Hubbard function with \(U=af_i\) | 204.975 | 65,743.94 | 320.74 |
| Hubbard function with \(U=bf_i\) | 353,640.05 | \(3.5494\times10^{11}\) | \(1.0037\times10^6\) |
| Hubbard function with \(U=Jf_i\) | 8,075.31 | \(1.3683\times10^8\) | 16,943.78 |
| Observed charged leptons, approximately | 206.77 | 3,477.23 | 16.82 |

The near-middle ratio in the \(a\) control is accompanied by a heavy mass
about nineteen times too large after light-mass normalization. A common
scale cannot correct it. A useful scale-independent shape statistic is

\[
Q=\frac{\log(m_3/m_2)}{\log(m_2/m_1)}.
\]

The controls give \(Q\) near 1.082–1.084, while the broad empirical windows
give approximately 0.529–0.531. Certified rational logarithm bounds make
these intervals disjoint. The wrong *shape*, not inadequate arithmetic
precision, is the issue.

The mechanism is understandable: the weak-coupling one-dimensional
Hubbard threshold contains an exponential in \(-1/U\), so using
\(U=sf_i\) amplifies the reciprocal eigenvalue spacings in a specific
pattern. The source theorem is for the half-filled one-dimensional
Hubbard chain, as discussed by
[Jeckelmann, Gebhard and Essler](https://arxiv.org/html/cond-mat/9911281).
Its use as a numerical texture does not derive an interacting 3+1D chiral
mass. These fixed tests do not exclude all nonlinear laws or all values
of a continuously varied coupling.

An arbitrary function \(f(F)\) is also not a prediction by itself. Since
\(F\) has three distinct eigenvalues, quadratic polynomial interpolation
can assign any three target values to them. Those polynomial coefficients
would encode the desired answer unless an independently specified source
law determines them.

## 5. Small normalization corrections cannot hide the discrepancy

Let \(r=\kappa_2(M)<3/2\) be the current mass-matrix condition number. For
new positive constant kinetic forms \(Z_L,Z_R\), canonical normalization
gives, in the usual Dirac convention,

\[
M_c=Z_L^{-1/2}MZ_R^{-1/2},\qquad
\kappa_2(M_c)\leq r\sqrt{\kappa_2(Z_L)\kappa_2(Z_R)}.
\]

The transpose in an all-left convention does not change this norm bound.
Reaching a heavy/light ratio of at least 3470 requires the product of
kinetic condition numbers to exceed \(5.35\times10^6\). If both
\(\|Z-I\|_2\leq\epsilon<1\), it requires
\(\epsilon>6937/6943\), approximately 0.999136.

Similarly, a constant additive mass correction with
\(\|\delta M\|_2=\delta\,s_{\min}(M)\), \(\delta<1\), obeys

\[
\kappa_2(M+\delta M)\leq\frac{r+\delta}{1-\delta}.
\]

The same target requires \(\delta>6937/6942\), approximately 0.999280.
These are necessary bounds, not sufficient constructions or claims about
the size of every possible quantum correction. They rule out explanations
that leave the current action almost unchanged in the specified norms.
Strong dynamics, additional fields and momentum-dependent self-energies
remain outside this argument.

## 6. Family structure must be recovered from more than one mean matrix

Every separate function of \(F\) shares its eigenvectors. Consequently
\(Y_uY_u^\dagger\) and \(Y_dY_d^\dagger\) commute when the Yukawa
matrices are functions of this one Hermitian feature. For nondegenerate
masses and canonical common kinetic terms their CKM matrix is trivial
up to phases and permutations. Choosing more elaborate eigenvalue
functions cannot supply physical quark mixing.

The new audit also computes the full connected occupation covariance

\[
C_{ij}=\mathbb E[n_in_j]-\mathbb E[n_i]\mathbb E[n_j]
\]

from all 66 configurations. Its full rank is 11, with the single constant
null direction required by fixed total load. Its compression to the actual
three-dimensional band is positive definite, has spread below \(5/4\),
and **commutes exactly with \(F\)**. Thus replacing the mean by this
equal-reference covariance does not produce either a large linear
hierarchy or an independent family basis. This is a computed negative
control, not a claim that every dynamical response is determined by these
two matrices.

There is a further symmetry restriction on using only this stationary
reference. Its native stabilizer contains a half-turn represented on the
family band by \(R=\operatorname{diag}(-1,+1,-1)\). Any covariant
endomorphism built solely from that reference and rotation-invariant
carrier data must commute with \(R\). For the same family action on both
fermion factors and a family-neutral Higgs, the allowed matrices have a
common two-dimensional block and an unmixed one-dimensional block.
They can support at most two-family mixing, with no three-family CKM
CP violation. This restriction is conditional on those inputs; it is not
a symmetry theorem for every full repair context.

### The existing update contains information that the stationary feature loses

The actual
[reference preparation](../../model/native_modular_reference/local_preparation.py)
already uses two distinct completion caps, \((0,1,5,7)\) and
\((0,1,2,8)\). The half-turn fixes the first and moves the second. Its full
update context therefore has less symmetry than its stationary reference.
The new
[native response audit](../../model/common_source/hierarchy_diagnosis/native_family_dynamics.py)
uses the existing finite selector and endpoint repair, without adding a
mass-dependent parameter.

It enumerates the reduced outcomes on every one of the 66 binary inputs,
including the selector's terminal timeout. Each accepted outcome retains
its seam, coin, current, output state and exact probability. The timeout
retains unchanged work with a separate terminal flag. The complete
one-proposal mean response is

\[
T=\sum_e w_e E_e+\tau I,
\qquad \sum_e w_e+\tau=1,
\]

where \(E_e\) averages the two endpoint occupations. Its diagonal native
reference is exactly stationary. With orthonormal band frame
\(U=B/\sqrt{10+2\sqrt5}\), define the compressed response
\(K=U^TTU\). The calculation finds

\[
[F,K]\ne0,\qquad [R,K]\ne0.
\]

All three off-diagonal links of \(K\) are nonzero. Since \(F\) has simple
eigenvalues, its spectral projectors \(P_i\) are polynomials in \(F\).
Each \(P_iKP_j/K_{ij}\) is a matrix unit. Therefore \(F,K\) generate
the full three-by-three matrix algebra, an exact statement rather than a
numerical rank guess.

This is a concrete positive result: **the existing update context supplies
independent family structure omitted by the mean-based mass ansatz**.
It does not select Yukawa matrices or predict mixing angles. These
retained matrices are real and do not by themselves derive a physical
CP phase. The caps and preparation protocol also remain declared
operational inputs, not a selected autonomous matter dynamics.

There is a second important limitation: the band is not closed under
this update. For \(V=(I-UU^T)TU\), the exact calculation gives

\[
V\ne0,\qquad U^TT^2U-K^2=V^TV\ne0.
\]

Thus iterating only the compressed matrix would lose the complement's
response. The displayed identity is an algebraic marginal comparison;
it is **not** a protocol that restarts after a terminal timeout. A continuing
reduction must retain the complementary modes or their derived memory,
along with the actual event and timeout rules. The parameter of this
calculation is one completion proposal, not physical time, modular time
or a calibrated count of elapsed native attempts. No eigenvalue of this
relaxing response has been identified with a Lorentzian particle mass.

## 7. What the next mass calculation must actually derive

The next target is the charged response of a specified continuing,
interacting repair construction. It must determine the effective mass
coupling and kinetic response together. Schematically its charged inverse
response has a kinetic part, a Higgs-dependent part and any interaction
self-energy; physical rest energies follow its singularities in the
observer's calibrated clock. An occupation statistic alone does not
supply that complete response.

Three checks make this a concrete research target:

1. Start from the existing update's independent family response and retain
   its complementary modes or derived memory. Derive the actual coupling to charged
   left and conjugated-right fields, rather than choose an eigenvalue
   function to match masses.
2. Compute the same law's kinetic normalization, background and charged
   response, with its state, regulator and clock assumptions explicit.
   Compare action parameters with observable poles or thresholds in the
   appropriate domain; native relaxation rates need their own bridge.
3. Test both independent charged-lepton ratios and nontrivial family
   mixing before moving to costly QCD or making physical-species claims.

This diagnosis does not require proving that one microscopic architecture
is uniquely possible. It requires a specified architecture to produce its
mass prediction without using the target masses to choose its couplings.
It also leaves intact the separate obligation to attach a complete chiral
field history to the observers' reconstructed spacetime; the current
frozen principal-symbol comparison is not that full attachment.
