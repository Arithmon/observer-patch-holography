# Charged finite response and an even observer apparatus

The executed four-mode ring has actual charged fermionic response. The CAR
operators connect its half-filled ground state to nonzero positive spectral
residues in the one-particle and three-particle sectors. These are excitations
of the declared continuing law, with an occupation charge and fermionic
algebra. Their identification as electrons, quarks or any other Standard
Model species requires additional results.

`response.report()` evaluates all finite Lehmann residues of the same
Hamiltonian used in `mass_response`. It consumes its native-reference
coefficient producer and the existing CAR implementation. No empirical mass
or coupling is read. The assignment of those coefficients to the continuing
SSH Hamiltonian remains the previously declared architecture choice.

## The exact finite identities

Write the unique two-particle ground as $|g\rangle$, with energy $E_g$.
For every three-particle eigenstate $|a\rangle$ and one-particle eigenstate
$|b\rangle$, define

```math
 R^+_{a,ij}=\langle g|c_i|a\rangle\langle a|c_j^\dagger|g\rangle,
 \qquad
 R^-_{b,ij}=\langle g|c_j^\dagger|b\rangle\langle b|c_i|g\rangle.
```

The addition and removal costs are $\Delta_a^+=E_a-E_g$ and
$\Delta_b^-=E_b-E_g$. Their signed spectral frequencies are respectively
$+\Delta_a^+$ and $-\Delta_b^-$. The retarded resolvent is

```math
 G^R_{ij}(z)=\sum_a\frac{R^+_{a,ij}}{z-\Delta_a^+}
           +\sum_b\frac{R^-_{b,ij}}{z+\Delta_b^-},\qquad \Im z>0.
```

Every residue is positive semidefinite, and CAR gives the full matrix sum
rule $\sum_aR_a^++\sum_bR_b^-=I_4$. The implementation also checks the
first-frequency-moment identity
$\sum_a\Delta_a^+R_a^+-\sum_b\Delta_b^-R_b^-=\langle\{[c_i,H],c_j^\dagger\}\rangle$.
These sum rules test weights as well as eigenvalues.

Centered densities $\delta n_i=n_i-\langle n_i\rangle$ have a separate
Lehmann decomposition entirely within the two-particle sector. Its integrated
residue is the density covariance. The sum over all densities has zero
response because total number is conserved. A density transition is neutral;
its gap is not a charged fermion mass merely because it is an eigenvalue.

All finite eigenspaces and dark residues are retained. Degenerate neutral
eigenspaces are summed to remove arbitrary eigenbasis choices. The displayed
arrays are floating evaluations of exact finite identities, with a stated
numerical visibility diagnostic; they are not outward arithmetic enclosures.

## Reading a charged response without odd charge coherence

A charged CAR operator alone is not an invariant apparatus interaction on a
fixed-total-charge history. The existing interferometer already supplies two
identical cells and two arms,

```math
 |r\rangle=|g_2\rangle_A|g_2\rangle_B,\qquad
 |q\rangle=|g_3\rangle_A|g_1\rangle_B.
```

Both arms have total charge four. The globally graded transfer
$T_{ij}=C_{Ai}^\dagger C_{Bj}$ is even and conserves total charge. With cell A
in the low bits, $C_A=I_B\otimes c_A$ and
$C_B=c_B\otimes(-1)^{N_A}$. Omitting that parity factor gives the wrong
operator on other charge sectors.

The transfer matrix element
$m_{ij}=\langle q|T_{ij}|r\rangle$ factorizes into a CAR addition amplitude
and a CAR removal amplitude, and is nonzero at the ground-arm frequency
$\omega_{\rm pair}=E_3+E_1-2E_2$ measured by the existing return protocol.
Thus the measured gap has a charged-CAR channel: it is not just an arbitrary
relaxation eigenvalue. For the unread two-arm state at model parameter $s$,

```math
 \langle X_{ij}\rangle=\Re(m_{ij}e^{+i\omega_{\rm pair}s}),\qquad
 \langle Y_{ij}\rangle=\Im(m_{ij}e^{+i\omega_{\rm pair}s}),
```

where $X=T+T^\dagger$ and $Y=-i(T-T^\dagger)$. A common chemical
potential in the convention $H(\mu)=H(0)+\mu N$ moves the addition cost by
$+\mu$ and the removal cost by $-\mu$, while this pair frequency remains
unchanged. Reading a single charged cost therefore also requires a chemical
potential convention; the symmetric half-gap uses particle-hole symmetry.
`transfer_visibility` admits that half-gap only when the supplied chemical
potential is zero and the two ground costs agree within its numerical check.
At nonzero chemical potential the pair frequency is retained and the
single-gap field is null. Pair interference alone cannot infer both costs.

## Inverse response and coordinate covariance

[kinetic_response.py](kinetic_response.py) forms the inverse and its analytic
frequency derivative from the complete same-law spectral matrix:

```math
 K(z)=G(z)^{-1},\qquad
 \partial_zG(z)=-\sum_r\frac{R_r}{(z-\omega_r)^2},\qquad
 \partial_zK(z)=-K(z)(\partial_zG(z))K(z).
```

The probe frequency remains strictly in the upper half-plane. For an
admitted packet the supplied residues are finite Hermitian positive matrices
whose actual sum agrees with the canonical CAR identity within the declared
numerical tolerance. The interface rejects missing/doubled nonzero residues
and inconsistent signed frequencies; it returns the actual transformed sum
and the expected coordinate Gram separately. This consistency validation
alone does not prove a supplied packet's Hamiltonian provenance. For an
invertible operator-coordinate change $c'=Ac$, $G'=AGA^\dagger$ and
$K'=A^{-\dagger}KA^{-1}$; the derivative transforms by the same congruence.
The free response has the pencil $K=zZ-M$. Its roots are the generalized
eigenvalues of $(M,Z)$, invariant under that simultaneous change. Ordinary
eigenvalues of $M$ alone can change and are not a repaired spectrum.

These identities use the full finite Hamiltonian response. The derivative
at a complex probe frequency is not asserted to be a positive physical
kinetic metric. In the interacting ring it depends on frequency, so a
static whitening alone does not remove the self-energy. Independent
finite-difference and free-action tests check the implementation. No detector
tomography, spatial derivative, physical chiral field, observer-clock matching
or measured mass follows from this finite calculation. The current
[particle receipt](results.json) includes this diagnostic separately from
the actually executed continuing record experiment.

## Local charge and apparatus scope

The invariant symmetry here is **global occupation U(1)**. This transfer is
not an established locally electromagnetic-gauge-invariant field. Gauging
the symmetry would require a link or dressing between distinct sites. The
new apparatus couples a declared pair of cell-local modes; emergent spacetime
locality of that apparatus pair has not been derived.

`transfer_read(state, site_a=0, site_b=0, quadrature="X")` executes a native
source instrument on the caller's existing work. It returns the existing
`history.Readout` type, with two unnormalized outputs. Bits 0–3 are cell A,
bits 4–7 are cell B, higher bits are spectators, and reference columns remain
untouched. The corresponding effects are $(I\pm X)/2$, or $(I\pm Y)/2$.
Since these observables also have a zero eigenspace, this is a binary POVM.

The actual compiler emits the Jordan–Wigner transfer rotation, uses the
existing diagonal pointer read for $n_A-n_B$, and undoes the rotation.
Native pulses implement $e^{i\theta P}$: the basis words are
$B_X=e^{+i\pi Y/4}$ and $B_Y=e^{-i\pi X/4}$, with
$B^\dagger(n_A-n_B)B=X$ or $Y$. Both positive Kraus maps preserve total
fermion charge at the completed-instrument boundary. This charge statement
does not apply to each intermediate processor pulse. Pointer preparation,
control and physical service accuracy remain declared apparatus services.

The instrument can be appended separately to every retained continuing
branch; earlier scratch aborts stay terminal. It adds no work preparation,
does not reset either cell, retains both new outcomes, and permits the
existing clocks and their untouched purification columns to continue.
A final expectation after other readouts includes their backreaction. It is
not an independent spectral pole measurement. The Lehmann matrix element and
the appended detector are separate, explicit checks of the same operator.

## What has and has not been identified

The charge labels are changes of encoded occupation, not established
electromagnetic charges. The CAR algebra supplies fermionic statistics in
this encoded finite model. The ring is spinless; no Lorentz spin-1/2
representation, three-dimensional interacting dispersion, infinite-volume
stable charged pole, gauge dressing, quark confinement, flavor identification
or mass in laboratory units follows from these finite residues. A discrete
finite spectral pole by itself is insufficient: every finite closed system
has discrete spectral lines.

The finite charged threshold can therefore be named a charged fermionic
excitation of this model. The calculation does not name a familiar physical
species or establish an interacting Lorentzian particle mass.

There is also a concrete obstruction to naming the existing elementary
charged poles a color triplet. In each of the N=1 and N=3 sectors the quartic
term vanishes and the four energies are
$-(a+b),-(a-b),a-b,a+b$, plus the sector's constant chemical-potential
shift. The exact ratios $e^b=8/5>1$ and $e^{a-b}=125/64>1$ prove their
strict ordering. Every eigenspace is one dimensional and the corresponding
CAR residue has rank one. The unitary commutant preserving this Hamiltonian
and charge in either charged sector is consequently $U(1)^4$. Relabeling
these poles cannot produce a nontrivial degenerate color triplet. This result
concerns the single ring's elementary charged sectors; it excludes no
composite, gauge-constrained, infinite-volume or refined construction.
