# Recovery with complete native records

The complete classical tape and final native load register recover every
initial load population exactly. They do not recover an arbitrary unknown
quantum state on the same carrier. The missing information is coherence
between different committed histories, not the classical currents needed
to reverse the load motion.

This conclusion has a constructive boundary. An explicit qubit code is
exactly correctable after one C or E event. Without intermediate quantum
recovery, suitable finite seam words rule out every nontrivial exact code.
The two record laws also differ under bounded classical control: E can
completely measure the initial carrier in seven binary or eight full-carrier
events; C cannot become entanglement breaking under any bounded
record-controlled seam policy on the existing graph. These are decisions
about the native candidate's passive and actively corrected memory, not a
programme-wide obstruction to quantum dynamics.

## The recorded process and allowed recovery

Use the existing [native update](native.py) and
[source C/E instruments](source_discrimination/DERIVATION.md). Inputs are
the 78 nonnegative integer total-two configurations, or their invariant
66-state binary subset. Write the chosen input set as `X`, with dimension
`D=|X|`. An event selects an oriented native seam `ab`, draws a fresh fair
coin `c`, and balances its endpoint loads by floor/ceiling placement. With
`d=x_a-x_b`, its current is

```text
j_c(d) = [d-(2c-1)(abs(d) mod 2)]/2
       = ceil(d/2)  if c=0,
         floor(d/2) if c=1.
```

C commits `(ab,c,j)`; E commits `(ab,c,d)`. E's record also determines `j`
by this identity. Both are **classical** records: distinct outcomes are
orthogonal density-matrix blocks. The coherent write isometry before
commitment is a different channel and can be inverted if its complete
quantum output remains accessible.

A seam policy may depend on earlier committed records and classical
randomness. It cannot inspect hidden loads or perform intervening quantum
encodings, resets or corrections. A final decoder may act on the entire
load register and tape, preserve the tape, use blank ancillas, and return a
state on the original carrier. Untouched references are allowed throughout.
This decoder access is stronger than the existing observer controller's
record-only access. Its mathematical existence does not supply a physical
load readout or coherent service calibration.

## Every history is a partial permutation followed by a projection

For a complete admitted history `h`, let `X_h` be the initial load states
compatible with it, and `P_h=sum_{x in X_h}|x><x|`. The branch probability
factor `w_h` is the product of its seam-selection probabilities and fair
coin factors. It is the same for every compatible input; the actual branch
probability is `w_h Tr(P_h rho)`, not necessarily `w_h`.

Every retained current determines a translation of the load vector.
Reversing those translations gives the initial vector uniquely from the
terminal vector and `h`. Therefore the history Kraus operator has the form

```text
K_h = sqrt(w_h) U_h P_h,          U_h* U_h = P_h,
sum_h w_h P_h = I,
Gamma(rho) = sum_h |h><h| tensor K_h rho K_h*.                     (1)
```

`U_h` is the injective basis map on `X_h`; it may be completed to a unitary
permutation of the carrier. A decoder controlled by the classical tape can
invert it while leaving `h` unchanged. Its output is exactly

```text
sum_h w_h |h><h| tensor P_h rho P_h.                              (2)
```

Thus all initial diagonal matrix units are recovered, while only the
coherences within each history fiber are restored. Discarding the tape
after (2) gives the Schur channel

```text
rho_xy -> H_xy rho_xy,
H_xy = sum_{h: x,y in X_h} w_h,           H_xx=1.                  (3)
```

These identities hold on all matrix units and therefore after tensoring
with any reference. They apply to the actual two-event observer as well as
the separate schedules below. For E, the existing `e`/`f` seam selector uses the
decoded **current** sign, not the contrast sign.

## Optimal recovery of the complete unknown input

Define entanglement fidelity for the maximally mixed input by the overlap
with `|Omega_D>=D^(-1/2)sum_x|x,x>` after applying a final recovery to one
half. Put `d_h=rank(P_h)`. Among all final CPTP decoders,

```text
F_opt = D^(-2) sum_h w_h d_h^2
      = D^(-2) sum_{x,y} H_xy.                                   (4)
```

The inverse in (2), followed by discarding the tape, attains (4). To prove
optimality, any decoder restricted to a classical branch is a CPTP map
with Kraus operators `R_ha`. In a basis `i` of `P_h`, Cauchy--Schwarz and
trace preservation give

```text
sum_a |Tr(R_ha U_h P_h)|^2
 = sum_a |sum_{i=1}^{d_h} <i|R_ha U_h|i>|^2
 <= d_h sum_{i,a} |<i|R_ha U_h|i>|^2
 <= d_h sum_{i,a} ||R_ha U_h|i>||^2 = d_h^2.
```

Multiplication by `w_h/D^2` and summation proves the upper bound. Mixing
different classical branches cannot improve it: the restriction of any
global decoder to each orthogonal input flag is itself CPTP. This is an
all-decoder optimum for this input convention, not only the performance
of one inverse circuit or a statement about an already encoded qubit.

There is also an exact two-state witness. For distinct `x,y`, let
`|psi_+/- >=(|x> +/- |y>)/sqrt(2)`. Orthogonal history flags and injectivity
inside each branch imply

```text
dist(Gamma(psi_+), Gamma(psi_-)) = H_xy,
dist(rho,sigma) = ||rho-sigma||_1/2.                              (5)
```

The inputs have distance one. Contractivity and the triangle inequality
therefore force every decoder to have error at least `(1-H_xy)/2` on at
least one of them. Here `psi_+/-` in a channel argument denotes its density
matrix; no outcome is postselected.

For the actual first seam `(0,1)`, take `x={0,2}` and `y={1,2}`.
Their contrasts are opposite. Their currents differ for each coin, so
both C and E have `H_xy=0` already after this first event. The two phase
states give identical complete recorded outputs, forcing recovery error
at least one half on one of them. Later processing cannot reverse that
loss. This concerns unencoded states; (8) supplies a distinct encoded
countercontrol.

For the existing two-event current-controlled C observer, and for E with
the same current-sign selector applied to its decoded current, exhaustive
history fibers give the following instances of (4):

| Input carrier | C optimal fidelity | E optimal fidelity | Full-input correctable block sizes, both laws |
| --- | --- | --- | --- |
| Binary, `D=66` | `5543/8712` | `769/1452` | `46,10,9,1` |
| Full, `D=78` | `1301/2028` | `3329/6084` | `56,10,9,1,1,1` |

The remaining quantum blocks are expected: this controller touches only
its declared seams and halts. These values do not assume that it explores
the entire graph.

## Full-input observable algebra versus encoded subspaces

On the entire input space, the maximal exactly correctable observable
algebra is

```text
A_correctable = {A: [A,P_h]=0 for every admitted h}.               (6)
```

Commutation is sufficient by (2) and `sum_h w_h P_h=I`. For necessity,
suppose a decoder `R` exactly corrects a unital algebra. For any element
`A`, the two Schwarz inequalities for `Gamma* R*` are equalities because
`A`, `A*` and `A*A` are all corrected. In a Stinespring representation `V`
of `Gamma`, this implies `(R*(A) tensor I_env)V=VA`, and the corresponding
adjoint identity. The environment projection onto history `h` commutes
with every output operator, so sandwiching it by these identities gives
`[A,V*(I tensor |h><h|)V]=0`. The sandwich is `w_h P_h`, proving (6).

Since the `P_h` are diagonal, (6) is the direct sum of full matrix algebras
on classes of basis states having identical history membership. Equivalently,
`x,y` share a class exactly when `H_xy=1`. With full-support random seam
choice on a connected graph, already one event makes this algebra precisely
the diagonal load algebra: any nonzero fixed-total difference `x-y` has
unequal values at the endpoints of some edge, and at least one coin then
gives different records.

This statement does **not** rule out an encoded quantum subspace. For a
code projector `Q`, the flagged Kraus operators in (1) obey the
[Knill--Laflamme condition, Theorem 3.2](https://arxiv.org/pdf/quant-ph/9604034)
exactly when

```text
Q P_h Q = alpha_h Q                 for every admitted h.        (7)
```

Cross-history products vanish because the output flags are orthogonal.
The condition allows noncommuting `Q`: each record need only have a
probability independent of the logical state. When `alpha_h>0`, the
normalized restriction `U_h P_h Q/sqrt(alpha_h)` is an isometry of the
code, so a record-conditioned inverse recovers it. Zero-probability
branches do not require postselection or a normalized update.

For a concrete one-event code on the retained graph, take

```text
|0_L> = (|{0,2}>+|{1,3}>)/sqrt(2),
|1_L> = (|{0,3}>+|{1,2}>)/sqrt(2).                               (8)
```

All four occupied pairs are nonedges. The logical supports are disjoint,
their one-port occupation means agree, and every native-edge product
`N_a N_b` has expectation zero. A binary C or E event effect is a linear
combination of `I,N_a,N_b,N_a N_b`; its off-diagonal logical matrix element
vanishes and its two diagonal logical elements agree. Thus (7) holds for
every seam, coin and record, for either law. One unknown logical qubit,
including its reference entanglement, can be recovered exactly after one
event. The native full-state indices are `(75,63)` and `(74,64)` for the
two logical supports; binary-only indices are `(64,53)` and `(63,54)`.

Repeated active correction and re-encoding could protect this code, but
that is an additional intervention between native events. The passive
history bounds below do not cover it, nor do they derive its physical
controls, cost or compatibility with the uniform geometry/clock protocol.

The separate [active encoded attachment](../common_source/ENCODED_MEMORY.md)
compiles that recovery and re-encoding inside a continuing common-source
history. Its matched passive control and paid interventions have their own
evidence; the final-recovery protocol and receipt here retain their stated
scope.

## Independent random seam attempts: exact all-integer rates

This section uses independent seam choices with probabilities `p_e>0` on
a connected graph and independent fair coins. It does not replace the
existing `e`/`f` seam controller with that schedule. For any equal-total integer
inputs put `z=x-y`. The ceiling/floor formula gives the same-record
probability on edge `ab` exactly:

```text
s_C(ab,z) = 1-min(abs(z_a-z_b),2)/2,
s_E(ab,z) = 1[z_a=z_b],
q_L(z) = sum_edges p_e s_L(e,z).                                 (9)
```

Two equal integers give equal currents for both coins; neighboring
integers do so for exactly one coin; integers separated by at least two
never do. This includes negative odd contrasts. On a matching-record
branch both configurations translate by the same current, so `z` remains
fixed. Conditional survival at the next event is again `q_L(z)`, giving

```text
H_xy(t) = q_L(x-y)^t,
F_opt(t) = D^(-2) sum_{x,y} q_L(x-y)^t.                           (10)
```

Here `t` counts seam attempts. Connectedness and equal total imply
`q_L(z)<1` for every nonzero admitted `z`. Thus `F_opt(t)` tends to `1/D`,
while all initial populations remain
perfectly recoverable. Retaining the tape preserves the classical input
even as the optimal full-state entanglement fidelity approaches its
complete-dephasing value.

On the existing uniform thirty-seam graph, the largest off-diagonal rates
are `q_C,max=5/6` and `q_E,max=7/10`, on both carriers. For binary pairs
with difference `e_u-e_v`, the C loss numerator is ten and the E numerator
is `2(10-1[uv is an edge])`, with denominator sixty. Binary disjoint pairs
have numerators `20-2(e_pp+e_mm)` and
`2[20-2(e_pp+e_mm)-e_pm]`: `e_pp,e_mm<=1` count edges within the two
positive/negative pairs, and `e_pm<=4` counts edges between them. These
give lower bounds sixteen and twenty-four. The extra full-carrier
differences are `2e_u-2e_v`, with both numerators at least eighteen, and
`2e_u-e_v-e_w`. For the latter, if `a` counts edges from `u` to `{v,w}` and
`b=1[vw is an edge]`, the numerators are `20-a-2b>=16` for C and
`2(15-a-2b)>=22` for E. These exhaust all total-two pair types and prove
the maxima directly from (9).

There is a uniform reference bound, not only a scalar fidelity statement.
Let `Delta` completely dephase the **initial** load basis. The channel
`Gamma_t Delta` is entanglement breaking. For every source/reference state,

```text
dist((Gamma_t tensor id)(rho), (Gamma_t Delta tensor id)(rho))
 <= min{1, (D-1) q_max^t/2}.                                    (11)
```

Write `rho=sum_xy |x><y| tensor rho_xy`, and `pi_x=Tr(rho_xx)`.
Positivity gives `||rho_xy||_1<=sqrt(pi_x pi_y)`. Each retained off-diagonal
matrix unit has trace norm `H_xy(t)`, so the triangle inequality bounds
twice (11) by `q_max^t sum_{x!=y}sqrt(pi_x pi_y) <= q_max^t(D-1)`.
The reference dimension is unrestricted. The comparison channel depends
on `t`; it includes the growing tape and does not erase its old records.

Closeness to an entanglement-breaking channel is not exact membership.
For an instrument (1), the normalized Choi state is a direct sum over
`h` of pure vectors with Schmidt rank `d_h`. Thus the complete instrument
is entanglement breaking **if and only if every positive-weight fiber has
rank one**. A rank-two fiber gives an entangled conditional Choi block;
local projection onto its flag proves that the unconditioned output is
entangled too. At every finite iid horizon, a repeated single-seam history
has a positive-weight fiber containing distinct spectator configurations.
Consequently neither retained C nor retained E is entanglement breaking
at any finite iid horizon, on either carrier, despite (11).

## Finite schedules and the classical-control boundary

For a fixed E seam word, matching histories are equivalent to `z` being
constant on every connected component of the graph of visited edges.
Necessity follows from the contrast record. Sufficiency follows because
matching currents preserve `z`. Alternatively, every measured contrast
can be translated back to its initial value using the preceding currents.
The history fiber therefore depends on this visited graph, not on the
coin word.

For full total-two inputs this identifies every basis state exactly if
and only if there is at most one visited component of size one or two.
Any nonzero constant on a component of size at least three would already
require more than two positive or negative loads. Two small components
permit a zero-total difference: sizes one and two balance with constants
two and minus one, while equal sizes balance with constants one and minus
one. On the binary carrier the constants have magnitude at most one, so
identification instead requires at most one singleton component and at
most one component of size two.

The full carrier thus needs at least eight distinct edges in a fixed E
word on twelve vertices; the binary carrier needs at least seven. The
following native forests attain those minima:

```text
Full:   (0,1),(0,5), (2,3),(2,4), (6,8),(8,9), (7,10),(10,11)
Binary:       (1,5), (2,3),(2,4), (6,8),(8,9), (7,10),(10,11)
```

Their component sizes are respectively `(3,3,3,3)` and `(1,2,3,3,3)`.
In the first case the maximum possible component count is four; in the
second it is five, proving the edge lower bounds. Repeated edges cannot
improve them. Each complete E history then has rank one, so the resulting
retained instrument is exactly entanglement breaking. These are minima
among fixed seam words, not a claim of optimality over adaptive policies.

C can learn an edge's contrast by using its current record twice. After
the first coin `c`, repeats with the same coin have zero current. The first
opposite coin gives `j_*=(2c-1)p`, where `p=abs(d) mod 2` was the initial
seam parity. Hence the original contrast is `d=2j_1+j_*`.

Double each edge in one of the displayed forests. All successive coin
pairs being opposite has input-independent probability `2^(-m)`, where
`m=8` or `7`. On those branches every initial basis point has its own
rank-one history effect. For a putative exact code of dimension greater
than one, (7) applied to each such `|x><x|` forces `Q|x>=0`: a rank-at-most-one
operator cannot be a nonzero scalar times `Q`. Every basis point occurs,
giving a contradiction. Thus the 16-event full or 14-event binary C word
admits no nontrivial exact code under final recovery. This is a sufficient
construction, not a minimum-length C theorem, and it is stronger than
showing loss of one particular unencoded phase pair.

C nevertheless cannot be entanglement breaking after any bounded
record-only adaptive schedule on this graph. Choose nonadjacent `u,v` and
a third site `w`, and compare `x={u,w}`, `y={v,w}`. Their difference is
`e_u-e_v`, whose endpoint difference on every admitted edge has magnitude
at most one. Whatever seam a common history selects, at least one fair
coin retains equal currents. For any horizon bounded by `T`, their shared
history probability is at least `2^(-T)`; halts may be padded by identities.
Some history fiber has rank at least two, proving non-entanglement-breaking
by the Choi argument above. Classical randomization of the seam policy
does not change this lower bound. Thus absence of exact codes after the
doubled forest does not imply entanglement breaking.

There is also an unbounded stopping policy. On each forest edge, repeat
until the coin first differs from its initial value, then advance. All
`m` contrasts are eventually learned with probability one. The stopping
time is `T=m+sum_{i=1}^m G_i`, with independent geometric variables of
success probability one half supported on positive integers:

```text
Pr(T=t) = binomial(t-m-1,m-1) / 2^(t-m),   t>=2m,
E[T] = 3m = 24 (full), 21 (binary).                              (12)
```

At stopping the complete instrument is entanglement breaking. There is
no deterministic event or tape-capacity bound, so this does not evade the
bounded-policy obstruction. These probing schedules are different from
the existing two-event observer and from the iid calibration.

## Source correspondence and candidate decision

[record_recovery.py](record_recovery.py) constructs the existing observer's
history fibers, inverts admitted currents and calculates the exact recovery
kernel and iid rates. The [schedule certificate](record_recovery_schedule.py)
decodes initial-load constraints from the actual tape, checks the two
forests and retained ambiguity witnesses, and reports incomplete stopped
histories without renormalizing them. [record_recovery_source.py](record_recovery_source.py)
uses the existing source C/E inverse repair words: it copies an already
classical archive into blank scratch, restores the source and clears the
scratch while retaining the archive. The single-seam check implements (2),
not recovery of the erased cross-record coherence. Only seam `(0,1)` is
executed there; the other native seams are mathematical port relabellings.

The relevant candidate decision for #1076/#1026 is precise. Native current
records suffice to reverse classical load transport and preserve declared
classical memory. They do not by themselves supply an exactly reversible
quantum response on the whole bare-load algebra. Encoded one-event recovery
is possible and gives a concrete active-control escape; coherent uncommitted
records, protected extra registers, or restricted exploration are other
changes of premises. None is silently counted as already present in the
passive native law. Their compatibility with the same geometry, clock,
physical source access and complete reversible current remains a separate
model requirement. No additional finite theorem is used as a substitute
for that common-model attachment or empirical law selection.
