# A bounded local repair family with a common faithful reference

One small change to the C/E repair law has a positive finite consequence:
an endpoint placement coin biased by two supplied local port weights preserves
a nonuniform faithful reference. Every complete C/E branch is covariant under
that reference's modular action. This constructs a compatible repair/reference
pair; it does not select the weights or execute a physical clock.

The [contract](campaign.json) freezes this finite alternative class separately
from the unchanged C/E campaigns. The [implementation](local_bias.py),
[tests](test_local_bias.py) and [exact receipt](results.json) retain the old
fair-coin law as a control. No predecessor receipt is edited.

## Compatibility first, then the update

For the existing twelve-port binary total-two carrier, independently choose
positive integer port weights `theta_i`. The compatibility target is

```text
pi(x) = product_i theta_i^x_i / Z,    sum_i x_i=2, x_i in{0,1}.
```

On a requested seam `(a,b)`, keep the native endpoint floor/ceiling map and
the complete C or E raw-record partition. Replace the fair coin by
`P(coin=1)=theta_a/(theta_a+theta_b)` and
`P(coin=0)=theta_b/(theta_a+theta_b)`. With one endpoint load, these place it
at a or b respectively. Empty and fully occupied binary endpoints stay fixed.
The bias reads no hidden source load or spectator: its two weights and seam
are declared control data.

The classical transition is precisely conditional resampling of pi at fixed
outside configuration and endpoint total. Thus the existing
[`FiniteConditionalRepair`](../../Lean/Thermodynamics/FiniteConditionalRepair.lean)
row-sum, detailed-balance, stationarity and fixed-mixture theorems apply to the
actual endpoint kernel. Exact enumeration checks every premise. Positive
fixed seam weights connect all 66 configurations; self-loops give aperiodicity
and the unique classical stationary law. No full quantum uniqueness claim is
needed for this construction.

Each quantum record map is `K_b=sqrt(p_coin) P_b`, with P_b the existing native
partial injection. Completeness is checked on both H78 and H66, including
zero raw slots. Pi is faithful only on H66. For decoded current J,

```text
pi(output)/pi(input) = (theta_b/theta_a)^J,
pi^it K_b pi^-it = (theta_b/theta_a)^(it J) K_b.
```

The multiplier is constant within each record block. It cancels in the CP
map, proving modular covariance of the complete instrument with an untouched
reference of any dimension. The sum of its outputs fixes pi. C and E retain
different coherence behavior despite their equal population transition.
This does not assert that the quantum channel is a conditional expectation,
that passive repair equals modular evolution, or that a modular phase service
has been compiled for the new reference.

The existing [regional factorization test](../native_modular_reference/common_caps.py)
applies: all bare two-port number-preserving algebras are invariant under this
product reference's modular action. The additional result here is an actual
changed local repair instrument preserving the same nonuniform state.

## One exact nontrivial member and one excluded fixed target

Set `theta_0=theta_1=25` and the ten other weights to8. This control chooses
the earlier target's ratio `w1/w0=25/8`; it is not a source-selected prediction
or a fit to a measured physical quantity. Crossing coins are
25/33 and8/33. The reference's core-occupation masses for0,1,2 loads are
576/1501,800/1501,125/1501. These differ from the earlier transported
three-charge target. Fixed nonuniform seam selection, identity holds and
tagged C/E mixtures preserve this state; record-dependent feedback requires
a separate joint source/controller analysis.

For the earlier target, individual configuration weights are
`w0=2/225,w1=1/36,w2=2/45`. Its non-product ratio is
`w1^2/(w0*w2)=125/64`. On any seam crossing a fixed two-port core, moving
outward with the second load inside the core requires forward/reverse rate
ratio5/8; the identical endpoint move with an exterior spectator requires
8/25. A positive endpoint-only reversible law cannot satisfy both. The
receipt gives actual native seams and two spectator configurations for all
66 cores.

This excludes that fixed target only for the stated reversible local class.
It does not exclude irreversible stationary currents, a larger local
neighborhood, an interaction, correlated controller state, or co-moving
reference metadata. The unchanged-law moving-core construction has a
different state space and retains its separate positive result. In particular,
changing a local law and extending its classical metadata are distinct options.

The [complete-branch algebra result](../completion/defect_algebra.md) also
survives this law adjustment. Every admitted seam and coin has strictly positive
probability, and each new Kraus map is a nonzero scalar multiple of the same
native partial injection. Each C/E endpoint therefore retains the full branch
star algebra on H66 and H78 and its scalar commutant. Tagged mixtures and
identity holds do not create a nontrivial Hermitian bare-source charge commuting
with every active branch. This says nothing stronger about nonreducing invariant
subspaces, arbitrary expectation-fixed observables, history/controller charges,
metastability, active recovery or refined carriers. Bias alone does not supply
exact bare-source particle protection.

## Frozen finite scope and the full CP-00 gate

The finite class permits primitive integer port weights1..32 with ratio at
most4, fixed positive seam weights1..16, rational C/E mixture and hold
probabilities with denominator at most256, and repair probability at least
1/16. Coin probabilities lie in[1/5,4/5]; every seam has probability at least
1/480 before holds. Complete records may be exported through the existing
classical record interface. No active recovery, adaptive source-conditioned
selector, new quantum register, phase primitive or refined carrier is included
in this certificate.

The weights are up to60 supplied bits before the primitive-gcd constraint;
the seam table is up to120 bits. Neither is derived from native symmetry or
observer records. Their inhomogeneity and the corresponding mathematical
`-sum_i log(theta_i) n_i` generator do not supply a physical potential, action
unit, defect mass or stress. A record store or controller using these data has
its own service and storage costs; this packet does not realize that hardware.

This completes a bounded source-law/reference decision under
[#1135](https://github.com/FloatingPragma/observer-patch-holography/issues/1135),
but cannot close the full issue. Four required portions remain unspecified:

1. The carrier/refinement class and how permitted memory, controls and services
   scale while describing a continuing physical regime.
2. The admissible interpretation and shared clock/length/energy dictionary,
   including which calibrations and external inputs are permitted.
3. The complete field, particle, neutrino and gravity target menu with exposure,
   uncertainty and covariance conventions.
4. The program-wide computational budgets, success thresholds and stopping
   coverage across that class.

The contract requires these omissions to remain explicit and rejects a
whole-program freeze or CP-00 acceptance flag. No physical comparison or
integrated gate is added. The original q=3 scalar history retains its
unexecuted status; this alternative does not retire it.

```sh
python -m model.local_repair_family.local_bias --check
python -m pytest -q -W error model/local_repair_family
```
