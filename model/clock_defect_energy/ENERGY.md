# Repair activity, rest energy and independent inertia

The mass-energy relation supplies a useful constraint on a repair action. It
does not identify that action from a count alone. This diagnostic tests one
simple activity ansatz on the same native carrier as the
[continuation](CONTINUATION.md) and [localized-pair](LOCALIZED.md) calculations,
then checks the kinetic coefficients that a relativistic candidate must have.
The [receipt](energy_results.json) retains exact rational results and source
bindings; [tests](test_energy.py) use independent Kraus effects, graph counts
and symbolic action derivatives.

## The tested native cost

Take all thirty seams uniformly and retain every C/E coin and record outcome.
The cost per attempted event is either one attempt, absolute current or
squared current. Binary current is zero or plus/minus one, so the last two
costs agree. This is a complete-instrument expectation, without conditioning
on a successful move. The localized target is independently defined as two
occupied neighboring ports, with thirty possible configurations.

| Background reference | Background movement per attempt | Adjacent-pair movement | Signed excess |
| --- | --- | --- | --- |
| Uniform binary66 | `5/33` | `2/15` | `-1/55` |
| Transported target, core `{0,1}` or `{1,5}` | `341/2250` | `2/15` | `-41/2250` |
| Transported target, core `{0,8}` or `{5,8}` | `409/2700` | `2/15` | `-49/2700` |

The last two rows use the exact three-charge preparation target transported
through the predecessor's actual histories. These are faithful full-shell
reference extensions, not normalized posterior corners or the unknown actual
finite-preparation state. Conditioning each target on adjacency defines an
offline ensemble; no adjacency measurement is added to the observer. Uniform
seam selection is common to all rows; correlated feedback has its own joint
reference and cannot inherit this subtraction unchanged.

Attempt cost has zero excess. Movement cost has negative excess for both C
and E and all these backgrounds. Multiplying it by a positive action scale
cannot give positive rest energy for this candidate. Changing the sign after
seeing the result or taking an absolute value would change the tested law.
This excludes this particular positive excess-activity ansatz; it neither
identifies a physical vacuum nor excludes other local maintenance actions.

Two further exact controls constrain the interpretation:

- A two-event native cycle returns the source load to its initial state while
  accumulating two unsigned movements. Its record archive changes. Unsigned
  activity therefore cannot be the increment of a source-only state energy
  on every branch. A history action or source-plus-environment energy remains
  possible.
- Two source configurations have identical attempt, absolute-current and
  squared-current expectations, but their declared reference weights have
  ratio `8/25`. No function of this instantaneous cost vector equals
  `-log(sigma)` up to a scalar shift on the whole carrier. The modular
  generator of that chosen target is not thereby identified as physical rest
  energy.

## What E₀ = mc² requires

Suppose the native process eventually supplies the proper-time action

```text
S = -E0 integral(d_tau)
  = integral[-E0 + E0*v^2/(2*c^2) + E0*v^4/(8*c^4) + ...] dt.
```

For independently obtained coefficients `L=-E0+A*v^2+B*v^4+...`, the inertial
mass measured by the quadratic velocity response is `m_kin=2*A`. The two
checks here are

```text
E0 / (2*A*c^2) = 1,
8*B*c^4 / E0 = 1.
```

The proper-time benchmark passes. Doubling its independently supplied inertia
fails the first ratio, giving `1/2`. Keeping its rest energy and inertia but
setting the quartic coefficient to zero passes `E0=m_kin*c^2` and fails the
next relativistic coefficient. These are necessary checks, not a derivation
of full dispersion, interactions or gravitational coupling. A common change
of energy, length and time units leaves both ratios unchanged.

A concrete repair-action hypothesis would be `S=-a*N`, with positive action
`a` per event and `N=r*integral(d_tau)` for a constant proper-time rate `r`.
It then gives `E0=a*r`. **Action per event times events per time has energy
units.** Energy per event times events per time instead has power units.
The sign, reference subtraction and proper-time attachment are substantive
parts of the hypothesis, not unit conversions. One may instead derive a
Hamiltonian directly; these diagnostics do not require a count action.

The current native pair has no independently constructed velocity response,
proper-time action or conserved gravitational stress. Neither its exit rate
nor a modular gap is used as a mass measurement. Assigning a wave speed to
SI `c` cannot repair a failed dimensionless ratio or supply those missing
observables. The next construction must obtain rest action, moving response
and stress from the same admitted process and calibrated clock dictionary.

## Replay

```sh
python -m model.clock_defect_energy.energy --check
python -m pytest -q -W error model/clock_defect_energy/test_energy.py
```
