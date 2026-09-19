# FZ-10 frozen target: the Koide conditional tau window

Registered 2026-07-28, before any comparison data beyond the embedded
reference existed. Source commit
`bbc2cebedf190af351fba1eb6c8887fb340eb6bb` in reverse-engineering-reality
carries the producing certificate
`code/particles/leptons/koide_balance_comparison_certificate.py` and the
frozen artifact whose byte copy sits beside this statement.

## The frozen stance

Under two declared premises, the balanced positive-chamber circulant
(`rho/a = 1/sqrt2`, supplied conditionally by the finite tracial
Gelfand-Naimark-Segal construction) and the mass ordering
`m_tau > m_mu`, the measured electron and muon masses determine the tau
mass through one quadratic. The outward-rounded 100-decimal-digit
enclosure is

    m_tau in [1776.968991, 1776.969063] MeV,   central 1776.969027 MeV.

The premise ancestry is declared: the balance condition was first
abstracted from the measured charged triple (Koide, 1981-1983). The
stance is therefore a conditional postdiction whose confirmation weight
is limited by that ancestry; its kill direction carries no such limit.

## Decision policy (frozen with the target)

Reference observable: the charged-tau pole-mass world average as
published by the Particle Data Group, or a dedicated single-experiment
measurement with a stated standard uncertainty `sigma_avg`, examined
after this registration.

- FAIL (kills the balanced-circulant premise): the central value differs
  from 1776.969027 MeV by more than `3 sigma_avg`.
- COMPATIBLE: the central value differs by at most `2 sigma_avg` and
  `sigma_avg <= 0.045 MeV` (half the 2026 world-average uncertainty).
- INCONCLUSIVE: every other outcome; no verdict is issued at
  insufficient precision.

The window half-width (36 eV) is negligible against any foreseeable
`sigma_avg` and the window midpoint is the comparison point; the full
window is recorded above. A kill is published with the same prominence
as a confirmation, after the audit / re-audit / repair ladder of the
kill-condition protocol.

## Custody

- `koide_balance_comparison_frozen_2026-07-28.json`: byte copy of the
  producing artifact at the source commit;
  sha256 `09efbad8de813790d10671b425daef69429831c58daa212f4be131c61f653898`.
- OpenTimestamps calendar stamps are submitted for this statement, the
  artifact copy, and the registration manifest; the Bitcoin upgrade
  follows on the standard `ots upgrade` pass.
