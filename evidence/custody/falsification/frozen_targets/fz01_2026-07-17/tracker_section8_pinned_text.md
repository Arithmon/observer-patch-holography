## 8. Freeze registry shortlist

The suite holds zero frozen prospective predictions (all four audits
agree). These are the rows to hash-pin in one registry artifact BEFORE the
next releases (PDG 2026, DESI DR3, O5, Hyper-K first light, LiteBIRD).
Rows marked (cond.) freeze together with their named condition.

| Row | Value | Compare against | Trigger |
|---|---|---|---|
| n_s = 1 - P/48 | 0.9660214 | SPT-3G, ACT DR6+, CMB-S4 combined tilt | next combined-likelihood release |
| Future W/Z physical readout; current chart diagnostics 80.330 / 91.119 / 0.22279 are ineligible | no physical value frozen | named W/Z release and convention | CL-03 full observable map plus theory covariance; otherwise `NOT_EVALUABLE` |
| m_H(m_t) criticality relation | 125.72 GeV at measured m_t | PDG m_H | PDG update |
| alpha_s(M_Z) | 0.11834 | world average | PDG update |
| Lambda_QCD(3) | 334.8 [319,350] MeV | FLAG | FLAG update |
| Ratio scoreboard (4 ratios) | CL-06 values | PDG/CODATA | PDG update |
| Glueball 0++ | 1.41-1.61 GeV | BESIII PWA f0(1500)/f0(1710) | next BESIII PWA |
| Capacity Lambda (EW branch) | Lambda l_P^2 = 2.668e-122 | Planck/DESI Lambda display | DESI DR3 |
| (w0, wa) (cond. fixed-N theorem) | (-1, 0) | DESI DR3 + SNe | DR3 release |
| Sigma m_nu (cond. rank lane) | 0.0588 eV, NO | DESI DR3 + CMB | DR3 release |
| m_bb (cond. rank lane) | 1.5-3.7 meV | LEGEND-1000/nEXO | first results |
| Proton decay | no signal (tau_gauge = infinity) | Hyper-K | first exposure years |
| Future ringdown transition (canonical integer-k Kerr family) | no value frozen; historical alpha = 4 branch is void because it implies k = e | O4/O5 stacked ringdowns | new source-derived transition and full-likelihood registration before data opening |
| Birefringence (cond. emission theorem) | 0.37501 deg | LiteBIRD, PR4 reanalyses | theorem, then release |
| Dark-siren H0 | 67.4 +/- 1 | LVK O5 dark sirens | O5 catalog |
| IR filter (cond. source emission) | (q_IR, t_IR), ell_IR ~ 32 | Planck low-ell likelihood | emission theorem |
| UHE coefficients (after VN-08) | real emission values | IceCube/Auger/LHAASO | VN-08 completion |
| Neutron-bottle stance | bottle correct | next beam remeasurement | publication |
| Exactly-3 block | N_nu = 3, no steriles/exotics | LHC, BEST resolution, N_eff | ongoing |
| No particle DM (cond. response branch) | all direct searches null | LZ full exposure, XLZD | ongoing |

Registry mechanics: one JSON + markdown artifact in
`reverse-engineering-reality/`, SHA-256 of inputs and code, named
acceptance band per row, scored on release, hits and kills recorded in
place. Pattern to copy: the boundary-scale candidate registry and the
`NO_UHE_DATA_USE` receipt.

---

