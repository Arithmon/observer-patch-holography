"""Exact rational checks of the analytic global threshold error budget.

Uses conservative rational boxes rather than a numerical optimizer or the
kinematic producer. The proof explains why these bounds control a global
three-dimensional minimum. No sampled scan is called a global certificate.
"""

from fractions import Fraction as F


def error_budget():
    m_lo, m_hi = F(500000), F(600000)
    d_lo, d_hi = F(5, 10**58), F(6, 10**58)
    lo, floor, soft_max = F(9*10**16), F(10**13), F(200)
    hi = F(10**20) + soft_max  # arbitrary incoming soft direction: |K n+s u|
    # All radicals are replaced by outward rational bounds.
    root3d_hi = F(43, 10**30)  # sqrt(3 d) < 4.3e-29
    h = (20*d_hi)**2 * hi**5 / 600
    t = (m_hi**4/floor**3 + 2*m_hi*m_hi*d_hi*hi + d_hi*d_hi*hi**5)/2
    photon = h + d_hi*d_hi*hi**5/2
    universal_gap_error = 2*(h+t) + photon
    photon_gap_error = 2*m_hi**4/lo**3 + photon
    # |D0'(Q)| for each branch, where D0 = 2 epsilon0.
    universal_lip = 2*m_hi*m_hi/lo**2 + F(9, 8)*d_hi*hi**2 + m_hi*root3d_hi
    photon_lip = 2*m_hi*m_hi/lo**2 + F(3, 2)*d_hi*hi**2
    radial_loss = 2*d_hi*hi**2 + d_hi*soft_max**2
    soft_energy_error = d_hi*soft_max**3
    universal = (universal_gap_error + soft_max*(universal_lip+radial_loss))/2 + soft_energy_error
    photon_only = (photon_gap_error + soft_max*(photon_lip+radial_loss))/2 + soft_energy_error
    # Energy-coordinate replacement K -> Omega(K), without changing domains.
    energy_relabel_error = max(universal_lip, photon_lip)*d_hi*hi**3/2
    small_leg_gap = m_lo*m_lo/(2*floor+m_hi) - 2*h
    feasible_gap = 2*m_hi*m_hi/lo + m_hi*root3d_hi*hi + 2*(h+t) + photon
    # A hard Cherenkov channel in the photon-only variant. This is a
    # strict existence witness, not an exact onset or an interaction rate.
    electron, emitted, recoil = F(3*10**17), F(2*10**17), F(10**17)
    vc_leading = d_lo*emitted**3/2 - m_hi*m_hi*emitted/(2*electron*recoil)
    vc_error = m_hi**4/(2*electron**3) + m_hi**4/(2*recoil**3)
    vc_error += (20*d_hi)**2*emitted**5/600 + d_hi*d_hi*emitted**5/2
    return dict(universal=universal, photon_only=photon_only,
                energy_relabel_error=energy_relabel_error,
                small_leg_gap=small_leg_gap, feasible_gap=feasible_gap,
                h=h, t=t, photon=photon,
                universal_lip=universal_lip, photon_lip=photon_lip,
                root3d_hi=root3d_hi, d_hi=d_hi, d_lo=d_lo,
                photon_only_cherenkov_margin=vc_leading-vc_error,
                asymmetric_leg_lower=m_lo/(root3d_hi*hi),
                symmetric_leg_lower=lo/2,
                threshold_lower=m_lo*m_lo/hi,
                threshold_upper=m_hi*m_hi/lo+d_hi*hi**3/4+max(universal,photon_only),
                compact_domain_parameter=20*d_hi*(5*hi)**2,
                compact_minimizer_ratio=(1+(soft_max+2*m_hi)/lo)**2,
                soft_domain_margin=F(10**17)-soft_max-lo,
                small_parameter=m_hi*m_hi/floor**2+d_hi*hi**2)


def verify_budget():
    b = error_budget()
    conditions = {
        "radical enclosure": b["root3d_hi"]**2 > 3*b["d_hi"],
        "small leg excluded": b["small_leg_gap"] > b["feasible_gap"],
        "leading asymmetric leg in expansion domain": b["asymmetric_leg_lower"] > 10**13,
        "leading symmetric leg in expansion domain": b["symmetric_leg_lower"] > 10**13,
        "positive lower threshold bracket": b["threshold_lower"] > max(b["universal"],b["photon_only"]),
        "upper bracket within soft domain": b["threshold_upper"] < 200,
        "compact domain inside cutoff": b["compact_domain_parameter"] < 1,
        "minimizer inside compact domain": b["compact_minimizer_ratio"] < F(19,5),
        "remaining pair momentum in proof range": b["soft_domain_margin"] > 0,
        "sqrt expansion valid": b["small_parameter"] < F(1, 2),
        "universal global threshold": b["universal"] + b["energy_relabel_error"] < F(1, 10**13),
        "photon-only global threshold": b["photon_only"] + b["energy_relabel_error"] < F(1, 10**13),
        "photon-only emission witness": b["photon_only_cherenkov_margin"] > F(7, 10**7),
    }
    for name, passed in conditions.items():
        if not passed:
            raise ArithmeticError(name)
    return b


if __name__ == "__main__":
    for name, value in verify_budget().items():
        print(f"{name}: {value} (approximately {float(value):.9g})")
