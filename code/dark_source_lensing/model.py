"""Closed-form producer; no natural observations or fitted parameters."""
from fractions import Fraction as F
from math import isqrt
import mpmath


def build():
    ctx = mpmath.mp.clone()
    ctx.dps = 70
    def mp(q):
        return ctx.mpf(q.numerator)/q.denominator
    def out(x):
        return ctx.nstr(x, 55)
    rows = []
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        lo, hi = u/(2*(1+u)), u/(1+2*u)
        for ratio in (2, 10, 100):
            for label, k in (('stiff_radial', lo), ('interior', (lo+hi)/2), ('zero_radial', hi)):
                theta = ctx.acos(ctx.mpf(ratio)**(mp(u)-1))
                sweep = 2*theta/((1-mp(u))*ctx.sqrt(1-2*mp(k)))
                alpha = sweep-2*theta
                inverse = (1-(2*theta/((1-mp(u))*(alpha+2*theta)))**2)/2
                rows.append(dict(u=str(u), R_over_r0=ratio, branch=label, kappa=str(k),
                                 radial_ratio=str(u/k-1-2*u),
                                 tangential_ratio=str(u*u*(1-2*k)/(2*k)),
                                 sweep=out(sweep), alpha=out(alpha), inferred_kappa=out(inverse)))
    cold = []
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        for eta in (F(0), F(1, 100), F(1, 2)):
            cold.append(dict(u=str(u), eta=str(eta),
                             lower=str(max(u/(2*(1+u)), u/(1+2*u+eta))),
                             upper=str(min(F(1, 2), u/(1+2*u-eta)))))
    kinetic = []
    def root_interval(q):
        scale = 10**18
        a = isqrt(q.numerator*scale*scale//q.denominator)
        return F(a, scale), F(a+1, scale)
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        for cap in (u, 2*u, F(1, 2)):
            lo = u*(1+u)/(1+2*u+2*u*u+cap)
            hi = u/(1+2*u)
            a0, a1 = root_interval(1/(1-2*lo))
            b0, b1 = root_interval(1+2*u)
            lower = (a0-(1-u))/(b1-(1-u))
            upper = min(F(1), (a1-(1-u))/(b0-(1-u)))
            kinetic.append(dict(u=str(u), speed_squared_cap=str(cap), lower_kappa=str(lo),
                                upper_kappa=str(hi), mass_fractional_width=str((hi-lo)/hi),
                                sqrt_low_kappa=[str(a0), str(a1)],
                                sqrt_high_kappa=[str(b0), str(b1)],
                                bending_ratio_lower=str(lower), bending_ratio_upper=str(upper)))
    controls = []
    # Q(x)-1 = sum c_p/x^p. All entries are transform controls, not matter sources.
    for terms in (((1, F(1, 3)),), ((2, F(2, 5)),), ((3, F(3, 7)),),
                  ((4, F(5, 9)),), ((1, F(1, 3)), (3, F(-1, 7)), (4, F(2, 9)))):
        for b in (F(1, 2), F(1), F(3)):
            angle = sum(mp(c)*ctx.sqrt(ctx.pi)*ctx.gamma(ctx.mpf(p+1)/2)/
                        ctx.gamma(ctx.mpf(p+2)/2)*mp(b)**(-p) for p, c in terms)
            controls.append(dict(terms=[[p, str(c)] for p, c in terms], x=str(b),
                                 alpha=out(angle), Q_minus_one=out(sum(mp(c)*mp(b)**(-p) for p, c in terms))))
    bounds = dict(density_lower=str(F(1, 27)-F(1, 160)),
                  radial_upper='1/2000',
                  tangential_upper=str(F(3, 196)+F(9, 1280)+(F(1, 160)+F(1, 2000))/196),
                  compactness_upper=str((2+F(1, 1280))/100),
                  bending_gap_per_epsilon=str(F(81, 3512320)),
                  bending_gap_at_epsilon_001=str(F(81, 351232000)))
    integral = 248*ctx.log(2)-ctx.mpf(1719)/10
    integral_lower = F(int(ctx.floor(integral*10**18)), 10**18)
    integral_upper = integral_lower+F(1, 10**18)
    enclosure = dict(integral_lower=str(integral_lower), integral_upper=str(integral_upper),
                     rows=[dict(epsilon=str(e), lower=str(e*integral_lower/5),
                                upper=str(e*integral_upper*F(539, 2000)))
                           for e in (F(1, 100), F(1, 10**6))])
    samples = []
    for r in (F(1), F(5, 4), F(3, 2), F(7, 4), F(2)):
        e = F(1, 100)
        m0 = e*r**3/(1+r**3)
        R0 = 3*e*r*r/(1+r**3)**2
        u = m0/(r-2*m0)
        du = (r*R0-m0)/(r-2*m0)**2
        h = (r-1)**4*(2-r)**4
        dh = 4*(r-1)**3*(2-r)**3*(3-2*r)
        C = (1+2*u)/r
        dC = 2*du/r-(1+2*u)/r**2
        for t in (F(-1, 10), F(0), F(1, 10)):
            R = R0+t*e*dh
            P = -t*e*h*C
            dP = -t*e*(dh*C+h*dC)
            T = (r*dP+(R+P)*u)/2
            samples.append(dict(r=str(r), t=str(t), m=str(m0+t*e*h),
                                u=str(u), density_scaled=str(R), radial_scaled=str(P),
                                tangential_scaled=str(T)))
    return dict(annulus=rows, cold_bounds=cold, kinetic_bounds=kinetic, abel=controls, global_bounds=bounds,
                global_samples=samples, global_ray_enclosure=enclosure,
                interpretation=dict(result='macroscopic_source_underdetermined',
                                    physical_promotion=False, natural_data_used=False,
                                    full_issue_751_closed=False,
                                    numeric_controls='high_precision_not_interval_proofs'))
