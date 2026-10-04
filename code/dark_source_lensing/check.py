"""Independent Einstein, geodesic and Abel replay. Imports no producer."""
from fractions import Fraction as F
from functools import lru_cache
import mpmath
import sympy as s
from . import global_ray
from .format import decimal, equal, keys, need, rational


@lru_cache(maxsize=1)
def symbolic_identities():
    r = s.Symbol('r', positive=True)
    u, k, e, t = s.symbols('u k e t', real=True)
    m, v = s.Function('m')(r), s.Function('v')(r)
    # Reconstruct angular Einstein component directly from m and nu'=v.
    angular = (-(s.diff(m, r)*r-m)*(1+r*v)/r +
               (1-2*m/r)*r*r*(v/r+v*v+s.diff(v, r)))/2
    radial = (r-2*m)*v-m/r
    tov = (r*s.diff(radial, r)+(s.diff(m, r)+radial)*r*v)/2
    need(s.cancel(angular-tov) == 0, 'independent angular Einstein/TOV identity')
    flat = {m:k*r, v:u/r}
    need(s.simplify(angular.subs(flat).doit()-u*u*(1-2*k)/2) == 0, 'flat angular stress')
    need(s.simplify(radial.subs(flat)-u+(1+2*u)*k) == 0, 'flat radial stress')
    need(s.diff(r*r*u/(1-u), r) == 2*r*u/(1-u), 'stable flat circular orbits')
    h = (r-1)**4*(2-r)**4
    for order in range(4):
        for endpoint in (1, 2):
            need(s.diff(h, r, order).subs(r, endpoint) == 0, 'bump junction derivative')
    need(s.expand(s.diff(h, r)-4*(r-1)**3*(2-r)**3*(3-2*r)) == 0, 'bump derivative')
    mass = e*r**3/(1+r**3)
    speed = mass/(r-2*mass)
    specific_l = r*r*speed/(1-speed)
    # Positive for 0<e<=1/100: numerator 4+r^3-6e*r^2 >0.
    target = e*r**3*(4+r**3-6*e*r**2)/(1+r**3-3*e*r**2)**2
    need(s.cancel(s.diff(specific_l, r)-target) == 0, 'global orbit stability polynomial')
    # Substitute after evaluating derivatives to avoid partial substitution artifacts.
    R = s.diff(mass+t*e*h, r)
    P = ((r-2*(mass+t*e*h))*speed/r-(mass+t*e*h)/r).cancel()
    V = speed/r
    T = (-(R*r-(mass+t*e*h))*(1+r*V)/r+
         (1-2*(mass+t*e*h)/r)*r*r*(V/r+V*V+s.diff(V, r)))/2
    return r, e, t, mass+t*e*h, speed, R, P, T


def close(ctx, actual, expected, label):
    value = decimal(ctx, actual)
    scale = abs(expected) if expected else ctx.mpf(1)
    need(abs(value-expected) <= ctx.mpf('1e-48')*scale, label)


def verify(packet):
    keys(packet, 'annulus cold_bounds kinetic_bounds abel global_bounds global_samples global_ray_enclosure winding_controls interpretation')
    equal(packet['interpretation'], dict(result='macroscopic_source_underdetermined',
          physical_promotion=False, natural_data_used=False, full_issue_751_closed=False,
          numeric_controls='high_precision_not_interval_proofs'), 'scientific boundary')
    ctx = mpmath.mp.clone()
    ctx.dps = 75
    def mp(q):
        return ctx.mpf(q.numerator)/q.denominator
    rows = packet['annulus']
    need(type(rows) is list and len(rows) == 27, 'complete 27-ray census')
    cursor = 0
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        low, high = u/(2+2*u), u/(1+2*u)
        for ratio in (2, 10, 100):
            for name, k in zip(('stiff_radial', 'interior', 'zero_radial'), (low, (low+high)/2, high)):
                row = rows[cursor]
                cursor += 1
                keys(row, 'u R_over_r0 branch kappa radial_ratio tangential_ratio sweep alpha inferred_kappa')
                equal([row['u'], row['R_over_r0'], row['branch'], row['kappa']],
                      [str(u), ratio, name, str(k)], 'fixed ray input and branch')
                # Reconstruct components before dividing by density.
                P, T = (1-2*k)*u-k, (1-2*k)*u*u/2
                need(0 <= P <= k and 0 <= T <= k, 'nonnegative dominant-energy source')
                equal(row['radial_ratio'], str(P/k), 'radial stress')
                equal(row['tangential_ratio'], str(T/k), 'tangential stress')
                # Independent nonsingular geodesic integral: r/r0=exp(z^2).
                U, K = mp(u), mp(k)
                q = 2*(1-U)
                def integrand(z):
                    if not z:
                        return 2/ctx.sqrt(q*(1-2*K))
                    return 2*z/ctx.sqrt((1-2*K)*ctx.expm1(q*z*z))
                sweep = 2*ctx.quad(integrand, [0, ctx.sqrt(ctx.log(ratio))])
                need(0 < sweep < 2*ctx.pi, 'retained annulus ray has no full winding')
                endpoint = ctx.asin(ctx.exp((U-1)*ctx.log(ratio)))
                alpha = sweep+2*endpoint-ctx.pi
                close(ctx, row['sweep'], sweep, 'integrated null sweep')
                close(ctx, row['alpha'], alpha, 'finite-endpoint bending')
                # Recover B directly from sweep and endpoint geometry.
                recovered_B = ((1-U)*sweep/(ctx.pi-2*endpoint))**2
                inferred = (1-1/recovered_B)/2
                close(ctx, row['inferred_kappa'], inferred, 'joint inverse')
                close(ctx, row['inferred_kappa'], K, 'original mass recovered')
    verify_winding(packet['winding_controls'])
    cold = packet['cold_bounds']
    need(type(cold) is list and len(cold) == 9, 'complete coldness census')
    cursor = 0
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        for eta in (F(0), F(1, 100), F(1, 2)):
            row = cold[cursor]
            cursor += 1
            keys(row, 'u eta lower upper')
            equal([row['u'], row['eta']], [str(u), str(eta)], 'fixed coldness input')
            lo, hi = rational(row['lower']), rational(row['upper'])
            # Solve the two pressure inequalities by their active equality.
            need((1+2*u+eta)*lo == u and (1+2*u-eta)*hi == u,
                 'sharp coldness endpoints')
            need(u/(2+2*u) <= lo <= hi < F(1, 2), 'coldness inside dominant energy')
    verify_kinetic(packet['kinetic_bounds'])
    r, e, t, mass, speed, R, P, T = symbolic_identities()
    samples = packet['global_samples']
    need(type(samples) is list and len(samples) == 15, 'complete global sample census')
    cursor = 0
    for radius in (F(1), F(5, 4), F(3, 2), F(7, 4), F(2)):
        for perturbation in (F(-1, 10), F(0), F(1, 10)):
            row = samples[cursor]
            cursor += 1
            keys(row, 'r t m u density_scaled radial_scaled tangential_scaled')
            equal([row['r'], row['t']], [str(radius), str(perturbation)], 'global sample inputs')
            sub = {r:s.Rational(radius.numerator, radius.denominator), e:s.Rational(1, 100),
                   t:s.Rational(perturbation.numerator, perturbation.denominator)}
            for field, expression in zip(('m', 'u', 'density_scaled', 'radial_scaled', 'tangential_scaled'),
                                         (mass, speed, R, P, T)):
                value = s.cancel(expression.subs(sub))
                equal(row[field], str(value), 'Einstein component '+field)
    b = packet['global_bounds']
    keys(b, 'density_lower radial_upper tangential_upper compactness_upper bending_gap_per_epsilon bending_gap_at_epsilon_001')
    values = {key:rational(value) for key, value in b.items()}
    density = F(1, 27)-F(1, 160)
    radial = F(1, 2000)
    tangential = F(3, 196)+F(9, 1280)+(F(1, 160)+radial)/196
    need(F(50, 49) < F(21, 20) and F(16, 100)+F(50, 49) < F(6, 5), 'coefficient bounds')
    need(F(21, 51200) < radial, 'radial bound from bump')
    need((F(1, 16)*F(21, 20)+F(1, 256)*F(6, 5))/10 == F(9, 1280), 'radial derivative bound')
    need(0 < radial < density and 0 < tangential < density, 'global strict dominant energy')
    compactness = (2+F(1, 1280))/100
    need(compactness < 1, 'no horizon')
    # Reconstruct every factor in the restricted positive ray integral.
    gap = 4*F(1, 2)*F(1, 10)*(F(3, 16)**4)*F(1, 2)/(F(7, 4)**3)
    expected = dict(density_lower=density, radial_upper=radial, tangential_upper=tangential,
                    compactness_upper=compactness, bending_gap_per_epsilon=gap,
                    bending_gap_at_epsilon_001=gap/100)
    equal({k:str(v) for k, v in values.items()}, {k:str(v) for k, v in expected.items()},
          'exact global certificate including nonzero gap')
    verify_ray_enclosure(packet['global_ray_enclosure'])
    global_ray.verify_enclosures(packet['global_ray_enclosure']['rows'])
    verify_abel(ctx, packet['abel'])


def verify_winding(rows):
    need(type(rows) is list and len(rows) == 3, 'complete winding catalogue')
    for n, row in enumerate(rows):
        keys(row, 'turns u r0_over_R endpoint_angle_over_pi reduced_sweep_over_pi sweep_over_pi alpha_over_pi kappa radial_ratio tangential_ratio')
        equal([row['turns'], row['u'], row['r0_over_R'], row['endpoint_angle_over_pi']],
              [n, '1/5', '2^(-5/4)', '1/6'], 'fixed winding geometry')
        k = rational(row['kappa'])
        sweep = rational(row['sweep_over_pi'])
        need(0 < k < F(1, 2) and sweep == 2*n+1, 'unwrapped ray sweep')
        # theta=pi/3; squaring the null-geodesic law is exact with positive roots.
        need(sweep*sweep*(1-F(1, 5))**2*(1-2*k) == F(4, 9),
             'winding branch null equation')
        equal(row['reduced_sweep_over_pi'], str(sweep-2*n), 'same endpoint azimuth modulo 2pi')
        equal(row['alpha_over_pi'], str(sweep+F(1, 3)-1), 'unwrapped finite bending')
        P = (1-2*k)/5-k
        T = (1-2*k)/50
        need(abs(P) <= k and abs(T) <= k, 'winding source obeys dominant energy')
        equal(row['radial_ratio'], str(P/k), 'winding radial stress')
        equal(row['tangential_ratio'], str(T/k), 'winding tangential stress')
    # On the nonnegative-pressure catalogue u<=1/5, even the maximal sweep
    # is <2pi: (1+2u)/(1-u)^2 <=35/16 <4.
    need(F(35, 16) < 4, 'weak-catalogue branch bound')


def verify_ray_enclosure(row):
    keys(row, 'integral_lower integral_upper rows')
    lo, hi = rational(row['integral_lower']), rational(row['integral_upper'])
    r = s.Symbol('r', positive=True)
    integral = s.integrate(s.expand((r-1)**4*(2-r)**4)/r**3, (r, 1, 2))
    need(s.expand(integral-248*s.log(2)+s.Rational(1719, 10)) == 0,
         'independently integrated ray weight')
    # Exact atanh series and positive geometric remainder, not mpmath.log.
    partial = sum(F(2, (2*j+1)*3**(2*j+1)) for j in range(32))
    tail = F(9, 4*65*3**65)
    lower, upper = 248*partial-F(1719, 10), 248*(partial+tail)-F(1719, 10)
    need(0 < lo <= lower < upper <= hi and hi-lo <= F(1, 10**18),
         'tight outward rational log-integral enclosure')
    # Upper kernel factors: lapse ratio, turning denominator, spatial derivative.
    lapse_factor = F(49, 48)
    need(lapse_factor*lapse_factor/4 < 1-F(25, 36), 'turning denominator at most 6/5')
    c = (2+F(1, 1280))/100
    need((1-c)**3*F(121, 100) >= 1, 'spatial derivative factor at most 11/10')
    need(lapse_factor*F(6, 5)*F(11, 10) == F(539, 400), 'full kernel upper factor')
    need(type(row['rows']) is list and len(row['rows']) == 2, 'two source strengths')
    for item, e in zip(row['rows'], (F(1, 100), F(1, 10**6))):
        keys(item, 'epsilon lower upper')
        equal(item['epsilon'], str(e), 'enclosed source strength')
        equal(item['lower'], str(e*lo/5), 'global bending lower enclosure')
        equal(item['upper'], str(e*hi*F(539, 2000)), 'global bending upper enclosure')


def verify_kinetic(rows):
    need(type(rows) is list and len(rows) == 9, 'complete kinetic census')
    cursor = 0
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        for cap in (u, 2*u, F(1, 2)):
            row = rows[cursor]
            cursor += 1
            keys(row, 'u speed_squared_cap lower_kappa upper_kappa mass_fractional_width sqrt_low_kappa sqrt_high_kappa bending_ratio_lower bending_ratio_upper')
            equal([row['u'], row['speed_squared_cap']], [str(u), str(cap)], 'fixed kinetic inputs')
            lo, hi = rational(row['lower_kappa']), rational(row['upper_kappa'])
            # Independent endpoint equalities: trace cap below, zero radial pressure above.
            need(u*(1+u)-(1+2*u+2*u*u)*lo == cap*lo, 'kinetic lower trace equality')
            need(u-(1+2*u)*hi == 0, 'kinetic upper radial equality')
            need(u/(2+2*u) <= lo <= hi < F(1, 2), 'kinetic bound inside dominant energy')
            width = (cap-u)/(1+2*u+2*u*u+cap)
            equal(row['mass_fractional_width'], str(width), 'exact kinetic mass width')
            need((hi-lo)/hi == width, 'mass width identity')
            intervals = []
            for name, target in (('sqrt_low_kappa', 1/(1-2*lo)), ('sqrt_high_kappa', 1+2*u)):
                values = row[name]
                need(type(values) is list and len(values) == 2, 'two root endpoints')
                a, b = map(rational, values)
                # No square root evaluation or producer's integer square root is used.
                need(0 < a <= b and a*a <= target <= b*b and b-a <= F(1, 10**18),
                     'outward exact squared-root bounds')
                intervals.append((a, b))
            (a0, a1), (b0, b1) = intervals
            need(a0 > 1-u and b0 > 1-u, 'strict positive bending denominators')
            equal(row['bending_ratio_lower'], str((a0-(1-u))/(b1-(1-u))), 'kinetic lensing lower enclosure')
            equal(row['bending_ratio_upper'], str(min(F(1), (a1-(1-u))/(b0-(1-u)))), 'kinetic lensing upper enclosure')


def verify_abel(ctx, rows):
    catalogue = (((1, F(1, 3)),), ((2, F(2, 5)),), ((3, F(3, 7)),), ((4, F(5, 9)),),
                 ((1, F(1, 3)), (3, F(-1, 7)), (4, F(2, 9))))
    need(type(rows) is list and len(rows) == 15, 'complete Abel census')
    cursor = 0
    for terms in catalogue:
        for radius in (F(1, 2), F(1), F(3)):
            row = rows[cursor]
            cursor += 1
            keys(row, 'terms x alpha Q_minus_one')
            equal([row['terms'], row['x']], [[[p, str(c)] for p, c in terms], str(radius)], 'fixed optical control')
            x = ctx.mpf(radius.numerator)/radius.denominator
            coeff = [(p, ctx.mpf(c.numerator)/c.denominator) for p, c in terms]
            # x_ray=b/cos(theta) transforms the full infinite radial integral.
            integrals = {p:ctx.quad(lambda z:ctx.cos(z)**p, [0, ctx.pi/2]) for p, _ in terms}
            alpha = sum(2*c*x**(-p)*integrals[p] for p, c in coeff)
            close(ctx, row['alpha'], alpha, 'direct Abel forward integral')
            # Differentiate the reconstructed alpha/b, then independently integrate its inverse.
            # b=x*cosh(z) cancels the square-root endpoint singularity.
            def integrand(z):
                b = x*ctx.cosh(z)
                return sum(-2*c*(p+1)*integrals[p]*b**(-p-2) for p, c in coeff)
            inverse = -x*x/ctx.pi*ctx.quad(integrand, [0, 1, ctx.inf])
            close(ctx, row['Q_minus_one'], inverse, 'independent Abel inverse')
            close(ctx, row['Q_minus_one'], sum(c*x**(-p) for p, c in coeff), 'input optical function recovered')
