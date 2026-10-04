"""Fresh simultaneous dimensionless solve; no imports from either producer.

The retained alpha certificates and measured tau are comparison inputs only.
All internal masses are ratios to v; no measured electroweak scale is used.
"""
from functools import lru_cache
import itertools

from mpmath import MPContext


def context():
    c = MPContext()
    c.dps = 65
    return c


def heat_table(c, group):
    # Collect degenerate Casimirs before summing, unlike the producer's
    # separate normalized probability for every representation.
    table = {}
    labels = ((n,) for n in range(121)) if group == 2 else itertools.product(range(91), repeat=2)
    for label in labels:
        if group == 2:
            n, = label
            d, casimir = n+1, c.mpf(n*(n+2))/4
        else:
            p, q = label
            d = (p+1)*(q+1)*(p+q+2)//2
            casimir = c.mpf(p*p+q*q+p*q+3*p+3*q)/3
        a, b = table.get(casimir, (c.mpf(0), c.mpf(0)))
        table[casimir] = a+d, b+d*c.log(d)
    return table


def mean_log(c, table, t):
    terms = [(c.exp(-t*k), a, b) for k, (a, b) in table.items()]
    return c.fsum(w*b for w, a, b in terms)/c.fsum(w*a for w, a, b in terms)


def stage5(c):
    roots = sorted(1+c.sqrt(2)*c.cos(c.mpf(2)/9+2*c.pi*k/3) for k in range(3))
    menu = []
    for electron in range(5, 12):
        ns = (electron, 4, 3)
        score = max(abs(2*c.log(roots[i]/roots[j])+(ns[i]-ns[j])*c.log(6))
                    for i, j in itertools.combinations(range(3), 2))
        menu.append((electron, score))
    chosen = min(menu, key=lambda x: x[1])[0]
    ns = (chosen, 4, 3)
    factor = c.root(2, 6)/c.root(c.fprod(r*r*c.sqrt(2)*6**n for r, n in zip(roots, ns)), 3)
    masses = dict(zip(('e', 'mu', 'tau'), (factor*r*r for r in roots)))
    masses.update({name: 1/(c.sqrt(2)*6**n) for name, n in
                   zip(('u', 'c', 't', 'd', 's', 'b'), (6, 3, 0, 6, 4, 2))})
    return menu, masses


def kernel(c, z, mass, asymptotic=False):
    ratio = z/mass
    if asymptotic:
        return (2*c.log(ratio)-c.mpf(5)/3)/(3*c.pi)
    a = ratio*ratio/4
    integral = -c.mpf(5)/18+1/(6*a)+(2*a-1)*c.sqrt(1+a)*c.asinh(c.sqrt(a))/(6*a**c.mpf('1.5'))
    return 2*integral/c.pi


def equations(c, tables, masses, P, g, z, mode):
    # Eliminate the Planck energy and tiny dimensionful scales before solving.
    log_ratio = -2*c.pi+c.mpf(2)/3*c.log(P)+c.pi/(2*g)-c.log(z)
    couplings = [1/(1/g+b*log_ratio/(2*c.pi)) for b in (c.mpf(33)/5, 1, -3)]
    a1, a2, a3 = couplings
    anchor = 1/a2+5/(3*a1)
    lepton = c.fsum(kernel(c, z, masses[k], mode == 'asymptotic') for k in ('e', 'mu', 'tau'))
    quark = c.fsum(kernel(c, z, masses[k], mode == 'asymptotic')*weight
                   for k, weight in (('u', c.mpf(4)/3), ('c', c.mpf(4)/3),
                                     ('d', c.mpf(1)/3), ('s', c.mpf(1)/3), ('b', c.mpf(1)/3)))
    readout = anchor+lepton+(1-3*a3/c.pi)*quark+(g if mode == 'gauge_width' else 0)
    phi = (1+c.sqrt(5))/2
    residuals = (z*z-c.pi*(a2+3*a1/5),
                 mean_log(c, tables[0], 4*c.pi*c.pi*a2)+mean_log(c, tables[1], 4*c.pi*c.pi*a3)-P/4,
                 (P-phi)*readout-c.sqrt(c.pi))
    return residuals, readout, couplings, anchor, lepton, quark


@lru_cache(maxsize=1)
def alpha_replay():
    c = context()
    tables = (heat_table(c, 2), heat_table(c, 3))
    menu, masses = stage5(c)
    fmt = lambda x: c.nstr(x, 45)
    rows = []
    for mode in ('structured', 'gauge_width', 'asymptotic'):
        f = lambda P, g, z: equations(c, tables, masses, P, g, z, mode)[0]
        P, g, z = c.findroot(f, (c.mpf('1.63'), c.mpf('.04'), c.mpf('.37')),
                            tol=c.mpf('1e-55'), maxsteps=30)
        residuals, inverse, couplings, anchor, lepton, quark = equations(c, tables, masses, P, g, z, mode)
        rows.append(dict(mode=mode, P=fmt(P), alpha_inverse=fmt(inverse), alpha_U=fmt(g),
                         mZ_over_v=fmt(z), couplings=list(map(fmt, couplings)),
                         anchor_inverse=fmt(anchor), lepton_transport=fmt(lepton),
                         unscreened_quark_transport=fmt(quark),
                         printed_residual_tolerance='1e-38'))
        if max(map(abs, residuals)) >= c.mpf('1e-50'):
            raise ValueError('independent simultaneous solve did not converge')
    return dict(rows=rows, electron_menu=[dict(exponent=n, score=fmt(s)) for n, s in menu],
                masses_over_v={k: fmt(v) for k, v in masses.items()})


def measured_pixel_diagnostic(inverse):
    c = context()
    P = (1+c.sqrt(5))/2+c.sqrt(c.pi)/c.mpf(inverse)
    tables = (heat_table(c, 2), heat_table(c, 3))
    _, masses = stage5(c)
    f = lambda g, z: equations(c, tables, masses, P, g, z, 'structured')[0][:2]
    g, z = c.findroot(f, (c.mpf('.04'), c.mpf('.37')), tol=c.mpf('1e-55'))
    base = c.mpf(alpha_replay()['rows'][0]['alpha_inverse'])
    return dict(P_measured=c.nstr(P, 45), alpha_U_measured_pixel=c.nstr(g, 45),
                mixed_inverse=c.nstr(base+g, 45))


def tau_replay(electron='0.51099895069', muon='105.6583755', *, ordered=True):
    c = context()
    e, m = c.mpf(electron), c.mpf(muon)
    # Solve the original ratio equation directly, not the producer's quadratic.
    Q = lambda t: (e+m+t)/(c.sqrt(e)+c.sqrt(m)+c.sqrt(t))**2-c.mpf(2)/3
    starts = (16*m, 18*m) if ordered else (m/40, m/20)
    return c.findroot(Q, starts, tol=c.mpf('1e-55'))


if __name__ == '__main__':
    import json
    print(json.dumps(alpha_replay(), indent=2))
    print('tau from electron and muon alone:', tau_replay())
