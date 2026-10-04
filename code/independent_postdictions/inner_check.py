"""Check the measured-pixel diagnostic using inverse running couplings."""
from functools import lru_cache
from mpmath import MPContext


@lru_cache(maxsize=4)
def measured_gauge(pixel):
    c = MPContext()
    c.dps = 60
    P = c.mpf(str(pixel))
    tables = []
    tables.append([(n+1, c.mpf(n*(n+2))/4) for n in range(121)])
    tables.append([((p+1)*(q+1)*(p+q+2)//2, c.mpf(p*p+q*q+p*q+3*p+3*q)/3)
                   for p in range(91) for q in range(91)])
    tables = [[(d, k, c.log(d)) for d, k in terms] for terms in tables]
    def entropy(index, inverse):
        weights = [(d*c.exp(-4*c.pi*c.pi*k/inverse), logd) for d, k, logd in tables[index]]
        return c.fsum(w*logd for w, logd in weights)/c.fsum(w for w, _ in weights)
    def equations(x, y):
        # x=alpha_2^-1, y=alpha_3^-1; solve for these instead of alpha_U,mZ/v.
        unified_inverse = (3*x+y)/4
        z = c.sqrt(c.pi*(1/x+3/(12*x-7*y)))
        scale = -2*c.pi+c.mpf(2)/3*c.log(P)+c.pi*unified_inverse/2-c.log(z)
        return entropy(0, x)+entropy(1, y)-P/4, c.pi*(x-y)/2-scale
    x, y = c.findroot(equations, (30, 8), tol=c.mpf('1e-50'), maxsteps=30)
    return 4/(3*x+y)
