#!/usr/bin/env python3
"""True implicit scale attachment, with declared finite representation cutoffs.

The finite-iteration receipt remains a separate object. Here a rectangle
self-map and contraction prove a unique exact scale before interval AD is
applied to its implicit derivative. No infinite heat-sum limit is asserted.
"""
from pathlib import Path
from types import ModuleType
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = HERE / "hierarchy_interval_witness.py"
w = ModuleType("hierarchy_interval_parent")
w.__file__ = str(PARENT)
exec(compile(PARENT.read_bytes(), str(PARENT), "exec"), w.__dict__)
iv, Dual = w.iv, w.Dual
dexp, dlog, dsqrt = w.dexp, w.dlog, w.dsqrt
N2, N3 = 128, 64
OUTPUT = HERE.parent / "certificates/implicit_scale_attachment_receipt.json"


def scale(a, mu):
    C = Dual.const if isinstance(a, Dual) or isinstance(mu, Dual) else iv.mpf
    P, pi = iv.mpf(w.P_STR), iv.pi
    MU = iv.exp(-2*pi)*P**(iv.mpf(1)/6)
    v = C(P**(-iv.mpf(1)/2))*dexp(C(-2*pi)/(C(4)*a))
    def alpha(b):
        return C(1)/(C(1)/a+C(b/(2*pi))*dlog(C(MU)/mu))
    return v/C(2)*dsqrt(C(4*pi)*alpha(iv.mpf(1))
                        +C(4*pi*iv.mpf(3)/5)*alpha(iv.mpf(33)/5))

def residual(a, P, b1, b2, b3, mu):
    """The declared hierarchy readback Phi(a): heat-kernel edge entropies minus P/4.

    `a` may be a Dual (interval dual number) or an interval. All constants are
    promoted through C() so no mpmath-side reflected operator ever fires on a
    Dual operand.
    """
    dual_mode = isinstance(a, Dual)

    def C(x):
        x = iv.mpf(x) if not hasattr(x, 'a') else x
        return Dual.const(x) if dual_mode else x

    one = iv.mpf(1)
    pi = iv.pi
    MU = iv.exp(-2 * pi) * iv.mpf(P) ** (one / 6)
    v = C(iv.mpf(P) ** (-one / 2)) * dexp(C(-2 * pi) / (C(4) * a))

    def alpha(mu, b):
        return C(1) / (C(1) / a + C(b / (2 * pi)) * dlog(C(MU) / mu))

    def fmu(mu):
        return v / C(2) * dsqrt(
            C(4 * pi) * alpha(mu, b2) + C(4 * pi) * C(iv.mpf(3) / 5) * alpha(mu, b1))

    t2 = C(4 * pi ** 2) * alpha(mu, b2)
    t3 = C(4 * pi ** 2) * alpha(mu, b3)

    def ell_su2(t):
        Z = C(0)
        S = C(0)
        for n in range(N2 + 1):
            j = iv.mpf(n) / 2
            d = 2 * j + 1
            Cas = j * (j + 1)
            w = C(d) * dexp(-(t * C(Cas)))
            Z = Z + w
            S = S + w * C(iv.log(d))
        return S / Z

    def ell_su3(t):
        Z = C(0)
        S = C(0)
        for p in range(N3 + 1):
            for q in range(N3 + 1):
                d = iv.mpf((p + 1) * (q + 1) * (p + q + 2)) / 2
                Cas = iv.mpf(p * p + q * q + p * q + 3 * p + 3 * q) / 3
                w = C(d) * dexp(-(t * C(Cas)))
                Z = Z + w
                S = S + w * C(iv.log(d))
        return S / Z

    return ell_su2(t2) + ell_su3(t3) - C(iv.mpf(P) / 4)



def _inside(value, bounds):
    outer = iv.mpf(bounds)
    if not (outer.a <= value.a and value.b <= outer.b):
        raise ValueError("computed enclosure exceeds the declared outer bound")


def _produce():
    a = iv.mpf([w.I_LO, w.I_HI])
    y = iv.mpf(["7e-18", "8e-18"])
    image = scale(a,y)
    slope = scale(Dual.const(a),Dual.var(y)).d
    _inside(image,["7.49e-18","7.52e-18"])
    _inside(slope,["0","0.006"])
    if not (y.a < image.a and image.b < y.b and slope.b < iv.mpf("0.01")):
        raise ValueError("implicit scale has no certified contraction rectangle")
    for _ in range(8):
        y=scale(a,y)
    tangent=scale(Dual.var(a),Dual.const(y)).d/(1-scale(Dual.const(a),Dual.var(y)).d)
    partial=residual(Dual.var(a),w.P_STR,iv.mpf(33)/5,iv.mpf(1),iv.mpf(-3),Dual(y,tangent))
    _inside(partial.d,["-12","-10"])
    for endpoint, bounds in [(w.I_LO,["0.0000109907","0.0000109908"]),
                              (w.I_HI,["-0.0000109904","-0.0000109902"])]:
        aa=iv.mpf(endpoint);yy=iv.mpf(["7e-18","8e-18"])
        for _ in range(16):
            yy=scale(aa,yy)
        _inside(residual(aa,w.P_STR,iv.mpf(33)/5,iv.mpf(1),iv.mpf(-3),yy),bounds)
    source_paths=[Path(__file__).resolve(),PARENT,HERE.parent/"validators/verify_implicit_scale_attachment.py",
                  HERE.parent/"test_interval_witness.py",ROOT/"code/interval_decimal.py"]
    return {
      "schema":"oph-hierarchy-implicit-scale-v1",
      "pixel_decimal":w.P_STR,
      "coupling_interval":[w.I_LO,w.I_HI],
      "scale_rectangle":["7e-18","8e-18"],
      "scale_image_outer":["7.49e-18","7.52e-18"],
      "scale_partial_mu_outer":["0","0.006"],
      "contraction_bound":"1/100",
      "endpoint_residual_outer":[["0.0000109907","0.0000109908"],["-0.0000109904","-0.0000109902"]],
      "implicit_residual_derivative_outer":["-12","-10"],
      "finite_representation_cutoffs":{"su2_max_n":128,"su3_max_p_and_q":64},
      "interval_image_refinements":{"whole_coupling_interval":8,"each_endpoint":16},
      "scope":{"true_implicit_scale":True,"unique_scale_on_declared_rectangle":True,
               "unique_hierarchy_root_on_declared_interval":True,
               "infinite_representation_sums":False,"physical_hierarchy_attachment":False,
               "pixel_target_independent":False},
      "source_pins":{p.relative_to(ROOT).as_posix():{"bytes":p.stat().st_size,
                     "sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in source_paths}
    }


def produce():
    previous=iv.dps
    try:
        iv.dps=60
        return _produce()
    finally:
        iv.dps=previous


if __name__=="__main__":
    OUTPUT.write_bytes((json.dumps(produce(),indent=2)+"\n").encode("utf-8"))
    print(OUTPUT)
