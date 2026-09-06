"""Exact dyadic interval Taylor enclosure of the finite charged action.

The phase-rotated canonical coordinates remove exponential coefficients.
All enclosure arithmetic uses integers with directed dyadic rounding.  The
time variable and action are supplied model inputs, not physical calibration.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/whitney_charged_enclosure_receipt.json"
BITS = 256
SCALE = 1 << BITS
POLY_BITS = 96
PINS = (
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "paper/tex_fragments/WHITNEY_CHARGED_MATTER.tex",
    "paper/tex_fragments/WHITNEY_CHARGED_EXECUTION.tex",
    "paper/tex_fragments/WHITNEY_EPHEMERIS_CLOCK.tex",
    "paper/tex_fragments/WHITNEY_CHARGED_ENCLOSURE.tex",
    "code/electromagnetism/whitney_charged_dynamics.py",
    "code/electromagnetism/verify_whitney_charged_dynamics.py",
    "code/electromagnetism/whitney_charged_enclosure.py",
    "code/electromagnetism/verify_whitney_charged_enclosure.py",
    "code/electromagnetism/test_whitney_charged_enclosure.py",
    "code/electromagnetism/runtime/whitney_charged_dynamics_receipt.json",
)


class Interval:
    """Closed interval [lo,hi]/2**BITS, outward rounded at every operation."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None, *, raw=False):
        if raw:
            self.lo, self.hi = lo, hi
        else:
            x = Fraction(lo)
            self.lo = (x.numerator*SCALE)//x.denominator
            self.hi = -((-x.numerator*SCALE)//x.denominator)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def cast(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        other = self.cast(other)
        return Interval(self.lo+other.lo, self.hi+other.hi, raw=True)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo, raw=True)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        v = (self.lo*other.lo, self.lo*other.hi,
             self.hi*other.lo, self.hi*other.hi)
        return Interval(min(v)//SCALE, -((-max(v))//SCALE), raw=True)

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError("interval division through zero")
        a, b = Fraction(SCALE*SCALE, self.lo), Fraction(SCALE*SCALE, self.hi)
        return Interval(min(a,b).__floor__(), max(a,b).__ceil__(), raw=True)

    def __truediv__(self, other):
        return self*self.cast(other).reciprocal()

    def __rtruediv__(self, other):
        return self.cast(other)*self.reciprocal()

    def square(self):
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo,self.hi*self.hi)
        upper = max(self.lo*self.lo,self.hi*self.hi)
        return Interval(lower//SCALE, -((-upper)//SCALE), raw=True)

    def contains(self, other):
        return self.lo <= other.lo and other.hi <= self.hi

    def interior_contains(self, other):
        return self.lo < other.lo and other.hi < self.hi

    def midpoint(self):
        m = (self.lo+self.hi)//2
        return Interval(m,m,raw=True)

    def inflate(self, factor=Fraction(9,8), minimum=1):
        m = (self.lo+self.hi)//2
        radius = max(m-self.lo,self.hi-m,minimum)
        radius = (radius*factor.numerator+factor.denominator-1)//factor.denominator
        return Interval(m-radius,m+radius,raw=True)

    def pair(self):
        return [str(self.lo),str(self.hi)]


ZERO = Interval(0)
ONE = Interval(1)


class Tape:
    """Truncated power-series arithmetic evaluated in topological order."""
    def __init__(self):
        self.nodes = []

    def node(self, op, *args):
        node = Node(self,len(self.nodes))
        self.nodes.append((op,args))
        return node

    def const(self, value):
        return self.node("constant",Interval.cast(value))

    def variable(self, index):
        return self.node("variable",index)

    def coefficients(self, initial, outputs, order):
        values = [[] for _ in self.nodes]
        solution = [[x] for x in initial]
        for n in range(order):
            for k,(op,args) in enumerate(self.nodes):
                if op == "constant":
                    v = args[0] if n == 0 else ZERO
                elif op == "variable":
                    v = solution[args[0]][n]
                elif op == "add":
                    v = values[args[0]][n]+values[args[1]][n]
                elif op == "negative":
                    v = -values[args[0]][n]
                elif op == "multiply":
                    a,b = values[args[0]],values[args[1]]
                    v = sum((a[j]*b[n-j] for j in range(n+1)),ZERO)
                elif op == "inverse":
                    a = values[args[0]]
                    v = a[0].reciprocal() if n == 0 else -sum(
                        (a[j]*values[k][n-j] for j in range(1,n+1)),ZERO)/a[0]
                else:
                    raise ValueError("unknown tape operation")
                values[k].append(v)
            for series, output in zip(solution,outputs,strict=True):
                series.append(values[output.index][n]/(n+1))
        return solution


class Node:
    def __init__(self,tape,index):
        self.tape,self.index=tape,index

    def cast(self,value):
        return value if isinstance(value,Node) else self.tape.const(value)

    def __add__(self,other):
        other=self.cast(other)
        return self.tape.node("add",self.index,other.index)

    __radd__=__add__

    def __neg__(self):
        return self.tape.node("negative",self.index)

    def __sub__(self,other):
        return self+-self.cast(other)

    def __rsub__(self,other):
        return self.cast(other)+-self

    def __mul__(self,other):
        other=self.cast(other)
        return self.tape.node("multiply",self.index,other.index)

    __rmul__=__mul__

    def __truediv__(self,other):
        other=self.cast(other)
        return self*self.tape.node("inverse",other.index)


def geometry_ratio():
    s = isqrt(5*SCALE*SCALE)
    sqrt5=Interval(s,s+1,raw=True)
    return 6/(7+3*sqrt5)


def vector_field(y,r):
    """Rational Hamilton equations, with alpha momentum fixed to zero."""
    _,cr,ci,br,bi,pr,pi,qr,qi=y
    e=Fraction(1,4)
    x=cr*cr+ci*ci
    z=br*br+bi*bi
    v=cr*br+ci*bi
    ur,ui=cr-br,ci-bi
    d=r+Fraction(2,525)*e*e*(ur*ur+ui*ui)
    w=[e*(3*ci+2*bi)/5,-e*(3*cr+2*br)/5,
       e*(bi-ci)/15,e*(cr-br)/15]
    p=[pr,pi,qr,qi]
    av=-sum((a*b for a,b in zip(w,p,strict=True)),0)/d
    hp=[8*pr-2*qr,8*pi-2*qi,-2*pr+Fraction(4,3)*qr,-2*pi+Fraction(4,3)*qi]
    xd=[a-b*av for a,b in zip(hp,w,strict=True)]
    # V=r|C-b|^2 + (1/2)<|sC+(1-s)b|^2>
    #                + (1/8)<|sC+(1-s)b|^4>.
    fx=r+Fraction(1,20)+(2*x+3*v+2*z)/280
    fv=-2*r+Fraction(3,20)+(3*x+8*v+10*z)/280
    fz=r+Fraction(3,10)+(2*x+10*v+30*z)/280
    grad=[2*cr*fx+br*fv,2*ci*fx+bi*fv,
          2*br*fz+cr*fv,2*bi*fz+ci*fv]
    dp=[e*(-3*pi/5+qi/15),e*(3*pr/5-qr/15),
        e*(-2*pi/5-qi/15),e*(2*pr/5+qr/15)]
    dd=[Fraction(4,525)*e*e*u for u in [ur,ui,-ur,-ui]]
    pd=[av*a+av*av*b/2-c for a,b,c in zip(dp,dd,grad,strict=True)]
    return [av]+xd+pd


def make_tape():
    tape=Tape()
    y=[tape.variable(i) for i in range(9)]
    rhs=vector_field(y,tape.const(geometry_ratio()))
    return tape,rhs


def initial():
    return [Interval(x) for x in [0,1,0,1,0,0,Fraction(3,10),0,Fraction(-3,10)]]


def horner(series,t):
    value=series[-1]
    for c in reversed(series[:-1]):
        value=c+t*value
    return value


def picard_box(y,h):
    time=Interval(0,h.hi,raw=True)
    field=vector_field(y,geometry_ratio())
    speed=max(max(abs(v.lo),abs(v.hi)) for v in field)
    pad=(2*h.hi*speed)//SCALE+SCALE//10**12
    box=[Interval(a.lo-pad,a.hi+pad,raw=True) for a in y]
    for _ in range(16):
        image=[a+time*b for a,b in zip(y,vector_field(box,geometry_ratio()),strict=True)]
        if all(a.interior_contains(b) for a,b in zip(box,image,strict=True)):
            return box
        box=[a if a.interior_contains(b) else
             Interval(min(a.lo,b.lo),max(a.hi,b.hi),raw=True).inflate(Fraction(3,2))
             for a,b in zip(box,image,strict=True)]
    raise ValueError("no Picard inclusion")


def advance(y,h,order,tape,outputs):
    box=picard_box(y,h)
    local=tape.coefficients(y,outputs,order)
    whole=tape.coefficients(box,outputs,order+1)
    power=ONE
    for _ in range(order+1):power=power*h
    rounding_margin=Interval(-(1 << (BITS-130)),1 << (BITS-130),raw=True)
    endpoint=[(horner(a,h)+b[-1]*power).inflate(Fraction(9,8))+rounding_margin
              for a,b in zip(local,whole,strict=True)]
    polynomial=[[str(((v.lo+v.hi)//2)//(1 << (BITS-POLY_BITS))) for v in a]
                for a in local]
    return endpoint,box,polynomial


def integrate(steps=80,order=32,progress=False):
    h=Interval(Fraction(2,steps))
    tape,outputs=make_tape()
    y=initial()
    rows=[{"state":[v.pair() for v in y]}]
    for i in range(steps):
        y,box,polynomial=advance(y,h,order,tape,outputs)
        rows[-1]["picard_box"]=[v.pair() for v in box]
        rows[-1]["polynomial"]=polynomial
        rows.append({"state":[v.pair() for v in y]})
        if progress and i%16==15:
            print(i+1, max(v.hi-v.lo for v in y)/SCALE,flush=True)
    return rows


def canonical(value):
    return (json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode("ascii")


def build():
    return {
        "schema":"oph.whitney_charged_enclosure.v1",
        "scope":"RIGOROUS_FINITE_CHARGED_ACTION_TIME_TUBE__SUPPLIED_MODEL_PARAMETER",
        "source_pins":{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in PINS},
        "parameters":{"charge":"1/4","mass_squared":"1/2","quartic":"1/4",
            "geometry_ratio":"6/(7+3*sqrt(5))","alpha_momentum":"0",
            "normalized_initial_energy":"49/40+9/(3200*r)"},
        "method":{"name":"Picard inclusion and interval Taylor remainder",
            "steps":80,"step":"1/40","time_span":["0","2"],"taylor_order":32,
            "interval_bits":BITS,"polynomial_bits":POLY_BITS,
            "coordinate_order":["alpha","Re(C)","Im(C)","Re(b)","Im(b)",
                "p_Re(C)","p_Im(C)","p_Re(b)","p_Im(b)"],
            "interval_encoding":"integer endpoint pairs divided by 2^256",
            "polynomial_encoding":"integer coefficients divided by 2^96, ascending powers of local time",
            "error_norm":"maximum absolute canonical coordinate difference",
            "approximant":"piecewise polynomial; adjacent polynomial endpoints need not agree"},
        "interpretation":{"units":"supplied dimensionless model units",
            "time":"exact supplied action parameter, not measured physical time",
            "action":"full dressed charged-scalar/Maxwell action on the fixed twenty-tetrahedron cone",
            "canonical_map":"C=exp(i*e*alpha)*c; momenta are for action divided by cone volume",
            "imports":["cone geometry","matter law and couplings","temporal gauge","exact initial data","action time"],
            "rigorous_time_enclosure":True,"whole_time_interval":True,
            "physical_clock_selected":False,"spatial_continuum_error_certified":False,
            "observer_history":False,"quantum_history":False,"empirical_comparison":False,
            "formalized_in_Lean":False,
            "historical_comparison":"old binary64 coordinates interpreted as exact dyadic numbers at nominal j/40; original float timestamps are metadata"},
        "bounds":{"uniform_canonical_polynomial_error_upper":"1/100000000000000000000",
            "historical_position_velocity_error_upper":"1/10000000000",
            "normalized_denominator_lower_bound":"2/5"},
        "states":integrate(),
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    args=parser.parse_args()
    packet=build()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt":str(args.output),"steps":80,"samples":81,
                      "bytes":args.output.stat().st_size},sort_keys=True))
