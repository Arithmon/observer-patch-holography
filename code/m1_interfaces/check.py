"""Independent integer balls, label meanings, unordered graphs and acceptance."""

from fractions import Fraction as F
from itertools import product
import hashlib
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def strict_load(path):
    from pathlib import Path
    def pairs(items):
        out = {}
        for key, value in items:
            require(key not in out, 'duplicate JSON key')
            out[key] = value
        return out
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs,
                      parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))


def root_floor(value, cap):
    low, high = 0, cap+1
    while high-low > 1:
        mid = (low+high)//2
        if mid*mid <= value: low = mid
        else: high = mid
    require(low*low <= value < (low+1)**2, 'integer square bracket')
    return low


def scale(t):
    require(type(t) is int and t in range(1,13), 'reference scale')
    base = 1 << t
    q, m, k, r = base**16, base**8, base**6, base**7
    require(m*k**4 == q*q and q == m*m and k**4 == m**3, 'coarse scaling')
    require(2*r <= m and r**5 >= m*m*k**3, 'repair and spectral window')
    require(m*q*base**4 == r**4 and r**4*base**4 == q*q, 'sharp scale factors')
    return dict(t=t, q=q, m=m, k=k, r=r, radius=m*k,
                tick_times_c_over_L=str(F(1,base**2)),
                compactness_error_factor=str(F(1,base**4)),
                added_area_factor=str(F(1,base**4)),
                alias_frequency_squared_scale=str(F(base**5)),
                spectral_saturation=str(F(base**4)))


def ball_columns(radius, q, spacing=1):
    require(type(radius) is int and 0<=radius<=128, 'ball reference scope')
    degree = flux = moment = fourth = axis2 = axis4 = cut = 0
    # Positive quadrant with multiplicities; polynomial sums independently
    # factor the producer's full signed-column enumeration.
    for x in range(radius+1):
        for y in range(radius+1):
            residual = radius*radius-x*x-y*y
            if residual < 0: continue
            z = root_floor(residual, radius)
            n = 2*z+1
            s2 = F(2*z**3,3)+z*z+F(z,3)
            s4 = F(2*z**5,5)+z**4+F(2*z**3,3)-F(z,15)
            mult = (2 if x else 1)*(2 if y else 1)
            degree += mult*n
            axis2 += mult*x*x*n
            axis4 += mult*x**4*n
            moment += int(mult*(n*(x*x+y*y)+s2))
            fourth += int(mult*(n*(x*x+y*y)**2+2*(x*x+y*y)*s2+s4))
            if x:
                my = 2 if y else 1
                flux += my*x*n
                cut += my*spacing*x*(q-spacing*y)*(n*q-spacing*z*(z+1))
    require(moment == 3*axis2, 'cubic moments')
    require(fourth <= radius*radius*moment, 'fourth-moment radius bound')
    if radius:
        require(F(3*(radius-1)**4,4)-F((2*radius+1)**2,8)<=flux<=F(22*(radius+1)**4,28),
                'spherical flux enclosure')
    return dict(degree=degree, flux=flux, moment=moment, fourth=fourth,
                axis2=axis2, axis4=axis4, coordinate_cut=cut)


def moment_case(kind):
    require(kind in ('raw','critical','repaired'), 'moment case')
    q, m, k = 65536, 256, 64
    r = dict(raw=0,critical=64,repaired=128)[kind]
    coarse, fine = ball_columns(k,q,m), ball_columns(r,q)
    remaining = [s for s in (1,2,4,8,16,32,64,128) if s>r]
    count = coarse['degree']+fine['degree']-1+6*len(remaining)
    second = m*m*coarse['moment']+fine['moment']+6*sum(s*s for s in remaining)
    fourth = m**4*coarse['fourth']+fine['fourth']+6*sum(s**4 for s in remaining)
    require(second <= 56*m*m*k**5, 'whole-symbol moment upper bound')
    require(64*second >= m*m*k**5, 'moment lower bound')
    coefficient = F(6*q*q,second)
    upper = coefficient*(F(968,49*m*m)*fine['axis2']+4*len(remaining))
    lower = max(F(0), coefficient*(F(18,m*m)*fine['axis2']-F(44**4,24*7**4*m**4)*fine['axis4']))
    if kind != 'raw': require(lower>0, 'alias lower enclosure must be positive')
    crossings = coarse['coordinate_cut']+fine['coordinate_cut']+q*q*sum(remaining)
    residue = (2*q**3//m)*(fine['flux']+sum(remaining))
    return dict(kind=kind,q=q,m=m,k=k,r=r,degree=count,
                second_norm_moment=second,fourth_norm_moment=fourth,
                coarse=coarse,fine=fine,remaining_axis_lengths=remaining,
                coordinate_cut=crossings,periodic_residue_cut=residue,
                pi_times_residue_entropy=str(F(4*residue,q**4)),
                uniform_alpha_L_squared=str(coefficient),
                one_alias_eigenvalue_L_squared_bounds=[str(lower),str(upper)])


CONFIGS = ((12,4,1,2),(16,4,1,2),(24,4,2,2))
LABELS = ('plane','residue','diagonal','checkerboard','ball','localized','connected','permuted_half')


def stencil(config, family):
    require(config in CONFIGS and family in ('T','U'), 'finite graph catalog')
    q,m,k,r = config
    radius = m*k
    require(2*radius<q, 'periodic read vectors must stay distinct')
    # Scan the actual finite displacement box and decide membership by integer
    # predicates rather than importing the generator's component union.
    result = []
    for v in product(range(-radius,radius+1),repeat=3):
        coarse = all(x%m==0 for x in v) and sum((x//m)**2 for x in v)<=k*k
        nonzero = [abs(x) for x in v if x]
        dyadic = len(nonzero)==1 and nonzero[0]<m and nonzero[0]&(nonzero[0]-1)==0
        short = family=='U' and sum(x*x for x in v)<=r*r
        if coarse or dyadic or short: result.append(v)
    return result


def labels(config):
    q,m,_,_ = config
    out = {name:[] for name in LABELS}
    for i,(x,y,z) in enumerate(product(range(q),repeat=3)):
        residue = (x//(m//2))%2==0
        ball = sum((4*p-2*q)**2 for p in (x,y,z))<q*q
        local = ball and residue
        values = (2*x<q,residue,(x+y+z)%q<q//2,(x+y+z)%2==0,
                  ball,local,local or (ball and 2*y==q and 2*z==q),
                  (i*137+29)%(q*q*q)<q*q*q//2)
        for name,value in zip(LABELS,values): out[name].append(int(value))
    return out


def graph_cuts(config, family, topology):
    require(topology in ('periodic','clipped'),'finite topology')
    q=config[0]
    vectors=stencil(config,family)
    masks=labels(config)
    cuts={name:0 for name in LABELS}
    reads=edges=0
    internal=boundary=missing=0
    for i,x in enumerate(product(range(q),repeat=3)):
        targets=set()
        field=3+x[0]+2*x[1]-x[2]
        for v in vectors:
            y=tuple(a+b for a,b in zip(x,v))
            if topology=='periodic': y=tuple(a%q for a in y)
            elif not all(0<=a<q for a in y): continue
            j=(y[0]*q+y[1])*q+y[2]
            require(j not in targets,'aliased duplicate parent')
            targets.add(j)
            reads+=1
            if i<j:
                edges+=1
                for name in LABELS: cuts[name]+=masks[name][i]!=masks[name][j]
                internal+=(field-(3+y[0]+2*y[1]-y[2]))**2
        absent=len(vectors)-len(targets)
        missing+=absent
        boundary+=absent*field*field
    require(reads==2*edges+q**3,'read/pair census mismatch')
    require(missing==q**3*len(vectors)-reads,'missing boundary slot census')
    require((missing==0)==(topology=='periodic'),'constant-mode boundary control')
    moment=sum(sum(a*a for a in v) for v in vectors)
    action=dict(boundary='periodic' if topology=='periodic' else 'fixed_zero_exterior',
                nonzero_diagonal_degree=len(vectors)-1,second_norm_moment=moment,
                missing_ordered_slots=missing,constant_quadratic=missing,
                affine_internal_quadratic=internal,affine_boundary_quadratic=boundary,
                affine_quadratic=internal+boundary,
                scaled_affine_quadratic=str(F(6*(internal+boundary),q*moment)))
    return masks,cuts,reads,action


def connected(config, mask):
    q=config[0]
    occupied={i for i,value in enumerate(mask) if value}
    if not occupied: return False
    seen={next(iter(occupied))};stack=list(seen)
    while stack:
        i=stack.pop();x=(i//(q*q),(i//q)%q,i%q)
        for axis,sign in product(range(3),(-1,1)):
            y=list(x);y[axis]+=sign
            if not all(0<=a<q for a in y): continue
            j=(y[0]*q+y[1])*q+y[2]
            if j in occupied and j not in seen: seen.add(j);stack.append(j)
    return seen==occupied
