"""Independent integer checks of the entropy-preserving sparse construction."""

from fractions import Fraction as F
from itertools import product

from . import check


def parameters(t):
    check.require(type(t) is int and 1<=t<=8,'balanced level')
    return 1<<(8*t),1<<(4*t),1<<(3*t)


def admissible(t,raw):
    v = check.vector(raw)
    _,m,k = parameters(t)
    if all(x%m==0 for x in v):
        return sum((x//m)**2 for x in v)<=k*k
    nonzero = [abs(x) for x in v if x]
    return len(nonzero)==1 and nonzero[0]<m and nonzero[0]&(nonzero[0]-1)==0


def route_targets():
    targets=[]
    for t in range(1,7):
        q,m,k=parameters(t)
        targets.extend((t,d) for d in ((0,0,0),(q-1,0,0),(q//2,-q//3,q//5),
                                       (-q+1,q//2,-1),(m//2,-m//2,1),(m*k+1,m-1,-2)))
    targets.extend((1,d) for d in product(range(-2,3),repeat=3))
    return targets


def check_route(t,d,steps):
    _,m,k=parameters(t)
    d=check.vector(d)
    check.require(type(steps) is list,'balanced route list')
    point=[0,0,0];points=[point.copy()]
    for raw in steps:
        v=check.vector(raw)
        check.require(admissible(t,v),'illegal balanced edge')
        point=[a+b for a,b in zip(point,v)];points.append(point)
    check.require(tuple(point)==d,'balanced endpoint')
    square=sum(x*x for x in d)
    check.require(square<=(len(steps)*m*k)**2,'balanced outer speed')
    left=(len(steps)-12*t-1)*(k-2)-2
    check.require(left<=0 or (left*m)**2<=square,'balanced charged time bound')
    for point in points:
        dot=sum(x*y for x,y in zip(point,d))
        if square and 0<dot<square:
            distance=F(sum(x*x for x in point))-F(dot*dot,square)
        else:
            endpoint=d if square and dot>=square else (0,0,0)
            distance=F(sum((x-y)**2 for x,y in zip(point,endpoint)))
        check.require(distance<=4*m*m,'balanced boundary tube')
    return dict(steps=len(steps),points_sha256=check.digest(points))


def certificate(t):
    check.require(type(t) is int and t in (1,2,3),'balanced cut level')
    q,m,k=parameters(t)
    coordinate=diagonal=moment=fourth=degree=flux=0
    # Octant enumeration and exact squared brackets independently reconstruct
    # the column endpoints; polynomial sums use a different factorization.
    for x in range(k+1):
        for y in range(k+1):
            residual=k*k-x*x-y*y
            if residual<0: continue
            z=check.floor_ratio_root(residual,1,k)
            s2=F(2*z**3,3)+z*z+F(z,3)
            s4=F(2*z**5,5)+z**4+F(2*z**3,3)-F(z,15)
            multiplicity=(2 if x else 1)*(2 if y else 1)
            count=2*z+1;xy=x*x+y*y
            degree+=multiplicity*count
            moment+=int(multiplicity*m*m*(count*xy+s2))
            fourth+=int(multiplicity*m**4*(count*xy*xy+2*xy*s2+s4))
            width=count*q-m*z*(z+1)
            if x:
                mult_y=2 if y else 1
                flux+=mult_y*x*count
                coordinate+=mult_y*m*x*(q-m*y)*width
            if t<=2:
                for sx,sy in product((-x,x) if x else (0,),(-y,y) if y else (0,)):
                    a,b=m*sx,m*sy
                    if a+b>0:
                        rectangle=check.ramp_rows(max(0,-a),min(q-1,q-1-a),
                                                  max(0,-b),min(q-1,q-1-b),q-a-b,q-1)
                        diagonal+=rectangle*width
    moment+=2*(m*m-1)
    fourth+=2*(m**4-1)//5
    degree+=24*t
    coordinate+=q*q*(m-1)
    expected=dict(q=q,spacing=m,coarse_radius=k,radius=m*k,degree=degree,
                  coordinate=coordinate,coarse_flux=flux,second_norm_moment=moment,
                  fourth_norm_moment=fourth)
    expected['stabilized']={
        'coarse_edge_weight_times_h_squared':str(F(3,moment)),
        'nearest_extra_weight_times_h_squared':str(F(1,2)),
        'fourth_error_coefficient_times_inverse_h_squared':str(F(moment+fourth,8*moment)),
        'naive_alias_modes':1<<(12*t),
        'naive_alias_eigenvalue_times_L_squared_upper':str(48*t*F(6*q*q,moment)),
        'stabilized_nonzero_alias_eigenvalue_times_L_squared_lower':8*(q//m)**2}
    if t<=2:
        expected['diagonal']=diagonal+2*q*(q*(m-1)-(m*m-1)//3)
    check.require(m*k**4==q*q,'balanced entropy scaling')
    check.require(m*(q-m*k)**2*flux<=coordinate<=m*q*q*flux+q*q*(m-1),'balanced cut bounds')
    # 3 < pi < 22/7 suffices for exact finite controls; the analytic limit
    # uses the actual spherical integral, not this rational enclosure.
    check.require(F(3*(k-1)**4,4)-F((2*k+1)**2,8)<=flux<=F(22*(k+1)**4,28),
                  'coarse spherical flux bounds')
    check.require(fourth<=(m*k)**2*moment,'balanced Taylor bound')
    check.require(moment<=(degree-1)*(m*k)**2,'balanced mean eigenvalue obstruction')
    check.require(64*moment>=m*m*k**5,'alias upper bound moment floor')
    check.require(F(moment+fourth,8*moment)<=F((m*k)**2,4),'stabilized Taylor bound')
    return expected


def verify(packet):
    check.require(type(packet) is dict and set(packet)=={'cuts','routes'},'balanced schema')
    check.require(type(packet['cuts']) is dict and set(packet['cuts'])=={'1','2','3'},'balanced cut catalog')
    for t in (1,2,3):
        check.require(check.canonical(packet['cuts'][str(t)])==check.canonical(certificate(t)),
                      'balanced cut or action certificate')
    targets=route_targets()
    rows=packet['routes']
    check.require(type(rows) is list and len(rows)==len(targets),'balanced route catalog')
    total=0
    for row,(t,d) in zip(rows,targets):
        check.require(type(row) is dict and set(row)=={'t','displacement','steps','points_sha256','instructions'},'balanced route schema')
        check.require(check.canonical([row['t'],row['displacement']])==check.canonical([t,d]),'balanced route target')
        check.require(type(row['instructions']) is list,'balanced instruction list')
        expanded=[]
        for item in row['instructions']:
            check.require(type(item) is list and len(item)==4 and all(type(x) is int for x in item),'balanced instruction')
            count,*v=item
            check.require(1<=count<=10000 and len(expanded)+count<=10000,'balanced repetition')
            expanded.extend([tuple(v)]*count)
        result=check_route(t,d,expanded)
        check.require(check.canonical(result)==check.canonical({key:row[key] for key in result}),'balanced path commitment')
        total+=result['steps']
    return dict(balanced_routes=len(rows),balanced_steps=total)
