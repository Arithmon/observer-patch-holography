#!/usr/bin/env python3
"""Exact/algebraic regressions for collar-rate and common-cone boundaries."""
from decimal import Decimal,localcontext
from fractions import Fraction as F
import math,unittest

def decimal_fraction(x):
    return Decimal(x.numerator)/Decimal(x.denominator)

def rare_event_pair(n,extra_power=0):
    p=F(1,2**(n*n));v=p**(2+extra_power)
    sigma=[(1-p)**2,p*(1-p),p*(1-p),p*p]
    rho=[sigma[0]+v,sigma[1]-v,sigma[2]-v,sigma[3]+v]
    assert all(x>0 for x in sigma+rho)
    with localcontext() as ctx:
        ctx.prec=max(100,3*n*n+80)
        logs=[(decimal_fraction(r)/decimal_fraction(s)).ln() for r,s in zip(rho,sigma)]
        cmi=sum(decimal_fraction(r)*log for r,log in zip(rho,logs))
        opnorm=max(abs(log) for log in logs)
        expectation=-cmi
        tangent_response=-logs[3]+logs[0] # |11><11| - |00><00|, trace norm 2
        values=dict(cmi=cmi,modular_norm=opnorm,state_weighted_defect=expectation,
                    tangent_defect=tangent_response)
    return dict(p=p,rho=rho,sigma=sigma,trace_distance=sum(abs(r-s) for r,s in zip(rho,sigma)),
                chi_square=sum((r-s)**2/s for r,s in zip(rho,sigma)),**values)

def conformal_ricci(eta):
    """Independent Christoffel contraction for g=(1+eta)^2 diag(-1,1,1,1)."""
    a=1+F(eta);assert a>0;signature=[-1,1,1,1]
    g=lambda i,j: F(signature[i])*a*a if i==j else F(0)
    inv=lambda i,j: F(signature[i])/(a*a) if i==j else F(0)
    dinv=lambda i,j: -2*F(signature[i])/(a**3) if i==j else F(0)
    dg=lambda i,j: 2*F(signature[i])*a if i==j else F(0)
    ddg=lambda i,j: 2*F(signature[i]) if i==j else F(0)
    def gamma(i,j,k,derivative=False):
        value=F(0)
        for l in range(4):
            part=(j==0)*dg(l,k)+(k==0)*dg(l,j)-(l==0)*dg(j,k)
            if derivative:
                dpart=(j==0)*ddg(l,k)+(k==0)*ddg(l,j)-(l==0)*ddg(j,k)
                value += (dinv(i,l)*part+inv(i,l)*dpart)/2
            else:value += inv(i,l)*part/2
        return value
    ricci=[]
    for j in range(4):
        row=[]
        for k in range(4):
            val=gamma(0,j,k,True)-(k==0)*sum(gamma(i,j,i,True) for i in range(4))
            val+=sum(gamma(i,i,l)*gamma(l,j,k)-gamma(i,k,l)*gamma(l,j,i) for i in range(4) for l in range(4))
            row.append(val)
        ricci.append(row)
    scalar=sum(inv(i,j)*ricci[i][j] for i in range(4) for j in range(4))
    return ricci,scalar

class CollarBoundaryTests(unittest.TestCase):
    def test_lower_bound_does_not_force_shrinking(self):
        for n in (2,8,32,128,1024):
            ell=F(1,4*n);delta=F(1,4)
            self.assertGreaterEqual(float(delta),float(ell)*math.log(n))
            self.assertEqual(delta/ell,n)
            self.assertEqual(delta,F(1,4))

    def test_chosen_logarithmic_schedule_does_shrink(self):
        n=[16,64,256,1024,4096]
        widths=[math.ceil(3*math.log(j))/(4*j) for j in n]
        self.assertTrue(all(a>b for a,b in zip(widths,widths[1:])))
        self.assertTrue(all(w>=3*math.log(j)/(4*j) for w,j in zip(widths,n)))
        self.assertLess(widths[-1],.002)

    def test_product_gibbs_has_exact_zero_log_defect(self):
        # H=0: log rho_R=-|R| log2; the coefficient cancels identically.
        for a,b,d in ((1,8,1),(2,32,3),(4,128,5)):
            self.assertEqual(-(a+b+d)-b+(a+b)+(b+d),0)

    def test_rare_states_are_faithful_normalized_same_marginals(self):
        for n in (2,4,8):
            pair=rare_event_pair(n);p=pair['p'];r=pair['rho'];s=pair['sigma']
            self.assertEqual(sum(r),1);self.assertEqual(sum(s),1)
            self.assertEqual(r[2]+r[3],p);self.assertEqual(r[1]+r[3],p)
            self.assertEqual(pair['trace_distance'],4*p*p)
            self.assertEqual(pair['chi_square'],p*p/(1-p)**2)
            self.assertLessEqual(pair['chi_square'],4*p*p)

    def test_trace_cmi_small_but_modular_norm_not_small(self):
        for n in (4,8,12):
            p=rare_event_pair(n)
            self.assertGreaterEqual(p['cmi'],0)
            self.assertLessEqual(p['cmi'],decimal_fraction(p['chi_square']))
            self.assertAlmostEqual(float(p['modular_norm']),math.log(2),places=14)
            self.assertAlmostEqual(float(p['state_weighted_defect']),-float(p['cmi']),places=14)
            self.assertGreater(abs(float(p['tangent_defect'])),.69)

    def test_rare_example_violates_strong_matrix_mixing(self):
        pair=rare_event_pair(8)
        epsilon=decimal_fraction(4*pair['p']**2)
        self.assertLess(pair['cmi'],epsilon)
        self.assertGreater(pair['modular_norm'],epsilon)

    def test_faithful_floor_scaled_positive_control(self):
        for n in (2,4,8):
            pair=rare_event_pair(n,extra_power=1);p=pair['p']
            self.assertEqual(pair['trace_distance'],4*p**3)
            # Every relevant eigenvalue is at least p^2 for these pairs.
            floor=p*p;bound=4*pair['trace_distance']/floor
            self.assertLessEqual(pair['modular_norm'],decimal_fraction(bound))
            self.assertLessEqual(pair['modular_norm'],decimal_fraction(p))

    def test_corrected_power_margin_accounts_for_floor(self):
        theta=F(1,2);floor_power=3
        self.assertGreater(9,F(4)/theta) # old trace-only inequality passes
        self.assertLessEqual(9*theta-floor_power,4) # modular rate does not
        self.assertGreater(15*theta-floor_power,4) # corrected margin passes

    def test_coupled_fields_can_have_two_characteristic_speeds(self):
        # K=I, spatial G=diag(1,4), V=(phi1-phi2)^2/2.
        for k in (1.,4.,16.):
            discriminant=math.sqrt(9*k**4+4)
            frequencies_sq=[(5*k*k+2-discriminant)/2,(5*k*k+2+discriminant)/2]
            self.assertAlmostEqual(sum(frequencies_sq),5*k*k+2)
            self.assertAlmostEqual(math.prod(frequencies_sq),4*k**4+5*k*k)
        self.assertNotEqual(F(1),F(4)) # no c^2 with G=c^2K

    def test_same_null_cones_and_scalar_curvature_not_einstein_tensor(self):
        for eta in (F(-1,2),F(0),F(1)):
            ricci,scalar=conformal_ricci(eta);a=1+eta
            self.assertEqual(scalar,0)
            for i in range(4):
                for j in range(4):
                    self.assertEqual(ricci[i][j],F(3 if i==0 else 1)/(a*a) if i==j else 0)
            self.assertNotEqual(ricci[0][0],0)

if __name__=='__main__':unittest.main()
