"""Exact rational filter controls and independent Fourier integrals."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent

def variance(lam,coefficients,rates,kappa=1):
    a=F(str(lam))/2;rates=list(map(F,rates));coefficients=list(map(F,coefficients))
    if a<=0 or any(g<=0 for g in rates): raise ValueError('Positive decay rates required')
    return F(kappa)*sum(c*d/(g+h)*(1/(g+a)+1/(h+a))
                         for c,g in zip(coefficients,rates,strict=True)
                         for d,h in zip(coefficients,rates,strict=True))

def fourier_variance(lam,coefficients,rates):
    a=lam/2
    def integrand(omega):
        W=sum(c/(g+1j*omega) for c,g in zip(coefficients,rates,strict=True))
        return abs(W)**2*a/(a*a+omega*omega)
    # Even integrand, convention dω/(2π) and native temporal PSD 2a/(a²+ω²).
    return 2/np.pi*quad(integrand,0,np.inf,epsabs=2e-12,epsrel=2e-11,limit=400)[0]

def execute():
    spec=json.loads((HERE/'spec.json').read_text());rows=[]
    for kernel in spec['filters']:
        vals=[]
        for lam in spec['lambdas']:
            v=variance(lam,kernel['coefficients'],kernel['rates'])
            numerical=fourier_variance(lam,kernel['coefficients'],kernel['rates'])
            assert abs(numerical-float(v))<2e-10*max(1,float(v))
            vals.append({'lambda':lam,'variance_exact':str(v),'lambda_variance_exact':str(F(str(lam))*v),
                         'fourier_variance':numerical})
        powers=[F(x['lambda_variance_exact']) for x in vals]
        assert all(b>=a for a,b in zip(powers,powers[1:]))
        rows.append({'name':kernel['name'],'values':vals,'monotone_exact':True})
    return {'schema':'oph.temporal-filter-certificate.v1',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'spec_sha256':hashlib.sha256((HERE/'spec.json').read_bytes()).hexdigest(),
            'scope':'finite rational examples independently checked by quadrature; general theorem proved analytically in REPORT.md',
            'filters':rows}

if __name__=='__main__':
    result=execute();(HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Three temporal kernels certified, including signed zero-DC memory')
