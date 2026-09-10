"""Exact algebraic checks for a conditional classical Cartan reduction."""
from pathlib import Path
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[2]


def derive():
    k1, k2 = s.symbols('k1 k2', positive=True)
    q = s.Matrix([s.Rational(1, 2)]*2)
    K = s.diag(k2, 2*k1)
    d = s.factor((q.T*K.inv()*q)[0])
    w = s.simplify(K.inv()*q/d)
    residual = list(K*w-q/d)+[(q.T*w)[0]-1, (w.T*K*w)[0]-1/d]
    code_q, code_K = s.Matrix([1, s.Rational(1, 2)]), s.diag(4*k2, 2*k1)
    residual.append((code_q.T*code_K.inv()*code_q)[0]-d)
    x, y, u, v = s.symbols('x y u v', real=True)
    H, DH = s.Matrix([x+s.I*y, 0]), s.Matrix([u+s.I*v, 0])
    generators = [s.I*s.Matrix([[0, 1], [1, 0]])/2,
                  s.Matrix([[0, 1], [-1, 0]])/2,
                  s.I*s.diag(1, -1)/2, s.I*s.eye(2)/2]
    currents = [s.expand(((t*H).conjugate().T*DH+DH.conjugate().T*t*H)[0])
                for t in generators]
    residual += currents[:2]+[currents[2]-currents[3]]
    # Every registered Yukawa monomial has two independent fermion factors.
    parent = json.loads((ROOT/'code/sm_local_action/local_action_receipt.json').read_text())
    ys = {k: v for k, v in parent['action_coefficients'].items() if k.startswith('Y_')}
    assert len(ys) == 54 and all(len(mask) == 2 for terms in ys.values() for mask, _ in terms)
    assert all(s.simplify(v) == 0 for v in residual)
    return {'standard_q': ['1/2', '1/2'], 'standard_K': ['k2', '2*k1'],
            'code_q': ['1', '1/2'], 'code_K': ['4*k2', '2*k1'],
            'd': str(d), 'w': [str(v) for v in w],
            'identities': [str(s.simplify(v)) for v in residual],
            'weak_and_hypercharge_currents': [str(v) for v in currents],
            'unrestricted_complex_yukawa_slots': len(ys)//2,
            'retained_yukawa_real_sectors': len(ys),
            'yukawa_monomial_degree': 2,
            'hypercharge_only_orthogonal_gauge_defect_at_k1_k2_1': '-4',
            'omitted_equations': ['color', 'weak_1', 'weak_2', 'lower_Higgs',
                                  'fermions_and_independent_conjugates', 'orthogonal_Cartan'],
            'local_log_plaquette_condition': 'abs(w_weak * curl(A)) < 2*pi',
            'general_coefficients': 'k1,k2 positive constants; m_squared,lambda real; Yukawas complex'}
