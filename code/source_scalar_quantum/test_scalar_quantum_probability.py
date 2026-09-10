"""Exact factors, immutable-parent and scientific-scope adversarial controls."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import copy
import importlib.util
import json
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('scalar_probability_independent',HERE/'verify_quantum_probability.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)

class ProbabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.p=v.load();cls.module,cls.parent=v.parent_context()

    def test_full_independent_replay(self):
        report=v.verify(self.p)
        self.assertTrue(report['parent_full_mathematical_replay']);self.assertEqual(report['parent_events_replayed'],5888)
        self.assertEqual(report['resolved_steps'],[5,6,7,8,9,10,11,13,14,15,16,17,18,20,21])
        self.assertGreater(F(report['step16']['split_response_abs_lower']),F(report['step16']['probability_error_upper']))

    def test_exact_single_mode_original_covariance(self):
        # Independent exact rational matrix powers; no trigonometric arithmetic.
        for omega in [F(1),F(2),F(3),F(7)]:
            tau=F(1,10);x=tau*tau*omega*omega/4;A=1-2*x;B=tau;C=-tau*omega*omega*(1-x)
            q,p=F(1),F(0) # top row of the evolving matrix
            for j in range(12):
                actual=q*q/(2*omega)+p*p*omega/2-F(1,2)/omega
                sin_squared=p*p*omega*omega*(1-x)
                expected=tau*tau*omega*sin_squared/(8*(1-x))
                self.assertEqual(actual,expected);self.assertGreaterEqual(actual,0)
                if j: self.assertNotEqual(actual,expected*2)
                q,p=q*A+p*C,q*B+p*A

    def test_probability_forgeries(self):
        mutations=[('variance zero',lambda p:p['covariance'].__setitem__('increase_upper','0')),
         ('original vacuum zero',lambda p:p['covariance'].__setitem__('original_variance_upper','0')),
         ('changed vacuum',lambda p:p['scope'].__setitem__('changed_vacuum',True)),
         ('physical clock',lambda p:p['scope'].__setitem__('native_clock',True)),
         ('q233 transfer',lambda p:p['scope'].__setitem__('spatial_continuum_error_transfer',True)),
         ('observed outcome',lambda p:p['scope'].__setitem__('observed_quantum_outcomes',True)),
         ('hbar',lambda p:p['scope'].__setitem__('hbar','2')),
         ('new preparation',lambda p:p['scope'].__setitem__('initial_state','modified vacuum')),
         ('drop time',lambda p:p['rows'].pop()),
         ('probability',lambda p:p['rows'][15].__setitem__('split_probability_interval',['1','1'])),
         ('zero error',lambda p:p['rows'][15].__setitem__('probability_error_upper','0')),
         ('rescale signal',lambda p:p['rows'][15].__setitem__('split_response_abs_lower','1/2')),
         ('mean sign',lambda p:p['rows'][15].__setitem__('signal_sign',1)),
         ('false early success',lambda p:p['rows'][0].__setitem__('resolved_above_total_bound',True)),
         ('moments',lambda p:p['moments'].__setitem__('g_energy_Qphi',['0','0'])),
         ('normalization',lambda p:p['constants'].__setitem__('tau_squared','1/49')),
         ('CFL',lambda p:p['constants'].__setitem__('cfl_upper','4')),
         ('boolean alias',lambda p:p['summary'].__setitem__('times',True)),
         ('source pins',lambda p:p['source_pins'].__setitem__('code/source_scalar_quantum/quantum_probability.py','0'*64)),
         ('extra scope',lambda p:p['scope'].__setitem__('whole_QFT',True))]
        for name,mutate in mutations:
            with self.subTest(name=name):
                bad=copy.deepcopy(self.p);mutate(bad)
                with self.assertRaises(ValueError):v.verify_arithmetic(bad,self.parent,self.module)

    def test_dependency_replay_is_required(self):
        with patch.object(self.module,'verify',side_effect=ValueError('forged parent')):
            with patch.object(v,'parent_context',return_value=(self.module,self.parent)):
                with self.assertRaisesRegex(ValueError,'forged parent'):v.verify(self.p)

    def test_parent_identity_not_rehashable(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/v.PARENT_PATH;p.parent.mkdir(parents=True);p.write_text('{}')
            with self.assertRaisesRegex(ValueError,'immutable parent receipt'):v.parent_context(root)
            p.write_bytes((v.RER/v.PARENT_PATH).read_bytes());q=root/v.VERIFIER_PATH;q.write_text('raise Exception("must not execute")')
            with self.assertRaisesRegex(ValueError,'immutable parent verifier'):v.parent_context(root)

    def test_detached_root_uses_its_own_source_pins(self):
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)
            paths=set(self.parent['source_pins']) | set(self.p['source_pins']) | {v.PARENT_PATH}
            for relative in paths:
                dest=target/relative;dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes((v.RER/relative).read_bytes())
            v.verify_arithmetic(self.p,self.parent,self.module,rer=target)
            forged=target/'code/source_scalar_quantum/quantum_probability.py'
            forged.write_bytes(forged.read_bytes()+b'\n# detached source change\n')
            with self.assertRaisesRegex(ValueError,'exact original-vacuum'):
                v.verify_arithmetic(self.p,self.parent,self.module,rer=target)

    def test_strict_json(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json'
            for raw in ['{"a":1,"a":2}','{"a":1.0}','{"a":NaN}','{"a":Infinity}']:
                p.write_text(raw)
                with self.assertRaises(ValueError):v.load(p)
            p.write_bytes(b' '*200001)
            with self.assertRaises(ValueError):v.load(p)

    def test_exact_negative_golden_coefficient_cancellation(self):
        R=self.module.R;x=self.module.parse(self.p['moments']['g_mass_squared_Qphi'])
        lo,hi=v.bracket(x,R);self.assertLess(F(61,1000),lo);self.assertLess(hi,F(62,1000))
        self.assertLessEqual((R(lo)-x).sign(),0);self.assertGreaterEqual((R(hi)-x).sign(),0)

if __name__=='__main__':unittest.main()
