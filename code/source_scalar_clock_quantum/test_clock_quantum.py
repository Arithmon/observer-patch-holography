"""Full parent joins, exact interval controls and coherent corruption rejection."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import importlib.util, sys, unittest

HERE=Path(__file__).resolve().parent
def module(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);m=importlib.util.module_from_spec(spec);sys.modules[name]=m
    exec(compile((HERE/file).read_bytes(),str(HERE/file),'exec'),m.__dict__);return m
V=module('_test_clock_quantum_verifier','verify_clock_quantum.py')
P=module('_test_clock_quantum_producer','clock_quantum.py')

class ClockQuantumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet=V.load();cls.q,cls.c,cls.v,cls.scalar,cls.quantum,cls.clock=V.dependency_context()
    def arithmetic(self,p=None,scalar=None,quantum=None,clock=None,root=V.RER):
        return V.verify_arithmetic(self.packet if p is None else p,self.scalar if scalar is None else scalar,
            self.quantum if quantum is None else quantum,self.clock if clock is None else clock,self.v,root)
    def test_fresh_full_both_parents(self):
        r=V.verify(self.packet)
        self.assertTrue(r['full_clock_replay'] and r['full_original_quantum_replay'])
        self.assertEqual(r['parent_events_replayed'],5888);self.assertEqual(r['times'],21)
        self.assertEqual(r['resolved_steps'],[5,6,7,8,9,10,11,13,14,15,16,17,18,20,21])
        self.assertEqual(r['unresolved_steps'],[1,2,3,4,12,19])
    def test_independent_producer_byte_parity(self):self.assertEqual(P.raw(P.produce()),V.OUTPUT.read_bytes())
    def test_exact_spectral_norm_lipschitz_and_all_endpoints(self):
        R=self.v.R;L=F(self.packet['moments']['continuous_probability_lipschitz_upper'])
        product=self.v.parse(self.packet['moments']['norm_product_squared_Qphi'])
        self.assertGreaterEqual((R(4*L*L)-product).sign(),0)
        # The next coarser downward grid value must not be represented as this bound.
        self.assertLess((R(4*(L-F(1,10**12))**2)-product).sign(),0)
        for row in self.packet['rows']:
            t=self.v.parse(row['nominal_time_Qphi']);delta=R(F(row['time_offset_upper']))
            for endpoint in row['recovered_elapsed_time_interval']:
                d=R(F(endpoint))-t
                self.assertGreaterEqual((delta-d).sign(),0);self.assertGreaterEqual((delta+d).sign(),0)
            old=F(row['original_probability_error_upper']);total=F(row['total_probability_error_upper'])
            self.assertGreaterEqual(total-old,L*F(row['time_offset_upper']))
    def test_all_reference_intervals_have_correct_baseline_sign(self):
        for row in self.packet['rows']:
            lower,upper=map(F,row['continuous_probability_interval_for_every_time'])
            self.assertTrue(0<=lower<=upper<=1)
            if row['uniformly_resolved']:
                if row['signal_sign']==1:self.assertGreater(lower,F(1,2))
                else:self.assertLess(upper,F(1,2))
            self.assertEqual([str(lower-F(1,2)),str(upper-F(1,2))],row['continuous_response_interval_for_every_time'])
    def test_consumer_forgeries(self):
        mutations=[
          ('omit unresolved',lambda p:p['rows'].pop(0)),
          ('omit endpoint',lambda p:p['rows'].pop()),
          ('wrong layer',lambda p:p['rows'][15].__setitem__('reconstructed_layer_index',16)),
          ('different field',lambda p:p['rows'][15].__setitem__('layer_value_sha256','0'*64)),
          ('point interval',lambda p:p['rows'][15].__setitem__('recovered_elapsed_time_interval',['1','1'])),
          ('other step time',lambda p:p['rows'][15].__setitem__('nominal_time_Qphi',p['rows'][14]['nominal_time_Qphi'])),
          ('clock error removed',lambda p:p['rows'][15].__setitem__('timing_probability_error_upper','0')),
          ('covariance error removed',lambda p:p['rows'][15].__setitem__('original_probability_error_upper','0')),
          ('false low total',lambda p:p['rows'][15].__setitem__('total_probability_error_upper','0')),
          ('false continuous sign',lambda p:p['rows'][15].__setitem__('continuous_response_interval_for_every_time',['1/10','1/5'])),
          ('half lipschitz',lambda p:p['moments'].__setitem__('continuous_probability_lipschitz_upper',str(F(p['moments']['continuous_probability_lipschitz_upper'])/2))),
          ('different preparation',lambda p:p['configuration_join'].__setitem__('preparation_velocity_sha256','0'*64)),
          ('wrong schedule join',lambda p:p['configuration_join']['schedule_layer_value_sha256']['descending_intervention'].__setitem__(17,'0'*64)),
          ('state noise promotion',lambda p:p['scope'].__setitem__('record_noise_certifies_quantum_state_or_field_error',True)),
          ('hidden existence promotion',lambda p:p['scope'].__setitem__('interval_overlap_proves_hidden_history_exists',True)),
          ('physical clock promotion',lambda p:p['scope'].__setitem__('physical_clock_or_count_volume_identified',True)),
          ('changed vacuum',lambda p:p['scope'].__setitem__('changed_vacuum',True)),
          ('false successful time',lambda p:p['rows'][0].__setitem__('uniformly_resolved',True)),
          ('bool time count',lambda p:p['summary'].__setitem__('times',True)),
          ('new producer source',lambda p:p['source_pins'].__setitem__('code/source_scalar_clock_quantum/clock_quantum.py','0'*64)),
          ('wrong parent pin',lambda p:p['dependencies']['reconstructed_clock'].__setitem__('sha256','0'*64))]
        for name,mutate in mutations:
            with self.subTest(name=name):
                bad=deepcopy(self.packet);mutate(bad)
                with self.assertRaises(ValueError):self.arithmetic(bad)
    def test_coherent_clock_uncertainty_erasure_rejected(self):
        p=deepcopy(self.packet)
        for r in p['rows']:
            r['time_offset_upper']='0';r['timing_probability_error_upper']='0'
            r['total_probability_error_upper']=r['original_probability_error_upper']
            e=F(r['total_probability_error_upper']);lo,hi=map(F,r['split_probability_interval'])
            interval=[max(F(0),lo-e),min(F(1),hi+e)]
            r['continuous_probability_interval_for_every_time']=list(map(str,interval))
            r['continuous_response_interval_for_every_time']=[str(x-F(1,2)) for x in interval]
            r['continuous_response_abs_lower']=str(max(F(0),F(r['split_response_abs_lower'])-e))
            r['uniformly_resolved']=F(r['split_response_abs_lower'])>e
        p['summary']['maximum_timing_probability_error_upper']='0'
        with self.assertRaisesRegex(ValueError,'comparison'):self.arithmetic(p)
    def test_coherent_parent_changes_rejected(self):
        for which in ('clock','quantum','scalar'):
            data=deepcopy(getattr(self,which))
            if which=='clock':data['stability']['schedules']['ascending_intervention']['elapsed_time_intervals'][16]=['1','2']
            elif which=='quantum':data['scope']['initial_state']='different coherent preparation'
            else:data['initial_velocity_Qphi'][0]=['1','0']
            with self.subTest(parent=which),self.assertRaisesRegex(ValueError,'immutable parent'):self.arithmetic(**{which:data})
    def test_wider_uncertainty_can_remove_resolution(self):
        L=F(self.packet['moments']['continuous_probability_lipschitz_upper'])
        row=V.row_bounds(self.quantum['rows'][15],16,['0','2'],L,self.v,'declared-test-layer')
        self.assertFalse(row['uniformly_resolved']);self.assertEqual(row['continuous_response_abs_lower'],'0')
    def test_rounding_domains_and_zero_width(self):
        with self.assertRaises(ValueError):V.upper(self.v.R(-1),self.v.R,True)
        with self.assertRaises(ValueError):V.row_bounds(self.quantum['rows'][15],True,['1','2'],F(1),self.v,'x')
        with self.assertRaises(ValueError):V.row_bounds(self.quantum['rows'][15],16,['1','1'],F(1),self.v,'x')
        with self.assertRaises(ValueError):V.rational('2/4')
    def test_arithmetic_only_dependency_success_is_rejected(self):
        context=(self.q,self.c,self.v,self.scalar,self.quantum,self.clock)
        goodq={'verdict':'PASS','parent_full_mathematical_replay':True,'parent_events_replayed':5888}
        goodc={'verdict':'PASS','full_parent_mathematical_replay':True,'events_authenticated':5888,'positive_gram_cases':40}
        for which in ('quantum','clock'):
            qr=deepcopy(goodq);cr=deepcopy(goodc)
            if which=='quantum':qr['parent_full_mathematical_replay']=False
            else:cr['full_parent_mathematical_replay']=False
            with self.subTest(parent=which),patch.object(V,'dependency_context',return_value=context),patch.object(self.q,'verify',return_value=qr),patch.object(self.c,'verify',return_value=cr):
                with self.assertRaisesRegex(ValueError,'full .* proof'):V.verify(self.packet)
    def test_detached_root_and_pinned_code_mutation(self):
        with TemporaryDirectory() as d:
            root=Path(d);paths=set(V.FILES)|set(self.scalar['source_pins'])|set(self.quantum['source_pins'])|set(self.clock['source_pins'])|{V.SCALAR,V.QUANTUM,V.CLOCK}
            for p in paths:
                target=root/p;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((V.RER/p).read_bytes())
            self.arithmetic(root=root);V.dependency_context(root)
            own=root/'code/source_scalar_clock_quantum/clock_quantum.py';own.write_bytes(own.read_bytes()+b'\n# source mutation\n')
            with self.assertRaisesRegex(ValueError,'comparison'):self.arithmetic(root=root)
            code=root/V.CVERIFY;code.write_bytes(code.read_bytes()+b'\n# changed pinned clock verifier\n')
            with self.assertRaisesRegex(ValueError,'pinned clock verifier'):V.dependency_context(root)
    def test_strict_parser(self):
        with TemporaryDirectory() as d:
            p=Path(d)/'packet.json'
            for blob in (b'{"a":1,"a":2}',b'{"a":1.0}',b'{"a":NaN}',b'{"a":Infinity}',b' '*200001):
                with self.subTest(size=len(blob)):
                    p.write_bytes(blob)
                    with self.assertRaises(ValueError):V.load(p)

if __name__=='__main__':unittest.main()
