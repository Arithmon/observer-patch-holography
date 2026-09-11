"""Exact source spacing, independent mode controls, custody and false promotions."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import importlib.util, json, shutil, subprocess, sys, unittest

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,HERE/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m
    spec.loader.exec_module(m);return m
P=module('_time_refinement_producer','source_scalar_time_refinement.py')
V=module('_time_refinement_verifier','verify_source_scalar_time_refinement.py')

class TimeRefinement(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent_verifier,cls.parent,cls.interval=V.context();cls.packet=V.load()
    def test_fresh_full_spatial_parent_and_all4levels(self):
        report=V.verify(self.packet);self.assertTrue(report['full_spatial_parent_replayed'])
        self.assertEqual(report['levels'],4);self.assertFalse(report['observer_event_log_executed'])
    def test_producer_is_byte_reproducible(self):self.assertEqual(P.raw(P.produce()),V.OUTPUT.read_bytes())
    def test_exact_golden_minimum_gaps_and_dyadic_steps(self):
        expected=[(13,[13,-8],9),(34,[34,-21],12),(89,[89,-55],14),(233,[233,-144],17)]
        for row,(q,pair,k) in zip(self.packet['levels'],expected):
            self.assertEqual(row['q'],q);self.assertEqual(row['minimum_gap_Qphi'],pair)
            self.assertEqual(row['step_exponent'],k)
            lo,hi=map(F,row['minimum_gap_interval']);tau=F(row['time_step'])
            self.assertLessEqual(tau,lo*lo);self.assertGreater(2*tau,hi*hi)
            self.assertLess(F(row['squared_step_operator_bound']),1)
    def test_coarse_comparisons_stay_unresolved(self):
        self.assertEqual([r['comparison_error_smaller_than_response'] for r in self.packet['levels']],
                         [False,False,False,True])
    def test_exact_small_mode_displacement_covariance_and_probability(self):
        # Rational KDK matrices, original-vacuum covariance and directed analytic
        # sine/exponential intervals; no trajectory or formula from the producer.
        I=self.interval.I;sin=self.interval.sin;exp=self.interval.exp;tau=F(1,16)
        count=0
        for omega in (1,2,4,8,16):
            lam=F(omega*omega);a=1-tau*tau*lam/2;b=tau;c=-tau*lam*(1-tau*tau*lam/4)
            aa,bb,cc,dd=F(1),F(0),F(0),F(1)
            for n in range(1,17):
                aa,bb,cc,dd=a*aa+b*cc,a*bb+b*dd,c*aa+a*cc,c*bb+a*dd
                self.assertEqual(aa*dd-bb*cc,1)
                t=n*tau;continuous=sin(I.of(t*omega))/omega
                displacement_error=(I.of(bb)-continuous).abs()
                classical_bound=tau*tau*(F(omega,6)+t*lam/18)
                self.assertLessEqual(F(displacement_error.hi,self.interval.SCALE),classical_bound)
                variance=aa*aa/F(2*omega)+bb*bb*F(omega,2)
                excess=variance-F(1,2*omega)
                self.assertEqual(excess,tau*tau*omega**3*bb*bb/8)
                self.assertGreaterEqual(excess,0);self.assertLessEqual(excess,tau*tau*omega/6)
                split_p=(1+exp(I.of(-variance/2))*sin(I.of(bb)))/2
                exact_p=(1+exp(I.of(F(-1,4*omega)))*sin(continuous))/2
                difference=(split_p-exact_p).abs()
                effect_bound=tau*tau*(F(omega,8)+t*lam/36)
                self.assertLessEqual(F(difference.hi,self.interval.SCALE),effect_bound)
                count+=1
        self.assertEqual(count,80)
    def test_grid_points_remain_inside_closed_window(self):
        for row in self.packet['levels']:
            tau=F(row['time_step']);first=row['first_window_step'];last=row['last_window_step']
            self.assertGreaterEqual(first*tau,F(19,20));self.assertEqual(last*tau,1)
            for t in (F(19,20),F(19,20)+tau/10,F(97,100),F(1)-tau/3,F(1)):
                shifted=t/tau+F(1,2)
                n=min(last,max(first,shifted.numerator//shifted.denominator))
                self.assertGreaterEqual(n,first);self.assertLessEqual(n,last)
                self.assertGreaterEqual(n*tau,F(19,20));self.assertLessEqual(n*tau,1)
                self.assertLessEqual(abs(n*tau-t),tau)
                for alternate in (max(first,n-1),min(last,n+1)):
                    self.assertLessEqual(abs(n*tau-t),abs(alternate*tau-t))
    def test_source_gaps_not_float_or_boolean_inputs(self):
        for q in (True,13.0,'13',0,1001):
            with self.subTest(value=q),self.assertRaises(ValueError):V.source_gaps(q)
    def test_coherent_formula_metadata_and_outcome_mutations(self):
        examples=[]
        packet=deepcopy(self.packet)
        for row in packet['levels']:
            row['vacuum_probability_error_upper']='0';row['temporal_probability_error_upper']='0'
            row['same_grid_total_error_upper']=row['spatial_probability_error_upper']
        examples.append(('erased_squeezing',packet))
        packet=deepcopy(self.packet);packet['levels'][-1]['time_step']='1/16';examples.append(('unstable_step',packet))
        packet=deepcopy(self.packet);packet['levels'][-1]['minimum_gap_Qphi']=[89,-55];examples.append(('larger_false_gap',packet))
        packet=deepcopy(self.packet);packet['preparation']['initial_mean_momentum']='q5 polynomial v';examples.append(('different_preparation',packet))
        packet=deepcopy(self.packet);packet['levels'][0]['comparison_error_smaller_than_response']=True;examples.append(('forced_coarse_success',packet))
        packet=deepcopy(self.packet);packet['levels'][-1]['first_window_step']-=1;examples.append(('outside_window',packet))
        packet=deepcopy(self.packet);packet['levels'][-1]['nearby_continuum_error_upper']=packet['levels'][-1]['same_grid_total_error_upper'];examples.append(('erased_time_alignment',packet))
        packet=deepcopy(self.packet);packet['scope']['reference_time_readout']='continuous split state at any time';examples.append(('invented_intermediate_state',packet))
        packet=deepcopy(self.packet);packet['scope']['continuum_probability_Lipschitz_upper']='0';examples.append(('erased_reference_motion',packet))
        for key in ('observer_event_log_executed','full_Fock_state_norm_convergence_claimed',
                    'velocity_norm_convergence_claimed','q5_clock_error_transferred','count_volume_for_substeps_established'):
            packet=deepcopy(self.packet);packet['scope'][key]=True;examples.append((key,packet))
        for name,packet in examples:
            with self.subTest(name=name),self.assertRaises(ValueError):V.verify_arithmetic(packet,self.parent)
    def test_strict_parser(self):
        for data in (b'{"x":1,"x":2}',b'{"x":NaN}',b'{"x":1.0}',b' '*200001):
            with self.subTest(size=len(data)),TemporaryDirectory() as directory:
                path=Path(directory)/'packet.json';path.write_bytes(data)
                with self.assertRaises(ValueError):V.load(path)
    def test_standalone_scientific_copy_and_source_pin_direction(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)/'scientific';copied=set()
            def copy(path):
                if path in copied:return
                copied.add(path);src=V.RER/path;dst=root/path;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
                if src.suffix=='.json':
                    data=json.loads(src.read_bytes())
                    for parent in data.get('source_pins',{}):copy(parent)
            copy(V.PARENT)
            for f in V.FILES:copy(f)
            copy('code/source_scalar_time_refinement/source_scalar_time_refinement_receipt.json')
            result=subprocess.run([sys.executable,str(root/'code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py'),
                                   '--root',str(root)],cwd=root,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertTrue(json.loads(result.stdout)['full_spatial_parent_replayed'])
            with (root/V.FILES[0]).open('a') as f:f.write('\n# altered detached source\n')
            with self.assertRaisesRegex(ValueError,'reconstruction'):V.verify(self.packet,root)

if __name__=='__main__':unittest.main()
