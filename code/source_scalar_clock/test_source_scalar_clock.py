"""Controls exercise actual graph inference, source binding and coherent forgeries."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import shutil, subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
import importlib.util, json, random, sys, unittest

HERE=Path(__file__).resolve().parent
def module(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;exec(compile((HERE/file).read_bytes(),str(HERE/file),'exec'),m.__dict__);return m
V=module('_dag_slice_test_verifier','verify_source_scalar_clock.py')
P=module('_dag_slice_test_producer','source_scalar_clock.py')

def reseal(trace,relabel=False,shuffle=False,refresh_reads=False):
    """Keep genuine writer/value references and the whole hash chain coherent."""
    t=deepcopy(trace);mapping={V.raw(e['id']):'opaque_node_'+str(i*73+19) for i,e in enumerate(t['events'])}
    ledger='0'*64;writes={}
    for e in t['events']:
        if relabel:
            e['id']=mapping[V.raw(e['id'])];e['write'][2]=e['id']
            for r in e['reads']:r[2]=mapping[V.raw(r[2])]
        if refresh_reads:
            for i,r in enumerate(e['reads']):
                if tuple(r[:2]) in writes:e['reads'][i]=deepcopy(writes[tuple(r[:2])])
        e['ledger_parent']=ledger;e['hash']=V.digest({k:e[k] for k in ('id','reads','write','ledger_parent')})
        ledger=e['hash'];writes[tuple(e['write'][:2])]=deepcopy(e['write'])
    t['final_hash']=ledger
    if shuffle:random.Random(918).shuffle(t['events'])
    return t

class DagSlices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.v,cls.parent=V.context();cls.m,cls.A,*_=cls.v.model(V.RER)
        cls.base=cls.parent['traces']['ascending_baseline'];cls.active=cls.parent['traces']['ascending_intervention']
    def recover(self,t):return V.recover(t,self.v,self.m,self.A)
    def test_full_independent_parent_and_clock_join(self):
        result=V.verify(V.load());self.assertEqual(result['events_authenticated'],5888)
        self.assertEqual(result['dynamic_reads_authenticated'],34944)
        self.assertTrue(result['full_parent_mathematical_replay'])
    def test_independent_producer_bytes(self):self.assertEqual(P.raw(P.produce()),V.OUTPUT.read_bytes())
    def test_opaque_id_relabel_storage_shuffle_and_no_adapter_copies(self):
        t=reseal(self.active,relabel=True,shuffle=True)
        for k in ('layers','schedule','seed_writes','update_writes','dynamic_field_reads'):t.pop(k)
        vectors,report=self.recover(t)
        self.assertEqual([[x.encode() for x in row] for row in vectors],self.active['layers'])
        self.assertEqual(report['defined_tick_centers'],20)
    def test_baseline_layer_recovery_does_not_make_clock(self):
        _,report=self.recover(reseal(self.base,relabel=True,shuffle=True))
        self.assertEqual(report['defined_tick_centers'],0);self.assertIsNone(report['recovered_squared_tick_Qphi'])
        self.assertEqual(report['undefined_center_indices'],list(range(21)))
    def test_both_initial_slices_have_zero_height(self):
        g,_,_=V.authenticate(self.base);layers,height,degree=V.infer_layers(g)
        self.assertEqual({height[v] for v in layers[0]|layers[1]},{0})
        self.assertEqual(len(layers[0]|layers[1]),128)
        self.assertEqual(sorted(set(degree.values())),[1,4,5,6,7])
    def test_ambiguous_two_root_scalar_graph_rejected(self):
        with self.assertRaisesRegex(ValueError,'ambiguous'):V.infer_layers({'a':set(),'b':set(),'v':{'a','b'}})
    def test_back_edge_cycle_rejected(self):
        with self.assertRaisesRegex(ValueError,'cycle'):V.infer_layers({'a':{'b'},'b':{'a'}})
    def test_same_layer_actual_read_rejected(self):
        t=deepcopy(self.base);t['events'][129]['reads'].append(deepcopy(t['events'][128]['write']))
        with self.assertRaises(ValueError):self.recover(reseal(t))
    def test_future_actual_writer_rejected_after_complete_rehash(self):
        t=deepcopy(self.base);t['events'][128]['reads'].append(deepcopy(t['events'][192]['write']))
        with self.assertRaisesRegex(ValueError,'backward'):self.recover(reseal(t))
    def test_delete_zero_read_preserving_height_rejected(self):
        t=deepcopy(self.base);t['events'][640]['reads'].pop()
        with self.assertRaisesRegex(ValueError,'stencil'):self.recover(reseal(t))
    def test_wrong_real_zero_writer_with_resealed_chain_rejected(self):
        t=deepcopy(self.base);t['events'][640]['reads'][0]=deepcopy(t['events'][512+1]['write'])
        with self.assertRaises(ValueError):self.recover(reseal(t))
    def test_audit_parent_substitution_rejected(self):
        t=deepcopy(self.base)
        for i,e in enumerate(t['events'][1:],1):e['reads']=[deepcopy(t['events'][i-1]['write'])]
        with self.assertRaises(ValueError):self.recover(reseal(t))
    def test_coherent_value_change_fails_full_action_residual(self):
        t=deepcopy(self.active);t['events'][10*64+17]['write'][3]=['1','0']
        forged=reseal(t,refresh_reads=True);forged.pop('layers');forged.pop('schedule')
        with self.assertRaisesRegex(ValueError,'residual'):self.recover(forged)
    def test_boolean_resource_alias_rejected(self):
        t=deepcopy(self.base);t['events'][64]['write'][0]=False
        with self.assertRaisesRegex(ValueError,'types'):self.recover(reseal(t))
    def test_wrong_layer_commitment_and_scope_rejected(self):
        for variant in ('late_commitment','physical_promotion'):
            forged=V.load()
            if variant=='late_commitment':forged['traces']['ascending_intervention']['layer_value_sha256'][-1]='0'*64
            else:forged['scope']['physical_or_regional_time_slice_claim']=True
            with self.subTest(variant=variant),self.assertRaisesRegex(ValueError,'canonical clock reconstruction'):V.verify(forged)
    def test_strict_loader(self):
        for blob in (b'{"x":1,"x":2}',b'{"x":NaN}',b'{"x":1.0}',b' '*500001):
            with self.subTest(size=len(blob)):
                with TemporaryDirectory() as d:
                    f=Path(d)/'packet.json';f.write_bytes(blob)
                    with self.assertRaises(ValueError):V.load(f)


producer=P
verifier=V
class ExactPerturbationControls(unittest.TestCase):
    def test_genuinely_nonuniform_noiseless_pair(self):
        b = producer.coefficient_box(2, -3, 1, 1, 0, 0, 0)
        self.assertEqual(b['previous_duration_interval'], [F(1), F(1)])
        self.assertEqual(b['next_duration_interval'], [F(2), F(2)])

    def test_729_exact_perturbed_nonuniform_pairs_and_residuals(self):
        # Independent construction: F=diag(1+dx,1+dy), c near (2,-3).
        # Perturb z directly, including a transverse-to-one-column residual.
        # Norm budgets use the rational l1 upper bound of the Euclidean norm.
        e = F(1, 1000)
        for dx, dy, da, db, r1, r2 in product((-e, F(0), e), repeat=6):
            a, b = 2+da, -3+db
            dz = ((1+dx)*a+r1-2, (1+dy)*b+r2+3)
            budget = producer.coefficient_box(2, -3, 1, 1, abs(dx), abs(dy),
                                             abs(dz[0])+abs(dz[1]), abs(r1)+abs(r2))
            self.assertLessEqual(abs(a-2), budget['alpha_radius'])
            self.assertLessEqual(abs(b+3), budget['beta_radius'])
            hm2, hp2 = -2*b/(a*(a+1)), -2*a*b/(a+1)
            for val, key in ((hm2, 'previous_duration_interval'), (hp2, 'next_duration_interval')):
                lo, hi = budget[key]
                self.assertLessEqual(lo*lo, val)
                self.assertLessEqual(val, hi*hi)

    def test_configuration_stationarity_for_unequal_positive_steps(self):
        # x=e1,y=e2; construct z from the original variable-step Euler equation.
        for hm, hp in product((F(1, 2), F(1), F(3, 2), F(2)), repeat=2):
            alpha, beta = hp/hm, -hp*(hp+hm)/2
            self.assertEqual(alpha/hp-1/hm, 0)
            self.assertEqual(beta/hp+(hp+hm)/2, 0)
            bounds = producer.coefficient_box(alpha, beta, 1, 1, 0, 0, 0)
            self.assertEqual(bounds['previous_duration_interval'], [hm, hm])
            self.assertEqual(bounds['next_duration_interval'], [hp, hp])

    def test_rounding_of_exact_and_nonexact_squares(self):
        for value in (F(0), F(1), F(1, 4), F(2), F(7, 31), F(10**8+1, 10**23)):
            lo, hi = producer.sqrt_interval(value)
            self.assertLessEqual(lo*lo, value)
            self.assertLessEqual(value, hi*hi)
            self.assertLessEqual(hi-lo, F(1, 10**12))

    def test_conditioning_and_sign_failure_are_not_certified(self):
        cases = [(1, -1, 1, 1, 1, 0, 0), (1, -1, 1, 1, 0, 0, 2),
                 (1, 0, 1, 1, 0, 0, 0), (0, -1, 1, 1, 0, 0, 0),
                 (1, -1, 0, 1, 0, 0, 0), (1, -1, 1, 1, -1, 0, 0)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                producer.coefficient_box(*args)
        for bad in (True, 0.01, '1/10', complex(1, 0)):
            with self.subTest(type=type(bad).__name__), self.assertRaises(ValueError):
                producer.coefficient_box(1, -1, 1, 1, bad, 0, 0)

    def test_nearly_degenerate_pair_requires_smaller_budget(self):
        # x=e1,y=e1+d e2 gives l_alpha=sqrt(1+d^-2), l_beta=d^-1.
        d = F(1, 1000)
        la, lb = F(1001), 1/d  # rigorous norm upper bounds
        with self.assertRaisesRegex(ValueError, 'conditioning'):
            producer.coefficient_box(1, -1, la, lb, d, d, 0)
        self.assertGreater(producer.coefficient_box(1, -1, la, lb, d/100, d/100, 0)['q'], 0)

    def test_overlapping_edges_are_intersections_not_unions(self):
        rows = [{'previous_duration_interval': [F(1), F(3)],
                 'next_duration_interval': [F(2), F(4)]} for _ in range(20)]
        edges = producer.intersect_edges(rows)
        self.assertEqual(edges[1:20], [[F(2), F(3)]]*19)
        rows[3]['previous_duration_interval'] = [F(5), F(6)]
        with self.assertRaisesRegex(ValueError, 'intersection'):
            producer.intersect_edges(rows)



class CanonicalClockControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.v,cls.parent=V.context();cls.packet=V.load()
    def test_all40_gram_and_stability_rows(self):
        report=V.verify_arithmetic(self.packet,self.v,self.parent)
        self.assertEqual(report['positive_gram_cases'],40)
        self.assertFalse(report['full_parent_mathematical_replay'])
        self.assertEqual(report['final_elapsed_interval'],['1341635903357/1000000000000','670822834819/500000000000'])
    def test_coherent_stability_and_scope_mutations(self):
        changes=[]
        p=deepcopy(self.packet)
        for row in p['stability']['records']:row['alpha_radius']=row['beta_radius']='0'
        p['stability']['summary']['maximum_alpha_radius']=p['stability']['summary']['maximum_beta_radius']='0'
        changes.append(('coherent_zero_error',p))
        p=deepcopy(self.packet);p['stability']['records'][0]['left_inverse_row_norm_upper']=['1','1'];changes.append(('wrong_Gram_norm',p))
        p=deepcopy(self.packet)
        for schedule in p['stability']['schedules'].values():schedule['elapsed_time_intervals'][-1]=['1','1']
        p['stability']['summary']['final_elapsed_interval']=['1','1'];p['stability']['summary']['final_elapsed_width']='0'
        changes.append(('coherent_final_time',p))
        p=deepcopy(self.packet);p['exact_duration']['positive_tick_Qphi']=['1','0'];changes.append(('wrong_exact_tick',p))
        p=deepcopy(self.packet);p['stability']['records'][0]['center_step']=True;changes.append(('boolean_center',p))
        p=deepcopy(self.packet);p['stability']['budget']['configuration_M_norm']='1/1000000';changes.append(('noise_budget',p))
        for key in ('arbitrary_site_dependent_lapse','interval_overlap_implies_history_existence',
                    'noisy_history_existence_certified','increment_residual_is_unscaled_EL_residual'):
            p=deepcopy(self.packet);p['scope'][key]=True;changes.append((key,p))
        for name,packet in changes:
            with self.subTest(name=name),self.assertRaises(ValueError):V.verify_arithmetic(packet,self.v,self.parent)
    def test_all_edges_and_elapsed_sums_contain_exact_duration(self):
        v=self.v
        for schedule in self.packet['stability']['schedules'].values():
            for lo,hi in schedule['edge_duration_intervals']:
                self.assertLessEqual((v.R(F(lo))-v.TAU).sign(),0);self.assertGreaterEqual((v.R(F(hi))-v.TAU).sign(),0)
            for j,(lo,hi) in enumerate(schedule['elapsed_time_intervals']):
                self.assertLessEqual((v.R(F(lo))-j*v.TAU).sign(),0);self.assertGreaterEqual((v.R(F(hi))-j*v.TAU).sign(),0)
    def test_full_vector_not_just_gram_projection(self):
        v=self.v;mass=[v.R(1)]*3;A=[[(i,v.R(1))] for i in range(3)]
        triple=[[v.R(-1),v.R(1),v.R()],[v.R(),v.R(1),v.R()],
                [v.R(1),v.R(F(244,245)),v.R(1)]]
        with self.assertRaisesRegex(ValueError,'full vector equation'):V.expected_row(v,mass,A,triple,'fixture',1)
    def test_nonuniform_positive_and_degenerate_controls(self):
        v=self.v
        out=P.infer(v,[v.R(1),v.R(1)],[v.R(1),v.R()],[v.R(),v.R(1)],[v.R(2),v.R(-3)])
        self.assertEqual(v.parse(out['previous_squared_duration_Qphi']),v.R(1))
        self.assertEqual(v.parse(out['next_squared_duration_Qphi']),v.R(4))
        with self.assertRaisesRegex(ValueError,'Gram'):
            P.infer(v,[v.R(1)],[v.R(1)],[v.R(1)],[v.R(0)])
        for hm,hp in [(F(1),F(1)),(F(1,2),F(7,2))]:self.assertEqual(-1/hm+(hm+hp)/2,0)
    def test_standalone_root_full_replay_and_rehashed_source_forgery(self):
        # Copy only the recursively pinned scientific inputs, without workspace
        # siblings, Git state, outside prototype files or installed package code.
        with TemporaryDirectory() as directory:
            root=Path(directory)/'standalone';copied=set()
            def copy(relative):
                if relative in copied:return
                source=V.RER/relative;target=root/relative;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,target);copied.add(relative)
                if source.suffix=='.json':
                    data=json.loads(source.read_bytes())
                    if isinstance(data,dict):
                        for name in data.get('source_pins',{}):
                            if (V.RER/name).is_file():copy(name)
            copy(V.PARENT);copy(V.VERIFIER)
            for f in V.FILES:copy(f)
            copy('code/source_scalar_clock/source_scalar_clock_receipt.json')
            result=subprocess.run([sys.executable,str(root/'code/source_scalar_clock/verify_source_scalar_clock.py'),'--root',str(root)],cwd=root,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertTrue(json.loads(result.stdout)['full_parent_mathematical_replay'])
            # An original-checkout fallback would incorrectly accept this.
            with (root/V.FILES[0]).open('a') as f:f.write('\n# altered detached source\n')
            with self.assertRaisesRegex(ValueError,'reconstruction'):V.verify(self.packet,root=root)

if __name__=='__main__':unittest.main()
