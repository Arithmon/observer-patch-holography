"""False-green controls for all-time rigorous continuous-pressure bounds."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('continuous_pressure_independent',HERE/'verify_continuous_pressure.py')
verifier=importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


class ContinuousPressureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet=json.loads((HERE/'continuous_pressure_receipt.json').read_text())

    def reject(self,mutate):
        changed=copy.deepcopy(self.packet)
        mutate(changed)
        result=verifier.verify(changed)
        self.assertFalse(result['receipt'],result)

    def test_complete_independent_replay(self):
        result=verifier.verify(self.packet)
        self.assertTrue(result['receipt'],result)
        self.assertTrue(result['rigorous_pressure_enclosures'])

    def test_omitted_time(self):
        self.reject(lambda r:r['samples'].pop())

    def test_native_promotion(self):
        self.reject(lambda r:r['scope'].update(native_repair_eos=True))

    def test_split_records_cannot_be_relabelled_continuous(self):
        self.reject(lambda r:r['scope'].update(recorded_split_values_identified_with_exact_continuous_values=True))

    def test_equilibrium_promotion(self):
        self.reject(lambda r:r['scope'].update(equilibrium_thermodynamics=True))

    def test_taylor_tail_cannot_be_removed(self):
        self.reject(lambda r:r['bounds'].update(v_tail='0',pressure_tail='0'))

    def test_later_time_tail_cannot_be_omitted(self):
        self.reject(lambda r:r['samples'][-1].update(pressure_tail_bound='0'))

    def test_regulator_cannot_drop_modes(self):
        self.reject(lambda r:r['polynomial'].update(all_operator_modes_retained=1))

    def test_false_pressure_interval(self):
        self.reject(lambda r:r['samples'][0].update(continuous_pressure_interval=['0','0']))

    def test_changed_velocity_polynomial(self):
        self.reject(lambda r:r['samples'][0]['polynomial_v_Qphi'][0].__setitem__(0,'99'))

    def test_parent_source_drift(self):
        self.reject(lambda r:r['source_sha256'].__setitem__(verifier.PARENT,'0'*64))

    def test_nonpositive_spectral_lower_bound(self):
        self.reject(lambda r:r['bounds'].update(operator_lower='0'))


if __name__=='__main__':
    unittest.main()
