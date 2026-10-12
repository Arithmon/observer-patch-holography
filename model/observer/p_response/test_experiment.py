"""Independent charge, Ramsey instrument and complete-history checks."""
from decimal import Decimal as D
from fractions import Fraction as Q
import cmath
from itertools import combinations, product
import json
import math
from pathlib import Path
import unittest

from model.observer.p_response import experiment as subject


def endpoint(loads, seam, coin):
    a, b = seam
    total = loads[a] + loads[b]
    low, high = total//2, (total+1)//2
    out = list(loads)
    out[a], out[b] = (low, high) if coin == 0 else (high, low)
    return tuple(out), (coin, loads[a]-out[a])


class ResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = json.loads(Path(__file__).with_name("results.json").read_text())
        _, cls.native = subject.sources()
        cls.states = cls.native.states()
        cls.g = float(cls.r["source_branch"]["g_sqrt_4pi_alpha"])

    def test_p_to_charge_and_low_energy_branch(self):
        branch = self.r["source_branch"]
        self.assertEqual(branch["mode"], "thomson_structured_running")
        self.assertAlmostEqual(self.g**2, 4*math.pi*float(branch["alpha"]), places=14)
        self.assertAlmostEqual(float(branch["inverse_alpha"]), 1/float(branch["alpha"]), places=11)
        self.assertAlmostEqual(float(branch["inverse_alpha"]),
                               float(branch["mz_anchor_inverse"])+float(branch["low_energy_inverse_shift"]), places=11)
        self.assertLess(abs(float(branch["outer_closure_residual"])), 1e-13)

    def test_independent_complex_ramsey_kraus_and_record_preservation(self):
        for j in (-1,0,1):
            for exposure in (0.0,0.3,1.0,2.0,-1.0):
                phi = self.g*j*exposure
                amps = [(cmath.exp(-0.5j*phi)-1j*y*cmath.exp(0.5j*phi))/2 for y in (-1,1)]
                self.assertAlmostEqual(sum(abs(a)**2 for a in amps),1,places=14)
                for y, amp in zip((-1,1),amps):
                    expected = float(subject.adapter_probability(j,y,D(str(self.g)),D(str(exposure))))
                    self.assertAlmostEqual(abs(amp)**2,expected,places=14)
        # Each Kraus is diagonal in old current; on a generic coherent input
        # its forgotten-outcome action leaves every diagonal population fixed.
        a = [complex(1,2),complex(2,-1),complex(-1,3)]
        for k,j in enumerate((-1,0,1)):
            mass = 0
            for y in (-1,1):
                phi = self.g*j
                m = (cmath.exp(-0.5j*phi)-1j*y*cmath.exp(0.5j*phi))/2
                mass += abs(m*a[k])**2
            self.assertAlmostEqual(mass,abs(a[k])**2,places=13)

    def test_exact_native_feedback_coefficients(self):
        # Symbolic affine probabilities: p(y|j)=(1+y*j*s)/2 for j=0,+/-1.
        mean = [Q(0),Q(0)]
        corr = [Q(0),Q(0)]
        bits = [Q(0),Q(0)]
        words = set()
        for occupied in combinations(range(12),2):
            x = tuple(int(k in occupied) for k in range(12))
            for c1 in (0,1):
                first,r1 = endpoint(x,(0,1),c1)
                for bit in (-1,1):
                    seam = (0,5) if bit==1 else (0,1)
                    for c2 in (0,1):
                        final,r2 = endpoint(first,seam,c2)
                        words.add((r1,bit,r2))
                        for k,w in enumerate((Q(1,528),Q(bit*r1[1],528))):
                            mean[k] += w*r2[1]
                            corr[k] += w*bit*r1[1]
                            bits[k] += w*bit
        self.assertEqual(mean,[Q(0),Q(5,264)])
        self.assertEqual(corr,[Q(0),Q(5,33)])
        self.assertEqual(bits,[Q(0),Q(0)])
        self.assertEqual(len(words),24)
        self.assertAlmostEqual(float(self.r["response"]["mean_second_native_current"]),
                               (5/264)*math.sin(self.g),places=14)

    def test_all_528_rows_independent_native_oracle(self):
        self.assertEqual(len(self.r["histories"]),528)
        mass = 0.0
        for row in self.r["histories"]:
            x = self.states[row["initial_source"]]
            z,r1 = endpoint(x,(0,1),row["first_native_record"][0])
            self.assertEqual(r1,tuple(row["first_native_record"]))
            bit = row["detector_bit"]
            seam = (0,5) if bit==1 else (0,1)
            final,r2 = endpoint(z,seam,row["second_native_record"][0])
            self.assertEqual(final,self.states[row["final_source"]])
            self.assertEqual(r2,tuple(row["second_native_record"]))
            probability = (1+bit*math.sin(self.g*r1[1]))/528
            self.assertAlmostEqual(float(row["history_probability"]),probability,places=16)
            mass += probability
        self.assertAlmostEqual(mass,1,places=13)

    def test_bit_only_actuator_preserves_bit_all_source_inputs(self):
        for x in self.states:
            for bit,coin in product((-1,1),(0,1)):
                final,retained,r = subject.controlled_native_step(self.native,x,bit,coin)
                expected,expected_r = endpoint(x,(0,5) if bit==1 else (0,1),coin)
                self.assertEqual((final,retained,r),(expected,bit,expected_r))
        # There is no current/first-record argument to either actuator interface.
        import inspect
        self.assertEqual(tuple(inspect.signature(subject.decoder_seam).parameters),("detector_bit",))
        self.assertEqual(tuple(inspect.signature(subject.controlled_native_step).parameters),
                         ("native","loads","detector_bit","coin"))

    def test_sample_and_sensitivity_controls(self):
        samples = self.r["sampled_execution"]
        self.assertEqual(samples["seed"],20261012)
        self.assertEqual(samples["trials"],200000)
        self.assertEqual(sum(row["count"] for row in samples["terminal_word_counts"]),samples["trials"])
        expected = float(self.r["response"]["mean_second_native_current"])
        self.assertLess(abs(samples["mean_second_native_current"]-expected),6*samples["standard_error_of_mean"])
        variants = self.r["sensitivity_controls"]
        self.assertNotEqual(variants[0]["detector_plus_given_j_plus"],variants[1]["detector_plus_given_j_plus"])
        for variant in variants:
            self.assertAlmostEqual(float(variant["final_occupation_0"]),1/6-5*math.sin(float(variant["g"]))/264,places=14)


if __name__ == "__main__":
    unittest.main()
