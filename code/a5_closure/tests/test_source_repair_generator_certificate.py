import copy
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P))
import source_repair_generator_certificate as cert
import verify_source_repair_generator_independent as independent

class GeneratorBoundaryTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.carrier=json.loads(cert.MANIFEST.read_text())
  cls.receipts={"a3":json.loads(cert.A3_PATH.read_text()),"directed":json.loads(cert.DIRECTED_PATH.read_text()),"repair":json.loads(cert.REPAIR_PATH.read_text())}

 def test_independent_graph_orbits_invariant_dimension_and_generator_pair(self):
  p=cert.classify();v=independent.audit()
  self.assertEqual((len(p["edges"]),len(p["group"]),p["invariant_dimension"]),(30,60,1))
  self.assertEqual(v["edge_orbit_sizes"],[30])
  self.assertTrue(v["uniform_and_biased_laws"]["nonproportional"])

 def test_conservation_and_exact_seam_support(self):
  e=cert.classify()["edges"];u=[Fraction(1)]*30;b=[Fraction(2 if i==0 else 1) for i in range(30)]
  Lu=cert.laplacian(12,e,u);Lb=cert.laplacian(12,e,b)
  self.assertTrue(all(sum(r)==0 for L in (Lu,Lb) for r in L))
  supp=lambda L:{(i,j) for i in range(12) for j in range(12) if i!=j and L[i][j]}
  expected={(a,b) for a,b in e}|{(b,a) for a,b in e}
  self.assertEqual(supp(Lu),expected);self.assertEqual(supp(Lb),expected)
  self.assertFalse(independent.proportional(Lb,Lu))

 def test_deleted_face_duplicate_face_and_reversed_face_are_rejected(self):
  m=copy.deepcopy(self.carrier);m["carrier"]["oriented_faces"].pop()
  with self.assertRaises(ValueError):cert.classify(m)
  m=copy.deepcopy(self.carrier);m["carrier"]["oriented_faces"][1]=m["carrier"]["oriented_faces"][0][:]
  with self.assertRaises(ValueError):cert.classify(m)
  m=copy.deepcopy(self.carrier);a,b,c=m["carrier"]["oriented_faces"][0];m["carrier"]["oriented_faces"][0]=[a,c,b]
  with self.assertRaises(ValueError):cert.classify(m)

 def test_deleted_duplicate_and_distance_two_serialized_edges_are_rejected(self):
  m=copy.deepcopy(self.carrier);m["carrier"]["edges"].pop()
  with self.assertRaises(ValueError):cert.classify(m)
  m=copy.deepcopy(self.carrier);m["carrier"]["edges"].append(m["carrier"]["edges"][0][:])
  with self.assertRaises(ValueError):cert.classify(m)
  m=copy.deepcopy(self.carrier);edges=m["carrier"]["edges"]
  adj={x:set() for x in m["carrier"]["ports"]}
  for a,b in edges:adj[a].add(b);adj[b].add(a)
  d2=next((a,b) for a in m["carrier"]["ports"] for b in m["carrier"]["ports"] if a<b and b not in adj[a] and adj[a]&adj[b])
  edges[-1]=list(d2)
  with self.assertRaises(ValueError):cert.classify(m)

 def test_receipt_pins_self_digests_and_scope_fail_closed(self):
  self.assertEqual(len(cert.read_receipts()["receipts"]),3)
  self.assertEqual(len(independent.verify_receipts()),3)
  m=copy.deepcopy(self.receipts);m["directed"]["uniform_s1_channel"]["sector"]="all integer states"
  self._redigest(m["directed"])
  with self.assertRaises(ValueError):cert.read_receipts(m)
  with self.assertRaises(AssertionError):independent.verify_receipts(m)
  m=copy.deepcopy(self.receipts);m["a3"]["a3_selection"]["selected_kernel_probability_per_seam"]="1/29"
  self._redigest(m["a3"])
  with self.assertRaises(ValueError):cert.read_receipts(m)
  with self.assertRaises(AssertionError):independent.verify_receipts(m)
  m=copy.deepcopy(self.receipts);m["directed"]["verdict"]["full_self_readback_and_universe_selection"]="exact"
  self._redigest(m["directed"])
  with self.assertRaises(ValueError):cert.read_receipts(m)
  with self.assertRaises(AssertionError):independent.verify_receipts(m)

 def _redigest(self,doc):
  doc["manifest_sha256"]=cert.canonical_digest({k:v for k,v in doc.items() if k!="manifest_sha256"})

if __name__=="__main__":unittest.main()
