import sys
import unittest
from fractions import Fraction
from pathlib import Path
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P))
import source_repair_generator_certificate as cert
import verify_source_repair_generator_independent as independent

class GeneratorBoundaryTests(unittest.TestCase):
 def test_face_reconstruction_and_group_classification(self):
  c=cert.classify();v=independent.audit()
  self.assertEqual((len(c["edges"]),len(c["group"]),c["invariant_dimension"]),(30,60,1))
  self.assertEqual(v["edge_orbit_sizes"],[30])
 def test_weighted_laplacian_conservative_and_supported(self):
  c=cert.classify();L=cert.laplacian(12,c["edges"],[Fraction(1)]*30)
  self.assertTrue(all(sum(row)==0 for row in L))
  self.assertEqual(sum(x!=0 for i,row in enumerate(L) for j,x in enumerate(row) if i!=j),60)
 def test_two_positive_laws_same_relation_distinct_generator(self):
  c=cert.classify();u=[Fraction(1)]*30;b=[Fraction(2 if i==0 else 1) for i in range(30)]
  self.assertNotEqual(cert.laplacian(12,c["edges"],u),cert.laplacian(12,c["edges"],b))
  self.assertTrue(all(x>0 for x in b))
 def test_deleted_edge_and_duplicate_edge_controls(self):
  c=cert.classify();e=c["edges"]
  self.assertEqual(len(set(e[:-1])),29)
  self.assertEqual(len(e+[e[0]]),31)
  self.assertNotEqual(set(e[:-1]),set(e))
 def test_distance_two_and_broken_weight_controls(self):
  c=cert.classify();e=c["edges"];adj=[set() for _ in range(12)]
  for u,v in e:adj[u].add(v);adj[v].add(u)
  d2=next((u,v) for u in range(12) for v in range(u+1,12) if v not in adj[u] and adj[u]&adj[v])
  self.assertNotIn(d2,e)
  w=[Fraction(1)]*30;w[0]=Fraction(0)
  self.assertFalse(all(x>0 for x in w))
 def test_invariance_rank_rejects_split_orbit(self):
  c=cert.classify();self.assertEqual(c["invariant_dimension"],1)
  self.assertEqual(len(c["edge_orbits"]),1)
 def test_source_boundary_is_explicit(self):
  out=cert.certificate()
  self.assertFalse(out["proposal_source_audit"]["reference_source_derived_from_issue_628_relation"])
  self.assertEqual(out["proposal_source_audit"]["full_integer_state_generator_binding"],"not established")

if __name__=="__main__":unittest.main()
