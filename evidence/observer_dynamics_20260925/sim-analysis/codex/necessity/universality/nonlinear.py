"""Exact density-projection witnesses for nonlinear native observables."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import certify as algebra

HERE=Path(__file__).resolve().parent


def build():
    spec=json.loads((HERE/"spec.json").read_text())
    row=next(g for g in spec["graphs"] if g["name"]=="path4")
    rule,pin=algebra.production_rule(spec["source_commit"])
    states,q,lap,incidence,z,kappa,pi=algebra.configuration_data(
        row["ports"],row["raised"],[(a,b,F(w)) for a,b,w in row["edges"]],rule)
    size=len(states)
    product=[state[0]*state[1] for state in z]
    mean=sum(product)/size
    even=[[value-mean] for value in product]
    linear=[[state[0]] for state in z]
    observations={"load":linear,"even_pair":even,"mixed":algebra.add(linear,even)}
    gp=algebra.pseudoinverse_connected_laplacian(algebra.scale(F(-1),q))
    lp=algebra.pseudoinverse_connected_laplacian(lap)
    rows={}
    for name,f in observations.items():
        b=algebra.scale(F(1)/(size*kappa),algebra.multiply(algebra.transpose(f),z))
        projected=algebra.multiply(z,algebra.transpose(b))
        g=algebra.add(f,algebra.scale(F(-1),projected))
        assert algebra.multiply(algebra.transpose(g),z)==algebra.zeros(1,len(z[0]))
        assert algebra.multiply(b,pi)==b
        cross=algebra.multiply(algebra.multiply(algebra.transpose(projected),gp),g)
        assert cross==[[F(0)]]
        total=algebra.scale(F(2,size),algebra.congruence(algebra.transpose(f),gp))
        density=algebra.scale(4*kappa,algebra.congruence(b,lp))
        residual=algebra.scale(F(2,size),algebra.congruence(algebra.transpose(g),gp))
        assert total==algebra.add(density,residual)
        rows[name]={"observable":algebra.strings(f),"B":algebra.strings(b),
                    "orthogonal_residual":algebra.strings(g),
                    "total_asymptotic_variance":str(total[0][0]),
                    "density_asymptotic_variance":str(density[0][0]),
                    "residual_asymptotic_variance":str(residual[0][0]),
                    "cross_poisson_covariance":"0","all_projection_checks_exact":True}
    assert rows["even_pair"]["B"]==[["0"]*4]
    assert F(rows["even_pair"]["total_asymptotic_variance"])>0
    assert rows["mixed"]["density_asymptotic_variance"]==rows["load"]["total_asymptotic_variance"]
    return {"schema":"oph.native-nonlinear-density-projection.receipt.v1",
            "inputs_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                ["nonlinear_spec.json","nonlinear.py","certify.py","spec.json"]},
            "production_rule_source":pin,"observables":rows,
            "scope":json.loads((HERE/"nonlinear_spec.json").read_text())["scope"]}


if __name__=="__main__":
    result=build()
    (HERE/"nonlinear_receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("Exact nonlinear density-projection witnesses passed")
