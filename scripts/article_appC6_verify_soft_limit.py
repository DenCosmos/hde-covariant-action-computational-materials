#!/usr/bin/env python3
"""V52-R1: article appC6 verify soft limit.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/verify_soft_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC6_verify_soft_limit.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article 6.10, eq:v34r-soft-scalar-vertex (195).
  article 6.10, eq:v33-soft-singularity (196).
  article C.6, eq:boundary-five-groups-leading (393).
  article C.6, eq:boundary-scalar-five-coefficients (394).
Method: Exact Laurent soft coefficients and five boundary leading groups.
Inputs: New radial pairs and contact records; r=1+mu*q before q->0.
Outputs below the selected results root: soft_v40/summary.json; soft_v40/checks.csv; soft_v40/five_groups.txt.
Provenance: GitHub_bundle/original_scripts/verify_soft_v40.py.
Scope limit: Does not prove existence of an unregulated unitary limit. Equal radii are not substituted before the joint soft limit.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import json,csv,itertools
import sympy as sp
import radial_coeff_algebra as a
from supplement_sec1_2_assemble_contacts import decode_file,load
from article_appC5_compute_angular_coefficients import transfer,f_boundary
B=Path(__file__).resolve().parents[1];D=RESULTS/'radial_generic';OUT=RESULTS/'soft_v40';R,Q=a.F.gens
CHECK=[];EXP=[]
def joint(v):
    def poly(p):return sum((co*(1+R*Q)**mon[0]*Q**mon[1] for mon,co in p.items()),a.F.zero)
    return a.Coeff({ph:poly(rat.numer)/poly(rat.denom) for ph,rat in v.d.items()})
def valuation(v):
    for rat in v.d.values():
        low=min(m[1] for m in rat.denom)
        assert all(m[0]==0 for m in rat.denom if m[1]==low), 'Nonconstant angular leading denominator'
    return min((min(m[1] for m in rat.numer)-min(m[1] for m in rat.denom) for rat in v.d.values()),default=999)
def check(name,test,detail):
    if not test:raise AssertionError(name+': '+str(detail))
    CHECK.append({'id':name,'status':'PASS','detail':str(detail)})
# V52: article eq:v34r-soft-scalar-vertex (195); article eq:v33-soft-singularity (196); article eq:boundary-five-groups-leading (393); article eq:boundary-scalar-five-coefficients (394).
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for sa,sb in itertools.product(range(3),repeat=2):
        key=f'02{sa}{sb}';d=decode_file(D/'pair'/(key+'.json'))
        f=transfer(f_boundary(0,2,sa,sb));cs=lambda s:a.O/a.SQ3 if s==0 else a.O
        om=cs(sa)-a.K.gens[0]*cs(sb)
        regA=a.la_add({p:transfer(v) for p,v in d['A'][0].items()},{-2:f,-1:a.I*om*f})
        regB=a.la_add({p:transfer(v) for p,v in d['B'][0].items()},{-1:-f})
        for name,vals in [('Areg',regA),('Breg',regB)]:
            for p,v in vals.items():
                order=valuation(joint(v));check(f'{key}.{name}.x{p}',order>=0,order)
        for I in [1,2]:
            for name in ['A','B']:
                for p,v in d[name][I].items():
                    val=joint(transfer(v));order=valuation(val)-(1 if I==2 else 0)
                    check(f'{key}.{name}.tensor{I}.x{p}',order>=0,order)
    for path in (D/'contacts').glob('*.json'):
        for p,v in load(json.loads(path.read_text())['full_H4']).items():
            # A full contact contains both transfer denominators. After mapping
            # the plus channel the opposite denominator stays regular at q=0.
            order=valuation(joint(transfer(v)))
            check(f'contact.{path.stem}.x{p}',order>=-2,order)
    ai,af,bi,bf=sp.symbols('a_i a_f b_i b_f');S=(af-ai)*(bf-bi)
    groups=[af*(bf-bi),-bi*(af-ai),af*bi,-af*bf/2,-ai*bi/2]
    reverse=[bf*(af-ai),-ai*(bf-bi),bf*ai,-bf*af/2,-bi*ai/2]
    check('five_boundary_groups_plus_reverse',sp.expand(sum(groups)+sum(reverse)-S)==0,0)
    qi,qf,q,mu=sp.symbols('x_i x_f q mu',positive=True);f=4*sp.sqrt(3)*sp.I*mu/q
    raw=24*sp.sqrt(3)*mu**2/q**3*(1/qi-1/qf)**2
    correction=(sum(groups)+sum(reverse)).subs({ai:f/qi,af:f/qf,bi:f/qi,bf:f/qf})/(2*q/sp.sqrt(3))
    check('scalar_q_minus3_cancellation',sp.simplify(sp.expand(raw+correction))==0,0)
    for n,(g,h) in enumerate(zip(groups,reverse),1):
        val=sp.factor((g+h).subs({ai:f/qi,af:f/qf,bi:f/qi,bf:f/qf})/(2*q/sp.sqrt(3)))
        EXP.append(f'boundary_group_{n} = {val}')
    with (OUT/'checks.csv').open('w',newline='') as out:w=csv.DictWriter(out,fieldnames=['id','status','detail']);w.writeheader();w.writerows(CHECK)
    (OUT/'five_groups.txt').write_text('\n'.join(EXP)+'\n')
    (OUT/'summary.json').write_text(json.dumps({'success':True,'checks':len(CHECK),'joint_limit':'r=1+mu*q, q->0; all nine hard pairs and all 81 contacts'},indent=2))
    print('PASS:',len(CHECK),'exact soft-limit checks')
if __name__=='__main__':main()
