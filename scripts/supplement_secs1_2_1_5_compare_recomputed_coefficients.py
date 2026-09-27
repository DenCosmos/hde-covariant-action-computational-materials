#!/usr/bin/env python3
"""Compare regenerated ADM coefficients with the unchanged V52 supplement.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics, V52; release V52-R1.
Supplement Sections 1.2, 1.3 and 1.5.2: contact representatives, 54 A/B pairs,
46 explicit crossed-constraint coefficients, homogeneous constraint sources.
The 81 full contacts are also compared with the PRINTED reconstruction map.
Method: parse printed R functions, restore t^2=(1-z)/(1+z), compare exact
expressions in Q(r,z,i,sqrt(3),sqrt(1-z^2)). Assumptions r>0,-1<z<1.
Inputs: unchanged mathematical_catalogue.tex and NEW radial_generic outputs.
Outputs: printed_vs_recomputed/checks.json; each test identifies source labels.
Run: python scripts/supplement_secs1_2_1_5_compare_recomputed_coefficients.py
Dependencies: SymPy and local parser/algebra modules. No TeX or network.
Provenance: NEW comparison routine; it does not generate or edit the catalogue.
The printed catalogue is the comparison object, not the generator's input.
This validates coefficient functions, not arbitrary-background scattering,
state choice, a regulator-free unitary limit or every analytical equation.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  supplement 1.2, sec:cat-contacts (1.2).
  supplement 1.3, sec:cat-pairs (1.3).
  supplement 1.2.14, sec:cat-component-table (1.2.14).
  supplement 1.5.2, sec:cat-contact-partitions (1.5.2).
  supplement 1.5.2, eq:cat-partition-homogeneous-polynomials (501).
  supplement 1.5.2, eq:cat-partition-auxiliary-minus (503).
Method: Exact restored r,z functions compared with printed R maps, reconstruction, 46 crossed auxiliary terms and homogeneous lapse polynomials.
Inputs: Unchanged V52 supplementary TeX; independently generated 54 pairs and 81 contacts.
Outputs below the selected results root: printed_vs_recomputed/checks.json.
Provenance: NEW parser/comparator; preserved algebra and pure half-angle transform reused.
Scope limit: Catalogue is never regenerated; tests carry exact per-coefficient source labels. No fictitious intermediate data.
"""
from __future__ import annotations
import argparse,json,re,sys,time
from pathlib import Path
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
import radial_coeff_algebra as a
from half_angle_transform import to_z
from release_paths import RESULTS
from supplement_sec1_2_assemble_contacts import load,decode_file
import supplement_sec1_verify_printed_catalogue as cat
ROOT=Path(__file__).resolve().parents[1]

def restored(v):
    even,odd=to_z(v)
    return (even.expr()+sp.sqrt(1-a.HALF_ANGLE**2)*odd.expr()).subs({a.RADIUS:cat.r,a.HALF_ANGLE:cat.z})

# V52: supplement sec:cat-contacts (1.2); supplement sec:cat-pairs (1.3).
def read_printed(path):
    text=path.read_text();eqs=[m[2] for m in re.finditer(r'\\begin\{(equation|align)\}([\s\S]*?)\\end\{\1\}',text)]
    bylabel={m[1]:e.split(r'\label')[0].strip() for e in eqs if (m:=re.search(r'\\label\{([^}]+)\}',e))};R={};N={};D={}
    for label,eq in bylabel.items():
        m=re.fullmatch(r'eq:cat-(rational|numerator|denominator)-(\d+)',label)
        if not m:continue
        kind,i=m[1],int(m[2]);rhs=eq.split('=',1)[1]
        if kind=='numerator':N[i]=cat.math(rhs)
        elif kind=='denominator':D[i]=cat.math(rhs)
        elif 'N_{' not in rhs:R[i]=cat.math(rhs)
    for i in N:R[i]=N[i]/D[i]
    if set(R)!=set(range(1,332)):raise ValueError('Printed R dictionary is incomplete')
    return text,bylabel,R

# V52: supplement sec:cat-component-table (1.2.14); supplement sec:cat-contact-partitions (1.5.2); supplement eq:cat-partition-homogeneous-polynomials (501); supplement eq:cat-partition-auxiliary-minus (503).
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,default=ROOT/'supplement/source/mathematical_catalogue.tex');ap.add_argument('--results',type=Path,default=RESULTS);args=ap.parse_args();start=time.monotonic();D=args.results/'radial_generic';out=args.results/'printed_vs_recomputed';out.mkdir(parents=True,exist_ok=True)
    pairs={p.stem:decode_file(p) for p in (D/'pair').glob('*.json')};contacts={p.stem:json.loads(p.read_text()) for p in (D/'contacts').glob('*.json')}
    if len(pairs)!=54 or len(contacts)!=81:raise FileNotFoundError('Requires 54 newly generated pair records and 81 assembled contacts')
    text,bylabel,R=read_printed(args.source);checks=[];bad=[]
    def check(identity,actual,expected,labels):
        residual=sp.cancel(sp.expand(actual-expected))
        if residual!=0:residual=sp.simplify(residual)
        row={'id':identity,'status':'PASS' if residual==0 else 'FAIL','residual':str(residual),'tolerance':'exact zero','publication_labels':labels,'function':'main'};checks.append(row)
        if residual!=0:bad.append(row)
    maps={};pairmaps={};maplabels={};pairlabels={}
    for label,eq in bylabel.items():
        for part,indices,power,idx in re.findall(r'([hbuv])\^\{(\d{4})\}_\{(-?\d+)\}=R_\{(\d+)\}',eq):maps[part,indices,int(power)]=R[int(idx)];maplabels[part,indices,int(power)]=[label,'eq:cat-rational-'+idx]
        for part,indices,I,power,idx in re.findall(r'([AB])\^\{(\d{4});(\d)\}_\{(-?\d+)\}=R_\{(\d+)\}',eq):pairmaps[part,indices,int(I),int(power)]=R[int(idx)];pairlabels[part,indices,int(I),int(power)]=[label,'eq:cat-rational-'+idx]
        for part,indices,I in re.findall(r'([AB])\^\{(\d{4});(\d)\}=0',eq):pairlabels[part,indices,int(I),'all']=[label]
    reps=sorted({key[1] for key in maps})
    if len(reps)!=13:raise ValueError('Expected 13 printed contact representatives')
    for key in reps:
        rec=contacts[key];parts={'h':load(rec['full_H4']),'b':load(rec['bare_H4']),'u':{},'v':{}}
        for pp in rec['partitions']:
            for dst,src in [('u','auxiliary'),('v','legendre')]:parts[dst]=a.la_add(parts[dst],load(pp[src]))
        for part,values in parts.items():
            for power in range(-4,1):
                check(f'representative.{key}.{part}.{power}',restored(values.get(power,a.Z)),maps.get((part,key,power),sp.S.Zero),maplabels.get((part,key,power),['sec:cat-contacts']))
        print('Compared contact',key,flush=True)
    rows=re.findall(r'^(\d{4}) & (\d{4}) & ([01]) & ([01]) & ([01])\\\\',text,re.M)
    if len(rows)!=81:raise ValueError('Missing printed reconstruction rows')
    for key,seed,flip,transpose,nonzero in rows:
        values=load(contacts[key]['full_H4'])
        for power in range(-4,1):
            expected=maps.get(('h',seed,power),sp.S.Zero) if nonzero=='1' else sp.S.Zero
            if flip=='1':expected=expected.subs(cat.z,-cat.z)
            if transpose=='1':expected=cat.r**(4+power)*expected.subs(cat.r,1/cat.r).xreplace({sp.I:-sp.I})
            check(f'full_contact.{key}.{power}',restored(values.get(power,a.Z)),expected,['sec:cat-component-table','eq:cat-contact-symmetries'])
    for key,rec in sorted(pairs.items()):
        for part,powers in [('A',range(-3,1)),('B',range(-2,1))]:
            for I,values in enumerate(rec[part]):
                if set(values)-set(powers):raise ValueError(('Unexpected source powers',key,part,I))
                for power in powers:
                    expected=pairmaps.get((part,key,I,power),sp.S.Zero);labels=pairlabels.get((part,key,I,power),pairlabels.get((part,key,I,'all'),['eq:cat-pair-polynomials']))
                    check(f'pair.{key}.{part}.{I}.{power}',restored(values.get(power,a.Z)),expected,labels)
        print('Compared pair',key,flush=True)
    plus={}
    for label,eq in bylabel.items():
        m=re.fullmatch(r'eq:cat-auxiliary-plus-(\d{4})-m(\d)',label)
        if not m:continue
        rhs=eq.split('=',1)[1].replace(r'q_{+}',r'\sqrt{1+r^{2}-2 r z}')
        plain=cat.plain_latex(rhs).replace('[','(').replace(']',')')
        plus[m[1],-int(m[2])]=(cat.parse_expr(plain,local_dict=cat.LOCALS,transformations=cat.TRANSFORMS),label)
    if len(plus)!=46:raise ValueError(f'Expected 46 printed crossed-constraint coefficients, found {len(plus)}')
    for key in reps:
        values=load(contacts[key]['partitions'][1]['auxiliary'])
        for power in range(-4,1):
            expected,label=plus.get((key,power),(sp.S.Zero,'eq:cat-partition-auxiliary-minus'))
            check(f'crossed_auxiliary_plus.{key}.{power}',restored(values.get(power,a.Z)),expected,[label])
    # Independent literal homogeneous polynomials printed in Sec. 1.5.2.
    rr=cat.r;ss=sp.sqrt(3)
    for prefix in ['01','23']:
        sign=-1 if prefix=='01' else 1;rad=1 if prefix=='01' else rr
        for aa in range(3):
            for bb in range(3):
                expected={-2:-18,-1:sign*20*ss*sp.I*rad/3,0:sp.Rational(4,3)*rad**2} if aa==bb==0 else ({-2:(-2 if aa==1 else 2),-1:(4*sign*sp.I*rad if aa==1 else -4*sign*sp.I*rad)} if aa==bb else {})
                values=pairs[f'{prefix}{aa}{bb}']['Jalpha']
                for power in range(-2,1):check(f'homogeneous_lapse.{prefix}{aa}{bb}.{power}',restored(values.get(power,a.Z)),expected.get(power,sp.S.Zero),['eq:cat-partition-homogeneous-polynomials'])
    report={'status':'PASS' if not bad else 'FAIL','checks':len(checks),'failed':len(bad),'seconds':time.monotonic()-start,'source_sha256':__import__('hashlib').sha256(args.source.read_bytes()).hexdigest(),'scope':'Independently generated ADM contacts and A/B sources versus unchanged printed V52 coefficient functions; 46 explicit crossed coefficients and homogeneous lapse polynomials included.','tests':checks}
    (out/'checks.json').write_text(json.dumps(report,indent=2)+'\n');print(report['status'],report['checks'],'tests;',report['failed'],'failures');return 0 if not bad else 1
if __name__=='__main__':raise SystemExit(main())
