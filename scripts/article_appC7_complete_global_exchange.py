#!/usr/bin/env python3
"""V52-R1: article appC7 complete global exchange.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/complete_angular_global_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC7_complete_global_exchange.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.5, eq:v34-angular-moments (371).
  article C.7, eq:v34r-global-free-flow (400).
  article C.7, eq:v34-global-source-polynomials (401).
  article C.7, eq:v34-global-contraction (402).
  article C.7, eq:v34-isolated-zero-distribution (403).
  article C.6, eq:v34r-boundary-table-mixed-dictionary (390).
  article C.6, eq:v34r-boundary-table-square-dictionary (391).
Method: Project AA,AB,BA,BB products for scalar and five homogeneous tensor directions; both time orderings.
Inputs: New radial source and 74 angular dictionaries; symbolic Gaussian covariance.
Outputs below the selected results root: angular_v40/j*.json, angular_v40/block_registry.csv, angular_v40/complete_summary.json.
Provenance: GitHub_bundle/original_scripts/complete_angular_global_v40.py; pure to_z extracted from editorial generator without executing it.
Scope limit: 12 nonhomogeneous plus 12 homogeneous ordered exchange blocks per position; no covariance or volume silently selected.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import argparse,json,csv,math,itertools
from pathlib import Path
import sympy as sp
import radial_coeff_algebra as a
from supplement_sec1_2_assemble_contacts import decode_file,load
from article_appC5_compute_angular_coefficients import channels,cweights,encoded_r
from half_angle_transform import to_z
B=Path(__file__).resolve().parents[1];D=RESULTS/'radial_generic';OUT=RESULTS/'angular_v40';R,Z=a.F.gens

def W(j,m,n):
    val=a.F.zero
    for b in range(2*j+1):
        ar=[j+n-b,b,m-n+b,j-m-b]
        if min(ar)<0:continue
        val+=sp.Rational((-1)**(m-n+b),2**j*math.prod(math.factorial(t) for t in ar))*(1+Z)**(j+(n-m)//2-b)*(1-Z)**((m-n)//2+b)
    return a.Coeff({(0,0):val})
# V52: article eq:v34-angular-moments (371).
def angular_moment(v):
    out=a.Z
    for phase,rat in v.d.items():
        assert all(mon[1]==0 for mon in rat.denom),'Global source is not a z polynomial'
        de=a.F(rat.denom)
        for mon,co in rat.numer.items():
            if mon[1]%2==0:out+=a.Coeff({phase:2*co*R**mon[0]/((mon[1]+1)*de)})
    return out

# V52: article eq:v34r-global-free-flow (400); article eq:v34-global-source-polynomials (401); article eq:v34-global-contraction (402); article eq:v34-isolated-zero-distribution (403); article eq:v34r-boundary-table-mixed-dictionary (390); article eq:v34r-boundary-table-square-dictionary (391).
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--j',type=int,nargs='+',default=[0,2,3,4],help='Nonnegative angular momenta; default reproduces the 74 published positions.')
    args=parser.parse_args()
    if any(j<0 for j in args.j):parser.error('j must be nonnegative')
    js=list(dict.fromkeys(args.j))
    pairs={p.stem:decode_file(p) for p in (D/'pair').glob('*.json') if p.stem[:2] in ['01','23']}
    # Products of the global Cartesian sources, before helicity projection.
    cache={}
    for keytuple in itertools.product(range(3),repeat=4):
        key=''.join(map(str,keytuple));cache[key]={}
        for late,early in [('01','23'),('23','01')]:
            da=pairs[late+key[int(late[0])]+key[int(late[1])]];db=pairs[early+key[int(early[0])]+key[int(early[1])]]
            for I in range(6):
              for form in ['AA','AB','BA','BB']:
                vals={}
                norm=1 if I==0 else da['norms'][I]
                for p,v in da[form[0]][I].items():
                    for s,w in db[form[1]][I].items():
                        if not v*w:continue
                        even,odd=to_z(v*w/norm)
                        assert not odd,(key,late,I,form,'global half-angle parity')
                        vals[p,s]=even
                cache[key][late,I,form]=vals
    rows=[];checks=0
    for j in js:
        ch=channels(j)
        for ia,inc in enumerate(ch):
          for ic,out in enumerate(ch):
            path=OUT/f'j{j}_{ic}_{ia}.json'
            if not path.exists():raise FileNotFoundError(path)
            obj=json.loads(path.read_text());hs=inc+out;m=inc[0]-inc[1];n=out[0]-out[1];wig=W(j,m,n)
            if 'global_exchange' not in obj:
                glob=[]
                for late,early in [('01','23'),('23','01')]:
                  for I in range(6):
                    rec={'late':list(map(int,late)),'early':list(map(int,early)),'internal':I,'kinetic_sign':-1 if I==0 else 1,'forms':{}}
                    for form in ['AA','AB','BA','BB']:
                        acc={}
                        for choices in itertools.product(*(cweights(h) for h in hs)):
                            key=''.join(str(t[0]) for t in choices);weight=a.O
                            for t in choices:weight*=t[1]
                            for deg,v in cache[key][late,I,form].items():acc[deg]=acc.get(deg,a.Z)+weight*v
                        rec['forms'][form]={f'{p},{s}':encoded_r(angular_moment(wig*v)) for (p,s),v in acc.items() if angular_moment(wig*v)}
                    glob.append(rec)
                obj['global_exchange']=glob
                obj['global_dictionary']='Each form is integral W(z)*A_a,p A_b,s etc / positive tensor norm. Multiply by p_jmn and the unchanged Gaussian X,Y covariance form; keep V normalization separate.'
            # Eliminate fixed endpoint-degree indices. The dictionaries in the
            # manuscript supply x_i^-1 and x_f^-1, exactly once.
            for rec in obj['boundary']:
                if 'fI_p_l' in rec:continue
                rec['fI_p_l']={','.join(k.split(',')[1:]):v for k,v in rec.pop('fI').items()}
                rec['If_p_l']={','.join([k.split(',')[0],k.split(',')[2]]):v for k,v in rec.pop('If').items()}
                rec['ff_l']={k.split(',')[2]:v for k,v in rec.pop('ff').items()}
            assert len(obj['exchange'])==12 and len(obj['global_exchange'])==12 and len(obj['boundary'])==4
            checks+=3
            path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
            rows.append([j,ic,ia,*hs,m,n,path.name])
            print(path.stem,'global and boundary dictionaries complete',flush=True)
    assert len(rows)==sum(len(channels(j))**2 for j in js)
    with (OUT/'block_registry.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['j','out_row','in_column','h0','h1','h2','h3','m','n','file']);w.writerows(rows)
    (OUT/'complete_summary.json').write_text(json.dumps({'success':True,'positions':len(rows),'bulk_positions':len(rows),'boundary_positions':len(rows),'ordered_exchanges_per_position':24,'structural_checks':checks},indent=2))
    print('COMPLETE:',len(rows),'bulk and boundary positions')
if __name__=='__main__':main()
