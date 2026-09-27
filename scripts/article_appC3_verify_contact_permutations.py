"""V52-R1: article appC3 verify contact permutations.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/check_catalogue_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC3_verify_contact_permutations.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.3, eq:radiation-contact-permutations (355).
  supplement 1.2, eq:cat-contact-symmetries (3).
  supplement 1.2.14, sec:cat-component-table (1.2.14).
Method: Within-pair swaps and pair transposition with r-dependent scaling, angle sign and complex conjugation.
Inputs: 81 newly assembled contact dictionaries.
Outputs below the selected results root: radial_generic/permutation_checks.json, radial_generic/reconstruction_map.csv.
Provenance: GitHub_bundle/original_scripts/check_catalogue_v40.py.
Scope limit: 243 full-contact identities plus 729 component identities = 972. Counts are distinct and not substituted for each other.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import itertools,json,csv
import radial_coeff_algebra as a
from supplement_sec1_2_assemble_contacts import parse,load
B=Path(__file__).resolve().parents[1];D=RESULTS/'radial_generic'
def reciprocal(v,index):
    """Exact variable inversion of a rational polynomial; no approximation."""
    def pp(p):
        m=max((mon[index] for mon in p),default=0)
        return a.F.ring.from_dict({tuple(m-e if j==index else e for j,e in enumerate(mon)):co for mon,co in p.items()}),m
    out={}
    for ph,rat in v.d.items():
        n,mn=pp(rat.numer);d,md=pp(rat.denom)
        out[ph]=a.F(n)/a.F(d)*a.F.gens[index]**(md-mn)
    return a.Coeff(out)
def conj(v):return a.Coeff({p:(-1)**p[0]*v for p,v in v.d.items()})
# V52: article eq:radiation-contact-permutations (355); supplement eq:cat-contact-symmetries (3); supplement sec:cat-component-table (1.2.14).
def main():
    values={p.stem:load(json.loads(p.read_text())['full_H4']) for p in (D/'contacts').glob('*.json')}
    assert len(values)==81, "Expected all 81 independently assembled contacts"
    checks=[];rows=[]
    for kind in ['bare_H4','auxiliary','legendre']:
        pieces={}
        for p in (D/'contacts').glob('*.json'):
            rec=json.loads(p.read_text())
            if kind=='bare_H4':v=load(rec[kind])
            else:
                v={}
                for part in rec['partitions']:v=a.la_add(v,load(part[kind]))
            pieces[p.stem]=v
        assert len(pieces)==81
        for name,h in pieces.items():
            aa,bb,cc,dd=name
            for perm in [bb+aa+cc+dd,aa+bb+dd+cc]:
                checks.append({'kind':kind,'component':name,'permutation':perm,'transformation':'t_inverse','status':'PASS'})
                for power in set(h)|set(pieces[perm]):
                    assert not(h.get(power,a.Z)-reciprocal(pieces[perm].get(power,a.Z),1)),(kind,name,perm,power)
            perm=cc+dd+aa+bb
            checks.append({'kind':kind,'component':name,'permutation':perm,'transformation':'r_inverse_conjugate_scale','status':'PASS'})
            for power in set(h)|set(pieces[perm]):
                assert not(h.get(power,a.Z)-a.K.gens[0]**(4+power)*conj(reciprocal(pieces[perm].get(power,a.Z),0))),(kind,name,perm,power)
    for name,h in values.items():
        A,Bb,C,Dd=map(int,name)
        for perm,flip in [(f'{Bb}{A}{C}{Dd}',True),(f'{A}{Bb}{Dd}{C}',True)]:
            checks.append({'kind':'full_H4','component':name,'permutation':perm,'transformation':'t_inverse','status':'PASS'})
            for power in set(h)|set(values[perm]):
                assert not(h.get(power,a.Z)-reciprocal(values[perm].get(power,a.Z),1)),(name,perm,power)
        perm=f'{C}{Dd}{A}{Bb}'
        checks.append({'kind':'full_H4','component':name,'permutation':perm,'transformation':'r_inverse_conjugate_scale','status':'PASS'})
        for power in set(h)|set(values[perm]):
            assert not(h.get(power,a.Z)-a.K.gens[0]**(4+power)*conj(reciprocal(values[perm].get(power,a.Z),0))),(name,perm,power,'transpose')
        key=[A,Bb,C,Dd];flip=False
        if key[0]>key[1]:key[0],key[1]=key[1],key[0];flip=not flip
        if key[2]>key[3]:key[2],key[3]=key[3],key[2];flip=not flip
        transpose=key[:2]>key[2:]
        if transpose:key=key[2:]+key[:2]
        rows.append((name,''.join(map(str,key)),int(flip),int(transpose),int(bool(h))))
    with (D/'reconstruction.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['component','seed','t_to_inverse_t','r_to_inverse_r_conjugate_scale','nonzero']);w.writerows(sorted(rows))
    assert len(checks)==972
    (D/'permutation_checks.json').write_text(json.dumps(checks,indent=2))
    (D/'symmetry_check.json').write_text(json.dumps({'components':81,'permutation_identities':len(checks),'full_contact_identities':sum(x['kind']=='full_H4' for x in checks),'nonzero_seeds':len({row[1] for row in rows if row[-1]}),'success':True},indent=2))
    print('All 972 exact permutation identities (full and three separate parts) passed; seeds=',len({row[1] for row in rows if row[-1]}))
if __name__=='__main__':main()
