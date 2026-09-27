#!/usr/bin/env python3
"""V52-R1: supplement sec1 3 generate radial sources.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/generate_radial_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/supplement_sec1_3_generate_radial_sources.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.5, eq:v34-linear-polarizations (364).
  article C.3, eq:radiation-vertex-coefficient-basis (353).
  article C.3, eq:radiation-cubic-source-definition (352).
  article C.2, eq:v34r-nonzero-constraints (345).
  article C.2, eq:v34r-global-reduction (348).
  supplement 1.3, eq:cat-pair-polynomials (56).
  article C.3, eq:radiation-catalogue-external-momenta (350).
  supplement 1.1, eq:cat-momenta (1).
Method: Exact rational field ADM expansion; 81 bare components, 54 unordered external-pair species combinations.
Inputs: Embedded symbolic r>0,t=tan(theta/2)>0; --task/--start/--stop/--resume control generation.
Outputs below the selected results root: radial_generic/bare/*.json, radial_generic/pair/*.json.
Provenance: GitHub_bundle/original_scripts/generate_radial_v40.py; preserved radiation_adm_v40.py.
Scope limit: TT norm 2 externally; internal cross norm 2q^2; homogeneous scalar negative kinetic sign preserved. Odd-reflection bare zeros are exact symmetry exclusions.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import argparse,itertools,json,os,time,traceback
from pathlib import Path
import sympy as sp
import radiation_adm_v40 as rad

r,t=rad.RADIUS,rad.HALF_ANGLE
Z,O=rad.Z,rad.O
SN=rad.alg(2*t/(1+t*t)); CO=rad.alg((1-t*t)/(1+t*t)); RR=rad.alg(r)
P=((Z,Z,O),(Z,Z,-O),(-RR*SN,Z,-RR*CO),(RR*SN,Z,RR*CO))
SIG=(1,1,-1,-1)
MAG=(O,O,RR,RR)
OUT=RESULTS/'radial_generic'


def atomic(path:Path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    os.replace(tmp,path)


def encode(d):
    return {str(k):str(rad.K.to_sympy(v)) for k,v in sorted(d.items()) if v}


# V52: article eq:v34-linear-polarizations (364); article eq:radiation-vertex-coefficient-basis (353).
def real_basis(p,mag=None):
    """Returns (Cartesian tensor, squared norm). No hidden sign reversal."""
    p2=rad.norm2(p)
    if not p2:
        return [(((O,Z,Z),(Z,-O,Z),(Z,Z,Z)),rad.alg(2)),
                (((O,Z,Z),(Z,O,Z),(Z,Z,-2*O)),rad.alg(6)),
                (((Z,O,Z),(O,Z,Z),(Z,Z,Z)),rad.alg(2)),
                (((Z,Z,O),(Z,Z,Z),(O,Z,Z)),rad.alg(2)),
                (((Z,Z,Z),(Z,Z,O),(Z,O,Z)),rad.alg(2))]
    e=(Z,O,Z); b=(-p[2],Z,p[0])
    if mag is not None:
        b=tuple(x/mag for x in b)
        E1=tuple(tuple(e[i]*e[j]-b[i]*b[j] for j in range(3)) for i in range(3))
        E2=tuple(tuple(e[i]*b[j]+b[i]*e[j] for j in range(3)) for i in range(3))
        return [(E1,rad.alg(2)),(E2,rad.alg(2))]
    E1=tuple(tuple(e[i]*e[j]-b[i]*b[j]/p2 for j in range(3)) for i in range(3))
    E2=tuple(tuple(e[i]*b[j]+b[i]*e[j] for j in range(3)) for i in range(3))
    return [(E1,rad.alg(2)),(E2,2*p2)]


def leg(i,species):
    om=SIG[i]*MAG[i]/(rad.SQ3 if species==0 else O)
    return rad.Leg(P[i],om,'w' if species==0 else 'v',None if species==0 else real_basis(P[i],MAG[i])[species-1][0])


def vertex(legs):
    return rad.Lvertex(tuple(legs))


# V52: article eq:radiation-cubic-source-definition (352); article eq:v34r-nonzero-constraints (345); article eq:v34r-global-reduction (348); supplement eq:cat-pair-polynomials (56).
def pair_data(i,j,a,b):
    ls=(leg(i,a),leg(j,b)); pa=tuple(P[i][k]+P[j][k] for k in range(3)); pp=tuple(-x for x in pa); p2=rad.norm2(pa)
    ja=vertex(ls+(rad.Leg(pp,0,'alpha'),))
    ji=[vertex(ls+(rad.Leg(pp,0,'beta'+str(k)),)) for k in range(3)]
    jb={}
    if p2:
        for k in range(3):jb=rad.la_add(jb,rad.la_scale(ji[k],-rad.I*pp[k]/p2))
        jt=[rad.la_add(ji[k],rad.la_scale(jb,-rad.I*pp[k])) for k in range(3)]
    else:
        if any(ji):raise AssertionError(f'Nonzero homogeneous shift source: {i,j,a,b}')
        jt=[{},{},{}]
    # Basis is chosen using pa, not the opposite probe momentum. Both endpoints
    # are transported into this common Cartesian basis during contractions.
    tensors=real_basis(pa)
    A=[vertex(ls+(rad.Leg(pp,0,'w',None,1,0),))]
    B=[vertex(ls+(rad.Leg(pp,0,'w',None,0,1),))]
    norms=[O if p2 else -O]
    for e,norm in tensors:
        A.append(vertex(ls+(rad.Leg(pp,0,'v',e,1,0),)))
        B.append(vertex(ls+(rad.Leg(pp,0,'v',e,0,1),)))
        norms.append(norm)
    return {'indices':[i,j,a,b],'p':[str(rad.K.to_sympy(x)) for x in pa],
            'p2':str(rad.K.to_sympy(p2)),'homogeneous':not bool(p2),
            'Jalpha':encode(ja),'Ji':[encode(d) for d in ji],'JB':encode(jb),'JT':[encode(d) for d in jt],
            'A':[encode(d) for d in A],'B':[encode(d) for d in B],
            'norms':[str(rad.K.to_sympy(x)) for x in norms]}


# V52: article eq:radiation-catalogue-external-momenta (350); supplement eq:cat-momenta (1).
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--task',choices=['bare','pair'],required=True)
    ap.add_argument('--resume',action='store_true',help='Reuse only checkpoints carrying this run fingerprint.');ap.add_argument('--start',type=int,default=0);ap.add_argument('--stop',type=int)
    args=ap.parse_args();started=time.monotonic()
    import hashlib,ast
    def semantic_hash(path):
        tree=ast.parse(path.read_text())
        for node in ast.walk(tree):
            if hasattr(node,'body') and isinstance(node.body,list):
                node.body[:]=[x for x in node.body if not (isinstance(x,ast.Expr) and isinstance(x.value,ast.Constant) and isinstance(x.value.value,str))]
        return hashlib.sha256(ast.dump(tree,include_attributes=False).encode()).hexdigest()
    fingerprint=hashlib.sha256((semantic_hash(Path(__file__))+semantic_hash(Path(rad.__file__))).encode()).hexdigest()
    atomic(OUT/'input_manifest.json',{'input_fingerprint':fingerprint,'configuration':'r>0,t>0 symbolic; external TT norm 2; internal cross norm 2*p^2; U=0; g_R=k=1','task':args.task})
    tasks=list(itertools.product(range(3),repeat=4)) if args.task=='bare' else [(i,j,a,b) for i,j in itertools.combinations(range(4),2) for a,b in itertools.product(range(3),repeat=2)]
    stop=args.stop if args.stop is not None else len(tasks)
    for idx,key in enumerate(tasks):
        if not args.start<=idx<stop:continue
        name=''.join(map(str,key));path=OUT/args.task/(name+'.json')
        if path.exists():
            if not args.resume:raise FileExistsError(f'{path.name} already exists; use --resume or a new output root')
            previous=json.loads(path.read_text())
            if previous.get('input_fingerprint')!=fingerprint:raise ValueError(f'Incompatible checkpoint: {path.name}')
            continue
        print(f'[{time.monotonic()-started:.1f}s] START {args.task} {idx+1}/{len(tasks)} {name}',flush=True)
        try:
            if args.task=='bare':
                # Reflection y -> -y kills an odd number of E2 factors exactly.
                data={'species':list(key),'parity_zero':sum(k==2 for k in key)%2==1}
                data['bare_L4']={} if data['parity_zero'] else encode(vertex(tuple(leg(i,a) for i,a in enumerate(key))))
            else:data=pair_data(*key)
            data['input_fingerprint']=fingerprint
            atomic(path,data)
            print(f'[{time.monotonic()-started:.1f}s] DONE {name}',flush=True)
        except Exception:
            print(traceback.format_exc(),flush=True);return 1
    return 0

if __name__=='__main__':raise SystemExit(main())
