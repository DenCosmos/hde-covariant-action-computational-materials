#!/usr/bin/env python3
"""V52-R1: article appC5 compute angular coefficients.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/angular_coefficients_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC5_compute_angular_coefficients.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.5, eq:v34-wigner-polynomial (373).
  supplement 1.5, eq:cat-angular-explicit-wigner-coefficients (492).
  article C.5, eq:v34-transfer-variable (370).
  supplement 1.5, eq:cat-angular-source-transfer (490).
  article C.5, eq:v34-helicity-transform (367).
  article C.5, eq:angular-block-complete-order (379).
  article C.5, eq:angular-exchange-finite-reconstruction (377).
  article C.5, eq:angular-contact-finite-reconstruction (378).
  article C.6, eq:boundary-finite-coefficient-reconstruction (392).
  supplement 1.5, eq:cat-angular-exchange-convolution (494).
  supplement 1.5, eq:cat-angular-contact-convolution (495).
  supplement 1.5, eq:cat-boundary-explicit-basis (496).
Method: Exact helicity conversion, Wigner polynomial, q-transfer Laurent coefficients; all scalar and tensor internal exchanges.
Inputs: New 54 pair and 81 contact records; default j=0,2,3,4.
Outputs below the selected results root: angular_v40/j*.json, angular_v40/block_registry_*.csv.
Provenance: GitHub_bundle/original_scripts/angular_coefficients_v40.py.
Scope limit: 74 coefficient positions, not 74 finite energy scattering eigenamplitudes. Finite cutoffs and state parameters remain explicit.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import argparse,itertools,json,time,math,csv
from pathlib import Path
import sympy as sp
import radial_coeff_algebra as a
from supplement_sec1_2_assemble_contacts import parse,load,decode_file,enc
B=Path(__file__).resolve().parents[1];D=RESULTS/'radial_generic';OUT=RESULTS/'angular_v40'
R,Q=a.F.gens
# In the target field, the second generator denotes q, not t_theta.
# V52: article eq:v34-transfer-variable (370); supplement eq:cat-angular-source-transfer (490).
def transfer(v,sgn=1):
    """Substitute t_theta**2=(1-z)/(1+z), z=sgn*(1+r*r-q*q)/(2*r)."""
    z=sgn*(1+R*R-Q*Q)/(2*R);u=(1-z)/(1+z)
    out={}
    def poly(p):
        parity={mon[1]%2 for mon in p}
        if len(parity)>1:raise ValueError('Mixed half-angle parity')
        odd=next(iter(parity),0)
        return sum((co*R**mon[0]*u**((mon[1]-odd)//2) for mon,co in p.items()),a.F.zero),odd
    for ph,rat in v.d.items():
        n,p=poly(rat.numer);d,k=poly(rat.denom)
        if p!=k:raise ValueError('Uncancelled half-angle square root')
        out[ph]=n/d
    return a.Coeff(out)

def target_laurent(v):
    """Exact Laurent expansion in q; reject a denominator with >1 q degree."""
    out={}
    for phase,rat in v.d.items():
        powers={m[1] for m in rat.denom}
        if len(powers)!=1:raise ValueError('Non-Laurent q denominator: '+str(rat.denom))
        shift=next(iter(powers));den=sum((co*R**m[0] for m,co in rat.denom.items()),a.F.zero)
        byq={}
        for mon,co in rat.numer.items():byq[mon[1]-shift]=byq.get(mon[1]-shift,a.F.zero)+co*R**mon[0]/den
        for l,rat1 in byq.items():out[l]=out.get(l,a.Z)+a.Coeff({phase:rat1})
    return out

# V52: article eq:v34-wigner-polynomial (373); supplement eq:cat-angular-explicit-wigner-coefficients (492).
def wigner(j,m,n,sgn):
    if max(abs(m),abs(n))>j:return a.Z
    z=sgn*(1+R*R-Q*Q)/(2*R);v=a.F.zero
    for b in range(2*j+1):
        arg=[j+n-b,b,m-n+b,j-m-b]
        if min(arg)<0:continue
        v+=sp.Rational((-1)**(m-n+b),2**j*math.prod(math.factorial(k) for k in arg))*(1+z)**(j+(n-m)//2-b)*(1-z)**((m-n)//2+b)
    return a.Coeff({(0,0):v})

def cweights(h):return [(0,a.O)] if h==0 else [(1,a.alg(sp.Rational(1,2))),(2,a.I*(1 if h>0 else -1)/2)]
def channels(j):
    pairs=[(0,0),(0,2),(0,-2),(2,2),(2,-2),(-2,-2)]
    return [p for p in pairs if j>=abs(p[0]-p[1]) and not(p[0]==p[1] and j%2)]
def hel_sum(fn,hels):
    out={}
    for choices in itertools.product(*(cweights(h) for h in hels)):
        key=''.join(str(c[0]) for c in choices);factor=a.O
        for c in choices:factor*=c[1]
        for deg,v in fn(key).items():out[deg]=out.get(deg,a.Z)+factor*v
    return {p:v for p,v in out.items() if v}
def f_boundary(i,j,sa,sb):
    if sa!=sb:return a.Z
    z=(1-a.K.gens[1]**2)/(1+a.K.gens[1]**2)
    # Only cross partitions have a soft propagating scalar line.
    if (i,j) in [(0,3),(1,2)]:z=-z
    rr=a.K.gens[0];q2=1+rr*rr-2*rr*z
    out=6*a.I*rr*(rr-1)*(1+z)/q2
    if sa==0:out/=a.SQ3
    elif sa==1:out*=1+z*z
    else:out*= -2*z
    # Simultaneous pi rotation of the two hard legs leaves this scalar invariant.
    return out

def encoded_r(v):
    ex=v.expr().subs(a.HALF_ANGLE,sp.Symbol('q'))
    return str(sp.factor(ex))
def save_coeff(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

# V52: article eq:v34-helicity-transform (367); article eq:angular-block-complete-order (379); article eq:angular-exchange-finite-reconstruction (377); article eq:angular-contact-finite-reconstruction (378); article eq:boundary-finite-coefficient-reconstruction (392); supplement eq:cat-angular-exchange-convolution (494); supplement eq:cat-angular-contact-convolution (495); supplement eq:cat-boundary-explicit-basis (496).
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--resume',action='store_true');ap.add_argument('--j',type=int,nargs='+',default=[0,2,3,4]);ap.add_argument('--start',type=int,default=0);ap.add_argument('--stop',type=int,default=1000000);args=ap.parse_args()
    if any(j<0 for j in args.j):raise ValueError('j must be nonnegative')
    started=time.monotonic();pairs={p.stem:decode_file(p) for p in (D/'pair').glob('*.json')}
    contacts={p.stem:json.loads(p.read_text()) for p in (D/'contacts').glob('*.json')}
    assert len(pairs)==54 and len(contacts)==81, 'Incomplete radial dependencies'
    import hashlib,ast
    def semantic(path):
        tree=ast.parse(path.read_text())
        for node in ast.walk(tree):
            if hasattr(node,'body') and isinstance(node.body,list):
                node.body[:]=[x for x in node.body if not (isinstance(x,ast.Expr) and isinstance(x.value,ast.Constant) and isinstance(x.value.value,str))]
        return ast.dump(tree,include_attributes=False).encode()
    digest=hashlib.sha256(sp.__version__.encode())
    for dependency in [Path(__file__),B/'lib/radial_coeff_algebra.py',B/'scripts/supplement_sec1_2_assemble_contacts.py']:
        digest.update(semantic(dependency))
    for folder in ['pair','contacts']:
        for source in sorted((D/folder).glob('*.json')):
            digest.update((folder+'/'+source.name).encode());digest.update(source.read_bytes())
    input_fingerprint=digest.hexdigest()

    # Cache the finite transfer polynomials before multiplying the Wigner factor.
    ptarget={}
    for p,d in pairs.items():
        if d['homogeneous']:continue
        sgn=1 if p[:2] in ['02','13'] else -1
        ptarget[p]={'A':[{power:transfer(v,sgn) for power,v in row.items()} for row in d['A']], 'B':[{power:transfer(v,sgn) for power,v in row.items()} for row in d['B']]}
    print('Sources mapped exactly to Q(r,q).',flush=True)
    hc={}
    for name,d in contacts.items():
        groups=[load(d['bare_H4'])]
        groups[0]=a.la_add(groups[0],a.la_add(load(d['partitions'][0]['auxiliary']),load(d['partitions'][0]['legendre'])))
        groups += [a.la_add(load(p['auxiliary']),load(p['legendre'])) for p in d['partitions'][1:]]
        hc[name]=[{power:transfer(v,sgn) for power,v in g.items()} for g,sgn in zip(groups,[1,1,-1])]
    print('Three contact groups mapped exactly.',flush=True)
    rows=[]
    for j in args.j:
      ch=channels(j)
      for ia,inc in enumerate(ch):
       for ic,out in enumerate(ch):
        position=ia*len(ch)+ic
        if not args.start<=position<args.stop:continue
        helicities=inc+out;m=inc[0]-inc[1];n=out[0]-out[1];label=f'j{j}_{ic}_{ia}'
        path=OUT/(label+'.json')
        if path.exists():
            if not args.resume:raise FileExistsError('Angular output exists; choose a new output root or --resume after validating radial inputs')
            if json.loads(path.read_text()).get('input_fingerprint') != input_fingerprint:raise ValueError('Incompatible angular checkpoint: '+path.name)
            rows.append([j,ic,ia,*helicities,path.name]);continue
        record={'input_fingerprint':input_fingerprint,'j':j,'row_out':ic,'column_in':ia,'helicities':list(helicities),'wigner_prefactor_squared':math.prod(math.factorial(t) for t in [j+m,j-m,j+n,j-n]),'bose_denominator_squared':(1+int(inc[0]==inc[1]))*(1+int(out[0]==out[1])),'contacts':[],'exchange':[],'boundary':[]}
        for group,sgn in enumerate([1,1,-1]):
            hh=hel_sum(lambda key:hc[key][group],helicities);coeff={}
            for power,v in hh.items():
                for ell,c in target_laurent(a.K.gens[1]*wigner(j,m,n,sgn)*v).items():coeff[f'{power},{ell}']=encoded_r(c)
            record['contacts'].append({'group':group,'sigma':sgn,'H_p_l':coeff})
        for partition,(left,right),sgn in [(1,((0,2),(1,3)),1),(2,((0,3),(1,2)),-1)]:
          for late,early in [(left,right),(right,left)]:
           for internal in range(3):
            c=a.O/a.SQ3 if internal==0 else a.O;nu=c*a.K.gens[1];norm=a.O if internal==0 else (2*a.O if internal==1 else 2*a.K.gens[1]**2);tau= -1 if internal==2 else 1
            def products(key,kind='bulk'):
                akey=f'{late[0]}{late[1]}{key[late[0]]}{key[late[1]]}';bkey=f'{early[0]}{early[1]}{key[early[0]]}{key[early[1]]}'
                da=ptarget[akey];db=ptarget[bkey]
                va=a.la_add(da['A'][internal],a.la_scale(da['B'][internal],-a.I*nu));vb=a.la_add(db['A'][internal],a.la_scale(db['B'][internal],a.I*nu))
                if kind!='bulk':
                    fa=transfer(f_boundary(*late,int(key[late[0]]),int(key[late[1]])),sgn);fb=transfer(f_boundary(*early,int(key[early[0]]),int(key[early[1]])),sgn)
                    if kind=='fI':va={-1:fa}
                    if kind=='If':vb={-1:fb}
                    if kind=='ff':va={-1:fa};vb={-1:fb}
                return {(p,s):tau*v*w/(2*c*norm) for p,v in va.items() for s,w in vb.items() if v*w}
            vals=hel_sum(products,helicities);coeff={}
            for (p,s),v in vals.items():
                for ell,c0 in target_laurent(wigner(j,m,n,sgn)*v).items():coeff[f'{p},{s},{ell}']=encoded_r(c0)
            record['exchange'].append({'partition':partition,'late':late,'early':early,'internal':internal,'sigma':sgn,'C_p_s_l':coeff})
            if internal==0:
                entry={'partition':partition,'late':late,'early':early,'sigma':sgn}
                for kind in ['fI','If','ff']:
                    vals=hel_sum(lambda key:products(key,kind),helicities);coeff={}
                    for (p,s),v in vals.items():
                        for ell,c0 in target_laurent(wigner(j,m,n,sgn)*v).items():coeff[f'{p},{s},{ell}']=encoded_r(c0)
                    entry[kind]=coeff
                record['boundary'].append(entry)
        save_coeff(path,record);rows.append([j,ic,ia,*helicities,path.name]);print(f'[{time.monotonic()-started:.1f}s] {label} complete',flush=True)
    with (OUT/f'block_registry_{args.j[0]}_{args.start}_{args.stop}.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['j','row_out','column_in','h1','h2','h3','h4','file']);w.writerows(rows)
    save_coeff(OUT/f'summary_{args.j[0]}_{args.start}_{args.stop}.json',{'j':args.j,'positions':len(rows),'bulk_and_boundary':True,'seconds':time.monotonic()-started,'success':True})
if __name__=='__main__':main()
