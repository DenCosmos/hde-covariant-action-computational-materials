#!/usr/bin/env python3
"""Exact mixed-scalar canonical normalization checks for article V52, V52-R1.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics. Appendix C.1, (326)-(335).
NEW small implementation directly from eq:v34-canonical-matrices through
 eq:v34-group-velocity; no recovered historical matrix calculation is claimed.
Method: generic two-by-two symbolic matrices, independent Euler differentiation,
polynomial reduction modulo the characteristic equation for simple poles.
Assumptions: positive K, real symmetric W, real antisymmetric Omega; the
residue needs a simple positive root, nonzero eigenvector and positive norm.
Inputs: generic symbolic matrices; outputs mixed_normalization/checks.json.
Run: python scripts/article_appC1_verify_mixed_normalization.py
Dependencies: SymPy. No supplied article, LaTeX, table, state or network needed.
This tests normalization identities, not mixed scattering eigenvalues, full
background gravitational corrections, positivity for all models or unitarity.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.1, eq:v34-canonical-matrices (326).
  article C.1, eq:v34-antisymmetric-matrix (327).
  article C.1, eq:v34-canonical-frequency (328).
  article C.1, eq:v34-canonical-action (329).
  article C.1, eq:v34-canonical-evolution (330).
  article C.1, eq:v34-symplectic-norm (331).
  article C.1, eq:v34-frequency-determinant (332).
  article C.1, eq:v34-frequency-roots (333).
  article C.1, eq:v34-spectral-residue (334).
  article C.1, eq:v34-group-velocity (335).
Method: Independent canonical action boundary check, Euler differentiation, symplectic identity, determinant and simple-pole polynomial remainders.
Inputs: Generic real symmetric W and real antisymmetric Omega; positive kinetic matrix; simple positive root, positive nonzero norm.
Outputs below the selected results root: mixed_normalization/checks.json.
Provenance: NEW implementation directly from unambiguous V52 equations; not found historical code.
Scope limit: 11 identities; no full independent-matter mixed scattering spectrum or all-background gravitational correction is computed.
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from release_paths import RESULTS

# V52: article eq:v34-canonical-matrices (326); article eq:v34-antisymmetric-matrix (327); article eq:v34-canonical-frequency (328); article eq:v34-canonical-action (329); article eq:v34-canonical-evolution (330); article eq:v34-symplectic-norm (331); article eq:v34-frequency-determinant (332); article eq:v34-frequency-roots (333); article eq:v34-spectral-residue (334); article eq:v34-group-velocity (335).
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=RESULTS/'mixed_normalization')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True);start=time.monotonic();tests=[]
    def check(name,expression,labels,method='exact algebra'):
        vals=list(expression) if isinstance(expression,s.MatrixBase) else [expression]
        residuals=[s.factor(s.cancel(s.expand(x))) for x in vals]
        if any(x!=0 for x in residuals):raise AssertionError((name,residuals))
        tests.append({'id':name,'status':'PASS','residuals':[str(x) for x in residuals],'method':method,'publication_labels':labels,'tolerance':'exact zero'})
    q=s.Matrix(s.symbols('q1 q2',real=True));v=s.Matrix(s.symbols('v1 v2',real=True));dd=s.Matrix(s.symbols('dd1 dd2',real=True))
    g1,g2,c11,c12,c21,c22,ds11,ds12,ds22,d11,d22,v11,v12,v22,p=s.symbols('g1 g2 c11 c12 c21 c22 ds11 ds12 ds22 d11 d22 v11 v12 v22 p',real=True)
    G=s.diag(g1,g2);C=s.Matrix([[c11,c12],[c21,c22]]);D=s.diag(d11,d22);V=s.Matrix([[v11,v12],[v12,v22]])
    O=(C-C.T)/2;S=(C+C.T)/2-G;Sd=s.Matrix([[ds11,ds12],[ds12,ds22]])
    W=p*p*D-V+G*C+C.T*G-G*G+Sd
    old=((v-G*q).dot(v-G*q)+2*((v-G*q).T*C*q)[0]-(q.T*(p*p*D-V)*q)[0])/2
    new=(v.dot(v)+2*(v.T*O*q)[0]-(q.T*W*q)[0])/2
    boundary=(v.T*S*q)[0]+(q.T*Sd*q)[0]/2
    check('canonical_action_up_to_boundary',old-new-boundary,['eq:v34-canonical-matrices','eq:v34-antisymmetric-matrix','eq:v34-canonical-frequency','eq:v34-canonical-action'])
    o,od,w,a,b,c=s.symbols('o od w a b c',real=True);O=s.Matrix([[0,o],[-o,0]]);Od=s.Matrix([[0,od],[-od,0]]);W=s.Matrix([[a,c],[c,b]])
    lag=(v.dot(v)+2*(v.T*O*q)[0]-(q.T*W*q)[0])/2
    momentum=s.Matrix([s.diff(lag,x) for x in v]);dtmom=momentum.jacobian(q)*v+momentum.jacobian(v)*dd+momentum.diff(o)*od
    euler=dtmom-s.Matrix([s.diff(lag,x) for x in q])
    check('canonical_momentum',momentum-v-O*q,['eq:v34-canonical-evolution'])
    check('Euler_equations',euler-dd-2*O*v-(Od+W)*q,['eq:v34-canonical-evolution'])
    fd=s.Matrix(1,2,s.symbols('fdag1 fdag2'));fvd=s.Matrix(1,2,s.symbols('fveldag1 fveldag2'))
    gg=s.Matrix(s.symbols('gsol1 gsol2'));gv=s.Matrix(s.symbols('gvel1 gvel2'))
    # Dagger evolution follows from Omega^T=-Omega, W^T=W, without freezing.
    fdd=2*fvd*O+fd*Od-fd*W;gdd=-2*O*gv-(Od+W)*gg
    derivative=s.I*(fd*gdd-fdd*gg+2*fvd*O*gg+2*fd*Od*gg+2*fd*O*gv)[0]
    check('symplectic_norm_conservation',derivative,['eq:v34-symplectic-norm'])
    F=W-w*w*s.eye(2)-2*s.I*w*O;det=s.expand(F.det());B=a+b+4*o*o;CC=a*b-c*c
    check('frozen_determinant',det-(w**4-B*w*w+CC),['eq:v34-frequency-determinant'])
    for sign in [-1,1]:
        root=(B+sign*s.sqrt(B*B-4*CC))/2
        check(f'squared_frequency_root_{sign}',s.expand(root*root-B*root+CC),['eq:v34-frequency-roots'])
    e=s.Matrix([-(c-2*s.I*w*o),a-w*w]);N=(e.conjugate().T*(2*w*s.eye(2)+2*s.I*O)*e)[0]
    check('eigenvector_on_characteristic_surface',s.Matrix([s.rem(s.expand(x),det,w) for x in F*e]),['eq:v34-frequency-determinant','eq:v34-spectral-residue'],'polynomial remainder modulo det F=0')
    numerator=F.adjugate()*N+(e*e.conjugate().T)*s.diff(det,w)
    check('simple_pole_residue',s.Matrix(2,2,[s.rem(s.expand(x),det,w) for x in numerator]),['eq:v34-spectral-residue'],'cleared denominators then polynomial remainder; requires simple pole, N nonzero and e nonzero')
    check('frequency_derivative_norm',(e.conjugate().T*F.diff(w)*e)[0]+N,['eq:v34-spectral-residue','eq:v34-group-velocity'])
    m1,m2,m12=s.symbols('m1 m2 m12',real=True);Dp=s.diag(d11,d22)
    Fp=p*p*Dp+s.Matrix([[m1,m12],[m12,m2]])-w*w*s.eye(2)-2*s.I*w*O
    check('momentum_derivative_for_group_velocity',(e.conjugate().T*Fp.diff(p)*e)[0]-2*p*(e.conjugate().T*Dp*e)[0],['eq:v34-group-velocity'])
    out={'status':'PASS','checks':len(tests),'scalar_residuals':sum(len(x['residuals']) for x in tests),'seconds':time.monotonic()-start,'tests':tests,'scope':'Canonical and simple-pole identities only; no independent mixed scattering matrix, chosen vacuum, stability domain or regulator-free unitary limit.'}
    (args.output/'checks.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',len(tests),'identities');return 0
if __name__=='__main__':raise SystemExit(main())
