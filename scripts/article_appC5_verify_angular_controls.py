#!/usr/bin/env python3
"""V52-R1: article appC5 verify angular controls.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/verify_angular_controls_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC5_verify_angular_controls.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.5, eq:v34-wigner-polynomial (373).
  article C.5, eq:v34-angular-contact-series (375).
Method: Wigner d comparison at theta=pi/3 for 74 positions plus three scalar contact moments.
Inputs: 74 new angular records and symbolic Wigner formula.
Outputs below the selected results root: angular_v40/independent_controls.csv.
Provenance: GitHub_bundle/original_scripts/verify_angular_controls_v40.py.
Scope limit: 77 checks. Wigner point tests alone are not an all-angle identity proof; coefficients were separately generated symbolically.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import json,csv,math
import sympy as s
from sympy.physics.wigner import wigner_d_small
import radial_coeff_algebra as a
from article_appC5_compute_angular_coefficients import channels,wigner,target_laurent
B=Path(__file__).resolve().parents[1];D=RESULTS/'angular_v40';CHECK=[]
def check(name,res):
    val=s.simplify(res)
    if val!=0:raise AssertionError((name,val))
    CHECK.append([name,'PASS','0'])
# V52: article eq:v34-wigner-polynomial (373); article eq:v34-angular-contact-series (375).
def main():
    for j in [0,2,3,4]:
        mat=wigner_d_small(j,s.pi/3);ch=channels(j)
        # r=1, q=1 gives z=1/2 independently of the half-angle parameter.
        for aa in ch:
            for bb in ch:
                m=aa[0]-aa[1];n=bb[0]-bb[1]
                p=s.sqrt(math.prod(math.factorial(k) for k in [j+m,j-m,j+n,j-n]))
                val=wigner(j,m,n,1).expr().subs({a.RADIUS:1,a.HALF_ANGLE:1})*p
                check(f'Wigner_matrix_j{j}_{aa}_{bb}',val-mat[j-m,j-n])
    for j,target in [(0,s.Rational(64,27)),(2,-s.Rational(16,135)),(4,s.S.Zero)]:
        rec=json.loads((D/f'j{j}_0_0.json').read_text());rad=s.Rational(3,2);low=rad-1;high=rad+1;total=0
        for group in rec['contacts']:
            for key,value in group['H_p_l'].items():
                p,ell=map(int,key.split(','))
                if p!=0:continue
                co=s.sympify(value).subs(s.Symbol('r'),rad)
                moment=s.log(high/low) if ell==-1 else (high**(ell+1)-low**(ell+1))/(ell+1)
                total+=co*moment
        total*=math.factorial(j)**2/rad
        check(f'scalar_h0_partial_integral_j{j}',total-target*rad**2)
    with (D/'independent_controls.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['id','status','residual']);w.writerows(CHECK)
    print('PASS',len(CHECK),'independent Wigner and contact projection controls')
if __name__=='__main__':main()
