#!/usr/bin/env python3
"""V52-R1: article appC3 verify radial controls.
Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (scientific version V52).
Method: exact rational ADM, source contraction, or coefficient validation,
as implemented below; derived from original_scripts/verify_radial_controls_v40.py unchanged in
scientific algebra. Only imports, paths, checkpoint guards and reporting adapted.
Assumptions: radiation background U=0, k=g_R=1, r>0, t=tan(theta/2)>0.
External real TT tensors have norm 2. Internal norms and homogeneous kinetic
signs are retained explicitly. No regulator-free unitary S-matrix is asserted.
Inputs and outputs: see SCRIPT_REFERENCE_MAP.md and REPRODUCIBILITY.md.
Run: python scripts/article_appC3_verify_radial_controls.py (see --help for generator options).
Dependencies: Python and SymPy; local lib modules. No network or LaTeX required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.3, eq:radiation-cubic-source-definition (352).
  article C.7, eq:v34-global-source-polynomials (401).
Method: 18 scalar/tensor source and homogeneous checks, including six sampled crossing controls.
Inputs: New pair records; explicit generic control formulas.
Outputs below the selected results root: radial_generic/independent_controls.csv.
Provenance: GitHub_bundle/original_scripts/verify_radial_controls_v40.py.
Scope limit: Some controls are exact substitutions at selected points, not identities at every angle.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"lib"))
from release_paths import RESULTS
import json,csv
import sympy as s
from supplement_sec1_2_assemble_contacts import load
import radial_coeff_algebra as a
B=Path(__file__).resolve().parents[1];D=RESULTS/'radial_generic';r,t=a.RADIUS,a.HALF_ANGLE;x=s.Symbol('x',positive=True);CHECK=[]
def check(name,ex):
    v=s.simplify(s.expand(ex))
    if v!=0:raise AssertionError((name,v))
    CHECK.append([name,'PASS','0'])
def pol(d):return sum(v.expr()*x**p for p,v in load(d).items())
# V52: article eq:radiation-cubic-source-definition (352); article eq:v34-global-source-polynomials (401).
def main():
    d=json.loads((D/'pair/0200.json').read_text())
    for z in [s.S.Zero,s.Rational(1,3),-s.Rational(1,3)]:
        sub={r:1,t:s.sqrt((1-z)/(1+z))};q=s.sqrt(2*(1-z))
        nu=q/s.sqrt(3)
        val=(pol(d['A'][0])-s.I*nu*pol(d['B'][0])).subs(sub)
        wanted=-2*s.I*nu*(1-z)/3+(-2*z*z+s.Rational(10,3)*z-s.Rational(8,3))/x-s.I*nu*(6*z+14)/x**2+(-6*z*z+6*z-28)/x**3
        check(f'scalar_cross_end_z={z}',val-wanted)
        ep=-(1+z)/2
        wanted=4*s.sqrt(3)*((z-s.Rational(2,3))/x-s.I*q/x**2+(3*z-4)/x**3)*ep
        val=(pol(d['A'][1])-s.I*q*pol(d['B'][1])).subs(sub)
        check(f'tensor_cross_end_z={z}',val-wanted)
    d=json.loads((D/'pair/0100.json').read_text())
    check('global_scalar_coordinate',pol(d['A'][0]))
    check('global_scalar_velocity',pol(d['B'][0])-(s.Rational(4,3)-16*s.I/s.sqrt(3)/x-8/x**2))
    eps=[0,-2,0,0,0]
    for I,ep in enumerate(eps,1):
        check(f'global_tensor_coordinate_{I}',pol(d['A'][I])-(4/s.sqrt(3)/x-4*s.I/x**2-4*s.sqrt(3)/x**3)*ep)
        check(f'global_tensor_velocity_{I}',pol(d['B'][I])-(4*s.I/x+4*s.sqrt(3)/x**2)*ep)
    with (D/'independent_controls.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['id','status','residual']);w.writerows(CHECK)
    print('PASS',len(CHECK),'independent cubic and homogeneous controls')
if __name__=='__main__':main()
