#!/usr/bin/env python3
"""V52 finite-time radiation checks and nonlinear identities, release V52-R1.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics. Scientific version V52.
Purpose: rerun the preserved exact ADM, Legendre, finite-window, scalar
angular, homogeneous-state and nonlinear tests; not a general S matrix.
Method: exact SymPy algebra; no source tables are used as computed answers.
Inputs: documented embedded test polynomials and momenta in the engine.
Outputs: engine_regressions/results.json and checks.json.
Run: python scripts/article_appC_verify_finite_time_and_nonlinear_identities.py
Dependencies: SymPy, lib/article_appC_noncom_engine.py.
Provenance: NEW explicit driver of the preserved noncom_original.py tests.
Scope: finite positive endpoints, radiation U=0, stated special nonlinear
families. Scalar j=0..8 projections are NOT the 74 mixed angular positions.
See SCRIPT_REFERENCE_MAP.md for per-function labels and printed numbers.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.2, eq:v34r-exact-scalar-density (342).
  article C.2, eq:v34r-exact-gravity-density (343).
  article C.2, eq:v34r-extrinsic-density (344).
  article C.2, eq:v34r-nonzero-constraints (345).
  article 6.10, eq:v33-legendre-contact (191).
  article C.3, eq:radiation-noncom-check-momenta (356).
  article C.4, eq:v34-elementary-single (359).
  article C.4, eq:v34-elementary-ordered (360).
  article C.4, eq:v34-single-time-series (361).
  article C.4, eq:v34-ordered-time-series (362).
  article C.4, eq:v34-time-order-identity (363).
  article C.2, eq:v34r-global-reduction (348).
  article C.7, eq:v34r-global-free-flow (400).
  article C.7, eq:v34-global-contraction (402).
  article C.4, eq:v34-scalar-leading-contact (357).
  article C.4, eq:v34-scalar-equal-radii-contact (358).
  article C.3, eq:radiation-contact-decomposition (351).
  article C.5, eq:v34-angular-exchange-series (374).
  article C.5, eq:v34-angular-contact-series (375).
  article C.3, eq:radiation-cubic-source-definition (352).
  article 6.10, eq:v33-channel-counts (194).
  article C.8, eq:v34-finite-time-leakage (411).
  article C.10, eq:v34-crossing-series-H (428).
  article C.10, eq:v34-crossing-series-L (429).
  article C.10, eq:v34r-Bianchi-monotonicity (437).
  article C.10, eq:v34r-simple-wave-equation (431).
  article C.10, eq:v34r-simple-wave-characteristics (432).
Method: 338 exact tests, including direct ADM points, finite-time/scalar angular/homogeneous and nonlinear identities.
Inputs: Embedded exact polynomials and special ADM configurations in lib/article_appC_noncom_engine.py.
Outputs below the selected results root: engine_regressions/checks.json, engine_regressions/results.json.
Provenance: NEW explicit driver of 13 preserved noncom_original.py test functions.
Scope limit: Per-function scope only. Scalar j=0..8 moments are separate from the 74 mixed angular positions. No regulator-free S matrix.
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from release_paths import RESULTS
import article_appC_noncom_engine as engine

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=RESULTS/'engine_regressions')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    names=['test_radiation_expansion','test_general_Legendre_transform',
           'test_complete_example','test_time_integrals','test_signed_Legendre_transform',
           'test_homogeneous_and_state','test_generic_scalar_contact','test_direct_com',
           'test_angular_projection','test_radial_kernel','test_channel_count_and_leakage',
           'test_nonlinear','test_simple_wave']
    engine.CHECKS.clear(); results={};records=[];start=time.monotonic()
    for index,name in enumerate(names):
        t=time.monotonic();before=len(engine.CHECKS)
        print(f'START {index+1}/{len(names)} {name}',flush=True)
        results[name]=getattr(engine,name)()
        for row in engine.CHECKS[before:]:row['function']=name
        records.append({'function':name,'checks':len(engine.CHECKS)-before,'seconds':time.monotonic()-t,'status':'PASS'})
        (args.output/'checkpoint.json').write_text(json.dumps({'status':'RUNNING','completed':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
        (args.output/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        print('DONE',name,records[-1]['checks'],'checks',flush=True)
    out={'status':'PASS','checks':len(engine.CHECKS),'seconds':time.monotonic()-start,'groups':records,'tests':engine.CHECKS,
         'scope':'Preserved exact identities and special radiation/scalar/nonlinear calculations; not a full arbitrary-background multichannel scattering calculation.'}
    (args.output/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    (args.output/'checkpoint.json').write_text(json.dumps({'status':'PASS','completed':records},indent=2)+'\n')
    print('PASS',len(engine.CHECKS),'exact checks');return 0
if __name__=='__main__':raise SystemExit(main())
