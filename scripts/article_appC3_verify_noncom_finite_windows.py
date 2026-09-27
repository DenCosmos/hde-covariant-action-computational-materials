#!/usr/bin/env python3
"""Evaluate the freshly generated V52 C.3 non-COM kernels on finite intervals.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics (V52); package V52-R1.
References: eq:radiation-noncom-check-momenta (356),
eq:v34-single-time-series, eq:v34-ordered-time-series (Appendix C.4),
eq:cat-ordered-exchange (supplement Sec. 1.4).
Method: independent tensor-product Gauss--Legendre quadrature on the ordered
triangle xi<y<x<xf. T_reduced=-i sum(h_p I_p)-sum(c_ps D_ps).
No high-frequency, zero-transfer or infinite-time limit is taken. Coupling,
external wave normalization and volume factors are omitted here; they are
nonzero common factors and do not change the zero/nonzero count.
Inputs: recomputed noncom/coefficients.json and parameters/noncom_windows.json.
Outputs: finite_windows.csv and finite_window_checks.json in noncom/.
Run: python scripts/article_appC3_verify_noncom_finite_windows.py
Dependencies: NumPy, SymPy and the exact-expression parser in the ADM engine.
Provenance: NEW independent numerical check of the existing coefficient
algorithm. Numerical nonzero witnesses complement exact structural zero
checks; point tests are not represented as proofs of functional identities.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article C.3, eq:radiation-noncom-check-momenta (356).
  article 6.10, eq:v33-radial-kernel (192).
  article C.4, eq:v34-elementary-single (359).
  article C.4, eq:v34-elementary-ordered (360).
Method: Independent 48/80-order nested quadrature over three positive finite intervals.
Inputs: Newly computed noncom/coefficients.json; parameters/noncom_windows.json read at runtime.
Outputs below the selected results root: noncom/finite_windows.csv, noncom/finite_window_checks.json.
Provenance: NEW independent nested Gauss-Legendre evaluation of the restored V52 finite-window kernel.
Scope limit: Nonzero witnesses and order-convergence comparison, not a rigorous quadrature bound or all-window identity proof.
"""
from __future__ import annotations
import sys,argparse,json,csv
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from article_appC_noncom_engine import parse_exact
ROOT=Path(__file__).resolve().parents[1]

def num(s):return complex(parse_exact(str(s)).evalf(18))

def evaluate(channels,xi,xf,n):
    if not 0<xi<xf:raise ValueError('Require 0<xi<xf')
    nodes,weights=np.polynomial.legendre.leggauss(n)
    x=xi+(xf-xi)*(nodes+1)/2;wx=weights*(xf-xi)/2
    y=xi+(x[:,None]-xi)*(nodes[None,:]+1)/2;wy=(x[:,None]-xi)*weights[None,:]/2
    ip={};dp={};out={}
    def I(p,A):
        key=(p,A)
        if key not in ip:ip[key]=np.sum(wx*x**p*np.exp(-1j*A*x))
        return ip[key]
    def D(p,s,A,B):
        key=(p,s,A,B)
        if key not in dp:dp[key]=np.sum(wx[:,None]*wy*x[:,None]**p*y**s*np.exp(-1j*(A*x[:,None]+B*y)))
        return dp[key]
    for key,rec in channels.items():
        omega=sum((1 if i<2 else -1)*2*(1/np.sqrt(3) if a==0 else 1) for i,a in enumerate(rec['external_species']))
        contact=-1j*sum(num(v)*I(int(p),omega) for p,v in rec['contact']['full_H4'].items())
        exchange=0j
        for e in rec['exchange']:
            A=num(e['A']).real;B=num(e['B']).real
            exchange-=sum(num(v)*D(*map(int,ps.split(',')),A,B) for ps,v in e['D_coefficients'].items())
        out[key]=(contact,exchange,contact+exchange)
    return out

# V52: article eq:radiation-noncom-check-momenta (356); article eq:v33-radial-kernel (192); article eq:v34-elementary-single (359); article eq:v34-elementary-ordered (360).
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,default=ROOT/'results/current/noncom/coefficients.json');ap.add_argument('--parameters',type=Path,default=ROOT/'parameters/noncom_windows.json');ap.add_argument('--output',type=Path,default=ROOT/'results/current/noncom');args=ap.parse_args()
    data=json.loads(args.input.read_text());channels=data['channels'];p=json.loads(args.parameters.read_text())
    if len(channels)!=81:raise ValueError('Full 81-component recomputation required')
    records=[];witnesses=set();max_delta=0.
    for xi,xf in p['windows']:
        low=evaluate(channels,xi,xf,p['orders'][0]);high=evaluate(channels,xi,xf,p['orders'][1]);print('Evaluated',xi,xf,flush=True)
        for key,(contact,exchange,total) in high.items():
            delta=abs(total-low[key][2]);threshold=max(p['absolute_zero_threshold'],p['difference_multiplier']*delta);max_delta=max(max_delta,delta)
            nonzero=abs(total)>threshold
            if nonzero:witnesses.add(key)
            records.append({'component':key,'x_i':xi,'x_f':xf,'contact_real':contact.real,'contact_imag':contact.imag,'exchange_real':exchange.real,'exchange_imag':exchange.imag,'total_real':total.real,'total_imag':total.imag,'abs_total':abs(total),'quadrature_order_difference':delta,'nonzero_threshold':threshold,'nonzero_witness':nonzero})
    exact_nonzero=set()
    for key,r in channels.items():
        # Aggregate exact exchange coefficients before testing structural zeros.
        from article_appC_noncom_engine import aggregate_exchange
        if r['contact']['full_H4'] or aggregate_exchange(r['exchange']):exact_nonzero.add(key)
    checks={'parameters':p,'exact_structurally_nonzero_count':len(exact_nonzero),'finite_window_nonzero_witness_count':len(witnesses),'same_support':witnesses==exact_nonzero,'expected_V52_count':57,'maximum_quadrature_order_difference':max_delta,'missing_witnesses':sorted(exact_nonzero-witnesses),'unexpected_numeric_nonzeros':sorted(witnesses-exact_nonzero),'scope':'Finite windows and two quadrature orders; numerical convergence comparison, not a rigorous quadrature error bound or a regulator-free limit.'}
    checks['status']='PASS' if witnesses==exact_nonzero and len(witnesses)==57 and max_delta<p['convergence_tolerance'] else 'FAIL'
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'finite_windows.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (args.output/'finite_window_checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2));return 0 if checks['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
