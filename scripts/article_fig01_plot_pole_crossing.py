#!/usr/bin/env python3
"""Recompute V52 Fig. 1: unreduced finite-k pole crossing, not a divergent ODE.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics, V52; release V52-R1.
Section 6.4, fig:pole-crossing-c12; eq:pole-denominator,
eq:unreduced-pole-action, eq:pole-A-explicit, eq:pole-unreduced-constraint,
eq:pole-constraint-coefficients, eq:pole-cauchy-determinant-general.
Method: DOP853 applied to Euler--Lagrange plus differentiated constraint;
at every evaluation solve a 2x2 system for Vddot and zetadot. Never divide
by D. All initial conditions and tolerances are read from the JSON below.
Inputs: parameters/figure01.json. Outputs: CSV, PDF, PNG and checks.json.
Run: python scripts/article_fig01_plot_pole_crossing.py --output results/fig01
Dependencies: NumPy, SciPy, Matplotlib. No LaTeX, network or article file.
Provenance: NEW implementation of the explicitly specified V52 calculation;
no original plotting program or numerical samples were supplied. No fit to
an image and no change to the article or its original figure are performed.
Scope: one stated trajectory and numerical checks, not a proof of general
regularity. All units and normalizations follow the article (H*=a*=M_Pl=1).

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article 6.4, eq:standard-pole-ratio (139).
  article 6.4, eq:unreduced-pole-action (136).
  article 6.4, eq:pole-A-explicit (137).
  article 6.4, eq:pole-unreduced-constraint (140).
  article 6.4, eq:pole-constraint-coefficients (141).
  article 6.4, eq:pole-cauchy-determinant (145).
  article 6.4, fig:pole-crossing-c12 (1).
Method: DOP853 integration of unreduced Euler equation plus differentiated constraint; determinant/residual/reversibility checks.
Inputs: parameters/figure01.json read at runtime.
Outputs below the selected results root: figure01/pole_crossing_c12.csv, figure01/pole_crossing_c12.pdf, figure01/checks.json.
Provenance: NEW implementation from V52 equations; original figure retained.
Scope limit: One specified regular trajectory. Scientific curves, not byte identity to original PDF.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def coefficients(t:float,c:float):
    # eq:pole-denominator and eq:pole-constraint-coefficients.
    eps=(c-1)/c; b=1+eps*t
    if b<=0:raise ValueError('Time lies outside the expanding standard solution')
    H=1/b; a=b**(1/eps); q=(1/3)/a**2; chi=-6*eps*H**2
    Hd=-eps*H**2;qd=-2*H*q;chid=-2*eps*H*chi
    D=2*c*q+(2-c)*chi;Dd=2*c*qd+(2-c)*chid
    p1=-2*q/(H*(c-2));p2=-2*q/c;p3=q*D/(3*H**2*(c-2))
    p1d=p1*(-2+eps)*H;p2d=-2*H*p2
    p3d=(qd*D+q*Dd-2*q*D*Hd/H)/(3*H**2*(c-2))
    return H,a,q,D,np.array([p1,p2,p3]),np.array([p1d,p2d,p3d])

def rhs(t,y,c):
    # Euler--Lagrange equation from eq:unreduced-pole-action with measure a^3 dt.
    # Second row is d/dt of eq:pole-unreduced-constraint, NOT zeta elimination.
    V,Vd,zeta=y;H,a,q,D,p,pd=coefficients(t,c);A=3/(c*(c-2))
    M=np.array([[2*A,p[0]],[p[0],p[2]]])
    source=np.array([-6*H*A*Vd-(pd[0]+3*H*p[0]-p[1])*zeta,
                     -(pd[0]+p[1])*Vd-pd[1]*V-pd[2]*zeta])
    Vdd,zetad=np.linalg.solve(M,source)
    return np.array([Vd,Vdd,zetad])

# V52: article fig:pole-crossing-c12 (1).
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parameters',type=Path,default=ROOT/'parameters/figure01.json')
    ap.add_argument('--output',type=Path,default=ROOT/'results/current/figure01')
    args=ap.parse_args();p=json.loads(args.parameters.read_text());args.output.mkdir(parents=True,exist_ok=True)
    c=p['c'];t0,t1=p['interval_tau'];V0=p['V_initial'];z0=p['zeta_initial']
    if c!=1.2 or p['q_star']!=1/3:raise ValueError('This V52 figure implementation supports the published c=1.2, q*=1/3 only')
    H,a,q,D,w,wd=coefficients(t0,c);vd0=-(w[1]*V0+w[2]*z0)/w[0];y0=np.array([V0,vd0,z0])
    kw=dict(method='DOP853',rtol=p['rtol'],atol=p['atol'],max_step=p['max_step'],dense_output=True)
    sol=solve_ivp(lambda t,y:rhs(t,y,c),(t0,t1),y0,**kw)
    if not sol.success:raise RuntimeError(sol.message)
    ts=np.linspace(t0,t1,p['samples']);ys=sol.sol(ts);diag=np.array([coefficients(t,c)[:4] for t in ts]);Hs,aa,qs,Ds=diag.T
    constraints=[];detres=[];dets=[]
    for t,y in zip(ts,ys.T):
        H,a,q,D,w,wd=coefficients(t,c);constraints.append(w[0]*y[1]+w[1]*y[0]+w[2]*y[2])
        det=np.linalg.det([[6/(c*(c-2)),w[0]],[w[0],w[2]]]);target=12*q*(c-1)/(c*c*(c-2))
        dets.append(det);detres.append(det-target)
    back=solve_ivp(lambda t,y:rhs(t,y,c),(t1,t0),ys[:,-1],**kw)
    if not back.success:raise RuntimeError(back.message)
    error=float(np.max(np.abs(back.y[:,-1]-y0)))
    Hl=(c-1)*ys[2]-ys[0]
    values=np.column_stack([ts,ys.T,Hl,Hs,aa,qs,Ds,constraints,dets])
    np.savetxt(args.output/'pole_crossing_c12.csv',values,delimiter=',',header='tau,V,V_dot,zeta,H_times_ell,H,a,q,D,constraint_residual,Cauchy_determinant',comments='')
    checks={'constraint_max_abs':float(np.max(np.abs(constraints))),'constraint_tolerance':1e-8,
            'Cauchy_determinant_min_abs':float(np.min(np.abs(dets))),
            'determinant_formula_max_abs_error':float(np.max(np.abs(detres))),
            'roundtrip_max_abs_error':error,'roundtrip_tolerance':2e-7,
            'finite_all_samples':bool(np.isfinite(values).all()),
            'pole_bracketed':bool(Ds[0]*Ds[-1]<0),
            'scope':'One specified initial condition; not a general proof of continuation',
            'parameters':p,'initial_vector':y0.tolist()}
    checks['status']='PASS' if checks['constraint_max_abs']<=1e-8 and error<=2e-7 and checks['finite_all_samples'] and checks['Cauchy_determinant_min_abs']>0 and checks['determinant_formula_max_abs_error']<1e-12 and checks['pole_bracketed'] else 'FAIL'
    (args.output/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(6.4,4.2),layout='constrained')
    ax.plot(ts,ys[0],label=r'$\mathcal{V}$');ax.plot(ts,Hl,label=r'$H\ell$')
    ax.axvline(0,linestyle='--',label=r'$\mathsf{D}=0$');ax.set_xlabel(r'$\tau=H_\star(t-t_\star)$');ax.set_ylabel('Dimensionless perturbation amplitude');ax.legend()
    for ext in ['pdf','png']:fig.savefig(args.output/('pole_crossing_c12.'+ext),dpi=160)
    plt.close(fig);print(json.dumps(checks,indent=2))
    return 0 if checks['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
