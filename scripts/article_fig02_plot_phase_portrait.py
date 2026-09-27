#!/usr/bin/env python3
"""Recompute V52 Fig. 2 and its numerical inputs without an old plotting cache.

Article: A minimal local covariant action for holographic dark energy:
constraints, perturbations, and nonlinear dynamics, V52; release V52-R1.
Section 7.6, fig:phase-portrait; eq:compact-phase-system,
eq:compact-fixed-points, eq:compact-eigenvalues, eq:early-regular-manifold.
Method: exact vector field and nullclines plus DOP853 trajectory initialized
on the leading early-time unstable-manifold expansion. Arrows are normalized.
Inputs: parameters/figure02.json (runtime configuration). Outputs: field,
nullcline and trajectory CSV files; PDF, PNG and checks.json.
Run: python scripts/article_fig02_plot_phase_portrait.py --output results/fig02
Dependencies: NumPy, SciPy, Matplotlib; no network or TeX.
Provenance: NEW implementation from V52 equations. c=1.2 and w_m=0 are in
the paper; grid, seed and tolerances are explicitly NEW numerical choices,
not recovered original plotting metadata. Reference figure is not replaced.
Scope: background dynamics, not perturbation stability. For this c the
late-time scalar perturbations have the gradient instability stated in V52.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article 7.6, eq:compact-phase-system (242).
  article 7.6, fig:phase-portrait (2).
  article 7.6, eq:compact-fixed-points (243).
  article 7.6, eq:compact-eigenvalues (244).
  article 7.6, eq:early-regular-manifold (245).
Method: DOP853 phase trajectory, fixed points, eigenvalues, nullclines, half-seed sensitivity test.
Inputs: parameters/figure02.json read at runtime; seed and plot mesh are explicit new reconstruction choices.
Outputs below the selected results root: figure02/phase_portrait_c12.pdf, figure02/phase_portrait_trajectory.csv, figure02/phase_portrait_field.csv, figure02/phase_portrait_nullclines.csv, figure02/checks.json.
Provenance: NEW implementation from V52 phase system; original figure retained.
Scope limit: Background phase portrait, not perturbation stability. Extended trajectory; original only displays an early segment.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def vector_field(N,y,c,wm):
    # eq:compact-phase-system; N=ln(a), not cosmic time.
    z,o=y;B=3*wm-1;C=c*(c+1)
    return np.array([z/2*(1+3*wm+2*z-2*C*z*z-B*o),(1-o)*(B*o+2*C*z*z)])

# V52: article fig:phase-portrait (2); article eq:compact-fixed-points (243); article eq:compact-eigenvalues (244); article eq:early-regular-manifold (245).
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--parameters',type=Path,default=ROOT/'parameters/figure02.json');ap.add_argument('--output',type=Path,default=ROOT/'results/current/figure02');args=ap.parse_args()
    p=json.loads(args.parameters.read_text());o=args.output;o.mkdir(parents=True,exist_ok=True);c=p['c'];wm=p['w_m'];C=c*(c+1)
    if wm!=0 or c!=1.2:raise ValueError('This plotting layout is for published c=1.2,w_m=0')
    def trajectory(seed):
        ini=[seed,C*seed**2];sol=solve_ivp(lambda N,y:vector_field(N,y,c,wm),[0,p['N_max']],ini,method='DOP853',rtol=p['rtol'],atol=p['atol'],max_step=p['max_step'],dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        Ns=np.linspace(0,p['N_max'],p['samples']);return Ns,sol.sol(Ns).T
    Ns,path=trajectory(p['seed_z']);_,fine=trajectory(p['seed_z']/2)
    # Compare Omega(z), not Omega at the same N: the seed changes time origin.
    zcompare=np.linspace(0.02,0.8,200)
    def graph(y):
        # Remove numerically repeated late-time z values before interpolation.
        z,ix=np.unique(y[:,0],return_index=True);return PchipInterpolator(z,y[ix,1])
    seed_difference=float(np.max(np.abs(graph(path)(zcompare)-graph(fine)(zcompare))))
    points=np.array([[0,0],[0,1],[1/c,1]]);residuals=np.array([vector_field(0,y,c,wm) for y in points])
    eig_expected=np.array([[.5,-1],[1,1],[-2-1/c,-1-2/c]])
    def jac(z,v):return np.array([[.5*(1+2*z-2*C*z*z+v)+z/2*(2-4*C*z),z/2],[4*C*z*(1-v),2*v-1-2*C*z*z]])
    eig_error=float(max(np.max(np.abs(np.sort(np.linalg.eigvals(jac(*pt)))-np.sort(ex))) for pt,ex in zip(points,eig_expected)))
    xx=np.linspace(*p['z_plot_range'],p['arrow_columns']);yy=np.linspace(0,1,p['arrow_rows']);zz,oo=np.meshgrid(xx,yy);f=vector_field(0,[zz,oo],c,wm);norm=np.hypot(*f);un=np.divide(f[0],norm,out=np.zeros_like(norm),where=norm>0);vn=np.divide(f[1],norm,out=np.zeros_like(norm),where=norm>0)
    np.savetxt(o/'phase_portrait_field.csv',np.column_stack([zz.ravel(),oo.ravel(),f[0].ravel(),f[1].ravel(),un.ravel(),vn.ravel()]),delimiter=',',header='z,Omega_H,z_prime,Omega_H_prime,normalized_z_prime,normalized_Omega_prime',comments='')
    np.savetxt(o/'phase_portrait_trajectory.csv',np.column_stack([Ns,path]),delimiter=',',header='N_from_seed,z,Omega_H',comments='')
    z=np.linspace(*p['z_plot_range'],2001);nc1=2*C*z*z-2*z-1;nc2=2*C*z*z
    np.savetxt(o/'phase_portrait_nullclines.csv',np.column_stack([z,nc1,nc2]),delimiter=',',header='z,Omega_for_z_prime_zero,Omega_for_Omega_prime_zero',comments='')
    checks={'fixed_point_max_residual':float(np.max(np.abs(residuals))),'eigenvalue_max_error':eig_error,'seed_halving_max_Omega_difference':seed_difference,'seed_halving_tolerance':2e-6,'physical_region_numerical':bool(np.min(path[:,0])>=0 and np.min(path[:,1])>=-1e-12 and np.max(path[:,1])<=1+1e-12),'attractor_endpoint_distance':float(np.linalg.norm(path[-1]-points[-1])),'parameters':p,'sampling_origin':'new numerical choices; original sampling metadata absent'}
    checks['status']='PASS' if checks['fixed_point_max_residual']<1e-12 and eig_error<1e-12 and seed_difference<2e-6 and checks['physical_region_numerical'] and checks['attractor_endpoint_distance']<1e-5 else 'FAIL'
    (o/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(6.4,4.8),layout='constrained');ax.quiver(zz,oo,un,vn,angles='xy',scale_units='xy',scale=25,width=.002)
    for nc,label in [(nc2,r'$\Omega_H^{\prime}=0$'),(nc1,r'$z^{\prime}=0$')]:
        mask=(nc>=0)&(nc<=1);ax.plot(z[mask],nc[mask],linestyle='--',label=label)
    ax.axvline(0,linestyle=':',linewidth=.8);ax.axhline(1,linestyle=':',linewidth=.8)
    ax.plot(path[:,0],path[:,1],linewidth=2,label='Regular early-time solution')
    ax.plot(points[:,0],points[:,1],linestyle='none',marker='o')
    for pt,label,offset in zip(points,['$P_m$','$P_v$','$P_H$'],[(7,7),(7,-16),(-30,-18)]):ax.annotate(label,pt,xytext=offset,textcoords='offset points')
    ax.set_xlim(p['z_plot_range']);ax.set_ylim(-.02,1.04);ax.set_xlabel(r'$z=(HL)^{-1}$');ax.set_ylabel(r'$\Omega_H$');ax.set_title(r'$c=1.2,\quad w_m=0$');ax.legend(loc='lower right',fontsize=8)
    for ext in ['pdf','png']:fig.savefig(o/('phase_portrait_c12.'+ext),dpi=160)
    plt.close(fig);print(json.dumps(checks,indent=2));return 0 if checks['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
