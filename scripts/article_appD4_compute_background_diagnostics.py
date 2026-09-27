#!/usr/bin/env python3
"""Numerical diagnostics for article V52 Appendix D.4; computational release V52-R1.
Preserved DOP853 algorithms; see the explicit references and parameter record below.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article D.1, eq:scalar-finite-k-complete-coefficients (442).
  article D.2, eq:power-canonical-mass-rational (450).
  article 7.4, eq:power-law-full-background (231).
  article 7.4, eq:power-law-full-eigenvalues (232).
  article 7.4, tab:v31-numerical-dynamics (3).
  article D.4, eq:numerical-crossing-point-v40 (460).
  article D.4, eq:numerical-crossing-initial-background-v40 (461).
  article D.4, tab:crossing-reproduction-v40 (5).
  article D.2, eq:power-canonical-mass-limits (451).
  article D.4, eq:numerical-late-initial-data-v40 (462).
  article D.4, eq:numerical-curvature-envelope-v40 (463).
  article D.4, tab:late-rate-convergence-v40 (6).
Method: DOP853 integration at stored tolerances; auxiliary linear algebra and asymptotic slopes.
Inputs: Embedded parameter grid; see parameters/embedded_parameters.md; runtime output argument only.
Outputs below the selected results root: n02_background_diagnostics/numerical_results.json, n02_background_diagnostics/saddle_rates.csv, n02_background_diagnostics/crossing_convergence.csv, n02_background_diagnostics/late_rates.csv.
Provenance: GitHub_bundle/scripts/n02_background_diagnostics.py.
Scope limit: Numerical examples and tolerances, not proofs of global dynamics.
"""
from __future__ import annotations
import csv
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/current/n02_background_diagnostics'
CHECKS: list[dict] = []
EXPRESSIONS: dict[str, str] = {}
R = sp.Rational

def check(name: str, residual: sp.Expr) -> None:
    value = sp.factor(sp.cancel(residual))
    if value != 0:
        raise AssertionError(f'{name}: {value}')
    CHECKS.append({'id': name, 'status': 'PASS', 'residual': '0'})

def save_csv(name: str, rows: list[dict]) -> None:
    with (OUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def integrate(fun, interval, initial, rtol=1e-12, atol=1e-14, **kwargs):
    sol=solve_ivp(fun, interval, np.asarray(initial,dtype=float), method='DOP853',
                  rtol=rtol, atol=atol, **kwargs)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol

# V52: article eq:scalar-finite-k-complete-coefficients (442); article eq:power-canonical-mass-rational (450).
def symbolic_reduction():
    M,H,L,q,chi,G1,G2 = sp.symbols('M H L q chi G1 G2', nonzero=True)
    alpha,B,r,pi,ell,v=sp.symbols('alpha B r pi ell v')
    lag=M**2*(-3*H**2*alpha**2-2*H*alpha*B)
    lag+=r*(-v+H*ell-alpha-L*B/3+L*q*pi/3)
    lag+=chi*ell*B/3-chi*q*pi**2/6-chi*ell*q*pi/3-H*G1*ell**2/2
    aux=[alpha,B,r,pi]
    mat=sp.hessian(lag,aux)
    source=sp.Matrix([sp.diff(lag,u).subs(dict.fromkeys(aux,0)) for u in aux])
    solution=mat.inv()*(-source)
    subs=dict(zip(aux,solution))
    red=sp.factor(lag.subs(subs))
    for i,u in enumerate(aux): check(f'auxiliary_{u}',sp.diff(lag,u).subs(subs))
    K=sp.factor(sp.diff(red,v,2)); C=sp.factor(sp.diff(red,v,ell)); V=sp.factor(sp.diff(red,ell,2))
    D=2*H*L*M**2*q+(2-H*L)*chi
    check('exact_K_ell',K+6*H*M**2*chi/(L*D))
    check('exact_C_ell',C+chi*(-6*H**2*M**2+2*H*L*M**2*q-H*L*chi+chi)/(L*D))
    check('auxiliary_determinant',mat.det()-R(2,9)*H*L*M**2*q*D)
    def dt(f):
        return sp.diff(f,H)*chi/(6*M**2)+sp.diff(f,L)*(H*L-1)-2*H*q*sp.diff(f,q)\
               +(G1*(H*L-1)-4*H*chi)*sp.diff(f,chi)+G2*(H*L-1)*sp.diff(f,G1)
    W=sp.factor(dt(C)+3*H*C-V)
    K1=-3*G1*(H*L-1)/(L**2*q)
    W0=G1*(3*H**2*L-3*H+L*q)/(L**2*q)
    check('Kdot_at_crossing',dt(K).subs(chi,0)-K1)
    check('W_at_crossing',W.subs(chi,0)-W0)
    check('main_Kdot_at_crossing',dt(K/H**2).subs(chi,0)-K1/H**2)
    nu=sp.symbols('nu')
    check('double_indicial_root',nu*(nu-1)+nu-nu**2)
    # Independent integration-by-parts reduction to Q=delta Phi.
    Q,Qdot=sp.symbols('Q Qdot')
    actionQ=sp.factor(lag.subs({r:chi*ell/L-Q/L**4,alpha:Q/(6*H*M**2*L**3),pi:-Q/(L**3*chi)}))
    cv=sp.diff(actionQ,v)
    actionQ=sp.factor(actionQ-cv*v-ell*(Qdot/L**4+dt(L**-4)*Q+3*H*Q/L**4)
                     +ell**2*(dt(chi/L)+3*H*chi/L)/2)
    ellsol=sp.solve(sp.diff(actionQ,ell),ell)[0]
    Delta=L*G1-chi
    check('Q_to_ell',ellsol+(Qdot+(4/L+chi/(6*M**2*H))*Q)/(L**2*Delta))
    actionQ=sp.factor(actionQ.subs(ell,ellsol))
    AQ=sp.diff(actionQ,Qdot,2)
    CQ=sp.diff(actionQ,Qdot,Q)
    VQ=sp.diff(actionQ,Q,2)
    cs2=-Delta/(3*chi)
    vratio=16/L**2+(Delta+4*chi)/(3*M**2*H*L)+chi**2/(36*M**4*H**2)-Delta/(6*M**2)
    check('Q_kinetic',AQ-1/(L**6*Delta))
    check('Q_mixed',CQ/AQ-(4/L+chi/(6*M**2*H)))
    check('Q_potential',VQ/AQ-(vratio-cs2*q))
    # Compact branch: spatial metric is positive definite, use an orthonormal frame.
    k1,k2,k3,k12,k13,k23=sp.symbols('k1 k2 k3 k12 k13 k23', real=True)
    km=sp.Matrix([[k1,k12,k13],[k12,k2,k23],[k13,k23,k3]])
    shear=km-sp.eye(3)/L
    check('compact_shear_identity',sp.trace(shear*shear)-(sp.trace(km*km)-2*sp.trace(km)/L+3/L**2))
    check('compact_trace_condition', (sp.trace(shear*shear)-(sp.trace(km*km)-3/L**2)).subs(k3,3/L-k1-k2))
    amp=sp.symbols('amplitude',real=True)
    s1=sp.Matrix([[k1,k12,k13],[k12,k2,k23],[k13,k23,-k1-k2]])
    check('second_order_obstruction',sp.expand(sp.trace((amp*s1)*(amp*s1))).coeff(amp,2)-sp.trace(s1*s1))
    # Reconstructed power law: arbitrary off-curve background variations allowed.
    n=sp.symbols('n',real=True)
    y,h=sp.symbols('y h',positive=True)
    g=3*((4-n)*y**(-n)+n*y**(-1-n/2))
    flow=sp.Matrix([h*y-1,(g-12*h*h)/6])
    jac=flow.jacobian([y,h]).subs({y:1,h:1})
    expected=sp.Matrix([[1,1],[n*(n-10)/4,-4]])
    for i in range(2):
        for j in range(2): check(f'jacobian_{i}{j}',jac[i,j]-expected[i,j])
    rr=sp.symbols('r')
    check('characteristic_polynomial',jac.charpoly(rr).as_expr()-(rr+(n-2)/2)*(rr-(n-8)/2))
    parallel=sp.Matrix([1,-n/2]); transverse=sp.Matrix([1,(n-10)/2])
    for i in range(2):
        check(f'tangent_eigenvector_{i}',(jac*parallel+(n-2)*parallel/2)[i])
        check(f'transverse_eigenvector_{i}',(jac*transverse-(n-8)*transverse/2)[i])
    mismatch=3*(h*h-y**(-n))
    dotm=sp.diff(mismatch,y)*flow[0]+sp.diff(mismatch,h)*flow[1]
    check('nonlinear_mismatch_evolution',dotm-(-4*h+n*y**(-n/2)/(y*(h+y**(-n/2))))*mismatch)
    # Radiation example is decelerating and not a dark-energy attractor.
    t,CL=sp.symbols('t CL',positive=True)
    Hrad=1/(2*t); Lrad=CL*sp.sqrt(t)-2*t
    check('radiation_L',sp.diff(Lrad,t)-(Hrad*Lrad-1))
    check('radiation_acceleration',sp.diff(Hrad,t)+2*Hrad**2)
    check('radiation_w',(-(-12*M**2*Hrad**2)/3-3*M**2*Hrad**2)/(3*M**2*Hrad**2)-R(1,3))
    check('late_curvature_exponent',-2-(n+6)/8+(n-2)/2-3*(n-10)/8)
    for key,expr in [('K_ell',K),('C_ell',C),('V_ell',V),('W_ell',W),('Kdot_crossing',K1),('W_crossing',W0),('jacobian',jac)]:
        EXPRESSIONS[key]=str(expr)
    return sp.lambdify((H,L,q,chi,G1,G2),[K.subs(M,1),W.subs(M,1)],'numpy',cse=True)

# V52: article eq:power-law-full-background (231); article eq:power-law-full-eigenvalues (232); article tab:v31-numerical-dynamics (3).
def saddle_numerics():
    print('Checking saddle eigenvalue convergence',flush=True)
    rows=[]
    for n in [10,12,14]:
        for direction,vec,rate in [('parallel',np.array([1.,-n/2]),-(n-2)/2),
                                    ('transverse',np.array([1.,(n-10)/2]),(n-8)/2)]:
            for amplitude in [1e-4,1e-5,1e-6]:
                # Integrate departures directly to avoid subtracting H=1 at each step.
                def rhs(tau,u):
                    dy,dh=u
                    ln=np.log1p(dy)
                    dg=3*((4-n)*np.expm1(-n*ln)+n*np.expm1((-1-n/2)*ln))
                    return [dy+dh+dy*dh, dg/6-4*dh-2*dh*dh]
                vals=[]
                for tol in [1e-9,1e-12]:
                    plus=integrate(rhs,(0,0.5),amplitude*vec,rtol=tol,atol=tol*amplitude*0.01).y[:,-1]
                    minus=integrate(rhs,(0,0.5),-amplitude*vec,rtol=tol,atol=tol*amplitude*0.01).y[:,-1]
                    projection=np.dot((plus-minus)/(2*amplitude),vec)/np.dot(vec,vec)
                    vals.append(float(np.log(abs(projection))/0.5))
                rows.append(dict(n=n,direction=direction,amplitude=amplitude,tau_end=0.5,
                                 expected=rate,measured=vals[-1],rtol_difference=abs(vals[-1]-vals[0]),
                                 error=abs(vals[-1]-rate)))
        assert max(z['error'] for z in rows if z['n']==n and z['amplitude']==1e-6)<1e-7
    save_csv('saddle_rates.csv',rows)
    return rows

# V52: article eq:numerical-crossing-point-v40 (460); article eq:numerical-crossing-initial-background-v40 (461); article tab:crossing-reproduction-v40 (5).
def crossing_numerics(kw):
    print('Checking logarithmic crossing solutions',flush=True)
    n=12; Lc=1.02; k=1.; chic_start=-1e-2
    def g(l): return 3*(-8*l**-12+12*l**-7)
    def g1(l): return 3*(96*l**-13-84*l**-8)
    def g2(l): return 3*(-1248*l**-14+672*l**-9)
    Hc=float(np.sqrt(g(Lc)/12)); cdotc=g1(Lc)*(Hc*Lc-1)
    K1=-3*cdotc/(Lc**2*k*k)
    assert cdotc>0
    rows=[]; histories=[]
    for tol in [1e-9,1e-12]:
        # Background parameter chi is monotone; its zero is imposed exactly.
        def background(ch, la):
            l,a=la; h=np.sqrt((g(l)-ch)/12)
            dotch=g1(l)*(h*l-1)-4*h*ch
            return [(h*l-1)/dotch,h*a/dotch]
        ini=integrate(background,(0,chic_start),[Lc,1.],rtol=tol,atol=tol*0.01).y[:,-1]
        l,a=ini; h=np.sqrt((g(l)-chic_start)/12)
        kin=kw(h,l,k*k/(a*a),chic_start,g1(l),g2(l))[0]
        P0=a**3*kin
        # Two fundamental initial conditions: (ell, dot ell)=(1,0),(0,1).
        init=[l,a,1.,0.,0.,P0]
        def rhs(u,state):
            ch=-np.exp(u);l,a,e1,j1,e2,j2=state
            h=np.sqrt((g(l)-ch)/12); qq=k*k/a**2
            cd=g1(l)*(h*l-1)-4*h*ch
            kin,w=kw(h,l,qq,ch,g1(l),g2(l))
            fac=ch/cd
            return [fac*(h*l-1),fac*h*a,fac*j1/(a**3*kin),-fac*a**3*w*e1,
                    fac*j2/(a**3*kin),-fac*a**3*w*e2]
        us=np.linspace(np.log(-chic_start),np.log(1e-12),1501)
        sol=integrate(rhs,(us[0],us[-1]),init,rtol=tol,atol=tol*0.01,t_eval=us)
        for scale in [1e-4,1e-8,1e-12]:
            idx=int(np.argmin(abs(sol.t-np.log(scale))))
            u=sol.t[idx]; state=sol.y[:,idx]; derivatives=rhs(u,state)
            for mode in [1,2]:
                ei,ji=(2,3) if mode==1 else (4,5)
                pred=sol.y[ji,-1]/K1
                rows.append(dict(rtol=tol,mode=mode,minus_chi=float(np.exp(u)),ell=float(state[ei]),
                                 d_ell_d_log_minus_chi=float(derivatives[ei]),
                                 flux_log_coefficient=float(pred),difference=abs(float(derivatives[ei]-pred))))
        if tol==1e-12:
            for i,u in enumerate(sol.t):
                histories.append(dict(minus_chi=float(np.exp(u)),L=sol.y[0,i],a=sol.y[1,i],
                                      ell1=sol.y[2,i],flux1=sol.y[3,i],ell2=sol.y[4,i],flux2=sol.y[5,i]))
            finalbg={'n':n,'dimensionless_equations':True,'normalization':'t*H0, L*H0, H/H0, k/H0, chi/(M^2*H0^2); not H0=M','Lc':Lc,'Hc':Hc,'k':k,'ac':1.,'dot_chi_c':cdotc,'K1':K1,
                     'initial_minus_chi':-chic_start,'initial_L':float(ini[0]),'initial_a':float(ini[1]),
                     'L_endpoint_error':abs(float(sol.y[0,-1]-Lc)),
                     'a_endpoint_error':abs(float(sol.y[1,-1]-1))}
    for mode in [1,2]:
        fine=[z for z in rows if z['rtol']==1e-12 and z['mode']==mode][-1]
        assert abs(fine['flux_log_coefficient'])>1e-5
        assert fine['difference']<1e-8
    save_csv('crossing_convergence.csv',rows);save_csv('crossing_history.csv',histories)
    return {'background':finalbg,'results':rows}

# V52: article eq:power-canonical-mass-limits (451).
def late_symbolic(n: int):
    d=sp.symbols('d',positive=True);x=1+d
    H=x**R(n,n-2);L=x**R(-2,n-2)
    eps=3*n*H**2*d/x
    # Differentiate G(L) first, then substitute the background.
    u=sp.symbols('L',positive=True)
    g=3*((4-n)*u**(-n)+n*u**R(-n-2,2))
    Delta=sp.factor((u*sp.diff(g,u)).subs(u,L)+eps)
    # A derivative with respect to log a, not partial differentiation of W.
    def DN(f):return sp.factor(-(n-2)*d*sp.diff(f,d)/2)
    A=1/(L**6*Delta)
    cr=4/L-eps/(6*H)
    vr=16/L**2+(Delta-4*eps)/(3*H*L)+eps**2/(36*H**2)-Delta/6
    nu=sp.factor(H*(R(3,2)+DN(A)/A/2))
    mass=sp.factor(-vr+2*nu*cr-nu**2+H*DN(cr)-H*DN(nu))
    cs2=sp.factor(Delta/(3*eps))
    expected= -R(1,4) if n==10 else -R((n-5)**2,4)
    check(f'mass_limit_n{n}',sp.limit(mass,d,0)-expected)
    if n==10:check('n10_exact_sound_speed',cs2-R(7,3))
    values=[H,L,eps,Delta,A,cr,nu,mass,cs2,DN(H)/H,DN(cs2)/cs2/2-1]
    EXPRESSIONS[f'late_n{n}_coefficients']=str(values)
    return sp.lambdify(d,values,'numpy',cse=True)

# V52: article eq:numerical-late-initial-data-v40 (462); article eq:numerical-curvature-envelope-v40 (463); article tab:late-rate-convergence-v40 (6).
def late_numerics():
    print('Checking late-time physical observables',flush=True)
    rows=[]
    for n in [10,12]:
        coeff=late_symbolic(n);d0=0.1;k=1.;end=8. if n==10 else 5.
        all_results=[]
        for tol in [1e-9,1e-12]:
            def rhs(N,state):
                d=d0*np.exp(-(n-2)*N/2)
                h,l,eps,delta,A,cr,nu,mass,cs2,hn,wn=coeff(d)
                w2=cs2*k*k*np.exp(-2*N)
                return [state[1],-hn*state[1]-(w2+mass)*state[0]/(h*h)]
            Ns=np.linspace(0,end,int(end*1200)+1)
            # Two fundamental real solutions, carried together; their norm has no zeros.
            def both(N,s):return rhs(N,s[:2])+rhs(N,s[2:])
            sol=integrate(both,(0,end),[1.,0.,0.,1.],rtol=tol,atol=tol*0.01,t_eval=Ns)
            observable=[]
            for i,N in enumerate(Ns):
                h,l,eps,delta,A,cr,nu,mass,cs2,hn,wn=coeff(d0*np.exp(-(n-2)*N/2))
                z=np.array([sol.y[0,i],sol.y[2,i]]);zp=np.array([sol.y[1,i],sol.y[3,i]])
                if n==10:
                    v=h*(h*zp+(cr-nu)*z)/(np.exp(1.5*N)*np.sqrt(A)*l*l*delta)
                    val=float(np.linalg.norm(v))
                else:
                    # Phase-space envelope for each real solution. Its corrections vanish
                    # in the adiabatic late limit; avoids taking derivatives through zeros.
                    omega=np.sqrt(cs2)*k*np.exp(-N)
                    zenv=np.sqrt(np.sum(z*z+(h*(zp+0.5*wn*z)/omega)**2))
                    val=4*h*k*k*np.exp(-2*N)*zenv/(np.exp(1.5*N)*np.sqrt(A)*eps*l**3)
                observable.append(val)
            obs=np.asarray(observable)
            for lo,hi in ([(4.,5.),(6.,7.),(7.,8.)] if n==10 else [(2.,3.),(3.,4.),(4.,5.)]):
                i0=np.argmin(abs(Ns-lo));i1=np.argmin(abs(Ns-hi))
                slope=float(np.log(obs[i1]/obs[i0])/(hi-lo))
                rows.append(dict(n=n,observable='norm_V' if n==10 else 'curvature_envelope',rtol=tol,
                                 log_a_start=lo,log_a_end=hi,expected=1. if n==10 else 3*(n-10)/8,
                                 measured=slope,absolute_asymptotic_error=abs(slope-(1. if n==10 else 3*(n-10)/8))))
            if tol==1e-12:
                save_csv(f'late_n{n}_history.csv',[dict(log_a=N,observable=obs[i],z1=sol.y[0,i],z1prime=sol.y[1,i],
                                                        z2=sol.y[2,i],z2prime=sol.y[3,i]) for i,N in enumerate(Ns)])
    assert abs(rows[-1]['measured']-0.75)<1e-5
    assert abs([z for z in rows if z['n']==10][-1]['measured']-1)<1e-5
    save_csv('late_rates.csv',rows)
    return rows

def main():
    import argparse
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/current')
    args=parser.parse_args()
    OUT=args.output/'n02_background_diagnostics'
    OUT.mkdir(parents=True,exist_ok=True)
    print('Deriving background and perturbation coefficients',flush=True)
    kw=symbolic_reduction()
    result={'environment':{'python':platform.python_version(),'sympy':sp.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
            'saddle':saddle_numerics(),'crossing':crossing_numerics(kw),'late':late_numerics()}
    result['symbolic_count']=len(CHECKS);result['symbolic_checks']=CHECKS
    (OUT/'numerical_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'coefficient_expressions.txt').write_text('\n\n'.join(k+'\n'+v for k,v in EXPRESSIONS.items()),encoding='utf-8')
    print('Saved numerical_results.json and convergence CSV files in '+str(OUT),flush=True)
    print(f'PASS: {len(CHECKS)} exact symbolic checks; all numerical assertions passed.')

if __name__=='__main__':main()
