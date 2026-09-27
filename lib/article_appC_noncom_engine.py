#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact ADM and finite-time verification engine for article V52, release V52-R1.
A minimal local covariant action for holographic dark energy: constraints,
perturbations, and nonlinear dynamics.
Preserved scientific algebra from original_scripts/noncom_original.py.
Appendices C.2-C.10; exact densities (342)-(349), fixed non-COM momenta (356),
finite-time integrals (359)-(363), scalar angular and nonlinear controls.
Method: multilinear nilpotent ADM polynomials over QQ(i,sqrt(2),sqrt(3));
exact SymPy verification for the explicitly named test functions below.
Inputs: Leg objects or embedded test configurations; no historical tables.
Outputs: caller-selected JSON objects; execute the themed scripts rather than
this legacy main entry point. See SCRIPT_REFERENCE_MAP.md for precise scope.
Assumptions: U=0 radiation, positive finite time endpoints, explicit tensor
norms, separately normalized homogeneous subsystem and state parameters.
No general IR-finite multichannel S matrix or regulator-free unitarity claim.
Dependencies: Python, SymPy. No LaTeX or network.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations, product
import sympy as sp

K=sp.QQ.algebraic_field(sp.I,sp.sqrt(2),sp.sqrt(3))
Z=K.zero; O=K.one

def alg(x):
    if isinstance(x,type(O)):return x
    return K.from_sympy(sp.sympify(x))
I=alg(sp.I); SQ3=alg(sp.sqrt(3))

class Ring:
    def __init__(self,momenta):
        self.n=len(momenta);self.all=(1<<self.n)-1
        self.p=[tuple(alg(x) for x in row) for row in momenta]
        self.ks={m:tuple(sum((self.p[j][i] for j in range(self.n) if m>>j&1),Z) for i in range(3)) for m in range(1<<self.n)}
    def c(self,x,power=0):
        x=alg(x)
        return Poly(self,{(0,power):x} if x else {})
    def leg(self,j,x=1,power=0):
        x=alg(x)
        return Poly(self,{(1<<j,power):x} if x else {})

class Poly:
    __slots__=('r','d')
    def __init__(self,r,d):self.r=r;self.d={k:v for k,v in d.items() if v}
    def __add__(self,b):
        if not isinstance(b,Poly):b=self.r.c(b)
        d=self.d.copy()
        for k,v in b.d.items():
            vv=d.get(k,Z)+v
            if vv:d[k]=vv
            else:d.pop(k,None)
        return Poly(self.r,d)
    __radd__=__add__
    def __neg__(self):return Poly(self.r,{k:-v for k,v in self.d.items()})
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return b+-self
    def __mul__(self,b):
        if not isinstance(b,Poly):
            b=alg(b)
            if not b:return self.r.c(0)
            return Poly(self.r,{k:v*b for k,v in self.d.items()})
        d={}
        for (a,p),v in self.d.items():
            for (bm,q),w in b.d.items():
                if a&bm:continue
                key=(a|bm,p+q)
                d[key]=d.get(key,Z)+v*w
        return Poly(self.r,d)
    __rmul__=__mul__
    def __truediv__(self,b):return self*(O/alg(b))
    def __pow__(self,n):
        if n<0:raise ValueError('Use series for inverse')
        out=self.r.c(1)
        for _ in range(n):out=out*self
        return out
    def dx(self,i):return Poly(self.r,{(m,p):v*I*self.r.ks[m][i] for (m,p),v in self.d.items()})
    def trunc(self,n):return Poly(self.r,{(m,p):v for (m,p),v in self.d.items() if m.bit_count()<=n})
    def take(self,mask=None):
        if mask is None:mask=self.r.all
        return {p:v for (m,p),v in self.d.items() if m==mask}
    def expr(self,mask=None):
        t=sp.Symbol('eta',positive=True)
        return sum((K.to_sympy(v)*t**p for p,v in self.take(mask).items()),sp.S.Zero)

def mat_zero(r):return [[r.c(0) for _ in range(3)] for __ in range(3)]
def mat_eye(r):return [[r.c(int(i==j)) for j in range(3)] for i in range(3)]
def mat_add(a,b):return [[a[i][j]+b[i][j] for j in range(3)] for i in range(3)]
def mat_scale(a,c):return [[a[i][j]*c for j in range(3)] for i in range(3)]
def mat_mul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(3)),a[0][0].r.c(0)) for j in range(3)] for i in range(3)]
def trace(a):return sum((a[i][i] for i in range(3)),a[0][0].r.c(0))

def norm2(p):return sum((x*x for x in p),Z)
def dot(p,q):return sum((a*b for a,b in zip(p,q)),Z)

@lru_cache(None)
def polarizations(p):
    """Real orthonormal TT tensors, with a deterministic transverse frame."""
    ps=sp.Matrix([K.to_sympy(x) if isinstance(x,type(O)) else sp.sympify(x) for x in p])
    k=sp.sqrt(sp.simplify(ps.dot(ps)))
    if k==0:raise ValueError('Homogeneous tensor is not a two-polarization wave')
    n=ps/k
    axis=next(sp.eye(3)[:,j] for j in range(3) if sp.simplify(1-n[j]**2)!=0)
    e=axis-n*(axis.dot(n));e=sp.simplify(e/sp.sqrt(sp.simplify(e.dot(e))))
    b=sp.simplify(n.cross(e))
    plus=sp.simplify((e*e.T-b*b.T)/sp.sqrt(2))
    cross=sp.simplify((e*b.T+b*e.T)/sp.sqrt(2))
    return [tuple(tuple(alg(m[i,j]) for j in range(3)) for i in range(3)) for m in (plus,cross)]

@dataclass(frozen=True)
class Leg:
    p:tuple
    omega:object
    kind:str='w'
    pol:tuple|None=None
    coordinate:object=1
    velocity:object|None=None


def lagrangian(legs:list[Leg],order:int|None=None):
    """Coefficient of epsilon_1 ... epsilon_n in the *unreduced* action,
    evaluated at the exact linear auxiliary solution for physical legs.
    Extra alpha/beta probe legs recover the second-order auxiliary sources.
    A coordinate=0, velocity=1 physical probe gives a Legendre source.
    """
    n=len(legs)
    if order is None:order=n
    r=Ring([leg.p for leg in legs]);zero=r.c(0);one=r.c(1)
    w=zero;u=zero;alpha=zero;beta=[zero for _ in range(3)]
    gamma=mat_zero(r);gammadot=mat_zero(r)
    for j,leg in enumerate(legs):
        c=alg(leg.coordinate);vel=alg(leg.velocity) if leg.velocity is not None else -I*alg(leg.omega)*c
        p=r.p[j];p2=norm2(p)
        if leg.kind=='w':
            w+=r.leg(j,c);u+=r.leg(j,vel)
            alpha+=r.leg(j,2*c,-1)
            if not p2:
                # Exact homogeneous linear lapse is 2 w_0', not 2 w_0/eta.
                alpha-=r.leg(j,2*c,-1)
                alpha+=r.leg(j,2*vel)
                continue
            for i in range(3):
                # -6 eta^-1 d_i Delta^-1 (w' - w/eta)
                beta[i]+=r.leg(j,6*I*p[i]*vel/p2,-1)+r.leg(j,-6*I*p[i]*c/p2,-2)
        elif leg.kind=='v':
            if leg.pol is None:raise ValueError('Missing tensor polarization')
            for a in range(3):
                for b in range(3):
                    e=alg(leg.pol[a][b])
                    gamma[a][b]+=r.leg(j,4*SQ3*e*c,-1)
                    gammadot[a][b]+=r.leg(j,4*SQ3*e*vel,-1)+r.leg(j,-4*SQ3*e*c,-2)
        elif leg.kind=='alpha':alpha+=r.leg(j,c)
        elif leg.kind.startswith('beta'):
            beta[int(leg.kind[-1])]+=r.leg(j,c)
        else:raise ValueError(leg.kind)
    powers=[mat_eye(r)]
    for m in range(1,order+1):powers.append(mat_mul(powers[-1],gamma))
    metric=mat_zero(r);inverse=mat_zero(r);metricdot=mat_zero(r)
    for m in range(order+1):
        fac=sp.factorial(m)
        metric=mat_add(metric,mat_scale(powers[m],sp.Rational(1,fac)))
        inverse=mat_add(inverse,mat_scale(powers[m],sp.Rational((-1)**m,fac)))
        if m:
            for j in range(m):
                term=mat_mul(mat_mul(powers[j],gammadot),powers[m-1-j])
                metricdot=mat_add(metricdot,mat_scale(term,sp.Rational(1,fac)))
    N=one+alpha
    invN=zero;invN3=zero
    for m in range(order+1):
        invN+=(-1)**m*alpha**m
        invN3+=(-1)**m*sp.binomial(m+2,2)*alpha**m
    # E_ij = Hconformal g_ij + 1/2(g'_ij-L_beta g_ij)
    covE=mat_zero(r)
    for i in range(3):
        for j in range(3):
            lie=sum((beta[k]*metric[i][j].dx(k)+metric[k][j]*beta[k].dx(i)+metric[i][k]*beta[k].dx(j) for k in range(3)),zero)
            covE[i][j]=metric[i][j]*r.c(1,-1)+(metricdot[i][j]-lie)/2
    E=mat_mul(inverse,covE)
    gravkin=trace(mat_mul(E,E))-trace(E)**2
    # Full spatial scalar curvature. No lapse-dependent integration by parts.
    Gamma=[[[zero for j in range(3)] for i in range(3)] for k in range(3)]
    for k,i,j in product(range(3),repeat=3):
        Gamma[k][i][j]=sum((inverse[k][l]*(metric[l][j].dx(i)+metric[l][i].dx(j)-metric[i][j].dx(l))/2 for l in range(3)),zero)
    curvature=zero
    for i,j in product(range(3),repeat=2):
        Ric=sum((Gamma[k][i][j].dx(k)-Gamma[k][i][k].dx(j) for k in range(3)),zero)
        Ric+=sum((Gamma[k][k][l]*Gamma[l][i][j]-Gamma[k][j][l]*Gamma[l][i][k] for k,l in product(range(3),repeat=2)),zero)
        curvature+=inverse[i][j]*Ric
    Lg=r.c(sp.Rational(1,24),2)*(N*curvature+invN*gravkin)
    grad=[w.dx(i) for i in range(3)]
    v=one+u-sum((beta[i]*grad[i] for i in range(3)),zero)
    grad2=sum((inverse[i][j]*grad[i]*grad[j] for i,j in product(range(3),repeat=2)),zero)
    Lm=(invN3*v**4-2*invN*v*v*grad2+N*grad2*grad2)/12
    return Lg+Lm



def scalar_variables(legs):
    r=Ring([leg.p for leg in legs]);z=r.c(0)
    w=z;u=z;rv=[z for _ in range(3)];T=mat_zero(r)
    for j,leg in enumerate(legs):
        if leg.kind!='w':continue
        c=alg(leg.coordinate);v=alg(leg.velocity) if leg.velocity is not None else -I*alg(leg.omega)*c
        p=r.p[j];p2=norm2(p)
        w+=r.leg(j,c);u+=r.leg(j,v)
        for a in range(3):
            rv[a]+=r.leg(j,-I*p[a]*v/p2)+r.leg(j,I*p[a]*c/p2,-1)
            for b in range(3):T[a][b]+=r.leg(j,p[a]*p[b]*v/p2)+r.leg(j,-p[a]*p[b]*c/p2,-1)
    g=[w.dx(i) for i in range(3)]
    g2=sum((x*x for x in g),z);rg=sum((x*y for x,y in zip(rv,g)),z);T2=sum((x*x for row in T for x in row),z)
    return r,w,u,rv,T,g,g2,rg,T2

def compact(legs,degree):
    r,w,u,rv,T,g,g2,rg,T2=scalar_variables(legs)
    if degree==3:
        return (u**3-u*g2)/3+r.c(1,-1)*(-3*w*T2+w*g2/3+6*u*rg)+r.c(1,-2)*(-12*w*rg-2*u*w*w)+r.c(sp.Rational(7,3),-3)*w**3
    if degree==4:
        return (u*u-g2)**2/12+r.c(1,-1)*(-2*g2*rg+sp.Rational(2,3)*g2*u*w+6*rg*u*u-2*u**3*w)+r.c(1,-2)*(6*T2*w*w-sp.Rational(2,3)*g2*w*w+18*rg**2-36*rg*u*w+6*u*u*w*w)+r.c(1,-3)*(48*rg*w*w-sp.Rational(20,3)*u*w**3)+r.c(2,-4)*w**4


import json,time
from pathlib import Path
from itertools import product

def la_add(a,b):
 d=a.copy()
 for p,v in b.items():
  v=d.get(p,Z)+v
  if v:d[p]=v
  else:d.pop(p,None)
 return d

def la_scale(a,c,power=0):
 c=alg(c)
 return {p+power:v*c for p,v in a.items() if v*c}

def la_mul(a,b):
 d={}
 for p,v in a.items():
  for q,w in b.items():d[p+q]=d.get(p+q,Z)+v*w
 return {p:v for p,v in d.items() if v}

def lc(a):return {str(p):str(sp.expand(K.to_sympy(v))) for p,v in sorted(a.items())}
def lp(a):
 t=sp.Symbol('eta',positive=True)
 return sp.simplify(sum(K.to_sympy(v)*t**p for p,v in a.items()))

VERTEX_CACHE={}
def Lvertex(legs):
 if legs not in VERTEX_CACHE:VERTEX_CACHE[legs]=lagrangian(list(legs)).take()
 return VERTEX_CACHE[legs]

@lru_cache(None)
def sources(la,lb):
 p=tuple(-alg(a)-alg(b) for a,b in zip(la.p,lb.p))
 sq=norm2(p)
 if not sq:raise ValueError('Homogeneous transfer requires independent constraints')
 j0=Lvertex((la,lb,Leg(p,0,'alpha')))
 ji=[Lvertex((la,lb,Leg(p,0,'beta'+str(i)))) for i in range(3)]
 jb={}
 for i in range(3):jb=la_add(jb,la_scale(ji[i],-I*p[i]/sq))
 jt=[la_add(ji[i],la_scale(jb,-I*p[i])) for i in range(3)]
 rw=Lvertex((la,lb,Leg(p,0,'w',None,0,1)))
 rv=[Lvertex((la,lb,Leg(p,0,'v',e,0,1))) for e in polarizations(p)]
 return j0,jb,jt,rw,rv

PAIRS=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))

def contact(legs):
 bare=Lvertex(tuple(legs));aux={};legendre={}
 parts=[]
 for ia,ib in PAIRS:
  aa,ab,jta,rwa,rva=sources(*(legs[i] for i in ia))
  ba,bb,jtb,rwb,rvb=sources(*(legs[i] for i in ib))
  pa=tuple(alg(legs[ia[0]].p[j])+alg(legs[ia[1]].p[j]) for j in range(3));sq=norm2(pa)
  one=la_scale(la_add(la_mul(aa,bb),la_mul(ba,ab)),-6,-1)
  one=la_add(one,la_scale(la_mul(ab,bb),-18,-2))
  for i in range(3):one=la_add(one,la_scale(la_mul(jta[i],jtb[i]),alg(24)/sq,-2))
  # Deterministic real frames differ by cross polarization sign under k -> -k.
  # Use one common frame at the two probes for contraction.
  es=polarizations(pa)
  rva=[Lvertex(tuple(legs[i] for i in ia)+(Leg(tuple(-x for x in pa),0,'v',e,0,1),)) for e in es]
  rvb=[Lvertex(tuple(legs[i] for i in ib)+(Leg(pa,0,'v',e,0,1),)) for e in es]
  lg=la_mul(rwa,rwb)
  for a,b in zip(rva,rvb):lg=la_add(lg,la_mul(a,b))
  aux=la_add(aux,one);legendre=la_add(legendre,lg)
  parts.append({'pair':ia,'auxiliary_H4':lc(one),'legendre_H4':lc(lg)})
 full=la_add(la_scale(bare,-1),la_add(aux,legendre))
 return {'bare_L4':lc(bare),'auxiliary_H4':lc(aux),'legendre_H4':lc(legendre),'full_H4':lc(full),'partitions':parts},full

def exchange(legs):
 result=[]
 for ia,ib in PAIRS:
  k=tuple(alg(legs[ia[0]].p[j])+alg(legs[ia[1]].p[j]) for j in range(3));sq=norm2(k)
  for A in range(3):
   nu=alg(sp.sqrt(K.to_sympy(sq)/(3 if A==0 else 1)))
   e=None if A==0 else polarizations(k)[A-1]
   kind='w' if A==0 else 'v'
   for late,early in ((ia,ib),(ib,ia)):
    pl=tuple(-alg(legs[late[0]].p[j])-alg(legs[late[1]].p[j]) for j in range(3))
    pe=tuple(-x for x in pl)
    vl=Lvertex(tuple(legs[i] for i in late)+(Leg(pl,nu,kind,e),))
    ve=Lvertex(tuple(legs[i] for i in early)+(Leg(pe,-nu,kind,e),))
    tensor={(r,s):v*w/(2*nu) for r,v in vl.items() for s,w in ve.items() if v*w}
    OmL=sum((alg(legs[i].omega) for i in late),Z)+nu
    OmE=sum((alg(legs[i].omega) for i in early),Z)-nu
    result.append({'pair_late':list(late),'pair_early':list(early),'internal':A,'frequency':str(K.to_sympy(nu)), 'A':str(sp.expand(K.to_sympy(OmL))),'B':str(sp.expand(K.to_sympy(OmE))), 'late_L3':lc(vl),'early_L3':lc(ve),'D_coefficients':{f'{r},{s}':str(sp.expand(K.to_sympy(v))) for (r,s),v in sorted(tensor.items())}})
 return result

P=[(0,1,sp.sqrt(3)),(0,1,-sp.sqrt(3)),(-sp.sqrt(3),-1,0),(sp.sqrt(3),-1,0)]

def sample_legs(species):
 return [Leg(tuple(alg(x) for x in p),alg((1 if i<2 else -1)*(2/sp.sqrt(3) if s==0 else 2)),'w' if s==0 else 'v',None if s==0 else polarizations(p)[s-1]) for i,(p,s) in enumerate(zip(P,species))]

def compute_all(out,start_index=0,stop_index=81):
 start=time.time();results=json.loads(out.read_text())["channels"] if out.exists() else {}

 for idx,s in enumerate(product(range(3),repeat=4)):
  if idx<start_index or idx>=stop_index:continue
  legs=sample_legs(s);ct,h4=contact(legs);ex=exchange(legs)
  results[''.join(map(str,s))]={'external_species':list(s),'contact':ct,'exchange':ex}
  out.write_text(json.dumps({'schema':'Finite-time connected tree kernels, not a partial-wave S matrix','units':'g_R=k0=1; eta_i, eta_f are parameters','momenta':[[str(x) for x in p] for p in P],'channels':results},ensure_ascii=False,indent=2),encoding='utf-8')
  print(len(results),s,'time',round(time.time()-start,2),'cache',len(VERTEX_CACHE),flush=True)
 print('DONE',time.time()-start)



# ----- Exact time integrals and independent checks -----
import ast
import argparse
import platform

CHECKS=[]

def check(name,expressions):
    if not isinstance(expressions,(list,tuple)):
        expressions=[expressions]
    residuals=[]
    for value in expressions:
        if isinstance(value,type(O)):
            value=K.to_sympy(value)
        value=sp.factor(sp.cancel(sp.expand(value)))
        if value!=0:value=sp.simplify(value)
        if value!=0:raise AssertionError(f'{name}: nonzero residual {value}')
        residuals.append('0')
    CHECKS.append({'name':name,'status':'PASS','residuals':residuals})


def parse_exact(text):
    """Parse only rational arithmetic, I and sqrt; never evaluate Python code."""
    def visit(node):
        if isinstance(node,ast.Expression):return visit(node.body)
        if isinstance(node,ast.Constant) and isinstance(node.value,int):return sp.Integer(node.value)
        if isinstance(node,ast.Name) and node.id=='I':return sp.I
        if isinstance(node,ast.UnaryOp):
            value=visit(node.operand)
            if isinstance(node.op,ast.USub):return -value
            if isinstance(node.op,ast.UAdd):return value
        if isinstance(node,ast.BinOp):
            a=visit(node.left);b=visit(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow):return a**b
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='sqrt' and len(node.args)==1 and not node.keywords:
            return sp.sqrt(visit(node.args[0]))
        raise ValueError(f'Unsupported expression in generated table: {text}')
    return sp.expand(visit(ast.parse(text,mode='eval')))


def L_power(z,T):
    z=sp.sympify(z)
    if z==0:return sp.log(T)
    if z.is_zero is False:return (T**z-1)/z
    return sp.Piecewise((sp.log(T),sp.Eq(z,0)),((T**z-1)/z,True))


def L_divided(u,v,T):
    u=sp.sympify(u);v=sp.sympify(v)
    derivative=sp.Piecewise((sp.log(T)**2/2,sp.Eq(u,0)),
        ((u*T**u*sp.log(T)-T**u+1)/u**2,True))
    if v==0:return derivative
    quotient=(L_power(u+v,T)-L_power(u,T))/v
    if v.is_zero is False:return quotient
    return sp.Piecewise((derivative,sp.Eq(v,0)),(quotient,True))


def exact_single(r,A,a,b):
    """Absolutely convergent exact series, not a numerical quadrature."""
    m=sp.Dummy('m',integer=True,nonnegative=True)
    return a**(r+1)*sp.Sum((-sp.I*A*a)**m/sp.factorial(m)*L_power(r+m+1,b/a),(m,0,sp.oo))


def exact_ordered(r,s,A,B,a,b):
    m=sp.Dummy('m',integer=True,nonnegative=True);n=sp.Dummy('n',integer=True,nonnegative=True)
    return a**(r+s+2)*sp.Sum((-sp.I*A*a)**m*(-sp.I*B*a)**n/sp.factorial(m)/sp.factorial(n)*L_divided(r+m+1,s+n+1,b/a),(m,0,sp.oo),(n,0,sp.oo))


def aggregate_exchange(rows):
    result={}
    for row in rows:
        key=(parse_exact(row['A']),parse_exact(row['B']))
        dest=result.setdefault(key,{})
        for pair,value in row['D_coefficients'].items():
            r,s=map(int,pair.split(','))
            dest[(r,s)]=sp.expand(dest.get((r,s),0)+parse_exact(value))
    return {key:{p:v for p,v in value.items() if v!=0} for key,value in result.items() if any(v!=0 for v in value.values())}


def inv_lap(poly):
    d={}
    for (mask,power),value in poly.d.items():
        sq=norm2(poly.r.ks[mask])
        if not sq:raise ValueError('Undefined homogeneous inverse Laplacian')
        d[(mask,power)]=-value/sq
    return Poly(poly.r,d)


# V52: article eq:v34r-exact-scalar-density (342); article eq:v34r-exact-gravity-density (343); article eq:v34r-extrinsic-density (344); article eq:v34r-nonzero-constraints (345).
def test_radiation_expansion():
    print("Running test_radiation_expansion", flush=True)
    t=sp.Symbol('eta',positive=True)
    p=(0,1,sp.sqrt(3));minus=tuple(-x for x in p)
    for kind,pol in [('w',None)]+[('v',e) for e in polarizations(p)]:
        ls=[Leg(p,3,kind,pol),Leg(minus,7,kind,pol)]
        expected=-21-sp.Rational(4,3)+20*sp.I/t+2/t**2 if kind=='w' else -25+10*sp.I/t+1/t**2
        check('exact_quadratic_'+kind+('_'+str(pol) if pol else ''),lagrangian(ls).expr()-expected)
    p1,p2,p3,p4=P;pm=tuple(-a-b for a,b in zip(p1,p2))
    for n,ls in [(3,[Leg(p1,3),Leg(p2,-2),Leg(pm,7)]),(4,[Leg(p1,3),Leg(p2,-2),Leg(p3,7),Leg(p4,4)])]:
        a=lagrangian(ls);b=compact(ls,n)
        check('independent_scalar_ADM_expansion_order_'+str(n),a.expr()-b.expr())
    alllegs=[Leg(p1,3),Leg(p2,-2),Leg(p3,7),Leg(p4,4)]
    for number,(ia,ib) in enumerate(PAIRS):
        ab=[alllegs[i] for i in ia];pm=tuple(-alg(ab[0].p[j])-alg(ab[1].p[j]) for j in range(3))
        r,w,u,rv,T,g,g2,rg,T2=scalar_variables(ab)
        ja=-sp.Rational(3,2)*T2+g2/6+r.c(1,-1)*(-6*rg+u*w)-r.c(sp.Rational(3,2),-2)*w*w
        ji=[-sum((T[i][j]*g[j] for j in range(3)),r.c(0))-r.c(sp.Rational(1,3),-1)*w*g[i] for i in range(3)]
        for kind,expected in zip(['alpha','beta0','beta1','beta2'],[ja]+ji):
            actual=lagrangian(ab+[Leg(pm,0,kind)]).expr()
            check(f'auxiliary_variation_{number}_{kind}',actual-expected.expr())
        rw=u*u-g2/3+r.c(6,-1)*rg-r.c(2,-2)*w*w
        for i in range(3):
            rw+=r.c(-6,-1)*inv_lap(u*g[i]).dx(i)+r.c(12,-2)*inv_lap(w*g[i]).dx(i)
            for j in range(3):rw+=r.c(-6,-1)*inv_lap(w*T[i][j]).dx(i).dx(j)
        actual=lagrangian(ab+[Leg(pm,0,'w',None,0,1)]).expr()
        check(f'scalar_Legendre_source_{number}',actual-rw.expr())
        for k,e in enumerate(polarizations(pm)):
            expected=-2*SQ3*r.c(1,-1)*sum((w*T[i][j]*e[i][j] for i in range(3) for j in range(3)),r.c(0))
            actual=lagrangian(ab+[Leg(pm,0,'v',e,0,1)]).expr()
            check(f'tensor_Legendre_source_{number}_{k}',actual-expected.expr())
            ls=ab+[Leg(pm,5,'v',e)]
            rr,ww,uu,rvv,TT,gg,gg2,rgg,TT2=scalar_variables(ls);z=rr.c(0)
            vv=[[rr.leg(2,e[i][j]) for j in range(3)] for i in range(3)]
            vp=[[rr.leg(2,-5*I*e[i][j]) for j in range(3)] for i in range(3)]
            expected=4*SQ3*sum((rr.c(1,-1)*vv[i][j]*gg[i]*gg[j]/6-rr.c(1,-1)*ww*TT[i][j]*vp[i][j]/2+rr.c(1,-2)*ww*TT[i][j]*vv[i][j]/2+sp.Rational(3,4)*rr.c(1,-1)*sum((vv[i][j].dx(a).dx(a) for a in range(3)),z)*rvv[i]*rvv[j] for i in range(3) for j in range(3)),z)
            check(f'complete_vww_without_time_integration_by_parts_{number}_{k}',lagrangian(ls).expr()-expected.expr())


# V52: article eq:v33-legendre-contact (191).
def test_general_Legendre_transform():
    print("Running test_general_Legendre_transform", flush=True)
    eps=sp.Symbol('eps');p,q=sp.symbols('p q');u,v=sp.symbols('u v')
    a,b,c,d,e,f=sp.symbols('a b c d e f')
    L3=a*u**3+b*u*u*v+c*u*v*v+d*v**3
    L4=e*u**4+f*u*u*v*v
    R=sp.Matrix([sp.diff(L3,u),sp.diff(L3,v)])
    sub={u:p,v:q};r=R.subs(sub)
    vel=sp.Matrix([p,q])-eps*r+eps**2*(R.jacobian([u,v]).subs(sub)*r-sp.Matrix([sp.diff(L4,u),sp.diff(L4,v)]).subs(sub))
    sv={u:vel[0],v:vel[1]}
    energy=p*vel[0]+q*vel[1]-(vel.dot(vel)/2+eps*L3.subs(sv,simultaneous=True)+eps**2*L4.subs(sv,simultaneous=True))
    actual=sp.expand(energy).coeff(eps,2)
    expected=-L4.subs(sub)+r.dot(r)/2
    check('two_field_Legendre_transform_independent',actual-expected)


# V52: article eq:radiation-noncom-check-momenta (356).
def test_complete_example():
    print("Running test_complete_example", flush=True)
    ls=sample_legs((0,0,0,0));ct,h=contact(ls)
    expected={'bare_L4':{0:sp.Rational(68,9),-2:sp.Rational(752,3),-4:504},
      'auxiliary_H4':{-2:sp.Rational(305,6),-4:sp.Rational(1275,8)},
      'legendre_H4':{0:24,-2:270,-4:sp.Rational(2955,4)},
      'full_H4':{0:sp.Rational(148,9),-2:sp.Rational(421,6),-4:sp.Rational(3153,8)}}
    for name,val in expected.items():
        actual={int(p):parse_exact(v) for p,v in ct[name].items()}
        check('qqqq_'+name,[actual.get(p,0)-val.get(p,0) for p in set(actual)|set(val)])
    rows=exchange(ls);raw=aggregate_exchange(rows)
    sq2=sp.sqrt(2);sq3=sp.sqrt(3);sq6=sp.sqrt(6);ii=sp.I
    polynomials=[
      (2*sq3,sq3/4,{0:40*sq3*ii/9,-1:34,-2:-34*sq3*ii,-3:-sp.Rational(53,2)}),
      (-2/sq3,sq3/4,{0:8*sq3*ii/9,-1:-6,-2:34*sq3*ii/3,-3:-sp.Rational(53,2)}),
      (sq2,sq2,{0:-2*sq2*ii,-1:-sp.Rational(47,6),-2:-31*sq2*ii/2,-3:-sp.Rational(215,8)}),
      (2+4/sq3,-sp.Rational(1,4),{-1:-5*sq6-6*sq2,-2:3*ii*(sq6+4*sq2),-3:15*sq6/4}),
      (2-4/sq3,-sp.Rational(1,4),{-1:5*sq6-6*sq2,-2:3*ii*(-sq6+4*sq2),-3:-15*sq6/4}),
      (sq6,-46/sq6,{-1:-5*sq6/12,-2:-3*ii/2,-3:-13*sq6/16})]
    wanted={}
    for A,factor,V in polynomials:
        key=(sp.expand(A),sp.expand(-A))
        wanted[key]={(r,s):sp.expand(factor*a*sp.conjugate(b)) for r,a in V.items() for s,b in V.items()}
    check('full_qqqq_exchange_six_polynomial_reduction',[
       raw.get(key,{}).get(p,0)-wanted.get(key,{}).get(p,0)
       for key in set(raw)|set(wanted) for p in set(raw.get(key,{}))|set(wanted.get(key,{}))])
    # Independent Feynman vs Hamiltonian stationary limit.
    feyn=Lvertex(tuple(ls)).get(0,Z);old=-h.get(0,Z)
    for ia,ib in PAIRS:
        k=tuple(alg(ls[ia[0]].p[j])+alg(ls[ia[1]].p[j]) for j in range(3))
        Om=sum((alg(ls[i].omega) for i in ia),Z);nu=alg(sp.sqrt(K.to_sympy(norm2(k))/3))
        ca=Lvertex(tuple(ls[i] for i in ia)+(Leg(tuple(-x for x in k),-Om),)).get(0,Z)
        cb=Lvertex(tuple(ls[i] for i in ib)+(Leg(k,Om),)).get(0,Z)
        feyn-=ca*cb/(Om*Om-nu*nu)
        for late,early in ((ia,ib),(ib,ia)):
            pl=tuple(-alg(ls[late[0]].p[j])-alg(ls[late[1]].p[j]) for j in range(3))
            va=Lvertex(tuple(ls[i] for i in late)+(Leg(pl,nu),)).get(0,Z)
            vb=Lvertex(tuple(ls[i] for i in early)+(Leg(tuple(-x for x in pl),-nu),)).get(0,Z)
            A=sum((alg(ls[i].omega) for i in late),Z)+nu
            old+=va*vb/(2*nu*A)
    check('local_Feynman_amplitude',feyn+alg(sp.Rational(52,27)))
    check('local_Hamiltonian_equals_Feynman',old-feyn)
    x,a,b=sp.symbols('x x_i x_f',positive=True)
    primitive=sp.Rational(148,9)*(b-a)+sp.Rational(421,6)*(1/a-1/b)+sp.Rational(1051,8)*(1/a**3-1/b**3)
    check('integrated_contact_upper_boundary',sp.diff(primitive,b)-(sp.Rational(148,9)+sp.Rational(421,6)/b**2+sp.Rational(3153,8)/b**4))
    check('integrated_contact_zero_interval',primitive.subs(b,a))
    p1=(sp.sqrt(2),0,1);p2=(-sp.sqrt(2),0,1);p3=(0,0,-2)
    for n,e in enumerate(polarizations(p3)):
        val=lagrangian([Leg(p1,1),Leg(p2,1),Leg(p3,-2,'v',e)]).expr()
        eta=sp.Symbol('eta',positive=True)
        wanted=4*sp.sqrt(6)/(3*eta)-4*sp.sqrt(6)*sp.I/eta**2-4*sp.sqrt(6)/eta**3 if n==0 else 0
        check('two_scalars_to_tensor_exact_'+str(n),val-wanted)
    primitive=4*sp.sqrt(6)/3*sp.log(b/a)+4*sp.sqrt(6)*sp.I*(1/b-1/a)+2*sp.sqrt(6)*(1/b**2-1/a**2)
    check('integrated_cubic_exact',sp.diff(primitive,b)-(4*sp.sqrt(6)/(3*b)-4*sp.sqrt(6)*sp.I/b**2-4*sp.sqrt(6)/b**3))
    pw=(sp.sqrt(3)/2,0,-sp.Rational(3,2));pa=(0,0,2);pb=(-sp.sqrt(3)/2,0,-sp.Rational(1,2))
    for i,e in enumerate(polarizations(pa)):
        for j,f in enumerate(polarizations(pb)):
            d=lagrangian([Leg(pw,-1),Leg(pa,2,'v',e),Leg(pb,-1,'v',f)]).take()
            check(f'on_shell_wvv_leading_cancellation_{i}{j}',d.get(-1,Z))
    return ct


# V52: article eq:v34-elementary-single (359); article eq:v34-elementary-ordered (360); article eq:v34-single-time-series (361); article eq:v34-ordered-time-series (362); article eq:v34-time-order-identity (363).
def test_time_integrals():
    print("Running test_time_integrals", flush=True)
    T=sp.Symbol('T',positive=True)
    for n in range(-6,5):
        L=L_power(n,T)
        check(f'exact_power_primitive_{n}',sp.diff(L,T)-T**(n-1))
        check(f'exact_power_lower_endpoint_{n}',L.subs(T,1))
    for r in range(-3,1):
        for s in range(-3,1):
            d=L_divided(r+1,s+1,T)
            check(f'ordered_moment_{r}_{s}',[sp.diff(d,T)-T**r*L_power(s+1,T),d.subs(T,1)])
            check(f'ordered_shuffle_{r}_{s}',d+L_divided(s+1,r+1,T)-L_power(r+1,T)*L_power(s+1,T))
    # The same calculation applies term by term to all Taylor indices.
    # Noninteger exponents test the divided-difference formula away from its removable poles.
    u,v=sp.symbols('u v',nonzero=True,real=True)
    d=( (T**(u+v)-1)/(u+v)-(T**u-1)/u )/v
    check('generic_ordered_moment_identity',sp.diff(d,T)-T**(u-1)*(T**v-1)/v)
    al,be=sp.symbols('alpha beta',complex=True);x,y,nu=sp.symbols('x y nu',real=True)
    left=(al*sp.exp(-sp.I*nu*x)+be*sp.exp(sp.I*nu*x))*(sp.conjugate(al)*sp.exp(sp.I*nu*y)+sp.conjugate(be)*sp.exp(-sp.I*nu*y))
    right=al*sp.conjugate(al)*sp.exp(-sp.I*nu*(x-y))+be*sp.conjugate(be)*sp.exp(sp.I*nu*(x-y))+al*sp.conjugate(be)*sp.exp(-sp.I*nu*(x+y))+be*sp.conjugate(al)*sp.exp(sp.I*nu*(x+y))
    check('Bogoliubov_four_frequency_weights',sp.expand(left-right,power_exp=True))


def audit_table(path):
    obj=json.loads(Path(path).read_text(encoding='utf-8'));channels=obj['channels']
    expected={''.join(map(str,s)) for s in product(range(3),repeat=4)}
    if set(channels)!=expected:raise AssertionError('The table does not contain exactly the 81 specified entries')
    count=0;groups=0
    for key,data in channels.items():
        ex=aggregate_exchange(data['exchange'])
        data['aggregated_exchange']=[{'A':str(A),'B':str(B),'D_coefficients':{f'{r},{s}':str(value) for (r,s),value in sorted(coeff.items())}} for (A,B),coeff in ex.items()]
        count+=bool(ex or data['contact']['full_H4']);groups+=len(ex)
        if len(data['exchange'])!=18:raise AssertionError('An exchange channel is missing')
        if not all(-4<=int(r)<=0 for r in data['contact']['full_H4']):raise AssertionError('Unexpected time degree')
        ss=list(map(int,key));Om=sum((1 if i<2 else -1)*(2/sp.sqrt(3) if s==0 else 2) for i,s in enumerate(ss))
        check('all_internal_frequencies_and_time_degrees_'+key,[parse_exact(e['A'])+parse_exact(e['B'])-Om for e in data['exchange']])
        rev=key[2:]+key[:2];phase=(-1)**sum(s!=0 for s in ss)
        a=data['contact']['full_H4'];b=channels[rev]['contact']['full_H4']
        check('contact_reflection_Hermiticity_'+key,[sp.conjugate(parse_exact(a.get(p,'0')))-phase*parse_exact(b.get(p,'0')) for p in set(a)|set(b)])
    check('number_of_entries',len(channels)-81)
    check('number_of_nonvanishing_entries_after_exact_aggregation',count-57)
    obj['counts']={'entries':81,'entries_not_vanishing_after_aggregation':count,'nonzero_ordered_frequency_groups':groups}
    obj['scope']='Reference conformal Fock basis; U=0 radiation, no independent matter; a fixed non-COM four-momentum configuration; connected tree 2-to-2 kernels only. Not an IR-finite partial-wave S matrix.'
    obj['units']='x=k0*eta. H4/(g_R^2*k0^4) and L3/(g_R*k0^3) have the listed Laurent coefficients. Full kernel and external normalization are given in the accompanying subsection.'
    obj['definitions']={'0':'canonical scalar w','1':'first real TT polarization','2':'second real TT polarization','key':'A1 A2 A3 A4, incoming slots 1,2, outgoing slots 3,4; all listed momenta are directed into vertices','F':'-i sum h_r I_r(Omega) - sum c_rs D_rs(A,B)','I_D':'the explicitly evaluated absolutely convergent series in the accompanying subsection, or exact_single/exact_ordered in this Python file'}
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
    return obj['counts']




@lru_cache(None)
def angular_polarizations(p):
    p=tuple(alg(x) for x in p)
    if not norm2(p):
        sq2=alg(sp.sqrt(2));sq6=alg(sp.sqrt(6))
        return [((O/sq2,Z,Z),(Z,-O/sq2,Z),(Z,Z,Z)),
                ((O/sq6,Z,Z),(Z,O/sq6,Z),(Z,Z,-2*O/sq6)),
                ((Z,O/sq2,Z),(O/sq2,Z,Z),(Z,Z,Z)),
                ((Z,Z,O/sq2),(Z,Z,Z),(O/sq2,Z,Z)),
                ((Z,Z,Z),(Z,Z,O/sq2),(Z,O/sq2,Z))]
    mag=MAGNITUDES[p]
    n=tuple(x/mag for x in p)
    e=(Z,O,Z);b=(-n[2],Z,n[0])
    sq2=alg(sp.sqrt(2))
    plus=tuple(tuple((e[i]*e[j]-b[i]*b[j])/sq2 for j in range(3)) for i in range(3))
    cross=tuple(tuple((e[i]*b[j]+b[i]*e[j])/sq2 for j in range(3)) for i in range(3))
    return (plus,cross)

@lru_cache(None)
def angular_sources(la,lb):
    p=tuple(-alg(a)-alg(b) for a,b in zip(la.p,lb.p))
    sq=norm2(p)
    j0=Lvertex((la,lb,Leg(p,0,'alpha')))
    ji=[Lvertex((la,lb,Leg(p,0,'beta'+str(i)))) for i in range(3)]
    if not sq:
        if any(ji):raise AssertionError('Nonzero total-momentum source on homogeneous shift')
        jb={};jt=[{}, {}, {}]
    else:
        jb={}
        for i in range(3):jb=la_add(jb,la_scale(ji[i],-I*p[i]/sq))
        jt=[la_add(ji[i],la_scale(jb,-I*p[i])) for i in range(3)]
    rw=Lvertex((la,lb,Leg(p,0,'w',None,0,1)))
    rv=[Lvertex((la,lb,Leg(p,0,'v',e,0,1))) for e in angular_polarizations(p)]
    return j0,jb,jt,rw,rv

def angular_contact(legs):
    bare=Lvertex(tuple(legs));aux={};legendre={};parts=[]
    for ia,ib in PAIRS:
        aa,ab,jta,rwa,rva=angular_sources(*(legs[i] for i in ia))
        ba,bb,jtb,rwb,rvb=angular_sources(*(legs[i] for i in ib))
        pa=tuple(alg(legs[ia[0]].p[j])+alg(legs[ia[1]].p[j]) for j in range(3));sq=norm2(pa)
        if sq:
            one=la_scale(la_add(la_mul(aa,bb),la_mul(ba,ab)),-6,-1)
            one=la_add(one,la_scale(la_mul(ab,bb),-18,-2))
            for i in range(3):one=la_add(one,la_scale(la_mul(jta[i],jtb[i]),alg(24)/sq,-2))
        else:
            # H_aa = 1/2, H4_aux = + J_0^2. Mixed pair coefficient is 2 Ja Jb.
            one=la_scale(la_mul(aa,ba),2)
        es=angular_polarizations(pa)
        rva=[Lvertex(tuple(legs[i] for i in ia)+(Leg(tuple(-x for x in pa),0,'v',e,0,1),)) for e in es]
        rvb=[Lvertex(tuple(legs[i] for i in ib)+(Leg(pa,0,'v',e,0,1),)) for e in es]
        lg=la_scale(la_mul(rwa,rwb),1 if sq else -1)
        for a,b in zip(rva,rvb):lg=la_add(lg,la_mul(a,b))
        aux=la_add(aux,one);legendre=la_add(legendre,lg)
        parts.append({'pair':ia,'homogeneous':not bool(sq),'auxiliary_H4':lc(one),'legendre_H4':lc(lg)})
    full=la_add(la_scale(bare,-1),la_add(aux,legendre))
    return {'bare_L4':lc(bare),'auxiliary_H4':lc(aux),'legendre_H4':lc(legendre),'full_H4':lc(full),'partitions':parts},full

def angular_exchange(legs):
    result=[];zero=[]
    for ia,ib in PAIRS:
        k=tuple(alg(legs[ia[0]].p[j])+alg(legs[ia[1]].p[j]) for j in range(3));sq=norm2(k)
        if not sq:
            for A in range(6):
                e=None if A==0 else angular_polarizations(k)[A-1]
                kind='w' if A==0 else 'v'
                for late,early in ((ia,ib),(ib,ia)):
                    al=Lvertex(tuple(legs[i] for i in late)+(Leg(k,0,kind,e,1,0),))
                    bl=Lvertex(tuple(legs[i] for i in late)+(Leg(k,0,kind,e,0,1),))
                    ae=Lvertex(tuple(legs[i] for i in early)+(Leg(k,0,kind,e,1,0),))
                    be=Lvertex(tuple(legs[i] for i in early)+(Leg(k,0,kind,e,0,1),))
                    zero.append({'late':list(late),'early':list(early),'global_index':A,'kinetic_sign':-1 if A==0 else 1,
                        'A':str(K.to_sympy(sum((legs[i].omega for i in late),Z))),
                        'B':str(K.to_sympy(sum((legs[i].omega for i in early),Z))),
                        'late_coordinate':lc(al),'late_velocity':lc(bl),'early_coordinate':lc(ae),'early_velocity':lc(be)})
            continue
        for A in range(3):
            nu=MAGNITUDES[k]/(SQ3 if A==0 else O)
            e=None if A==0 else angular_polarizations(k)[A-1];kind='w' if A==0 else 'v'
            for late,early in ((ia,ib),(ib,ia)):
                pl=tuple(-legs[late[0]].p[j]-legs[late[1]].p[j] for j in range(3));pe=tuple(-x for x in pl)
                vl=Lvertex(tuple(legs[i] for i in late)+(Leg(pl,nu,kind,e),))
                ve=Lvertex(tuple(legs[i] for i in early)+(Leg(pe,-nu,kind,e),))
                tensor={(r,s):v*w/(2*nu) for r,v in vl.items() for s,w in ve.items() if v*w}
                OmL=sum((legs[i].omega for i in late),Z)+nu;OmE=sum((legs[i].omega for i in early),Z)-nu
                result.append({'late':list(late),'early':list(early),'internal':A,'frequency':str(K.to_sympy(nu)),
                    'A':str(K.to_sympy(OmL)),'B':str(K.to_sympy(OmE)),
                    'late_L3':lc(vl),'early_L3':lc(ve),
                    'D_coefficients':{f'{r},{s}':str(K.to_sympy(v)) for (r,s),v in sorted(tensor.items())}})
    return result,zero



# The homogeneous lapse solution used here is N_1=2*g_R*w_0', not the
# nonzero-momentum solution 2*g_R*w/eta. These are distinct reductions.
# No call to a homogeneous inverse Laplacian is made.

def initialize_com(z):
    """Exact algebraic COM configuration; all four momenta point inward."""
    global MAGNITUDES
    z=sp.sympify(z);sn=sp.sqrt(1-z*z)
    momenta=[(0,0,1),(0,0,-1),(-sn,0,-z),(sn,0,z)]
    MAGNITUDES={}
    for p in momenta:
        p=tuple(alg(x) for x in p)
        MAGNITUDES[p]=O;MAGNITUDES[tuple(-x for x in p)]=O
    for ia,ib in PAIRS:
        p=tuple(alg(momenta[ia[0]][j]+momenta[ia[1]][j]) for j in range(3))
        if norm2(p):
            mag=alg(sp.sqrt(K.to_sympy(norm2(p))))
            MAGNITUDES[p]=mag;MAGNITUDES[tuple(-x for x in p)]=mag
    return [Leg(tuple(alg(x) for x in p),alg((1 if j<2 else -1)/sp.sqrt(3)))
            for j,p in enumerate(momenta)]


# V52: article eq:radiation-contact-decomposition (351); article eq:v34-scalar-equal-radii-contact (358).
def test_direct_com():
    print("Running test_direct_com", flush=True)
    """Independent full ADM evaluations, not tests against a supplied table."""
    t=sp.Symbol('eta',positive=True);points=[]
    for z in [sp.S.Zero,sp.Rational(1,3),-sp.Rational(1,3)]:
        legs=initialize_com(z)
        ct,h=angular_contact(legs)
        wanted=sp.Rational(4,3)-4*z*z/9+(64+16*z*z-32/(1-z*z))/t**2+(368+24*z*z)/t**4
        check('new_full_COM_contact_'+str(z),lp(h)-wanted)
        for j,pair in enumerate(PAIRS[0]):
            so=angular_sources(*(legs[i] for i in pair));sign=1 if j==0 else -1
            check(f'new_global_lapse_{z}_{j}',lp(so[0])-(sp.Rational(4,3)-sign*20*sp.I/sp.sqrt(3)/t-18/t**2))
            check(f'new_global_scalar_momentum_{z}_{j}',lp(so[3])-(sp.Rational(4,3)-sign*16*sp.I/sp.sqrt(3)/t-8/t**2))
        ex,ze=angular_exchange(legs)
        for row in ex:
            late=row['late'];zz=z if set(late) in [{0,2},{1,3}] else -z
            nu=parse_exact(row['frequency'])
            val=sum(parse_exact(v)*t**int(r) for r,v in row['late_L3'].items())
            if row['internal']==0:
                target=-2*sp.I*nu*(1-zz)/3+(-2*zz**2+sp.Rational(10,3)*zz-sp.Rational(8,3))/t-sp.I*nu*(6*zz+14)/t**2+(-6*zz**2+6*zz-28)/t**3
            else:
                p=tuple(legs[late[0]].p[j]+legs[late[1]].p[j] for j in range(3))
                pol=angular_polarizations(p)[row['internal']-1]
                a=next(j for j in late if j<2);v=legs[a].p
                ep=K.to_sympy(sum((pol[i][j]*v[i]*v[j] for i in range(3) for j in range(3)),Z))
                target=4*sp.sqrt(3)*((zz-sp.Rational(2,3))/t-sp.I*nu/t**2+(3*zz-4)/t**3)*ep
            check(f'new_nonzero_vertex_{z}_{late}_{row["internal"]}',val-target)
        for row in ze:
            if row['global_index']==0:
                check(f'new_global_scalar_coordinate_{z}_{row["late"]}',[parse_exact(v) for v in row['late_coordinate'].values()])
                continue
            late=row['late'];v=legs[late[0]].p
            pol=angular_polarizations((0,0,0))[row['global_index']-1]
            ep=K.to_sympy(sum((pol[i][j]*v[i]*v[j] for i in range(3) for j in range(3)),Z))
            sign=1 if late==[0,1] else -1
            A=4/sp.sqrt(3)/t-sign*4*sp.I/t**2-4*sp.sqrt(3)/t**3
            B=sign*4*sp.I/t+4*sp.sqrt(3)/t**2
            for key,target in [('late_coordinate',A*ep),('late_velocity',B*ep)]:
                val=sum(parse_exact(v)*t**int(r) for r,v in row[key].items())
                check(f'new_global_tensor_{z}_{late}_{row["global_index"]}_{key}',val-target)
        points.append({'cos_theta':str(z),'contact':ct})
    return points


def invariant_scalar_cubic(velocities,squared_momenta,gram,h):
    """Multilinear coefficient of the exact compact scalar cubic action."""
    r=[a-h for a in velocities]
    out=2*sp.prod(velocities)
    for a,b,c in ((0,1,2),(1,0,2),(2,0,1)):
        out+=sp.Rational(2,3)*velocities[a]*gram[b,c]
        out-=6*h*r[b]*r[c]*gram[b,c]**2/(squared_momenta[b]*squared_momenta[c])
        out-=sp.Rational(2,3)*h*gram[b,c]
    for a,b,c in permutations(range(3)):
        out+=6*h*velocities[a]*r[b]*gram[b,c]/squared_momenta[b]
        out-=12*h*h*r[b]*gram[b,c]/squared_momenta[b]
    return sp.factor(out-4*h*h*sum(velocities)+14*h**3)


# V52: article eq:v34-scalar-leading-contact (357); article eq:v34-scalar-equal-radii-contact (358).
def test_generic_scalar_contact():
    print("Running test_generic_scalar_contact", flush=True)
    """All-angle derivation by Gram matrices, independent of direct ADM points."""
    h,z=sp.symbols('h z',real=True);c=1/sp.sqrt(3);ii=sp.I
    D=sp.Matrix([[1,-1,-z,z],[-1,1,z,-z],[-z,z,1,-1],[z,-z,-1,1]])
    u=[-ii*c,-ii*c,ii*c,ii*c];r=[a-h for a in u]
    bare=2*sp.prod(u)
    for a,b,d,e in permutations(range(4)):
        bare+=u[a]*u[b]*D[d,e]/6+D[a,b]*D[d,e]/12
        bare+=h*(2*D[a,b]*r[d]*D[d,e]-sp.Rational(2,3)*D[a,b]*u[d]+6*r[a]*D[a,b]*u[d]*u[e]-2*u[a]*u[b]*u[d])
        bare+=h*h*(6*r[a]*r[b]*D[a,b]**2+sp.Rational(2,3)*D[a,b]+18*r[a]*D[a,b]*r[d]*D[d,e]-36*r[a]*D[a,b]*u[d]+6*u[a]*u[b])
        bare+=h**3*(48*r[a]*D[a,b]-sp.Rational(20,3)*u[a])
    bare=sp.factor(bare+48*h**4)
    check('new_generic_bare_L4',bare-(sp.Rational(4,3)*(1+z*z)+sp.Rational(320,3)*h*h+(624+384*z*z)*h**4))
    ua=-ii*c;ub=ii*c;ra=ua-h;rb=ub-h;Dab=-z;q2=2*(1-z)
    Ja=sp.factor(-3*ra*rb*Dab**2-Dab/3+h*(-6*Dab*(ra+rb)+ua+ub)-3*h*h)
    JB=sp.factor(Dab*(ra+rb)*(1+Dab)/q2+h/3)
    JTprod=-sp.Rational(2,3)*z*z*(1+z)
    auxiliary=sp.factor(-12*h*Ja*JB-18*h*h*JB**2+24*h*h*JTprod/q2)
    ui=sp.Symbol('ui');nu=sp.Symbol('nu',positive=True)
    Di=sp.Matrix([[1,-z,z-1],[-z,1,z-1],[z-1,z-1,q2]])
    cubic=invariant_scalar_cubic([ua,ub,ui],[1,1,q2],Di,h)
    rw=sp.diff(cubic,ui)
    legendre=sp.factor(rw**2+6*h**4*(1+z)**2)
    j0=sp.Rational(4,3)-20*ii*c*h-18*h*h
    r0=sp.Rational(4,3)-16*ii*c*h-8*h*h
    aux0=2*j0*sp.conjugate(j0)
    leg0=-r0*sp.conjugate(r0)+(16*h*h+48*h**4)*(z*z-sp.Rational(1,3))
    total=sp.factor(-bare+auxiliary+auxiliary.subs(z,-z)+legendre+legendre.subs(z,-z)+aux0+leg0)
    target=sp.Rational(4,3)-4*z*z/9+(64+16*z*z-32/(1-z*z))*h*h+(368+24*z*z)*h**4
    check('new_generic_full_H4',total-target)
    vs=-2*ii*nu*(1-z)/3+(-2*z*z+sp.Rational(10,3)*z-sp.Rational(8,3))*h-ii*nu*(6*z+14)*h*h+(-6*z*z+6*z-28)*h**3
    check('new_generic_scalar_exchange',cubic.subs(ui,-ii*nu)-vs)
    vt=4*sp.sqrt(3)*(h/3-ii*nu*h*h-h**3-sp.Rational(3,2)*q2*h*ra*rb)
    check('new_generic_tensor_exchange',vt-4*sp.sqrt(3)*((z-sp.Rational(2,3))*h-ii*nu*h*h+(3*z-4)*h**3))
    lead=-target.subs(h,0)+sp.Rational(4,9)*((1-z)**2+(1+z)**2)
    check('new_COM_leading_amplitude',lead-sp.Rational(4,9)*(3*z*z-1))
    check('new_leading_covariant_independent',sp.Rational(4,3)*(1+z*z)-sp.Rational(16,9)-lead)
    return {'bare_L4':str(sp.expand(bare)),'full_H4':str(target),'scalar_vertex':str(vs),'tensor_vertex':str(vt),'leading_amplitude':str(sp.factor(lead))}


# V52: article eq:v33-legendre-contact (191); article eq:v34r-global-reduction (348).
def test_signed_Legendre_transform():
    """Independent Legendre inversion with one negative global kinetic direction."""
    eps=sp.Symbol('eps');p,q,u,v=sp.symbols('p q u v')
    a,b,c,d,e,f=sp.symbols('a b c d e f')
    Kinv=sp.diag(-1,1);vel0=Kinv*sp.Matrix([p,q])
    L3=a*u**3+b*u*u*v+c*u*v*v+d*v**3;L4=e*u**4+f*u*u*v*v
    grad=sp.Matrix([sp.diff(L3,u),sp.diff(L3,v)]);sub={u:vel0[0],v:vel0[1]}
    R=grad.subs(sub);vel1=-Kinv*R
    grad4=sp.Matrix([sp.diff(L4,u),sp.diff(L4,v)]).subs(sub)
    vel2=-Kinv*(grad.jacobian([u,v]).subs(sub)*vel1+grad4)
    vel=vel0+eps*vel1+eps**2*vel2
    sv={u:vel[0],v:vel[1]}
    energy=sp.Matrix([p,q]).dot(vel)-(vel.T*Kinv*vel)[0]/2-eps*L3.subs(sv,simultaneous=True)-eps**2*L4.subs(sv,simultaneous=True)
    expected=-L4.subs(sub)+(R.T*Kinv*R)[0]/2
    check('new_signed_Legendre_transform',sp.expand(energy).coeff(eps,2)-expected)


# V52: article eq:v34r-global-free-flow (400); article eq:v34-global-contraction (402).
def test_homogeneous_and_state():
    print("Running test_homogeneous_and_state", flush=True)
    n,g,u,b=sp.symbols('N g u b',positive=True);a=1+g*u
    lag=(a**4/n**3-3*b/n)/(12*g*g);nsol=a*a/sp.sqrt(b)
    red=-b**sp.Rational(3,2)/(6*g*g*a*a)
    check('new_exact_homogeneous_lapse',sp.diff(lag,n).subs(n,nsol))
    check('new_exact_homogeneous_reduction',lag.subs(n,nsol)-red)
    target=-1/(6*g*g)+u/(3*g)-u*u/2+2*g*u**3/3-5*g*g*u**4/6
    check('new_homogeneous_scalar_quartic_expansion',sp.series(red.subs(b,1),u,0,5).removeO()-target)
    check('new_homogeneous_lapse_Hessian',sp.diff(lag,n,2).subs({n:1,u:0,b:1})-1/(2*g*g))
    pols=angular_polarizations((0,0,0))
    check('new_five_global_tensors_orthonormal',[
       K.to_sympy(sum((e[i][j]*f[i][j] for i in range(3) for j in range(3)),Z))-int(a==b)
       for a,e in enumerate(pols) for b,f in enumerate(pols)])
    projector=[]
    for i,j,k,l in product(range(3),repeat=4):
        target=sp.Rational(1,2)*(int(i==k)*int(j==l)+int(i==l)*int(j==k))-sp.Rational(1,3)*int(i==j)*int(k==l)
        projector.append(K.to_sympy(sum((e[i][j]*e[k][l] for e in pols),Z))-target)
    check('new_global_traceless_projector',projector)
    x,y,xi,sg=sp.symbols('x y x_i sigma',real=True)
    cq,cp,cm=sp.symbols('Cqq Cpp Cqp',real=True)
    W=cq+sg*(x+y-2*xi)*cm+(x-xi)*(y-xi)*cp-sp.I*sg*(x-y)/2
    check('new_free_global_correlator_equation',sp.diff(W,x,2))
    check('new_free_global_correlator_Hermiticity',W-sp.conjugate(W.subs({x:y,y:x},simultaneous=True)))
    check('new_global_equal_time_commutator',
          (sp.diff(W,y)-sp.diff(W,x)).subs(y,x)-sp.I*sg)
    # Two time orderings. The covariance terms fill a rectangle; the remaining
    # commutator is i Re(E). Verify its coefficient algebra independently.
    fx,bx,fy,by=sp.symbols('fx bx fy by',complex=True)
    E=fx*sp.conjugate(by)-bx*sp.conjugate(fy)
    second=sp.conjugate(fx)*by-sp.conjugate(bx)*fy
    check('new_two_orderings_commutator',sp.expand_complex(sp.I*(E+second)/2-sp.I*sp.re(E)))
    return {'homogeneous_lapse':str(nsol),'homogeneous_reduced_density':str(red),'homogeneous_Wightman':str(W)}


def polynomial_integral(poly,z):
    p=sp.Poly(sp.expand(poly),z)
    return sp.expand(sum(coef*(1-(-1)**(power[0]+1))/(power[0]+1) for power,coef in p.terms()))


def angular_polynomial_coefficients(j):
    """Computed finite polynomials for every integer j >= 0, no differentiation of action left."""
    if not isinstance(j,int) or j<0:raise ValueError('j must be a nonnegative integer')
    q=sp.Symbol('q',real=True)
    v={0:-sp.I*q**3/(3*sp.sqrt(3)),
       -1:-sp.Rational(4,3)+q*q/3-q**4/2,
       -2:-sp.I*(20*q-3*q**3)/sp.sqrt(3),
       -3:-28+3*q*q-sp.Rational(3,2)*q**4}
    t={-1:4*sp.sqrt(3)/3-2*sp.sqrt(3)*q*q,
       -2:-4*sp.I*sp.sqrt(3)*q,
       -3:-4*sp.sqrt(3)-6*sp.sqrt(3)*q*q}
    records=[]
    for name,V,weight in [('scalar',v,sp.sqrt(3)),('tensor',t,(2-q*q/2)**2/8)]:
        for r,a in V.items():
            for s,b in V.items():
                pol=sp.Poly(sp.expand(sp.legendre(j,1-q*q/2)*weight*a*sp.conjugate(b)),q)
                records.append({'internal':name,'r':r,'s':s,'coefficients':{str(p[0]):str(c) for p,c in pol.terms() if c!=0}})
    return records


# V52: article eq:v34-angular-exchange-series (374); article eq:v34-angular-contact-series (375).
def test_angular_projection():
    print("Running test_angular_projection", flush=True)
    z,q,h=sp.symbols('z q h',real=True);xx,yy=sp.symbols('x y',positive=True)
    contact={};coefficients={}
    for j in range(9):
        Pj=sp.legendre(j,z);m0=polynomial_integral(Pj,z);m2=polynomial_integral(z*z*Pj,z)
        check(f'new_contact_angular_moments_{j}',[m0-2*int(j==0),m2-sp.Rational(2,3)*int(j==0)-sp.Rational(4,15)*int(j==2)])
        remainder=sp.cancel((Pj-1)/(1-z))
        check(f'new_endpoint_log_finite_part_{j}',polynomial_integral(remainder,z)+2*sp.harmonic(j))
        h0=sp.Rational(4,3)*m0-sp.Rational(4,9)*m2
        h2=64*m0+16*m2+32*(1+(-1)**j)*sp.harmonic(j)
        h4=368*m0+24*m2
        contact[str(j)]={'h0':str(h0),'h2_finite_part':str(h2),'h4':str(h4),'log_coefficient_H4':str(-16*(1+(-1)**j))}
        records=angular_polynomial_coefficients(j);coefficients[str(j)]=records
        # Polynomial integration is independently evaluated term-by-term and by
        # an antiderivative. All angular and Taylor coefficients are exact.
        residuals=[]
        for row in records:
            pol=sum(parse_exact(c)*q**int(m) for m,c in row['coefficients'].items())
            for n in range(3):
                formula=sum(parse_exact(c)*2**(int(m)+n+1)/sp.Integer(int(m)+n+1) for m,c in row['coefficients'].items())
                primitive=sp.Poly(pol*q**n,q).integrate().as_expr()
                residuals.append(primitive.subs(q,2)-primitive.subs(q,0)-formula)
        check(f'new_all_exchange_angular_moments_j{j}_n0_to2',residuals)
    for n in range(8):
        check(f'new_double_time_binomial_{n}',(xx-yy)**n-sum((-1)**l*sp.binomial(n,l)*xx**(n-l)*yy**l for l in range(n+1)))
    check('new_contact_j0_values',[parse_exact(contact['0'][key])-v for key,v in [('h0',sp.Rational(64,27)),('h2_finite_part',sp.Rational(416,3)),('h4',752)]])
    check('new_contact_j2_values',[parse_exact(contact['2'][key])-v for key,v in [('h0',-sp.Rational(16,135)),('h2_finite_part',sp.Rational(1504,15)),('h4',sp.Rational(32,5))]])
    return contact,coefficients


def legendre_Q_explicit(j,u):
    return sp.legendre(j,u)*sp.log((u+1)/(u-1))/2-sum(sp.legendre(r-1,u)*sp.legendre(j-r,u)/r for r in range(1,j+1))


def radial_transverse_moment(j,u):
    Q=legendre_Q_explicit(j,u)
    m0=2*int(j==0);m1=sp.Rational(2,3)*int(j==1)
    m2=sp.Rational(2,3)*int(j==0)+sp.Rational(4,15)*int(j==2)
    return -2*((2*u-4*u**3)*Q+(u*u-u**4)*sp.diff(Q,u))-m2-2*u*m1-(3*u*u-1)*m0


# V52: article eq:radiation-cubic-source-definition (352).
def test_radial_kernel():
    print("Running test_radial_kernel", flush=True)
    k,l,z,h=sp.symbols('k l z h',real=True,positive=True);c=1/sp.sqrt(3)
    a=c*k*l*z-sp.I*h*(l*z+k/3);b=c*k*l*z+sp.I*h*(k*z+l/3)
    q2=k*k+l*l-2*k*l*z
    j2=a*a+b*b+2*z*a*b;qj=k*a-l*b+z*(k*b-l*a)
    transverse=sp.factor(j2-qj*qj/q2)
    target=(k+l)**2*z*z*(1-z*z)*(c*k*l+sp.I*h*(k-l))**2/q2
    check('new_general_radial_transverse_projection',transverse-target)
    H=-24*h*h*transverse/q2
    check('new_equal_radius_transverse_contact',H.subs(l,k)+8*h*h*k*k*z*z*(1+z)/(1-z))
    check('new_forward_residue',sp.limit((1-z)*H.subs(l,k),z,1)+16*h*h*k*k)
    u=sp.Symbol('upsilon',positive=True);z=sp.Symbol('z',real=True)
    data={}
    for j in range(9):
        P=sp.legendre(j,z);Pu=sp.legendre(j,u)
        Q=legendre_Q_explicit(j,u)
        direct=Pu*sp.log((u+1)/(u-1))/2+polynomial_integral(sp.cancel((P-Pu)/(u-z)),z)/2
        check(f'new_Legendre_Q_finite_formula_{j}',Q-direct)
        F=z*z*(1-z*z)*P;Fu=F.subs(z,u)
        cauchy=Fu*sp.log((u+1)/(u-1))+polynomial_integral(sp.cancel((F-Fu)/(u-z)),z)
        R=radial_transverse_moment(j,u)
        check(f'new_radial_angular_integral_{j}',R+sp.diff(cauchy,u))
        data[str(j)]=str(R)
    return data


# V52: article eq:v33-channel-counts (194); article eq:v34-finite-time-leakage (411).
def test_channel_count_and_leakage():
    print("Running test_channel_count_and_leakage", flush=True)
    def count(nscalar,j):
        helicities=[0]*nscalar+[2,-2];result=[]
        for a in range(len(helicities)):
            for b in range(a,len(helicities)):
                if abs(helicities[a]-helicities[b])>j:continue
                if a==b and j%2:continue
                result.append((a,b))
        return result
    data={}
    for n in [1,2]:
        data[str(n)]={}
        for j in range(10):
            pairs=count(n,j)
            expected=(3*int(j%2==0)+2*int(j>=2)+int(j>=4)) if n==1 else (4*int(j%2==0)+1+4*int(j>=2)+int(j>=4))
            check(f'new_helicity_channel_count_{n}_{j}',len(pairs)-expected)
            data[str(n)][str(j)]=pairs
    count_residuals=[]
    for n_global,wanted in [(0,[1,3,5]),(1,[0,2,4])]:
        local_degree=3-n_global
        allowed=sorted({2+local_degree-2*n_annihilate for n_annihilate in range(min(2,local_degree)+1)})
        count_residuals.extend(a-b for a,b in zip(allowed,wanted))
    check('new_cubic_local_particle_numbers_with_global_variables',count_residuals)
    aa,bb,g,k,V=sp.symbols('x_i x_f g k V',positive=True)
    R=sp.log(bb/aa)/3+(1/bb**2-1/aa**2)/2
    J=1/bb-1/aa
    amp=sp.I*sp.sqrt(6)*g*sp.sqrt(k)/sp.sqrt(V)*(R+sp.I*J)
    check('new_one_particle_leakage_probability',sp.expand(amp*sp.conjugate(amp))-6*g*g*k/V*(R*R+J*J))
    check('new_leakage_zero_interval',(R*R+J*J).subs(bb,aa))
    return data


# V52: article eq:v34-crossing-series-H (428); article eq:v34-crossing-series-L (429); article eq:v34r-Bianchi-monotonicity (437).
def test_nonlinear():
    print("Running test_nonlinear", flush=True)
    L,chi=sp.symbols('L chi',positive=True);rho,U,Up,Gp=sp.symbols('rho U Up Gp')
    # D_i Phi=0. The independent spatial derivative is D_i L.
    Dr=chi/L;Dchi=Gp-4*chi/L;Delta=L*Gp-chi
    check('new_spatial_Phi_constraint',
       (L**4*(Up-Dr)+4*L**3*(U-rho)).subs(rho,U+(L*Up-chi)/4))
    check('new_spatial_chi_constraint',Gp-4*Dr-Dchi)
    check('new_spatial_pressure_constraint',(Dr-Gp)/3+Delta/(3*L))
    acc=sp.Symbol('a_i');DL=-L*chi*acc/Delta
    check('new_nonlinear_Euler_equation',-chi*acc/3-Delta*DL/(3*L))
    # Exact FLRW vector field; Taylor coefficients obtained by differentiating it.
    HH,LL,CC,M=sp.symbols('H L chi M');G1,G2=sp.symbols('G1 G2')
    JJ=HH*LL-1;fc=G1*JJ-4*HH*CC;fh=CC/(6*M*M);fl=JJ
    def flow(expr):
        return sp.diff(expr,HH)*fh+sp.diff(expr,LL)*fl+sp.diff(expr,CC)*fc+sp.diff(expr,G1)*G2*fl
    A=G1*JJ;B=G2*JJ**2-3*HH*G1*JJ
    check('new_transverse_chi_second_derivative',flow(fc).subs(CC,0)-B)
    check('new_transverse_H_second_derivative',flow(fh).subs(CC,0)-A/(6*M*M))
    check('new_transverse_H_third_derivative',flow(flow(fh)).subs(CC,0)-B/(6*M*M))
    check('new_transverse_L_second_derivative',flow(fl).subs(CC,0)-HH*JJ)
    check('new_transverse_L_third_derivative',flow(flow(fl)).subs(CC,0)-(LL*A/(6*M*M)+HH*HH*JJ))
    tau=sp.Symbol('tau',real=True);AA=sp.Function('A')(tau);BB=sp.Function('B')(tau);Li=sp.Symbol('L_i')
    sol=AA*(Li-BB);theta=sp.Symbol('theta')
    check('new_exact_Einstein_branch_transport',sp.diff(sol,tau).subs({sp.diff(AA,tau):theta*AA/3,sp.diff(BB,tau):1/AA})-theta*sol/3+1)
    return {'D_rho':str(Dr),'D_chi':str(Dchi),'D_L':str(DL),'Taylor_A':str(A),'Taylor_B':str(B)}


# V52: article eq:v34r-simple-wave-equation (431); article eq:v34r-simple-wave-characteristics (432).
def test_simple_wave():
    print("Running test_simple_wave", flush=True)
    # Conformal invariance of sqrt(-g)*X_psi^2 in four dimensions.
    a,X=sp.symbols('a X',positive=True)
    check('new_radiation_scalar_conformal_identity',a**4*sp.Rational(4,3)*(X/a**2)**2-sp.Rational(4,3)*X**2)
    r,Q=sp.symbols('r Q',real=True,positive=True);c=1/sp.sqrt(3)
    # r=e^y avoids relying on assumptions in hyperbolic simplification.
    ch=(r+1/r)/2;sh=(r-1/r)/2
    q=Q*r**c;v=sh/ch;speed=(v+c)/(1+c*v)
    d=lambda f:sp.simplify(r*sp.diff(f,r))
    tau=-q*ch;kappa=q*sh
    check('new_simple_wave_integrability',-speed*d(kappa)-d(tau))
    check('new_simple_wave_field_equation',-speed*d(q*q*tau)-d(q*q*kappa))
    check('new_characteristic_speed_derivative',d(speed)-(1-speed*speed))
    aa,psi=sp.symbols('a psi',positive=True)
    Lrestore=aa*psi/q;chirestore=-4*q**4/aa**4
    check('new_simple_wave_original_fields',[
        (-q**4/aa**4)*Lrestore**4+psi**4,
        -Lrestore**3*chirestore-4*psi**3*q/aa])
    check('new_simple_wave_original_transport',
          -speed*(c*ch+sh/3)+(c*sh+ch/3))
    dt,l1,y1=sp.symbols('dt lambda_prime y_prime')
    jac=1+dt*l1*y1
    check('new_characteristic_jacobian',sp.diff(sp.Symbol('xi')+sp.Symbol('lambda')*dt,dt)-sp.Symbol('lambda'))
    check('new_gradient_catastrophe_denominator',jac.subs(dt,-1/(l1*y1)))
    return {'rapidity_parameter':'r=exp(y)','speed':str(sp.factor(speed)),'scalar_norm':str(q),'Jacobian':str(jac)}


def main_new():
    parser=argparse.ArgumentParser(description='Exact symbolic continuation of V31. No article or JSON input files are required.')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--recompute-fixed-81',action='store_true',help='Additionally recompute all 81 entries of the previous fixed-momentum table from the action.')
    parser.add_argument('--audit-fixed-reference',type=Path,help='Additionally audit the supplied previous 81-entry JSON table.')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    start=time.time()
    test_radiation_expansion();test_general_Legendre_transform();test_complete_example();test_time_integrals()
    previous_count=len(CHECKS)
    results={'scope':'State-parameterized finite-time radiation four-scalar angular kernel; evaluated scalar angular projections and one full radial contact contribution; nonlinear necessary conditions and exact special families. NOT the general IR-finite physical multichannel S matrix.'}
    test_signed_Legendre_transform()
    results['homogeneous']=test_homogeneous_and_state()
    results['all_angle_scalar']=test_generic_scalar_contact()
    results['direct_ADM_points']=test_direct_com()
    angular,coeff=test_angular_projection();results['contact_angular_moments']=angular
    results['radial_transverse_moments']=test_radial_kernel()
    results['channel_indices']=test_channel_count_and_leakage()
    results['nonlinear']=test_nonlinear();results['simple_wave']=test_simple_wave()
    new_check_count=len(CHECKS)-previous_count
    reference_before=len(CHECKS)
    if args.recompute_fixed_81:
        reference_path=args.output/'reference_fixed_momentum_81.json'
        compute_all(reference_path)
        audit_table(reference_path)
    if args.audit_fixed_reference:
        audit_table(args.audit_fixed_reference)
    reference_check_count=len(CHECKS)-reference_before
    (args.output/'symbolic_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (args.output/'angular_coefficients.json').write_text(json.dumps({'definition':'Polynomial coefficients d^A_{jrs,m} in the subsection; the infinite time/angular series is explicit. These are NOT all mixed scalar-tensor channel amplitudes.','j_values':list(range(9)),'data':coeff},ensure_ascii=False,indent=2),encoding='utf-8')
    payload={'python':platform.python_version(),'sympy':sp.__version__,'checks_passed':len(CHECKS),'previous_regression_checks':previous_count,'new_checks':new_check_count,'additional_fixed_table_checks':reference_check_count,'checks_failed':0,'duration_seconds':round(time.time()-start,3),'all_residuals_exact':True,'scope':results['scope'],'checks':CHECKS}
    (args.output/'verification_report.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'PASS {len(CHECKS)} exact checks: {previous_count} prior regressions, {new_check_count} new checks, {reference_check_count} optional fixed-table checks.',flush=True)
    print('General IR-finite physical multichannel S matrix: not computed or claimed.',flush=True)

if __name__=='__main__':
    main_new()
