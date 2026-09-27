"""Sparse exact ADM expansion used by Appendices C.2–C.7.

The nilpotent leg variables select one factor of every external field. Laurent
coefficients in conformal time are rational functions of r and t=tan(theta/2)
over QQ(i,sqrt(3)). V40 uses unnormalized real polarization tensors and
contracts them with their explicit Gram norms.
This is the algebraic engine, not a stored table and not an S-matrix solver.
Its scope is the constant-potential radiation background. General finite-time
and state data remain explicit in the time-integral and global-sector modules.

Adapted from provenance/original/verify_scattering_extension.py. Test functions
are relocated to the matching sections; no archived script is imported.
"""

from __future__ import annotations


from dataclasses import dataclass


from functools import lru_cache


from itertools import permutations, product


import sympy as sp


RADIUS, HALF_ANGLE = sp.symbols('r t', positive=True)
K=sp.QQ.algebraic_field(sp.I,sp.sqrt(3)).frac_field(RADIUS, HALF_ANGLE)


Z=K.zero; O=K.one


Z=K.zero; O=K.one


def alg(x):
    if isinstance(x,type(O)):return x
    return K.from_sympy(sp.sympify(x))


I=alg(sp.I); SQ3=alg(sp.sqrt(3))


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


import ast


import argparse


import platform


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


def legendre_Q_explicit(j,u):
    return sp.legendre(j,u)*sp.log((u+1)/(u-1))/2-sum(sp.legendre(r-1,u)*sp.legendre(j-r,u)/r for r in range(1,j+1))


def radial_transverse_moment(j,u):
    Q=legendre_Q_explicit(j,u)
    m0=2*int(j==0);m1=sp.Rational(2,3)*int(j==1)
    m2=sp.Rational(2,3)*int(j==0)+sp.Rational(4,15)*int(j==2)
    return -2*((2*u-4*u**3)*Q+(u*u-u**4)*sp.diff(Q,u))-m2-2*u*m1-(3*u*u-1)*m0

