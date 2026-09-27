"""Exact sparse coefficients in Q(r,t)[i,sqrt(3)].

Every polynomial gcd is computed over Q, not over a primitive algebraic
extension. This removes a severe normalization cost in the final contraction.
The four basis elements are 1, i, sqrt(3), i*sqrt(3).
"""
from __future__ import annotations
import sympy as sp
from sympy.polys.fields import field
F,fr,ft=field('r,t',sp.QQ)
RADIUS,HALF_ANGLE=sp.symbols('r t',positive=True)
class Coeff:
    __slots__=('d',)
    def __init__(self,d=None):self.d={k:v for k,v in (d or {}).items() if v}
    def __bool__(self):return bool(self.d)
    def __add__(self,other):
        b=alg(other);d=self.d.copy()
        for k,v in b.d.items():
            vv=d.get(k,F.zero)+v
            if vv:d[k]=vv
            else:d.pop(k,None)
        return Coeff(d)
    __radd__=__add__
    def __neg__(self):return Coeff({k:-v for k,v in self.d.items()})
    def __sub__(self,b):return self+-alg(b)
    def __rsub__(self,b):return alg(b)+-self
    def __mul__(self,other):
        b=alg(other);d={}
        for (i,j),v in self.d.items():
            for (k,l),w in b.d.items():
                phase=(i+k)%2,(j+l)%2
                factor=(-1)**((i+k)//2)*3**((j+l)//2)
                d[phase]=d.get(phase,F.zero)+factor*v*w
        return Coeff(d)
    __rmul__=__mul__
    def inverse(self):
        if len(self.d)!=1:raise ValueError('Denominator must be a single radical phase')
        (i,j),v=next(iter(self.d.items()))
        return Coeff({(i,j):F.one/(v*((-1)**i)*3**j)})
    def __truediv__(self,b):return self*alg(b).inverse()
    def __rtruediv__(self,b):return alg(b)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        out=O
        for _ in range(n):out=out*self
        return out
    def expr(self,factor=True):
        out=sp.S.Zero
        for (i,j),v in self.d.items():
            lc=v.denom.LC;n=v.numer.quo_ground(lc).as_expr();d=v.denom.quo_ground(lc).as_expr()
            if factor:n=sp.factor(n);d=sp.factor(d)
            out+=sp.I**i*sp.sqrt(3)**j*n/d
        return out.subs({sp.Symbol('r'):RADIUS,sp.Symbol('t'):HALF_ANGLE})

def alg(x):
    if isinstance(x,Coeff):return x
    if x==0:return Z
    if x==sp.I:return I
    if x==sp.sqrt(3):return SQ3
    if isinstance(x,(int,sp.Integer,sp.Rational)):return Coeff({(0,0):F(sp.Rational(x))})
    ex=sp.expand(sp.cancel(x));d={}
    # Used only by small independent control formulas, never by table parsing.
    for term in sp.Add.make_args(ex):
        imag=1 if term.has(sp.I) else 0
        root=1 if term.has(sp.sqrt(3)) else 0
        rat=sp.cancel(term/(sp.I**imag*sp.sqrt(3)**root))
        d[imag,root]=d.get((imag,root),F.zero)+F.from_expr(rat.subs({RADIUS:sp.Symbol('r'),HALF_ANGLE:sp.Symbol('t')}))
    return Coeff(d)
Z=Coeff();O=Coeff({(0,0):F.one});I=Coeff({(1,0):F.one});SQ3=Coeff({(0,1):F.one})
class Domain:
    gens=(Coeff({(0,0):fr}),Coeff({(0,0):ft}))
    @staticmethod
    def convert(x):return alg(x)
K=Domain()
def la_add(a,b):
    d=a.copy()
    for k,v in b.items():
        vv=d.get(k,Z)+v
        if vv:d[k]=vv
        else:d.pop(k,None)
    return d
def la_scale(a,c,power=0):
    c=alg(c)
    return {k+power:v*c for k,v in a.items() if v*c}
def la_mul(a,b):
    d={}
    for p,v in a.items():
        for q,w in b.items():d[p+q]=d.get(p+q,Z)+v*w
    return {p:v for p,v in d.items() if v}
