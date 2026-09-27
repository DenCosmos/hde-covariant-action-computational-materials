"""V52-R1 exact half-angle conversion used by supplement Sec. 1.5.
Pure to_z() extracted verbatim from build_mathematical_catalogue_v40.py.
No catalogue generation or manuscript editing is imported or executed.
Method t^2=(1-z)/(1+z); even/odd radical sectors remain separate.
"""
import radial_coeff_algebra as a
R,Z=a.F.gens
def to_z(v):
    u=(1-Z)/(1+Z)
    def pol(p):
        parity={m[1]%2 for m in p}
        if len(parity)>1:raise ValueError('Mixed parity')
        b=next(iter(parity),0)
        return sum((c*R**m[0]*u**((m[1]-b)//2) for m,c in p.items()),a.F.zero),b
    d={};odd={}
    for ph,rat in v.d.items():
        n,p=pol(rat.numer);de,b=pol(rat.denom)
        if p==b:d[ph]=n/de
        elif p-b==1:odd[ph]=n/de/(1+Z)
        else:odd[ph]=n/de/(1-Z)
    return a.Coeff(d),a.Coeff(odd)
