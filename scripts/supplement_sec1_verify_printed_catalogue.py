#!/usr/bin/env python3
"""Read the printed catalogue, not a historical coefficient cache.

Checks its 331 functions, 65 h=b+u+v identities, 81-row reconstruction map,
and agreement of 27 printed W derivatives / five vertices with an independently
regenerated local-vertex calculation. This is NOT a rederivation of the full
mixed gravitational exchange amplitude.

Run from the package root: python scripts/verify_catalogue.py

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  supplement 1.2, sec:cat-contacts (1.2).
  supplement 1.2.14, sec:cat-component-table (1.2.14).
  supplement 1.3, sec:cat-pairs (1.3).
  supplement 1.6, sec:cat-local-vertices (1.6).
Method: Parse 331 rational functions and reconstruction/vertex tables; exact identities in printed coefficients.
Inputs: supplement/source; newly calculated local_vertices_v40.
Outputs below the selected results root: catalogue_checks.json.
Provenance: GitHub_bundle/scripts/verify_catalogue.py.
Scope limit: Printed-coefficient checks, not independent ADM derivation; the separate strict comparator supplies that comparison.
"""
from __future__ import annotations
import argparse, itertools, json, re, time
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication

TRANSFORMS = standard_transformations + (implicit_multiplication,)
NAMES = 'r z L A f2 f3 f4 f5 s phi phi_dot grad_phi_squared'.split() + [f'W{i}{j}' for i in range(7) for j in range(7)]
LOCALS = {name: sp.Symbol(name) for name in NAMES}
LOCALS.update(i=sp.I, sqrt=sp.sqrt)
r,z = LOCALS['r'], LOCALS['z']

def group(text: str, pos: int) -> tuple[str,int]:
    while text[pos].isspace(): pos += 1
    if text[pos] != '{': raise ValueError('Expected braced argument: '+text[pos:pos+80])
    depth, start = 1, pos+1
    pos += 1
    while pos < len(text):
        depth += (text[pos] == '{') - (text[pos] == '}')
        if depth == 0: return text[start:pos], pos+1
        pos += 1
    raise ValueError('Unbalanced argument')

def plain_latex(text: str) -> str:
    # Flatten display-only wrappers; retain every algebraic term and sign.
    text = text.replace(r'\shortstack[l]', '').replace(r'\vcenter','').replace(r'\hbox','')
    text = text.replace(r'\displaystyle','').replace(r'\notag','').replace('&','').replace('$','').replace(r'\\',' ')
    text = text.replace(r'\left','').replace(r'\right','').replace(r'\,',' ').strip().rstrip('.')
    text = re.sub(r'\\mathcal W_\{(\d)(\d)\}',lambda m:'W'+m[1]+m[2],text)
    text = text.replace(r'A_{L}','A').replace(r'G_{\phi}','grad_phi_squared')
    text = text.replace(r'\phi','phi')
    text = re.sub(r'f_\{([2-5])\}',r'f\1',text)
    text = re.sub(r'\bu\b','phi_dot',text)
    for cmd in [r'\frac',r'\sqrt']:
        while cmd in text:
            a=text.index(cmd); arg,b=group(text,a+len(cmd))
            if cmd == r'\frac':
                arg2,b=group(text,b); new='(('+plain_latex(arg)+')/('+plain_latex(arg2)+'))'
            else: new='sqrt('+plain_latex(arg)+')'
            text=text[:a]+new+text[b:]
    text=text.replace('{}','').replace('{','(').replace('}',')').replace('^','**')
    text=' '.join(text.split())
    if '\\' in text: raise ValueError('Unsupported command: '+text)
    return text

def math(text: str) -> sp.Expr:
    return parse_expr(plain_latex(text),local_dict=LOCALS,transformations=TRANSFORMS)

# V52: supplement sec:cat-contacts (1.2); supplement sec:cat-component-table (1.2.14); supplement sec:cat-pairs (1.3); supplement sec:cat-local-vertices (1.6).
def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1]/'supplement/source')
    parser.add_argument('--results',type=Path,default=Path(__file__).resolve().parents[1]/'results/current')
    args=parser.parse_args(); args.results.mkdir(exist_ok=True,parents=True)
    text=(args.source/'mathematical_catalogue.tex').read_text()
    out=args.results/'catalogue_checks.json'; checks=[]; start=time.monotonic()
    def check(ok: bool, desc: str) -> None:
        checks.append({'test':desc,'passed':bool(ok)})
        if not ok: raise AssertionError(desc)
    equations=[m[2] for m in re.finditer(r'\\begin\{(equation|align)\}([\s\S]*?)\\end\{\1\}',text)]
    bylabel={m[1]:e.split(r'\label')[0].strip() for e in equations if (m:=re.search(r'\\label\{([^}]+)\}',e))}
    R={}; Ns={}; Ds={}
    for label,eq in bylabel.items():
        m=re.fullmatch(r'eq:cat-(rational|numerator|denominator)-(\d+)',label)
        if not m:continue
        kind,index=m[1],int(m[2]); rhs=eq.split('=',1)[1]
        if kind=='numerator': Ns[index]=math(rhs)
        elif kind=='denominator': Ds[index]=math(rhs)
        elif 'N_{' not in rhs:R[index]=math(rhs)
    for index in Ns:R[index]=Ns[index]/Ds[index]
    check(set(R)==set(range(1,332)), 'All 331 R functions are defined exactly once and parsed')
    referenced=set(map(int,re.findall(r'R_\{(\d+)\}',text)))
    check(referenced==set(R),'Every R index resolves to a definition')
    check(set(Ns)==set(Ds),'Numerators and denominators are paired')
    print(f'Parsed 331 algebraic functions in {time.monotonic()-start:.1f}s',flush=True)
    maps={(part,indices,int(power)):R[int(idx)] for part,indices,power,idx in re.findall(r'([hbuv])\^\{(\d{4})\}_\{(-?\d+)\}=R_\{(\d+)\}',text)}
    reps=sorted({key[1] for key in maps})
    check(len(reps)==13,'Thirteen contact representatives')
    for indices in reps:
        for power in range(-4,1):
            expr=maps.get(('h',indices,power),0)-sum(maps.get((part,indices,power),0) for part in 'buv')
            check(sp.cancel(sp.expand(expr))==0, f'h=b+u+v: {indices}, power {power}')
        print(f'Contact {indices}: exact identities passed',flush=True)
    rows=re.findall(r'^(\d{4}) & (\d{4}) & ([01]) & ([01]) & ([01])\\\\',text,re.M)
    check(len(rows)==81,'The reconstruction table has 81 rows')
    check(len({row[0] for row in rows})==81,'Unique reconstruction-table component indices')
    for comp,rep,sflag,tflag,nonzero in rows:
        a,b,c,d=map(int,comp); swaps=0
        if a>b:a,b=b,a;swaps+=1
        if c>d:c,d=d,c;swaps+=1
        pair_swap=int((a,b)>(c,d))
        if pair_swap:a,b,c,d=c,d,a,b
        expected=''.join(map(str,(a,b,c,d)))
        check(rep==expected and int(sflag)==swaps%2 and int(tflag)==pair_swap and int(nonzero)==(comp.count('2')%2==0),f'Reconstruction indexing {comp}')
    check(sum(int(row[4]) for row in rows)==41,'41 even-polarization components; 40 odd components vanish')
    # Local derivatives / vertices from a separately regenerated calculation.
    raw=json.loads((args.results/'local_vertices_v40/results.json').read_text())
    for key,expected in raw['derivatives'].items():
        i,j=key.split(','); rhs=bylabel[f'eq:cat-local-W-{i}-{j}'].split('=',1)[1]
        check(sp.cancel(math(rhs)-sp.sympify(expected,locals=LOCALS))==0,f'Printed local derivative W{i}{j}')
    for order,expected in raw['vertices'].items():
        rhs=bylabel[f'eq:cat-local-L-{order}'].split('=',1)[1]
        check(sp.expand(math(rhs)-sp.sympify(expected,locals=LOCALS))==0,f'Printed local vertex L{order}')
    check(len(re.findall(r'\\subsubsection\{Pair \d{2}:\d{2}\}',text))==54,'All 54 pair-species headings are present')
    summary={'passed':True,'checks':len(checks),'runtime_s':time.monotonic()-start,'sympy':sp.__version__,
             'scope':'Printed-coefficient algebra and indexing, NOT full rederivation of mixed gravitational amplitudes','tests':checks}
    out.write_text(json.dumps(summary,indent=2))
    print(f'PASS: {len(checks)} tests; {summary["runtime_s"]:.1f}s',flush=True)

if __name__=='__main__':main()
