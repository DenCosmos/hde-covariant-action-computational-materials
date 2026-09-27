#!/usr/bin/env python3
"""Check source connectivity and printed numerical tables against rerun CSV files.

This validates references, input files and numerical transcription; it does not
establish all analytical claims in the manuscript. No network is required.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article D.4, app:numerical-protocols (D.4).
  article D.4, tab:crossing-reproduction-v40 (5).
  article D.4, tab:late-rate-convergence-v40 (6).
  article 7.4, tab:v31-numerical-dynamics (3).
Method: Resolve references; parse 18 table values; compare to displayed-decimal tolerances.
Inputs: Explicit --article-source main.tex or folder, references.bib, three figures; supplementary source; newly computed numerical CSVs.
Outputs below the selected results root: release_checks.json.
Provenance: GitHub_bundle/scripts/verify_release.py.
Scope limit: No article redistributed. Saddle maxima are expected test constants; other 18 values parsed. No .aux needed.
"""
from __future__ import annotations
import argparse, collections, csv, json, math, re
from pathlib import Path


# V52: article app:numerical-protocols (D.4); article tab:crossing-reproduction-v40 (5); article tab:late-rate-convergence-v40 (6); article tab:v31-numerical-dynamics (3).
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--article-source', '--source', dest='source', type=Path, required=True,
                    help='V52 main.tex or its directory; references.bib and the three figure files must be alongside it. No .aux or LaTeX installation is needed.')
    ap.add_argument('--supplement-source', type=Path,
                    default=Path(__file__).resolve().parents[1]/'supplement/source')
    ap.add_argument('--report', type=Path,
                    default=Path(__file__).resolve().parents[1]/'results/current/release_checks.json')
    ap.add_argument('--results', type=Path, default=Path(__file__).resolve().parents[1]/'results/current')
    args = ap.parse_args()
    if args.source.is_file():
        if args.source.name != 'main.tex':
            ap.error('--article-source must name the V52 main.tex file or its directory')
        args.source = args.source.parent
    for required in ('main.tex', 'references.bib'):
        if not (args.source / required).is_file():
            ap.error(f'Missing author-supplied article input: {required}')
    tests = []
    def check(condition: bool, description: str) -> None:
        tests.append({'test': description, 'passed': bool(condition)})
        if not condition:
            raise AssertionError(description)
    texts = {'main.tex': (args.source/'main.tex').read_text(encoding='utf-8')}
    catalogue = (args.supplement_source/'mathematical_catalogue.tex').read_text(encoding='utf-8')
    clean = '\n'.join(re.sub(r'(?<!\\)%[^\n]*', '', t) for t in texts.values())
    labels = re.findall(r'\\label\{([^}]+)\}', clean)
    refs = re.findall(r'\\(?:eqref|ref|pageref|autoref)\*?\{([^}]+)\}', clean)
    check(len(labels)==len(set(labels)), 'All LaTeX labels are unique')
    check(set(refs)<=set(labels), 'All LaTeX cross-references resolve')
    bib = (args.source/'references.bib').read_text()
    keys = re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', bib)
    cited = set(k.strip() for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', clean) for k in group.split(','))
    check(len(keys)==len(set(keys))==68, '68 unique bibliography records, including the V52 supplement')
    check(cited==set(keys), 'Every citation resolves; every bibliography entry is cited')
    inputs = re.findall(r'\\input\{([^}]+)\}', clean)
    check(all((args.source/(p if Path(p).suffix else p+'.tex')).is_file() for p in inputs), 'Every LaTeX input file is present')
    figures = re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', clean)
    check(len(figures)==3 and all((args.source/p).is_file() for p in figures), 'All three included figure files are present')
    check(sum(t.count(r'\documentclass') for t in texts.values())==1, 'A single top-level document')
    check('V40MathematicalCatalogue' not in clean+bib, 'No obsolete external catalogue citation')
    main_count = len(re.findall(r'\\begin\{equation\}', texts['main.tex']))
    cat_count = len(re.findall(r'\\begin\{(?:equation|align)\}', catalogue))
    check(main_count==477 and cat_count==1145, '477 article/appendix equations and 1145 catalogue equations')
    cat_clean = re.sub(r'(?<!\\)%[^\n]*', '', catalogue)
    cat_labels = re.findall(r'\\label\{([^}]+)\}', cat_clean)
    cat_refs = re.findall(r'\\(?:eqref|ref|pageref|autoref)\*?\{([^}]+)\}', cat_clean)
    check(len(cat_labels) == len(set(cat_labels)), 'All supplement labels are unique')
    check(set(cat_refs) <= set(cat_labels), 'Supplement references resolve without the article')
    cat_bib = (args.supplement_source/'references.bib').read_text(encoding='utf-8')
    cat_keys = re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', cat_bib)
    cat_cited = set(k.strip() for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', cat_clean) for k in group.split(','))
    check(len(cat_keys) == len(set(cat_keys)) == 1 and cat_cited == set(cat_keys),
          'Supplement bibliography contains the cited V52 main article only')
    check(not re.search(r'\\(?:input|include|externaldocument)\s*\{', cat_clean),
          'Supplement has no external source or auxiliary-file inputs')
    appendix = texts['main.tex']
    def rows(name: str) -> list[dict]:
        with (args.results/'n02_background_diagnostics'/name).open(newline='') as f:
            return list(csv.DictReader(f))
    late = rows('late_rates.csv')
    printed_late = re.findall(r'\$(10|12)\$\s*&\s*\$\[(\d+),(\d+)\]\$\s*&\s*\$(-?[\d.]+)\$\s*&', appendix)
    check(len(printed_late)==6, 'Six late-time numerical rows located in the manuscript')
    comparisons=[]
    for n, lo, hi, printed in printed_late:
        candidates=[r for r in late if int(r['n'])==int(n) and float(r['rtol'])==1e-12 and float(r['log_a_start'])==int(lo) and float(r['log_a_end'])==int(hi)]
        check(len(candidates)==1, f'Unique late-time CSV match n={n}, interval=[{lo},{hi}]')
        actual=float(candidates[0]['measured']); tol=0.51*10**(-len(printed.split('.')[1]))
        check(abs(float(printed)-actual)<=tol, f'Printed late-time value agrees to all displayed decimals: {printed}')
        comparisons.append({'kind':'late','printed':printed,'rerun':actual,'tolerance':tol})
    crossing=rows('crossing_convergence.csv')
    printed_cross=re.findall(r'\$10\^\{-(\d+)\}\$\s*&\s*([12])\s*&\s*\$(-?[\d.]+)\$\s*&\s*\$(-?[\d.]+)\$', appendix)
    check(len(printed_cross)==6, 'Six crossing numerical rows located in the manuscript')
    for power, mode, printed, flux in printed_cross:
        candidates=[r for r in crossing if float(r['rtol'])==1e-12 and int(r['mode'])==int(mode) and math.isclose(float(r['minus_chi']),10**(-int(power)), rel_tol=1e-10)]
        check(len(candidates)==1, f'Unique crossing CSV match mode={mode}, power={power}')
        for text, col in ((printed,'d_ell_d_log_minus_chi'), (flux,'flux_log_coefficient')):
            actual=float(candidates[0][col]); tol=0.51*10**(-len(text.split('.')[1]))
            check(abs(float(text)-actual)<=tol, f'Printed crossing value agrees to all displayed decimals: {text}')
            comparisons.append({'kind':col,'printed':text,'rerun':actual,'tolerance':tol})
    saddle=rows('saddle_rates.csv')
    check(len(saddle)==18, 'All 18 saddle parameter combinations reproduced')
    smallest=[r for r in saddle if float(r['amplitude'])==1e-6]
    err=max(float(r['error']) for r in smallest); change=max(float(r['rtol_difference']) for r in smallest)
    check(abs(err-6.5163e-9)<=0.51e-13, 'Printed maximum saddle error matches the rerun')
    check(abs(change-3.2515e-10)<=0.51e-14, 'Printed saddle tolerance-change estimate matches the rerun')
    out={'passed':True,'checks':len(tests),'labels':len(labels),'reference_uses':len(refs),'bibliography_records':len(keys),'figures':figures,'article_equations':main_count,'catalogue_equations':cat_count,'tests':tests,'numerical_comparisons':comparisons,'numerical_input_directory':str(args.results),'scope':'Source connectivity and transcription of numerical tables, not a complete independent proof.'}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(f'PASS: {len(tests)} checks; {len(labels)} labels; {len(keys)} references; {len(comparisons)} printed numerical values.')

if __name__=='__main__':
    main()
