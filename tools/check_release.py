#!/usr/bin/env python3
"""Static local release preflight; not a substitute for run_all.py --mode full.

Checks checksums, Python syntax, local Markdown links, forbidden files and
candidate metadata. --for-publication adds author decision gates; no remote
access or official CFF schema validation is performed.
"""
from __future__ import annotations
import argparse, ast, json, re, subprocess, sys
from pathlib import Path
from release_common import ROOT, payload_files, check_checksums, write_checksums

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--update-checksums', action='store_true')
    parser.add_argument('--for-publication', action='store_true')
    parser.add_argument('--output', type=Path, help='Use an external path to avoid changing the checked payload')
    args = parser.parse_args()
    if args.update_checksums:
        write_checksums()
    errors = check_checksums()
    files = payload_files()
    checked_links, python_files = 0, 0
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        if path.suffix.lower() in {'.ttf','.otf','.woff','.woff2','.pfx','.pem','.key'}:
            errors.append('Forbidden distributable type: '+relative)
        if relative.startswith('article/') or path.name in {'V52_complete_package.zip','main.tex'}:
            errors.append('Main manuscript/private input unexpectedly included: '+relative)
        if path.suffix == '.py':
            try:
                ast.parse(path.read_text()); python_files += 1
            except SyntaxError as exc:
                errors.append(f'Python syntax: {relative}:{exc.lineno}')
        if path.suffix == '.md':
            text = path.read_text()
            for target in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\s]+)\)', text):
                if target.startswith(('https://','http://','mailto:','#')):
                    continue
                target = target.split('#',1)[0]
                if target:
                    checked_links += 1
                    if not (path.parent/target).exists():
                        errors.append('Broken local link: '+relative+' -> '+target)
    required = ['LICENSE','LICENSE_SCOPE.md','README.md','REPRODUCIBILITY.md','KNOWN_LIMITATIONS.md','SCRIPT_REFERENCE_MAP.md',
                'PUBLICATION_GUIDE_RU.md','CITATION.cff','CITATION.bib','requirements.txt',
                'publication/REFERENCE_MAP.json','publication/REFERENCE_MAP.csv',
                'validation/fresh_numbering.json','reference_results/v52_r1/run_summary.json']
    errors += ['Required file missing: '+x for x in required if not (ROOT/x).is_file()]
    frozen = ROOT/'reference_results/v52_r1/run_summary.json'
    if frozen.is_file():
        run = json.loads(frozen.read_text())
        if run.get('mode') != 'full' or run.get('status') != 'PASS':
            errors.append('Bundled actual full run is not PASS')
    code = 1 if errors else 0
    if args.for_publication:
        result = subprocess.run([sys.executable,str(ROOT/'tools/finalize_metadata.py')],check=False)
        if result.returncode:
            code = code or 2
    report = {'status': 'FAIL' if errors else ('BLOCKED' if code else 'PASS'),
              'scope':'Static release preflight, not a fresh scientific calculation or official CFF schema validation',
              'files':len(files),'python_files':python_files,'local_links':checked_links,
              'publication_author_gates_checked':args.for_publication,'errors':errors}
    text = json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(text)
    print(text)
    return code
if __name__ == '__main__':
    raise SystemExit(main())
