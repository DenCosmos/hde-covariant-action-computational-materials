#!/usr/bin/env python3
"""Compile a separate copy of the unchanged V52 supplement, never its generator.

Requires pdflatex and BibTeX locally. No network is used and the mathematical
source remains one main .tex with the provided bibliography. The scientific
Python checks do not require this optional editorial operation.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, subprocess, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/supplement_build')
    parser.add_argument('--latex', default='pdflatex')
    parser.add_argument('--bibtex', default='bibtex')
    args = parser.parse_args()
    output = args.output.resolve()
    if output == ROOT or ROOT/'supplement' in output.parents or output == ROOT/'supplement':
        parser.error('Refusing to overwrite the supplement or repository')
    if output.exists() and any(output.iterdir()):
        parser.error('Use a new empty build directory')
    output.mkdir(parents=True, exist_ok=True)
    source = ROOT/'supplement/source'
    copied = {}
    for path in source.iterdir():
        if path.is_file():
            shutil.copy2(path, output/path.name)
            copied[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    records = []
    commands = [[args.latex, '-interaction=nonstopmode', '-halt-on-error', 'mathematical_catalogue.tex'],
                [args.bibtex, 'mathematical_catalogue'],
                [args.latex, '-interaction=nonstopmode', '-halt-on-error', 'mathematical_catalogue.tex'],
                [args.latex, '-interaction=nonstopmode', '-halt-on-error', 'mathematical_catalogue.tex']]
    status, code = 'PASS', 0
    for index, command in enumerate(commands, 1):
        if not shutil.which(command[0]):
            status, code = 'BLOCKED', 2
            records.append({'command': command, 'status': status, 'reason': 'Executable not installed'})
            break
        start = time.monotonic()
        with (output/f'pass{index}.txt').open('w') as log:
            result = subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=False)
        records.append({'command': command, 'returncode': result.returncode, 'seconds': time.monotonic()-start})
        if result.returncode:
            status, code = 'FAIL', 1
            break
    if status == 'PASS':
        log = (output/'mathematical_catalogue.log').read_text(errors='replace')
        if 'There were undefined references' in log or 'Citation' in log and 'undefined' in log:
            status, code = 'FAIL', 1
    report = {'status': status, 'source_sha256': copied, 'commands': records,
              'scope': 'Independent LaTeX build of unchanged supplement; no historical editorial generators'}
    (output/'build_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(status, output)
    return code
if __name__ == '__main__':
    raise SystemExit(main())
