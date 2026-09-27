#!/usr/bin/env python3
"""Local, explicitly limited CITATION.cff validation; no network.

This checks YAML syntax, required CFF 1.2 fields, release consistency and
identifier syntax/check digits. It is NOT the official CFF JSON schema or
cffconvert. Install cffconvert separately and run `cffconvert --validate`
for the official recommended validation before publication.
"""
from __future__ import annotations
import argparse, datetime, json, re, sys
from pathlib import Path
import yaml
ROOT = Path(__file__).resolve().parents[1]

def validate(path: Path) -> dict:
    tests = []
    def check(name: str, condition: bool) -> None:
        tests.append({'id': name, 'status': 'PASS' if condition else 'FAIL'})
    data = yaml.safe_load(path.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError('CITATION.cff must be a YAML mapping')
    for field in ('cff-version', 'message', 'title', 'authors'):
        check('required.' + field, bool(data.get(field)))
    check('cff_version', data.get('cff-version') == '1.2.0')
    check('software_type', data.get('type') == 'software')
    check('release_version', data.get('version') == (ROOT/'VERSION').read_text().strip())
    for i, author in enumerate(data.get('authors', [])):
        check(f'author.{i}.name', bool(author.get('family-names') and author.get('given-names')))
        if 'orcid' in author:
            value = str(author['orcid'])
            valid = bool(re.fullmatch(r'https://orcid.org/\d{4}-\d{4}-\d{4}-\d{3}[\dX]', value))
            if valid:
                digits = value.rsplit('/', 1)[1].replace('-', '')
                total = 0
                for digit in digits[:15]:
                    total = (total + int(digit)) * 2
                remainder = (12 - total % 11) % 11
                valid = digits[-1] == ('X' if remainder == 10 else str(remainder))
            check(f'author.{i}.orcid', valid)
    for field in ('repository-code', 'url'):
        if field in data:
            check(field + '.syntax', bool(re.fullmatch(r'https://[^\s<>]+', str(data[field]))))
    if 'doi' in data:
        check('doi.syntax', bool(re.fullmatch(r'10\.\d{4,9}/[^\s<>]+', str(data['doi']))))
    if 'date-released' in data:
        value = data['date-released']
        try:
            datetime.date.fromisoformat(str(value)); valid = True
        except ValueError:
            valid = False
        check('release_date.syntax', valid)
    if 'license' in data:
        check('license.nonempty', bool(data['license']))
    return {'status': 'PASS' if all(t['status']=='PASS' for t in tests) else 'FAIL',
            'level': 'Local YAML/required-field/identifier/release consistency checks only; NOT official CFF schema validation',
            'official_cffconvert_status': 'NOT RUN', 'tests': tests}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--file', type=Path, default=ROOT/'CITATION.cff')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = validate(args.file)
    except Exception as exc:
        result = {'status': 'FAIL', 'reason': str(exc), 'official_cffconvert_status': 'NOT RUN'}
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text)
    return 0 if result['status'] == 'PASS' else 1
if __name__ == '__main__':
    raise SystemExit(main())
