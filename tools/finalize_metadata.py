#!/usr/bin/env python3
"""Apply ONLY author-completed local metadata; no remote actions or licence choice.

Edit publication/metadata_pending.json first. --apply requires a real GitHub
URL, actual planned release date, licence descriptions and existing licence
texts plus three author confirmations. For publication_target="github", DOI
is optional. For "github-and-zenodo", the actual reserved version DOI is required.
The optional custom license_url must be an actual absolute HTTPS URL.
No commit hash is written into its own commit. DOI reservation is NOT publication.
"""
from __future__ import annotations
import argparse, datetime, json, re
from pathlib import Path
import yaml
ROOT = Path(__file__).resolve().parents[1]

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    path = ROOT/'publication/metadata_pending.json'
    meta = json.loads(path.read_text())
    missing = []
    for field in ('author_metadata_confirmed','public_access_confirmed','licenses_confirmed'):
        if meta.get(field) is not True:
            missing.append(field)
    if not re.fullmatch(r'https://github\.com/[^/\s<>]+/[^/\s<>]+/?', meta.get('repository_url') or ''):
        missing.append('repository_url')
    target = meta.get('publication_target', 'github-and-zenodo')
    if target not in ('github', 'github-and-zenodo'):
        missing.append('publication_target (github or github-and-zenodo)')
    doi = meta.get('version_doi')
    if (target != 'github' or doi) and not re.fullmatch(r'10\.\d{4,9}/[^\s<>]+', doi or ''):
        missing.append('version_doi')
    license_url = meta.get('license_url')
    if license_url and not re.fullmatch(r'https://[^\s<>]+', license_url):
        missing.append('license_url')
    try:
        datetime.date.fromisoformat(meta.get('actual_release_date') or '')
    except ValueError:
        missing.append('actual_release_date')
    for kind in ('code','text','data'):
        if not meta.get(kind+'_license'):
            missing.append(kind+'_license')
        value = meta.get(kind+'_license_file')
        candidate = (ROOT/value).resolve() if value else None
        if not candidate or ROOT not in candidate.parents or not candidate.is_file() or candidate.stat().st_size == 0:
            missing.append(kind+'_license_file')
    if missing:
        print('BLOCKED: author decisions or verified values missing:', ', '.join(missing))
        return 2
    print('PASS: local author-controlled fields present; identifiers were not resolved over the network.')
    if not args.apply:
        print('Dry run only. Use --apply after author review.')
        return 0
    cff_path = ROOT/'CITATION.cff'
    cff = yaml.safe_load(cff_path.read_text())
    message = ('Please cite computational version V52-R1 and the associated scientific publication. '
               'Copyright 2026 Danylo Yerokhin. All rights reserved; '
               'limited scholarly verification permission only. See LICENSE.')
    cff.update({'message': message, 'repository-code': meta['repository_url'],
                'date-released': meta['actual_release_date']})
    if doi:
        cff['doi'] = doi
    else:
        cff.pop('doi', None)
    if license_url:
        cff['license-url'] = license_url
    else:
        cff.pop('license-url', None)
    # CFF has one release-wide licence field, not file-scope maps. Avoid claiming
    # that distinct code/text/data licences are alternative grants for every file.
    cff.pop('license', None)
    cff_path.write_text(yaml.safe_dump(cff, sort_keys=False, allow_unicode=True, width=100))
    title = meta['title']
    doi_line = f'  doi = {{{doi}}},\n' if doi else ''
    bib = ('@software{YerokhinV52R1Computational,\n'
           '  author = {Yerokhin, Danylo},\n'
           f'  title = {{{title}}},\n'
           '  version = {V52-R1},\n'
           f'  year = {{{meta["actual_release_date"][:4]}}},\n'
           f'  date = {{{meta["actual_release_date"]}}},\n'
           f'{doi_line}'
           f'  url = {{{meta["repository_url"]}}},\n'
           '  note = {Computational release for scientific article V52; proposed tag v52-r1}\n}\n')
    (ROOT/'CITATION.bib').write_text(bib)
    # Preserve the author-approved scope and platform caveats; do not replace
    # the existing licence notice as a side effect of citation finalization.
    if not (ROOT/'LICENSE_SCOPE.md').exists():
        (ROOT/'LICENSE_SCOPE.md').write_text('# Author-approved licence scopes\n\n'+
            '\n'.join(f'* {kind}: {meta[kind+"_license"]}; [{meta[kind+"_license_file"]}]({meta[kind+"_license_file"]}).'
                        for kind in ('code','text','data'))+'\n\nThird-party rights remain applicable.\n')
    meta['status'] = 'AUTHOR_METADATA_PREPARED_NOT_PUBLISHED'
    path.write_text(json.dumps(meta, indent=2, ensure_ascii=False)+'\n')
    (ROOT/'publication/RELEASE_IDENTITY.md').write_text(
        '# Author-prepared release identity\n\nVersion: V52-R1. Article: V52. Tag: `v52-r1`.\n\n'
        f'Repository: {meta["repository_url"]}\n\nVersion DOI: {doi or "not assigned; GitHub-only route"}\n\n'
        f'Intended actual release date: {meta["actual_release_date"]}\n\n'
        'No publication has been performed by this helper. Correct the date and rerun before committing if it changes. '
        'Keep final commit SHA and ZIP SHA-256 in the external release record.\n')
    print('Written CITATION files and scoped author metadata. Revalidate and regenerate SHA256SUMS before committing.')
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
