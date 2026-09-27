#!/usr/bin/env python3
"""Build local ZIP candidates; --github-only omits Zenodo. Never upload/publish.

Refresh checksum manifest, run preflight, then write outside the repository.
Metadata changes require rebuilding before the final commit. After a clean
final commit, the same command is deterministic and must not change its files.
--allow-candidate is required while author-controlled publication gates remain
unfilled. ZIP SHA-256 values belong to an external release record, not the ZIP.
"""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys, zipfile
from pathlib import Path
from release_common import ROOT,payload_files,sha256,write_checksums

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--allow-candidate',action='store_true')
    parser.add_argument('--github-only',action='store_true',help='Build only the GitHub ZIP; do not create a Zenodo candidate')
    args=parser.parse_args(); destination=args.output_dir.resolve()
    if destination==ROOT or ROOT in destination.parents:
        parser.error('Archive output must be outside the repository')
    write_checksums()
    command=[sys.executable,str(ROOT/'tools/check_release.py')]
    if not args.allow_candidate:
        command.append('--for-publication')
    result=subprocess.run(command,check=False)
    if result.returncode:
        return result.returncode
    destination.mkdir(parents=True,exist_ok=True)
    github=destination/('V52_R1_GitHub_ready_All_Rights_Reserved.zip' if args.github_only else 'V52_R1_GitHub_ready.zip')
    with zipfile.ZipFile(github,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in payload_files():
            name='V52_R1_repository/'+path.relative_to(ROOT).as_posix()
            info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o100644<<16
            archive.writestr(info,path.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    outputs = [github]
    if not args.github_only:
        zenodo=destination/'V52_R1_Zenodo_deposit.zip'
        shutil.copyfile(github,zenodo)
        outputs.append(zenodo)
    report={'status':'LOCAL_CANDIDATE_NOT_PUBLISHED','release':'V52-R1','scientific_version':'V52',
            'github_only':args.github_only,
            'identical_zip_bytes':None if args.github_only else sha256(github)==sha256(zenodo),
            'files':[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)} for p in outputs],
            'final_commit_sha':None,'instruction':'Author records final verified commit SHA here, outside the repository. No remote actions performed.'}
    (destination/'release_record.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return 0
if __name__=='__main__':
    raise SystemExit(main())
