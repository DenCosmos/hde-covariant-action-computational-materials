#!/usr/bin/env python3
"""Copy real heavy checkpoints to a NEW output directory for explicit resume.

This is only staging, NOT scientific verification. The following full runner
must execute with --resume. It rechecks each generator's input fingerprint,
reassembles contacts and reruns all checks; derived homogeneous exchange is
cleared from staged angular records so it is recalculated, not only counted; it never treats success=true as
verification. Preserve the fresh-source computation logs separately.
"""
from __future__ import annotations
import argparse, importlib.util, json, shutil
from pathlib import Path
from release_common import ROOT,sha256

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-results',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();source=args.from_results.resolve();target=args.output.resolve()
    if not source.is_dir():
        parser.error('Source results directory does not exist')
    if target==ROOT or ROOT/'reference_results' in target.parents or target==ROOT/'reference_results':
        parser.error('Refusing to overwrite a protected reference directory')
    if target.exists() and any(target.iterdir()):
        parser.error('The target must be empty')
    target.mkdir(parents=True,exist_ok=True)
    records=[]
    for name in ('radial_generic','angular_v40','noncom'):
        folder=source/name
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob('*')):
            if not path.is_file() or path.suffix not in {'.json','.csv','.txt'}:
                continue
            relative=path.relative_to(source);out=target/relative
            out.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,out)
            source_hash=sha256(path)
            cleared=[]
            if name=='angular_v40' and path.name.startswith('j') and path.suffix=='.json':
                obj=json.loads(out.read_text())
                for key in ('global_exchange','global_dictionary'):
                    if key in obj:
                        obj.pop(key);cleared.append(key)
                if cleared:
                    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
            records.append({'file':relative.as_posix(),'source_sha256':source_hash,'staged_sha256':sha256(out),
                            'cleared_derived_fields_for_recomputation':cleared})
    if not records:
        parser.error('No heavy checkpoint files found; nothing staged')
    spec=importlib.util.spec_from_file_location('release_runner',ROOT/'run_all.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    fingerprint=module.fingerprint()
    context={'release':'V52-R1','scientific_version':'V52','input_fingerprint':fingerprint,
             'mode':'full','resume':True,'seeded_checkpoint_copy':True,
             'scientific_verification_status':'NOT RUN; invoke run_all.py --mode full --resume'}
    (target/'run_context.json').write_text(json.dumps(context,indent=2)+'\n')
    report={'operation':'STAGED_REAL_CHECKPOINTS','verification_status':'NOT RUN',
            'source_label':source.name,'copied_files':records,
            'instruction':'Only subsequent actually executed full-run tests establish PASS. Per-generator fingerprints are checked there.'}
    (target/'seed_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('STAGED',len(records),'files. Scientific checks NOT RUN. Next: run_all.py --mode full --resume --output',target)
    return 0
if __name__=='__main__':
    raise SystemExit(main())
