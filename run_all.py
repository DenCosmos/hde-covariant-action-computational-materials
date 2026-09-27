#!/usr/bin/env python3
"""Local V52-R1 computation runner. No network, remote writes or LaTeX required.

Quick executes the documented short identity, numerical and figure checks;
full also generates all 81 non-COM kernels, 81 radial bare contacts, 54 sources,
243 full-contact permutation identities and 74 mixed angular dictionaries.
The separately supplied V52 article is mandatory for full transcription tests.
PASS means a command actually returned zero, not a saved success field.
Resume validates algorithm/input fingerprints before using heavy checkpoints.
Use a NEW writable output folder; reference_results is never overwritten.
"""
from __future__ import annotations
import argparse,ast,datetime,hashlib,importlib.metadata,json,os,platform,subprocess,sys,time,traceback,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def semantic_bytes(path:Path)->bytes:
    if path.suffix!='.py':return path.read_bytes()
    tree=ast.parse(path.read_text())
    for node in ast.walk(tree):
        if hasattr(node,'body') and isinstance(node.body,list):
            node.body[:]=[x for x in node.body if not (isinstance(x,ast.Expr) and isinstance(x.value,ast.Constant) and isinstance(x.value.value,str))]
    return ast.dump(tree,include_attributes=False).encode()


def fingerprint()->str:
    digest=hashlib.sha256()
    for folder in ['scripts','lib','parameters','supplement/source']:
        for p in sorted((ROOT/folder).rglob('*')):
            if p.is_file() and p.suffix not in ['.pyc','.aux','.log','.out','.toc']:
                digest.update(p.relative_to(ROOT).as_posix().encode());digest.update(semantic_bytes(p))
    for dep in ['sympy','numpy','scipy','matplotlib']:
        digest.update((dep+'='+importlib.metadata.version(dep)).encode())
    return digest.hexdigest()


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=['quick','full'],default='quick')
    ap.add_argument('--output',type=Path,default=ROOT/'results/quick')
    ap.add_argument('--article-source',type=Path,help='Author-supplied V52 main.tex or its directory, with references.bib and figures')
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--resume',action='store_true',help='Resume fingerprinted heavy checkpoints; short checks still rerun')
    args=ap.parse_args()
    if not 1<=args.workers<=16:ap.error('--workers must be between 1 and 16; 4 was tested')
    out=args.output.resolve()
    for protected in ['reference_results','supplement','data','publication','provenance','scripts','lib']:
        if out==ROOT/protected or ROOT/protected in out.parents:ap.error('Refusing a protected output directory')
    if out==ROOT:ap.error('The repository root cannot be the output directory')
    fp=fingerprint();manifest=out/'run_context.json'
    if out.exists() and any(out.iterdir()):
        if not args.resume:ap.error('Output is not empty. Choose a new directory, or --resume for this same computation.')
        if not manifest.is_file() or json.loads(manifest.read_text()).get('input_fingerprint')!=fp:
            ap.error('No matching runner fingerprint. Existing files cannot be treated as compatible checkpoints.')
    out.mkdir(parents=True,exist_ok=True);logs=out/'logs';logs.mkdir(exist_ok=True)
    manifest.write_text(json.dumps({'release':'V52-R1','scientific_version':'V52','input_fingerprint':fp,'mode':args.mode,'resume':args.resume},indent=2)+'\n')
    env=os.environ.copy();env['V52_RESULTS_ROOT']=str(out);env['V52_OFFLINE']='1'
    env['PYTHONPATH']=str(ROOT/'tools/offline_guard')
    env['MPLCONFIGDIR']=str(out/'.matplotlib-cache')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    environment={'python':platform.python_version(),'platform':platform.platform(),'executable':'<active Python>','dependencies':{k:importlib.metadata.version(k) for k in ['sympy','mpmath','numpy','scipy','matplotlib','PyYAML']},'offline_guard':True,'workers':args.workers,'numeric_threads':1}
    (out/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    records=[];lock=threading.Lock();started=time.monotonic()
    def save():
        with lock:
            (out/'progress.json').write_text(json.dumps({'mode':args.mode,'elapsed_seconds':time.monotonic()-started,'records':records},indent=2)+'\n')
    def record(name,status,**kw):
        row={'stage':name,'status':status,**kw}
        with lock:records.append(row)
        save();return row
    def run(name,file,arguments=()):
        command=[sys.executable,'-u',str(ROOT/'scripts'/file),*map(str,arguments)]
        visible=['python','scripts/'+file,*[str(x).replace(str(out),'<OUTPUT>').replace(str(ROOT),'<REPO>') for x in arguments]]
        t=time.monotonic();print('START',name,flush=True)
        try:
            with (logs/(name+'.log')).open('w',encoding='utf-8') as stream:
                proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
                beat=t
                while proc.poll() is None:
                    time.sleep(1)
                    if time.monotonic()-beat>=60:
                        elapsed=time.monotonic()-t
                        (logs/(name+'.checkpoint.json')).write_text(json.dumps({'status':'RUNNING','elapsed_seconds':elapsed,'command':visible},indent=2)+'\n')
                        print('RUNNING',name,f'{elapsed:.0f}s',flush=True);beat=time.monotonic()
                code=proc.returncode
            status='PASS' if code==0 else ('BLOCKED' if code==2 else 'FAIL')
            row=record(name,status,command=visible,returncode=code,seconds=time.monotonic()-t,log='logs/'+name+'.log')
            (logs/(name+'.run.json')).write_text(json.dumps(row,indent=2)+'\n')
            print(status,name,f'{row["seconds"]:.1f}s',flush=True);return status=='PASS'
        except Exception:
            (logs/(name+'.error.txt')).write_text(traceback.format_exc())
            record(name,'FAIL',command=visible,returncode=1,seconds=time.monotonic()-t,reason='Runner exception; see error log')
            return False
    def block(name,reason):record(name,'BLOCKED',returncode=2,reason=reason)
    jobs=[
      ('analytic','article_appsDE_verify_reduced_identities.py',['--output',out/'analytic_v40']),
      ('local_vertices','supplement_sec1_6_verify_local_vertices.py',['--output',out/'local_vertices_v40']),
      ('entropy','article_sec7_3_verify_entropy.py',[]),
      ('printed_catalogue','supplement_sec1_verify_printed_catalogue.py',['--results',out]),
      ('background','article_appD4_compute_background_diagnostics.py',['--output',out]),
      ('figure01','article_fig01_plot_pole_crossing.py',['--output',out/'figure01']),
      ('figure02','article_fig02_plot_phase_portrait.py',['--output',out/'figure02']),
      ('figure03','article_fig03_plot_strong_coupling.py',['--output',out/'figures','--data-output',out/'figures']),
      ('mixed_normalization','article_appC1_verify_mixed_normalization.py',['--output',out/'mixed_normalization']),
      ('finite_time_identities','article_appC_verify_finite_time_and_nonlinear_identities.py',['--output',out/'engine_regressions'])]
    for job in jobs:run(*job)
    if args.article_source:
        if args.article_source.exists():run('source_transcription','article_V52_verify_source_transcription.py',['--article-source',args.article_source.resolve(),'--results',out,'--report',out/'release_checks.json'])
        else:block('source_transcription','The --article-source path does not exist')
    elif args.mode=='full':block('source_transcription','Full mode requires --article-source. Main article intentionally not redistributed.')
    else:record('source_transcription','NOT RUN',required_in_mode=False,reason='External source check; supply --article-source. Mandatory in full mode.')
    heavy_names=['noncom','radial_generation','contacts','permutations','radial_controls','printed_vs_recomputed','angular','global_exchange','angular_controls','soft','noncom_windows']
    if args.mode=='full':
        flags=['--resume'] if args.resume else []
        parallel=[('noncom','article_appC3_compute_noncom_matrix_elements.py',['--output',out/'noncom',*flags])]
        for task,n,step in [('pair',54,9),('bare',81,27)]:
            for i in range(0,n,step):parallel.append((f'radial_{task}_{i}_{min(i+step,n)}','supplement_sec1_3_generate_radial_sources.py',['--task',task,'--start',i,'--stop',min(i+step,n),*flags]))
        outcomes={}
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures={executor.submit(run,*j):j[0] for j in parallel}
            for f in as_completed(futures):outcomes[futures[f]]=f.result()
        radial_ok=all(ok for name,ok in outcomes.items() if name.startswith('radial_'))
        if radial_ok:
            contacts_ok=run('contacts','supplement_sec1_2_assemble_contacts.py')
            if contacts_ok:
                for name,file in [('permutations','article_appC3_verify_contact_permutations.py'),('radial_controls','article_appC3_verify_radial_controls.py'),('printed_vs_recomputed','supplement_secs1_2_1_5_compare_recomputed_coefficients.py')]:run(name,file)
                angular_ok=run('angular','article_appC5_compute_angular_coefficients.py',flags)
                if angular_ok:
                    glob_ok=run('global_exchange','article_appC7_complete_global_exchange.py')
                    if glob_ok:run('angular_controls','article_appC5_verify_angular_controls.py')
                    else:block('angular_controls','Global exchange step failed')
                else:
                    block('global_exchange','Angular generation failed');block('angular_controls','Angular generation failed')
                run('soft','article_appC6_verify_soft_limit.py')
            else:
                for name in ['permutations','radial_controls','printed_vs_recomputed','angular','global_exchange','angular_controls','soft']:block(name,'Contact assembly failed')
        else:
            for name in ['contacts','permutations','radial_controls','printed_vs_recomputed','angular','global_exchange','angular_controls','soft']:block(name,'Radial generation incomplete or failed')
        if outcomes['noncom']:
            run('noncom_windows','article_appC3_verify_noncom_finite_windows.py',['--input',out/'noncom/coefficients.json','--output',out/'noncom'])
        else:block('noncom_windows','Non-COM generation failed')
    else:
        for name in heavy_names:record(name,'NOT RUN',required_in_mode=False,reason='Full-only calculation; quick PASS does not cover this stage.')
    states=[r['status'] for r in records]
    code=1 if 'FAIL' in states else (2 if 'BLOCKED' in states else 0)
    final={'release':'V52-R1','scientific_version':'V52','mode':args.mode,'status':'PASS' if code==0 else ('FAIL' if code==1 else 'BLOCKED'),'returncode':code,'seconds':time.monotonic()-started,'input_fingerprint':fp,'records':records,'scope':'Full documented computational scope' if args.mode=='full' else 'Quick subset only; heavy computations are NOT RUN','no_regulator_free_unitarity_claim':True}
    (out/'run_summary.json').write_text(json.dumps(final,indent=2)+'\n')
    print(final['status'],args.mode,'return code',code,flush=True);return code
if __name__=='__main__':
    raise SystemExit(main())
