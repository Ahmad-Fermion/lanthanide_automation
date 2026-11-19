#!/usr/bin/env python3
"""Conservative, stage-aware recovery for ordinary and Screwer TURBOMOLE runs.

This utility *never* invokes DEFINE or Screwer, never changes an electronic
configuration, and never changes calculation files at submission time. It only
submits dependency-linked SLURM jobs for stages justified by existing evidence.
Default is a read-only plan; --submit is explicit.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import sys

METALS = ('la', 'ce', 'pr', 'nd', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu')
SCF_DONE = re.compile(r'convergence\s+criteria\s+satisfied\s+after\s+\d+\s+iterations', re.I)
SCF_ENERGY = re.compile(r'ENERGY\s+CONVERGED\s*!', re.I)
SCF_END = re.compile(r'ridft\s*:\s*all\s+done', re.I)

class RecoveryError(RuntimeError):
    pass


def text(path):
    try:
        return path.read_text(errors='replace')
    except OSError:
        return ''


def scf_ok(path):
    txt = text(path)
    return bool(SCF_DONE.search(txt) and SCF_ENERGY.search(txt) and SCF_END.search(txt))


def geometry_ok(folder):
    return (not (folder / 'not.converged').exists() and
            not (folder / 'GEO_OPT_FAILED').exists() and
            ((folder / 'GEO_OPT_CONVERGED').is_file() or (folder / 'converged').is_file()) and
            (folder / 'gradient').is_file() and (folder / 'gradient').stat().st_size > 0)


def hessian_ok(folder, method):
    if method == 'none':
        return False
    hess = folder / 'hessian'
    spec = folder / 'vibspectrum'
    if not all(p.is_file() and p.stat().st_size for p in (hess, spec)):
        return False
    opt = folder / 'jobex.out'
    if opt.is_file() and spec.stat().st_mtime + 1 < opt.stat().st_mtime:
        return False  # possible inherited pre-Screwer spectrum
    if method == 'aoforce':
        output = text(folder / 'aoforce.out')
        # Match the existing tm_workflow.py acceptance rule: some TM output
        # variants provide files without the literal 'aoforce : all done'.
        return bool(output.strip())
    return (folder / 'numforce.out').is_file() and (folder / 'numforce.out').stat().st_size > 0


def identify(path, base=None, metal=None):
    path = Path(path).expanduser().resolve()
    if metal and base is None:
        base = path  # EF-root convenience: --plan --metal nd[/screwer_nd]
    if base:
        root = Path(base).expanduser().resolve()
        if metal is None:
            path = root
        elif '/' in metal:
            outer, inner = metal.split('/', 1)
            if outer not in METALS:
                raise RecoveryError('--metal path must begin with a lanthanide symbol')
            path = root / metal
        elif metal in METALS or metal in ('la_closedshell', 'lu_closedshell'):
            path = root / metal
        else:
            raise RecoveryError('invalid --metal')
    # Correction cycles are nested to preserve each preceding calculation.
    ancestor = path.parent
    while ancestor.name not in METALS and ancestor.parent != ancestor:
        ancestor = ancestor.parent
    if ancestor.name in METALS and path != ancestor:
        metal = ancestor.name
        parts = path.relative_to(ancestor).parts
        if parts and parts[0] == f'screwer_{metal}' and all(
                part == f'screwer_{metal}_cycle{i}' for i, part in enumerate(parts[1:], 2)):
            return path, metal, 'screwer'
    if path.name in METALS:
        return path, path.name, 'standard'
    if path.name.endswith('_closedshell') and path.name[:-12] in ('la', 'lu'):
        return path, path.name[:-12], 'standard'
    raise RecoveryError('run inside a metal/ or metal/screwer_metal/ directory, or use --base PROJECT --metal SYMBOL[/screwer_SYMBOL]')


def method_from_chain(folder):
    plan = folder / 'chain_plan.txt'
    if plan.exists():
        lines = plan.read_text().splitlines()
        methods = [x for x in lines if x in ('aoforce', 'numforce')]
        if methods:
            return methods[-1]
        if 'jobex_eiger' in lines:
            return 'none'
    records = [folder / 'chain_submission.tsv'] + list(folder.glob('recovery_submission.*.tsv'))
    for path in sorted((p for p in records if p.is_file()), key=lambda p: p.stat().st_mtime, reverse=True):
        entries = [x.split('\t')[0] for x in path.read_text().splitlines()[1:]]
        methods = [x for x in entries if x in ('aoforce', 'numforce')]
        if methods:
            return methods[-1]
    return None


def selected_method(folder, kind, requested):
    from_plan = method_from_chain(folder)
    if kind == 'screwer':
        manifest = folder / 'screwer_prepared.json'
        if not manifest.is_file():
            raise RecoveryError('Screwer manifest missing; cannot automatically recover an unverified correction')
        try:
            info = json.loads(manifest.read_text())
        except (ValueError, OSError) as exc:
            raise RecoveryError(f'invalid Screwer manifest: {exc}') from exc
        expected = info.get('method')
        if expected not in ('aoforce', 'numforce'):
            raise RecoveryError(f'invalid Screwer Hessian method: {expected!r}')
        if from_plan is not None and from_plan != expected:
            raise RecoveryError(f'Screwer method conflicts with chain: {expected} versus {from_plan}')
    else:
        expected = from_plan
    if requested == 'auto':
        if expected is None:
            raise RecoveryError('frequency method cannot be determined from submission history; supply --frequency none|ao|num explicitly')
        return expected
    chosen = {'ao': 'aoforce', 'num': 'numforce'}.get(requested, requested)
    if expected is not None and expected != chosen:
        raise RecoveryError(f'requested {chosen} conflicts with recorded method {expected}')
    if kind == 'screwer' and chosen != expected:
        raise RecoveryError(f'Screwer was prepared for {expected}; recovery cannot change method')
    return chosen


def integer_scf_prepared(folder):
    """Recognize the already-completed occupation fix after ridft_eiger_fix.

    An initial run of ridft_eiger_fix removes $fermi before finishing. Recovery
    must not run that non-idempotent step again. Require the saved pre-removal
    control as evidence rather than inferring readiness from a missing group.
    This is read-only: never overwrite manually edited SCF settings.
    """
    control = text(folder / 'control')
    if re.search(r'(?m)^[ \t]*\$fermi\b', control, re.I):
        return False
    previous = folder / 'control.with_fermi'
    if not previous.is_file() or not re.search(r'(?m)^[ \t]*\$fermi\b', text(previous), re.I):
        raise RecoveryError('control has no $fermi, but control.with_fermi does not verify '
                            'the occupation fix; inspect the control before resuming integer SCF')
    return True


def verified_legacy_resume(folder):
    marker = folder / 'ridft.resume.ok'
    if not marker.is_file():
        return False
    match = re.search(r'(?m)^log=(ridft\.resume\.[A-Za-z0-9_.-]+\.out)$', text(marker))
    return bool(match and scf_ok(folder / match.group(1)))


def stage_plan(folder, kind, method):
    if not folder.is_dir() or folder.is_symlink():
        raise RecoveryError('calculation folder missing or symbolic link')
    for name in ('control', 'coord', 'tm_workflow.py'):
        if not (folder / name).is_file() or (folder / name).is_symlink():
            raise RecoveryError(f'missing or unsafe required file {name}')
    if (folder / 'control.template').exists():
        raise RecoveryError('unexpected control.template')
    manual = not (folder / 'metal.json').is_file()
    if (folder / 'metal.json').is_symlink():
        raise RecoveryError('metal.json symlink is not safe')
    if manual and kind == 'screwer':
        raise RecoveryError('Screwer child missing required metal.json')
    if manual and not geometry_ok(folder):
        raise RecoveryError('manual calculation without metal.json: only an already-converged '
                            'geometry may use Hessian-only recovery without charge/spin metadata')
    if kind == 'screwer':
        if not (folder / 'screwer_prepared.json').is_file():
            raise RecoveryError('Screwer manifest is missing')
    is_restricted = (folder / 'restricted_prepared.json').is_file()
    if folder.name.endswith('_closedshell') and not is_restricted:
        raise RecoveryError('closed-shell directory lacks restricted preparation manifest')
    if is_restricted and kind != 'standard':
        raise RecoveryError('unexpected restricted-Screwer combination')
    if geometry_ok(folder):
        if method == 'none':
            return [], 'geometry converged; no Hessian requested'
        if hessian_ok(folder, method):
            if kind == 'screwer':
                try:
                    from screwer_workflow import spectrum_modes
                    negatives = [(mode, val) for mode, val in spectrum_modes(folder / 'vibspectrum') if val < 0]
                    if negatives:
                        return [], ('current Hessian completed with %d negative mode(s); '
                                    'use lan-screwer --cycle to prepare another correction, not lan-recover' % len(negatives))
                except (RuntimeError, OSError, ValueError):
                    pass
            return [], 'geometry converged; requested final stage already complete'
        return [method], 'geometry converged; no verified current Hessian; submit Hessian only'
    if is_restricted:
        ctrl = text(folder / 'control')
        if re.search(r'(?m)^\s*\$uhf\b', ctrl):
            raise RecoveryError('restricted calculation contains $uhf; investigate manually')
        if scf_ok(folder / 'ridft.out') or any(scf_ok(p) for p in folder.glob('ridft.resume.*.out')):
            return ['jobex_eiger'] + ([method] if method != 'none' else []), 'restricted SCF converged; restart JOBEX'
        if not ((folder / 'define.out').is_file() and re.search(r'define\s+ended\s+normally', text(folder / 'define.out'), re.I)):
            raise RecoveryError('no successful restricted DEFINE; do not automatically regenerate state')
        return ['restricted_ridft_no_define', 'jobex_eiger'] + ([method] if method != 'none' else []), 'restricted DEFINE complete; restart SCF without DEFINE'
    if kind == 'screwer':
        good_define = any(re.search(r'define\s+ended\s+normally', text(folder / n), re.I)
                          for n in ('define.out', 'define.fallback.out'))
        if not good_define:
            raise RecoveryError('Screwer DEFINE not verified; do not automatically re-run DEFINE')
        if scf_ok(folder / 'ridft.out') or any(scf_ok(p) for p in folder.glob('ridft.resume.*.out')):
            stages = ['jobex_eiger']
            reason = 'Screwer SCF converged; restart JOBEX without DEFINE'
        else:
            stages = ['screwer_ridft_jobex_eiger']
            reason = 'Screwer DEFINE complete; restart SCF + JOBEX without DEFINE'
        return stages + ([method] if method != 'none' else []), reason
    if scf_ok(folder / 'ridft.integer.out'):
        return ['jobex_eiger'] + ([method] if method != 'none' else []), 'integer SCF converged; restart JOBEX'
    if not (folder / 'initial_occupations').is_file():
        raise RecoveryError('initial_occupations missing; cannot safely restore the original open-shell state')
    if (scf_ok(folder / 'ridft.out') or scf_ok(folder / 'ridft.retry.out') or
            any(scf_ok(p) for p in folder.glob('ridft.resume.*.out')) or verified_legacy_resume(folder)):
        if integer_scf_prepared(folder):
            stages = ['integer_ridft', 'jobex_eiger']
            reason = 'initial SCF converged; $fermi already removed; resume integer SCF without repeating EIGER/occupation fix'
        else:
            stages = ['eiger_fix_integer_ridft', 'jobex_eiger']
            reason = 'initial SCF converged; perform occupation fix then integer SCF'
        return stages + ([method] if method != 'none' else []), reason
    if not (folder / 'define_setup.ok').is_file():
        raise RecoveryError('initial DEFINE not verified; manual inspection required')
    return ['ridft-recovery', 'eiger_fix_integer_ridft', 'jobex_eiger'] + ([method] if method != 'none' else []), 'restart initial SCF without rerunning DEFINE'


def scheduler_clear(folder):
    if not shutil_which('squeue'):
        raise RecoveryError('squeue unavailable; cannot rule out active jobs')
    # Also match recorded job IDs; some SLURM configurations report N/A for
    # %Z while a dependent job has not been assigned a work directory.
    recorded_ids = set()
    journals = [folder / 'chain_submission.tsv'] + list(folder.glob('recovery_submission.*.tsv'))
    for journal in journals:
        if not journal.is_file():
            continue
        for row in journal.read_text().splitlines()[1:]:
            cells = row.split('\t')
            if len(cells) >= 2 and cells[1].isdecimal():
                recorded_ids.add(cells[1])
    r = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i|%Z|%r'],
                       text=True, capture_output=True)
    if r.returncode:
        raise RecoveryError(f'squeue failed: {r.stderr.strip()}')
    for line in r.stdout.splitlines():
        fields = line.split('|', 2)
        if len(fields) == 3 and (fields[1].strip() == str(folder) or fields[0] in recorded_ids):
            raise RecoveryError(f'active/pending job {fields[0]} for {folder} ({fields[2]}); '
                                'inspect and cancel obsolete dependencies before retrying')


def shutil_which(command):
    from shutil import which
    return which(command)


def submit(folder, metal, kind, method, stages, script_root):
    if not shutil_which('sbatch'):
        raise RecoveryError('sbatch is unavailable')
    scheduler_clear(folder)
    if (folder / '.lan_submission_pending').exists() or (folder / '.lan_submission_lock').exists():
        raise RecoveryError('original submission lock/pending marker is present')
    lock = folder / '.lan_recovery_lock'
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise RecoveryError('.lan_recovery_lock exists; inspect another recovery submission') from exc
    try:
        scheduler_clear(folder)
        marker = folder / '.lan_recovery_pending'
        if marker.exists():
            raise RecoveryError('.lan_recovery_pending exists; previous submission may be ambiguous; inspect scheduler')
        run_stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        journal = folder / f'recovery_submission.{run_stamp}.{os.getpid()}.tsv'
        with journal.open('x') as handle:
            handle.write('stage\tjob_id\tpartition\twalltime\tmemory\tdependency\n')
        dep = None
        for stage in stages:
            cpus = os.environ.get('TM_CPUS', '20')
            if not cpus.isdecimal() or int(cpus) < 1:
                raise RecoveryError('invalid TM_CPUS')
            memory = os.environ.get('HESS_MEM' if stage in ('numforce', 'aoforce') else 'RIDFT_MEM' if stage in ('ridft-recovery', 'integer_ridft', 'eiger_fix_integer_ridft', 'restricted_ridft_no_define') else 'JOBEX_MEM',
                                    '100G' if stage in ('numforce', 'aoforce') else '20G')
            script = 'turbo.resume.slurm' if stage == 'ridft-recovery' else 'turbo.pipeline.slurm'
            script_path = script_root / script
            if not script_path.is_file():
                raise RecoveryError(f'development recovery script missing: {script_path}')
            if stage == 'screwer_ridft_jobex_eiger' and 'screwer_ridft_jobex_eiger)' not in text(script_path):
                raise RecoveryError('development pipeline does not include Screwer recovery')
            exports = ['ALL', f'TM_STAGE={stage}', 'TM_RECOVERY=1']
            if kind == 'screwer':
                exports.append('TM_SCREWER_CHILD=1')
                exports.append('TM_RECOVERY_SCREWER_TOOL=' + str((script_root / 'screwer_workflow.py').resolve()))
            if not (folder / 'metal.json').is_file():
                exports.extend(['TM_MANUAL_RECOVERY=1', f'TM_MANUAL_METAL={metal}'])
            if method == 'none' and stage == 'jobex_eiger':
                exports.append('TM_FINAL_NO_HESSIAN=1')
            opts = ['sbatch', '--parsable', '--partition=general', '--time=2-00:00:00',
                    f'--mem={memory}', f'--cpus-per-task={cpus}',
                    f'--job-name=tm-{metal}-recover-{stage}', '--export=' + ','.join(exports)]
            if dep:
                opts.append(f'--dependency=afterok:{dep}')
            opts.append(str(script_path))
            marker.write_text(f'job_stage={stage}\n')
            result = subprocess.run(opts, cwd=folder, text=True, capture_output=True)
            if result.returncode:
                marker.unlink(missing_ok=True)
                raise RecoveryError(f'sbatch failed for {stage}: {result.stderr.strip()}')
            jobid = result.stdout.strip().split(';')[0]
            if not jobid.isdecimal():
                raise RecoveryError(f'ambiguous sbatch output: {result.stdout.strip()} (pending marker retained)')
            with journal.open('a') as handle:
                handle.write(f'{stage}\t{jobid}\tgeneral\t2-00:00:00\t{memory}\t{("afterok:"+dep) if dep else ""}\n')
            marker.unlink()
            print(f'  submitted {stage}: {jobid}' + (f' (afterok:{dep})' if dep else ''))
            dep = jobid
        return journal
    finally:
        lock.rmdir()


def recovery_warnings(folder, kind):
    warnings = []
    if re.search(r'(?m)^\s*\$fermi\b', text(folder / 'control')):
        if kind == 'screwer':
            warnings.append('Screwer contains $fermi: retained as user-configured; verify intended smearing/occupations before submitting')
        elif scf_ok(folder / 'ridft.integer.out'):
            warnings.append('Integer SCF completed but $fermi remains: verify intended occupation treatment before submitting')
    return warnings


def main(argv=None):
    parser = argparse.ArgumentParser(description='Universal, non-destructive SLURM recovery for normal or Screwer calculations')
    parser.add_argument('method', nargs='?', choices=['ao', 'num', 'aoforce', 'numforce'],
                        help='optional Hessian method, as in legacy lan-recover')
    parser.add_argument('--base', help='EF project root when combined with --metal')
    parser.add_argument('--metal', help='metal symbol or metal/screwer_metal path when using --base')
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--plan', action='store_true', help='read-only recovery plan (default)')
    choice.add_argument('--submit', action='store_true', help='submit validated recovery plan; never runs DEFINE')
    parser.add_argument('--frequency', default='auto', choices=['auto', 'none', 'ao', 'aoforce', 'num', 'numforce'],
                        help='auto uses recorded original stage; explicitly choose when missing history')
    parser.add_argument('--json', action='store_true', help='JSON plan on stdout; read-only (cannot combine with --submit)')
    args = parser.parse_args(argv)
    try:
        if args.json and args.submit:
            raise RecoveryError('--json cannot be used with --submit')
        folder, metal, kind = identify(Path.cwd(), args.base, args.metal)
        if args.method and args.frequency != 'auto':
            raise RecoveryError('choose either positional method or --frequency, not both')
        method = selected_method(folder, kind, args.method or args.frequency)
        stages, reason = stage_plan(folder, kind, method)
        warnings = recovery_warnings(folder, kind)
        output = {'directory': str(folder), 'metal': metal, 'kind': kind, 'frequency': method,
                  'stages': stages, 'reason': reason, 'warnings': warnings,
                  'writes_inputs_on_submission': False}
        if args.json:
            print(json.dumps(output, indent=2))
        else:
            print(f'RECOVERY PLAN: {folder}\n  kind: {kind}; Hessian: {method}\n  reason: {reason}')
            print('  stages: ' + (' -> '.join(stages) if stages else '(none)'))
            for warning in warnings:
                print('  WARNING: ' + warning)
        if args.submit:
            if not stages:
                print('Nothing to recover; no jobs submitted')
            else:
                journal = submit(folder, metal, kind, method, stages, Path(__file__).resolve().parent)
                print(f'Submitted recovery; journal: {journal}')
        return 0
    except RecoveryError as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2

if __name__ == '__main__':
    sys.exit(main())
