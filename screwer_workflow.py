#!/usr/bin/env python3
"""Opt-in one-pass imaginary-mode correction, never changing parent calculations.

TURBOMOLE 7.9 VACC.  All work is in a NEW <metal>/screwer_<metal> directory.
Templates screwer.inp and def_screwer.inp belong to the project owner.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from setup_lanthanides import SPIN, make_define, read_charge

METALS = ('la', 'ce', 'pr', 'nd', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu')
STAGES = ('screwer_jobex_eiger',)
NEWCOORD_MODE = re.compile(r'shifted\s+along\s+normal\s+mode\s+(\d+)', re.I)


def correction_dir(parent, metal):
    return parent / f'screwer_{metal}'


def legacy_dir(parent):
    return parent / 'screwer'



def execute(args, *, cwd, input_text=None, timeout=180):
    return subprocess.run(args, cwd=cwd, input=input_text, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          timeout=timeout, check=False)


def spectrum_modes(path):
    """Parse at most eight *mode rows* of a standard TURBOMOLE vibspectrum."""
    if not path.is_file():
        raise RuntimeError(f'missing {path}')
    modes = []
    with path.open(errors='replace') as file:
        heading = file.readline()
        if not heading.lstrip().startswith('$vibrational spectrum'):
            raise RuntimeError(f'{path}: invalid vibrational spectrum header')
        for line in file:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            if line.startswith('$') or len(modes) == 8:
                break
            cols = line.split()
            if not cols or not cols[0].isdigit():
                raise RuntimeError(f'{path}: invalid mode row: {line.strip()}')
            number = int(cols[0])
            if number != len(modes) + 1:
                raise RuntimeError(f'{path}: expected mode {len(modes)+1}, got {number}')
            if len(cols) < 4:
                raise RuntimeError(f'{path}: incomplete mode {number}')
            value = cols[1] if cols[1][0] in '+-.0123456789' else cols[2]
            try:
                freq = float(value.replace('D', 'E').replace('d', 'e'))
            except ValueError as error:
                raise RuntimeError(f'{path}: invalid wave number in mode {number}') from error
            if not (-float('inf') < freq < float('inf')):
                raise RuntimeError(f'{path}: non-finite wave number in mode {number}')
            modes.append((number, freq))
    if not modes:
        raise RuntimeError(f'{path}: no mode rows')
    return modes



def latest_submission(parent):
    records = [p for p in ([parent / 'chain_submission.tsv'] +
                           list(parent.glob('recovery_submission.*.tsv'))) if p.is_file()]
    return max(records, key=lambda p: p.stat().st_mtime_ns) if records else None


def verify_parent(parent):
    """Don't trigger on an unfinished run or a previous spectrum."""
    state = parent / 'workflow_status.json'
    try:
        data = json.loads(state.read_text())
    except (FileNotFoundError, ValueError) as error:
        raise RuntimeError(f'{parent}: no valid workflow_status.json') from error
    if data.get('stage') != 'DONE' or data.get('status') != 'SUCCESS':
        raise RuntimeError(f'{parent}: frequency job has not recorded DONE/SUCCESS')
    record = latest_submission(parent)
    if record is None:
        raise RuntimeError(f'{parent}: no tracked submission record')
    fields = record.read_text().splitlines()
    if not any(line.split('\t')[0] in ('numforce', 'aoforce') for line in fields[1:]):
        raise RuntimeError(f'{parent}: latest submission did not request a Hessian')
    spec = parent / 'vibspectrum'
    if not spec.is_file() or spec.stat().st_mtime + 1 < record.stat().st_mtime:
        raise RuntimeError(f'{parent}: missing/stale vibspectrum for current attempt')
    # An active job in this calculation directory may still be mutating inputs.
    if shutil.which('squeue') is None:
        raise RuntimeError('squeue unavailable: cannot rule out concurrent calculations')
    queue = execute(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i|%Z'], cwd=parent)
    if queue.returncode:
        raise RuntimeError(f'squeue failed: {queue.stdout.strip()}')
    current_gate = os.environ.get('SLURM_JOB_ID', '')
    if any((row.partition('|')[2].strip() == str(parent.resolve()) and
            row.partition('|')[0].strip() != current_gate)
           for row in queue.stdout.splitlines() if '|' in row):
        raise RuntimeError(f'{parent}: an active/pending job still uses this directory')
    return record


def parse_coordinates(output, parent_coord):
    """Validate sdg -b newcoord output before making a new coordinate file."""
    rows = []
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 4:
            raise RuntimeError(f'newcoord has non-coordinate row: {line!r}')
        try:
            for field in fields[:3]:
                value = float(field.replace('D', 'E'))
                if not (-float('inf') < value < float('inf')):
                    raise ValueError('non-finite coordinate')
        except ValueError as error:
            raise RuntimeError(f'newcoord has invalid Cartesian values: {line!r}') from error
        if not re.fullmatch(r'[A-Za-z]{1,3}', fields[3]):
            raise RuntimeError(f'newcoord has invalid atom label: {line!r}')
        rows.append((fields[3].lower(), line))
    orig = []
    for line in parent_coord.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] != '$coord' and not line.lstrip().startswith('$'):
            orig.append(fields[3].lower())
    if not orig or [row[0] for row in rows] != orig:
        raise RuntimeError('displaced coordinates do not match the original atom count and order')
    return '$coord\n' + '\n'.join(line for _, line in rows) + '\n$end\n'


def prepare(base, metal, method, *, min_negative=0):
    parent = base / metal
    if not parent.is_dir():
        raise RuntimeError(f'{parent}: no metal directory')
    record = verify_parent(parent)
    modes = spectrum_modes(parent / 'vibspectrum')
    negative = [(n, freq) for n, freq in modes if freq < -min_negative]
    if not negative:
        print(f'{metal}: no frequencies below -{min_negative:g} cm^-1 in first 8 modes; no screwer folder created')
        return None
    template_name = 'def_screwer_gd.inp' if metal == 'gd' else 'def_screwer.inp'
    for file in ('screwer.inp', template_name):
        if not (base / file).is_file() or not (base / file).stat().st_size:
            raise RuntimeError(f'{base / file} required (user-authored template)')
    try:
        meta = json.loads((parent / 'metal.json').read_text())
        generated_define = make_define((base / template_name).read_text(), metal, SPIN[metal])
        screwer_charge = read_charge(generated_define)
    except (OSError, ValueError, KeyError) as error:
        raise RuntimeError(f'{metal}: could not validate screwer DEFINE template: {error}') from error
    if screwer_charge is None or screwer_charge != meta.get('charge'):
        raise RuntimeError(f'{metal}: screwer DEFINE EHT charge {screwer_charge} does not match the original project charge {meta.get("charge")}; refusing')
    # Require that the requested correction uses the existing Hessian method.
    content = record.read_text().splitlines()
    recorded_methods = [line.split('\t')[0] for line in content[1:] if line.split('\t')[0] in ('numforce', 'aoforce')]
    if not recorded_methods or recorded_methods[-1] != method:
        raise RuntimeError(f'{metal}: current Hessian is {recorded_methods[-1] if recorded_methods else "unknown"}, not {method}')
    child = correction_dir(parent, metal)
    if child.exists() or legacy_dir(parent).exists():
        raise RuntimeError(f'{parent}: screwer correction directory already exists; left untouched')
    if not shutil.which('cpc') or not shutil.which('screwer') or not shutil.which('sdg') or not shutil.which('kdg'):
        raise RuntimeError('cpc, screwer, sdg and kdg must be in PATH (load turbomole/7.9)')
    print(f'{metal}: negative modes {negative}; creating {child}', flush=True)
    copy = execute(['cpc', child.name], cwd=parent)
    if copy.returncode or not child.is_dir() or child.is_symlink():
        raise RuntimeError(f'cpc {child.name} failed ({copy.returncode}): {copy.stdout.strip()}; inspect folder manually')
    # Never risk following a link back to the completed parent calculation.
    for name in ('control', 'coord', 'metal.json', 'alpha', 'beta', 'mos',
                 'screwer.inp', 'def_screwer.inp'):
        if (child / name).is_symlink():
            raise RuntimeError(f'{child / name}: cpc left a symbolic link; refusing to modify linked files')
    shutil.copy2(base / 'screwer.inp', child / 'screwer.inp')
    (child / 'def_screwer.inp').write_text(generated_define)
    # cpc copies TURBOMOLE data, not necessarily our custom metadata.
    # The correction retains the parent's Ln/charge/spin; never invent values.
    if not (parent / 'metal.json').is_file():
        raise RuntimeError(f'{parent}: missing metal.json')
    shutil.copy2(parent / 'metal.json', child / 'metal.json')
    script_dir = Path(__file__).resolve().parent
    for name in ('tm_workflow.py', 'turbo.pipeline.slurm',
                 'turbo.resume.slurm', 'screwer_workflow.py', 'setup_lanthanides.py'):
        shutil.copy2(script_dir / name, child / name)
    screwer_input = (child / 'screwer.inp').read_text()
    out = execute(['screwer'], cwd=child, input_text=screwer_input)
    (child / 'screwer.out').write_text(out.stdout)
    if out.returncode:
        raise RuntimeError(f'{child}: screwer failed ({out.returncode}); see screwer.out')
    control_text = (child / 'control').read_text()
    match = re.search(r'(?ims)^\s*\$newcoord\b(.*?)(?=^\s*\$|\Z)', control_text)
    if not match:
        raise RuntimeError(f'{child}: screwer produced no $newcoord (see screwer.out)')
    mode_match = NEWCOORD_MODE.search(match.group(1))
    if not mode_match:
        raise RuntimeError(f'{child}: cannot verify displaced mode from $newcoord comment')
    mode = int(mode_match[1])
    if mode not in dict(negative):
        raise RuntimeError(f'{child}: screwer displaced mode {mode}, which is not among verified negative modes {negative}')
    coords = execute(['sdg', '-b', 'newcoord'], cwd=child)
    if coords.returncode:
        raise RuntimeError(f'{child}: sdg -b newcoord failed: {coords.stdout.strip()}')
    new_coord = parse_coordinates(coords.stdout, parent / 'coord')
    (child / 'coord.screwer.tmp').write_text(new_coord)
    os.replace(child / 'coord.screwer.tmp', child / 'coord')
    # Perform exactly the four user-specified kdg operations, in that order.
    for name in ('newcoord', 'coord', 'symmetry', 'optimize'):
        result = execute(['kdg', name], cwd=child)
        if result.returncode:
            raise RuntimeError(f'{child}: kdg {name} failed: {result.stdout.strip()}')
    # Prevent copying completed parent results into the child run's history.
    # This is a NEW directory; the original parent files remain untouched.
    copied_results = ('GEO_OPT_CONVERGED',
                      'GEO_OPT_FAILED', 'converged', 'not.converged', 'ridft.resume.ok',
                      'workflow_status.json', 'workflow_history.jsonl',
                      'initial_occupations', 'control.with_fermi', 'chain_submission.tsv',
                      'chain_plan.txt', '.lan_submission_pending', 'occupation_pre.json',
                      'occupation_post.json', 'define_setup.ok')
    archive = child / 'parent_results'
    for filename in copied_results:
        p = child / filename
        if p.is_file():
            archive.mkdir(exist_ok=True)
            os.replace(p, archive / filename)
    if (child / '.lan_submission_lock').exists():
        raise RuntimeError(f'{child}: a copied submission lock was found; inspect manually')
    # The screwer continuation must start from the inherited converged orbitals.
    # DEFINE is checked for changes to those files and to shell occupations.
    import tm_workflow
    inherited = tm_workflow.screwer_inheritance(child)
    info = {'workflow_version': 2, 'inherited_state': inherited,
            'parent': str(parent), 'metal': metal, 'mode': mode,
            'define_source': template_name, 'charge': screwer_charge,
            'frequency_cm1': dict(negative)[mode], 'negative_modes_first8': negative,
            'method': method, 'prepared_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'one_pass_only': True}
    (child / 'screwer_prepared.json').write_text(json.dumps(info, indent=2) + '\n')
    print(f'{metal}: screwer prepared mode {mode} ({dict(negative)[mode]:g} cm^-1) in {child}', flush=True)
    return child


def submit(child, method):
    """Submit one 2-allocation JOBEX/EIGER + Hessian chain; no automatic repeat."""
    manifest = child / 'screwer_prepared.json'
    if not manifest.is_file():
        raise RuntimeError(f'{child}: no screwer preparation record')
    info = json.loads(manifest.read_text())
    if info['method'] != method:
        raise RuntimeError(f'{child}: prepared for {info["method"]}, not {method}')
    if info.get('workflow_version') != 2 or not info.get('inherited_state'):
        raise RuntimeError(f'{child}: old screwer preparation; review manually (never resubmit under a changed stage plan)')
    if any((child / name).exists() for name in ('chain_submission.tsv', 'chain_plan.txt', '.lan_submission_pending', '.lan_submission_lock')):
        raise RuntimeError(f'{child}: already submitted or ambiguous; refusing duplicate submission')
    for name in ('control', 'coord', 'def_screwer.inp', 'turbo.pipeline.slurm',
                 'tm_workflow.py', 'setup_lanthanides.py'):
        if not (child / name).is_file():
            raise RuntimeError(f'{child}: missing {name}')
    # Prepared child folders keep their own script copy. Refuse an older
    # copy that would skip the requested pre-JOBEX RIDFT; never overwrite it.
    script = (child / 'turbo.pipeline.slurm').read_text()
    if "screwer continuation SCF before JOBEX" not in script:
        raise RuntimeError(f'{child}: older Screwer SLURM script lacks pre-JOBEX RIDFT; inspect/update it manually before submitting')
    if shutil.which('sbatch') is None:
        raise RuntimeError('sbatch not found; submit on VACC')
    stages = (*STAGES, method)
    # Lock the child folder, including when submit is invoked automatically.
    lock = child / '.lan_submission_lock'
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise RuntimeError(f'{child}: submission lock already exists') from error
    dep = None
    submitted = 0
    try:
        for i, stage in enumerate(stages):
            mem = os.environ.get('HESS_MEM' if stage in ('numforce', 'aoforce') else 'JOBEX_MEM' if stage == 'screwer_jobex_eiger' else 'RIDFT_MEM',
                                 '100G' if stage in ('numforce', 'aoforce') else '20G')
            cpus = os.environ.get('TM_CPUS', '20')
            if not cpus.isdecimal() or int(cpus) < 1:
                raise RuntimeError('invalid TM_CPUS')
            opts = ['sbatch', '--parsable', '--partition=general', '--time=2-00:00:00',
                    f'--mem={mem}', f'--cpus-per-task={cpus}',
                    f'--job-name=tm-{info["metal"]}-screwer-{stage}',
                    f'--export=ALL,TM_STAGE={stage},TM_SCREWER_CHILD=1']
            if dep:
                opts.append(f'--dependency=afterok:{dep}')
            opts.append('turbo.pipeline.slurm')
            (child / '.lan_submission_pending').write_text(f'stage={stage}\n')
            result = execute(opts, cwd=child)
            if result.returncode:
                (child / '.lan_submission_pending').unlink(missing_ok=True)
                raise RuntimeError(f'{child}: sbatch failed at {stage}: {result.stdout.strip()}')
            jobid = result.stdout.strip().split(';')[0]
            if not jobid.isdecimal():
                raise RuntimeError(f'{child}: ambiguous sbatch output for {stage}; pending marker retained')
            if i == 0:
                (child / 'chain_plan.txt').write_text('\n'.join(stages) + '\n')
                (child / 'chain_submission.tsv').write_text('stage\tjob_id\tpartition\twalltime\tmemory\n')
            with (child / 'chain_submission.tsv').open('a') as log:
                log.write(f'{stage}\t{jobid}\tgeneral\t2-00:00:00\t{mem}\n')
            (child / '.lan_submission_pending').unlink()
            dep = jobid
            submitted += 1
            print(f'{child}: submitted {stage} job {jobid}', flush=True)
    finally:
        lock.rmdir()
    return submitted


def clean_inherited_frequencies(child):
    """Remove stale frequency files from the fresh correction directory only.

    Delay until Hessian stage: DEFINE and screwer may need original spectra
    and control's file references while preparing the corrected calculation.
    """
    if not (child / 'screwer_prepared.json').is_file():
        raise RuntimeError('refusing to clear spectra outside a prepared screwer correction')
    archive = child / 'parent_results'
    for name in ('vibspectrum', 'hessian', 'vib_normal_modes', 'nprhessian',
                 'nprvibspectrum', 'npr_vib_normal_modes'):
        p = child / name
        if p.is_file():
            archive.mkdir(exist_ok=True)
            # Archived inherited spectra must not overwrite an existing backup.
            target = archive / name
            if target.exists():
                raise RuntimeError(f'parent frequency backup already exists: {target}')
            os.replace(p, target)
    print(f'{child}: inherited frequency files archived before new Hessian', flush=True)


def extract_noproj(child):
    """Extract absent $noproj data files from control, without replacing data."""
    if not (child / 'control').is_file():
        raise RuntimeError(f'{child}: no control file')
    names = (('nprhessian', 'nprhessian'),
             ('nprvibrational spectrum', 'nprvibspectrum'),
             ('nprvibrational normal modes', 'npr_vib_normal_modes'))
    found = []
    for group, filename in names:
        target = child / filename
        if target.is_file() and target.stat().st_size:
            found.append(filename)
            continue
        if shutil.which('sdg') is None:
            raise RuntimeError('sdg is not available; load turbomole/7.9')
        result = execute(['sdg', group], cwd=child)
        if result.returncode:
            print(f'{child}: no {group} group available', flush=True)
            continue
        txt = result.stdout
        first_line = txt.splitlines()[0] if txt.splitlines() else ''
        if (not txt.lstrip().lower().startswith('$' + group.lower()) or
                len(txt.splitlines()) < 3 or 'file=' in first_line.lower()):
            print(f'{child}: {group} unavailable or only a file-reference; not writing {filename}', flush=True)
            continue
        temp = child / (filename + '.tmp.' + str(os.getpid()))
        temp.write_text(txt)
        os.replace(temp, target)
        found.append(filename)
    print(f'{child}: $noproj files available: {", ".join(found) if found else "none"}', flush=True)
    return found


def main():
    parser = argparse.ArgumentParser(
        description='Opt-in screwer reoptimization of completed imaginary-frequency calculations (one pass only).',
        epilog=('Examples: lan-screwer num --metals gd --prepare-only; '
                'lan-screwer num --metals gd; lan-screwer --export-only --metals gd'))
    parser.add_argument('method', nargs='?', choices=['num', 'ao', 'numforce', 'aoforce'],
                        help='Hessian method from the completed original calculation')
    parser.add_argument('--base', type=Path, default=Path.cwd(), help='project directory')
    parser.add_argument('--metals', nargs='+', choices=METALS, required=True, help='metals to process')
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--prepare-only', action='store_true', help='create screwer_<metal>/ inputs but do not submit')
    action.add_argument('--submit-prepared', action='store_true', help='submit an already prepared screwer_<metal>/ directory')
    action.add_argument('--export-only', action='store_true', help='extract $noproj data for completed screwer corrections')
    action.add_argument('--clean-frequency', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--min-imag', type=float, default=0.0,
                        help='only trigger on frequencies strictly below -VALUE cm^-1 (default 0)')
    args = parser.parse_args()
    if not (0 <= args.min_imag < float('inf')):
        parser.error('--min-imag must be a finite nonnegative number')
    method = {'num': 'numforce', 'ao': 'aoforce', 'numforce': 'numforce', 'aoforce': 'aoforce'}.get(args.method)
    if not (args.export_only or args.clean_frequency) and not method:
        parser.error('choose num or ao (unless using --export-only)')
    base = args.base.resolve()
    failures = 0
    for metal in args.metals:
        try:
            if args.clean_frequency:
                clean_inherited_frequencies(correction_dir(base / metal, metal))
            elif args.export_only:
                new_child = correction_dir(base / metal, metal)
                extract_noproj(new_child if new_child.is_dir() else legacy_dir(base / metal))
            elif args.submit_prepared:
                submit(correction_dir(base / metal, metal), method)
            else:
                dest = prepare(base, metal, method, min_negative=args.min_imag)
                if dest is not None and not args.prepare_only:
                    submit(dest, method)
        except (RuntimeError, OSError, subprocess.TimeoutExpired, ValueError) as error:
            print(f'ERROR {metal}: {error}', file=sys.stderr)
            failures += 1
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
