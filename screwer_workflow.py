#!/usr/bin/env python3
"""Opt-in, bounded imaginary-mode correction, never changing parent calculations.

TURBOMOLE 7.9 VACC. Each cycle has its own nested directory.
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
import tempfile


METALS = ('la', 'ce', 'pr', 'nd', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu')
STAGES = ('screwer_jobex_eiger',)
NEWCOORD_MODE = re.compile(r'shifted\s+along\s+normal\s+mode\s+(\d+)', re.I)


def correction_dir(parent, metal, cycle=1):
    return parent / (f'screwer_{metal}' if cycle == 1 else f'screwer_{metal}_cycle{cycle}')


def cycle_directory(base, metal, cycle):
    if not isinstance(cycle, int) or not 1 <= cycle <= 10:
        raise ValueError('cycle must be an integer from 1 to 10')
    path = base / metal
    for n in range(1, cycle + 1):
        path = correction_dir(path, metal, n)
    return path


def cycle_parent(base, metal, cycle):
    return (base / metal) if cycle == 1 else cycle_directory(base, metal, cycle - 1)


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


def verify_parent(parent, method):
    """Check a completed Hessian calculation, with or without lan-run history.

    Missing BOTH automation records and workflow status permits a conservative
    manual-results check.  Partial/contradictory automation history is never
    silently reclassified as a manual calculation.
    """
    state = parent / 'workflow_status.json'
    record = latest_submission(parent)
    automation_markers = ('workflow_status.json', 'workflow_history.jsonl',
                          'chain_submission.tsv', 'chain_plan.txt')
    has_history = record is not None or any((parent / name).exists()
                                            for name in automation_markers)
    if has_history:
        try:
            data = json.loads(state.read_text())
        except (FileNotFoundError, ValueError) as error:
            raise RuntimeError(f'{parent}: incomplete automation history; review manually') from error
        if data.get('stage') != 'DONE' or data.get('status') != 'SUCCESS':
            raise RuntimeError(f'{parent}: frequency job has not recorded DONE/SUCCESS')
        if record is None:
            raise RuntimeError(f'{parent}: no tracked submission record for automated result')
        fields = record.read_text().splitlines()
        recorded_methods = [line.split('\t')[0] for line in fields[1:]
                            if line.split('\t')[0] in ('numforce', 'aoforce')]
        if not recorded_methods or recorded_methods[-1] != method:
            raise RuntimeError(f'{parent}: current Hessian is {recorded_methods[-1] if recorded_methods else "unknown"}, not {method}')
        spec = parent / 'vibspectrum'
        if not spec.is_file() or spec.stat().st_mtime + 1 < record.stat().st_mtime:
            raise RuntimeError(f'{parent}: missing/stale vibspectrum for current attempt')
        # A DONE flag alone is insufficient: reject an incomplete or inherited
        # frequency result before starting another computationally costly cycle.
        for name in ('control', 'coord', 'gradient', 'jobex.out', 'hessian',
                     'vibspectrum', 'aoforce.out' if method == 'aoforce' else 'numforce.out'):
            file = parent / name
            if not file.is_file() or not file.stat().st_size:
                raise RuntimeError(f'{parent}: verified history but missing/empty {name}')
        if (parent / 'vibspectrum').stat().st_mtime + 1 < (parent / 'jobex.out').stat().st_mtime:
            raise RuntimeError(f'{parent}: vibspectrum predates the latest JOBEX')
        source = 'automation'
    else:
        # A manual calculation does not have a workflow_status.json or
        # chain_submission.tsv. Verify tangible TURBOMOLE results instead.
        if (parent / 'GEO_OPT_FAILED').exists() or (parent / 'not.converged').exists():
            raise RuntimeError(f'{parent}: failure marker present; verify geometry optimization')
        if not (parent / 'GEO_OPT_CONVERGED').is_file():
            raise RuntimeError(f'{parent}: missing GEO_OPT_CONVERGED for manual calculation')
        required = ('control', 'coord', 'basis', 'auxbasis', 'gradient',
                    'jobex.out', 'hessian', 'vibspectrum',
                    'aoforce.out' if method == 'aoforce' else 'numforce.out')
        for name in required:
            file = parent / name
            if not file.is_file() or file.stat().st_size == 0:
                raise RuntimeError(f'{parent}: missing/empty {name} for manual {method} calculation')
        if not ((parent / 'alpha').is_file() and (parent / 'beta').is_file()) and not (parent / 'mos').is_file():
            raise RuntimeError(f'{parent}: missing inherited alpha/beta or mos orbitals')
        if (parent / 'vibspectrum').stat().st_mtime + 1 < (parent / 'jobex.out').stat().st_mtime:
            raise RuntimeError(f'{parent}: vibspectrum predates latest JOBEX; review stale Hessian')
        source = 'manual'

    # Read and validate the same spectrum under both paths, before creating
    # anything.  An active job may still be changing the source calculation.
    spectrum_modes(parent / 'vibspectrum')
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
    if source == 'manual':
        print(f'WARNING {parent.name}: no lan-run history; accepted manually completed '
              f'{method} result using TURBOMOLE optimization and Hessian files', flush=True)
    return source


def parse_coordinates(output, parent_coord):
    """Validate sdg -b newcoord output before making a new coordinate file."""
    rows = []
    for line in output.splitlines():
        fields = line.split()
        if not fields or line.lstrip().startswith('#'):
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
    # A TURBOMOLE coord file may also contain $redundant and other data groups.
    # Only atom records inside $coord belong in the atom-order comparison.
    orig = []
    in_coord = False
    for line in parent_coord.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith('$coord'):
            in_coord = True
            continue
        if in_coord and stripped.startswith('$'):
            break
        if in_coord and stripped and not stripped.startswith('#'):
            fields = stripped.split()
            if len(fields) >= 4:
                orig.append(fields[3].lower())
    if not orig or [row[0] for row in rows] != orig:
        raise RuntimeError('displaced coordinates do not match the original atom count and order')
    return '$coord\n' + '\n'.join(line for _, line in rows) + '\n$end\n'


def prepare(base, metal, method, *, min_negative=0, define_backend="legacy", cycle=1, max_cycles=1):
    if not (1 <= cycle <= max_cycles <= 10):
        raise RuntimeError('require 1 <= cycle <= max_cycles <= 10')
    parent = cycle_parent(base, metal, cycle)
    if not parent.is_dir():
        raise RuntimeError(f'{parent}: no metal directory')
    parent_source = verify_parent(parent, method)
    modes = spectrum_modes(parent / 'vibspectrum')
    negative = [(n, freq) for n, freq in modes if freq < -min_negative]
    if not negative:
        print(f'{metal}: no frequencies below -{min_negative:g} cm^-1 in first 8 modes; no screwer folder created')
        return None
    template_name = 'def_screwer_gd.inp' if metal == 'gd' else 'def_screwer.inp'
    for file in ('screwer.inp', template_name):
        if not (base / file).is_file() or not (base / file).stat().st_size:
            raise RuntimeError(f'{base / file} required (user-authored template)')
    # Screwer DEFINE only changes the displaced coordinates.  It does not
    # establish EHT charge, basis, or orbital occupations again.
    generated_define = (base / template_name).read_text()
    metal_meta = parent / 'metal.json'
    if metal_meta.is_file():
        try:
            meta = json.loads(metal_meta.read_text())
        except ValueError as error:
            raise RuntimeError(f'{metal_meta}: invalid metal metadata') from error
        if not isinstance(meta, dict) or meta.get('metal') != metal:
            raise RuntimeError(f'{metal_meta}: incorrect metal identity')
    else:
        # Metadata is for automation reporting, NOT an input to the Screwer
        # calculation.  Do not guess the charge or spin of a manual project.
        meta = {'metal': metal, 'charge': None, 'unpaired': None,
                'source': 'manual TURBOMOLE control and orbitals'}
        print(f'WARNING {metal}: no metal.json; child metadata will not infer '
              'charge or spin from the Screwer DEFINE template', flush=True)
    child = correction_dir(parent, metal, cycle)
    if child.exists():
        raise RuntimeError(f'{child}: correction directory already exists; left untouched')
    if legacy_dir(parent).exists():
        print(f'NOTE {metal}: preserving existing legacy {legacy_dir(parent)}', flush=True)
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
    # Optional DEFINE fallback: a separate no-ired input, copied now so queued
    # jobs use the tested project template even if the master later changes.
    fallback = 'def_screwer_gd_ln.inp' if metal == 'gd' else 'def_screwer_ln.inp'
    if (base / fallback).is_file():
        shutil.copy2(base / fallback, child / fallback)
    else:
        print(f'NOTE {metal}: no {fallback}; DEFINE will stop if its primary input fails', flush=True)
    # cpc copies TURBOMOLE data, not necessarily our custom metadata.
    # The correction retains the parent's Ln/charge/spin; never invent values.
    if metal_meta.is_file():
        shutil.copy2(metal_meta, child / 'metal.json')
    else:
        (child / 'metal.json').write_text(json.dumps(meta, indent=2) + '\n')
    script_dir = Path(__file__).resolve().parent
    for name in ('tm_workflow.py', 'turbo.pipeline.slurm', 'adaptive_define.py',
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
            'parent': str(parent), 'parent_source': parent_source, 'metal': metal, 'mode': mode,
            'define_source': template_name, 'charge': meta.get('charge'),
            'frequency_cm1': dict(negative)[mode], 'negative_modes_first8': negative,
            'method': method, 'prepared_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'one_pass_only': cycle == 1 and max_cycles == 1, 'define_backend': define_backend,
            'cycle': cycle, 'max_cycles': max_cycles, 'project_root': str(base),
            'min_negative_cm1': float(min_negative)}
    (child / 'screwer_prepared.json').write_text(json.dumps(info, indent=2) + '\n')
    print(f'{metal}: screwer prepared mode {mode} ({dict(negative)[mode]:g} cm^-1) in {child}', flush=True)
    return child



# Files DEFINE may update while processing the displaced geometry.  Preserve
# the pre-DEFINE inputs for a clean fallback attempt, *only inside the child*.
DEFINE_INPUTS = ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta',
                 'mos', 'statistics', 'intdef', 'wherefrom', 'gradient',
                 'energy', 'hessapprox', 'job.start', 'job.last')


def validate_screwer_define(child, result):
    """DEFINE can return 0 after a stopped menu; verify normal completion.

    The geometry and atom data groups can be present inline or through
    file= references.  This intentionally does NOT compare orbital contents,
    shell occupations, or the inherited control to their original values.
    """
    if result.returncode:
        raise RuntimeError(f'DEFINE exited with code {result.returncode}')
    if not re.search(r'(?i)define\s*:\s*all\s+done|define\s+ended\s+normally',
                     result.stdout):
        raise RuntimeError('DEFINE did not report normal completion')
    control = child / 'control'
    if not control.is_file() or not control.stat().st_size:
        raise RuntimeError('DEFINE did not produce a nonempty control file')
    data = control.read_text(errors='replace')
    for group in ('coord', 'atoms'):
        if not re.search(r'(?im)^[ \t]*\$' + group + r'(?=\s|$)', data):
            raise RuntimeError(f'DEFINE left control without ${group}')
    coord = child / 'coord'
    if not coord.is_file() or not coord.stat().st_size:
        raise RuntimeError('DEFINE left the coordinate file missing/empty')
    # The basis and MO files are required to continue RIDFT, but may be
    # rewritten by DEFINE (a legitimate continuation behavior).
    for name in ('basis', 'auxbasis'):
        if not (child / name).is_file() or not (child / name).stat().st_size:
            raise RuntimeError(f'DEFINE left {name} missing/empty')
    if not ((child / 'alpha').is_file() and (child / 'beta').is_file()) and not (child / 'mos').is_file():
        raise RuntimeError('DEFINE left no alpha/beta or mos orbitals for continuation')


def run_define_with_fallback(child):
    """Run Screwer DEFINE once; retry without ired from a clean snapshot.

    Unlike initial lan-run, Screwer DEFINE may update control and MOs.  Its
    only validation is normal program completion and usable required inputs.
    """
    manifest = child / 'screwer_prepared.json'
    if not manifest.is_file():
        raise RuntimeError(f'{child}: no Screwer preparation manifest')
    info = json.loads(manifest.read_text())
    metal = info.get('metal')
    backend = info.get('define_backend', 'legacy')
    if backend not in ('legacy', 'adaptive'):
        raise RuntimeError(f'{manifest}: unrecognized DEFINE backend {backend!r}')
    if metal not in METALS:
        raise RuntimeError(f'{manifest}: invalid metal')
    primary = child / 'def_screwer.inp'
    fallback_name = 'def_screwer_gd_ln.inp' if metal == 'gd' else 'def_screwer_ln.inp'
    fallback = child / fallback_name
    if not primary.is_file() or not primary.stat().st_size:
        raise RuntimeError(f'{primary}: missing/empty primary DEFINE input')
    with tempfile.TemporaryDirectory(prefix='.screwer_define_snapshot_', dir=child) as scratch:
        backups = Path(scratch)
        existed = set()
        for name in DEFINE_INPUTS:
            path = child / name
            if path.is_symlink():
                raise RuntimeError(f'{path}: refusing to snapshot a symlink')
            if path.is_file():
                shutil.copy2(path, backups / name)
                existed.add(name)

        def restore():
            for name in DEFINE_INPUTS:
                target = child / name
                if name in existed:
                    if target.is_symlink():
                        target.unlink()
                    shutil.copy2(backups / name, target)
                elif target.exists() and target.is_file():
                    target.unlink()  # newly generated DEFINE input/output only

        def attempt(template, logfile):
            if backend == 'adaptive':
                from adaptive_define import define_supervised
                result = define_supervised(child, template.read_text())
            else:
                result = execute(['define'], cwd=child, input_text=template.read_text(), timeout=None)
            (child / logfile).write_text(result.stdout)
            try:
                validate_screwer_define(child, result)
                return None
            except RuntimeError as error:
                return str(error)

        first_error = attempt(primary, 'define.out')
        if first_error is None:
            print(f'{child}: DEFINE succeeded using {primary.name} (backend={backend})', flush=True)
            return
        print(f'WARNING {child}: {primary.name} failed: {first_error}', flush=True)
        restore()
        if not fallback.is_file() or not fallback.stat().st_size:
            raise RuntimeError(f'{child}: no {fallback_name} to retry DEFINE; prepared inputs restored; see define.out')
        if re.search(r'(?im)^\s*ired(?:\s|$)', fallback.read_text()):
            raise RuntimeError(f'{fallback}: fallback still requests ired; prepared inputs restored')
        print(f'{child}: retrying DEFINE with {fallback_name} from restored inputs', flush=True)
        second_error = attempt(fallback, 'define.fallback.out')
        if second_error is not None:
            restore()
            raise RuntimeError(f'{child}: DEFINE fallback also failed: {second_error}; '
                               'prepared inputs restored; see define.out and define.fallback.out')
        print(f'{child}: DEFINE fallback succeeded using {fallback_name}', flush=True)


def submit(child, method):
    """Submit one 2-allocation Screwer cycle; next cycle uses a fresh directory."""
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
    if '--define-fallback' not in script:
        raise RuntimeError(f'{child}: copied SLURM script lacks DEFINE fallback; inspect/update it manually')
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


def advance_cycle(folder):
    """Advance from a *completed* correction, subject to its pinned cap.

    Intended for the completed Hessian allocation; never repeats the old Hessian.
    """
    folder = folder.resolve()
    manifest = folder / 'screwer_prepared.json'
    if not manifest.is_file():
        raise RuntimeError(f'{folder}: not a prepared Screwer correction')
    info = json.loads(manifest.read_text())
    cycle = int(info.get('cycle', 1))
    maximum = int(info.get('max_cycles', 1))
    if not (1 <= cycle <= maximum <= 10):
        raise RuntimeError('invalid cycle bounds in Screwer manifest')
    if cycle >= maximum:
        print(f'{folder}: reached configured limit of {maximum} Screwer cycle(s)')
        return None
    metal = info['metal']
    base = Path(info.get('project_root') or folder.parents[1]).resolve()
    if cycle_directory(base, metal, cycle) != folder:
        raise RuntimeError('Screwer cycle directory does not match recorded project root')
    next_folder = cycle_directory(base, metal, cycle + 1)
    if next_folder.exists():
        raise RuntimeError(f'{next_folder}: next correction already exists; refusing duplicate submission')
    method = info['method']
    # Do not re-run an incomplete Hessian, or mistake an inherited spectrum for a fresh result.
    verify_parent(folder, method)
    negatives = [(mode, val) for mode, val in spectrum_modes(folder / 'vibspectrum') if val < 0]
    threshold = float(info.get('min_negative_cm1', 0))
    selected = [(mode, val) for mode, val in negatives if val < -threshold]
    if not selected:
        print(f'{folder}: no modes below -{threshold:g} cm^-1; no additional cycle')
        return None
    next_child = prepare(base, metal, method, min_negative=threshold,
                         define_backend=info.get('define_backend', 'legacy'),
                         cycle=cycle + 1, max_cycles=maximum)
    if next_child is not None:
        submit(next_child, method)
    return next_child


def main():
    parser = argparse.ArgumentParser(
        description='Opt-in Screwer correction; optionally chain new cycles after verified Hessians.',
        epilog=('Examples: lan-screwer num --metals gd --prepare-only; '
                'lan-screwer num --metals gd; lan-screwer --export-only --metals gd'))
    parser.add_argument('method', nargs='?', choices=['num', 'ao', 'numforce', 'aoforce'],
                        help='Hessian method from the completed original calculation')
    parser.add_argument('--base', type=Path, default=Path.cwd(), help='project directory')
    parser.add_argument('--metals', nargs='+', choices=METALS, help='metals to process')
    parser.add_argument('--cycle', type=int, default=1, help='correction cycle to prepare or submit (default 1)')
    parser.add_argument('--max-cycles', type=int, default=None, help='maximum automatically chained cycles (1-10); opt-in')
    parser.add_argument('--folder', type=Path, help=argparse.SUPPRESS)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--prepare-only', action='store_true', help='create screwer_<metal>/ inputs but do not submit')
    action.add_argument('--submit-prepared', action='store_true', help='submit an already prepared screwer_<metal>/ directory')
    action.add_argument('--export-only', action='store_true', help='extract $noproj data for completed screwer corrections')
    action.add_argument('--clean-frequency', action='store_true', help=argparse.SUPPRESS)
    action.add_argument('--define-fallback', action='store_true', help=argparse.SUPPRESS)
    action.add_argument('--advance-cycle', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--define-backend', choices=('legacy', 'adaptive'), default='legacy',
                        help='experimental adaptive pexpect supervisor; pinned at preparation, legacy is default')
    parser.add_argument('--min-imag', type=float, default=0.0,
                        help='only trigger on frequencies strictly below -VALUE cm^-1 (default 0)')
    args = parser.parse_args()
    if not (0 <= args.min_imag < float('inf')):
        parser.error('--min-imag must be a finite nonnegative number')
    if args.define_fallback:
        try:
            run_define_with_fallback(Path.cwd())
            return 0
        except (RuntimeError, OSError, ValueError) as error:
            print(f'ERROR Screwer DEFINE: {error}', file=sys.stderr)
            return 1
    if args.advance_cycle:
        if args.folder is None:
            parser.error('--advance-cycle requires --folder')
        try:
            advance_cycle(args.folder)
            return 0
        except (RuntimeError, OSError, ValueError, subprocess.TimeoutExpired) as error:
            print(f'ERROR cycle advance: {error}', file=sys.stderr)
            return 1
    if args.folder and (args.clean_frequency or args.export_only):
        try:
            if args.clean_frequency:
                clean_inherited_frequencies(args.folder.resolve())
            else:
                extract_noproj(args.folder.resolve())
            return 0
        except (RuntimeError, OSError, ValueError) as error:
            print(f'ERROR {args.folder}: {error}', file=sys.stderr)
            return 1
    if not args.metals:
        parser.error('--metals is required')
    if not (1 <= args.cycle <= 10):
        parser.error('--cycle must be 1..10')
    if args.max_cycles is not None and not (args.cycle <= args.max_cycles <= 10):
        parser.error('--max-cycles must be >= --cycle and <= 10')
    method = {'num': 'numforce', 'ao': 'aoforce', 'numforce': 'numforce', 'aoforce': 'aoforce'}.get(args.method)
    if not (args.export_only or args.clean_frequency) and not method:
        parser.error('choose num or ao (unless using --export-only)')
    base = args.base.resolve()
    maximum = args.max_cycles if args.max_cycles is not None else args.cycle
    failures = 0
    for metal in args.metals:
        try:
            if args.clean_frequency:
                clean_inherited_frequencies(cycle_directory(base, metal, args.cycle))
            elif args.export_only:
                new_child = cycle_directory(base, metal, args.cycle)
                extract_noproj(new_child if new_child.is_dir() else legacy_dir(base / metal))
            elif args.submit_prepared:
                submit(cycle_directory(base, metal, args.cycle), method)
            else:
                dest = prepare(base, metal, method, min_negative=args.min_imag,
                               define_backend=args.define_backend, cycle=args.cycle,
                               max_cycles=maximum)
                if dest is not None and not args.prepare_only:
                    submit(dest, method)
        except (RuntimeError, OSError, subprocess.TimeoutExpired, ValueError) as error:
            print(f'ERROR {metal}: {error}', file=sys.stderr)
            failures += 1
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
