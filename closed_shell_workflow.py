#!/usr/bin/env python3
"""Opt-in restricted La/Lu input preparation and SLURM submission.

A user-authored def_closedshell.inp is required; automatic conversion of a UHF
EHT dialogue to a closed-shell dialogue is not scientifically dependable.
This never overwrites existing paths and never changes the project master.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from setup_lanthanides import ROOT, make_coord, read_charge

SUPPORTED = ('la', 'lu')
# Atomic numbers for the elements ordinarily encountered in these models;
# reject unknown atoms rather than assume even electron parity.
SYMBOLS = 'h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni cu zn ga ge as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs ba la ce pr nd pm sm eu gd tb dy ho er tm yb lu'.split()
ATOMIC_Z = {x: i + 1 for i, x in enumerate(SYMBOLS[:56])}
ATOMIC_Z.update({'la': 57, 'lu': 71})


def sha(s):
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


def restricted_control(source):
    # Remove *entire groups*, not just their first lines. All master DFT,
    # ECP/X2CSR, RI, dispersion and basis references are otherwise untouched.
    remove = {'uhf', 'uhfmo_alpha', 'uhfmo_beta', 'alpha shells', 'beta shells',
              'closed shells', 'scfmo', 'fermi'}
    marker = re.compile(r'(?m)^\s*\$([A-Za-z][^\n]*?)(?=\s*(?:\n|$))')
    # Group header with optional arguments, e.g. "$alpha shells".
    lines = source.splitlines(keepends=True)
    out = []
    dropping = False
    for line in lines:
        m = re.match(r'^\s*\$([\w-]+)(?:\s+(shells))?\b', line, re.I)
        if m:
            group = m.group(1).lower() + (' shells' if m.group(2) else '')
            dropping = group in remove
        if not dropping:
            out.append(line)
    changed = ''.join(out)
    if not re.search(r'(?m)^\s*\$end\s*$', changed):
        raise ValueError('master control has no $end')
    if re.search(r'(?m)^\s*\$(?:uhf\b|fermi\b|alpha\s+shells\b|beta\s+shells\b)', changed, re.I):
        raise ValueError('restricted control still contains open-shell initialization')
    return changed


def restricted_define(source, metal):
    if re.search(r'(?im)^\s*u\s+\d+\s*$', source):
        raise ValueError('def_closedshell.inp contains a u N (unrestricted occupation) command')
    updated, count = re.subn(r'(?im)^(\s*b\s+)["\'](?:la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu)["\'](\s+\S+\s*)$',
                             lambda m: m.group(1) + '"' + metal + '"' + m.group(2), source)
    if count != 1:
        raise ValueError(f'def_closedshell.inp must contain exactly one quoted Ln basis selector (found {count})')
    if read_charge(updated) is None:
        raise ValueError('cannot read molecular charge after eht in def_closedshell.inp')
    return updated


def electron_count(coord, charge):
    atoms = []
    active = False
    for line in coord.splitlines():
        stripped = line.strip()
        if stripped.startswith('$coord'):
            active = True
            continue
        if active and stripped.startswith('$'):
            break
        if active and stripped and not stripped.startswith('#'):
            fields = stripped.split()
            if len(fields) < 4:
                raise ValueError('malformed coord atom row')
            symbol = fields[3].lower()
            if symbol not in ATOMIC_Z:
                raise ValueError(f'unknown atomic number for {symbol}; cannot verify electron parity')
            atoms.append(symbol)
    if not atoms:
        raise ValueError('coord has no atom rows')
    ne = sum(ATOMIC_Z[a] for a in atoms) - charge
    if ne <= 0 or ne % 2:
        raise ValueError(f'molecular electron count {ne} is not a positive even integer; restricted closed-shell setup refused')
    return ne


def inputs(base, metal):
    for name in ('control', 'coord', 'def_closedshell.inp'):
        if not (base / name).is_file():
            raise ValueError(f'missing project input {base / name}')
    coord = make_coord((base / 'coord').read_text(), metal)
    define = restricted_define((base / 'def_closedshell.inp').read_text(), metal)
    charge = read_charge(define)
    count = electron_count(coord, charge)
    control = restricted_control((base / 'control').read_text())
    return {'coord': coord, 'control': control, 'def.inp': define}, charge, count


def prepare(base, metal):
    folder = base / f'{metal}_closedshell'
    if folder.exists():
        print(f'SKIP {folder.name}: existing path preserved')
        return False
    files, charge, electrons = inputs(base, metal)
    folder.mkdir()
    try:
        for name, content in files.items():
            (folder / name).write_text(content)
        for name in ('tm_workflow.py', 'turbo.pipeline.slurm', 'turbo.resume.slurm',
                     'setup_lanthanides.py', 'screwer_workflow.py'):
            shutil.copy2(ROOT / name, folder / name)
        meta = {'metal': metal, 'unpaired': 0, 'charge': charge, 'restricted': True,
                'electron_count_all_electron': electrons}
        (folder / 'metal.json').write_text(json.dumps(meta, indent=2) + '\n')
        manifest = {'workflow_version': 1, 'metal': metal, 'variant': 'restricted',
                    'files_sha256': {name: sha(content) for name, content in files.items()},
                    'metal_sha256': sha((folder / 'metal.json').read_text()),
                    'scientific_master': str(base), 'electrons_all_electron': electrons}
        (folder / 'restricted_prepared.json').write_text(json.dumps(manifest, indent=2) + '\n')
    except BaseException:
        shutil.rmtree(folder)
        raise
    print(f'CREATED {folder.name}: restricted, charge={charge}, electrons={electrons}')
    return True


def pristine(folder, base, metal):
    """Validate an unsubmitted restricted folder without forbidding DEFINE edits.

    The prepared manifest is provenance, not a lock on the user's scientific
    DEFINE transcript. We still reject missing files, wrong metal/charge,
    unrestricted initialization, and folders that have already been run.
    """
    allowed = {'coord', 'control', 'def.inp', 'metal.json', 'restricted_prepared.json',
               'tm_workflow.py', 'turbo.pipeline.slurm', 'turbo.resume.slurm',
               'setup_lanthanides.py', 'screwer_workflow.py'}
    if folder.is_symlink() or not folder.is_dir() or set(p.name for p in folder.iterdir()) != allowed:
        raise ValueError(f'{folder} is not a pristine prepared restricted folder')
    info = json.loads((folder / 'restricted_prepared.json').read_text())
    if info.get('variant') != 'restricted' or info.get('metal') != metal:
        raise ValueError('invalid restricted provenance')

    # Keep the existing guards for prepared geometry and control. These are not
    # the manually editable DEFINE transcript, and drift may invalidate the
    # restricted initial electronic state.
    sources = {'coord': make_coord((base / 'coord').read_text(), metal),
               'control': restricted_control((base / 'control').read_text())}
    for filename, content in sources.items():
        p = folder / filename
        if p.is_symlink() or not p.is_file() or p.read_text() != content:
            raise ValueError(f'{folder}/{filename} was changed after preparation')
        if info['files_sha256'].get(filename) != sha(content):
            raise ValueError(f'{filename} manifest hash mismatch')

    meta_path = folder / 'metal.json'
    if meta_path.is_symlink() or not meta_path.is_file():
        raise ValueError('missing or linked metal.json')
    if info.get('metal_sha256') != sha(meta_path.read_text()):
        raise ValueError('metal.json altered after preparation')
    meta = json.loads(meta_path.read_text())
    if meta.get('metal') != metal or meta.get('restricted') is not True or meta.get('unpaired') != 0:
        raise ValueError('invalid restricted metal metadata')

    define_path = folder / 'def.inp'
    if define_path.is_symlink() or not define_path.is_file():
        raise ValueError('restricted def.inp missing or linked')
    define_text = define_path.read_text()
    if info.get('files_sha256', {}).get('def.inp') != sha(define_text):
        print(f'WARNING {folder}: def.inp was edited after preparation; preserving user-defined DEFINE input',
              flush=True)

    # Validate only the essentials instead of demanding byte-for-byte equality
    # with the master transcript. Dialogue, basis details and other manual
    # changes remain under the user's control.
    if re.search(r'(?im)^\s*u\s+\d+\s*$', define_text):
        raise ValueError('edited restricted def.inp contains a u N unrestricted occupation command')
    selectors = re.findall(
        r'(?im)^\s*b\s+["\'](la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu)["\']\s+\S+\s*$',
        define_text)
    if [symbol.lower() for symbol in selectors] != [metal]:
        raise ValueError(f'edited restricted def.inp must specify exactly one basis for {metal}')
    charge = read_charge(define_text)
    if charge is None or charge != meta.get('charge'):
        raise ValueError('edited restricted def.inp EHT charge disagrees with metal.json; reconcile charge before submission')


def submit(base, metal, method):
    folder = base / f'{metal}_closedshell'
    pristine(folder, base, metal)
    if not shutil.which('sbatch') or not shutil.which('squeue'):
        raise ValueError('SLURM submission/status commands missing')
    query = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i|%Z'],
                           capture_output=True, text=True)
    if query.returncode:
        raise ValueError('could not verify scheduler status')
    if any(x.split('|', 1)[-1].strip() == str(folder) for x in query.stdout.splitlines() if '|' in x):
        raise ValueError('existing active/pending SLURM job in prepared directory')
    stages = ['restricted_ridft', 'jobex_eiger']
    if method not in ('aoforce', 'numforce', 'none'):
        raise ValueError('invalid LAN_CLOSED_FREQ')
    if method != 'none':
        stages.append(method)
    # Maintain anti-double-submit guard even if an interrupted sbatch call is
    # accepted by SLURM before its job ID reaches the journal.
    if any((folder / x).exists() for x in ('chain_plan.txt', 'chain_submission.tsv', '.lan_submission_pending')):
        raise ValueError('already submitted/ambiguous; inspect scheduler; never submit twice')
    lock = folder / '.lan_submission_lock'
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise ValueError('another submission is in progress') from exc
    try:
        # Recheck after acquiring the local lock to avoid concurrent submissions.
        query = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i|%Z'],
                               capture_output=True, text=True)
        if query.returncode or any(x.split('|', 1)[-1].strip() == str(folder)
                                   for x in query.stdout.splitlines() if '|' in x):
            raise ValueError('scheduler failed or a job appeared while acquiring lock')
        dep = None
        for stage in stages:
            mem = os.environ.get('HESS_MEM' if stage in ('aoforce', 'numforce') else 'JOBEX_MEM' if stage == 'jobex_eiger' else 'RIDFT_MEM',
                                 '100G' if stage in ('aoforce', 'numforce') else '20G')
            threads = os.environ.get('TM_CPUS', '20')
            if not threads.isdigit() or int(threads) < 1:
                raise ValueError('invalid TM_CPUS')
            cmd = ['sbatch', '--parsable', '--partition=general', '--time=2-00:00:00',
                   f'--cpus-per-task={threads}', f'--mem={mem}',
                   f'--job-name=tm-{metal}-restricted-{stage}',
                   f'--export=ALL,TM_STAGE={stage}' + (',TM_FINAL_NO_HESSIAN=1' if stage == 'jobex_eiger' and method == 'none' else '')]
            if dep:
                cmd.append(f'--dependency=afterok:{dep}')
            cmd.append('turbo.pipeline.slurm')
            (folder / '.lan_submission_pending').write_text(f'stage={stage}\n')
            result = subprocess.run(cmd, cwd=folder, capture_output=True, text=True)
            if result.returncode:
                (folder / '.lan_submission_pending').unlink(missing_ok=True)
                raise ValueError(f'sbatch failed for {stage}: {result.stderr.strip()}')
            jid = result.stdout.split(';')[0].strip()
            if not jid.isdigit():
                raise ValueError('ambiguous sbatch response; inspect .lan_submission_pending')
            if dep is None:
                (folder / 'chain_plan.txt').write_text('\n'.join(stages) + '\n')
                (folder / 'chain_submission.tsv').write_text('stage\tjob_id\tpartition\twalltime\tmemory\n')
            with (folder / 'chain_submission.tsv').open('a') as f:
                f.write(f'{stage}\t{jid}\tgeneral\t2-00:00:00\t{mem}\n')
            (folder / '.lan_submission_pending').unlink()
            print(f'SUBMITTED {folder.name}: {stage} job {jid}' + (f' afterok:{dep}' if dep else ''))
            dep = jid
    finally:
        lock.rmdir()


def main(argv=None):
    parser = argparse.ArgumentParser(description='Opt-in restricted La/Lu DFT; creates separate *_closedshell directories')
    parser.add_argument('method', nargs='?', choices=['ao', 'num', 'aoforce', 'numforce'],
                        help='optional Hessian method (default geometry only)')
    parser.add_argument('--base', type=Path, default=Path.cwd())
    parser.add_argument('--metals', nargs='+', choices=list(SUPPORTED), required=True,
                        help='explicit La and/or Lu; never converts existing metal folders')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--include-existing', action='store_true')
    args = parser.parse_args(argv)
    base = args.base.resolve()
    made = []
    errors = 0
    method = {'ao': 'aoforce', 'num': 'numforce', 'aoforce': 'aoforce', 'numforce': 'numforce'}.get(args.method, 'none')
    for metal in args.metals:
        folder = base / f'{metal}_closedshell'
        was_existing = folder.exists()
        try:
            if not was_existing:
                prepare(base, metal)
                made.append(metal)
            if not args.prepare_only and (not was_existing or args.include_existing):
                submit(base, metal, method)
        except (ValueError, OSError, json.JSONDecodeError) as err:
            print(f'ERROR {metal}: {err}', file=sys.stderr)
            errors += 1
    return 1 if errors else 0

if __name__ == '__main__':
    sys.exit(main())
