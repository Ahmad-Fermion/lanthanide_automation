#!/usr/bin/env python3
"""Fail-closed checks for prepared inputs and incomplete SLURM submissions.

Only reads project/TURBOMOLE files. No scientific inputs are ever edited here.
"""
import argparse
import json
import getpass
import os
from pathlib import Path
import subprocess
import sys

from setup_lanthanides import SPIN, make_control, make_coord, make_define, read_charge

# setup_lanthanides.py creates precisely these files, before any calculations.
PREPARED_FILES = {
    'control', 'coord', 'def.inp', 'metal.json', 'tm_workflow.py',
    'turbo.pipeline.optimized.slurm', 'turbo.resume.slurm',
}


class UnsafeSubmission(Exception):
    pass


def verify_prepared(base: Path, metal: str, check_queue: bool = False):
    """Only permit a pristine lan-setup directory under --include-existing.

    Compare inputs to the *current* master templates. If masters changed after
    preparation, fail closed rather than silently accepting an edited control.
    """
    directory = base / metal
    if not directory.is_dir() or directory.is_symlink():
        raise UnsafeSubmission('not a real prepared metal directory')
    children = {item.name for item in directory.iterdir()}
    # lan-run holds this temporary directory lock while doing the verification.
    unexpected = sorted(children - PREPARED_FILES - {'.lan_submission_lock'})
    missing = sorted(PREPARED_FILES - children)
    if unexpected or missing:
        raise UnsafeSubmission(
            f'folder is not pristine (extra files: {unexpected or "none"}; '
            f'missing files: {missing or "none"}). '
            'Existing calculation outputs or unrecognized files require manual review.'
        )
    for name in PREPARED_FILES:
        p = directory / name
        if p.is_symlink() or not p.is_file():
            raise UnsafeSubmission(f'{name} must be a regular file, not a symlink or directory')
    source = base / ('def_gd.inp' if metal == 'gd' else 'def.inp')
    if not source.is_file():
        raise UnsafeSubmission(f'master {source.name} unavailable for verifying the original preparation')
    try:
        expected = {
            'control': make_control((base / 'control').read_text(), SPIN[metal]),
            'coord': make_coord((base / 'coord').read_text(), metal),
            'def.inp': make_define(source.read_text(), metal, SPIN[metal]),
        }
        for name, text in expected.items():
            if (directory / name).read_text() != text:
                raise UnsafeSubmission(
                    f'{name} differs from freshly prepared inputs. '
                    'It may have been edited or the project master changed; refusing to run define.'
                )
        meta = json.loads((directory / 'metal.json').read_text())
        if meta != {'metal': metal, 'unpaired': SPIN[metal],
                    'charge': read_charge(expected['def.inp'])}:
            raise UnsafeSubmission('metal.json does not match prepared metal/charge/spin')
    except (OSError, UnicodeError, ValueError, KeyError) as exc:
        raise UnsafeSubmission(f'could not verify prepared inputs: {exc}') from exc
    if check_queue:
        # A manually queued job may have no output files or chain record yet.
        # Do not submit DEFINE to a directory with an existing queued job.
        queued = run_command([
            'squeue', '-h', '-u', os.environ.get('USER') or getpass.getuser(),
            '-o', '%i|%Z'
        ])
        for line in queued.splitlines():
            fields = line.split('|', 1)
            if len(fields) == 2 and fields[1] and Path(fields[1]).resolve() == directory:
                raise UnsafeSubmission(f'active or pending job {fields[0]} already uses {directory}')
    print('PREPARED')


def run_command(args):
    try:
        completed = subprocess.run(args, check=False, capture_output=True,
                                   text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise UnsafeSubmission(f'cannot verify scheduler state ({args[0]}): {exc}') from exc
    if completed.returncode:
        raise UnsafeSubmission(f'{args[0]} failed while checking previous jobs: {completed.stderr.strip()}')
    return completed.stdout


def check_job_states(jids, directory: Path):
    """Check recorded job IDs using one squeue call and at most one sacct call.

    A job that has left squeue is not assumed successful: its accounting record
    must confirm COMPLETED/0:0 and the expected work directory.
    """
    queued = run_command([
        'squeue', '-h', '-u', os.environ.get('USER') or getpass.getuser(),
        '-o', '%i|%T|%Z|%r'
    ])
    active = {}
    requested = set(jids)
    for line in queued.splitlines():
        fields = line.split('|', 3)
        if len(fields) != 4 or fields[0] not in requested:
            continue
        jid, state, workdir, reason = fields
        if not workdir or Path(workdir).resolve() != directory.resolve():
            raise UnsafeSubmission(f'job {jid} belongs to {workdir}, not {directory}')
        if state not in {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING',
                         'SUSPENDED', 'RESIZING', 'STAGE_OUT'}:
            raise UnsafeSubmission(f'job {jid} has unsafe queue state {state}')
        if 'DependencyNeverSatisfied' in reason:
            raise UnsafeSubmission(f'job {jid} has an unsatisfiable dependency')
        active[jid] = 'ACTIVE'

    missing = [jid for jid in jids if jid not in active]
    if not missing:
        return active
    # Query completed/failed jobs once, instead of sending repeated RPCs to
    # the SLURM accounting daemon. Exact JobIDRaw excludes .batch/step rows.
    accounted = run_command([
        'sacct', '-n', '-P', '-j', ','.join(missing),
        '--format=JobIDRaw,State,ExitCode,WorkDir'
    ])
    for line in accounted.splitlines():
        fields = line.split('|')
        if len(fields) < 4 or fields[0] not in missing:
            continue
        jid, state, exit_code, workdir = fields[:4]
        if jid in active:
            raise UnsafeSubmission(f'multiple accounting records for job {jid}; manual review required')
        if not workdir or Path(workdir).resolve() != directory.resolve():
            raise UnsafeSubmission(f'cannot verify workdir of completed job {jid}: {workdir!r}')
        parsed_state = state.split()[0].split('+')[0] if state.strip() else 'UNKNOWN'
        if parsed_state != 'COMPLETED' or exit_code != '0:0':
            raise UnsafeSubmission(
                f'job {jid} ended in state {state} with exit code {exit_code}; '
                'use lan-recover if appropriate'
            )
        active[jid] = 'COMPLETED'
    for jid in missing:
        if jid not in active:
            raise UnsafeSubmission(f'job {jid} absent from squeue and sacct; unable to verify success')
    return active


def verify_resume(directory: Path, stages):
    record = directory / 'chain_submission.tsv'
    plan = directory / 'chain_plan.txt'
    pending = directory / '.lan_submission_pending'
    if pending.exists():
        raise UnsafeSubmission(
            f'{pending.name} exists: submission may have been interrupted after sbatch '
            'accepted a job. Inspect squeue/sacct and reconcile manually; do NOT delete blindly.'
        )
    if not record.is_file():
        raise UnsafeSubmission('chain_submission.tsv missing; no submitted stage to resume')
    if not plan.is_file():
        raise UnsafeSubmission(
            'chain_plan.txt missing (legacy submission). Original Hessian choice cannot '
            'be verified; inspect this partial chain manually rather than guessing.'
        )
    saved_stages = plan.read_text().splitlines()
    if saved_stages != stages:
        raise UnsafeSubmission(f'original chain plan {saved_stages} differs from requested {stages}')
    lines = record.read_text().splitlines()
    if not lines or lines[0] != 'stage\tjob_id\tpartition\twalltime\tmemory':
        raise UnsafeSubmission('chain_submission.tsv has an unrecognized header')
    rows = [line.split('\t') for line in lines[1:]]
    if not rows or any(len(row) != 5 for row in rows):
        raise UnsafeSubmission('chain record is empty/malformed; manual reconciliation required')
    if len(rows) > len(stages) or [r[0] for r in rows] != stages[:len(rows)]:
        raise UnsafeSubmission('recorded stages are not a valid prefix of the requested chain')
    seen = set()
    for row in rows:
        jid = row[1]
        if not jid.isdecimal() or jid in seen:
            raise UnsafeSubmission(f'bad or repeated SLURM job ID {jid!r}')
        seen.add(jid)
    if len(rows) == len(stages):
        # Nothing can be resubmitted: a complete submission record is enough
        # to refuse duplicates, without requiring historical accounting data.
        print('FULL')
        return
    states = check_job_states([row[1] for row in rows], directory)
    ordered = [states[row[1]] for row in rows]
    if any(earlier == 'ACTIVE' and later == 'COMPLETED'
           for i, earlier in enumerate(ordered) for later in ordered[i + 1:]):
        raise UnsafeSubmission('inconsistent SLURM states: a later job completed while an earlier dependency is active')
    last_state = ordered[-1]
    # afterok is necessary when the last recorded job is still active.
    # When the job is already COMPLETED, submitting with a historical afterok
    # dependency may never release on SLURM after MinJobAge has expired.
    last_dep = rows[-1][1] if last_state == 'ACTIVE' else '-'
    print(f'PARTIAL\t{len(rows)}\t{last_dep}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='operation', required=True)
    prepared = sub.add_parser('prepared')
    prepared.add_argument('--base', type=Path, required=True)
    prepared.add_argument('--metal', choices=tuple(SPIN), required=True)
    prepared.add_argument('--check-queue', action='store_true')
    resumed = sub.add_parser('resume')
    resumed.add_argument('--dir', type=Path, required=True)
    resumed.add_argument('--stages', nargs='+', required=True)
    args = p.parse_args()
    try:
        if args.operation == 'prepared':
            verify_prepared(args.base.resolve(), args.metal, args.check_queue)
        else:
            verify_resume(args.dir.resolve(), args.stages)
    except UnsafeSubmission as exc:
        print(f'UNSAFE: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
