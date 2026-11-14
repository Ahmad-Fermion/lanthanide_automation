#!/usr/bin/env python3
"""Summarize live SLURM jobs and recorded TURBOMOLE results separately."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess

from setup_lanthanides import SPIN


RUNNING_STATES = {'RUNNING', 'COMPLETING', 'CONFIGURING', 'RESIZING', 'STAGE_OUT',
                  'SIGNALING', 'SUSPENDED', 'STOPPED'}
PENDING_STATES = {'PENDING', 'REQUEUED', 'REQUEUE_HOLD', 'REQUEUE_FED', 'SPECIAL_EXIT'}
FAILURE_STATES = {'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL',
                  'BOOT_FAIL', 'DEADLINE', 'PREEMPTED', 'REVOKED'}


def read_json(path):
    try:
        parsed = json.loads(path.read_text()) if path.is_file() else {}
        return parsed if isinstance(parsed, dict) else {}
    except (ValueError, OSError):
        return {}


def read_history(path):
    if not path.is_file():
        return []
    try:
        lines = path.read_text(errors='replace').splitlines()
    except OSError:
        return []
    entries = []
    for line in lines:
        try:
            obj = json.loads(line)
            if isinstance(obj, dict):
                entries.append(obj)
        except ValueError:
            continue
    return entries


def recorded_result(status):
    """A historical RUNNING/RETRY record is not proof a job is still executing."""
    return 'INCOMPLETE' if status in ('RUNNING', 'RETRY') else (status or '-')


def stage_status(history, stage):
    for event in reversed(history):
        if event.get('stage') == stage:
            return recorded_result(event.get('status'))
    return '-'


def hole_status(directory, stage):
    data = read_json(directory / f'occupation_{stage}.json')
    if not data:
        return '-'
    status = data.get('status', 'UNKNOWN')
    if status == 'NO_HOLES_IN_PRINTED_WINDOW':
        status = 'NO_HOLES_SHOWN'  # Display only; retain the full diagnostic in JSON.
    spinholes = len(data.get('holes_within_each_spin', []))
    fractional = len(data.get('fractional_occupations', []))
    return f'{status} (h={spinholes}, f={fractional})'


def read_submission(path):
    """Return the recorded job stages; flag malformed records rather than guessing."""
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return [], 'unreadable submission record'
    if not lines or not lines[0].startswith('stage\tjob_id\t'):
        return [], 'invalid submission record header'
    jobs = []
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) < 2 or not fields[1].isdecimal():
            return jobs, 'malformed submission record'
        jobs.append((fields[0], fields[1]))
    if not jobs:
        return [], 'empty submission record'
    if len({jid for _, jid in jobs}) != len(jobs):
        return jobs, 'repeated job IDs in submission record'
    return jobs, ''


def latest_attempt(directory):
    """Prefer the most recently written recovery/normal submission record."""
    paths = ([directory / 'chain_submission.tsv'] +
             list(directory.glob('recovery_submission.*.tsv')))
    paths = [p for p in paths if p.is_file()]
    if not paths:
        return None
    return max(paths, key=lambda p: (p.stat().st_mtime_ns, p.name))


def run_scheduler(command):
    try:
        res = subprocess.run(command, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, f'{command[0]} unavailable: {exc}'
    if res.returncode:
        return None, f'{command[0]} failed: {res.stderr.strip() or res.returncode}'
    return res.stdout, ''


def normalize_slurm_state(raw):
    return raw.strip().split(' ', 1)[0].split('+', 1)[0].upper() if raw.strip() else 'UNKNOWN'


def matches_directory(raw, expected):
    if not raw or raw in ('Unknown', 'N/A', 'None'):
        return False
    try:
        return Path(raw).resolve() == expected.resolve()
    except (OSError, ValueError):
        return False


def scheduler_snapshot(attempts):
    """One squeue query and one batched sacct query for the entire project."""
    ids = {jid for attempt in attempts.values() if attempt
           for _, jid in attempt['jobs']}
    user = os.environ.get('USER') or os.environ.get('LOGNAME') or ''
    queue_text, queue_error = run_scheduler(
        ['squeue', '-h', '-u', user, '-o', '%i|%T|%Z'])
    queue = {}
    if queue_text is not None:
        for line in queue_text.splitlines():
            fields = line.split('|', 2)
            if len(fields) == 3 and fields[0] in ids:
                queue[fields[0]] = (normalize_slurm_state(fields[1]), fields[2])
    accounting = {}
    accounting_error = ''
    if ids:
        # -X returns job allocations, not the .batch/.extern task-step records.
        accounting_text, accounting_error = run_scheduler(
            ['sacct', '-n', '-P', '-X', '-j', ','.join(sorted(ids)),
             '--format=JobIDRaw,State,ExitCode,WorkDir'])
        if accounting_text is not None:
            for line in accounting_text.splitlines():
                fields = line.split('|')
                if len(fields) < 4 or fields[0] not in ids:
                    continue
                jid = fields[0]
                record = (normalize_slurm_state(fields[1]), fields[2], fields[3])
                # Multiple distinct allocation rows for the same job ID are ambiguous.
                if jid in accounting and accounting[jid] != record:
                    accounting[jid] = ('AMBIGUOUS', '-', '')
                else:
                    accounting[jid] = record
    return queue, queue_error, accounting, accounting_error


def stage_label(stage):
    return {
        'ridft_eiger_fix': 'RIDFT', 'ridft': 'RIDFT',
        'ridft-recovery': 'RIDFT', 'eiger_pre': 'EIGER',
        'occupation_fix': 'OCCUPATIONS', 'integer_ridft': 'INTEGER_SCF',
        'eiger_fix_integer_ridft': 'INTEGER_SCF',
        'jobex_eiger': 'JOBEX', 'jobex': 'JOBEX', 'eiger_post': 'EIGER',
        'numforce': 'FREQ', 'aoforce': 'FREQ',
    }.get(stage, stage.upper() if stage else '-')


def event_is_current(timestamp, attempt):
    """Check timestamps without treating a previous attempt as current."""
    if not attempt or not timestamp:
        return False
    try:
        stamp = datetime.fromisoformat(timestamp)
        if stamp.tzinfo is None:
            return False
        # Recorded workflow timestamps have whole-second resolution.
        return stamp.timestamp() + 1 >= attempt['path'].stat().st_mtime
    except (ValueError, OSError, TypeError):
        return False


def final_event_for_attempt(history, attempt):
    """Require a current DONE/SUCCESS event, not an old success marker."""
    for event in reversed(history):
        if event.get('stage') == 'DONE' and event.get('status') == 'SUCCESS':
            return event_is_current(event.get('updated_utc'), attempt)
    return False


def classify(directory, prepared, attempt, history, state, queue, queue_error,
             accounting, accounting_error):
    if not directory.exists():
        return 'NOT_PREPARED', '-', 'INACTIVE' if not queue_error else 'UNKNOWN', '-', '-', 'No directory'
    if not prepared:
        return 'NEEDS_REVIEW', '-', 'INACTIVE' if not queue_error else 'UNKNOWN', '-', '-', 'Existing unrecognized directory'
    if attempt is None:
        untracked = bool(history or (state and state.get('stage') not in ('PREPARED',))
                         or (directory / 'define_setup.ok').exists())
        return ('NEEDS_REVIEW' if untracked else 'PREPARED', '-',
                'INACTIVE' if not queue_error else 'UNKNOWN', '-', '-',
                'No tracked job; calculation evidence exists' if untracked else 'Not submitted')

    jobs = attempt['jobs']
    note = []
    blockers = []

    def require_review(reason):
        note.append(reason)
        blockers.append(reason)

    if attempt['error']:
        require_review(attempt['error'])
    if (directory / '.lan_submission_pending').exists():
        require_review('ambiguous pending submission marker')
    if (directory / '.lan_submission_lock').exists():
        require_review('submission lock present')
    if attempt['path'].name == 'chain_submission.tsv':
        plan_file = directory / 'chain_plan.txt'
        if plan_file.is_file():
            try:
                plan = plan_file.read_text().splitlines()
                if not plan or [s for s, _ in jobs] != plan[:len(jobs)] or len(jobs) > len(plan):
                    require_review('submission stages do not match saved chain plan')
                elif len(jobs) < len(plan):
                    require_review(f'partial chain: {len(jobs)}/{len(plan)} jobs submitted')
            except OSError:
                require_review('unreadable chain plan')
        else:
            require_review('legacy chain: intended final stage unverified')

    live = []
    historical = []
    uncertain = []
    failed = []
    for stage, jid in jobs:
        if jid in queue:
            status, workdir = queue[jid]
            if not matches_directory(workdir, directory):
                uncertain.append(f'job {jid} has mismatched/unknown workdir in squeue')
            elif status in RUNNING_STATES | PENDING_STATES:
                live.append((stage, jid, status))
            elif status in FAILURE_STATES:
                failed.append((stage, jid, status))
            else:
                uncertain.append(f'job {jid}: unexpected live state {status}')
        elif queue_error:
            uncertain.append('live queue unavailable')
        else:
            if jid not in accounting:
                uncertain.append(f'no accounting record for job {jid}')
            else:
                status, exit_code, workdir = accounting[jid]
                if not matches_directory(workdir, directory):
                    uncertain.append(f'job {jid} has mismatched/unknown workdir in sacct')
                elif status == 'COMPLETED' and exit_code == '0:0':
                    historical.append((stage, jid, status))
                elif status in FAILURE_STATES or (status == 'COMPLETED' and exit_code != '0:0'):
                    failed.append((stage, jid, status if status != 'COMPLETED' else 'FAILED'))
                else:
                    uncertain.append(f'job {jid}: accounting state {status} / {exit_code}')

    if accounting_error and any(jid not in queue for _, jid in jobs):
        note.append(accounting_error)
    if queue_error:
        require_review(queue_error)
    note.extend(uncertain)
    blockers.extend(issue for issue in uncertain
                    if not issue.startswith('no accounting record'))
    failure = failed[0] if failed else None
    current = next(((stage, jid, status) for stage, jid, status in live
                    if status in RUNNING_STATES), None)
    if current is None:
        current = next(((stage, jid, status) for stage, jid, status in live
                        if status in PENDING_STATES), None)
    selected = current or failure or (jobs[-1] + ('-',) if jobs else ('-', '-', '-'))
    last_stage = stage_label(selected[0])
    jobid = selected[1]
    slurm = selected[2] if current else ('UNKNOWN' if queue_error else 'INACTIVE')
    # Prefer the actual failed stage for explaining unsuccessful chains.
    final_slurm = failure[2] if failure else (
        'COMPLETED' if historical and len(historical) == len(jobs) else '-')

    if failure:
        overall = 'FAILED'
        note.insert(0, f'job {failure[1]} ended {failure[2]}')
        last_stage = stage_label(failure[0])
        jobid = failure[1]
    elif current and not blockers:
        # A live job remains RUNNING/PENDING when sacct has not yet recorded
        # predecessors; accounting gaps remain visible in the notes and CSV.
        overall = 'RUNNING' if current[2] in RUNNING_STATES else 'PENDING'
    elif note or uncertain:
        overall = 'NEEDS_REVIEW'
    elif len(historical) == len(jobs) and final_event_for_attempt(history, attempt):
        overall = 'COMPLETED'
        # The grouped JOBEX allocation ends with a real post-JOBEX EIGER check.
        if jobs and jobs[-1][0] in ('jobex_eiger', 'eiger_post'):
            last_stage = 'EIGER'
    else:
        overall = 'NEEDS_REVIEW'
        note.append('no verified final DONE/SUCCESS for this attempt')
    return overall, last_stage, slurm, final_slurm, jobid, '; '.join(dict.fromkeys(note)) or '-'


def main():
    parser = argparse.ArgumentParser(
        description='Show overall status and scheduler history for each metal; write workflow_status.csv.',
        epilog='Examples: lan-status    |    lan-status --details')
    parser.add_argument('--base', type=Path, default=Path.cwd(),
                        help='project directory (default: current working directory)')
    parser.add_argument('--details', action='store_true',
                        help='also print SCF/geometry/frequency and orbital occupation details')
    args = parser.parse_args()
    base = args.base.resolve()

    attempts = {}
    for metal in SPIN:
        directory = base / metal
        path = latest_attempt(directory) if directory.is_dir() else None
        if path:
            jobs, error = read_submission(path)
            attempts[metal] = {'path': path, 'jobs': jobs, 'error': error}
        else:
            attempts[metal] = None
    queue, queue_error, accounting, accounting_error = scheduler_snapshot(attempts)

    rows = []
    for metal, spin in SPIN.items():
        directory = base / metal
        prepared = (directory / 'metal.json').is_file()
        history = read_history(directory / 'workflow_history.jsonl')
        state = read_json(directory / 'workflow_status.json')
        attempt = attempts[metal]
        jobs = attempt['jobs'] if attempt else []
        overall, last_stage, slurm, slurm_final, jid, reason = classify(
            directory, prepared, attempt, history, state, queue, queue_error,
            accounting, accounting_error)
        if not directory.exists():
            stage = 'NOT_PREPARED'
            result = '-'
        elif not prepared:
            stage = 'EXISTING_PRESERVED'
            result = recorded_result(state.get('status'))
        elif attempt and not event_is_current(state.get('updated_utc'), attempt):
            # Never display a previous SUCCESS as the present run's result.
            stage = 'NO_CURRENT_RECORD'
            result = 'INCOMPLETE'
        else:
            stage = state.get('stage', 'PREPARED')
            result = recorded_result(state.get('status'))
        running = [(stage_name, job_id, queue[job_id][0]) for stage_name, job_id in jobs
                   if job_id in queue]
        active = next((r for r in running if r[2] in RUNNING_STATES), None)
        if active is None:
            active = next((r for r in running if r[2] in PENDING_STATES), None)
        frequency_method = next((s for s, _ in reversed(jobs)
                                 if s in ('numforce', 'aoforce')), '-')
        rows.append({
            'metal': metal, 'overall': overall, 'last_stage': last_stage,
            'unpaired': spin, 'jobid': jid, 'slurm': slurm,
            'slurm_final': slurm_final, 'stage': stage,
            'result': result,
            'scf': stage_status(history, 'RIDFT'),
            'integer_scf': stage_status(history, 'INTEGER_SCF'),
            'geometry': stage_status(history, 'JOBEX'),
            'frequency': (stage_status(history, 'NUMFORCE') if frequency_method == 'numforce'
                          else stage_status(history, 'AOFORCE') if frequency_method == 'aoforce'
                          else '-'),
            'frequency_method': frequency_method,
            'holes_pre': hole_status(directory, 'pre'),
            'holes_post': hole_status(directory, 'post'),
            'submitted_stages': len(jobs),
            'queue_stage': active[0] if active else '-',
            'attempt': attempt['path'].name if attempt else '-',
            'status_note': reason,
        })

    outpath = base / 'workflow_status.csv'
    with outpath.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'{"METAL":<6} {"OVERALL":<15} {"LAST STAGE":<13} {"LIVE SLURM":<14} {"SLURM FINAL":<13}')
    for row in rows:
        print(f'{row["metal"]:<6} {row["overall"]:<15} {row["last_stage"]:<13} '
              f'{row["slurm"]:<14} {row["slurm_final"]:<13}')
    if args.details:
        print(f'\n{"METAL":<6} {"STAGE":<15} {"RESULT":<12} {"SCF":<12} '
              f'{"INT SCF":<12} {"GEO":<12} {"FREQ":<12} {"METHOD":<9} HOLES PRE / POST')
        for row in rows:
            print(f'{row["metal"]:<6} {row["stage"]:<15} {row["result"]:<12} '
                  f'{row["scf"]:<12} {row["integer_scf"]:<12} {row["geometry"]:<12} {row["frequency"]:<12} '
                  f'{row["frequency_method"]:<9} {row["holes_pre"]} / {row["holes_post"]}')
        print('Stage results without a current record show INCOMPLETE; per-stage history is retained.')
    needs_review = [row for row in rows if row['overall'] in ('FAILED', 'NEEDS_REVIEW')]
    if needs_review:
        print('\nAttention:')
        for row in needs_review:
            print(f'  {row["metal"]}: {row["overall"]} — {row["status_note"]}')
    print(f'\nWrote: {outpath}')


if __name__ == '__main__':
    main()
