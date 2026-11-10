#!/usr/bin/env python3
"""Read-only status summary for the dependent TURBOMOLE lanthanide chains."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

from setup_lanthanides import SPIN

parser = argparse.ArgumentParser(
    description='Show SLURM/job progress for all lanthanides; write workflow_status.csv.',
    epilog='Example: lan-status    (run from any EF project directory)',
)
parser.add_argument('--base', type=Path, default=Path.cwd(),
                    help='project directory (default: current working directory)')
args = parser.parse_args()
base = args.base.resolve()

queue = {}
try:
    res = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i|%T'],
                         text=True, capture_output=True, timeout=10)
    if res.returncode == 0:
        for line in res.stdout.splitlines():
            if '|' in line:
                jobid, status = line.split('|', 1)
                queue[jobid] = status
except (OSError, subprocess.TimeoutExpired):
    pass


def read_json(path):
    try:
        return json.loads(path.read_text()) if path.is_file() else {}
    except (ValueError, OSError):
        return {}


def stage_status(history, stage):
    for event in reversed(history):
        if event.get('stage') == stage:
            return event.get('status', '-')
    return '-'


def hole_status(directory, stage):
    data = read_json(directory / f'occupation_{stage}.json')
    if not data:
        return '-'
    status = data.get('status', 'UNKNOWN')
    spinholes = len(data.get('holes_within_each_spin', []))
    fractional = len(data.get('fractional_occupations', []))
    return f'{status} (h={spinholes}, f={fractional})'


rows = []
for metal, spin in SPIN.items():
    d = base / metal
    prepared = (d / 'metal.json').is_file()
    history = []
    hist_file = d / 'workflow_history.jsonl'
    if hist_file.is_file():
        for line in hist_file.read_text(errors='replace').splitlines():
            try:
                history.append(json.loads(line))
            except ValueError:
                pass
    state = read_json(d / 'workflow_status.json')
    jobs = []
    chain_file = d / 'chain_submission.tsv'
    if chain_file.is_file():
        for line in chain_file.read_text().splitlines()[1:]:
            parts = line.split('\t')
            if len(parts) >= 2:
                jobs.append((parts[0], parts[1]))
    running = [(stage, jid, queue[jid]) for stage, jid in jobs if jid in queue]
    if running:
        priority = {'RUNNING': 0, 'COMPLETING': 1, 'PENDING': 2}
        active_stage, jid, slurm = min(running, key=lambda x: priority.get(x[2], 3))
    else:
        active_stage, jid, slurm = ('-', jobs[-1][1], 'NOT_QUEUED') if jobs else ('-', '-', '-')
    if not d.exists():
        stage = 'NOT_PREPARED'
    elif not prepared:
        stage = 'EXISTING_PRESERVED'
    else:
        stage = state.get('stage', 'PREPARED')
    rows.append({
        'metal': metal, 'unpaired': spin, 'jobid': jid, 'slurm': slurm,
        'stage': stage, 'result': state.get('status', '-'),
        'scf': stage_status(history, 'RIDFT'),
        'integer_scf': stage_status(history, 'INTEGER_SCF'),
        'geometry': stage_status(history, 'JOBEX'),
        'frequency': stage_status(history, 'NUMFORCE') if any(s == 'numforce' for s, _ in jobs)
                     else stage_status(history, 'AOFORCE'),
        'frequency_method': jobs[-1][0] if jobs and jobs[-1][0] in ('numforce', 'aoforce') else '-',
        'holes_pre': hole_status(d, 'pre'),
        'holes_post': hole_status(d, 'post'),
        'submitted_stages': len(jobs), 'queue_stage': active_stage,
    })

outpath = base / 'workflow_status.csv'
with outpath.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f'{"METAL":<6} {"SLURM":<12} {"STAGE":<15} {"RESULT":<12} {"SCF":<8} {"GEO":<8} {"FREQ":<9} {"METHOD":<9} HOLES PRE / POST')
for r in rows:
    print(f'{r["metal"]:<6} {r["slurm"]:<12} {r["stage"]:<15} {r["result"]:<12} '
          f'{r["scf"]:<8} {r["geometry"]:<8} {r["frequency"]:<9} '
          f'{r["frequency_method"]:<9} {r["holes_pre"]} / {r["holes_post"]}')
print(f'\nWrote: {outpath}')
