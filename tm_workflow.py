#!/usr/bin/env python3
"""TURBOMOLE helper: validate DEFINE output, check SCF, inspect eiger, track stage."""
import argparse
import hashlib
import datetime as dt
import json
import os
from pathlib import Path
import re
import sys

HERE = Path.cwd()
META = HERE / 'metal.json'
STATE = HERE / 'workflow_status.json'
OCC_REFERENCE = HERE / 'initial_occupations'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def metadata():
    return json.loads(META.read_text())


def record(stage, status, message=''):
    obj = json.loads(STATE.read_text()) if STATE.exists() else {}
    obj.update({'metal': metadata()['metal'], 'stage': stage, 'status': status,
                'message': message, 'updated_utc': now()})
    tmp = STATE.with_name(STATE.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(obj, indent=2) + '\n')
    tmp.replace(STATE)
    with (HERE / 'workflow_history.jsonl').open('a') as log:
        log.write(json.dumps({'updated_utc': obj['updated_utc'], 'stage': stage, 'status': status, 'message': message}) + '\n')
    print(f'[{now()}] {stage}: {status} {message}', flush=True)


def parse_groups(text):
    """Return groups with original lines, respecting multi-line $dft etc."""
    out = []
    for line in text.splitlines(keepends=True):
        if re.match(r'^\s*\$[A-Za-z]', line):
            out.append([re.match(r'^\s*\$([\w-]+)', line).group(1).lower(), [line]])
        elif out:
            out[-1][1].append(line)
    return out



# Identify data-group boundaries by the next $keyword, not by $ricore:
# valid controls can contain $energy, $grad, and other groups before $ricore.
DATA_GROUP_RE = re.compile(r'(?m)^[ \t]*\$[^\r\n]*')
SHELL_OCC_RE = re.compile(r'\(\s*([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)\s*\)')


def group_spans(text):
    """Yield (name, start, end) preserving every byte outside each data group."""
    heads = list(DATA_GROUP_RE.finditer(text))
    for i, hit in enumerate(heads):
        header = hit.group().lstrip()[1:].strip().lower()
        if re.match(r'^alpha\s+shells(?:\s|$)', header):
            name = 'alpha'
        elif re.match(r'^beta\s+shells(?:\s|$)', header):
            name = 'beta'
        elif re.match(r'^fermi(?:\s|$)', header):
            name = 'fermi'
        else:
            name = None
        yield name, hit.start(), heads[i + 1].start() if i + 1 < len(heads) else len(text)


def shell_blocks(text, label):
    """Find exactly one alpha and beta shell block, irrespective of their lengths."""
    result = {}
    for name, start, end in group_spans(text):
        if name in ('alpha', 'beta'):
            if name in result:
                raise RuntimeError(f'{label}: multiple ${name} shells groups')
            result[name] = (start, end, text[start:end])
    if set(result) != {'alpha', 'beta'}:
        raise RuntimeError(f'{label}: expected one $alpha shells and one $beta shells group; found {sorted(result)}')
    return result


def save_initial_occupations():
    """Freeze the integer shell assignments produced by define, before first ridft."""
    groups = shell_blocks((HERE / 'control').read_text(), 'control after define')
    for name in ('alpha', 'beta'):
        lines = groups[name][2].splitlines()[1:]
        rows = [line for line in lines if line.strip()]
        if not rows:
            raise RuntimeError(f'define produced an empty ${name} shells group')
        for line in rows:
            match = SHELL_OCC_RE.search(line)
            if not match or abs(float(match.group(1)) - 1.0) > 1e-8:
                raise RuntimeError(f'define produced a non-integer or unrecognized ${name} shell: {line.strip()}')
    # Preserve the exact two blocks. No model- or metal-specific tables needed.
    reference = groups['alpha'][2] + groups['beta'][2]
    tmp = OCC_REFERENCE.with_name(OCC_REFERENCE.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(reference)
    tmp.replace(OCC_REFERENCE)
    print(f'Saved define occupations in {OCC_REFERENCE.name}', flush=True)


def restore_shells(text):
    """Replace only $alpha/$beta shell groups; preserve all other control groups."""
    if not OCC_REFERENCE.is_file():
        raise RuntimeError(f'{OCC_REFERENCE.name} missing; it must be saved immediately after define')
    expected = shell_blocks(OCC_REFERENCE.read_text(), OCC_REFERENCE.name)
    current = shell_blocks(text, 'current control')
    changed = any(current[name][2] != expected[name][2] for name in ('alpha', 'beta'))
    # Reverse-order editing preserves the offsets, even if beta precedes alpha.
    for name, (start, end, _) in sorted(current.items(), key=lambda item: item[1][0], reverse=True):
        text = text[:start] + expected[name][2] + text[end:]
    return text, changed


def restore_after_integer_scf():
    """Ensure step 4 leaves the original integer occupations for jobex."""
    control = HERE / 'control'
    before = control.read_text()
    after, changed = restore_shells(before)
    if changed:
        control.write_text(after)
        with (HERE / 'occupation_warnings.log').open('a') as log:
            log.write(f'{now()} integer_scf: restored initial $alpha/$beta shells after ridft\n')
    record('INTEGER_SCF', 'OK',
           'initial integer alpha/beta shells restored after ridft' if changed
           else 'initial integer alpha/beta shells confirmed after ridft')


def screwer_inheritance(directory):
    """Fingerprint parent MOs and shell groups; Screwer must reuse this state.

    This is recorded *after* generating the displaced coord, before DEFINE.
    It does not guarantee electronic-state correctness, only continuity.
    """
    orb_files = ['alpha', 'beta'] if all((directory / f).is_file() for f in ('alpha', 'beta')) else (
        ['mos'] if (directory / 'mos').is_file() else [])
    if not orb_files:
        raise RuntimeError(f'{directory}: inherited alpha/beta or mos orbital files missing; refusing screwer continuation')
    orbitals = {}
    for name in orb_files:
        checksum = hashlib.sha256()
        with (directory / name).open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                checksum.update(block)
        orbitals[name] = checksum.hexdigest()
    control = (directory / 'control').read_text()
    # Find whole shell blocks with explicit data-group boundaries.
    pattern = re.compile(r'(?ims)^[ \t]*\$(?:alpha|beta|closed)\s+shells\b.*?(?=^[ \t]*\$|\Z)')
    shells = sorted(m.group(0).strip() for m in pattern.finditer(control))
    if not shells:
        raise RuntimeError(f'{directory}: missing inherited shell occupations')
    if re.search(r'(?im)^[ \t]*\$fermi\b', control):
        raise RuntimeError(f'{directory}: $fermi unexpectedly present; inspect before screwer continuation')
    return {'orbitals_sha256': orbitals, 'shell_groups': shells}


def check_screwer_define():
    """Validate DEFINE without resetting the inherited SCF electronic state."""
    manifest = HERE / 'screwer_prepared.json'
    if not manifest.is_file():
        raise RuntimeError('missing screwer_prepared.json')
    info = json.loads(manifest.read_text())
    if info.get('workflow_version') != 2:
        raise RuntimeError('outdated screwer manifest; refusing incompatible restart')
    if (HERE / 'control.template').exists():
        raise RuntimeError('unexpected control.template in screwer continuation')
    for name in ('control', 'basis', 'auxbasis'):
        if not (HERE / name).is_file():
            raise RuntimeError(f'screwer DEFINE missing {name}')
    expected = info.get('inherited_state')
    current = screwer_inheritance(HERE)
    if not expected or current != expected:
        raise RuntimeError('DEFINE changed inherited orbital files or shell occupations; stop and inspect def_screwer.inp')
    (HERE / 'define_setup.ok').write_text(now() + '\n')
    record('DEFINE', 'OK', 'screwer DEFINE verified; inherited orbitals and occupations unchanged; no Fermi smearing')


def check_define():
    """Verify DEFINE ran against the sole control file; never rewrite its output."""
    meta = metadata()
    if (HERE / 'control.template').exists():
        raise RuntimeError('unexpected control.template; each metal directory must contain only control')
    for name in ('control', 'basis', 'auxbasis'):
        if not (HERE / name).is_file():
            raise RuntimeError(f'define did not produce {name}; inspect SLURM log')
    text = (HERE / 'control').read_text()
    matches = re.findall(r'(?im)^\s*\$fermi\b[^\n]*?\bnue\s*=\s*([+\-]?\d+(?:\.\d+)?)', text)
    if len(matches) != 1:
        raise RuntimeError(f'expected exactly one $fermi nue= in control after define; found {len(matches)}')
    if abs(float(matches[0]) - meta['unpaired']) > 1e-5:
        raise RuntimeError(f'control has nue={matches[0]}, but {meta["metal"]} needs nue={meta["unpaired"]}')
    save_initial_occupations()
    (HERE / 'define_setup.ok').write_text(now() + '\n')
    record('DEFINE', 'OK', f'control generated by define; charge={meta.get("charge") if meta.get("charge") is not None else "from def.inp"}, nue={meta["unpaired"]}; no post-define overlay')


def remove_fermi():
    p = HERE / 'control'
    original = p.read_text()
    # Restore the occupations saved directly after define before integer SCF.
    # Other groups ($energy, $grad, $ricore, etc.) are left untouched.
    restored, changed = restore_shells(original)
    fermi = [(start, end) for name, start, end in group_spans(restored) if name == 'fermi']
    if len(fermi) != 1:
        raise RuntimeError('expected exactly one $fermi group before jobex')
    start, end = fermi[0]
    edited = restored[:start] + restored[end:]
    (HERE / 'control.with_fermi').write_text(original)
    p.write_text(edited)
    record('INTEGER_SCF', 'RUNNING',
           '$fermi removed; reference alpha/beta occupations restored' if changed else
           '$fermi removed; reference alpha/beta occupations already matched')


def check_ridft(filename):
    p = HERE / filename
    text = p.read_text(errors='replace') if p.exists() else ''
    criteria = (bool(re.search(r'ENERGY\s+CONVERGED\s*!', text, re.I)) and
                bool(re.search(r'convergence\s+criteria\s+satisfied\s+after\s+\d+\s+iterations', text, re.I)) and
                bool(re.search(r'ridft\s*:\s*all\s+done', text, re.I)))
    if not criteria:
        record('RIDFT', 'FAILED', f'missing convergence signatures in {filename}')
        return 1
    match = re.search(r'convergence\s+criteria\s+satisfied\s+after\s+(\d+)\s+iterations', text, re.I)
    record('RIDFT', 'OK', f'{filename}; {match.group(1)} iterations')
    return 0


# In TURBOMOLE's eiger output a missing Occupation column value means an empty MO.
# The printed spectrum is a limited window, not a guarantee about all 1204 orbitals.
ORBITAL_RE = re.compile(
    r'^\s*(\d+)\.\s+([ab])\s+(\d+)\s+([A-Za-z0-9_+\-]+)\s+'
    r'(?:(\d+(?:\.\d+)?)\s+)?([+\-]?\d+(?:\.\d+)?)\s+H\s*=', re.M)


def orbitals(filename, stage):
    p = HERE / filename
    data = p.read_text(errors='replace') if p.exists() else ''
    rows = []
    for m in ORBITAL_RE.finditer(data):
        rows.append({'nr': int(m[1]), 'spin': m[2], 'index': int(m[3]),
                     'occupation': float(m[5]) if m[5] is not None else 0.0,
                     'energy_hartree': float(m[6])})
    eps = 1e-3
    occ = sorted((r for r in rows if r['occupation'] >= 1 - eps), key=lambda x: x['energy_hartree'])
    unocc = sorted((r for r in rows if r['occupation'] <= eps), key=lambda x: x['energy_hartree'])
    partial = [r for r in rows if eps < r['occupation'] < 1 - eps]
    top = max((r['energy_hartree'] for r in occ), default=None)
    # Combined spin channels can contain expected inversions for an open-shell state.
    global_holes = [r for r in unocc if top is not None and r['energy_hartree'] < top - 1e-6]
    spin_holes = []
    for spin in ('a', 'b'):
        highest = max((r['energy_hartree'] for r in occ if r['spin'] == spin), default=None)
        if highest is not None:
            spin_holes.extend(r for r in unocc if r['spin'] == spin and r['energy_hartree'] < highest - 1e-6)
    above = [r for r in unocc if top is not None and r['energy_hartree'] > top]
    above.sort(key=lambda r: r['energy_hartree'])
    result = {
        'metal': metadata()['metal'], 'stage': stage, 'source': filename,
        'parsed_orbitals': len(rows), 'frontier_homo_hartree': top,
        'first_five_above_homo': above[:5],
        'holes_global_spin_mixed': global_holes,
        'holes_within_each_spin': spin_holes,
        'fractional_occupations': partial,
        'coverage': 'Only orbitals printed by eiger (frontier window). Not a full-spectrum proof.',
        'status': 'UNKNOWN' if not rows or not re.search(r'eiger\s*:\s*all\s+done', data, re.I)
                  else ('WARNING' if global_holes or spin_holes or partial else 'NO_HOLES_IN_PRINTED_WINDOW'),
    }
    out = HERE / ('occupation_' + stage + '.json')
    out.write_text(json.dumps(result, indent=2) + '\n')
    line = (f'{now()} {stage} {result["status"]}: global={len(global_holes)}, '
            f'spin={len(spin_holes)}, fractional={len(partial)}, '
            f'parsed={len(rows)}; see {out.name}')
    with (HERE / 'occupation_warnings.log').open('a') as f:
        f.write(line + '\n')
    print(line, flush=True)
    record('EIGER_' + stage.upper(), result['status'], 'non-blocking occupation check')
    return 0


def check_numforce():
    # NumForce should complete normally AND create frequency data.
    p = HERE / 'numforce.out'
    output = p.read_text(errors='replace') if p.exists() else ''
    ok = (bool(output.strip()) and (HERE / 'hessian').is_file() and
          (HERE / 'hessian').stat().st_size > 0 and
          (HERE / 'vibspectrum').is_file() and
          (HERE / 'vibspectrum').stat().st_size > 0)
    if not ok:
        record('NUMFORCE', 'FAILED', 'missing/empty numforce.out, hessian or vibspectrum')
        return 1
    record('NUMFORCE', 'OK', 'numerical frequencies produced; inspect imaginary modes manually')
    return 0


def check_aoforce():
    text = (HERE / 'aoforce.out').read_text(errors='replace') if (HERE / 'aoforce.out').exists() else ''
    ok = bool(re.search(r'aoforce\s*:\s*all\s+done', text, re.I)) or (
        (HERE / 'hessian').exists() and (HERE / 'vibspectrum').exists() and bool(text.strip()))
    if not ok:
        record('AOFORCE', 'FAILED', 'no completion marker or hessian+vibspectrum outputs')
        return 1
    record('AOFORCE', 'OK', 'frequency output generated; inspect imaginary modes manually')
    return 0


def main():
    ap = argparse.ArgumentParser(
        description='Internal helper: records stage results and shell occupations.',
        epilog='Called automatically by the SLURM stages; normally do not run directly.',
    )
    sub = ap.add_subparsers(dest='command', required=True, help='stage helper action')
    sub.add_parser('check-define', help='verify files and $fermi nue after define without rewriting control')
    sub.add_parser('check-screwer-define', help='verify continuation DEFINE kept inherited orbitals, occupations, and no $fermi')
    sub.add_parser('remove-fermi', help='save control.with_fermi, restore define occupations and remove $fermi before integer SCF')
    sub.add_parser('restore-occupations', help='restore define occupations after integer SCF if changed (before jobex)')
    r = sub.add_parser('check-ridft', help='verify SCF convergence signatures in an output file'); r.add_argument('filename', help='ridft output file to check')
    e = sub.add_parser('eiger', help='log frontier orbital occupations and holes (warnings only)'); e.add_argument('filename', help='eiger output file'); e.add_argument('stage', choices=('pre', 'integer', 'post'), help='occupation check stage')
    s = sub.add_parser('record', help='update the per-metal workflow status and history'); s.add_argument('stage', help='calculation stage'); s.add_argument('status', help='stage state'); s.add_argument('message', nargs='?', default='', help='optional status details')
    sub.add_parser('check-aoforce', help='verify aoforce finished and generated expected results')
    sub.add_parser('check-numforce', help='verify NumForce generated numerical Hessian and frequencies')
    args = ap.parse_args()
    try:
        if args.command == 'check-define': check_define()
        elif args.command == 'check-screwer-define': check_screwer_define()
        elif args.command == 'remove-fermi': remove_fermi()
        elif args.command == 'restore-occupations': restore_after_integer_scf()
        elif args.command == 'check-ridft': return check_ridft(args.filename)
        elif args.command == 'eiger': return orbitals(args.filename, args.stage)
        elif args.command == 'record': record(args.stage, args.status, args.message)
        elif args.command == 'check-aoforce': return check_aoforce()
        elif args.command == 'check-numforce': return check_numforce()
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
