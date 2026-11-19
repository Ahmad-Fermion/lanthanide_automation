#!/usr/bin/env python3
"""TURBOMOLE helper: validate DEFINE output, check SCF, inspect eiger, track stage."""
import argparse
import hashlib
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path.cwd()
META = HERE / 'metal.json'
STATE = HERE / 'workflow_status.json'
OCC_REFERENCE = HERE / 'initial_occupations'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def metadata():
    if META.is_file():
        return json.loads(META.read_text())
    if os.environ.get('TM_MANUAL_RECOVERY') == '1':
        # Metadata-free legacy Hessian recovery only; no charge or spin guess.
        metal = os.environ.get('TM_MANUAL_METAL', '')
        if metal not in {'la', 'ce', 'pr', 'nd', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu'}:
            raise RuntimeError('manual recovery requires a verified metal identity')
        return {'metal': metal, 'unpaired': None, 'charge': None,
                'source': 'manual Hessian-only recovery'}
    raise RuntimeError('metal.json missing (not an authorized manual Hessian recovery)')


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


def _orbital_data_sha256(path):
    """Hash orbital data, not DEFINE's formatting and comments.

    TURBOMOLE can remove SCF-energy comments or rewrite padded MO headers.
    Preserve the numerical characters and orbital ordering; only normalize
    whitespace and discard full-line comments.
    """
    digest = hashlib.sha256()
    rows = 0
    with path.open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            data = line.strip()
            if not data or data.startswith('#'):
                continue
            if data.startswith('$') or 'eigenvalue=' in data.lower():
                data = ' '.join(data.split())
            else:
                # Packed D-format coefficients may be separated by optional
                # spaces; keep all numerical characters and their order.
                data = ''.join(data.split())
            digest.update(data.encode('utf-8') + b'\n')
            rows += 1
    if not rows:
        raise RuntimeError(f'{path}: empty orbital data')
    return digest.hexdigest()


def _raw_sha256(path):
    checksum = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def screwer_inheritance(directory):
    """Fingerprint parent MOs and shell groups; Screwer must reuse this state.

    This is recorded *after* generating the displaced coord, before DEFINE.
    It does not guarantee electronic-state correctness, only continuity.
    """
    orb_files = ['alpha', 'beta'] if all((directory / f).is_file() for f in ('alpha', 'beta')) else (
        ['mos'] if (directory / 'mos').is_file() else [])
    if not orb_files:
        raise RuntimeError(f'{directory}: inherited alpha/beta or mos orbital files missing; refusing screwer continuation')
    orbitals = {name: _orbital_data_sha256(directory / name) for name in orb_files}
    control = (directory / 'control').read_text()
    # Find whole shell blocks with explicit data-group boundaries.
    pattern = re.compile(r'(?ims)^[ \t]*\$(?:alpha|beta|closed)\s+shells\b.*?(?=^[ \t]*\$|\Z)')
    shells = sorted(m.group(0).strip() for m in pattern.finditer(control))
    if not shells:
        raise RuntimeError(f'{directory}: missing inherited shell occupations')
    if re.search(r'(?im)^[ \t]*\$fermi\b', control):
        print(f'WARNING {directory}: inherited $fermi retained as configured; '
              'verify intended SCF smearing and occupations', flush=True)
    return {'orbitals_sha256': orbitals, 'shell_groups': shells}


def check_screwer_define():
    """Legacy command for existing child scripts; do NOT compare MO/control data.

    Screwer DEFINE can legitimately transform inherited orbitals and update
    control, especially when the input uses `use ../control`. The actual
    DEFINE exit status and subsequent RIDFT/JOBEX checks determine success.
    """
    record('DEFINE', 'OK', 'Screwer DEFINE finished; no orbital/control comparison')



# Initial Ln DEFINE may fail when ired tries to construct redundant internal
# coordinates for a changed metal. Unlike Screwer, these are fresh inputs;
# retain the metal-specific orbital, EHT, basis and relativistic instructions.
INITIAL_DEFINE_OUTPUTS = ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta',
                          'mos', 'statistics', 'intdef', 'wherefrom', 'gradient',
                          'energy', 'hessapprox', 'job.start', 'job.last')
IREDS_COORD_MENU = re.compile(
    r'(?im)^([ \t]*a[ \t]+coord[ \t]*\r?\n)[ \t]*ired[ \t]*\r?\n[ \t]*\*[ \t]*(?:\r?\n|$)')


def initial_noired_transcript(script):
    """Replace precisely the geometry menu path; never rewrite chemistry menus.

    DEFINE's no-ired route needs BOTH the geometry terminator '*' and the
    subsequent 'no'. Simply deleting 'ired' shifts all subsequent answers.
    """
    revised, count = IREDS_COORD_MENU.subn(lambda m: m.group(1) + '*\nno\n', script)
    if count > 1:
        raise RuntimeError('multiple a coord / ired / * sections; cannot choose a safe retry')
    return revised if count == 1 else None


def run_initial_define():
    """Run an initial DEFINE transcript; retry cleanly without ired if it fails.

    The original def.inp stays authoritative and is never changed. A retry
    starts from a snapshot of the PRE-DEFINE inputs, never the failed control.
    Successful DEFINE is still checked by `check-define` before SCF.
    """
    source = HERE / 'def.inp'
    if not source.is_file() or not source.stat().st_size:
        raise RuntimeError(f'{source}: missing or empty DEFINE transcript')
    script = source.read_text()
    noired = initial_noired_transcript(script)

    def attempt(input_text, filename):
        result = subprocess.run(['define'], cwd=HERE, input=input_text,
                                text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        (HERE / filename).write_text(result.stdout)
        success = (result.returncode == 0 and
                   re.search(r'(?i)define\s*:\s*all\s+done|define\s+ended\s+normally',
                             result.stdout) is not None and
                   all((HERE / n).is_file() and (HERE / n).stat().st_size
                       for n in ('control', 'coord', 'basis', 'auxbasis')))
        return success, result.returncode

    # Templates already using 'a coord / * / no' run once as before.
    if noired is None:
        ok, code = attempt(script, 'define.out')
        if not ok:
            raise RuntimeError(f'initial DEFINE failed (exit {code}); no recognized '
                               'a coord / ired / * menu path for retry; see define.out')
        return

    with tempfile.TemporaryDirectory(prefix='.initial_define_snapshot_', dir=HERE) as temp:
        backup = Path(temp)
        present = set()
        for name in INITIAL_DEFINE_OUTPUTS:
            item = HERE / name
            if item.is_symlink():
                raise RuntimeError(f'{item}: cannot safely snapshot a symlink')
            if item.is_file():
                shutil.copy2(item, backup / name)
                present.add(name)

        def restore():
            for name in INITIAL_DEFINE_OUTPUTS:
                target = HERE / name
                if target.is_symlink():
                    raise RuntimeError(f'{target}: cannot safely restore a symlink')
                if name in present:
                    shutil.copy2(backup / name, target)
                elif target.is_file():
                    target.unlink()

        ok, code = attempt(script, 'define.out')
        if ok:
            print('Initial DEFINE succeeded with original input', flush=True)
            return
        print(f'WARNING initial DEFINE failed (exit {code}); '
              'retrying without ired from original inputs', flush=True)
        restore()
        (HERE / 'def.noired.generated.inp').write_text(noired)
        ok, code = attempt(noired, 'define.noired.out')
        if not ok:
            restore()
            raise RuntimeError(f'initial DEFINE and no-ired fallback both failed '
                               f'(fallback exit {code}); prepared inputs restored; '
                               'see define.out and define.noired.out')
        print('Initial DEFINE fallback succeeded (a coord / * / no)', flush=True)

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


def check_restricted_define():
    """Verify DEFINE really generated a restricted closed-shell input.

    Does not infer correct *physical* spin multiplicity, only that the user-
    requested RHF/RKS-style orbital representation was actually prepared.
    """
    meta = metadata()
    if not meta.get('restricted') or meta.get('unpaired') != 0:
        raise RuntimeError('restricted validation requires opt-in metadata with zero unpaired electrons')
    define_log = (HERE / 'define.out').read_text(errors='replace') if (HERE / 'define.out').is_file() else ''
    if not re.search(r'(?i)define\s+ended\s+normally', define_log):
        raise RuntimeError('restricted DEFINE did not end normally; see define.out')
    if re.search(r'(?i)ERROR IN CARTESIAN COORDINATES\s*\?', define_log):
        raise RuntimeError('restricted DEFINE reported a geometry error; see define.out')
    ctl = (HERE / 'control').read_text()
    for group in ('uhf', 'uhfmo_alpha', 'uhfmo_beta', 'alpha shells', 'beta shells', 'fermi'):
        if re.search(r'(?im)^\s*\$' + re.escape(group) + r'\b', ctl):
            raise RuntimeError(f'RESTRICTED DEFINE unexpectedly produced ${group}')
    if not re.search(r'(?im)^\s*\$closed\s+shells\b', ctl):
        raise RuntimeError('DEFINE did not produce $closed shells')
    if not re.search(r'(?im)^\s*\$scfmo\b', ctl):
        raise RuntimeError('DEFINE did not produce $scfmo; expected restricted orbital reference')
    for group in ('coord', 'atoms', 'dft', 'scfconv', 'basis'):
        if not re.search(r'(?im)^\s*\$' + group + r'\b', ctl):
            raise RuntimeError(f'restricted control lacks ${group}')
    if not (HERE / 'mos').is_file() or not (HERE / 'mos').stat().st_size:
        raise RuntimeError('DEFINE did not produce nonempty restricted mos file')
    for name in ('basis', 'auxbasis'):
        if not (HERE / name).is_file() or not (HERE / name).stat().st_size:
            raise RuntimeError(f'restricted DEFINE did not generate {name}')
    record('DEFINE', 'OK', 'restricted orbitals and closed-shell occupations verified')


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
    sub.add_parser('check-restricted-define', help='verify restricted molecular orbitals and $closed shells')
    sub.add_parser('check-define', help='verify files and $fermi nue after define without rewriting control')
    sub.add_parser('run-initial-define', help='run initial DEFINE with a clean no-ired retry when needed')
    sub.add_parser('check-screwer-define', help='legacy nonblocking Screwer DEFINE completion record')
    sub.add_parser('remove-fermi', help='save control.with_fermi, restore define occupations and remove $fermi before integer SCF')
    sub.add_parser('restore-occupations', help='restore define occupations after integer SCF if changed (before jobex)')
    r = sub.add_parser('check-ridft', help='verify SCF convergence signatures in an output file'); r.add_argument('filename', help='ridft output file to check')
    e = sub.add_parser('eiger', help='log frontier orbital occupations and holes (warnings only)'); e.add_argument('filename', help='eiger output file'); e.add_argument('stage', choices=('pre', 'integer', 'post'), help='occupation check stage')
    s = sub.add_parser('record', help='update the per-metal workflow status and history'); s.add_argument('stage', help='calculation stage'); s.add_argument('status', help='stage state'); s.add_argument('message', nargs='?', default='', help='optional status details')
    sub.add_parser('check-aoforce', help='verify aoforce finished and generated expected results')
    sub.add_parser('check-numforce', help='verify NumForce generated numerical Hessian and frequencies')
    args = ap.parse_args()
    try:
        if args.command == 'check-restricted-define': check_restricted_define()
        elif args.command == 'check-define': check_define()
        elif args.command == 'run-initial-define': run_initial_define()
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
