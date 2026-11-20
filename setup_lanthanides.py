#!/usr/bin/env python3
"""Prepare Ln inputs using one control file; preserve all existing folders."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

# 4f electron high-spin Ln3+ (exclude Pm entirely)
SPIN = {'la': 0, 'ce': 1, 'pr': 2, 'nd': 3, 'sm': 5, 'eu': 6, 'gd': 7,
        'tb': 6, 'dy': 5, 'ho': 4, 'er': 3, 'tm': 2, 'yb': 1, 'lu': 0}
ROOT = Path(__file__).resolve().parent
RE_LN = re.compile(r'^(la|ce|pr|nd|pm|sm|eu|gd|tb|dy|ho|er|tm|yb|lu)$', re.I)


def make_coord(text, metal):
    """Replace the element field of exactly one Ln atom within $coord.

    TURBOMOLE permits extra per-atom fields, especially ``f`` for a fixed
    atom.  These fields, original whitespace, and other data groups must
    survive preparation unchanged.
    """
    changed = 0
    in_coord = False
    out = []
    for line in text.splitlines(keepends=True):
        stripped = line.lstrip()
        if re.match(r'(?i)^\$coord(?=\s|$)', stripped):
            in_coord = True
        elif stripped.startswith('$'):
            in_coord = False
        elif in_coord and not stripped.startswith('#'):
            parts = line.split()
            if len(parts) >= 4 and RE_LN.fullmatch(parts[3]):
                # Locate the fourth *token*, not the end of the row: a fixed
                # metal may appear as ``... gd f`` rather than ``... gd``.
                hit = re.match(r'^(\s*\S+\s+\S+\s+\S+\s+)(\S+)', line)
                if hit is None or hit.group(2) != parts[3]:
                    raise ValueError('cannot locate lanthanide element in coord row')
                line = line[:hit.start(2)] + metal + line[hit.end(2):]
                changed += 1
        out.append(line)
    if changed != 1:
        raise ValueError(f'coord must contain exactly one lanthanide atom (found {changed})')
    return ''.join(out)


def make_define(text, metal, spin):
    # Keep whichever basis the master define input uses (including X2C bases).
    text, n1 = re.subn(r'(?im)^(\s*b\s+)["\'](?:la|ce|pr|nd|pm|sm|eu|gd|tb|dy|ho|er|tm|yb|lu)["\'](\s+\S+\s*)$',
                      lambda m: m.group(1) + '"' + metal + '"' + m.group(2), text)
    text, n2 = re.subn(r'(?im)^(\s*u\s+)\d+(\s*)$',
                      lambda m: m.group(1) + str(spin) + m.group(2), text)
    if n1 != 1 or n2 != 1:
        raise ValueError(f'def.inp has unexpected format (metal_basis={n1}; unpaired={n2})')
    return text


def read_charge(define_text):
    """Read the project-defined EHT charge for reporting only; never constrain it."""
    match = re.search(r'(?im)^[ \t]*eht[ \t]*\r?\n(?:[ \t]*\r?\n)*[ \t]*([+-]?\d+)[ \t]*(?:\r?\n|$)', define_text)
    return int(match.group(1)) if match else None


def make_control(text, spin):
    """Keep the user's original control, changing only $fermi nue for this Ln3+."""
    pattern = r'(?im)^(\s*\$fermi\b[^\n]*?\bnue\s*=\s*)([+\-]?\d+(?:\.\d+)?)'
    updated, count = re.subn(pattern, lambda m: m.group(1) + f'{spin:.2f}', text)
    if count != 1:
        raise ValueError(f'control must have exactly one $fermi nue= setting (found {count})')
    return updated


def setup_one(base, metal):
    dest = base / metal
    if dest.exists():
        print(f'SKIP {metal}: directory already exists (left completely unchanged)')
        return False
    # Gd's EHT menu differs from the other metals. Its master input has
    # the extra response needed by define; retain its contents and just
    # normalize the Gd basis selector and spin, as for the other metals.
    define_source = base / ('def_gd.inp' if metal == 'gd' else 'def.inp')
    if not define_source.is_file():
        raise FileNotFoundError(f'missing master {define_source.name} in {base}; required when creating {metal}/')
    created_by_us = False
    try:
        # Request cpc first, as in the original workflow. Missing template-linked output
        # files can cause cpc to fail; in that case create the same clean directory.
        cpc = shutil.which('cpc')
        if cpc:
            cmd = subprocess.run([cpc, metal], cwd=base, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if dest.is_dir():
                created_by_us = True
            if cmd.returncode:
                print(f'NOTE {metal}: cpc returned {cmd.returncode}; using clean inputs ({cmd.stdout.strip()[:200]})')
        if not dest.exists():
            dest.mkdir()
            created_by_us = True
        if not dest.is_dir():
            raise RuntimeError(f'{dest} is not a directory')
        # Remove stale copies made by cpc. DEFINE must generate new basis/MOs for this metal.
        for p in dest.iterdir():
            if p.is_file() or p.is_symlink():
                p.unlink()
            elif p.is_dir():
                shutil.rmtree(p)
        (dest / 'coord').write_text(make_coord((base / 'coord').read_text(), metal))
        define_text = make_define(define_source.read_text(), metal, SPIN[metal])
        (dest / 'def.inp').write_text(define_text)
        charge = read_charge(define_text)
        (dest / 'control').write_text(make_control((base / 'control').read_text(), SPIN[metal]))
        shutil.copy2(ROOT / 'turbo.pipeline.slurm', dest / 'turbo.pipeline.slurm')
        shutil.copy2(ROOT / 'turbo.resume.slurm', dest / 'turbo.resume.slurm')
        shutil.copy2(ROOT / 'turbo.screwer.slurm', dest / 'turbo.screwer.slurm')
        # The in-folder screwer gate imports make_define/SPIN from this helper.
        shutil.copy2(ROOT / 'setup_lanthanides.py', dest / 'setup_lanthanides.py')
        shutil.copy2(ROOT / 'screwer_workflow.py', dest / 'screwer_workflow.py')
        shutil.copy2(ROOT / 'adaptive_define.py', dest / 'adaptive_define.py')
        shutil.copy2(ROOT / 'tm_workflow.py', dest / 'tm_workflow.py')
        (dest / 'metal.json').write_text(json.dumps({'metal': metal, 'unpaired': SPIN[metal],
                                                     'charge': charge}, indent=2) + '\n')
        print(f'CREATED {metal}: charge={charge if charge is not None else "from def.inp"}, unpaired={SPIN[metal]}')
        return True
    except Exception:
        # A pre-existing user folder never reaches this block because we skip it above.
        if created_by_us and dest.exists():
            shutil.rmtree(dest)
        raise



def prepare_default_restricted(base, metal):
    """Prepare complementary restricted La/Lu input without submitting it.

    Standard unrestricted projects must keep working without a restricted
    DEFINE transcript. Existing folders are never overwritten.
    """
    if metal not in ('la', 'lu'):
        return
    if not (base / 'def_closedshell.inp').is_file():
        print(f'NOTE {metal}_closedshell: not prepared (no def_closedshell.inp in project root)')
        return
    from closed_shell_workflow import prepare
    try:
        prepare(base, metal)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        # An optional restricted template must not block standard UHF runs.
        print(f'WARNING {metal}_closedshell: not prepared ({exc})', file=sys.stderr)



def main():
    if '--closed-shell' in sys.argv[1:]:
        raw = [x for x in sys.argv[1:] if x != '--closed-shell']
        if '--prepare-only' not in raw:
            raw.append('--prepare-only')
        from closed_shell_workflow import main as closed_main
        return closed_main(raw)
    ap = argparse.ArgumentParser(
        description='Prepare lanthanide folders (14 metals, excluding Pm); do not submit jobs.',
        epilog=(
            'Run from a project folder with master control, coord, def.inp and '
            'def_gd.inp (needed only for new Gd). Existing metal folders are untouched. '
            'New folders receive metal-specific inputs, control and SLURM helpers.\n\n'
            'Examples:\n'
            '  lan-setup --metals sm\n'
            '  lan-setup --metals gd tb\n'
            '  lan-setup\n'
            'To prepare and submit in one command, use: lan-run num --metals sm'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--base', type=Path, default=Path.cwd(),
                    help='project directory containing master control, coord and def.inp; also def_gd.inp if creating Gd (default: current directory)')
    ap.add_argument('--metals', nargs='+', choices=list(SPIN), default=list(SPIN),
                    metavar='METAL', help='one or more symbols to prepare (e.g. sm tb er); default: all 14, excluding Pm')
    args = ap.parse_args()
    base = args.base.resolve()
    for name in ('control', 'coord', 'def.inp'):
        if not (base / name).is_file():
            sys.exit(f'ERROR: missing master {name} in {base}')
    created = []
    for metal in args.metals:
        try:
            if setup_one(base, metal): created.append(metal)
            prepare_default_restricted(base, metal)
        except Exception as e:
            print(f'ERROR {metal}: {e}', file=sys.stderr)
            return 1
    print('To submit already prepared folders: lan-run num --include-existing')
    return 0


if __name__ == '__main__':
    sys.exit(main())
