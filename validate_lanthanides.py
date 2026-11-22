#!/usr/bin/env python3
"""Read-only TURBOMOLE input preflight (turbomoleio-backed).

Does not invoke DEFINE, TURBOMOLE, SLURM, or write any files. It inspects
project masters, prepared calculations, and Screwer continuation directories.
A passing report means only that the specified *input* invariants passed.
"""
import argparse
import json
import math
from pathlib import Path
import re
import sys

from setup_lanthanides import SPIN, make_coord, make_control, make_define, read_charge

GROUP_RE = re.compile(r'(?m)^\s*\$([A-Za-z][\w-]*)\b')
# TURBOMOLE may append repeated result groups during optimizations/Hessians.
# Duplicate output records are not invalid initialization inputs.
REPEATABLE_OUTPUT = {'last', 'vibrational', 'nprvibrational'}
FILE_RE = re.compile(r'\bfile\s*=\s*([^\s]+)', re.I)
NUMERIC_RE = re.compile(r'^[+-]?\d+$')
SCIENTIFIC_FLAGS = ('disp4', 'noproj', 'rij', 'scfdiis', 'uhf')
REFERENCE_GROUPS = ('basis', 'ecp', 'jbas', 'jkbas', 'cbas', 'coord',
                    'uhfmo_alpha', 'uhfmo_beta', 'scfmo')
RELATIVISTIC_RE = re.compile(r'^(?:x2c(?:[-_].*)?|dkh(?:[-_].*)?|bss(?:[-_].*)?|soghf|soc)$', re.I)


class Report:
    def __init__(self, label):
        self.label = label
        self.issues = []

    def add(self, level, message):
        self.issues.append({'level': level, 'message': message})

    def error(self, message):
        self.add('FAIL', message)

    def warn(self, message):
        self.add('WARN', message)

    @property
    def status(self):
        if any(x['level'] == 'FAIL' for x in self.issues):
            return 'FAIL'
        if any(x['level'] == 'WARN' for x in self.issues):
            return 'WARN'
        return 'PASS'

    def serialize(self):
        return {'target': self.label, 'status': self.status, 'issues': self.issues}


def read_control(directory, report, Control):
    path = directory / 'control'
    if not path.is_file():
        report.error('missing control')
        return None
    try:
        source = path.read_text()
        if not re.search(r'(?m)^\s*\$end\s*$', source):
            report.error('control has no terminating $end')
        names = [m.group(1).lower() for m in GROUP_RE.finditer(source)]
        for name in set(names):
            if name not in REPEATABLE_OUTPUT and name != 'end' and names.count(name) > 1:
                report.error(f'duplicate ${name} data group')
        return Control.from_file(str(path))
    except Exception as exc:
        report.error(f'turbomoleio could not parse control: {exc}')
        return None


def data(control, name, report):
    if control is None:
        return None
    try:
        return control.show_data_group(name, show_from_subfile=False)
    except Exception as exc:
        report.error(f'cannot inspect ${name}: {exc}')
        return None


def atom_symbols(path, report):
    if not path.is_file():
        report.error('missing coord')
        return None
    try:
        lines = path.read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        report.error(f'cannot read coord: {exc}')
        return None
    active = False
    atoms = []
    end_found = False
    for line in lines:
        s = line.strip()
        if not active:
            if s.lower().startswith('$coord'):
                active = True
            continue
        if s.startswith('$'):
            end_found = True
            break
        if not s or s.startswith('#'):
            continue
        parts = s.split()
        if len(parts) < 4:
            report.error('coord has a malformed atom row')
            return None
        try:
            xyz = [float(v.replace('D', 'E').replace('d', 'e')) for v in parts[:3]]
            if not all(math.isfinite(v) for v in xyz):
                raise ValueError('nonfinite')
        except ValueError:
            report.error('coord has a non-finite or invalid Cartesian position')
            return None
        if not re.fullmatch(r'[A-Za-z]{1,3}', parts[3]):
            report.error('coord has an invalid atom symbol')
            return None
        atoms.append(parts[3].lower())
    if not active or not end_found or not atoms:
        report.error('coord needs a nonempty $coord group and a closing group marker')
        return None
    return atoms


def check_refs(directory, ctl, report, *, require_generated):
    if ctl is None:
        return
    for group in REFERENCE_GROUPS:
        value = data(ctl, group, report)
        if value is None:
            continue
        matches = FILE_RE.findall(value)
        if len(matches) > 1:
            report.error(f'${group} declares multiple file references')
            continue
        if not matches:
            continue
        filename = matches[0].strip('"\'')
        target = directory / filename
        # Project masters intentionally refer to basis/auxbasis and orbitals
        # before DEFINE has generated them. Existing calculations may not.
        if require_generated and (not target.is_file() or target.stat().st_size == 0):
            report.error(f'${group}: missing/empty referenced file {filename}')
        if not require_generated and not target.exists() and group not in (
                'basis', 'ecp', 'jbas', 'jkbas', 'cbas',
                'uhfmo_alpha', 'uhfmo_beta', 'scfmo'):
            report.error(f'${group}: missing referenced file {filename}')


def method_details(ctl, report):
    if ctl is None:
        return None
    dft = data(ctl, 'dft', report)
    if dft is None:
        report.error('missing $dft')
        return None
    parsed = {}
    for key in ('functional', 'gridsize'):
        found = re.search(r'(?im)^\s*' + key + r'\s+(\S+)', dft)
        parsed[key] = found.group(1).lower() if found else None
        if parsed[key] is None:
            report.error(f'$dft has no explicit {key} (cannot verify intended settings)')
    parsed['weight_derivatives'] = bool(re.search(r'(?im)^\s*weight\s+derivatives\s*$', dft))
    parsed['scfconv'] = (data(ctl, 'scfconv', report) or '').strip()
    if not parsed['scfconv']:
        report.error('missing explicit $scfconv')
    return parsed


def compare_science(original, candidate, report, *, source):
    """Compare stable method settings; do not compare orbitals or SCF results."""
    if original is None or candidate is None:
        return
    master = method_details(original, report)
    child = method_details(candidate, report)
    if master and child:
        for key in ('functional', 'gridsize', 'weight_derivatives', 'scfconv'):
            if master[key] != child[key]:
                report.error(f'{key} differs from {source}: {child[key]!r} != {master[key]!r}')
    for name in SCIENTIFIC_FLAGS:
        a = data(original, name, report) is not None
        b = data(candidate, name, report) is not None
        if a != b:
            report.error(f'${name} presence differs from {source}')
    a = data(original, 'ecp', report) is not None
    b = data(candidate, 'ecp', report) is not None
    if a != b:
        report.error(f'$ecp presence differs from {source} (ECP/all-electron method mismatch)')
    # Include X2C / DKH / related data groups for method identity.
    a_names = {m.group(1).lower() for m in GROUP_RE.finditer(str(original)) if RELATIVISTIC_RE.match(m.group(1))}
    b_names = {m.group(1).lower() for m in GROUP_RE.finditer(str(candidate)) if RELATIVISTIC_RE.match(m.group(1))}
    if a_names != b_names:
        report.error(f'relativistic data-group names differ from {source}: {sorted(b_names)} vs {sorted(a_names)}')


def check_geometry(directory, ctl, report, metal=None, *, prepared=False):
    atoms = atom_symbols(directory / 'coord', report)
    if atoms is None:
        return
    if metal:
        present = [x for x in atoms if x in SPIN or x == 'pm']
        if present != [metal]:
            report.error(f'expected exactly one {metal} atom; found lanthanides {present}')
    if prepared and data(ctl, 'coord', report) is None:
        report.error('control is missing $coord after DEFINE')
    if prepared and data(ctl, 'atoms', report) is None:
        report.error('control is missing $atoms after DEFINE')
    rundimensions = data(ctl, 'rundimensions', report)
    if rundimensions:
        found = re.search(r'\bnatoms\s*=\s*(\d+)', rundimensions)
        if found and int(found.group(1)) != len(atoms):
            report.error(f'$rundimensions natoms={found.group(1)} but coord has {len(atoms)} atoms')


def check_project(base, Control, metals):
    report = Report('project master')
    for filename in ('control', 'coord', 'def.inp'):
        if not (base / filename).is_file():
            report.error(f'missing project master {filename}')
    ctl = read_control(base, report, Control)
    check_geometry(base, ctl, report)
    check_refs(base, ctl, report, require_generated=False)
    if ctl is not None:
        method_details(ctl, report)
    try:
        master_control = (base / 'control').read_text()
        coord = (base / 'coord').read_text()
    except OSError:
        return report, ctl
    for metal in metals:
        template = base / ('def_gd.inp' if metal == 'gd' else 'def.inp')
        if not template.is_file():
            report.error(f'{metal}: missing project template {template.name}')
            continue
        try:
            make_define(template.read_text(), metal, SPIN[metal])
            make_control(master_control, SPIN[metal])
            make_coord(coord, metal)
        except Exception as exc:
            report.error(f'{metal}: cannot generate prepared inputs: {exc}')
        if read_charge(template.read_text()) is None:
            report.warn(f'{metal}: charge not readable from {template.name}; do not infer it from $charge from ridft')
    return report, ctl


def defined_inputs_present(folder, ctl, report, *, screwer=False):
    """Distinguish not-yet-defined directories from incomplete DEFINE output.

    A project master and an auto-prepared metal do not have $atoms/$coord in
    their control yet. Screwer preparation may intentionally remove $coord.
    """
    has_coord = data(ctl, 'coord', report) is not None
    has_atoms = data(ctl, 'atoms', report) is not None
    if has_coord and has_atoms:
        return True
    if screwer:
        # The preparation stage has its own manifest and has not started DEFINE.
        # Old inherited computation outputs may already be present via cpc.
        fresh = (folder / 'screwer_prepared.json').is_file() and not any(
            (folder / name).exists() for name in ('define.out', 'define.fallback.out'))
    else:
        # Genuine fresh lan-setup inputs have no existing computational results.
        fresh = (folder / 'metal.json').is_file() and (folder / 'def.inp').is_file() and not any(
            (folder / name).exists() for name in ('define.out', 'ridft.out', 'jobex.out',
                                                   'define_setup.ok', 'GEO_OPT_CONVERGED',
                                                   'workflow_status.json', 'basis', 'alpha', 'beta', 'mos'))
    if fresh:
        report.add('INFO', 'prepared before DEFINE; generated $coord/$atoms and basis/orbital files are not yet required')
    else:
        missing = [x for x, present in (('coord', has_coord), ('atoms', has_atoms)) if not present]
        report.error('incomplete existing DEFINE/control: missing ' + ', '.join('$' + x for x in missing))
    return False


def check_metal(base, metal, master, Control):
    folder = base / metal
    report = Report(metal)
    if not folder.is_dir():
        report.warn('metal folder not present; no existing inputs to inspect')
        return report, None
    ctl = read_control(folder, report, Control)
    ready = defined_inputs_present(folder, ctl, report) if ctl is not None else False
    check_geometry(folder, ctl, report, metal, prepared=ready)
    check_refs(folder, ctl, report, require_generated=ready)
    compare_science(master, ctl, report, source='project master')
    meta = folder / 'metal.json'
    if meta.is_file():
        try:
            record = json.loads(meta.read_text())
            if record.get('metal') != metal:
                report.error(f'metal.json metal={record.get("metal")!r} instead of {metal!r}')
            if record.get('unpaired') is not None and record['unpaired'] != SPIN[metal]:
                report.warn(f'metal.json unpaired={record["unpaired"]} differs from initial Ln3+ high-spin expectation {SPIN[metal]}')
            template = folder / 'def.inp'
            if template.is_file():
                expected = read_charge(template.read_text())
                if expected is not None and record.get('charge') != expected:
                    report.error(f'metal.json charge={record.get("charge")} differs from local def.inp EHT charge={expected}')
        except (ValueError, TypeError, AttributeError, OSError) as exc:
            report.error(f'cannot read metal.json: {exc}')
    else:
        report.warn('manual/legacy metal folder: metal.json absent; charge/spin provenance not verified')
    if not (folder / 'def.inp').is_file() and meta.is_file():
        report.warn('no local def.inp to verify initialization charge')
    return report, ctl


def check_screwer(base, metal, parent, parent_ctl, Control):
    child = parent / f'screwer_{metal}'
    report = Report(f'{metal}/screwer_{metal}')
    if not child.is_dir():
        report.warn('no Screwer correction folder')
        return report
    ctl = read_control(child, report, Control)
    ready = defined_inputs_present(child, ctl, report, screwer=True) if ctl is not None else False
    check_geometry(child, ctl, report, metal, prepared=ready)
    check_refs(child, ctl, report, require_generated=True)
    if parent_ctl is None:
        report.warn('parent control unavailable: inheritance settings could not be checked')
    else:
        compare_science(parent_ctl, ctl, report, source=f'{metal}/control')
    if not (child / 'screwer_prepared.json').is_file():
        report.warn('missing screwer_prepared.json provenance')
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description='Read-only turbomoleio-backed scientific input preflight; never submits or edits jobs.')
    ap.add_argument('--base', type=Path, default=Path.cwd(), help='EF project directory (default current directory)')
    ap.add_argument('--metals', nargs='+', choices=tuple(SPIN), help='specific metals; default existing metal folders')
    ap.add_argument('--project-only', action='store_true', help='inspect only project master inputs')
    ap.add_argument('--include-screwer', action='store_true', help='also inspect screwer_<metal> folders')
    ap.add_argument('--json', action='store_true', help='print a JSON report to stdout; no files written')
    args = ap.parse_args(argv)
    try:
        from turbomoleio.core.control import Control
    except ImportError as exc:
        print('ERROR: turbomoleio is not installed in this Python environment. '
              'Activate an isolated venv and run: python -m pip install "turbomoleio==1.6.0" '
              f'(details: {exc})', file=sys.stderr)
        return 2
    base = args.base.resolve()
    # From a metal/ folder, use the containing EF project automatically.
    if (base.name in SPIN or base.name.endswith('_closedshell')) and (base.parent / 'def.inp').is_file():
        base = base.parent
    metals = args.metals if args.metals is not None else [m for m in SPIN if (base / m).is_dir()]
    if args.project_only:
        metals = args.metals if args.metals is not None else list(SPIN)
    root_report, master = check_project(base, Control, metals)
    results = [root_report]
    if not args.project_only:
        for metal in metals:
            report, ctl = check_metal(base, metal, master, Control)
            results.append(report)
            if args.include_screwer and (base / metal / f'screwer_{metal}').is_dir():
                results.append(check_screwer(base, metal, base / metal, ctl, Control))
    if args.json:
        print(json.dumps({'base': str(base), 'read_only': True,
                          'results': [x.serialize() for x in results]}, indent=2))
    else:
        print(f'Read-only input validation: {base}')
        for result in results:
            print(f'  {result.status:4s}  {result.label}')
            for issue in result.issues:
                print(f'      {issue["level"]}: {issue["message"]}')
        print('NOTE: This does not certify SCF/geometry convergence, the electronic ground state, or job status.')
    return 1 if any(r.status == 'FAIL' for r in results) else 0


if __name__ == '__main__':
    sys.exit(main())
