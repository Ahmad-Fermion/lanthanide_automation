"""Offline contract tests of the read-only validator using turbomoleio's public API shape.

The real turbomoleio-1.6.0 parser must also be tested on VACC (integration test).
"""
import contextlib
import hashlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import validate_lanthanides as val


class FakeControl:
    """Minimal stand-in for the documented Control.from_file/show_data_group API."""
    def __init__(self, text):
        self.text = text

    @classmethod
    def from_file(cls, path):
        return cls(Path(path).read_text())

    def __str__(self):
        return self.text

    def show_data_group(self, name, show_from_subfile=False):
        groups = re.findall(r'(?ms)^\$([\w-]+)([^$]*)', self.text)
        found = [x[1] for x in groups if x[0].lower() == name.lower()]
        if len(found) > 1:
            raise RuntimeError('duplicate group')
        return found[0] if found else None


MASTER = '''$title
$noproj
$scfiterlimit 3000
$basis file=basis
$ecp file=basis
$uhfmo_alpha file=alpha
$uhfmo_beta file=beta
$uhf
$scfdiis
$dft
    functional tpssh
    gridsize 4
    weight derivatives
$scfconv 7
$disp4
$rij
$jbas file=auxbasis
$fermi tmstrt=1000 nue= 6.00
$end
'''

COORD = '''$coord
 0.0 0.0 0.0 la
 0.0 0.0 1.0 o
$end
'''
DEFINE = '''a coord
ired
*
no
b "la" def2-SVP
*
eht
-1
u 0
*
'''


def generate_fixture(root):
    root = Path(root)
    for name, contents in [('control', MASTER), ('coord', COORD), ('def.inp', DEFINE),
                            ('def_gd.inp', DEFINE)]:
        (root / name).write_text(contents)
    # Model a finished metal control, with proper TurboMole group references.
    metal = root / 'la'
    metal.mkdir()
    ctl = MASTER.replace('$title', '$title\n$coord file=coord\n$atoms\n   la 1 \n   o 2')
    ctl = ctl.replace('$fermi tmstrt=1000 nue= 6.00\n', '')
    (metal / 'control').write_text(ctl)
    (metal / 'coord').write_text(COORD)
    (metal / 'basis').write_text('$basis\n*\n$ecp\n*\n$end\n')
    (metal / 'auxbasis').write_text('$jbas\n*\n$end\n')
    (metal / 'alpha').write_text('$uhfmo_alpha\n$end\n')
    (metal / 'beta').write_text('$uhfmo_beta\n$end\n')
    (metal / 'def.inp').write_text(DEFINE)
    (metal / 'metal.json').write_text(json.dumps({'metal': 'la', 'unpaired': 0, 'charge': -1}))
    return metal


def digest_tree(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in Path(root).rglob('*') if p.is_file()}


class ValidateTests(unittest.TestCase):
    def test_master_does_not_require_uncreated_basis_files(self):
        with tempfile.TemporaryDirectory() as d:
            generate_fixture(d)
            report, ctl = val.check_project(Path(d), FakeControl, ['la'])
            self.assertEqual(report.status, 'PASS', report.issues)
            self.assertIsNotNone(ctl)

    def test_valid_finished_metal(self):
        with tempfile.TemporaryDirectory() as d:
            generate_fixture(d)
            root = Path(d)
            project, master = val.check_project(root, FakeControl, ['la'])
            report, child_ctl = val.check_metal(root, 'la', master, FakeControl)
            self.assertEqual(project.status, 'PASS', project.issues)
            self.assertEqual(report.status, 'PASS', report.issues)
            self.assertIsNotNone(child_ctl)

    def test_missing_basis_reference_is_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            (metal / 'basis').unlink()
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'FAIL')
            self.assertIn('basis', str(report.issues))

    def test_wrong_functional_or_dispersion_fails(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            p = metal / 'control'
            p.write_text(p.read_text().replace('functional tpssh', 'functional pbe').replace('$disp4\n', ''))
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'FAIL')
            self.assertIn('functional', str(report.issues))
            self.assertIn('$disp4', str(report.issues))

    def test_wrong_metal_and_charge_fail(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            (metal / 'coord').write_text(COORD.replace('la\n', 'pr\n'))
            (metal / 'metal.json').write_text(json.dumps({'metal': 'pr', 'unpaired': 2, 'charge': 7}))
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'FAIL')
            self.assertIn('expected exactly one la', str(report.issues))

    def test_screwer_does_not_compare_mo_coefficients(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            child = metal / 'screwer_la'
            child.mkdir()
            for name in ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta'):
                (child / name).write_bytes((metal / name).read_bytes())
            (child / 'alpha').write_text('DEFINE legitimately changed the orbitals')
            (child / 'screwer_prepared.json').write_text('{}')
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            parent_report, parent_ctl = val.check_metal(Path(d), 'la', master, FakeControl)
            report = val.check_screwer(Path(d), 'la', metal, parent_ctl, FakeControl)
            self.assertEqual(parent_report.status, 'PASS', parent_report.issues)
            self.assertEqual(report.status, 'PASS', report.issues)

    def test_screwer_detects_method_drift(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            child = metal / 'screwer_la'
            child.mkdir()
            for name in ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta'):
                (child / name).write_bytes((metal / name).read_bytes())
            (child / 'control').write_text((child / 'control').read_text().replace('functional tpssh', 'functional pbe0'))
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            _, parent_ctl = val.check_metal(Path(d), 'la', master, FakeControl)
            report = val.check_screwer(Path(d), 'la', metal, parent_ctl, FakeControl)
            self.assertEqual(report.status, 'FAIL')

    def test_no_writes_anywhere(self):
        with tempfile.TemporaryDirectory() as d:
            generate_fixture(d)
            before = digest_tree(d)
            root = Path(d)
            _, master = val.check_project(root, FakeControl, ['la'])
            _, parent_ctl = val.check_metal(root, 'la', master, FakeControl)
            val.check_screwer(root, 'la', root / 'la', parent_ctl, FakeControl)
            after = digest_tree(d)
            self.assertEqual(before, after)

    def test_cartesian_nan_fails(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            (metal / 'coord').write_text(COORD.replace('0.0 0.0 1.0', 'nan 0.0 1.0'))
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'FAIL')


# Additional cases for the actual lifecycle boundaries.
class PreparationPhaseTests(unittest.TestCase):
    def test_fresh_setup_is_not_incomplete_define(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            (metal / 'control').write_text(MASTER)
            for f in ('basis', 'auxbasis', 'alpha', 'beta'):
                (metal / f).unlink()
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'PASS', report.issues)
            self.assertIn('prepared before DEFINE', str(report.issues))

    def test_incomplete_define_is_not_prepared(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            (metal / 'control').write_text(MASTER)
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            report, _ = val.check_metal(Path(d), 'la', master, FakeControl)
            self.assertEqual(report.status, 'FAIL')
            self.assertIn('incomplete existing DEFINE', str(report.issues))

    def test_screwer_preparation_is_not_incomplete_define(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            child = metal / 'screwer_la'
            child.mkdir()
            for name in ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta'):
                (child / name).write_bytes((metal / name).read_bytes())
            (child / 'control').write_text((child / 'control').read_text().replace('$coord file=coord\n', ''))
            (child / 'screwer_prepared.json').write_text('{}')
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            _, parent_ctl = val.check_metal(Path(d), 'la', master, FakeControl)
            report = val.check_screwer(Path(d), 'la', metal, parent_ctl, FakeControl)
            self.assertEqual(report.status, 'PASS', report.issues)

    def test_screwer_failed_define_missing_coord_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            metal = generate_fixture(d)
            child = metal / 'screwer_la'
            child.mkdir()
            for name in ('control', 'coord', 'basis', 'auxbasis', 'alpha', 'beta'):
                (child / name).write_bytes((metal / name).read_bytes())
            (child / 'control').write_text((child / 'control').read_text().replace('$coord file=coord\n', ''))
            (child / 'screwer_prepared.json').write_text('{}')
            (child / 'define.out').write_text('error in cartesian coordinates')
            _, master = val.check_project(Path(d), FakeControl, ['la'])
            _, parent_ctl = val.check_metal(Path(d), 'la', master, FakeControl)
            report = val.check_screwer(Path(d), 'la', metal, parent_ctl, FakeControl)
            self.assertEqual(report.status, 'FAIL')


if __name__ == '__main__':
    unittest.main()
