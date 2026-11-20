"""Regression tests for fixed-atom TURBOMOLE coordinates (no cluster needed)."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import setup_lanthanides as setup
import submission_guard as guard
import closed_shell_workflow as closed


COORD = ("$coord\n"
         "   -2.12915379847462 -11.65202040293635 4.19916256966488 h\n"
         "   -1.91839326993634 -10.84402329081377 2.29230392958542 c f\n"
         "    5.85765020891686 -7.11363400093588 4.26891021279086 c f\n"
         "    0.98944548104506 -7.73053131486772 4.23397295619878 n f\n"
         "    1.49251514186006 -0.41838347426602 -1.81193554870555 gd f\n"
         "$end\n")


class ConstrainedMasterCoordinates(unittest.TestCase):
    def test_replace_fixed_metal_preserves_exact_non_element_bytes(self):
        expected = COORD.replace('gd f', 'la f')
        self.assertEqual(expected, setup.make_coord(COORD, 'la'))
        self.assertEqual(4, setup.make_coord(COORD, 'ce').count(' f\n'))
        self.assertEqual(COORD.replace('gd f', 'yb f'), setup.make_coord(COORD, 'yb'))

    def test_unconstrained_coordinates_unchanged_beyond_metal_name(self):
        plain = COORD.replace(' f\n', '\n')
        self.assertEqual(plain.replace('gd\n', 'nd\n'), setup.make_coord(plain, 'nd'))

    def test_flags_comments_and_other_groups_preserved(self):
        source = ('$coord angs\r\n'
                  '  1.0\t2.0 3.0\tGD\tf   \r\n'
                  '  2 3 4 h f\r\n'
                  '$redundant\r\n'
                  '  0 0 0 la f\r\n'
                  '$end\r\n')
        expected = source.replace('GD\tf', 'lu\tf')
        self.assertEqual(expected, setup.make_coord(source, 'lu'))

    def test_reject_missing_or_two_lanthanide_atoms(self):
        with self.assertRaisesRegex(ValueError, 'found 0'):
            setup.make_coord('$coord\n 0 0 0 h f\n$end\n', 'la')
        with self.assertRaisesRegex(ValueError, 'found 2'):
            setup.make_coord('$coord\n 0 0 0 gd f\n 1 0 0 lu f\n$end\n', 'la')

    def test_real_setup_and_submission_guard_accept_fixed_metal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'coord').write_text(COORD)
            (root / 'control').write_text('$title\n$fermi tmstrt=1000 nue=7.00\n$end\n')
            (root / 'def.inp').write_text('a coord\n*\nno\nb "gd" def2-TZVP\n*\neht\n-1\nno\nu 7\n*\nri\nq\n')
            with mock.patch.object(setup.shutil, 'which', return_value=None):
                self.assertTrue(setup.setup_one(root, 'la'))
            generated = (root / 'la' / 'coord').read_text()
            self.assertEqual(COORD.replace('gd f', 'la f'), generated)
            self.assertEqual(-1, json.loads((root / 'la' / 'metal.json').read_text())['charge'])
            # Existing guard must recognize the unchanged f flags as pristine.
            guard.verify_prepared(root, 'la')

    def test_restricted_preparation_uses_same_fixed_atom_logic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'coord').write_text('$coord\n0 0 0 gd f\n0 0 1 h\n$end\n')
            (root / 'control').write_text('$title\n$uhf\n$fermi nue=0.00\n$dft\n functional tpssh\n$end\n')
            (root / 'def_closedshell.inp').write_text('a coord\n*\nno\nb "gd" def2-TZVP\n*\neht\n0\nri\nq\n')
            files, charge, electrons = closed.inputs(root, 'la')
            self.assertEqual(0, charge)
            self.assertEqual(58, electrons)
            self.assertIn('la f\n', files['coord'])



if __name__ == '__main__':
    unittest.main()
