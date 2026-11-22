"""Regression tests for user-edited restricted DEFINE transcripts."""
import contextlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import closed_shell_workflow as closed


class RestrictedManualDefine(unittest.TestCase):
    @staticmethod
    def setup(root):
        (root / 'coord').write_text('$coord\n 0 0 0 gd f\n 0 0 1 h\n$end\n')
        (root / 'control').write_text('$title\n$uhf\n$fermi nue=7.0\n$dft\n functional tpssh\n$end\n')
        (root / 'def_closedshell.inp').write_text(
            'a coord\n*\nno\nb all def2-SVP\nb "gd" def2-TZVP\n*\neht\n0\n\nri\n\nq\n')
        closed.prepare(root, 'la')
        return root / 'la_closedshell'

    def test_edited_define_warns_but_submits_without_overwriting(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder = self.setup(base)
            def_file = folder / 'def.inp'
            user_edit = def_file.read_text().replace('def2-TZVP', 'def2-QZVP')
            def_file.write_text(user_edit)
            calls = []

            def runner(cmd, **kw):
                calls.append(cmd)
                if cmd[0] == 'squeue':
                    return subprocess.CompletedProcess(cmd, 0, '', '')
                return subprocess.CompletedProcess(cmd, 0, str(900 + sum(x[0] == 'sbatch' for x in calls)) + '\n', '')

            output = io.StringIO()
            with mock.patch.object(closed.shutil, 'which', return_value='/usr/bin/mock'), \
                 mock.patch.object(closed.subprocess, 'run', side_effect=runner), \
                 contextlib.redirect_stdout(output):
                closed.submit(base, 'la', 'none')
            self.assertIn('WARNING', output.getvalue())
            self.assertIn('def.inp was edited', output.getvalue())
            self.assertEqual(user_edit, def_file.read_text())
            self.assertEqual(2, sum(x[0] == 'sbatch' for x in calls))
            self.assertEqual(['restricted_ridft', 'jobex_eiger'],
                             (folder / 'chain_plan.txt').read_text().splitlines())

    def test_invalid_unrestricted_command_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder = self.setup(base)
            d = folder / 'def.inp'
            d.write_text(d.read_text() + 'u 0\n')
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'unrestricted occupation'):
                    closed.pristine(folder, base, 'la')

    def test_changed_charge_or_wrong_metal_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder = self.setup(base)
            d = folder / 'def.inp'
            original = d.read_text()
            d.write_text(original.replace('eht\n0', 'eht\n-2'))
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'charge disagrees'):
                    closed.pristine(folder, base, 'la')
            d.write_text(original.replace('b "la"', 'b "lu"'))
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'exactly one basis for la'):
                    closed.pristine(folder, base, 'la')

    def test_master_define_drift_does_not_reject_prepared_transcript(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder = self.setup(base)
            (base / 'def_closedshell.inp').write_text('changed master template (not used for existing prepared input)')
            closed.pristine(folder, base, 'la')


if __name__ == '__main__':
    unittest.main()
