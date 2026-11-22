"""Offline behavioral contracts; no TURBOMOLE, VACC or SLURM dependencies."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import adaptive_define as adaptive
import closed_shell_workflow as closed
import recover_universal as recovery
import validate_lanthanides as validator
from test_validate_lanthanides import FakeControl, MASTER, COORD, DEFINE

SCF_OK = '''energy converged !
convergence criteria satisfied after 22 iterations
ridft : all done
'''


def populate(parent, metal='nd', *, screwer=False, method='aoforce'):
    path = parent / metal
    if screwer:
        path = path / f'screwer_{metal}'
    path.mkdir(parents=True)
    (path / 'control').write_text('$dft\n functional tpssh\n$end\n')
    (path / 'coord').write_text('$coord\n 0 0 0 nd\n$end\n')
    (path / 'tm_workflow.py').write_text('# supplied helper\n')
    (path / 'metal.json').write_text(json.dumps({'metal': metal, 'unpaired': 3, 'charge': -1}))
    (path / 'chain_plan.txt').write_text('screwer_jobex_eiger\n' + method + '\n' if screwer else
                                         'ridft_eiger_fix\ninteger_ridft\njobex_eiger\n' + method + '\n')
    if screwer:
        (path / 'screwer_prepared.json').write_text(json.dumps({'method': method, 'metal': metal}))
    return path


class AdaptiveDefineTests(unittest.TestCase):
    def test_pexpect_supervisor_detects_ired_error(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            executable = root / 'define'
            executable.write_text('#!/usr/bin/env python3\nimport sys\nfor line in sys.stdin:\n if line.strip()=="ired":\n  print("SINCE THERE ARE THREE COLLINEAR BOND PARTNERS",flush=True)\n  break\n')
            executable.chmod(0o755)
            with mock.patch.dict(os.environ, {'PATH': str(root) + os.pathsep + os.environ['PATH']}):
                result = adaptive.define_supervised(root, 'a coord\nired\n*\nq\n')
            self.assertNotEqual(0, result.returncode)
            self.assertIn('COLLINEAR', result.stdout)

    def test_pexpect_supervisor_accepts_normal_define(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            executable = root / 'define'
            executable.write_text('#!/usr/bin/env python3\nimport sys\nfor line in sys.stdin:\n if line.strip()=="q":\n  print(" **** define : all done ****", flush=True)\n  print("define ended normally",flush=True)\n  break\n')
            executable.chmod(0o755)
            with mock.patch.dict(os.environ, {'PATH': str(root) + os.pathsep + os.environ['PATH']}):
                result = adaptive.define_supervised(root, 'a coord\n*\nno\nq\n')
            self.assertEqual(0, result.returncode, result.stdout)
            self.assertIn('define ended normally', result.stdout)


class ValidatorRegressionTests(unittest.TestCase):
    def test_repeated_turbomole_result_groups_are_not_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            for name, value in [('coord', COORD), ('def.inp', DEFINE), ('def_gd.inp', DEFINE)]:
                (root / name).write_text(value)
            original = MASTER.replace('$end\n', '$last step ridft\n$last step vibrations\n'
                                      '$vibrational spectrum\n 1 0.0\n$vibrational spectrum\n 1 0.0\n'
                                      '$nprvibrational spectrum\n 1 0.0\n$nprvibrational spectrum\n 1 0.0\n$end\n')
            (root / 'control').write_text(original)
            report, _ = validator.check_project(root, FakeControl, ['la'])
            self.assertEqual('PASS', report.status, report.issues)


class ClosedShellTests(unittest.TestCase):
    @staticmethod
    def master(root):
        (root / 'control').write_text(MASTER + '$alpha shells\n a 1-123 (1)\n$beta shells\n a 1-118 (1)\n')
        (root / 'coord').write_text('$coord\n 0 0 0 la\n 0 0 1 h\n$end\n')
        (root / 'def_closedshell.inp').write_text('a coord\n*\nno\nb "la" def2-TZVP\n*\neht\n0\ny\n*\nq\n')

    def test_restricted_removes_uhf_groups_without_touching_dft(self):
        out = closed.restricted_control(MASTER + '$alpha shells\n a 1-124\n$beta shells\n a 1-117\n')
        for group in ('$uhf\n', '$uhfmo_alpha', '$uhfmo_beta', '$alpha shells', '$beta shells', '$fermi'):
            self.assertNotIn(group, out)
        for group in ('$dft', '$disp4', '$scfconv', '$ecp'):
            self.assertIn(group, out)
        self.assertIn('functional tpssh', out)

    def test_requires_user_template_without_unrestricted_occupation(self):
        with self.assertRaisesRegex(ValueError, 'unrestricted'):
            closed.restricted_define('b "la" def2-SVP\neht\n0\nu 0\n', 'lu')
        text = closed.restricted_define('b "la" def2-SVP\neht\n0\ny\n', 'lu')
        self.assertIn('b "lu"', text)

    def test_electron_parity_rejects_odd(self):
        with self.assertRaisesRegex(ValueError, 'not a positive even'):
            closed.electron_count('$coord\n 0 0 0 la\n$end\n', 0)
        self.assertEqual(58, closed.electron_count('$coord\n 0 0 0 la\n 0 0 1 h\n$end\n', 0))

    def test_fresh_closed_shell_is_preserved_and_verified(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.master(root)
            self.assertTrue(closed.prepare(root, 'lu'))
            folder = root / 'lu_closedshell'
            self.assertIn('$dft', (folder / 'control').read_text())
            self.assertNotIn('$uhf', (folder / 'control').read_text())
            self.assertEqual(72, json.loads((folder / 'metal.json').read_text())['electron_count_all_electron'])
            closed.pristine(folder, root, 'lu')
            self.assertFalse(closed.prepare(root, 'lu'))
            (folder / 'control').write_text((folder / 'control').read_text() + 'tampered')
            with self.assertRaisesRegex(ValueError, 'changed after preparation'):
                closed.pristine(folder, root, 'lu')

    def test_restricted_slurm_dependency_chain(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.master(root)
            closed.prepare(root, 'la')
            calls = []
            def fake_run(command, **kwargs):
                calls.append(command)
                if command[0] == 'squeue':
                    return subprocess.CompletedProcess(command, 0, '', '')
                return subprocess.CompletedProcess(command, 0, f'{701 + len([x for x in calls if x[0] == "sbatch"])}\n', '')
            with mock.patch.object(closed.shutil, 'which', return_value='/usr/bin/fake'), mock.patch.object(closed.subprocess, 'run', side_effect=fake_run):
                closed.submit(root, 'la', 'aoforce')
            stages = (root / 'la_closedshell' / 'chain_plan.txt').read_text().splitlines()
            self.assertEqual(['restricted_ridft', 'jobex_eiger', 'aoforce'], stages)
            sbatches = [cmd for cmd in calls if cmd[0] == 'sbatch']
            self.assertIn('--dependency=afterok:702', sbatches[1])
            self.assertIn('--dependency=afterok:703', sbatches[2])
            self.assertFalse((root / 'la_closedshell' / '.lan_submission_pending').exists())


class UniversalRecoveryTests(unittest.TestCase):
    def test_resume_post_occupation_fix_starts_integer_ridft_not_second_remove_fermi(self):
        # Reproduces EF1741: ridft_eiger_fix already removed $fermi,
        # subsequent recovery mistakenly called remove-fermi a second time.
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), metal='la', method='numforce')
            (folder / 'ridft.out').write_text(SCF_OK)
            (folder / 'initial_occupations').write_text('$alpha shells\n a 1-2 (1)\n$beta shells\n a 1-2 (1)\n')
            (folder / 'control').write_text('$dft\n functional tpssh\n$alpha shells\n a 1-2 (1)\n$beta shells\n a 1-2 (1)\n$end\n')
            (folder / 'control.with_fermi').write_text('$dft\n functional tpssh\n$fermi nue=0.00\n$alpha shells\n a 1-2 (1)\n$beta shells\n a 1-2 (1)\n$end\n')
            control_before = (folder / 'control').read_bytes()
            stages, reason = recovery.stage_plan(folder, 'standard', 'numforce')
            self.assertEqual(['integer_ridft', 'jobex_eiger', 'numforce'], stages)
            self.assertIn('$fermi already removed', reason)
            self.assertEqual(control_before, (folder / 'control').read_bytes())

    def test_fermi_present_still_runs_occupation_fix(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), method='numforce')
            (folder / 'initial_occupations').write_text('saved occupied shells')
            (folder / 'ridft.out').write_text(SCF_OK)
            (folder / 'control').write_text('$fermi nue=3.00\n$end\n')
            stages, _ = recovery.stage_plan(folder, 'standard', 'numforce')
            self.assertEqual(['eiger_fix_integer_ridft', 'jobex_eiger', 'numforce'], stages)

    def test_missing_fermi_without_saved_prior_control_requires_review(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), method='numforce')
            (folder / 'initial_occupations').write_text('saved occupied shells')
            (folder / 'ridft.out').write_text(SCF_OK)
            with self.assertRaisesRegex(recovery.RecoveryError, 'control.with_fermi'):
                recovery.stage_plan(folder, 'standard', 'numforce')

    def test_standard_jobex_failure_skips_completed_scf(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t))
            (folder / 'ridft.integer.out').write_text(SCF_OK)
            (folder / 'not.converged').write_text('stale failure')
            method = recovery.selected_method(folder, 'standard', 'auto')
            stages, reason = recovery.stage_plan(folder, 'standard', method)
            self.assertEqual(['jobex_eiger', 'aoforce'], stages)
            self.assertIn('JOBEX', reason)
            self.assertTrue((folder / 'not.converged').exists(), 'planning must not delete input/output files')

    def test_screwer_jobex_failure_has_no_define_stage(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            (folder / 'define.fallback.out').write_text('define ended normally\n')
            (folder / 'ridft.out').write_text(SCF_OK)
            (folder / 'not.converged').write_text('old error')
            method = recovery.selected_method(folder, 'screwer', 'auto')
            stages, _ = recovery.stage_plan(folder, 'screwer', method)
            self.assertEqual(['jobex_eiger', 'aoforce'], stages)
            self.assertNotIn('screwer_jobex_eiger', stages)

    def test_screwer_scf_failure_restarts_without_define(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            (folder / 'define.out').write_text('define ended normally\n')
            stages, _ = recovery.stage_plan(folder, 'screwer', 'aoforce')
            self.assertEqual(['screwer_ridft_jobex_eiger', 'aoforce'], stages)
            (folder / 'ridft.resume.101.out').write_text(SCF_OK)
            stages, _ = recovery.stage_plan(folder, 'screwer', 'aoforce')
            self.assertEqual(['jobex_eiger', 'aoforce'], stages)

    def test_screwer_failed_define_is_not_repeated(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            (folder / 'define.out').write_text('ERROR IN CARTESIAN COORDINATES ?\n')
            with self.assertRaisesRegex(recovery.RecoveryError, 'not verified'):
                recovery.stage_plan(folder, 'screwer', 'aoforce')

    def test_hessian_failure_restarts_hessian_only(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            (folder / 'GEO_OPT_CONVERGED').touch()
            (folder / 'gradient').write_text('good gradient')
            (folder / 'hessian').write_text('old hessian')
            (folder / 'vibspectrum').write_text('old spectrum')
            stages, reason = recovery.stage_plan(folder, 'screwer', 'aoforce')
            self.assertEqual(['aoforce'], stages)
            self.assertIn('Hessian', reason)

    def test_completed_hessian_is_not_scheduled_again(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            (folder / 'GEO_OPT_CONVERGED').touch()
            (folder / 'gradient').write_text('converged gradient')
            (folder / 'jobex.out').write_text('optimization done')
            (folder / 'hessian').write_text('generated Hessian')
            (folder / 'vibspectrum').write_text('$vibrational spectrum\n 1 0.00\n$end')
            (folder / 'aoforce.out').write_text('frequency analysis completed with results')
            stages, _ = recovery.stage_plan(folder, 'screwer', 'aoforce')
            self.assertEqual([], stages)

    def test_manual_converged_geometry_can_recover_hessian_without_guessing_charge(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), metal='gd')
            (folder / 'metal.json').unlink()
            (folder / 'GEO_OPT_CONVERGED').touch()
            (folder / 'gradient').write_text('existing gradient')
            method = recovery.selected_method(folder, 'standard', 'ao')
            self.assertEqual('aoforce', method)
            stages, _ = recovery.stage_plan(folder, 'standard', method)
            self.assertEqual(['aoforce'], stages)
            (folder / 'GEO_OPT_CONVERGED').unlink()
            with self.assertRaisesRegex(recovery.RecoveryError, 'Hessian-only'):
                recovery.stage_plan(folder, 'standard', method)

    def test_bad_frequency_override_does_not_change_screwer_method(self):
        with tempfile.TemporaryDirectory() as t:
            folder = populate(Path(t), screwer=True)
            with self.assertRaisesRegex(recovery.RecoveryError, 'method'):
                recovery.selected_method(folder, 'screwer', 'num')

    def test_detects_paths_in_both_stages_and_restricted(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.assertEqual('standard', recovery.identify(root / 'nd')[2])
            self.assertEqual('screwer', recovery.identify(root / 'nd' / 'screwer_nd')[2])
            folder, metal, kind = recovery.identify(root / 'lu_closedshell')
            self.assertEqual(('lu', 'standard'), (metal, kind))

    def test_submit_recovery_only_new_stages_with_afterok(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            folder = populate(root, screwer=True)
            calls = []
            def fake_run(command, **kwargs):
                calls.append(command)
                if command[0] == 'squeue':
                    return subprocess.CompletedProcess(command, 0, '', '')
                return subprocess.CompletedProcess(command, 0, f'{900 + len([x for x in calls if x[0] == "sbatch"])}\n', '')
            with mock.patch.object(recovery, 'shutil_which', return_value='/usr/bin/fake'), mock.patch.object(recovery.subprocess, 'run', side_effect=fake_run):
                journal = recovery.submit(folder, 'nd', 'screwer', 'aoforce', ['jobex_eiger', 'aoforce'], Path(recovery.__file__).parent)
            self.assertTrue(journal.exists())
            self.assertFalse((folder / '.lan_recovery_pending').exists())
            self.assertFalse((folder / '.lan_recovery_lock').exists())
            self.assertIn('afterok:901', journal.read_text())
            sbatches = [x for x in calls if x[0] == 'sbatch']
            self.assertTrue(any(x.startswith('--export=ALL,TM_STAGE=jobex_eiger,TM_RECOVERY=1,TM_SCREWER_CHILD=1') for x in sbatches[0]))
            self.assertTrue(any(x == '--dependency=afterok:901' for x in sbatches[1]))

    def test_active_scheduler_job_blocks_submission(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            folder = populate(root, screwer=True)
            fake = subprocess.CompletedProcess(['squeue'], 0, f'923|{folder}|DependencyNeverSatisfied\n', '')
            with mock.patch.object(recovery, 'shutil_which', return_value='/bin/squeue'), mock.patch.object(recovery.subprocess, 'run', return_value=fake):
                with self.assertRaisesRegex(recovery.RecoveryError, 'cancel obsolete'):
                    recovery.scheduler_clear(folder)

if __name__ == '__main__':
    unittest.main()
