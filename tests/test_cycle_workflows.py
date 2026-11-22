"""Development regression checks; mocked scheduling, no real TURBOMOLE."""
import contextlib
import io
import json
import sys
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest import mock

import recover_universal as recovery
import screwer_workflow as screwer
import tm_workflow as tm
from test_development_workflows import populate


class RecoveryPlanTests(unittest.TestCase):
    def test_fermi_is_warning_not_fatal_to_plan_and_mock_submission(self):
        with tempfile.TemporaryDirectory() as t:
            child = populate(Path(t), screwer=True)
            (child / 'define.fallback.out').write_text('define ended normally\n')
            (child / 'ridft.out').write_text('convergence criteria satisfied after 1 iterations\nENERGY CONVERGED !\nridft : all done')
            (child / 'control').write_text((child / 'control').read_text() + '\n$fermi tmstrt=1000 nue= 5.00\n')
            stages, _ = recovery.stage_plan(child, 'screwer', 'aoforce')
            self.assertEqual(['jobex_eiger', 'aoforce'], stages)
            warns = recovery.recovery_warnings(child, 'screwer')
            self.assertEqual(1, len(warns))
            self.assertIn('$fermi', warns[0])
            self.assertIn('$fermi', (child / 'control').read_text())
            with mock.patch.object(recovery.Path, 'cwd', return_value=child), \
                 mock.patch.object(recovery, 'submit', return_value=child / 'mock-journal') as submit:
                with contextlib.redirect_stdout(io.StringIO()) as stdout:
                    self.assertEqual(0, recovery.main(['--submit']))
                self.assertIn('WARNING: Screwer contains $fermi', stdout.getvalue())
                submit.assert_called_once()

    def test_completed_hessian_negative_does_not_resubmit_same_hessian(self):
        with tempfile.TemporaryDirectory() as t:
            child = populate(Path(t), screwer=True)
            (child / 'GEO_OPT_CONVERGED').touch()
            (child / 'gradient').write_text('gradient')
            (child / 'jobex.out').write_text('optimized')
            (child / 'hessian').write_text('computed')
            (child / 'vibspectrum').write_text('$vibrational spectrum\n 1 a -2.50 0.0 0.0\n$end\n')
            (child / 'aoforce.out').write_text('done')
            stages, reason = recovery.stage_plan(child, 'screwer', 'aoforce')
            self.assertEqual([], stages)
            self.assertIn('negative mode', reason)
            self.assertIn('--cycle', reason)

    def test_nested_cycle_path_identified(self):
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            path = screwer.cycle_directory(base, 'lu', 3)
            self.assertEqual('screwer_lu_cycle3', path.name)
            self.assertEqual(('lu', 'screwer'), recovery.identify(path)[1:])
            self.assertEqual(path.resolve(), recovery.identify(base, base, 'lu/screwer_lu/screwer_lu_cycle2/screwer_lu_cycle3')[0])


class ScrewerCoordinateTests(unittest.TestCase):
    def test_parent_coord_screwer_comment_is_not_an_atom(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / 'coord'
            atoms = [f'    {i}.0 0.0 0.0  c' for i in range(56)]
            parent.write_text('$coord\n# cartesian coordinates shifted along normal mode 1 by *******\n'
                              + '\n'.join(atoms) + '\n$end\n')
            displaced = '# cartesian coordinates shifted along normal mode 1 by *******\n' + '\n'.join(atoms)
            result = screwer.parse_coordinates(displaced, parent)
            self.assertTrue(result.startswith('$coord\n'))
            self.assertEqual(56, len(result.splitlines()) - 2)
            self.assertEqual('c', result.splitlines()[1].split()[3])


    def test_parent_coord_redundant_group_is_not_atom_data(self):
        # TURBOMOLE's coord can contain thousands of non-Cartesian lines
        # after $coord; these must not become apparent extra "atoms".
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / 'coord'
            atoms = ['    0.0 0.0 0.0  h', '    1.0 2.0 3.0  nd']
            parent.write_text('$coord\n' + '\n'.join(atoms) +
                              '\n$user-defined bonds\n$redundant\n'
                              '    1 k 1.00000 stre 1 2 val=1.0\n'
                              '    2 k 1.00000 tors 1 2 1 2 val=2.0\n$end\n')
            displaced = '# cartesian coordinates shifted\n' + '\n'.join(atoms)
            result = screwer.parse_coordinates(displaced, parent)
            self.assertEqual(2, len(result.splitlines()) - 2)
            self.assertIn(' nd\n', result)

    def test_parent_coord_redundant_group_still_rejects_wrong_atom_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / 'coord'
            parent.write_text('$coord\n 0 0 0 h\n 1 1 1 nd\n'
                              '$redundant\n 1 k 1.0 tors 1 2 1 2\n$end\n')
            with self.assertRaisesRegex(RuntimeError, 'atom count and order'):
                screwer.parse_coordinates(' 0 0 0 nd\n 1 1 1 h\n', parent)


class ScrewerFermiTests(unittest.TestCase):
    def test_inherited_fermi_warns_but_retains_orbitals_and_shells(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            (root / 'alpha').write_text('$uhfmo_alpha\n 1 a eigenvalue=-1.0\n$end\n')
            (root / 'beta').write_text('$uhfmo_beta\n 1 a eigenvalue=-1.0\n$end\n')
            control = ('$uhf\n$alpha shells\n a 1-124 (1)\n'
                       '$beta shells\n a 1-117 (1)\n'
                       '$fermi tmstrt= 1000.00 tmend=100.00 nue= 3.00\n$end\n')
            (root / 'control').write_text(control)
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                result = tm.screwer_inheritance(root)
            self.assertIn('WARNING', stdout.getvalue())
            self.assertIn('$fermi', stdout.getvalue())
            self.assertEqual({'alpha', 'beta'}, set(result['orbitals_sha256']))
            self.assertEqual(2, len(result['shell_groups']))
            self.assertEqual(control, (root / 'control').read_text())

    def test_no_fermi_emits_no_warning(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            (root / 'mos').write_text('$scfmo\n 1 a eigenvalue=-1.0\n$end\n')
            (root / 'control').write_text('$closed shells\n a 1-3 (2)\n$end\n')
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                result = tm.screwer_inheritance(root)
            self.assertNotIn('WARNING', stdout.getvalue())
            self.assertEqual({'mos'}, set(result['orbitals_sha256']))

    def test_second_cycle_preparation_with_user_fermi_writes_manifest(self):
        # Exercise the actual prepare -> cpc -> screwer -> sdg -> kdg ->
        # inheritance -> preparation-record path using stand-in executables.
        from types import SimpleNamespace
        import shutil
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            parent = screwer.cycle_directory(base, 'nd', 1)
            parent.mkdir(parents=True)
            (base / 'def_screwer.inp').write_text('a coord\n*\nno\nq\n')
            (base / 'def_screwer_ln.inp').write_text('a coord\n*\nno\nq\n')
            (base / 'screwer.inp').write_text('1\n')
            (parent / 'metal.json').write_text(json.dumps(dict(metal='nd', charge=-1, unpaired=3)))
            (parent / 'control').write_text(
                '$alpha shells\n a 1-124 (1)\n$beta shells\n a 1-117 (1)\n'
                '$fermi tmstrt=1000 nue=3.00\n$end\n')
            (parent / 'coord').write_text('$coord\n 0 0 0 nd\n 1 1 1 o\n$end\n')
            (parent / 'alpha').write_text('$uhfmo_alpha\n1 a eigenvalue=-1.0\n$end\n')
            (parent / 'beta').write_text('$uhfmo_beta\n1 a eigenvalue=-1.0\n$end\n')

            def simulate(command, *, cwd, input_text=None, timeout=180):
                if command[0] == 'cpc':
                    child = cwd / command[1]
                    child.mkdir()
                    for path in cwd.iterdir():
                        if path.is_file():
                            shutil.copy2(path, child / path.name)
                    return SimpleNamespace(returncode=0, stdout='')
                if command[0] == 'screwer':
                    control = cwd / 'control'
                    ctl = control.read_text()
                    control.write_text(ctl + '\n$newcoord\n'
                                       '# cartesian coordinates shifted along normal mode 1 by 0.1\n'
                                       ' 0 0 0 nd\n 1 1 1 o\n$end\n')
                    return SimpleNamespace(returncode=0, stdout='screwer complete')
                if command[:3] == ['sdg', '-b', 'newcoord']:
                    return SimpleNamespace(returncode=0,
                        stdout='# shifted along normal mode 1\n 0 0 0 nd\n 1 1 1 o\n')
                if command[0] == 'kdg':
                    return SimpleNamespace(returncode=0, stdout='')
                raise AssertionError(f'unexpected command: {command}')

            with mock.patch.object(screwer, 'verify_parent', return_value='automation'), \
                 mock.patch.object(screwer, 'spectrum_modes', return_value=[(1, -15.93)]), \
                 mock.patch.object(screwer, 'execute', side_effect=simulate), \
                 mock.patch.object(screwer.shutil, 'which', return_value='/fake/bin'), \
                 contextlib.redirect_stdout(io.StringIO()) as stdout:
                child = screwer.prepare(base, 'nd', 'aoforce', cycle=2, max_cycles=2)
            self.assertTrue((child / 'screwer_prepared.json').is_file())
            self.assertIn('$fermi', (child / 'control').read_text())
            self.assertIn('WARNING', stdout.getvalue())
            self.assertEqual(2, json.loads((child / 'screwer_prepared.json').read_text())['cycle'])

    def test_fermi_does_not_bypass_missing_orbitals(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            (root / 'control').write_text('$alpha shells\n a 1-3 (1)\n$fermi nue=3\n$end\n')
            with self.assertRaisesRegex(RuntimeError, 'orbital files missing'):
                tm.screwer_inheritance(root)


class AdaptiveComputePythonTests(unittest.TestCase):
    def prepared_child(self, base, backend='adaptive'):
        child = base / 'pr' / 'screwer_pr'
        child.mkdir(parents=True)
        (child / 'screwer_prepared.json').write_text(json.dumps(dict(
            metal='pr', method='aoforce', workflow_version=2,
            inherited_state=True, define_backend=backend,
        )))
        for name in ('control', 'coord', 'def_screwer.inp',
                     'tm_workflow.py', 'setup_lanthanides.py'):
            (child / name).write_text('test\n')
        (child / 'turbo.pipeline.slurm').write_text(
            'screwer continuation SCF before JOBEX\n--define-fallback\n'
        )
        return child

    def test_adaptive_submission_exports_absolute_python_to_both_jobs(self):
        with tempfile.TemporaryDirectory() as t:
            child = self.prepared_child(Path(t))
            calls = []

            def fake_execute(command, *, cwd, **kwargs):
                self.assertEqual(child, cwd)
                calls.append(command)
                if command[0] == 'sbatch':
                    return SimpleNamespace(returncode=0, stdout=str(7700000 + len(calls)))
                return SimpleNamespace(returncode=0, stdout='')

            with mock.patch.object(screwer, 'execute', side_effect=fake_execute), \
                 mock.patch.object(screwer.shutil, 'which', return_value='/usr/bin/sbatch'), \
                 mock.patch.dict(screwer.os.environ, {'TM_WORKFLOW_PYTHON': sys.executable}):
                self.assertEqual(2, screwer.submit(child, 'aoforce'))
            self.assertEqual([sys.executable, '-c', 'import pexpect'], calls[0])
            self.assertEqual(3, len(calls))
            for command in calls[1:]:
                export = next(x for x in command if x.startswith('--export='))
                self.assertIn(f'TM_WORKFLOW_PYTHON={sys.executable}', export)
                self.assertIn('TM_SCREWER_CHILD=1', export)
            self.assertTrue((child / 'chain_submission.tsv').is_file())

    def test_missing_pexpect_fails_before_submitting(self):
        with tempfile.TemporaryDirectory() as t:
            child = self.prepared_child(Path(t))
            calls = []

            def fake_execute(command, *, cwd, **kwargs):
                calls.append(command)
                return SimpleNamespace(returncode=1, stdout='ModuleNotFoundError: pexpect')

            with mock.patch.object(screwer, 'execute', side_effect=fake_execute), \
                 mock.patch.object(screwer.shutil, 'which', return_value='/usr/bin/sbatch'), \
                 mock.patch.dict(screwer.os.environ, {'TM_WORKFLOW_PYTHON': sys.executable}):
                with self.assertRaisesRegex(RuntimeError, 'cannot import pexpect'):
                    screwer.submit(child, 'aoforce')
            self.assertEqual(1, len(calls))
            self.assertFalse((child / 'chain_submission.tsv').exists())

    def test_legacy_backend_does_not_require_or_pin_pexpect(self):
        with tempfile.TemporaryDirectory() as t:
            child = self.prepared_child(Path(t), backend='legacy')
            commands = []

            def fake_execute(command, *, cwd, **kwargs):
                commands.append(command)
                return SimpleNamespace(returncode=0, stdout=str(7800000 + len(commands)))

            with mock.patch.object(screwer, 'execute', side_effect=fake_execute), \
                 mock.patch.object(screwer.shutil, 'which', return_value='/usr/bin/sbatch'), \
                 mock.patch.dict(screwer.os.environ, {}, clear=True):
                self.assertEqual(2, screwer.submit(child, 'aoforce'))
            self.assertEqual(2, len(commands))
            for command in commands:
                self.assertEqual('sbatch', command[0])
                self.assertNotIn('TM_WORKFLOW_PYTHON=', ' '.join(command))

    def test_batch_script_uses_pinned_python_for_all_helpers(self):
        script = (Path(screwer.__file__).parent / 'turbo.pipeline.slurm').read_text()
        self.assertIn('DEFINE_PYTHON_BIN="${TM_WORKFLOW_PYTHON:-python3}"', script)
        self.assertIn('"$DEFINE_PYTHON_BIN" screwer_workflow.py --define-fallback', script)
        self.assertIn('python3 tm_workflow.py', script)
        self.assertNotIn('python3 screwer_workflow.py --define-fallback', script)


class AdditionalCyclesTests(unittest.TestCase):
    def test_automatic_cap_and_no_negative_modes_are_noops(self):
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            child = screwer.cycle_directory(base, 'lu', 1)
            child.mkdir(parents=True)
            info = dict(metal='lu', method='aoforce', cycle=1, max_cycles=1, project_root=str(base))
            (child / 'screwer_prepared.json').write_text(json.dumps(info))
            self.assertIsNone(screwer.advance_cycle(child))
            info['max_cycles'] = 3
            (child / 'screwer_prepared.json').write_text(json.dumps(info))
            (child / 'vibspectrum').write_text('$vibrational spectrum\n 1 0.00 0 0\n$end\n')
            with mock.patch.object(screwer, 'verify_parent'), mock.patch.object(screwer, 'prepare') as prepare:
                self.assertIsNone(screwer.advance_cycle(child))
                prepare.assert_not_called()

    def test_negative_spectrum_uses_new_cycle_and_pinned_cap(self):
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            child = screwer.cycle_directory(base, 'lu', 1)
            child.mkdir(parents=True)
            info = dict(metal='lu', method='aoforce', cycle=1, max_cycles=3,
                        project_root=str(base), define_backend='legacy', min_negative_cm1=1.0)
            (child / 'screwer_prepared.json').write_text(json.dumps(info))
            (child / 'vibspectrum').write_text('$vibrational spectrum\n 1 -12.50 0 0\n$end\n')
            new_path = screwer.cycle_directory(base, 'lu', 2)
            with mock.patch.object(screwer, 'verify_parent'), mock.patch.object(screwer, 'prepare', return_value=new_path) as prepare, mock.patch.object(screwer, 'submit') as submit:
                self.assertEqual(new_path, screwer.advance_cycle(child))
                self.assertEqual(2, prepare.call_args.kwargs['cycle'])
                self.assertEqual(3, prepare.call_args.kwargs['max_cycles'])
                submit.assert_called_once_with(new_path, 'aoforce')
            new_path.mkdir(parents=True)
            with mock.patch.object(screwer, 'verify_parent') as gate:
                with self.assertRaisesRegex(RuntimeError, 'already exists'):
                    screwer.advance_cycle(child)
                gate.assert_not_called()

    def test_manually_requested_cycle_2_reuses_parent_not_original(self):
        with tempfile.TemporaryDirectory() as t:
            base = Path(t)
            parent = screwer.cycle_directory(base, 'lu', 1)
            parent.mkdir(parents=True)
            new = parent / 'screwer_lu_cycle2'
            with mock.patch.object(screwer, 'verify_parent', return_value='automation'), \
                 mock.patch.object(screwer, 'spectrum_modes', return_value=[(1, -12.0)]), \
                 mock.patch.object(screwer, 'execute') as execute, \
                 mock.patch.object(screwer.shutil, 'which', return_value=None):
                # The parent was found and its spectrum inspected, before external cpc verification.
                (parent / 'vibspectrum').write_text('$vibrational spectrum\n 1 -12.0 0 0\n$end')
                with self.assertRaisesRegex(RuntimeError, 'required'):
                    screwer.prepare(base, 'lu', 'aoforce', cycle=2, max_cycles=2)
            self.assertFalse(new.exists())

if __name__ == '__main__':
    unittest.main()
