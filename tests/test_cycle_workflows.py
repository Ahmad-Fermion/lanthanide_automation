"""Development regression checks; mocked scheduling, no real TURBOMOLE."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import recover_universal as recovery
import screwer_workflow as screwer
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
