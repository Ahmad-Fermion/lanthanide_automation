"""Verified legacy SCF can continue without a fabricated DEFINE snapshot."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import recover_universal as recovery
import tm_workflow as tm

SCF_OK = 'ENERGY CONVERGED !\nconvergence criteria satisfied after 12 iterations\nridft : all done\n'
CTRL = ("$symmetry c1\n$fermi tmstrt= 1000.00 nue= 6.00\n"
        "$alpha shells\n a 1-157 ( 1 )\n$beta shells\n a 1-151 ( 1 )\n$end\n")


def sample(folder, *, resume=True, fermi=True):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'tm_workflow.py').write_text('# helper')
    (folder / 'coord').write_text('$coord\n0 0 0 eu\n$end\n')
    (folder / 'metal.json').write_text(json.dumps({'metal':'eu','charge':-2,'unpaired':6}))
    (folder / 'control').write_text(CTRL if fermi else CTRL.replace('$fermi tmstrt= 1000.00 nue= 6.00\n',''))
    (folder / 'alpha').write_text('a orbital data')
    (folder / 'beta').write_text('b orbital data')
    (folder / 'chain_plan.txt').write_text('ridft_eiger_fix\ninteger_ridft\njobex_eiger\nnumforce\n')
    if resume:
        (folder / 'ridft.resume.ok').write_text('job_id=123\nlog=ridft.resume.123.out\n')
        (folder / 'ridft.resume.123.out').write_text(SCF_OK)


class LegacyOccupationReference(unittest.TestCase):
    def test_eu_plan_is_read_only_and_starts_with_occupation_fix(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/'eu'
            sample(folder)
            stages, reason = recovery.stage_plan(folder,'standard','numforce')
            self.assertEqual(['eiger_fix_integer_ridft','jobex_eiger','numforce'],stages)
            self.assertIn('legacy',reason)
            self.assertFalse((folder/'recovery_occupations').exists())
            self.assertEqual(CTRL,(folder/'control').read_text())

    def test_capture_recovery_shells_then_restore_after_integer_scf(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/'eu'
            sample(folder)
            with (mock.patch.object(tm,'HERE',folder),
                  mock.patch.object(tm,'META',folder/'metal.json'),
                  mock.patch.object(tm,'STATE',folder/'workflow_status.json'),
                  mock.patch.object(tm,'OCC_REFERENCE',folder/'initial_occupations'),
                  mock.patch.object(tm,'RECOVERY_OCC_REFERENCE',folder/'recovery_occupations')):
                with contextlib.redirect_stdout(io.StringIO()) as stdout:
                    tm.remove_fermi()
                self.assertIn('NOT initial_occupations',stdout.getvalue())
                self.assertFalse((folder/'initial_occupations').exists())
                self.assertIn('$alpha shells',(folder/'recovery_occupations').read_text())
                self.assertNotIn('$fermi',(folder/'control').read_text())
                self.assertIn('$fermi',(folder/'control.with_fermi').read_text())
                (folder/'control').write_text((folder/'control').read_text().replace('1-157','1-156'))
                with contextlib.redirect_stdout(io.StringIO()):
                    tm.restore_after_integer_scf()
                self.assertIn('1-157',(folder/'control').read_text())
                self.assertFalse((folder/'initial_occupations').exists())

    def test_reject_unsupported_or_unverified_state_without_mutating_control(self):
        for variant in ('no_resume','wrong_spin','wrong_nue','noninteger','missing_orbital'):
            with self.subTest(variant=variant), tempfile.TemporaryDirectory() as t:
                folder=Path(t)/'eu'
                sample(folder,resume=(variant!='no_resume'))
                control=(folder/'control')
                if variant=='wrong_spin': control.write_text(control.read_text().replace('1-151','1-152'))
                if variant=='wrong_nue': control.write_text(control.read_text().replace('nue= 6.00','nue= 5.00'))
                if variant=='noninteger': control.write_text(control.read_text().replace('( 1 )','( 0.5 )',1))
                if variant=='missing_orbital': (folder/'alpha').unlink()
                before=control.read_bytes()
                with self.assertRaises(recovery.RecoveryError):
                    recovery.stage_plan(folder,'standard','numforce')
                self.assertEqual(before,control.read_bytes())
                self.assertFalse((folder/'recovery_occupations').exists())

    def test_normal_initial_reference_remains_preferred(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/'eu'
            sample(folder,resume=False)
            (folder/'initial_occupations').write_text('$alpha shells\n a 1-157 ( 1 )\n$beta shells\n a 1-151 ( 1 )\n')
            stages,_=recovery.stage_plan(folder,'standard','numforce') if (folder/'ridft.out').exists() else ([], '')
            with (mock.patch.object(tm,'HERE',folder),
                  mock.patch.object(tm,'META',folder/'metal.json'),
                  mock.patch.object(tm,'STATE',folder/'workflow_status.json'),
                  mock.patch.object(tm,'OCC_REFERENCE',folder/'initial_occupations'),
                  mock.patch.object(tm,'RECOVERY_OCC_REFERENCE',folder/'recovery_occupations')):
                with contextlib.redirect_stdout(io.StringIO()): tm.remove_fermi()
            self.assertFalse((folder/'recovery_occupations').exists())
