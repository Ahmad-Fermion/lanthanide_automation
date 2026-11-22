"""Initial DEFINE fallback retains metal-specific input and restores inputs."""
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tm_workflow as tm


PRIMARY = '''a coord
ired
*
b all def2-SVP
b "gd" def2-TZVP
*
use ../control

y
ri
on
q
'''


class InitialDefineFallbackTests(unittest.TestCase):
    def sample_dir(self, root, transcript=PRIMARY):
        (root / 'def.inp').write_text(transcript)
        control = ('$dft\n functional tpssh\n$fermi tmstrt=1000.0 nue=3.00\n'
                   '$scfdamp start=6.0 step=0.05 min=0.10\n$end\n')
        (root / 'control').write_text(control)
        (root / 'coord').write_text('$coord\n 0 0 0 gd\n$end\n')
        return control

    def test_switches_complete_menu_sequence_and_preserves_gd_tail(self):
        derived = tm.initial_noired_transcript(PRIMARY)
        self.assertIn('a coord\n*\nno\n', derived)
        self.assertNotIn('\nired\n', derived)
        self.assertEqual(PRIMARY.split('b all def2-SVP')[1], derived.split('b all def2-SVP')[1])
        self.assertIsNone(tm.initial_noired_transcript('a coord\n*\nno\nq\n'))

    def test_failed_ired_retry_restores_control_and_keeps_original_transcript(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            control = self.sample_dir(root)
            calls = []
            def simulated_define(command, **kwargs):
                transcript = kwargs['input']
                calls.append(transcript)
                self.assertEqual(control, (root / 'control').read_text())
                if '\nired\n' in transcript:
                    (root / 'control').write_text('$corrupted\n$end\n')
                    (root / 'basis').write_text('bad basis')
                    return subprocess.CompletedProcess(command, 1, 'SINCE THERE ARE THREE COLLINEAR BOND PARTNERS\n')
                (root / 'basis').write_text('$basis\n')
                (root / 'auxbasis').write_text('$jbas\n')
                return subprocess.CompletedProcess(command, 0, 'define ended normally\n')
            with mock.patch.object(tm, 'HERE', root), mock.patch.object(tm.subprocess, 'run', side_effect=simulated_define), contextlib.redirect_stdout(io.StringIO()) as output:
                tm.run_initial_define()
            self.assertEqual(2, len(calls))
            self.assertEqual(PRIMARY, (root / 'def.inp').read_text())
            self.assertEqual(control, (root / 'control').read_text())
            self.assertIn('no\n', (root / 'def.noired.generated.inp').read_text())
            self.assertIn('fallback succeeded', output.getvalue())
            self.assertIn('COLLINEAR', (root / 'define.out').read_text())
            self.assertIn('ended normally', (root / 'define.noired.out').read_text())

    def test_unmodified_no_ired_succeeds_only_once(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            original = self.sample_dir(root, transcript='a coord\n*\nno\nq\n')
            def fake(command, **kwargs):
                (root / 'basis').write_text('basis')
                (root / 'auxbasis').write_text('auxbasis')
                return subprocess.CompletedProcess(command, 0, 'define ended normally\n')
            with mock.patch.object(tm, 'HERE', root), mock.patch.object(tm.subprocess, 'run', side_effect=fake) as call:
                tm.run_initial_define()
            call.assert_called_once()
            self.assertEqual(original, (root / 'control').read_text())
            self.assertFalse((root / 'def.noired.generated.inp').exists())

    def test_both_failed_restores_original_and_never_hides_failure(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            control = self.sample_dir(root)
            def fail(command, **kwargs):
                (root / 'control').write_text('corrupt')
                (root / 'basis').write_text('partial output')
                return subprocess.CompletedProcess(command, 1, 'define ended abnormally\n')
            with mock.patch.object(tm, 'HERE', root), mock.patch.object(tm.subprocess, 'run', side_effect=fail), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, 'both failed'):
                    tm.run_initial_define()
            self.assertEqual(control, (root / 'control').read_text())
            self.assertFalse((root / 'basis').exists())
            self.assertTrue((root / 'define.noired.out').exists())

    def test_original_script_not_matching_known_menu_is_not_rewritten(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            self.sample_dir(root, transcript='a coord\nired\ninvalid-answer\nq\n')
            def fail(command, **kwargs):
                return subprocess.CompletedProcess(command, 1, 'error\n')
            with mock.patch.object(tm, 'HERE', root), mock.patch.object(tm.subprocess, 'run', side_effect=fail) as call:
                with self.assertRaisesRegex(RuntimeError, 'no recognized'):
                    tm.run_initial_define()
            call.assert_called_once()

    def test_pipeline_calls_helper_and_check(self):
        batch = (Path(tm.__file__).parent / 'turbo.pipeline.slurm').read_text()
        self.assertIn('python3 tm_workflow.py run-initial-define', batch)
        self.assertIn('python3 tm_workflow.py check-define', batch)


if __name__ == '__main__':
    unittest.main()

class PreparedOverwriteTests(unittest.TestCase):
    def test_force_replaces_only_child_when_no_jobs_are_active(self):
        import screwer_workflow as sw
        with tempfile.TemporaryDirectory() as t:
            parent = Path(t) / 'nd' / 'screwer_nd'
            child = parent / 'screwer_nd_cycle2'
            child.mkdir(parents=True)
            (parent / 'control').write_text('parent unchanged')
            (child / 'broken_input').write_text('stale')
            with mock.patch.object(sw.shutil, 'which', return_value='/usr/bin/squeue'), \
                 mock.patch.object(sw, 'execute', return_value=subprocess.CompletedProcess([], 0, '', '')):
                sw.replace_prepared_child(child, parent, force=True)
            self.assertFalse(child.exists())
            self.assertEqual('parent unchanged', (parent / 'control').read_text())

    def test_force_rejects_active_job(self):
        import screwer_workflow as sw
        with tempfile.TemporaryDirectory() as t:
            parent = Path(t) / 'lu'
            child = parent / 'screwer_lu'
            child.mkdir(parents=True)
            with mock.patch.object(sw.shutil, 'which', return_value='/usr/bin/squeue'), \
                 mock.patch.object(sw, 'execute', return_value=subprocess.CompletedProcess([], 0, '6606465|' + str(child.resolve()) + '\n', '')):
                with self.assertRaisesRegex(RuntimeError, 'active/pending'):
                    sw.replace_prepared_child(child, parent, force=True)
            self.assertTrue(child.exists())

    def test_noninteractive_requires_opt_in_and_preserves_dir(self):
        import screwer_workflow as sw
        with tempfile.TemporaryDirectory() as t:
            parent = Path(t) / 'nd'
            child = parent / 'screwer_nd'
            child.mkdir(parents=True)
            with mock.patch.object(sw.sys.stdin, 'isatty', return_value=False):
                with self.assertRaisesRegex(RuntimeError, '--overwrite-prepared'):
                    sw.replace_prepared_child(child, parent)
            self.assertTrue(child.exists())

    def test_rejects_symlink_child(self):
        import screwer_workflow as sw
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            parent = root / 'nd'
            parent.mkdir()
            (root / 'important').mkdir()
            child = parent / 'screwer_nd'
            child.symlink_to(root / 'important', target_is_directory=True)
            with self.assertRaisesRegex(RuntimeError, 'not a direct, real directory'):
                sw.replace_prepared_child(child, parent, force=True)
            self.assertTrue((root / 'important').exists())
