"""Regression coverage for automatic La/Lu dual preparation and submission."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def masters(base, *, restricted=True):
    (base / 'coord').write_text('$coord\n 0 0 0 gd f\n 0 0 1 h\n$end\n')
    (base / 'control').write_text('$title\n$uhf\n$fermi tmstrt=1000.00 nue=7.00\n$dft\n functional tpssh\n$end\n')
    (base / 'def.inp').write_text('a coord\n*\nno\nb all def2-SVP\nb "gd" def2-TZVP\n*\neht\n0\nno\nu 7\n*\nri\nq\n')
    if restricted:
        (base / 'def_closedshell.inp').write_text('a coord\n*\nno\nb all def2-SVP\nb "gd" def2-TZVP\n*\neht\n0\n\nri\n\nq\n')


def cmd(*args, cwd):
    p = subprocess.run(args, cwd=cwd, text=True, capture_output=True, env=os.environ.copy())
    if p.returncode:
        raise AssertionError(f'failed {args}: {p.stdout}\nSTDERR: {p.stderr}')
    return p.stdout + p.stderr


class DualLaLuDefaults(unittest.TestCase):
    def test_bash3_nounset_portability_for_empty_optional_arrays(self):
        """Bash 3.2 treats zero-element arrays as unset with `set -u`."""
        shell = (ROOT / 'submit-turbo-chain.sh').read_text()
        self.assertNotIn('${1,,}', shell)  # Bash 4+ only
        for name in ('args', 'created', 'submit_metals', 'closed_method'):
            # Every direct array expansion must occur inside the guarded form.
            guarded = '${' + name + '[@]+"${' + name + '[@]}"}'
            bare = '${' + name + '[@]}'
            self.assertGreater(shell.count(guarded), 0, name)
            self.assertEqual(shell.count(bare), shell.count(guarded), name)
        self.assertNotIn('${#created[@]}', shell)
        self.assertNotIn('${#selected_metals[@]}', shell)

    def test_setup_creates_both_formulations_and_preserves_frozen_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            msg = cmd(sys.executable, str(ROOT / 'setup_lanthanides.py'), '--base', str(base), '--metals', 'la', 'lu', cwd=base)
            for metal in ('la', 'lu'):
                self.assertIn(f'CREATED {metal}:', msg)
                self.assertIn(f'CREATED {metal}_closedshell:', msg)
                u = base / metal
                r = base / f'{metal}_closedshell'
                self.assertTrue(u.is_dir() and r.is_dir())
                self.assertIn(f'{metal} f', (u / 'coord').read_text())
                self.assertIn(f'{metal} f', (r / 'coord').read_text())
                self.assertIn(f'u {0}', (u / 'def.inp').read_text())
                self.assertNotIn('u 0', (r / 'def.inp').read_text())
                self.assertNotIn('$uhf', (r / 'control').read_text())
                self.assertTrue(json.loads((r / 'metal.json').read_text())['restricted'])
            # Repeated setup must not touch existing scientific folders.
            (base / 'la' / 'manual_notes.txt').write_text('untouched')
            cmd(sys.executable, str(ROOT / 'setup_lanthanides.py'), '--base', str(base), '--metals', 'la', 'lu', cwd=base)
            self.assertEqual('untouched', (base / 'la' / 'manual_notes.txt').read_text())

    def test_no_restricted_template_keeps_unrestricted_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base, restricted=False)
            msg = cmd(sys.executable, str(ROOT / 'setup_lanthanides.py'), '--base', str(base), '--metals', 'la', 'lu', cwd=base)
            self.assertTrue((base / 'la').exists())
            self.assertTrue((base / 'lu').exists())
            self.assertFalse((base / 'la_closedshell').exists())
            self.assertIn('no def_closedshell.inp', msg)

    def test_lan_run_prepare_only_prepares_both_in_new_and_existing_projects(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            cmd('bash', str(ROOT / 'submit-turbo-chain.sh'), '--base', str(base), '--metals', 'la', 'lu', '--prepare-only', cwd=base)
            self.assertTrue((base / 'la_closedshell').exists())
            self.assertTrue((base / 'lu_closedshell').exists())
            # If a restricted folder is missing, it can be recreated without
            # modifying the existing unrestricted calculation.
            import shutil
            shutil.rmtree(base / 'la_closedshell')
            (base / 'la' / 'test.marker').write_text('remain')
            cmd('bash', str(ROOT / 'submit-turbo-chain.sh'), '--base', str(base), '--metals', 'la', '--prepare-only', cwd=base)
            self.assertTrue((base / 'la_closedshell').is_dir())
            self.assertEqual('remain', (base / 'la' / 'test.marker').read_text())

    def test_normal_run_submits_both_formulations(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            fakebin = base / 'fakebin'
            fakebin.mkdir()
            sbatch = fakebin / 'sbatch'
            sbatch.write_text('#!/bin/sh\nprintf "%s|%s\\n" "$PWD" "$*" >> "$FAKE_SBATCH_LOG"\necho 76543\n')
            sbatch.chmod(0o755)
            squeue = fakebin / 'squeue'
            squeue.write_text('#!/bin/sh\nexit 0\n')
            squeue.chmod(0o755)
            env = os.environ.copy()
            env['PATH'] = str(fakebin) + os.pathsep + env['PATH']
            env['FAKE_SBATCH_LOG'] = str(base / 'sbatch.log')
            run = subprocess.run(['bash', str(ROOT / 'submit-turbo-chain.sh'),
                                  '--base', str(base), '--metals', 'la', 'lu'],
                                 cwd=base, env=env, text=True, capture_output=True)
            self.assertEqual(0, run.returncode, run.stdout + run.stderr)
            self.assertTrue((base / 'la_closedshell' / 'restricted_prepared.json').exists())
            self.assertTrue((base / 'lu_closedshell' / 'restricted_prepared.json').exists())
            calls = (base / 'sbatch.log').read_text().splitlines()
            self.assertEqual(10, len(calls), calls)  # 3+2 stages per metal
            for metal in ('la', 'lu'):
                self.assertEqual(3, sum(f'/{metal}|' in item for item in calls))
                self.assertEqual(2, sum(f'/{metal}_closedshell|' in item for item in calls))
                restricted = base / f'{metal}_closedshell'
                self.assertEqual(['restricted_ridft', 'jobex_eiger'],
                                 (restricted / 'chain_plan.txt').read_text().splitlines())
            # A repeat invocation must never double-submit existing chains.
            repeat = subprocess.run(['bash', str(ROOT / 'submit-turbo-chain.sh'),
                                     '--base', str(base), '--metals', 'la', 'lu'],
                                    cwd=base, env=env, text=True, capture_output=True)
            self.assertEqual(0, repeat.returncode, repeat.stdout + repeat.stderr)
            self.assertEqual(10, len((base / 'sbatch.log').read_text().splitlines()))

    def test_hessian_selection_applies_to_both_formulations(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            fakebin = base / 'fakebin'
            fakebin.mkdir()
            for name, body in {
                'sbatch': '#!/bin/sh\necho 12345\n',
                'squeue': '#!/bin/sh\nexit 0\n',
            }.items():
                f = fakebin / name
                f.write_text(body)
                f.chmod(0o755)
            env = os.environ.copy()
            env['PATH'] = str(fakebin) + os.pathsep + env['PATH']
            run = subprocess.run(['bash', str(ROOT / 'submit-turbo-chain.sh'),
                                  'ao', '--base', str(base), '--metals', 'la', 'lu'],
                                 cwd=base, env=env, text=True, capture_output=True)
            self.assertEqual(0, run.returncode, run.stdout + run.stderr)
            for metal in ('la', 'lu'):
                self.assertEqual('aoforce', (base / metal / 'chain_plan.txt').read_text().splitlines()[-1])
                self.assertEqual('aoforce', (base / f'{metal}_closedshell' / 'chain_plan.txt').read_text().splitlines()[-1])

    def test_prepared_restricted_folder_is_submitted_in_normal_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            cmd('bash', str(ROOT / 'submit-turbo-chain.sh'), '--base', str(base),
                '--metals', 'la', '--prepare-only', cwd=base)
            fakebin = base / 'fakebin'
            fakebin.mkdir()
            f = fakebin / 'sbatch'
            f.write_text('#!/bin/sh\necho 12345\n')
            f.chmod(0o755)
            f = fakebin / 'squeue'
            f.write_text('#!/bin/sh\nexit 0\n')
            f.chmod(0o755)
            env = os.environ.copy()
            env['PATH'] = str(fakebin) + os.pathsep + env['PATH']
            run = subprocess.run(['bash', str(ROOT / 'submit-turbo-chain.sh'),
                                  '--base', str(base), '--metals', 'la', '--include-existing'],
                                 cwd=base, env=env, text=True, capture_output=True)
            self.assertEqual(0, run.returncode, run.stdout + run.stderr)
            self.assertTrue((base / 'la' / 'chain_submission.tsv').is_file())
            self.assertTrue((base / 'la_closedshell' / 'chain_submission.tsv').is_file())

    def test_missing_restricted_template_keeps_unrestricted_submission(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base, restricted=False)
            fakebin = base / 'fakebin'
            fakebin.mkdir()
            f = fakebin / 'sbatch'
            f.write_text('#!/bin/sh\necho 12345\n')
            f.chmod(0o755)
            env = os.environ.copy()
            env['PATH'] = str(fakebin) + os.pathsep + env['PATH']
            run = subprocess.run(['bash', str(ROOT / 'submit-turbo-chain.sh'),
                                  '--base', str(base), '--metals', 'la'],
                                 cwd=base, env=env, text=True, capture_output=True)
            self.assertEqual(0, run.returncode, run.stdout + run.stderr)
            self.assertFalse((base / 'la_closedshell').exists())
            self.assertIn('SKIP la_closedshell', run.stdout)
            self.assertTrue((base / 'la' / 'chain_submission.tsv').exists())

    def test_explicit_restricted_only_mode_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            masters(base)
            cmd(sys.executable, str(ROOT / 'setup_lanthanides.py'), '--base', str(base), '--metals', 'la', '--closed-shell', cwd=base)
            self.assertTrue((base / 'la_closedshell').is_dir())
            self.assertFalse((base / 'la').exists())


if __name__ == '__main__':
    unittest.main()
