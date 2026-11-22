"""Regression tests for documented command-line entry points (no HPC required)."""
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def help_text(*command):
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                            timeout=15)
    if result.returncode:
        raise AssertionError(f'help command failed: {command}\n{result.stderr}')
    return result.stdout


class HelpMenuTests(unittest.TestCase):
    def test_recovery_separates_universal_and_legacy_modes(self):
        output = help_text('bash', str(ROOT / 'submit-turbo-recovery.sh'), '--help')
        for text in ('Universal recovery (preferred)', '--plan', '--submit',
                     '--frequency', '--json', 'Legacy SCF-only recovery',
                     '--continue-only', '--full-chain', '--ridft-only'):
            self.assertIn(text, output)
        self.assertIn('NOT universal', output)

    def test_run_describes_dual_formulations_and_prepare_only(self):
        output = help_text('bash', str(ROOT / 'submit-turbo-chain.sh'), '--help')
        for text in ('def_closedshell.inp', '--closed-shell', '--prepare-only',
                     '--include-existing', '--resume-partial', '--screwer'):
            self.assertIn(text, output)

    def test_setup_describes_automatic_la_lu(self):
        output = help_text(sys.executable, str(ROOT / 'setup_lanthanides.py'), '--help')
        self.assertIn('def_closedshell.inp', output)
        self.assertIn('lan-setup --closed-shell --metals la', output)

    def test_screwer_documents_cycle_and_overwrite(self):
        output = help_text(sys.executable, str(ROOT / 'screwer_workflow.py'), '--help')
        for text in ('--cycle', '--max-cycles', '--overwrite-prepared',
                     '--define-backend', '--submit-prepared', '--min-imag'):
            self.assertIn(text, output)

    def test_restricted_help_still_available(self):
        output = help_text('bash', str(ROOT / 'submit-turbo-chain.sh'),
                           '--closed-shell', '--help')
        self.assertIn('--metals', output)
        self.assertIn('--include-existing', output)


if __name__ == '__main__':
    unittest.main()
