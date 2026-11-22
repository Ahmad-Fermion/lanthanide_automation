#!/usr/bin/env python3
"""Opt-in supervised DEFINE backend for inherited Screwer continuations.

Uses pexpect to catch known geometry failures immediately rather than waiting
for a fixed answer script to reach EOF. The existing, user-authored DEFINE
scripts remain authoritative: this is NOT yet a replacement for fresh ECP or
X2CSR initialization. Caller owns snapshot/restore and final input validation.
"""
import io
import os
from pathlib import Path
import re
import subprocess

IREDS = (
    r'SINCE THERE ARE THREE COLLINEAR BOND PARTNERS',
    r'ERROR IN CARTESIAN COORDINATES\s*\?',
    r'THE B\*m\*Bt-MATRIX IS SINGULAR',
    r'CAUTION, within .*? steps NO convergence',
)


def define_supervised(cwd, script):
    """Run a pinned DEFINE transcript; return CompletedProcess-like diagnostic.

    The expect loop treats unexpected early EOF, stuck prompts and known ired
    errors as failures. Nothing is retried or repaired here: screwer_workflow
    restores the prepared snapshot and selects the no-ired template.
    """
    try:
        import pexpect
    except ImportError as exc:
        raise RuntimeError('adaptive DEFINE requires pexpect in the SLURM Python environment') from exc
    timeout = int(os.environ.get('LAN_DEFINE_TIMEOUT', '300'))
    if timeout < 10 or timeout > 3600:
        raise RuntimeError('LAN_DEFINE_TIMEOUT must be between 10 and 3600 seconds')
    log = io.StringIO()
    child = pexpect.spawn('define', cwd=str(Path(cwd).resolve()),
                          encoding='utf-8', timeout=timeout, echo=False)
    child.logfile_read = log
    code = 1
    try:
        child.send(script)
        patterns = [*IREDS, r'define\s+ended\s+normally',
                    r'define\s+ended\s+abnormally', pexpect.EOF, pexpect.TIMEOUT]
        idx = child.expect(patterns, timeout=timeout)
        if idx == len(IREDS):
            # Normal termination, not just an exit code of 0.
            follow = child.expect([pexpect.EOF, pexpect.TIMEOUT], timeout=30)
            if follow == 0:
                child.close()
                code = child.exitstatus if child.exitstatus is not None else 1
        elif idx in (len(IREDS) + 1, len(IREDS) + 2, len(IREDS) + 3):
            code = 1
        else:
            # Known ired error; do not send guessed answers to an unknown menu.
            code = 1
    finally:
        if child.isalive():
            child.terminate(force=True)
        child.close(force=True)
    return subprocess.CompletedProcess(['define'], code, log.getvalue())
