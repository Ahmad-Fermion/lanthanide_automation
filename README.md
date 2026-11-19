# Lanthanide automation · TURBOMOLE 7.9 on VACC

Prepare, submit, monitor, and recover **14 Ln3+ calculations (La–Lu, excluding Pm)** using a shared command-line toolset. Developed for **UVM's Vermont Advanced Computing Center (VACC)** with TURBOMOLE 7.9 and SLURM, the multi-stage workflow can be adapted to other HPC clusters by adjusting their SLURM resources, software modules, and scratch configuration.

A representative calculation is first configured interactively using `tee def.inp | define`, which records the setup instructions in `def.inp`. A  master `control` file needs to prepared. The automation then reuses these templates, followed by SCF calculations (`ridft`), geometry optimization (`jobex`), and optional Hessian calculations (`aoforce` or `NumForce`). Although developed for lanthanides, the approach is readily generalizable to other elements and systematic DFT benchmarking studies involving different basis sets, functionals, and molecular systems.

Run the commands from a project directory containing `control`, `coord`, and `def.inp` (and `def_gd.inp` when preparing Gd). Each project's input files define the molecular system and theory.

## Install once on VACC (from GitLab)

Copy the **SSH clone URL** from this repository's GitLab page, then on VACC:

```bash
mkdir -p ~/src ~/bin
git clone <SSH-CLONE-URL> ~/src/lanthanide_automation
bash ~/src/lanthanide_automation/install.sh
export PATH="$HOME/bin:$PATH"   # add to ~/.bashrc only if ~/bin is not already in PATH
```

The installer adds `lan-run`, `lan-recover`, `lan-setup`, and `lan-status` to `~/bin`, including the submission safety checker used by `lan-run`. To update an existing installation:

```bash
cd ~/src/lanthanide_automation
git pull --ff-only
bash install.sh
```

**Script filenames:** The central submission launcher is `submit-turbo-chain.sh`, and the calculation-stage batch script is `turbo.pipeline.slurm`. The recovery and Screwer scripts remain `submit-turbo-recovery.sh`, `turbo.resume.slurm`, and `turbo.screwer.slurm`.

**Existing metal directories:** Previously prepared metal directories contain their own copied batch scripts. Renaming the central scripts does **not** update those copies. Do not rename or replace an in-use batch script while jobs are pending or running; first inspect outstanding jobs, then manually migrate the per-metal script when safe. Existing metal directories are never changed automatically.

**Requirements:** VACC account, Python 3, TURBOMOLE 7.9, and SLURM. 

## Four commands to remember

| Command | Purpose |
| --- | --- |
| `lan-run` | Prepare missing metal folders and submit the normal SLURM chain |
| `lan-recover` | Restart failed RIDFT with your manual `control` edits; optionally continue |
| `lan-setup` | Prepare folders **only**, without submitting |
| `lan-status` | Show progress; write `workflow_status.csv` |
| `lan-validate` | Read-only input validation using optional `turbomoleio` |

By default, `lan-run` stops after geometry optimization and the post-optimization EIGER check, with **no Hessian** (useful for large models). Choose `lan-run num` = **NumForce** (for X2CSR) or `lan-run ao` = **aoforce** to include a Hessian. Full names `numforce`/`aoforce` also work.

## Common examples

Run from a single EF project directory, e.g.:

```bash
cd ~/lanM_shared/small-model/x2csr-d4-wd-scf7/EF1738
lan-run num --metals sm --prepare-only      # create only Sm for inspection
lan-run num --metals sm --include-existing  # submit that prepared Sm
lan-run                                     # submit all other *missing* metals, no Hessian
lan-run ao                                  # same workflow, ending with aoforce
lan-run ao --metals gd tb                   # only new Gd/Tb (Gd uses def_gd.inp)
lan-status                                  # progress + workflow_status.csv
```

`--include-existing` is only for a **pristine, never-run** directory created by `lan-setup`. Before submission, the tool now checks that the metal-specific `control`, `coord`, `def.inp`, and `metal.json` match freshly prepared inputs from the current project masters, and that no calculation outputs, queued jobs, or other unexpected files are present. If a project master was changed after setup, the directory is conservatively refused for manual review. It will **not** run `define` in an edited or previously calculated directory. Without `--include-existing`, existing folders are skipped. No previously existing scientific files are overwritten during these checks.

### Interrupted SLURM submissions

If a `sbatch` submission fails for one metal, `lan-run` records the partial chain, reports an error, **continues submitting other metals**, and exits nonzero at the end if there were failures. It does not restart already submitted jobs. To resume only the **missing submissions**, explicitly use the same `num`/`ao` selector as the original command:

```bash
lan-run num --metals gd --resume-partial  # only unsent stages of an interrupted NumForce chain
lan-run --metals sm --resume-partial      # only unsent stages of a no-Hessian chain
```

The new `chain_plan.txt` records the original full 3- or 4-job stage plan; `chain_submission.tsv` records each accepted SLURM job ID. Resumption checks that the recorded stages are a valid prefix of that plan and queries `squeue` and `sacct` for the previous jobs' states and working directory. It depends on the last job if still active, or submits the next stage without a stale dependency only after confirming the last job finished successfully. The calculations themselves are **not rerun**. Failed, canceled, missing, mismatched, or ambiguous job records stop resumption and require manual review or `lan-recover` as appropriate. `--resume-partial` requires `--metals` and cannot be mixed with `--include-existing`.

**Important:** Previously submitted chains from the older tool do not have `chain_plan.txt`, so they are **not automatically resumable** by this new command. The original Hessian choice cannot safely be inferred from an incomplete record. A `.lan_submission_pending` marker means the tool may have been interrupted after SLURM accepted a job but before its ID was saved; inspect SLURM carefully before reconciliation, and do not blindly delete this marker or `chain_submission.tsv`. A `.lan_submission_lock` prevents concurrent submissions to the same metal; stale locks also require manual verification.

## Milestone 1: read-only input validation (`lan-validate`)

This optional command uses the **installed `turbomoleio` 1.6.0 library** to parse
TURBOMOLE `control` data groups. It does not call DEFINE, RIDFT, JOBEX, SLURM,
or change any inputs/results. **No existing execution command calls it:** the
production `lan-run`/`lan-screwer`/`lan-recover` pipelines are unchanged.

Install the dependency in a separate Python virtual environment (do not modify
an environment being used by running jobs):

```bash
python3 -m venv ~/venvs/lan-input-validation
source ~/venvs/lan-input-validation/bin/activate
python -m pip install "turbomoleio==1.6.0"
# From this source repository: bash install.sh
```

The `lan-validate` launcher uses `python3` on the current `PATH`; activate the
above environment whenever using the command. For the installed command, run:

```bash
cd /gpfs1/pi/crreinha/lanM_shared/small-model/svp-d4-wd-scf7/EF1738
lan-validate --project-only                 # check masters for all 14 Ln
lan-validate --metals pr eu er tm --include-screwer
lan-validate --metals la gd                # check specific metal input folders
lan-validate --metals pr --include-screwer --json  # JSON to stdout; nothing written
```

The validator distinguishes **project masters**, **fresh prepared inputs before
DEFINE**, and **existing DEFINE-generated input files** (including Screwer
corrections). Project masters may refer to basis/MO files not generated yet.
Checks include readable data groups, duplicate groups, coordinate atom counts and
finite positions, metal identity, referenced-file presence after DEFINE, and
matching *stable* scientific settings from the project master or parent:
functional, DFT grid, weight derivatives, SCF convergence, D4, ECP presence,
selected relativistic flags, and related method switches. It checks metal.json
identity and charge against the local initial `def.inp` when both are available.
The initial charge is **never guessed** from `$charge from ridft` and the code
does not compare or restore Screwer orbital coefficients or occupations.

**Limitations:** Milestone 1 is an *input preflight*, not validation of the
physically correct ground state, exact basis/ECP assignment for every atom,
or SCF/geometry/frequency convergence. It does not yet parse the output logs.
The Sci settings checks are scoped to those listed above; unusual X2CSR flags,
auxiliary basis content, electronic occupations and all method-specific choices
still need independent scientific review. Because actively running jobs change
files, an input snapshot taken while they run can produce transient warnings.
The validation exit code is 1 for any FAIL, 0 for PASS/WARN, and 2 if the
`turbomoleio` dependency is unavailable. It is not connected to automatic
submission blocking. Offline unit tests run with
`PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v`.
A real TURBOMOLE 7.9/VACC integration test is required before relying on its
PASS result for production submission decisions.


### Adaptive Screwer DEFINE: Python on compute nodes

`lan-screwer --define-backend adaptive` requires `pexpect`. On VACC,
`module load turbomole/7.9` selects a separate TURBOMOLE Python interpreter,
which may not have `pexpect` installed. When submitting an *adaptive* Screwer
cycle, `lan-screwer` now checks the Python interpreter used by the submitter,
then passes its absolute path as `TM_WORKFLOW_PYTHON` to the SLURM jobs. The
batch script uses that interpreter for adaptive DEFINE after module loading;
other workflow helpers and normal/legacy submissions continue to use the
existing module Python.

Activate your development environment when submitting adaptive Screwer jobs:

```bash
source ~/venvs/lan-dev/bin/activate
~/bin/lan-dev/lan-screwer ao --metals pr --define-backend adaptive --prepare-only
~/bin/lan-dev/lan-screwer ao --metals pr --submit-prepared
```

For an alternative environment, set `TM_WORKFLOW_PYTHON` to an absolute path
to an executable Python with `pexpect` installed. Submission checks that the
selected interpreter can import `pexpect` **before** scheduling either job.
This does not modify `control`, SCF settings, or DEFINE templates.

## Monitoring calculations

```bash
lan-status               # one clear status per metal + workflow_status.csv
lan-status --details     # also show SCF, geometry, frequency, and occupation history
                           # separate screwer_<metal> rows appear when present
lan-status --base DIR    # check a different EF project directory
```

The summary separates **OVERALL**, **LAST STAGE**, **LIVE SLURM**, and **SLURM FINAL**. `OVERALL` is `COMPLETED` only when the current attempt has a recorded `DONE/SUCCESS` **and** all tracked SLURM jobs are verified as successfully completed. `RUNNING`/`PENDING` require live queue evidence; `FAILED` includes failed, cancelled, timed-out, and out-of-memory jobs (the exact cause remains visible under `SLURM FINAL` and in the CSV). `NEEDS_REVIEW` means the evidence is incomplete, inconsistent, or a submission is partial. `PREPARED` and `NOT_PREPARED` identify jobs that have not been submitted. `LAST STAGE` describes the current or last submitted computational stage; the combined JOBEX job also runs post-optimization EIGER.

`LIVE SLURM` comes from `squeue`; `INACTIVE` means no tracked jobs were found in the live queue, **not** that they succeeded. Historical states come from a batched `sacct` query. The script does not claim success when accounting data, the expected working directory, or a current final workflow record cannot be verified. `INCOMPLETE` in the recorded workflow result means there is no final outcome recorded for this attempt; it does **not** claim the calculation is presently running. Legacy chains without an original `chain_plan.txt` may show `NEEDS_REVIEW` because their intended final stage is uncertain.

The CSV retains detailed SCF, integer-SCF, geometry, frequency, and orbital-hole fields, plus `overall`, `last_stage`, `slurm_final`, and `status_note`. The `NEG FREQ` column and CSV `negative_freqs` field count strictly negative wave numbers in the first eight mode rows of a current, completed `vibspectrum` file (excluding `-0.00`); `-` means unavailable or unverified. Detailed per-stage history can include earlier attempts; inspect the per-metal `workflow_history.jsonl` when investigating a recovery. The status command is read-only with respect to metal folders: it only updates the project-level `workflow_status.csv`. If `squeue`/`sacct` is unavailable or a job is too old for accounting retention, results may require manual review.

## When RIDFT needs manual intervention

After the original failed/blocked jobs are stopped, edit e.g., `tb/control` manually. From the EF project directory:

```bash
lan-recover --metal tb                      # run only RIDFT; inspect state
lan-recover --metal tb --continue-only      # later, if recovery converged; no Hessian
```

Or, if confident in your changes:

```bash
lan-recover --metal tb --full-chain         # resume RIDFT + later stages; no Hessian
lan-recover ao --metal tb --full-chain      # same, ending with aoforce
lan-recover num --metal tb --continue-only  # continue after recovery, ending with NumForce
```

Like `lan-run`, `lan-recover` only requests a Hessian when `ao` or `num` is specified. `--ridft-only` remains the default recovery mode and never schedules later stages.

## What happens automatically

`define < def.inp` → `ridft` → `eiger` → occupation fix → second `ridft` → `jobex` → `eiger` → **stop** (default `lan-run`), or → `NumForce`/**`aoforce`** (`lan-run num`/`lan-run ao`). The **scientific steps and their order are unchanged**, but lightweight steps run inside existing compute allocations:

- **Job 1:** `define` + initial RIDFT + pre-optimization EIGER + occupation fix.
- **Job 2:** integer-occupation RIDFT.
- **Job 3:** JOBEX + post-optimization EIGER; finish here for `lan-run` without Hessian.
- **Job 4 (optional):** NumForce or aoforce.

Thus `lan-run` submits **3 dependent jobs** without a Hessian, or **4** with one. The three expensive SCF/optimization stages and optional Hessian keep separate **48-hour `general`** allocations with their existing CPU/memory settings. EIGER and occupation editing execute inside those jobs; no extra EIGER/`short` allocations are submitted. Jobs use `afterok` dependencies.

`lan-recover --full-chain` similarly uses **3 jobs** without Hessian, or **4** with it: recovered RIDFT (without `define`), combined EIGER + occupation fix + integer RIDFT, then JOBEX + EIGER, then optional Hessian. `--continue-only` skips the recovered RIDFT job (**2 or 3 jobs**); `--ridft-only` remains **1 job**. Pre-existing metal folders are not modified automatically: older copies of the per-metal SLURM script must be reviewed and updated manually before using these consolidated stages.

After `define`, the initial integer alpha/beta shell assignments are saved in `initial_occupations`. The occupation-fix stage restores those blocks and **removes `$fermi`** before the second RIDFT (existing behavior). Orbital holes are **warnings**; true SCF or geometry failure blocks downstream jobs. `control` is the only TURBOMOLE control file; no `control.template`. The validated NumForce OpenMP environment and `PARNODES=1 NumForce -ri -central` remain unchanged.

**Requirements:** Python 3, SLURM, TURBOMOLE 7.9 on VACC. For new `gd/`, provide `def_gd.inp` in the project root. Defaults: 20 CPUs/20 GB for RIDFT/JOBEX, 100 GB for frequency; overrides via `TM_CPUS`, `RIDFT_MEM`, `JOBEX_MEM`, `HESS_MEM`, `EIGER_CPUS`, `EIGER_MEM` (legacy standalone-stage settings; not used by consolidated submissions).

**Help:** `lan-run --help`, `lan-recover --help`, `lan-setup --help`, `lan-status --help`. The per-metal files (including `chain_submission.tsv`, `initial_occupations`, stage logs and warnings) remain in that EF project. Verify orbitals/occupations and frequencies scientifically.

## Optional imaginary-frequency correction (`lan-screwer`)

This is an **opt-in, one-pass** screwer workflow. It does not change existing
`lan-run` or `lan-recover` behavior unless selected. `lan-status` also
shows separate rows for any correction directories that exist.
Place your **own tested** `screwer.inp` (mode/displacement prompts) and
`def_screwer.inp` (DEFINE menu responses for displaced coordinates) in each
project root alongside the master `def.inp`. For **Gd**, additionally provide
`def_screwer_gd.inp` containing the Gd-specific continuation menu responses.
The Screwer DEFINE template is copied **verbatim** (not modified using the
initial setup's `b`/`u` basis and spin selectors), because the continuation
inherits orbitals and shell occupations from the completed metal calculation.
For Gd, `def_screwer_gd.inp` is selected and copied as `def_screwer.inp`.
For non-Gd metals, optionally place **`def_screwer_ln.inp`** next to the other
project templates: it is a complete alternate DEFINE input that skips `ired`.
If the primary DEFINE does not finish normally or leaves `control` missing
`$coord` or `$atoms`, the new Screwer correction job restores its prepared
inputs and retries once with `def_screwer_ln.inp`. The fallback is copied into
the correction directory during preparation. For Gd, the generic Ln template
is **never** substituted: a distinct optional `def_screwer_gd_ln.inp` is needed
for automatic fallback if its dedicated DEFINE input fails.
The Screwer template does **not** need `eht`, a charge, `b`, or `u` commands:
those settings are inherited from the completed calculation's `control` and
orbital files. When available, the charge recorded in `screwer_prepared.json`
comes from the original metal's `metal.json`, which was prepared using
`def.inp` or `def_gd.inp`; it is **not** inferred from `def_screwer.inp`.
For a manual parent without `metal.json`, the charge is recorded as unknown
rather than guessed; the existing `control` is inherited. After DEFINE,
Screwer DEFINE is allowed to change inherited orbital files, occupations,
and `control`; it does not perform the initial Fermi/occupation-reset stages.
Only DEFINE's successful exit and required resulting input files/data groups
are checked before RIDFT.

**Manually completed TURBOMOLE calculations:** `lan-screwer` also accepts a
metal folder without `workflow_status.json`, `chain_submission.tsv`, or
`metal.json` when the existing results independently demonstrate a converged
JOBEX and a completed Hessian. It checks `GEO_OPT_CONVERGED`, nonempty
`jobex.out`, `gradient`, `hessian`, `vibspectrum`, method-specific
`aoforce.out` or `numforce.out`, `control`, `coord`, `basis`, `auxbasis`,
existing `alpha`/`beta` or `mos`, and an inactive SLURM working directory.
Incomplete or contradictory automation records are **not** treated as manual.
The requested `ao`/`num` method must match the available frequency output.
If no parent `metal.json` exists, only the correction child receives minimal
metadata (`metal` plus unknown charge/spin); these values are **not** invented
or used to reinitialize the electronic state. The original parent remains
untouched. A pre-existing legacy `<metal>/screwer/` is left unchanged and
**does not prevent** creating a separate `<metal>/screwer_<metal>/`.
Use the project root for `--base`; if you are inside `gd/`, use `--base ..`.
For example, from inside a manually completed `gd/` folder:

```bash
lan-screwer ao --base .. --metals gd --prepare-only
```

For a completed, verified Hessian calculation:

```bash
lan-screwer num --metals gd --prepare-only  # inspect displaced inputs before submitting
lan-screwer num --metals gd --submit-prepared # submit one correction chain after inspection
# Or in one step when confident and no screwer_<metal>/ directory exists:
lan-screwer num --metals gd
lan-screwer ao --metals ce
```

For *future* calculations, opt in when you submit:

```bash
lan-run num --metals gd --screwer
lan-run ao --metals sm --screwer
```

`--screwer` requires `num` or `ao` (there is no spectrum for bare
`lan-run`). It adds **one** dependent `short` gate job (1 CPU, 1 GB, 10 minutes)
after the fourth Hessian job. Once the Hessian has recorded DONE/SUCCESS, the
gate checks the first eight mode rows of `vibspectrum`. If none is strictly
negative, it does nothing. If a negative mode exists, it runs `cpc screwer_<metal>`,
copies the screwer input and a validated metal-specific DEFINE template to
the **new** `<metal>/screwer_<metal>/` folder,
runs `screwer < screwer.inp`, checks that `$newcoord` describes a genuinely
negative mode, generates `$coord` using `sdg -b newcoord`, then invokes
`kdg newcoord`, `kdg coord`, `kdg symmetry`, `kdg optimize` in that order.
It prepares `screwer_prepared.json` and submits a **separate two-job**
correction chain: **(1) primary Screwer DEFINE, with one clean fallback
attempt if necessary → `actual -r` → RIDFT → JOBEX → EIGER;
(2) AOFORCE or NumForce**. The first DEFINE output is saved as `define.out`;
any retry output is saved separately as `define.fallback.out`. A successful
DEFINE must report normal completion and leave usable `$coord` and `$atoms`
groups in `control`, as well as basis and orbital inputs. If both attempts fail,
RIDFT is **never** started. The single pre-JOBEX
RIDFT runs on the displaced structure, with no separate SLURM allocation
or automatic SCF retry; convergence is checked before JOBEX. JOBEX also
performs its own RI-DFT SCF/gradient cycles as required. The original
Fermi-smearing initialization, pre-JOBEX EIGER, occupation restoration,
and integer-RIDFT stages are deliberately **not** repeated.

Unlike the normal initial DEFINE check, the Screwer DEFINE completion check
does **not** compare MO coefficients, shell occupations, or control contents
against the prior state. This is an intentional continuation; a successful
DEFINE must produce structurally usable input for the subsequent RIDFT.
`def_screwer_ln.inp` is used only if the primary DEFINE fails; it is not used
for the initial `lan-run` configuration.
The scientific correction runs **only once**; no
recursive automatic re-distortion is scheduled, even if the new Hessian has
remaining imaginary frequencies.

The supplied `screwer.inp` determines which mode is displaced and the
step/temperature. This program verifies the produced `$newcoord` comment and
refuses to submit if the selected mode is not negative. It never invents a
mode or displacement setting. Use **`lan-screwer ... --prepare-only` first**
for an unfamiliar project or input template. An already existing
`<metal>/screwer_<metal>/` is never overwritten. An existing legacy
`<metal>/screwer/` is preserved and does not block the new folder. If preparation fails after `cpc`,
inspect the directory manually rather than rerunning destructively.

**Small imaginary frequencies:** A negative value such as -8.65 cm⁻¹ in
one of the first six modes can be a numerical translation/rotation rather
than a chemically meaningful unstable coordinate. The default is to trigger
on *strictly negative* values as requested, but for an existing completed job
you can choose `--min-imag 20` to require a mode below -20 cm⁻¹:

```bash
lan-screwer num --metals gd --min-imag 20 --prepare-only
```

The original completed metal folder is never edited by screwer preparation.
Old copied status/scheduler files in the newly created child are archived
under `screwer_<metal>/parent_results/` (including old initial occupations
and Fermi-era control backups); inherited frequency files are archived
before the new Hessian. The child Hessian step attempts to extract unavailable
`$noproj` data groups with `sdg` to `nprhessian`, `nprvibspectrum`,
`npr_vib_normal_modes` without overwriting existing nonempty files. You can
repeat just this cheap extraction on the login node with:

```bash
lan-screwer --export-only --metals gd
```

Caveats: `def_screwer.inp` / `def_screwer_gd.inp` must preserve the inherited
MO files and occupations and must **not reintroduce `$fermi`**. The child
stages are separate jobs with their own chain records, not the original
parent's history. `lan-status` and `lan-status --details` include an
additional `<metal>/screwer_<metal>` row whenever that folder exists, with
independent SLURM, frequency and negative-mode results; historical
`<metal>/screwer` directories are also recognized for monitoring.
Review `<metal>/screwer_<metal>/chain_submission.tsv` and `screwer.out`
separately. An old four-job screwer preparation cannot be submitted under
the new two-job stage plan without manual review. Likewise, a previously
prepared child containing an older SLURM script without the pre-JOBEX RIDFT
is rejected; inspect and update its copied script manually if appropriate.
SLURM submissions launched from the short gate job and actual TURBOMOLE
`screwer`/`sdg`/`kdg`/DEFINE behavior require a controlled **VACC validation**.


## Development branch: turbomoleio integration and stage-aware recovery

The experimental additions below are **not validated on VACC/TURBOMOLE 7.9**.
Use them only on independent copies of calculations and keep your production
installation untouched until the tests have passed. Installing into
`$HOME/bin/lan-dev` keeps the commands separate from production.

### Input inspection

`lan-validate --project-only` inspects the **EF project root**; it does not
scan metal folders. `lan-validate --include-screwer` scans existing metal
folders and correction children. When invoked inside a metal folder it
recognizes the enclosing EF project. Repeated result data groups like `$last`,
`$vibrational`, and `$nprvibrational` are legitimate TURBOMOLE results and no
longer cause false FAILs. It is strictly read-only and is **not a proof** that
one has the correct electronic ground state.

### Universal recovery (experimental)

`lan-recover` retains its original `--metal ... --ridft-only|--full-chain|--continue-only`
behavior. New **read-only by default** recovery planning works from inside a
metal folder or its `screwer_<metal>` correction child:

```
lan-recover --plan                    # no modifications and no submissions
lan-recover --submit                  # submit the verified stages
lan-recover --plan --frequency ao     # for an untracked frequency method
```

From the EF root, target a folder explicitly:

```
lan-recover --plan --base "$PWD" --metal nd
lan-recover --plan --base "$PWD" --metal nd/screwer_nd
```

The planner uses the recorded Hessian method (or requires the user to name
it when unknown), checks SCF signatures, geometry convergence markers and
frequency evidence, and refuses to submit if SLURM still has queued or
running jobs in the same directory. In particular, a failed Screwer JOBEX
resumes with `jobex_eiger` **without DEFINE**. An SCF failure after
successful Screwer DEFINE resumes with `screwer_ridft_jobex_eiger`, which
also skips DEFINE. Hessian failures resume the selected Hessian alone.

When a JOBEX recovery starts on a compute node, it first archives existing
`jobex.out`/failure markers in `recovery_attempts/<jobid>/`, then removes
`not.converged` and stale geometry markers before rerunning JOBEX. Frequency
recovery similarly archives stale frequency outputs before rerunning. No
science files are overwritten merely by asking for a plan. If old dependent
jobs remain PENDING after a failed parent job, inspect/cancel those obsolete
jobs separately before using `--submit`—the tool will never cancel jobs on
your behalf. All submitted stages are journaled, with afterok dependencies.

For manual calculations without `metal.json`, **Hessian-only recovery from a
verified converged geometry** is supported with an explicitly selected method
(e.g., `lan-recover --plan --frequency ao`). The resumed Hessian uses the
existing control/orbitals, with `TM_MANUAL_RECOVERY` constrained to Hessian
stages only; the tool does not invent charge or spin metadata. Manual SCF or
JOBEX recovery without provenance is intentionally rejected.

### Opt-in restricted (closed-shell) La/Lu

The whole **complex** must admit an even-electron closed-shell state; a La(III)
or Lu(III) metal center alone is insufficient proof of that state. Provide a
user-authored `def_closedshell.inp` at the EF project root. It must implement
the actual *restricted* EHT dialogue in TURBOMOLE 7.9, contain one quoted
lanthanide basis assignment to replace, and set the EF-specific molecular
charge. It must not contain an unrestricted `u N` command. Do not simply
rename `def.inp`—the menus differ. No such template is generated here.

```
lan-setup --closed-shell --metals la lu
lan-run ao --closed-shell --metals lu --include-existing
# or one-step (on previously unprepared directories)
lan-run ao --closed-shell --metals la lu
```

This creates only `la_closedshell` and/or `lu_closedshell`, leaving the
existing `la/` and `lu/` untouched. The original control's DFT, RI, ECP/X2CSR
and dispersion settings remain in place, but open-shell and Fermi groups are
removed from the **new copy**. A positive even all-electron molecular count
is required; this is only a necessary precondition, not a multiplicity proof.
The new pipeline has `restricted_ridft` → `jobex_eiger` → optional Hessian,
without the open-shell occupation-removal stage. DEFINE must terminate
normally and produce `$closed shells`, `$scfmo`, and `mos`; otherwise RIDFT
will not start. These paths require a real VACC DEFINE test.

### Opt-in adaptive Screwer supervisor

The existing tested `def_screwer.inp` → `def_screwer_ln.inp` fallback remains
the default. For **new** correction directories, one may instead prepare with:

```
lan-screwer ao --metals pr --prepare-only --define-backend adaptive
```

The prepared manifest pins the backend; submission uses the prepared choice.
The experimental `pexpect` backend recognizes collinearity/geometry failures
while DEFINE is running, rather than waiting for an exit status. On an error,
the existing Screwer workflow restores the prepared files and retries the
no-`ired` template. This is a supervised implementation of the *existing*
user-defined menus, **not** a general replacement for initial ECP/X2CSR
DEFINE. The `turbomoleio` library is a separately installed dependency for
`lan-validate`; it is not vendored into this repository.

### Development acceptance tests

Run `python -m unittest discover -s tests -v` to exercise offline contracts.
Before considering a merge into production, also test in fresh copies on
VACC: La/Gd open-shell initial ECP and X2CSR, one restricted La/Lu system
with a verified charge, Pr Screwer fallback, Screwer RIDFT/JOBEX failure,
and Hessian-only recovery. Compare basis/ECP, occupations, method data
groups, convergence evidence and resulting energies with reference jobs.
No automated test can establish a unique physical electronic ground state.


### Development: bounded Screwer cycles and recovery plans

`lan-recover --plan` never modifies files. A `$fermi` entry in a Screwer
`control` produces a prominent warning, not a planning failure. Submission
retains the user's existing smearing choice; verify it before using `--submit`.
A verified Hessian with negative modes is **not** resubmitted by recovery; use
Screwer to create a separate correction cycle.

Screwer cycle 1 keeps its existing `METAL/screwer_METAL` directory. Cycle 2
is `METAL/screwer_METAL/screwer_METAL_cycle2`; cycle 3 is nested inside
cycle 2, and so on. Every cycle preserves the preceding results. A cycle
is only eligible after the previous cycle has completed its Hessian.

From the EF project root, manually start a further correction after cycle 1:

```bash
lan-screwer ao --metals lu --cycle 2 --prepare-only
lan-screwer ao --metals lu --cycle 2 --submit-prepared
```

For a **fresh cycle 1**, opt into bounded automatic continuation:

```bash
lan-screwer ao --metals lu --max-cycles 3 --prepare-only
lan-screwer ao --metals lu --max-cycles 3 --submit-prepared
```

Upon successful Hessian verification, each cycle checks the first eight
`vibspectrum` mode rows and, if a mode is below `--min-imag`, creates and
submits the next cycle. The cap is immutable in the prepared manifest. No
looping occurs by default (`--max-cycles` omitted means only the selected
cycle). Automatic scheduling errors are reported in `cycle-advance.log` and
do not invalidate an otherwise successful Hessian. This feature is
**development-only** until tested with TURBOMOLE 7.9 and SLURM on VACC.

Universal recovery now recognizes nested Screwer cycles as well; it still
never re-runs Screwer displacement or DEFINE automatically.
