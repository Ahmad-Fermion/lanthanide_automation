# Lanthanide automation · TURBOMOLE 7.9 on VACC

Prepare, submit, monitor, and recover **14 Ln³⁺ calculations (La–Lu, excluding Pm)** using a shared command-line toolset. Developed for **UVM's Vermont Advanced Computing Center (VACC)** with TURBOMOLE 7.9 and SLURM. The multi-stage workflow can be adapted to other HPC clusters by adjusting their SLURM resources, software modules, and scratch configuration.

Run the commands from a project directory containing `control`, `coord`, and `def.inp` (and `def_gd.inp` when preparing Gd). Each project's input files define the molecular system and theory; the tools can be reused across multiple models, methods (ECP/X2CSR), and binding sites.

## Install once on VACC (from GitLab)

Copy the **SSH clone URL** from this repository's GitLab page, then on VACC:

```bash
mkdir -p ~/src ~/bin
git clone <SSH-CLONE-URL> ~/src/lanthanide_automation
bash ~/src/lanthanide_automation/install.sh
export PATH="$HOME/bin:$PATH"   # add to ~/.bashrc only if ~/bin is not already in PATH
```

For a private repository, the VACC account also needs GitLab SSH-key access (separate from your Mac). The installer adds `lan-run`, `lan-recover`, `lan-setup`, and `lan-status` to `~/bin`. You don't need to copy the tools into every calculation project. To update an existing installation:

```bash
cd ~/src/lanthanide_automation
git pull --ff-only
bash install.sh
```

**Requirements:** VACC account, Git access, Python 3, TURBOMOLE 7.9, and SLURM. Installation does not submit jobs.

## Four commands to remember

| Command | Purpose |
| --- | --- |
| `lan-run` | Prepare missing metal folders and submit the normal SLURM chain |
| `lan-recover` | Restart failed RIDFT with your manual `control` edits; optionally continue |
| `lan-setup` | Prepare folders **only**, without submitting |
| `lan-status` | Show progress; write `workflow_status.csv` |

Choose `num` = **NumForce** (for X2CSR) or `ao` = **aoforce**. Full names `numforce`/`aoforce` also work. Default is `num`.

## Common examples

Run from a single EF project directory, e.g.:

```bash
cd ~/lanM_shared/small-model/x2csr-d4-wd-scf7/EF1738
lan-run num --metals sm --prepare-only      # create only Sm for inspection
lan-run num --metals sm --include-existing  # submit that prepared Sm
lan-run num                                 # submit all other *missing* metals
lan-run ao                                  # same workflow, ending with aoforce
lan-run ao --metals gd tb                   # only new Gd/Tb (Gd uses def_gd.inp)
lan-status                                  # progress + workflow_status.csv
```

`--include-existing` is for a prepared folder that has **not** already been chain-submitted; `chain_submission.tsv` prevents duplicate chain submission. Without this flag, existing folders are **skipped and never modified**.

## When RIDFT needs manual intervention

After the original failed/blocked jobs are stopped, edit e.g., `tb/control` manually. From the EF project directory:

```bash
lan-recover num --metal tb                  # run only RIDFT; inspect state
lan-recover num --metal tb --continue-only  # later, if recovery converged
```

Or, if confident in your changes:

```bash
lan-recover num --metal tb --full-chain     # resume RIDFT + all later stages
```

## What happens automatically

`define < def.inp` → `ridft` → `eiger` → occupation fix → second `ridft` → `jobex` → `eiger` → `NumForce`/**`aoforce`**. Each expensive stage has its own **48-hour `general`** allocation; EIGER takes 10 minutes; the occupation edit gets **10 minutes, 1 CPU, 1 GB on `short`**. Jobs use `afterok` dependencies.

After `define`, the initial integer alpha/beta shell assignments are saved in `initial_occupations`. The occupation-fix stage restores those blocks and **removes `$fermi`** before the second RIDFT (existing behavior). Orbital holes are **warnings**; true SCF or geometry failure blocks downstream jobs. `control` is the only TURBOMOLE control file; no `control.template`. The validated NumForce OpenMP environment and `PARNODES=1 NumForce -ri -central` remain unchanged.

**Requirements:** Python 3, SLURM, TURBOMOLE 7.9 on VACC. For new `gd/`, provide `def_gd.inp` in the project root. Defaults: 20 CPUs/20 GB for RIDFT/JOBEX, 100 GB for frequency; overrides via `TM_CPUS`, `RIDFT_MEM`, `JOBEX_MEM`, `HESS_MEM`, `EIGER_CPUS`, `EIGER_MEM`.

**Help:** `lan-run --help`, `lan-recover --help`, `lan-setup --help`, `lan-status --help`. The per-metal files (including `chain_submission.tsv`, `initial_occupations`, stage logs and warnings) remain in that EF project. Verify orbitals/occupations and frequencies scientifically.
