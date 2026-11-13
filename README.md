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

The installer adds `lan-run`, `lan-recover`, `lan-setup`, and `lan-status` to `~/bin`. To update an existing installation:

```bash
cd ~/src/lanthanide_automation
git pull --ff-only
bash install.sh
```

**Requirements:** VACC account, Python 3, TURBOMOLE 7.9, and SLURM. 

## Four commands to remember

| Command | Purpose |
| --- | --- |
| `lan-run` | Prepare missing metal folders and submit the normal SLURM chain |
| `lan-recover` | Restart failed RIDFT with your manual `control` edits; optionally continue |
| `lan-setup` | Prepare folders **only**, without submitting |
| `lan-status` | Show progress; write `workflow_status.csv` |

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

`--include-existing` is for a prepared folder that has **not** already been chain-submitted; `chain_submission.tsv` prevents duplicate chain submission. Without this flag, existing folders are **skipped and never modified**.

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
