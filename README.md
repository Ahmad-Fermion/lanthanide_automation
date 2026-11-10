# Lanthanide automation · TURBOMOLE 7.9 on VACC

One shared toolset for **14 Ln3+ metals (La–Lu, excluding Pm)**. 
Run it from **any one of your 16 EF project directories** (small/large model × ECP/X2CSR × EF1738–EF1741). 
The project's own `control`, `coord`, `def.inp`, and Gd `def_gd.inp` determine the level of theory.

## Install once in `~/bin`

Upload the ZIP to your home directory, then:

```bash
mkdir -p ~/bin
unzip -o ~/lanthanide_automation.zip -d ~/bin
bash ~/bin/lanthanide_automation/install.sh
export PATH="$HOME/bin:$PATH"   # add to ~/.bashrc if needed for future logins
```

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
