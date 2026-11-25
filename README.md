# Lanthanide Automation

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A Python and Bash toolkit for automating lanthanide density functional theory (DFT) calculations with TURBOMOLE and SLURM on high-performance computing (HPC) clusters.

The workflow generates metal-specific inputs from a shared molecular template and runs calculations across 14 trivalent lanthanides (La–Lu, excluding Pm). It automates structure preparation, calculation setup, dependent job submission, output parsing, monitoring, and recovery of failed calculations. It is designed for reproducible studies of lanthanide complexes across molecular models and computational methods.

## Features

- Generate metal-specific coordinates, basis selections, spin settings, and calculation directories.
- Automate SCF calculations, orbital occupation handling, geometry optimization, and optional vibrational frequency calculations.
- Submit dependent SLURM jobs, track job IDs and convergence, and export status summaries to CSV.
- Resume interrupted submissions or recover failed calculation stages without repeating verified completed stages unnecessarily.
- Validate TURBOMOLE inputs and optionally correct imaginary frequencies in separate calculation directories.

## Installation

Requires Python 3, Bash, TURBOMOLE 7.9, and a SLURM-based HPC environment.

```bash
git clone https://github.com/Ahmad-Fermion/lanthanide_automation.git
cd lanthanide_automation
bash install.sh
export PATH="$HOME/bin:$PATH"
```

The installer adds six `lan-*` commands to `~/bin`. Add the `export` line to your shell startup file to make them available in future sessions.

The workflow was developed for the Vermont Advanced Computing Center (VACC) at the University of Vermont. For other HPC environments, adjust the SLURM account, partitions, resource settings, TURBOMOLE module, and scratch paths before installation and submission.

## Usage

### 1. Prepare a molecular template

Set up a representative lanthanide complex in TURBOMOLE. Record the interactive DEFINE session using `tee def.inp | define` and prepare the master `control` file. Your project directory should contain:

```text
project/
├── coord          # Molecular coordinates with one lanthanide atom
├── control        # Calculation settings, including $fermi nue
├── def.inp        # DEFINE input with metal basis, charge, and spin settings
└── def_gd.inp     # Gd-specific DEFINE input (required for Gd)
```

The workflow uses the supplied scientific settings and adapts the metal-specific inputs. Project templates are user-provided and should be checked for the chosen molecular system and method.

### 2. Prepare and submit calculations

Run from the project directory:

```bash
lan-setup --metals sm eu                     # Prepare inputs without submitting
lan-run ao --metals sm eu --include-existing # Submit inspected, pristine inputs
lan-run num --metals gd tb                   # Prepare and submit with NumForce
lan-run                                      # Run all remaining new metals, without a Hessian
```

`lan-run` defaults to geometry optimization without a Hessian. Use `ao` for analytical frequencies (`aoforce`) or `num` for numerical frequencies (`NumForce`). Existing calculation directories are skipped by default.

### 3. Monitor and recover calculations

```bash
lan-status                 # Summarize progress and write workflow_status.csv
lan-status --details       # Inspect SCF, geometry, and frequency status

cd sm
lan-recover --plan         # Inspect the recovery plan without changing files
lan-recover --submit       # Submit the remaining stages after reviewing the plan
```

### 4. Optional input validation and frequency correction

From the project directory, with `turbomoleio` installed in your Python environment:

```bash
lan-validate --project-only
lan-validate --metals sm eu
```

For a completed analytical Hessian calculation, provide `screwer.inp` and `def_screwer.inp` in the project directory, then run:

```bash
lan-screwer ao --metals sm --prepare-only
lan-screwer ao --metals sm --submit-prepared
```

Gd requires its own `def_screwer_gd.inp` continuation template. Review displaced inputs before submission.

## Testing

Run the offline unit tests from the repository root:

```bash
python3 -m unittest discover -s tests -v
```

Input validation and convergence checks assist workflow execution but do not replace scientific review of electronic states, geometries, and vibrational frequencies. Actual calculations require TURBOMOLE and a configured SLURM environment.
