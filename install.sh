#!/usr/bin/env bash
# Install one shared copy of the automation and short launch commands.
set -Eeuo pipefail
source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
bin_dir="${1:-$HOME/bin}"
mkdir -p -- "$bin_dir"
bin_dir="$(cd -- "$bin_dir" && pwd -P)"
dest="$bin_dir/lanthanide_automation"
mkdir -p -- "$dest"
for name in README.md setup_lanthanides.py submission_guard.py monitor_lanthanides.py tm_workflow.py screwer_workflow.py validate_lanthanides.py recover_universal.py closed_shell_workflow.py adaptive_define.py \
            submit-turbo-chain.sh submit-turbo-recovery.sh \
            turbo.pipeline.slurm turbo.resume.slurm turbo.screwer.slurm install.sh; do
    if [[ ! "$source_dir/$name" -ef "$dest/$name" ]]; then
        cp -- "$source_dir/$name" "$dest/$name"
    fi
done
make_launcher() {
    local name="$1" runner="$2" file="$3" path="$bin_dir/$1"
    printf '#!/usr/bin/env bash\nexec %s %q "$@"\n' "$runner" "$dest/$file" > "$path"
    chmod 755 "$path"
}
make_launcher lan-run bash submit-turbo-chain.sh
make_launcher lan-recover bash submit-turbo-recovery.sh
make_launcher lan-setup python3 setup_lanthanides.py
make_launcher lan-status python3 monitor_lanthanides.py
make_launcher lan-screwer python3 screwer_workflow.py
make_launcher lan-validate python3 validate_lanthanides.py
printf 'Installed one shared copy: %s\nCommands: lan-run | lan-recover | lan-setup | lan-status | lan-screwer | lan-validate\n' "$dest"
if [[ ":$PATH:" != *":$bin_dir:"* ]]; then
    printf 'Add to your shell startup file (~/.bashrc), then open a new shell:\n  export PATH="%s:$PATH"\n' "$bin_dir"
fi
