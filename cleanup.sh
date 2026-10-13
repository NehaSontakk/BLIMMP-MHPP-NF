#!/bin/bash
# cleanup.sh — remove Nextflow run artifacts from the current directory
# Targets: work/, results dirs, .nextflow/, .nextflow.log*, nextflow-*.out/err

set -euo pipefail

# Collect items to delete
to_delete=()

# Nextflow work cache
[[ -d work ]]        && to_delete+=(work/)
[[ -d .nextflow ]]   && to_delete+=(.nextflow/)

# Nextflow log files
for f in .nextflow.log .nextflow.log.*; do
    [[ -e "$f" ]] && to_delete+=("$f")
done

# SLURM output/error files from nextflow head jobs
for f in nextflow-*.out nextflow-*.err; do
    [[ -e "$f" ]] && to_delete+=("$f")
done

# Results directories (any directory matching results*)
for d in results*/; do
    [[ -d "$d" ]] && to_delete+=("$d")
done

if [[ ${#to_delete[@]} -eq 0 ]]; then
    echo "Nothing to clean up."
    exit 0
fi

echo "The following will be deleted:"
echo ""
for item in "${to_delete[@]}"; do
    du -sh "$item" 2>/dev/null | awk '{printf "  %-10s %s\n", $1, $2}' || echo "  ?          $item"
done
echo ""

read -r -p "Delete all of the above? [y/N] " answer
case "$answer" in
    [yY][eE][sS]|[yY])
        for item in "${to_delete[@]}"; do
            echo "Removing $item ..."
            rm -rf "$item"
        done
        echo "Done."
        ;;
    *)
        echo "Aborted."
        exit 0
        ;;
esac
