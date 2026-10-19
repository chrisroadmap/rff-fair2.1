#!/usr/bin/env bash
# Zip stochastic output NC files in batches of 1000:
#   output0.zip = run00001.nc .. run01000.nc
#   output1.zip = run01001.nc .. run02000.nc, etc.
set -euo pipefail

# Load DATADIR from .env in the git repo root, if present
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)"
if [[ -f "${REPO_ROOT}/.env" ]]; then
    set -a
    # shellcheck disable=SC1091
    source "${REPO_ROOT}/.env"
    set +a
fi

: "${DATADIR:?DATADIR is not set (set it in the environment or in ${REPO_ROOT}/.env)}"

INDIR="${DATADIR}/data_output/stochastic"
OUTDIR="${1:-$INDIR}"   # optional: directory for zip files (default: same as input)
NFILES=10000
BATCH=1000

cd "$INDIR"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

for ((b = 0; b < NFILES / BATCH; b++)); do
    start=$((b * BATCH + 1))
    end=$(((b + 1) * BATCH))
    files=()
    for ((i = start; i <= end; i++)); do
        f=$(printf "run%05d.nc" "$i")
        [[ -f "$f" ]] || { echo "Missing file: $f" >&2; exit 1; }
        files+=("$f")
    done
    echo "Creating output${b}.zip (run$(printf '%05d' $start).nc - run$(printf '%05d' $end).nc)"
    zip -q -j "${OUTDIR}/output${b}.zip" "${files[@]}"
done

echo "Done."
