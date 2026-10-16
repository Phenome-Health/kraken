#!/bin/bash
# Download the archived build inputs from S3 and unpack them, ready for `build.py`.
#
#   scripts/restore_inputs.sh [BASE_DIR]
#
# Puts the files in BASE_DIR/input_data (default BASE_DIR: the kraken repo root). The build reads
# inputs from `base_path` in config/build_config.yaml, so set that to BASE_DIR afterwards (or to
# null when BASE_DIR is the repo root). Uses the AWS profile in $AWS_PROFILE (default kraken-admin).
#
# The archive (s3://kraken-kg-archive/input_data) stores large text files compressed to save
# transfer time: babel and translator-kg-open as .zst, robokop and kg2 as their original .gz
# downloads. This puts back the uncompressed files the build reads, keeping the compressed copies.
# Needs the AWS CLI and zstd (`brew install zstd`). Safe to re-run: files already downloaded or
# unpacked are skipped. Expect ~250 GB of free disk space.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/.." && pwd)"
BASE_DIR="${1:-${REPO}}"
DEST="${BASE_DIR}/input_data"
PROFILE="${AWS_PROFILE:-kraken-admin}"

if ! command -v zstd &> /dev/null; then
    echo "zstd not found; install it first (e.g. brew install zstd)." >&2
    exit 1
fi

echo "Downloading s3://kraken-kg-archive/input_data to ${DEST} ..."
aws s3 sync s3://kraken-kg-archive/input_data "${DEST}" --profile "${PROFILE}" --region us-west-2

echo "Unpacking .zst files ..."
find "${DEST}" -name '*.zst' | while read -r f; do
    out="${f%.zst}"
    if [[ ! -f "${out}" ]]; then
        echo "  ${out#"${DEST}"/}"
        zstd -d -q -T0 "${f}" -o "${out}"
    fi
done

echo "Unpacking robokop and kg2 .gz files ..."
find "${DEST}/robokop" "${DEST}/kg2" -name '*.jsonl.gz' 2> /dev/null | while read -r f; do
    out="${f%.gz}"
    if [[ ! -f "${out}" ]]; then
        echo "  ${out#"${DEST}"/}"
        gunzip -c "${f}" > "${out}.tmp" && mv "${out}.tmp" "${out}"
    fi
done

# Babel reads NCBI's taxonomy from its own release directory; the archive keeps a single copy,
# under ncbigene/, so link each Babel release to it.
for release in "${DEST}"/babel/*/; do
    if [[ ! -e "${release}taxdump.tar.gz" && -f "${DEST}/ncbigene/taxdump.tar.gz" ]]; then
        echo "Linking ${release#"${DEST}"/}taxdump.tar.gz -> ncbigene/taxdump.tar.gz"
        ln -s ../../ncbigene/taxdump.tar.gz "${release}taxdump.tar.gz"
    fi
done

echo ""
echo "Done. Inputs are in ${DEST}."
if [[ "${BASE_DIR}" == "${REPO}" ]]; then
    echo "Set base_path: null in config/build_config.yaml, then run: uv run python build.py"
else
    echo "Set base_path: \"${BASE_DIR}\" in config/build_config.yaml, then run: uv run python build.py"
fi
