#!/usr/bin/env bash
set -euo pipefail

# Merge only evt_tree from an explicit list of ROOT files.
#
# Usage:
#   scripts/merge_evt_tree.sh OUTPUT_FILE [--dry-run] [--max-files N] [--force] FILE1 [FILE2 ...]
#
# Example:
#   scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --dry-run --max-files 10 prepped_data/*.root
# srun --partition=production --account=clas12 --mem-per-CPU=5000 --nodes=1 --chdir=/w/hallb-scshelf2102/clas12/users/tjhellst/parnassus_data_prep --output=/farm_out/%u/%x-%j-%N.out --error=/farm_out/%u/%x-%j-%N.err --job-name==merge_files_train --time=2:00:00 

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 OUTPUT_FILE [--dry-run] [--max-files N] [--force] FILE1 [FILE2 ...]" >&2
  exit 1
fi

OUTPUT_FILE="$1"
shift 1

FORCE=false
DRY_RUN=false
MAX_FILES=-1

INPUT_FILES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-files)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --max-files requires a numeric value." >&2
        exit 1
      fi
      MAX_FILES="$2"
      shift 2
      ;;
    --max-files=*)
      MAX_FILES="${1#--max-files=}"
      shift 1
      ;;
    --force)
      FORCE=true
      shift 1
      ;;
    --dry-run)
      DRY_RUN=true
      shift 1
      ;;
    --)
      shift 1
      while [[ $# -gt 0 ]]; do
        INPUT_FILES+=("$1")
        shift 1
      done
      ;;
    *)
      INPUT_FILES+=("$1")
      shift 1
      ;;
  esac
done

if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
  echo "ERROR: Provide at least one input ROOT file." >&2
  exit 1
fi

ROOT_FORCE=false
if [[ "$FORCE" == true ]]; then
  ROOT_FORCE=true
fi

if [[ "$MAX_FILES" != "-1" ]]; then
  if ! [[ "$MAX_FILES" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --max-files must be a non-negative integer." >&2
    exit 1
  fi
  if [[ "$MAX_FILES" -ge 0 && ${#INPUT_FILES[@]} -gt "$MAX_FILES" ]]; then
    INPUT_FILES=("${INPUT_FILES[@]:0:$MAX_FILES}")
  fi
fi

if [[ "$DRY_RUN" == true ]]; then
  echo "Dry run: would merge the following ${#INPUT_FILES[@]} file(s) into $OUTPUT_FILE:"
  for input_file in "${INPUT_FILES[@]}"; do
    echo "$input_file"
  done
  exit 0
fi

if [[ -e "$OUTPUT_FILE" && "$FORCE" != true ]]; then
  echo "ERROR: Output file exists. Re-run with --force to overwrite: $OUTPUT_FILE" >&2
  exit 1
fi

MANIFEST_FILE="$(mktemp)"
trap 'rm -f "$MANIFEST_FILE"' EXIT

for input_file in "${INPUT_FILES[@]}"; do
  if [[ ! -f "$input_file" ]]; then
    echo "ERROR: Not a regular file: $input_file" >&2
    exit 1
  fi
  printf '%s\n' "$input_file" >> "$MANIFEST_FILE"
done

root -l -b -q "macros/merge_evt_tree.C(\"${MANIFEST_FILE}\",\"${OUTPUT_FILE}\",false,${ROOT_FORCE})"

# The raw eventNumber copied from each input file is only unique within that
# file, so the merged evt_tree contains duplicate event numbers. Reassign it
# sequentially (writing to a temp file, then moving it over the merged output).
REINDEX_FILE="${OUTPUT_FILE}.reindex.tmp"
rm -f "$REINDEX_FILE"
root -l -b -q "macros/reassign_event_number.C(\"${OUTPUT_FILE}\",\"${REINDEX_FILE}\")"
mv -f "$REINDEX_FILE" "$OUTPUT_FILE"
