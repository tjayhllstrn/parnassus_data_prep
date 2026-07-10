# parnassus_data_prep

This repository contains a CLAS12 data-preparation pipeline for turning raw `.hipo` event files into flat ROOT `TTree`s for downstream analysis and ML workflows, especially the Parnassus particle-flow pipeline.

The code is organized as ROOT macros in `macros/` that `#include` the analysis implementation from `src/`. There is no separate build system; the macros are run directly with `root -l -b -q`, typically inside the CLAS12 ROOT/apptainer environment.

## Repository purpose

The main pipeline reads reconstructed and Monte Carlo truth particle banks from CLAS12 HIPO files, computes DIS/SIDIS kinematics, applies detector and PID cuts, and writes ROOT trees for later plotting, validation, and machine-learning use.

Key outputs include:
- `reco_EventTree` for reconstructed particles
- `mc_EventTree` for Monte Carlo truth particles
- `evt_tree` for the ML-facing Parnassus tree

Supporting macros in this repo also provide quick validation and plotting of the prepped ROOT files.

## Merge script

The script [scripts/merge_evt_tree.sh](scripts/merge_evt_tree.sh) merges only `evt_tree` from an explicit list of ROOT files into a single output ROOT file.

It is intended for cases where you already know exactly which files should be merged and want a simple, safe wrapper around the ROOT macro.

### Syntax

```bash
scripts/merge_evt_tree.sh OUTPUT_FILE [--dry-run] [--max-files N] [--force] FILE1 [FILE2 ...]
```

### Arguments

- `OUTPUT_FILE`: destination ROOT file to create.
- `--dry-run`: print the resolved input file list and exit without creating output.
- `--max-files N`: only merge the first `N` input files after shell expansion.
- `--force`: overwrite an existing output file.
- `FILE1 [FILE2 ...]`: explicit input ROOT files to merge. Shell globbing is supported, so `prepped_data/*.root` works.

### Examples

Preview the exact files that would be merged:

```bash
scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --dry-run prepped_data/*.root
```

Merge the first 10 files from a globbed list:

```bash
scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --max-files 10 --force prepped_data/*.root
```

Merge a hand-picked list of files:

```bash
scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --force \
  prepped_data/file1.root \
  prepped_data/file2.root \
  prepped_data/file3.root
```

## ROOT environment

Run ROOT macros inside the CLAS12 analysis container or an environment that provides `root` and the CLAS12 libraries. The repository was written to be used with the clas12root/apptainer setup described in `CLAUDE.md`.

## Other useful macros

- `macros/hipo2tree_parnassus.C` converts a HIPO file into the standard prepped ROOT trees.
- `macros/check_classID.C` checks `evt_tree` for unknown particle class assignments.
- `macros/plot_prepped_data.C` compares reconstructed and truth kinematics from a prepped file.

## Notes

- Generated data products under `prepped_data/` and `plots/` are treated as reproducible artifacts.
- ACLiC build products such as `.d`, `.pcm`, and `.so` files are normal ROOT byproducts and can be regenerated.
