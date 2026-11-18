# parnassus_data_prep

ROOT macros that turn CLAS12 HIPO files into the flat ROOT trees used to train and evaluate the
Parnassus fast simulation of CLAS12 RG-A ([core repository](https://github.com/omphalos-mc/parnassus-clas12)). They read the reconstructed
particles and the generator-level (`MC::Lund`) particles of `clasdis` Monte Carlo, apply the
particle cuts listed below, and write one `evt_tree` entry per event.

The tag `paper-v1` is the version that goes with the paper. No data files are included; the
`clasdis` HIPO files are available to CLAS Collaboration members on the JLab farm.

## Layout

```text
macros/               ROOT entry points; they #include the analysis code in src/
src/                  particle loading (CLAS12Ana), cuts (CutManager), kinematics, tree writers
scripts/              SLURM jobs for the JLab farm, the merge wrapper, a plotting script
paper_plotting.ipynb  reconstructed vs generated theta and phi figure
```

## Environment

The macros need ROOT with the CLAS12 HIPO reader (`clas12root`). The v8 samples were made on the
JLab farm (ifarm) with clas12root (v6.38.04). There is no build step: the macros are
run with  `clas12root -l -b -q` and compiled on the fly.

The SLURM scripts contain JLab farm paths, the `clas12` account and the author's user directories;
edit them before running anywhere else.

## How the v8 samples were made

The input is `clasdis` Monte Carlo for RG-A Fall 2018 inbending, pass-2 reconstruction, 45 nA
background merging, 10.604 GeV beam:

```text
/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/clasdis_rga_fa18_inb_45nA_10604MeV-*.hipo
```

1. **Convert each HIPO file.** `macros/hipo2tree_parnassus.C` reads one HIPO file and writes
   `reco_EventTree`, `mc_EventTree` and `evt_tree` (the tree Parnassus reads) to
   `<output dir>/<file stem>.root`. `scripts/submit_hipo2tree_parnassus.sbatch` runs it as a
   SLURM array with one task per file; set `--array` to the number of files:

   ```bash
   sbatch scripts/submit_hipo2tree_parnassus.sbatch <directory with clasdis*.hipo> v8
   ```

   For a single file:

   ```bash
   clas12root -l -b -q 'macros/hipo2tree_parnassus.C("<file>.hipo","prepped_data",-1)'
   ```

   The last argument is the number of events (`-1` for all). A file whose path contains `reco` is
   treated as real data, and no generator-level particles are written for it.

2. **Merge.** `scripts/merge_evt_tree.sh` (see [Merge script](#merge-script)) concatenates
   `evt_tree` from a list of files and renumbers `eventNumber` sequentially with
   `macros/reassign_event_number.C`, since the raw event numbers are only unique within one HIPO
   file. With `--filterScatE` it keeps only events whose reconstructed scattered electron passes
   the cuts (`pflow_scatE == 1`). `scripts/submit_merge.sbatch` runs the merge for a train and a
   test file list on the farm.

| File | Events | Merge | Hipo Files
|---|---|---|---|
| `v8_train.root` | 59.0M | `--filterScatE`: every event has a reconstructed scattered electron | files 2-31 |
| `v8_test.root` | 29.91M | no filter| files 1,32,33 |
| `pre-filter_train_data_v8.root` | 19.91M | no filter | files 2,3 |


## Cuts

`CutManager::filter_particles` (`src/CutManager.C`) cuts the reconstructed particles before they
are written.

Every particle must satisfy:

- 5° < θ < 35° (Forward Detector acceptance)
- p ≤ 10.6 GeV
- transverse vertex √(v<sub>x</sub>² + v<sub>y</sub>²) < 20 cm

On top of that, per species, where "v<sub>z</sub>" means −40 cm < v<sub>z</sub> < 30 cm and "β" means 0 < β < 2:

| Particle | Cuts |
|---|---|
| e⁻ | DC fiducial (pass-2 edge distances), pass-2 sampling fraction, PCAL edges (l<sub>v</sub>, l<sub>w</sub> > 9 cm; l<sub>v</sub> > 13.5 cm in sector 4), β, v<sub>z</sub> |
| e⁺ | β, v<sub>z</sub> |
| π⁺, π⁻ | DC fiducial, pion χ²<sub>pid</sub> cut (−3C < χ²<sub>pid</sub> < 3C, tighter above 2.44 GeV; C = 0.88 for π⁺, 0.93 for π⁻), v<sub>z</sub> |
| K⁺, K⁻, p | DC fiducial, abs(χ²<sub>pid</sub>) < 3, v<sub>z</sub> |
| p̄ | abs(χ²<sub>pid</sub>) < 3, v<sub>z</sub> |
| γ | PCAL edges (l<sub>v</sub>, l<sub>w</sub> > 14 cm), PCAL energy > 0, β, v<sub>z</sub> |
| n | v<sub>z</sub> |
| PID 0 and 45 | removed |
| any other PID | only the cuts for every particle |

The reconstructed scattered electron is the highest-energy electron, chosen before the cuts. If
the cuts remove it, the event's DIS variables are reset to −999 and `pflow_scatE` is 0.

At generator level (`CLAS12Ana::load_mc_particles`), the final-state `MC::Lund` particles are kept,
except electron neutrinos and particles with θ < 1°.

## Merge script

[scripts/merge_evt_tree.sh](scripts/merge_evt_tree.sh) merges only `evt_tree` from an explicit
list of ROOT files into one output file:

```bash
scripts/merge_evt_tree.sh OUTPUT_FILE [--dry-run] [--max-files N] [--force] [--filterScatE] FILE1 [FILE2 ...]
```

- `OUTPUT_FILE`: the ROOT file to create.
- `--dry-run`: print the input file list and exit.
- `--max-files N`: merge only the first `N` input files.
- `--force`: overwrite an existing output file.
- `--filterScatE`: keep only events with `pflow_scatE == 1`.
- `FILE1 [FILE2 ...]`: the input ROOT files; shell globs such as `prepped_data/*.root` work.

For example:

```bash
scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --dry-run prepped_data/*.root
scripts/merge_evt_tree.sh prepped_data/merged_evt_tree.root --force --filterScatE prepped_data/*.root
```

## Other tools

- `macros/plot_prepped_data.C`: reconstructed vs generated kinematics of one prepped file, written
  to `plots/<file stem>/` as PNG, PDF and the histograms in ROOT files.
- `paper_plotting.ipynb`: restyles the θ and φ comparison from
  `plots/<sample>/compare_theta.root` and `compare_phi.root` (written by `plot_prepped_data.C`)
  for the paper.
- `scripts/plot_kinematics_per_species.py`: reconstructed vs generated p, θ, φ, v<sub>z</sub> and β per
  species (needs `uproot`, `awkward` and `matplotlib`).
- `macros/filter_electrons.C`: drops every event with an electron or positron in `pflow_pdgId`.

Output under `prepped_data/` and `plots/` is not tracked, and neither are ROOT's compiled-macro
files (`*.so`, `*.d`, `*.pcm`).

## License

MIT; see [LICENSE](LICENSE).
