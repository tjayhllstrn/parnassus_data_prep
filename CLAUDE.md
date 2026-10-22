# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository does

This is a CLAS12 (Jefferson Lab Hall B) data-preparation pipeline. It converts raw `.hipo` event files
(read via `clas12root`/`hipo` banks) into flat ROOT `TTree`s, with a specific focus on producing input
trees for "Parnassus" — a particle-flow-style ML pipeline. For every event it loads reconstructed (reco)
and Monte Carlo (MC/truth) particles, computes DIS/SIDIS kinematics, applies detector fiducial/PID cuts,
and writes everything out as ROOT trees for downstream plotting/ML use.

There is no build system (no Makefile/CMake) and no test suite — this is a ROOT/CINT macro-based analysis
codebase, run interactively or in batch via `root -l -b -q`.

## Running the code

All code in this repo must be run inside the clas12root apptainer container — the CLAS12/ROOT/hipo
environment is not available on the bare cluster filesystem. Enter the container with:
```
apptainer exec /hpc/group/vossenlab/niw5/apptainer_images_new/analysis_latest.sif
```
Run this first (or prefix commands with it) before invoking any `root` command below; ROOT/`clas12root`
are only on PATH inside the container.

All entry points are ROOT macros in `macros/`. They `#include` the corresponding `.C` sources from `src/`
directly (ROOT's ACLiC/CLING interprets these), so there is nothing to "build" separately — just invoke
the macro with `root`.

Convert a hipo file to ROOT trees (main pipeline):
```
root -l -b -q 'macros/hipo2tree_parnassus.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001_small_version.root",-1)'
```
The third argument is `maxEvents` (`-1` means process the whole file).

Plot/compare reco vs. truth kinematics from a prepped file:
```
root -l -b -q 'macros/plot_prepped_data.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001_small_version.root")'
```

Sanity-check the particle-class assignment (`pflow_class`/`truth_class`) in a prepped file for unknown PIDs:
```
root -l -b -q 'macros/check_classID.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001_small_version.root")'
```

Measure what each particle-level cut does on its own — for every cut in the registry, apply only that cut
and count how many events end up with no surviving particles (`npflow==0`):
```
root -l -b -q 'macros/cut_study.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","plots/cut_study",20000)'
```
This one must read a raw `.hipo`, not a prepped `.root` — the prepped trees are filled *after*
`filter_particles`, so the pre-cut particles only exist in the hipo. The set of cuts studied is a single
registry block at the top of `macros/cut_study.C` (one line per cut); the counting machinery lives in
`src/CutStudy.h`. Results go to stdout plus `<outDir>/<hipo stem>/summary.csv` and `removed_by_pid.csv`.

`hipo2tree_parnassus.C` writes the v7 `evt_tree` layout directly. To convert processed files made before that
(e.g. v5) to v7 without rerunning the hipo step, build it from their `reco_EventTree`/`mc_EventTree`
(adds `pflow_p`/`truth_p`, unfolds `truth_theta`, applies the e+/e- beta cut).
`"train"` keeps only events with a reco scattered electron; `"test"` keeps every event:
```
root -l -b -q 'macros/make_v7.C+("<processed dir>","v7_train.root","train",-1,50000000)'   // nFiles=-1 (all), nEvents
root -l -b -q 'macros/make_v7.C+("<processed dir>","v7_test.root","test",-1,5000000,200000000)'  // last: first eventNumber
```

These all require running inside the apptainer container described above.

ACLiC compilation artifacts (`macros/*.d`, `macros/*.pcm`, `macros/*.so`) are build byproducts of running a
macro once (e.g. `plot_prepped_data.C`) — they are gitignored and safe to delete/regenerate.

## Architecture

### Pipeline shape

`hipo2tree_parnassus.C` is the main driver. Per event it:
1. `CLAS12Ana::get_event_info` — pulls run/event/torus/helicity/target/polarization metadata (via `Constants.h` run-range lookup tables).
2. `CLAS12Ana::load_reco_particles` — builds a `std::vector<part>` from the REC banks, enriched with calorimeter/DC/Cherenkov info via `HipoBankInterface`.
3. `CLAS12Ana::find_reco_scattered_electron` — identifies the scattered electron in the reco particle list.
4. `CLAS12Ana::fill_reco_event_variables` — computes DIS kinematics (x, y, Q2, W, nu) into an `EVENT` struct.
5. `CutManager::filter_particles` — applies per-particle PID/fiducial/momentum cuts, returning a filtered `std::vector<part>`.
6. `CLAS12Ana::load_mc_particles` / `fill_mc_event_variables` — same as above but for the Lundified MC truth bank.
7. Three tree writers (`EventTree` for reco, `EventTree` for MC, `ParnassusTree` combining both) are filled and written to the output ROOT file.

### Key classes (`src/`)

- **`CLAS12Ana`** (`CLAS12Ana.h/.C`) — the core per-event analysis driver: loads reco/MC particles, finds the
  scattered electron, fills DIS and dihadron kinematic variables, and matches MC truth particles to reco
  particles (`match_mc_to_reco`, `match_mc_to_reco_mcmatch`). Also contains dihadron combinatorics helpers
  (`dihadron_idxs`, `unique_combinations`, `remove_duplicates`) for building hadron pairs from a particle list.
- **`CutManager`** (`CutManager.h/.C`) — encapsulates all detector/PID fiducial cuts (DC fiducial cuts in
  XY and theta/phi parameterizations, calorimeter sampling-fraction cuts, chi2pid, vertex cuts, etc.),
  keyed off `RUN_PERIOD` (RGA/RGB/RGC/RGK/MC) and torus bending, determined from the run number. Entry
  point is `filter_particles(vector<part>)`.
- **`HipoBankInterface`** (`HipoBankInterface.h/.C`) — low-level glue that reads CLAS12 hipo banks
  (`REC::Calorimeter`, `REC::Cherenkov`, `REC::Traj`, etc.) by `pindex` and copies calorimeter/drift-chamber/
  Cherenkov info into a `part` struct (`loadBankData`).
- **`Kinematics`** (`Kinematics.h/.C`) — stateless static-method library for standard DIS/SIDIS kinematic
  formulas (Q2, x, y, W, xF, phi_h, phi_R, depolarization factors, etc.), built on `TLorentzVector`.
- **`TreeManager.C`** — defines the output tree schemas: `BaseTree` (abstract), `EventTree` (per-particle
  reco or MC branches, fixed-size arrays up to `kNmax=100`), `DihadronTree` (event-level dihadron kinematics),
  and `ParnassusTree` (the ML-facing tree: `pflow_*` reco particles + `truth_*` MC particles, each assigned a
  coarse `class` id via `GetClassID` — 0=charged hadron, 1=electron, 2=muon, 3=neutral hadron, 4=photon, -1=unknown).
  `ParnassusTree` writes the v7 layout: `pflow_pt`/`pflow_p`, `truth_pt`/`truth_p`, and `truth_theta` as the
  true polar angle in [0, pi] (`atan2(truept, truepz)`; `Kinematics::th` is `atan2` too).
- **`Constants.h`** — header-only constants and lookup tables: run-period classification (`runPeriod`),
  beam polarization/energy/torus-bending by run range (`runPolarization`, `runBeamEnergy`, `runTorusBending`),
  and large hardcoded fit-parameter tables for DC fiducial cuts and calorimeter sampling-fraction cuts,
  split by inbending/outbending torus setting.
- **`Structs.h`** — plain data structs shared across the codebase: `part` (one reconstructed or MC particle,
  with parallel `true*` fields for MC truth and detector-subsystem fields), `EVENT` (per-event DIS/dihadron
  kinematics, with parallel `true*` fields), `EVENT_INFO` (run/target/polarization metadata), and `FS`
  (final-state hadron pair specification, with pi0-aware logic in `get_FS`).
- **`style.C`** — ROOT plotting style helpers (colors, `TStyle` setup) used only by `plot_prepped_data.C`.

### Data flow / naming conventions

- Every kinematic/particle field has an MC-truth counterpart prefixed with `true` (e.g. `Q2`/`trueQ2`,
  `px`/`truepx`). This convention runs through `part`, `EVENT`, and the tree branches.
- `part.pindex` is the index into the original hipo REC::Particle bank; MC truth matching is done by
  `mcindex`. Dihadron code refers to particle pairs by index pairs (`i`,`ii`,`j`,`jj` in `EVENT`).
- Target codes (`EVENT_INFO.target`) are fixed integers: 0=NH3, 1=ND3, 2=C, 3=CH2, 4=CD2, 5=Empty.
- Run-number ranges are the primary key used throughout `Constants.h`/`CutManager` to select run-period-
  and torus-bending-dependent physics constants and cuts — when adding support for a new run period, add
  the run range to `runPeriod`, `runPolarization`, `runBeamEnergy`, `runTorusBending`, and (if fiducial
  cuts differ) the `CutManager` cut logic.
- Individual cuts in `CutManager` are `protected` on purpose — production code goes through
  `filter_particles` only. `src/CutStudy.h` defines a `ProbeableCutManager` subclass that re-exposes them
  for the cut study, so a cut you want to be able to study must be a named method (not an inline condition
  inside `filter_particles`) and needs a `using` line in that subclass.
- `prepped_data/` and `plots/` hold generated `.hipo`/`.root`/`.png` outputs and are gitignored — treat
  files there as regenerable artifacts, not source.
