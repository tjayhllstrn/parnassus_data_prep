# Parnassus CLAS12 v7 dataset: build task (p instead of pt, unfolded truth theta)

## What this is

`parnassus_data_prep` turns CLAS12 RG-A clasdis Monte Carlo `.hipo` files into flat ROOT trees.
For every processed file it writes `reco_EventTree`, `mc_EventTree` and the ML-facing `evt_tree`
(see `CLAUDE.md` and `parnassus_data_structure.md`). The merged `evt_tree` files are the training
and test data for Parnassus, a flow-matching fast simulation trained on NERSC Perlmutter. The
current dataset is v6.

The task is to produce **v7 = v6 with exactly two changes**, built on ifarm from the processed
per-hipo ROOT files:

1. Store the momentum magnitude `p` in place of the transverse momentum `pt`, for both collections.
2. Store the true polar angle of every truth particle in [0, pi]. Do not fold backward-going
   particles into the forward hemisphere.

The event set, the event order, the particle lists and every other branch must be identical to v6.

## Why

- **Momentum variable.** The particle model reproduces the pt and theta residuals of full
  simulation per species, but not the |p| residuals. CLAS12 tracking measures p and theta with
  nearly independent errors: |corr(dp/p, dtheta)| <= 0.15 for e, pi, K and p. In the (pt, theta)
  basis the same errors are strongly correlated, corr(dpt/pt, cot(theta) dtheta) = 0.81 to 0.88
  for charged hadrons. The model only partly learns that correlation, so |p| comes out 1.1x to
  1.6x too wide. Training directly on p removes the problem.
- **Folded truth theta.** `Kinematics::th` in `src/Kinematics.C` returns `abs(atan(Pt/Pz))`, so a
  particle with pz < 0 is stored at pi - theta_true. A particle at 160 deg is stored at 20 deg,
  inside the Forward Detector acceptance, where the model cannot tell it from a real forward
  particle. In `prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.root` (100k events), 4.12% of
  truth particles have pz < 0, and 25.6% of events contain at least one.

## Scope rules

- **Do not change any selection from the v6 selection** That covers event selection, particle selection, the scatE
  filtering, the electron beta cut and the ordering. v7 inherits all of them from v6 unchanged.
- Do not rerun `hipo2tree_parnassus.C` unless the processed files are gone (see Fallback).

## Target: what the v7 files must look like

| v7 file | Built from | Entries |
|---|---|---|
| `v7_train.root` | `v5/train_files` | 50,000,000 |
| `v7_test.root` | `v5/test_files` | 5,000,000 |
| `README.md` | (new) | n/a |

**Out of scope:** `v7_test_filtered.root` and `v7_test_filtered_truth.root`. They are produced
later on Perlmutter from `v7_test.root` by the parnassus-core script `scatE_filter.py`.

**Layout:** the tree is named `evt_tree`. Each v7 file has the same entries, in the same order, with
the same `eventNumber` per entry as its v6 counterpart.

| Branch (leaflist) | v7 content |
|---|---|
| `eventNumber/I`, `ntruth/I`, `npflow/I`, `pflow_scatE/I`, `truth_scatE/I` | identical to v6 |
| `pflow_p[npflow]/D` | replaces `pflow_pt`; equals `pflow_pt / sin(pflow_theta)`, in GeV |
| `truth_p[ntruth]/D` | replaces `truth_pt`; equals `truth_pt / sin(truth_theta_v6)`, in GeV. This is exact even for folded particles, because sin(pi - x) = sin(x) |
| `truth_theta[ntruth]/D` | unfolded: `atan2(pt, pz)` in [0, pi]. Equivalently `theta_v6` if `truepz >= 0`, else `pi - theta_v6` |
| `pflow_theta[npflow]/D` | identical to v6; every reco particle has pz > 0 (theta 5 to 35 deg), so it already equals `atan2(pt, pz)` |
| `*_phi`, `*_vx`, `*_vy`, `*_vz`, `*_beta`, `*_pdgId` | identical to v6 |

- **No `*_pt` branches in v7.** Downstream code will derive pt as p sin(theta).
- **Branch order:** consumers read by name, so order does not matter. Keep v6's order, with each
  `*_pt` replaced in place by `*_p`.
- **Leaflists** follow `ParnassusTree` in `src/TreeManager.C`: `/D` and `/I`, with `npflow` and
  `ntruth` as the counters.
- **Compression:** ZLIB level 1 (ROOT setting 101), as in v6.

## Inputs on ifarm

- **Processed per-hipo files**, which contain `evt_tree`, `reco_EventTree` and `mc_EventTree`:
  - `/volatile/clas12/users/tjhellst/parnassus_data/v5/test_files/`: 16 files.
  - `/volatile/clas12/users/tjhellst/parnassus_data/v5/train_files/`: 148 files. Only the first
    two, in `LC_ALL=C sort` order, fed v6.
  - Beta branches were added to these files in place by `macros/add_beta_branches.C`.
  - `/volatile` is auto-cleaned, so confirm the files still exist before planning around them.
- **v6 files:** look for them under `/volatile/clas12/users/tjhellst/parnassus_data/v6/`. 
  (`/global/cfs/cdirs/m3246/tjhellst/parnassus_data/v6/`: `v6_train.root` is 17 GB, `v6_test.root`
  is 1.6 GB).
- **Fallback, only if the processed files are gone:** regenerate them with the unchanged prep code
  via `scripts/submit_hipo2tree_parnassus.sbatch`, from
  `/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/`.

## The only new information needed: the sign of truepz

Everything in v7 except the theta unfolding can be computed from v6 alone. The sign of pz for each
truth particle is not stored in v6. It must come from `mc_EventTree.truepz` in the processed files.

The facts below were verified on `prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.root`
(cycles `evt_tree;2`, `mc_EventTree;3`, `reco_EventTree;2`; 100k entries each):

- **Alignment.** The three trees line up entry by entry and particle index by particle index:
  - `mc_EventTree.Nmax == evt_tree.ntruth` and `reco_EventTree.Nmax == evt_tree.npflow` for every
    entry.
  - `truepid == truth_pdgId`, `hypot(truepx, truepy) == truth_pt` and `pid == pflow_pdgId` for
    every particle.
  - This is the same premise `macros/add_beta_branches.C` relies on.
- **p and theta to machine precision.** `pt / sin(theta)` matches the stored `truep` (truth) and
  `p` (reco) to a relative difference of at most 5e-16. The unfolded theta,
  `atan2(hypot(truepx, truepy), truepz)`, matches `theta_v6` or `pi - theta_v6` to at most 5e-16.
- **pz exactly zero.** It occurs for 73 particles in 100k events, mostly photons, all with
  pt > 0. Theta is pi/2 either way.
- **No NaN.** There is no NaN truth_theta in the first 5M v6_train entries.

## Mapping v6 events to processed-file entries

v6 `eventNumber` is not the hipo event number; it was renumbered with
`macros/reassign_event_number.C`. The evidence:

- **v6_train:** `eventNumber` runs 0 to 41,767,869 with gaps, over 41,738,530 entries. That fits a
  renumbering followed by `apply_lepton_beta_cut.C` with `dropNoElectron = true`, which dropped
  29,340 events. Every event has `pflow_scatE == 1`.
- **v6_test:** `eventNumber` runs 200,000,000 to 204,999,999 with no gaps, i.e. a 5M slice taken
  with `copy_slice.C` from a renumbered merged file. Only 18.8% of events have `pflow_scatE == 1`
  (unfiltered). The lepton beta cut was applied without dropping events.
- `parnassus_data/v6/README.md` (on Perlmutter) adds: "All of this data came from the test_files
  folder of v5 and the first two files of the train_files folder of v5. The test dataset was taken
  from entry 200M onward from the resulting file. Then the training data came from filtering out
  the non-electron reconstructed events from the rest of this data and appending the two files from
  train_files to make up for more statistics."

The exact v6 recipe is not recorded anywhere, so **do not rely on reproducing it.** Use a content
join instead:

1. **Key.** Build a per-event key that `evt_tree` processing never changes: a 64-bit hash over the
   raw bytes of `(ntruth, truth_pt[], truth_pdgId[])`. Every v6 post-processing step (merging, the
   scatE filter, the electron beta cut, renumbering, slicing) touches only pflow arrays or
   `eventNumber`. The truth arrays in v6 are therefore bit-identical to those in the processed
   `evt_tree`.
   - On the sample, the key was unique across all 100k events.
   - Check for duplicate keys across the full input set. If any occur, extend the key with
     `truth_phi[]` and `truth_vz[]`.
2. **Sign table.** Build a map from key to the per-particle sign pattern of `truepz`, taken from
   `mc_EventTree` in the 16 test files and the first 2 train files.
   - For each entry, assert the alignment checks above before using it.
   - Store `truepz < 0` per particle. A bitmask works: `ntruth <= 25` in v6, but processed files
     allow up to 100, so size the mask accordingly.
3. **Write v7.** For every v6 entry, in order: look up the key, flip theta where `truepz < 0`,
   compute `p = pt / sin(theta_v6)` for both collections, and write the v7 entry.
   - Every v6 entry must match exactly once. Report the number of unmatched entries, which must be
     0. If any are unmatched, extend the sign table to more `train_files` before going further.
   - About 250M input events at 12 bytes per entry is roughly 3 GB, which fits on one node.
   - A compiled ROOT macro (ACLiC, `+`) modeled on `add_beta_branches.C` is the most portable
     choice inside the CLAS12 container. PyROOT or RDataFrame also work if available.
4. **Optional shortcut.** If a `TChain` over the sorted `test_files` at global entry
   200,000,000 + k reproduces v6_test entry k exactly (same truth arrays), a direct index mapping
   can build `v7_test.root`. Run the full validation either way.

## Validation (acceptance criteria)

Run these on the complete v7 files and report the numbers.

1. **Entries.** v7 entry counts equal v6 (41,738,530 and 5,000,000), and `eventNumber` is
   identical entry by entry.
2. **Unchanged branches.** Every branch except `*_pt`/`*_p` and `truth_theta` is bit-identical to
   v6, entry by entry.
3. **Momentum.** `|p * sin(theta_v7) - pt_v6| / pt_v6 < 1e-12` for every particle in both
   collections.
4. **Angles.**
   - `truth_theta_v7` is `theta_v6` or `pi - theta_v6` for every particle.
   - `pflow_theta` is unchanged; in v6 it ranges from 0.08727 to 0.61086 rad.
   - `truth_theta` lies within (1 deg, 179 deg), and its maximum exceeds pi/2. In v6 it is exactly
     pi/2.
5. **Backward fraction.** Report the fraction of truth particles with theta > pi/2 (expect about
   4%; 4.12% in the -0001 sample) and the fraction of events with at least one (expect about 26%;
   25.6% in the sample).
6. **Physics closure.** Use E = p / beta (stored `truth_beta` is the generator's p/E) and
   pz = p cos(theta). The sum over truth particles of (E - pz) should equal the proton mass,
   0.938272 GeV, within 3e-4 GeV in about 83% of events.
   - In the sample, with the true pz signs, the rate is 82.6%.
   - With the v6 folded theta the same test passes in only 61.9% of events (50k v6_test events).
   - 100% is not expected. Even with true signs about 17% of events fail, because the stored
     truth list is incomplete: particles are dropped by the 1 deg theta cut, and neutrinos are
     removed.
7. **Finite values.** No NaN or inf appears in any `*_p` or `*_theta`.
8. **Spot check.** For 1000 random v7 entries, go back to the matched `mc_EventTree` entry and
   confirm `truth_p == truep` and `truth_theta == atan2(hypot(truepx, truepy), truepz)` to 1e-12.

## Deliverables

- **Data on ifarm:** `/volatile/clas12/users/tjhellst/parnassus_data/v7/` containing
  `v7_train.root`, `v7_test.root` and `README.md`. The final destination is
  `/global/cfs/cdirs/m3246/tjhellst/parnassus_data/v7/` on Perlmutter; the transfer is a separate
  step.
- **`README.md` in the v7 directory**, covering:
  - that v7 is v6 with identical events and order;
  - the two changes;
  - the macro and exact commands used;
  - the validation numbers.
- **The new macro or script** in `macros/`, with a usage header in the style of
  `add_beta_branches.C`. Existing files stay untouched.
- **A short report:** commands run, input files used, entry counts, unmatched or duplicate keys,
  and the output of every validation item.

## Findings to carry forward

- **Truth particles beyond 179 deg are lost.** `CLAS12Ana::load_mc_particles` applies its 1 deg
  cut to the folded theta, so particles with true theta > 179 deg were dropped. They are absent
  from the processed files and cannot be recovered here; they stay absent in v7. They should be
  very rare. Recovering them needs the prep-code fix plus a hipo rerun, which is out of scope.
- **`mc_EventTree.truetheta` is folded too.** It goes through the same `Kinematics::th`, with a
  maximum of exactly 90 deg. Take the sign from `truepz`, not from `truetheta`.
- **Generator beam energy.** In complete events, the sum of truth E minus the proton mass comes
  out near 10.600 GeV, not 10.604. This is irrelevant to v7; it is noted only so it is not mistaken
  for a bug introduced by the conversion.

## Rejected alternatives (do not re-propose without new evidence)

- **Recovering pz signs from v6 alone through per-event conservation of (E - pz).** Tested on 50k
  v6_test events:
  - 61.9% need no flips.
  - 20.6% have a unique sign assignment.
  - 17.3% show a deficit that no assignment explains. Their median |sum pT| imbalance is 26 MeV,
    against 0.1 MeV in complete events, which points to particles missing from the truth list.
  
  It is not exact, so the processed files are used instead.
- **Inverting the fold from beta and mass.** These give |p| only, never the sign of pz.

## Deliberately left alone

- **The `Kinematics::th` fix in `src/` and a hipo rerun.** Separate change.
- **v7 filtered test files.** Produced on Perlmutter with `scatE_filter.py`.
- **parnassus-core changes needed to consume `*_p`.** Pt-based computed variables (ptrel, ht,
  Pxh, Pyh, DIS variables) will derive pt as p sin(theta) on the Perlmutter side.

## Status

| Item | Status to confirm |
|---|---|
| v7 data files | Not produced. Done when both files exist under the v7 directory and pass every validation item. |
| v7 build macro | Not written. Done when a new file in `macros/` builds v7 and `git status` shows no changes to existing files. |
| `Kinematics::th` fix | Out of scope. `src/Kinematics.C` still uses `abs(atan(Pt/Pz))`. |
