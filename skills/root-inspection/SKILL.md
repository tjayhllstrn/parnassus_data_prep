---
name: root-inspection
description: Inspect, slice and plot segments of large ROOT TTrees from this CLAS12/Parnassus repository. Use whenever the user wants to look at a subset of events in a prepped ROOT file — e.g. "look at all events within the cut eventNumber > 3392000 && eventNumber < 3395200 and run the plotting macro on them", "plot pflow_pt vs truth_pt for these events", "dump the particles in event N", "what branches/trees are in this file", "how many events pass this cut". Handles evt_tree (Parnassus pflow/truth tree), reco_EventTree and mc_EventTree.
---

# ROOT inspection: slicing and plotting segments of a prepped TTree

Prepped files in `prepped_data/` are 100 MB–GB scale ROOT files with 10^5–10^7 entries.
This skill exists so a *segment* of such a file (defined by a TTree selection) can be
looked at, skimmed to a small file, and plotted — without ever mutating the input.

## Environment (do this first, every time)

ROOT is **not** on the bare cluster PATH. It comes from the parnassus-core conda env:

```bash
source ../parnassus-core/setup.sh      # module load conda; conda activate parnassus-core_env
```

Do **not** try the apptainer container mentioned in `CLAUDE.md` — it does not exist on this
cluster (NERSC). The wrapper script below sources the setup automatically, so prefer it over
raw `root` calls. Override the setup path with `PARNASSUS_CORE_SETUP=/path/to/setup.sh` if needed.

## The one entry point

```bash
skills/root-inspection/scripts/root_inspect.sh <subcommand> [options]
```

| subcommand | what it does |
| --- | --- |
| `info FILE [--values]` | trees, **cycles**, entry counts, branches (grouped by leaf counter), scalar min/max |
| `count FILE --cut C` | how many entries pass a selection (cheap; do this before skimming) |
| `select FILE --cut C` | write a *new* small ROOT file with only the passing entries |
| `dump FILE --cut C` | print the selected events: event scalars + a per-particle table |
| `plot FILE --expr E [--cut C]` | one histogram (1D, `y:x` 2D colz, or `--compare` pflow-vs-truth overlay) |
| `plot-all FILE --cut C` | skim by the cut, then run `macros/plot_prepped_data.C` on the skim |

All subcommands take `--tree NAME` (default `evt_tree`) and `--cycle N`.
Every subcommand prints exactly which files it wrote.

## Standard workflows

### "Look at all events in this cut and run the plotting macro on them"

```bash
S=skills/root-inspection/scripts/root_inspect.sh
F=prepped_data/merged_evt_tree.root
CUT="eventNumber > 3392000 && eventNumber < 3395200"

$S count    $F --cut "$CUT"                          # 1. how big is the segment?
$S plot-all $F --cut "$CUT" --label seg_3392k        # 2. skim + full plot set
```

`plot-all` writes the skim to `prepped_data/skims/<base>__seg_3392k.root` (kept by default;
`--no-keep-skim` deletes it after plotting) and the plots to `plots/<skim base>/`, because
`plot_prepped_data.C` derives its output directory from the input filename. Report the printed
paths back to the user.

### "Just one plot within that selection"

```bash
$S plot $F --expr pt --compare --cut "$CUT" --label seg_3392k     # pflow_pt vs truth_pt overlay
$S plot $F --expr pflow_pt --cut "$CUT" --bins 60 --xmin 0 --xmax 3 --label seg_3392k
$S plot $F --expr "pflow_eta:pflow_phi" --cut "$CUT" --bins 40 --logz --label seg_3392k
```

- `--compare` takes the bare stem (`pt`, `eta`, `vz`, `beta`, `charge`) and overlays
  `pflow_<stem>` against `truth_<stem>` in the repo's own style (`src/style.C`, same blue/magenta
  as `plot_prepped_data.C`), normalized to unit area.
- Single histograms are raw counts by default; add `--norm` to normalize, `--no-norm` to force counts.
- Plots go to `plots/<file base>/<label or cut tag>/<expr>.png`. `--save-root` also writes the
  histograms to a sidecar `.root`. `--out` sets an explicit path (must stay under `plots/`).
- Output stats (entries, mean, RMS, under/overflow) are printed — quote them, they are often the
  actual answer.

### "Show me the actual events / particles"

```bash
$S dump $F --cut "eventNumber == 3392001" --max-events 5 --max-particles 30
$S dump $F --cut "$CUT" --vars eventNumber,npflow,ntruth --no-particles   # event-level only
```

`dump` prints one block per entry with every `pflow_*` and `truth_*` row, and annotates
`*_class` with its meaning (0 charged hadron, 1 electron, 2 muon, 3 neutral hadron, 4 photon,
-1 unknown).

### "Make me a small file I can work with"

```bash
$S select $F --cut "$CUT" --label seg_3392k --dry-run   # count only, writes nothing
$S select $F --cut "$CUT" --label seg_3392k             # writes prepped_data/skims/...
```

The skim carries a `skim_info` TNamed recording the source file, tree, cycle and cut, so it is
never a mystery file later. `--also-trees reco_EventTree,mc_EventTree` skims extra trees with the
same cut, but only if the cut's branches exist there (see the `evnum` note below).

## Cut and expression syntax

Cuts are ordinary `TTree::Draw` selection strings, validated before use:

- Every identifier must be a real branch of the tree or a whitelisted math symbol; typos are
  rejected with a suggestion (`eventnumber` → "Did you mean: eventNumber?"). ROOT would otherwise
  silently treat an unknown name as 0 and return an empty or wrong selection.
- Shell/CINT injection is impossible: strings never touch a shell, and `;`, quotes, `>>`,
  `gSystem`, `ProcessLine`, backticks etc. are refused.
- `Sum$`, `Max$`, `Min$`, `Length$`, `Entry$`, `Alt$`, `abs`, `sqrt`, `TMath::*` … are available.

**Per-particle cuts:** `pflow_pt > 1` selects an entry if *any* element passes, and then the
whole event is copied/plotted. For explicit event-level logic write it out:
`Sum$(pflow_class==1) > 0` (event has a reco electron), `Max$(pflow_pt) > 2`,
`Sum$(pflow_pdgId==22) >= 2` (at least two reco photons).

## Facts about these files that will bite you

- **`eventNumber` is only unique in merged files.** In a single per-hipo-file prepped ROOT file it
  restarts (e.g. 1…10000 repeated across 200k entries), so `eventNumber == 3393` can match 20
  entries. `scripts/merge_evt_tree.sh` reassigns it sequentially at merge time
  (`macros/reassign_event_number.C`), which is where large values like 3392000 come from. Run
  `info --values` to see the real range before trusting a cut, and say so if the cut looks like
  it targets a merged file but the file is a single-run one.
- **Duplicate tree cycles.** Some prepped files contain `evt_tree;1` *and* `evt_tree;2` (repeated
  `Write()` calls) with different entry counts. Bare `Get("evt_tree")` silently takes the highest
  cycle. `info` lists all cycles, every subcommand warns when there is more than one, and
  `--cycle N` pins a specific one. Mention which cycle was used.
- **Older files lack `*_beta` and `*_charge`.** `plot_prepped_data.C` requires them, so `plot-all`
  checks first: if they are missing it warns and produces an equivalent built-in comparison set
  plus printed `class`/`pdgId` composition tables instead. `--strict` makes that a hard failure.
- **`evt_tree` vs the EventTrees.** `evt_tree` (ParnassusTree) has `eventNumber`, `npflow`/`ntruth`
  and the `pflow_*`/`truth_*` arrays. `reco_EventTree`/`mc_EventTree` are a different schema
  (~105 branches) and use **`evnum`**, not `eventNumber` — translate cuts accordingly and pass
  `--tree reco_EventTree`.
- `prepped_data/` and `plots/` are gitignored; skims and plots are regenerable artifacts.

## Safety rules this skill enforces (don't work around them)

- Inputs are opened `READ` only. Nothing is ever modified in place — `select` always writes a new
  file, and no subcommand deletes anything except a skim it just created under `--no-keep-skim`.
- Outputs must live under `prepped_data/skims/` (ROOT files) or `plots/` (images), must have the
  right extension, may never equal an input, and never overwrite an existing file without `--force`.
- `select`/`plot-all` refuse to copy more than `--max-entries` (default 500,000) entries. If a cut
  selects more, tighten the cut rather than raising the cap by reflex.
- A cut selecting 0 entries is an error, not a silently empty plot.

Prefer these scripts over ad-hoc `root -l -b -q` one-liners: a hand-written macro that opens a
300 MB file with the wrong tree cycle, a typo'd branch name, or `"UPDATE"` mode is exactly the
failure mode this skill removes. If something genuinely isn't covered, add a subcommand to
`scripts/root_tools.py` rather than bypassing it.
