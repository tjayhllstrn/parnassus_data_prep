# The Parnassus data format

A guide to `parnassus_data_sample.root` for readers who know physics but have never
opened a ROOT file, and who are new to CLAS12.

---

## 1. What this file is, in one paragraph

CLAS12 is a large-acceptance spectrometer in Hall B at Jefferson Lab. A 10.6 GeV electron beam hits a
proton target; the scattered electron and the hadrons produced along with it fly out through a stack of
tracking chambers, Cherenkov counters, time-of-flight walls and calorimeters. Reconstruction software
turns the raw detector hits into a list of *candidate particles* per event — each with a guessed
identity, a momentum and a vertex. This file is a Monte Carlo (simulated) sample, so for every event we
also know the **truth**: the exact list of particles the generator actually produced, before the
detector had a chance to smear, lose or misidentify anything.

Each entry of this file is therefore **one event containing two lists of particles**:

| | what it is | how many, on average |
|---|---|---|
| `pflow_*` | what the detector + reconstruction *thinks* was there | 2.07 particles/event |
| `truth_*` | what was *actually* there | 7.15 particles/event |

The name "pflow" is short for *particle flow* — the reconstruction-level particle list. "Parnassus" is
the machine-learning pipeline these trees feed: the task is to learn the map from the messy `pflow`
list to the clean `truth` list. The large gap between 2.07 and 7.15 is the physics of the problem: most
produced particles are never reconstructed, because they miss the detector's angular acceptance, fall
below thresholds, or fail quality cuts.

**Provenance.** This sample was produced by `macros/hipo2tree_parnassus.C` from `clasdis` Monte Carlo
generated for the RG-A Fall 2018 inbending configuration (beam energy **10.6041 GeV**, fixed proton
target, torus magnet inbending). Several per-file outputs were then concatenated with
`macros/merge_evt_tree.C` and renumbered with `macros/reassign_event_number.C`.

**Size.** One tree, `evt_tree`, with **1,000,000 entries**, 19 branches, ~265 MiB on disk. It holds
2,072,701 reconstructed particles and 7,147,997 truth particles in total.

---

## 2. What a ROOT file and a TTree are

If you have used pandas, HDF5 or a FITS table, the mental model transfers directly.

A **ROOT file** (`.root`) is a container, a bit like a zip archive or an HDF5 file. It holds named
objects. Listing them:

```
$ root -l parnassus_data_sample.root
root [0] .ls
TFile**         parnassus_data_sample.root
 TFile*         parnassus_data_sample.root
  KEY: TTree    evt_tree;1      evt_tree
```

or, without ROOT at all:

```python
>>> import uproot; uproot.open("parnassus_data_sample.root").keys()
['evt_tree;1']
```

This file contains exactly one object: a **TTree** called `evt_tree`. (Sibling files produced directly
by the pipeline, e.g. `clasdis_..._small_version.root`, also contain `reco_EventTree` and
`mc_EventTree` — much wider trees with full detector-level detail. Those are *not* in this file.)

A TTree is a **columnar table**:

- one **entry** (row) per event — here, 1,000,000 of them;
- one **branch** (column) per quantity.

The `;1` after the name is a *cycle number*, ROOT's version counter for objects with the same name.
Ignore it; `f["evt_tree"]` or `f->Get("evt_tree")` always gives you the latest.

### Jagged branches and their counters

Most columns here are not scalars. Event 4 has 6 reconstructed particles and 13 truth particles;
event 5 has different numbers. So `pflow_pt` is not one number per event, it is a **variable-length
array** per event ("jagged" or "ragged" in the Python world).

ROOT stores the length in a separate integer branch, and the array branch declares its dependence on
it. In this file:

- every `pflow_*` array has length `npflow`;
- every `truth_*` array has length `ntruth`.

You will see this written as `pflow_pt[npflow]/D` — "an array of doubles, `npflow` of them". `npflow`
and `ntruth` are called *leaf counters*. Modern readers (uproot, RDataFrame, `TTree::Draw`) follow the
counter automatically and hand you a correctly-sized array; you only need to know about it if you write
raw C++ against the tree.

> The writer that produced this tree uses fixed C arrays of size 100 internally, so `npflow` and
> `ntruth` can never exceed 100. Observed maxima in this sample are 29 and 23.

### Type codes

ROOT spells types with one-letter codes. Only two appear here:

| code | C++ | meaning |
|---|---|---|
| `/I` | `int32_t` | 32-bit signed integer |
| `/D` | `double` | 64-bit float |

---

## 3. Conventions: units, angles, coordinates

Fixed for every branch in the file.

**Coordinate system.** The CLAS12 lab frame. The *z* axis points along the beam, downstream. The origin
is the nominal centre of the CLAS12 detector, *not* the target — for this run configuration the target
sits near *z* ≈ −3 cm. *x* and *y* are transverse to the beam.

**Units.**

| quantity | unit |
|---|---|
| momentum (`*_pt`) | GeV/*c*, written GeV throughout |
| angles (`*_theta`, `*_phi`) | **radians** |
| vertex positions (`*_vx`, `*_vy`, `*_vz`) | **centimetres** |

Radians and centimetres are both easy to get wrong: the CLAS12 literature quotes angles in degrees and
some CLAS12 banks use millimetres. This file does neither. Convert with `theta*180/pi` when comparing
to a paper.

**Angle definitions.** `theta` is the polar angle from the beam axis: θ = 0 is straight downstream,
θ = π/2 is transverse. `phi` is the azimuth about the beam axis, `atan2(py, px)`, so it runs on
(−π, +π]. `pt` is the momentum component transverse to the beam, √(p<sub>x</sub>² + p<sub>y</sub>²).

Note that θ and φ and p<sub>T</sub> specify a momentum vector only up to the *sign* of
p<sub>z</sub>. See the folding caveat in §6.

**Particle identity.** `*_pdgId` uses the [PDG Monte Carlo numbering
scheme](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-monte-carlo-numbering.pdf). Sign encodes
particle vs. antiparticle. The codes that actually occur in this file are tabulated in §5.

**Missing / sentinel values.** There are none in this file — every array slot up to the counter holds a
real particle. (Other trees in this codebase use `-999` as a sentinel; `evt_tree` does not.)

---

## 4. Reading the file

### What you need installed

Nothing CLAS12-specific. This is a plain ROOT file: any of the following is enough.

*note: these commands are installing software on the cluster, so should be adapted to whichever environment you are working in. 
pip install should usually be done in a python virtual environment and conda install is meant for instead a conda environemnt. 
Both of these environment tools are ways to control where your software is contained and to keep your space from getting too cluttered.

```bash
pip install uproot awkward numpy      # pure Python, no ROOT needed — easiest
conda install -c conda-forge root     # full ROOT, gives you the C++ prompt and PyROOT
```

Worth stressing, because it surprises people: *producing* these files needs the `clas12root`
environment (a container with the CLAS12 and `hipo` libraries), since the pipeline reads raw `.hipo`
input. *Reading* a finished `.root` file needs none of that. If a `parnassus_data_sample.root` has been
handed to you, you can analyse it anywhere.

In the examples below, substitute your own path for `parnassus_data_sample.root`.

### Python (uproot) — recommended

`uproot` reads ROOT files with no ROOT installation at all, returning numpy/awkward arrays.

```python
import uproot, awkward as ak, numpy as np

f = uproot.open("parnassus_data_sample.root")
t = f["evt_tree"]

t.num_entries          # 1000000
t.keys()               # list of branch names
t.show()               # branch names, types, and their leaf counters

# read three branches for the first 1000 events
a = t.arrays(["npflow", "pflow_pt", "pflow_pdgId"], entry_stop=1000)

a["pflow_pt"][4]          # array of pt for the 6 particles in event 4
ak.num(a["pflow_pt"])     # per-event multiplicity; equals a["npflow"]
ak.flatten(a["pflow_pt"]) # all particles from all events in one flat array
```

`ak.flatten` is the workhorse: single-particle distributions (a p<sub>T</sub> spectrum, a θ
distribution) are almost always made from flattened arrays, while per-event quantities
(multiplicities, sums) keep the jagged structure.



### C++ / interactive ROOT

```
$ root -l parnassus_data_sample.root
root [0] evt_tree->Print()                                        // full branch listing
root [1] evt_tree->Scan("npflow:ntruth:pflow_pdgId", "", "", 5)   // first 5 events
root [2] evt_tree->Draw("pflow_pt", "", "hist")                   // flattens automatically
root [3] evt_tree->Draw("truth_theta*TMath::RadToDeg()", "truth_pdgId==22")
```

Opening the file on the command line puts the tree in scope as a global named `evt_tree`, which is why
the snippets above can use the bare name.

`TTree::Draw` is the quickest way to look at anything: the first argument is an expression over branch
names, the second a per-particle or per-event selection.

### RDataFrame (C++ or PyROOT)

```python
import ROOT
df = ROOT.RDataFrame("evt_tree", "parnassus_data_sample.root")
df.Count().GetValue()
h = df.Filter("pflow_scatE == 1").Histo1D("pflow_pt")
```

---

## 5. Branch reference

19 branches. Five are per-event scalars; the other fourteen are per-particle arrays — seven describing
reconstructed particles and the same seven repeated for truth particles.

All statistics below were measured on this file. Counts, fractions, minima and maxima are over all
1,000,000 entries; percentiles (p1/p50/p99) are over the first 200,000 entries.

### 5.1 Event-level branches

#### `eventNumber` — `int32_t`

A unique sequential index for the entry, `0 … 999999`, strictly increasing.
Treat it as a row label: useful for referring to a specific event, meaningless as physics.

#### `npflow` — `int32_t`

Number of reconstructed particles in this event, and the length of every `pflow_*` array.


#### `ntruth` — `int32_t`

Number of truth (generator-level) particles in this event, and the length of every `truth_*` array.

Range 2–23, mean 7.15. Never zero — every event has at least the scattered electron plus something
else. Peaks around 6–7.

#### `pflow_scatE` — `int32_t`

Flag, 0 or 1: does the reconstructed particle list contain the identified **scattered electron**?

Measured: **`pflow_scatE == 1` in 18.7% of events.** So four events in five have no usable
reconstructed electron, even though (see below) every event contains one in truth.


#### `truth_scatE` — `int32_t`

The same flag for the truth list: is one of the truth particles the scattered electron?

**Measured value: 1 for all 1,000,000 entries.** The `clasdis` generator produces DIS events by
construction, so a scattered electron is always present. The branch is constant in this file and
carries no information; it exists for symmetry with `pflow_scatE`, and would be informative for a real
data or non-DIS sample.


### 5.2 Reconstructed-particle branches (`pflow_*`, length `npflow`)


#### `pflow_pdgId[npflow]` — `int32_t`

Reconstructed particle identity, as a PDG code. This is the detector's **hypothesis**, not the truth —
misidentification is exactly one of the effects the Parnassus model has to learn.

Eleven distinct values occur:

| code | particle | count | fraction |
|---|---|---:|---:|
| `22` | photon (γ) | 881,485 | 42.53% |
| `211` | π⁺ | 362,391 | 17.48% |
| `2112` | neutron | 307,255 | 14.82% |
| `11` | electron (e⁻) | 189,565 | 9.15% |
| `2212` | proton | 148,788 | 7.18% |
| `-211` | π⁻ | 105,271 | 5.08% |
| `-11` | positron (e⁺) | 22,561 | 1.09% |
| `321` | K⁺ | 22,239 | 1.07% |
| `-321` | K⁻ | 13,650 | 0.66% |
| `45` | deuteron | 12,540 | 0.61% |
| `-2212` | antiproton | 6,956 | 0.35% |



#### `pflow_pt[npflow]` — `double`, GeV

Transverse momentum, √(p<sub>x</sub>² + p<sub>y</sub>²).

Range 0.0007 – 2.4974 GeV; p1 = 0.022, median = 0.332, p99 = 1.362.


Note that p<sub>T</sub> alone does not give you the momentum magnitude; combine with `pflow_theta`:
p = p<sub>T</sub> / sin θ.

#### `pflow_theta[npflow]` — `double`, radians

Polar angle from the beam axis.

Range **exactly [0.087266, 0.610865] rad = [5.000°, 35.000°]**, median 0.320 rad (18.3°).


#### `pflow_phi[npflow]` — `double`, radians

Azimuthal angle about the beam axis, `atan2(py, px)`, on (−π, +π].

Range −3.1416 – 3.1416, essentially flat, median ≈ 0. The φ distribution is *not* perfectly uniform:
the Forward Detector is built from six sectors with dead material between them, and the fiducial cuts
carve out the sector edges, so a fine-binned φ histogram shows six-fold structure.

#### `pflow_vx[npflow]`, `pflow_vy[npflow]` — `double`, centimetres

Transverse position of the particle's reconstructed origin (its vertex).

Both are sharply peaked on the beam axis with long tails. For `pflow_vx`: median 0.000, interquartile
range ±0.23, p5/p95 = ±1.7, p1/p99 = ±5.7, full range −57.1 to +85.3. `pflow_vy` is statistically
identical. 70% of particles have |v| < 0.5 cm.

The peak is the beam spot; the tails are reconstruction artefacts. Neutrals have no track, so their
"vertex" is assigned rather than measured, and mis-tracked charged particles can be extrapolated far
off axis. Charged particles alone show the same core width, so the tails are not purely a neutral
effect.

#### `pflow_vz[npflow]` — `double`, centimetres

Longitudinal position of the reconstructed vertex, along the beam.

Median −3.05, interquartile range [−4.53, −1.55], p5/p95 = [−7.33, +1.08], p1/p99 = [−17.8, +10.7],
full range −166.5 to +172.4.

The median at **−3 cm** is the physical target position for this run configuration. The interquartile
spread of ~3 cm reflects the target length plus vertex resolution; the extreme values are again
artefacts of assigned neutral vertices and mis-tracking. Note this is the one coordinate whose
distribution is *not* centred on zero — see the truth/reco offset noted in §5.3.

In about 40% of multi-particle events, *all* particles share bit-identical `vx`, `vy` and `vz` — these
are events where the Event Builder assigned every particle the same common event vertex. Treat the
per-particle vertex as weakly informative for neutrals.

### 5.3 Truth-particle branches (`truth_*`, length `ntruth`)

#### `truth_pdgId[ntruth]` — `int32_t`

True particle identity, a genuine PDG code (no CLAS12 conventions here). Fourteen distinct values:

| code | particle | count | fraction |
|---|---|---:|---:|
| `22` | photon (γ) | 2,612,263 | 36.55% |
| `211` | π⁺ | 1,330,874 | 18.62% |
| `11` | electron (e⁻) | 1,015,739 | 14.21% |
| `-211` | π⁻ | 964,949 | 13.50% |
| `2212` | proton | 568,858 | 7.96% |
| `2112` | neutron | 439,467 | 6.15% |
| `321` | K⁺ | 101,685 | 1.42% |
| `130` | K<sub>L</sub>⁰ | 54,738 | 0.77% |
| `-321` | K⁻ | 32,865 | 0.46% |
| `-11` | positron (e⁺) | 15,640 | 0.22% |
| `-2212` | antiproton | 6,129 | 0.09% |
| `-2112` | antineutron | 4,706 | 0.07% |
| `13` | μ⁻ | 42 | 0.0006% |
| `-13` | μ⁺ | 42 | 0.0006% |

The photon fraction is dominated by π⁰ → γγ: the π⁰ itself is not a final-state particle, so it never
appears — only its two daughter photons do. Electrons are 14.2% of truth particles: one per event is the
scattered electron, the rest come from conversions and Dalitz decays.

The truth composition is close to charge-symmetric (π⁺ 18.6% vs. π⁻ 13.5%, the excess coming from the
proton target's valence content), whereas the reconstructed composition is strongly asymmetric
(π⁺ 17.5% vs. π⁻ 5.1%). That asymmetry is the inbending torus field, which focuses negatives out of and
positives into the Forward Detector acceptance — a detector effect, not physics.

#### `truth_pt[ntruth]` — `double`, GeV

True transverse momentum. Range 0.0001 – 2.0809 GeV; p1 = 0.015, median = 0.285, p99 = 1.290.

Unlike `pflow_pt` there is **no p<sub>T</sub> cut here**; the 2.08 GeV maximum is simply where the
statistics run out. If you compare the two spectra, remember the reco side is truncated at 2.5 GeV and
the truth side is not.

#### `truth_theta[ntruth]` — `double`, radians

True polar angle — **but folded into the forward hemisphere. Read §6 before using this branch.**

Range [0.017454, 1.570796] rad = [1.000°, **90.000°**]. Median 0.363 rad (20.8°). The lower edge is the
1° load-time filter. The upper edge of exactly π/2 is the folding artefact: every backward-going
particle is mapped into (0, π/2], so the region just below 90° holds a mixture of genuinely
near-transverse particles and mirrored backward ones.

#### `truth_phi[ntruth]` — `double`, radians

True azimuth, `atan2(py, px)`, on (−π, +π]. Flat, median ≈ 0. Unlike `pflow_phi` there is no sector
structure — no detector geometry has been imposed. Unlike `truth_theta`, this branch is *not* folded.

#### `truth_vx[ntruth]`, `truth_vy[ntruth]`, `truth_vz[ntruth]` — `double`, centimetres

True production vertex of the particle.

For **primary** particles — those produced at the interaction point — the vertex is the beam spot:
transverse spread of order ±0.01 cm and `vz` exactly 0. The scattered electron has `truth_vz == 0` in
99.9% of events, and K<sub>L</sub>⁰ in 100%.

For **secondary** particles — decay products and interaction products — the vertex is displaced
downstream, sometimes far downstream. Measured `truth_vz` percentiles by species:

| species | p50 | p99 | max | fraction with `vz == 0` |
|---|---:|---:|---:|---:|
| e⁻ | 0.00 | 0.00 | 66.6 | 99.9% |
| π⁺ | 0.00 | 8.36 | 152.9 | 96.8% |
| γ | 0.00 | 18.2 | 257.3 | 94.1% |
| n | 0.00 | 28.1 | 257.3 | 89.8% |
| p | 0.00 | 32.0 | 262.7 | 87.6% |

Full range −2.15 to +280.4 cm. A vertex at *z* = 250 cm is a hadronic interaction or decay deep inside
the detector volume, not at the target. So `truth_vz != 0` is a practical (if imperfect) tag for "this
particle is a secondary".

**A systematic offset to be aware of:** the truth interaction point is at *z* = 0, while the
reconstructed vertex distribution is centred at *z* ≈ −3 cm (the physical target position). Truth and
reco vertices are therefore displaced by ~3 cm relative to each other in this sample. Do not treat
`truth_vz` and `pflow_vz` as living on a common zero.

---

## 6. Caveats

Ranked roughly by how likely they are to bite you.

### `pflow` and `truth` are two unordered sets, not aligned lists

There is **no index correspondence**. `pflow_pt[3]` and `truth_pt[3]` describe unrelated particles.
The two lists have different lengths in 98.4% of events (`npflow < ntruth` in 96.5%, `npflow > ntruth`
in 1.9%), and in the 1.6% where the lengths coincide, position *i* in one list still has nothing to do
with position *i* in the other.

This is by design — it is the whole point of a particle-flow learning task. The model must map an
unordered set of detector-level objects onto an unordered set of true particles, and the loss has to be
permutation-invariant (or use an explicit assignment such as Hungarian matching). If you want a
per-particle truth label attached to a reco particle, `evt_tree` cannot give it to you.

The ordering *within* each list is also worth knowing:

- `truth` is generator order, and the scattered electron is always element 0.
- `pflow` is Event Builder order with cut-failing entries deleted. `pflow_pdgId[0] == 11` in only 29% of
  non-empty events, because the electron is so often removed. **Never assume `pflow[0]` is the
  electron** — use `pflow_scatE` and search for `pflow_pdgId == 11`.



### `truth_scatE` is constant

It is 1 for every entry. Do not use it as a feature or a selection; it will silently do nothing.

### Cuts have already been applied, and they are not reversible

The reconstructed list is post-cut. `CutManager::filter_particles` applied, per particle:

- **Forward Detector acceptance** (`is_in_FD`): 5° ≤ θ ≤ 35°, applied to every particle.
- **p<sub>T</sub> ceiling** (`maxPtCut`): p<sub>T</sub> ≤ 2.5 GeV, applied to every particle.
- **Per-species cuts**, keyed off the reconstructed PDG code: drift-chamber fiducial cuts, calorimeter
  sampling-fraction and minimum-energy cuts, calorimeter edge cuts, χ²<sub>PID</sub> and vertex-*z* cuts
  for electrons; pion and photon cuts referencing the scattered electron; a minimum-energy cut for
  neutrons. Particles with `pdgId == 0` are always rejected.
- Species not named explicitly (kaons, protons, deuterons) pass with only the θ and p<sub>T</sub> cuts.


### No DIS kinematic variables

`Q2`, `x`, `W` and ν are *computed* by the pipeline but they are absent here. For a `pflow_scatE == 1`
event you can recompute them yourself from the scattered electron: with beam energy E = 10.6041 GeV,
electron energy E′ ≈ p<sub>T</sub>/sin θ (the electron mass is negligible), and proton mass
M<sub>p</sub>:

```
ν  = E − E′
y  = ν / E
Q² = 4 E E′ sin²(θ/2)
W  = sqrt(M_p² + 2 M_p ν − Q²)
x  = Q² / (2 M_p ν)
```


---

## 7. Recipes

### Reconstructed p<sub>T</sub> spectrum for identified pions

```python
import uproot, awkward as ak
t = uproot.open("parnassus_data_sample.root")["evt_tree"]
a = t.arrays(["pflow_pt", "pflow_pdgId"], entry_stop=200000)
pions = a["pflow_pt"][abs(a["pflow_pdgId"]) == 211]
pt = ak.to_numpy(ak.flatten(pions))
```

### Reconstruction efficiency vs. truth θ, inside the Forward Detector window

Because the two lists are unmatched, "efficiency" here means a ratio of *counts*, not a per-particle
success rate:

```python
import numpy as np
a = t.arrays(["pflow_theta", "pflow_pdgId", "truth_theta", "truth_pdgId"], entry_stop=200000)
bins = np.deg2rad(np.arange(5, 36, 2.5))

def spectrum(theta, pid, want):
    v = ak.to_numpy(ak.flatten(theta[abs(pid) == want]))
    return np.histogram(v, bins=bins)[0]

reco  = spectrum(a["pflow_theta"], a["pflow_pdgId"], 211)
truth = spectrum(a["truth_theta"], a["truth_pdgId"], 211)
ratio = reco / np.maximum(truth, 1)
```

The truth histogram is restricted to the same 5°–35° window by the bin edges, so the ratio is
interpretable — with the §6 folding caveat: a few mirrored backward pions land in these bins.

### Multiplicity correlation, events with a reconstructed electron

```python
a = t.arrays(["npflow", "ntruth", "pflow_scatE"], entry_stop=200000)
sel = a["pflow_scatE"] == 1
np_, nt = ak.to_numpy(a["npflow"][sel]), ak.to_numpy(a["ntruth"][sel])
```

### Dump one event

```python
i = 4
a = t.arrays(entry_stop=i + 1)
for b in ["npflow", "pflow_pdgId", "pflow_pt", "pflow_theta", "pflow_phi", "pflow_scatE",
          "ntruth", "truth_pdgId", "truth_pt", "truth_theta", "truth_phi", "truth_scatE"]:
    print(f"{b:14s}", a[b][i].tolist() if a[b].ndim > 1 else a[b][i])
```

Event 4 of this file, for reference — 6 reconstructed particles (a neutron and five photons, no
electron, so `pflow_scatE == 0`) against 13 truth particles led by the scattered electron:

```
npflow         6
pflow_pdgId    [2112, 22, 22, 22, 22, 22]
pflow_pt       [0.301, 0.071, 0.047, 0.085, 0.045, 0.043]
pflow_theta    [0.395, 0.147, 0.151, 0.302, 0.312, 0.321]
pflow_scatE    0
ntruth         13
truth_pdgId    [11, 2212, -211, 321, 211, -211, 211, 22, -211, 211, -211, 22, 22]
truth_pt       [0.573, 0.287, 0.079, 0.207, 0.092, 0.203, 0.209, 0.101, 0.147, 0.258, 0.277, 0.133, 0.375]
truth_theta    [0.325, 0.158, 0.259, 0.134, 0.366, 0.220, 0.726, 0.260, 0.146, 0.314, 0.396, 0.379, 0.409]
truth_scatE    1
```

Note in passing how hard the task is: none of the six reconstructed objects is obviously any one of the
thirteen true particles.

### Larger-than-memory scans

```python
total = 0
for chunk in t.iterate(["npflow"], step_size="200 MB"):
    total += int(ak.sum(chunk["npflow"]))
```
