#!/usr/bin/env python3
"""Plot reco (pflow) vs truth kinematics per particle species from a prepped Parnassus file.

Reads evt_tree from one ROOT file written by macros/hipo2tree_parnassus.C (v7 layout, with pflow_p/truth_p)
and draws one row per species and one column per variable (p, theta, phi, vz, beta). Reco and truth are each
normalised to unit area. Values outside an axis range are piled into its edge bins, so outliers stay visible.
evt_tree is filled after CutManager::filter_particles, so the reco side shows particles after the analysis cuts.

Usage:
  python3 scripts/plot_kinematics_per_species.py prepped_data/v8_data_sample.root          # all events
  python3 scripts/plot_kinematics_per_species.py prepped_data/v8_data_sample.root 100000   # first 100k events
  python3 scripts/plot_kinematics_per_species.py in.root 100000 -o plots/my_plot.png

Default output: plots/<file stem>/kinematics_per_species.png in the repository.
Requires numpy, matplotlib and uproot (on ifarm: python3 -m pip install --user uproot matplotlib).
"""
import argparse
import os

import awkward as ak
import matplotlib
import numpy as np
import uproot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SPECIES = [(11, "e⁻"), (-11, "e⁺"), (211, "π⁺"), (-211, "π⁻"), (321, "K⁺"), (-321, "K⁻"),
           (2212, "p"), (22, "γ"), (2112, "n")]
# (branch suffix, axis label, axis range, log y, stored in radians)
VARIABLES = [
    ("p", "p [GeV]", (0, 11), False, False),
    ("theta", "θ [deg]", (0, 180), True, True),
    ("phi", "φ [deg]", (-180, 180), False, True),
    ("vz", "v_z [cm]", (-60, 60), True, False),
    ("beta", "β", (-0.5, 2.5), True, False),
]
NBINS = 80
# (branch prefix, legend label, color); colors are slots 1-2 of the dataviz reference palette
COLLECTIONS = [("truth", "truth", "#eb6834"), ("pflow", "reco (pflow)", "#2a78d6")]
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def fill_histograms(path, n_events):
    """Accumulate fixed-binning histograms chunk by chunk, so memory does not grow with the file size."""
    edges = {v: np.linspace(*rng, NBINS + 1) for v, _, rng, _, _ in VARIABLES}
    hists = {(c, pid, v): np.zeros(NBINS) for c, _, _ in COLLECTIONS for pid, _ in SPECIES for v, *_ in VARIABLES}
    counts = {(c, pid): 0 for c, _, _ in COLLECTIONS for pid, _ in SPECIES}
    nonfinite = {}
    branches = [f"{c}_{v}" for c, _, _ in COLLECTIONS for v in ["pdgId"] + [v for v, *_ in VARIABLES]]

    tree = uproot.open(path)["evt_tree"]
    stop = None if n_events is None or n_events < 0 else min(n_events, tree.num_entries)
    n_read = 0
    for arrays in tree.iterate(branches, entry_stop=stop, step_size="200 MB", library="ak"):
        n_read += len(arrays)
        for c, _, _ in COLLECTIONS:
            pdg = ak.to_numpy(ak.flatten(arrays[f"{c}_pdgId"]))
            values = {v: ak.to_numpy(ak.flatten(arrays[f"{c}_{v}"])) for v, *_ in VARIABLES}
            for pid, _ in SPECIES:
                m = pdg == pid
                counts[(c, pid)] += int(m.sum())
                for v, _, rng, _, in_rad in VARIABLES:
                    x = values[v][m]
                    finite = np.isfinite(x)
                    if not finite.all():
                        nonfinite[(c, pid, v)] = nonfinite.get((c, pid, v), 0) + int((~finite).sum())
                    x = x[finite]
                    if in_rad:
                        x = np.degrees(x)
                    x = np.clip(x, rng[0], np.nextafter(rng[1], rng[0]))  # overflow into the edge bins
                    hists[(c, pid, v)] += np.histogram(x, bins=edges[v])[0]
    return edges, hists, counts, nonfinite, n_read


def plot(path, edges, hists, counts, n_read, out):
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "text.color": INK, "axes.facecolor": SURFACE,
                         "figure.facecolor": SURFACE})
    fig, axs = plt.subplots(len(SPECIES), len(VARIABLES), figsize=(15, 2.0 * len(SPECIES)), constrained_layout=True)
    for i, (pid, name) in enumerate(SPECIES):
        for j, (v, label, _, logy, _) in enumerate(VARIABLES):
            ax = axs[i, j]
            e = edges[v]
            drawn = False
            for c, legend, color in COLLECTIONS:
                h = hists[(c, pid, v)]
                if h.sum() == 0:
                    continue
                density = h / (h.sum() * np.diff(e))
                # weights on bin starts reproduce the histogram; works with old matplotlib (no ax.stairs)
                ax.hist(e[:-1], bins=e, weights=density, histtype="step", lw=1.5, color=color,
                        label=f"{legend} ({counts[(c, pid)]:,})")
                drawn = True
            if logy and drawn:
                ax.set_yscale("log")
            ax.grid(color=GRID, lw=0.6)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if i == len(SPECIES) - 1:
                ax.set_xlabel(label)
            if j == 0:
                ax.set_ylabel(f"{name}\ndensity")
                if drawn:
                    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.suptitle(f"{os.path.basename(path)} ({n_read:,} events): reco vs truth per species. "
                 "Edge bins hold overflow; θ, v_z, β on log scale", color=INK)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=110)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root_file", help="prepped ROOT file with an evt_tree")
    parser.add_argument("n_events", nargs="?", type=int, default=-1, help="number of events to read (-1: all)")
    parser.add_argument("-o", "--output", help="output image (default: plots/<file stem>/kinematics_per_species.png)")
    args = parser.parse_args()

    out = args.output
    if out is None:
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        stem = os.path.splitext(os.path.basename(args.root_file))[0]
        out = os.path.join(repo, "plots", stem, "kinematics_per_species.png")

    edges, hists, counts, nonfinite, n_read = fill_histograms(args.root_file, args.n_events)
    print(f"Read {n_read:,} events from {args.root_file}")
    print(f"{'species':8s} {'reco':>10s} {'truth':>10s}")
    for pid, name in SPECIES:
        print(f"{name:8s} {counts[('pflow', pid)]:10,d} {counts[('truth', pid)]:10,d}")
    for (c, pid, v), n in sorted(nonfinite.items()):
        print(f"WARNING: {n} non-finite {c}_{v} values for pdgId {pid} (left out of the plot)")
    plot(args.root_file, edges, hists, counts, n_read, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
