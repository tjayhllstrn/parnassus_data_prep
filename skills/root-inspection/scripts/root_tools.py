#!/usr/bin/env python3
"""Inspect, slice and plot segments of large ROOT TTrees from this repository.

Invoke through the wrapper (it activates the ROOT conda environment):

    skills/root-inspection/scripts/root_inspect.sh <subcommand> [options]

Subcommands
    info      list trees / cycles / branches (and scalar ranges) in a ROOT file
    count     count how many entries pass a TTree selection
    select    write a *new* ROOT file holding only the entries passing a cut
    dump      print the selected events (event scalars + per-particle rows)
    plot      one histogram of an arbitrary expression, optionally within a cut
    plot-all  skim by cut, then run macros/plot_prepped_data.C on the skim

Safety model (deliberate, do not weaken):
    * input files are only ever opened READ; nothing is written in place
    * cut/expression strings are validated against the tree's real branch names
      and a math whitelist, and are passed to ROOT as Python strings (never
      through a shell), so they cannot inject shell or CINT code
    * outputs go under prepped_data/skims/ or plots/ (both gitignored) and an
      existing file is never clobbered without --force
    * an output path may never equal the input path
    * `select` refuses to copy more than --max-entries entries (default 500k)
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import os
import re
import subprocess
import sys

# --------------------------------------------------------------------------- #
# repo layout
# --------------------------------------------------------------------------- #

_HERE = os.path.dirname(os.path.realpath(__file__))
# realpath (not abspath): this repo is reachable through symlinked paths
# (e.g. ~/link_to_work/...), and mixing the two breaks path comparisons.
REPO_ROOT = os.path.realpath(
    os.environ.get("PARNASSUS_DATA_PREP_ROOT", os.path.join(_HERE, "..", "..", ".."))
)
SKIM_DIR = os.path.join(REPO_ROOT, "prepped_data", "skims")
PLOT_DIR = os.path.join(REPO_ROOT, "plots")
STYLE_MACRO = os.path.join(REPO_ROOT, "src", "style.C")
PLOT_MACRO = os.path.join(REPO_ROOT, "macros", "plot_prepped_data.C")

DEFAULT_TREE = "evt_tree"

# Variables plot_prepped_data.C compares, with its fixed axis ranges. Used by
# the built-in fallback when a file predates some of those branches.
COMPARE_VARS = [
    # (var, x label, xmin, xmax, nbins)
    ("pt", "p_{T} [GeV/c]", -0.1, 2.0, 100),
    ("phi", "#phi [rad]", -3.14159, 3.14159, 100),
    ("eta", "#eta", -1.0, 5.0, 100),
    ("vx", "v_{x} [cm]", -4.0, 4.0, 100),
    ("vy", "v_{y} [cm]", -4.0, 4.0, 100),
    ("vz", "v_{z} [cm]", -10.0, 10.0, 100),
    ("beta", "#beta", 0.0, 1.1, 100),
    ("charge", "charge", -1.0, 2.0, 3),
]
CATEGORY_VARS = ["class", "pdgId"]

CLASS_NAMES = {
    0: "charged hadron",
    1: "electron",
    2: "muon",
    3: "neutral hadron",
    4: "photon",
    -1: "unknown",
}

# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #


class ToolError(Exception):
    """User-facing error: printed without a traceback."""


def info(msg: str) -> None:
    print(msg, flush=True)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr, flush=True)


def import_root():
    try:
        import ROOT  # noqa: PLC0415
    except ImportError as exc:  # pragma: no cover
        raise ToolError(
            "PyROOT is not importable. Run this tool through "
            "skills/root-inspection/scripts/root_inspect.sh so the "
            "parnassus-core conda environment is activated."
        ) from exc
    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    return ROOT


def open_read(path: str):
    """Open a ROOT file strictly read-only."""
    ROOT = import_root()
    if not os.path.isfile(path):
        raise ToolError(f"not a readable file: {path}")
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise ToolError(f"could not open ROOT file: {path}")
    return f


def tree_inventory(f):
    """Return [(name, cycle, entries)] for every TTree key, newest cycle first."""
    out = []
    for key in f.GetListOfKeys():
        if not key.GetClassName().startswith("TTree"):
            continue
        obj = key.ReadObj()
        out.append((key.GetName(), key.GetCycle(), int(obj.GetEntries())))
    out.sort(key=lambda t: (t[0], -t[1]))
    return out


def get_tree(f, name: str, cycle: int | None = None):
    """Fetch a tree by name (and optional cycle).

    Prepped files in this repo sometimes contain several *cycles* of the same
    tree (e.g. evt_tree;1 and evt_tree;2) from repeated Write() calls. Bare
    Get() silently returns the highest cycle, which is a real trap, so cycles
    are reported and can be pinned with --cycle.
    """
    key = f"{name};{cycle}" if cycle else name
    tree = f.Get(key)
    if not tree or not tree.InheritsFrom("TTree"):
        available = sorted({n for n, _, _ in tree_inventory(f)})
        raise ToolError(
            f"no TTree named {key!r} in the file. Available trees: "
            + (", ".join(available) if available else "(none)")
        )
    cycles = [c for n, c, _ in tree_inventory(f) if n == name]
    if cycle is None and len(cycles) > 1:
        warn(
            f"'{name}' exists in {len(cycles)} cycles {sorted(cycles)} with "
            f"different entry counts; using the highest ({max(cycles)}). "
            "Pin one with --cycle N."
        )
    return tree


def branch_map(tree) -> dict[str, dict]:
    """name -> {'type': root type char, 'counter': leaf-count branch or None}."""
    out = {}
    for br in tree.GetListOfBranches():
        name = br.GetName()
        leaf = br.GetLeaf(name) or (br.GetListOfLeaves().At(0) if br.GetListOfLeaves().GetEntries() else None)
        if leaf is None:
            out[name] = {"type": "?", "counter": None}
            continue
        counter = leaf.GetLeafCount()
        out[name] = {
            "type": leaf.GetTypeName(),
            "counter": counter.GetName() if counter else None,
        }
    return out


# --------------------------------------------------------------------------- #
# expression / cut validation
# --------------------------------------------------------------------------- #

# TTreeFormula pseudo-functions and math available in ROOT selections.
_ALLOWED_IDENTIFIERS = {
    "abs", "fabs", "sqrt", "exp", "log", "log10", "log2", "pow", "sin", "cos",
    "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh", "tanh", "floor",
    "ceil", "round", "int", "fmod", "min", "max", "sign", "TMath", "Pi", "Abs",
    "Sqrt", "Power", "Hypot", "Sq", "Entry$", "Entries$", "Length$", "Iteration$",
    "Sum$", "Min$", "Max$", "MinIf$", "MaxIf$", "Alt$", "Length", "kTRUE",
    "kFALSE", "true", "false", "pi", "e",
}

_IDENT_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\$?")
# Characters a legitimate TTree formula needs. Note: no ';', no backtick, no
# quotes, no braces, no backslash, no '#'.
_ALLOWED_CHARS_RE = re.compile(r"^[A-Za-z0-9_$.:+\-*/%()\[\]<>=!&|^~?, \t]*$")
_FORBIDDEN_SUBSTRINGS = (">>", "<<", "::system", "system(", "gSystem", "gROOT",
                         "ProcessLine", "Exec", "TFile", "unlink", "remove(")


def validate_formula(expr: str, tree, kind: str = "cut", allow_colon: bool = False) -> str:
    """Validate a TTree selection/expression against the tree's real branches.

    Rejects anything that is not plain formula syntax, and rejects identifiers
    that are neither a branch/leaf of the tree nor a whitelisted math symbol —
    which also catches typos (e.g. 'eventnumber' vs 'eventNumber') before ROOT
    silently evaluates them as 0.
    """
    expr = (expr or "").strip()
    if not expr:
        raise ToolError(f"empty {kind} expression")
    if len(expr) > 2000:
        raise ToolError(f"{kind} expression is unreasonably long ({len(expr)} chars)")
    if not _ALLOWED_CHARS_RE.match(expr):
        bad = sorted({c for c in expr if not _ALLOWED_CHARS_RE.match(c)})
        raise ToolError(
            f"{kind} contains characters that are not valid in a TTree formula: "
            f"{bad!r}. Only branch names, numbers and math operators are allowed."
        )
    lowered = expr.lower()
    for token in _FORBIDDEN_SUBSTRINGS:
        if token.lower() in lowered:
            raise ToolError(f"{kind} contains a forbidden token {token!r}")
    if ":" in expr and not allow_colon:
        raise ToolError(f"{kind} may not contain ':' (that is 2D-draw syntax)")
    if expr.count("(") != expr.count(")") or expr.count("[") != expr.count("]"):
        raise ToolError(f"unbalanced brackets in {kind}: {expr}")

    known = set(branch_map(tree))
    for leaf in tree.GetListOfLeaves():
        known.add(leaf.GetName())
    for alias in (tree.GetListOfAliases() or []):
        known.add(alias.GetName())

    pos = 0
    for match in _IDENT_RE.finditer(expr):
        ident = match.group(0)
        pos = match.end()
        # skip a member access such as TMath::Pi -> 'Pi' handled by whitelist
        if ident in known or ident in _ALLOWED_IDENTIFIERS:
            continue
        # a numeric suffix like 1e5 lexes as 'e5'
        if re.fullmatch(r"[eE]\d+", ident) and match.start() > 0 and expr[match.start() - 1].isdigit():
            continue
        suggestion = difflib.get_close_matches(ident, sorted(known), n=3, cutoff=0.6)
        hint = f" Did you mean: {', '.join(suggestion)}?" if suggestion else ""
        raise ToolError(
            f"unknown identifier {ident!r} in {kind}; it is not a branch of "
            f"'{tree.GetName()}'.{hint}"
        )
    del pos
    return expr


def cut_tag(cut: str | None) -> str:
    """Short, filesystem-safe tag for a cut string."""
    if not cut:
        return "all"
    slug = re.sub(r"[^A-Za-z0-9]+", "_", cut).strip("_")[:60]
    digest = hashlib.sha1(cut.encode()).hexdigest()[:6]
    return f"{slug}_{digest}" if slug else digest


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "x"


def file_base(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


# --------------------------------------------------------------------------- #
# output-path guards
# --------------------------------------------------------------------------- #


def _under(path: str, root: str) -> bool:
    path = os.path.realpath(path)
    root = os.path.realpath(root)
    return path == root or path.startswith(root + os.sep)


def resolve_output(out: str, inputs: list[str], allowed_roots: list[str],
                   suffixes: tuple[str, ...], force: bool) -> str:
    """Validate an output path: extension, location, no-clobber, not an input."""
    out = os.path.realpath(out)
    if not out.endswith(suffixes):
        raise ToolError(f"output must end in one of {suffixes}: {out}")
    for src in inputs:
        if os.path.exists(src) and os.path.exists(out) and os.path.samefile(src, out):
            raise ToolError(f"output would overwrite the input file: {out}")
        if os.path.realpath(src) == out:
            raise ToolError(f"output would overwrite the input file: {out}")
    if not any(_under(out, root) for root in allowed_roots) and not force:
        pretty = ", ".join(os.path.relpath(r, REPO_ROOT) for r in allowed_roots)
        raise ToolError(
            f"refusing to write outside {pretty} without --force: {out}"
        )
    if os.path.exists(out) and not force:
        raise ToolError(f"output already exists (use --force to overwrite): {out}")
    if os.path.isdir(out):
        raise ToolError(f"output path is a directory: {out}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    return out


def rel(path: str) -> str:
    """Repo-relative path for display (falls back to the absolute path)."""
    try:
        return os.path.relpath(os.path.realpath(path), REPO_ROOT)
    except ValueError:  # pragma: no cover
        return path


# --------------------------------------------------------------------------- #
# selection
# --------------------------------------------------------------------------- #


def selected_entries(tree, cut: str | None, limit: int | None = None) -> list[int]:
    """Entry numbers passing `cut` (all entries when cut is None).

    Uses a TEntryList, i.e. the same TTreeFormula semantics as TTree::CopyTree.
    NOTE: for per-particle array branches an entry is selected if the formula is
    true for *any* element; wrap per-particle conditions in Sum$()/Max$() if you
    want explicit event-level logic.
    """
    ROOT = import_root()
    n_total = int(tree.GetEntries())
    if not cut:
        entries = range(n_total)
        return list(entries)[:limit] if limit else list(entries)

    ROOT.gROOT.cd()  # keep the temporary entry list out of the input file
    listname = "ri_elist"
    tree.SetEntryList(ROOT.nullptr)
    tree.Draw(f">>{listname}", cut, "entrylist")
    elist = ROOT.gROOT.FindObject(listname) or ROOT.gDirectory.Get(listname)
    if not elist:
        raise ToolError(f"ROOT could not evaluate the selection: {cut!r}")
    n_sel = int(elist.GetN())
    out = []
    for i in range(n_sel if limit is None else min(n_sel, limit)):
        out.append(int(elist.GetEntry(i)))
    return out


def count_selection(tree, cut: str | None) -> int:
    if not cut:
        return int(tree.GetEntries())
    return len(selected_entries(tree, cut))


# --------------------------------------------------------------------------- #
# subcommand: info
# --------------------------------------------------------------------------- #


def cmd_info(args) -> int:
    f = open_read(args.file)
    inv = tree_inventory(f)
    size_mb = os.path.getsize(args.file) / 1e6
    info(f"file    : {rel(os.path.abspath(args.file))}  ({size_mb:.1f} MB)")
    if not inv:
        info("no TTrees found in this file")
        return 0
    info("trees   :")
    for name, cycle, entries in inv:
        info(f"  {name};{cycle:<3d} {entries:>12,d} entries")

    tree = get_tree(f, args.tree, args.cycle)
    bmap = branch_map(tree)
    info("")
    info(f"tree '{tree.GetName()}' : {int(tree.GetEntries()):,d} entries, "
         f"{len(bmap)} branches")
    scalars = [n for n, m in bmap.items() if m["counter"] is None]
    arrays = [n for n, m in bmap.items() if m["counter"] is not None]

    info("  event-level branches:")
    for name in scalars:
        line = f"    {name:<24s} {bmap[name]['type']:<10s}"
        if args.values:
            line += f" min={tree.GetMinimum(name):<14.6g} max={tree.GetMaximum(name):<14.6g}"
        info(line)
    if arrays:
        info("  per-particle branches (variable length):")
        by_counter: dict[str, list[str]] = {}
        for name in arrays:
            by_counter.setdefault(bmap[name]["counter"], []).append(name)
        for counter, names in by_counter.items():
            info(f"    [{counter}]  " + ", ".join(sorted(names)))
    return 0


# --------------------------------------------------------------------------- #
# subcommand: count
# --------------------------------------------------------------------------- #


def cmd_count(args) -> int:
    f = open_read(args.file)
    tree = get_tree(f, args.tree, args.cycle)
    cut = validate_formula(args.cut, tree, "cut") if args.cut else None
    total = int(tree.GetEntries())
    n = count_selection(tree, cut)
    frac = 100.0 * n / total if total else 0.0
    info(f"tree '{tree.GetName()}': {total:,d} entries")
    info(f"cut  : {cut or '(none)'}")
    info(f"pass : {n:,d} entries ({frac:.3f}%)")
    return 0


# --------------------------------------------------------------------------- #
# subcommand: select
# --------------------------------------------------------------------------- #


def default_skim_path(src: str, cut: str | None, label: str | None) -> str:
    tag = safe_name(label) if label else cut_tag(cut)
    return os.path.join(SKIM_DIR, f"{file_base(src)}__{tag}.root")


def do_select(src: str, tree_name: str, cycle: int | None, cut: str | None,
              out: str | None, label: str | None, max_entries: int,
              dry_run: bool, force: bool, extra_trees: list[str] | None = None) -> str | None:
    """Copy the entries passing `cut` into a brand-new ROOT file."""
    ROOT = import_root()
    fin = open_read(src)
    tree = get_tree(fin, tree_name, cycle)
    vcut = validate_formula(cut, tree, "cut") if cut else None

    total = int(tree.GetEntries())
    n_sel = count_selection(tree, vcut)
    info(f"source   : {rel(os.path.abspath(src))}  (tree '{tree_name}', {total:,d} entries)")
    info(f"cut      : {vcut or '(none)'}")
    info(f"selected : {n_sel:,d} entries")

    if n_sel == 0:
        raise ToolError(
            "the cut selects 0 entries. Check the branch ranges with the "
            "`info --values` subcommand (e.g. eventNumber is per-file and is "
            "only unique/large in merged files)."
        )
    if n_sel > max_entries:
        raise ToolError(
            f"selection has {n_sel:,d} entries, above the safety cap of "
            f"{max_entries:,d}. Tighten the cut or raise --max-entries "
            "deliberately."
        )

    out_path = os.path.abspath(out) if out else default_skim_path(src, cut, label)
    if dry_run:
        info(f"dry run  : would write {rel(out_path)}")
        return None

    out_path = resolve_output(
        out_path, [src], [SKIM_DIR, os.path.join(REPO_ROOT, "prepped_data")],
        (".root",), force
    )

    fout = ROOT.TFile(out_path, "RECREATE")
    if fout.IsZombie():
        raise ToolError(f"could not create output file: {out_path}")
    fout.cd()
    copied = tree.CopyTree(vcut or "")
    if not copied:
        raise ToolError("CopyTree returned null; nothing was written")
    copied.SetName(tree_name)
    n_out = int(copied.GetEntries())
    copied.Write(tree_name, ROOT.TObject.kOverwrite)

    for name in extra_trees or []:
        other = fin.Get(name)
        if not other or not other.InheritsFrom("TTree"):
            warn(f"skipping extra tree '{name}': not found in {rel(src)}")
            continue
        try:
            validate_formula(vcut, other, "cut") if vcut else None
        except ToolError as exc:
            warn(f"skipping extra tree '{name}': {exc}")
            continue
        fout.cd()
        clone = other.CopyTree(vcut or "")
        if clone:
            clone.SetName(name)
            clone.Write(name, ROOT.TObject.kOverwrite)
            info(f"           + tree '{name}': {int(clone.GetEntries()):,d} entries")

    # provenance, so a skim is never a mystery file later on
    prov = ROOT.TNamed(
        "skim_info",
        f"source={os.path.abspath(src)}; tree={tree_name};"
        f"cycle={cycle if cycle else 'latest'}; cut={vcut or 'none'}",
    )
    prov.Write()
    fout.Close()

    if n_out != n_sel:
        warn(f"CopyTree wrote {n_out:,d} entries but the entry list had {n_sel:,d}")
    info(f"wrote    : {rel(out_path)}  ({n_out:,d} entries, "
         f"{os.path.getsize(out_path) / 1e6:.1f} MB)")
    return out_path


def cmd_select(args) -> int:
    do_select(args.file, args.tree, args.cycle, args.cut, args.out, args.label,
              args.max_entries, args.dry_run, args.force,
              extra_trees=split_list(args.also_trees))
    return 0


def split_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [v.strip() for v in value.split(",") if v.strip()]


# --------------------------------------------------------------------------- #
# subcommand: dump
# --------------------------------------------------------------------------- #


def _fmt(value) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def cmd_dump(args) -> int:
    f = open_read(args.file)
    tree = get_tree(f, args.tree, args.cycle)
    cut = validate_formula(args.cut, tree, "cut") if args.cut else None
    bmap = branch_map(tree)

    n_sel_total = count_selection(tree, cut)
    entries = selected_entries(tree, cut, limit=args.max_events)
    info(f"tree '{tree.GetName()}' : {int(tree.GetEntries()):,d} entries")
    info(f"cut  : {cut or '(none)'}")
    info(f"pass : {n_sel_total:,d} entries; showing {len(entries)}"
         + (f" (raise with --max-events)" if n_sel_total > len(entries) else ""))

    wanted = split_list(args.vars)
    for name in wanted:
        if name not in bmap:
            raise ToolError(f"--vars: '{name}' is not a branch of {tree.GetName()}")

    scalars = [n for n, m in bmap.items() if m["counter"] is None]
    if wanted:
        scalars = [n for n in scalars if n in wanted]
    groups: dict[str, list[str]] = {}
    for name, meta in bmap.items():
        if meta["counter"] is None:
            continue
        if wanted and name not in wanted:
            continue
        groups.setdefault(meta["counter"], []).append(name)
    for counter in groups:
        groups[counter].sort()

    for entry in entries:
        tree.GetEntry(entry)
        head = "  ".join(f"{n}={_fmt(getattr(tree, n))}" for n in scalars)
        info("")
        info(f"--- entry {entry}  {head}")
        if args.no_particles:
            continue
        for counter, names in groups.items():
            n_obj = int(getattr(tree, counter))
            if n_obj <= 0:
                info(f"    [{counter}] empty")
                continue
            header = "    " + f"{'idx':>4s} " + " ".join(f"{n:>14s}" for n in names)
            info(header)
            for j in range(min(n_obj, args.max_particles)):
                row = " ".join(f"{_fmt(getattr(tree, n)[j]):>14s}" for n in names)
                extra = ""
                for n in names:
                    if n.endswith("_class"):
                        extra = f"  # {CLASS_NAMES.get(int(getattr(tree, n)[j]), '?')}"
                info(f"    {j:>4d} {row}{extra}")
            if n_obj > args.max_particles:
                info(f"    ... {n_obj - args.max_particles} more (raise --max-particles)")
    return 0


# --------------------------------------------------------------------------- #
# plotting helpers
# --------------------------------------------------------------------------- #


def load_style():
    """Apply the repository's plotting style (src/style.C) if available."""
    ROOT = import_root()
    colors = {"blue": ROOT.kBlue + 1, "magenta": ROOT.kMagenta + 1}
    if os.path.isfile(STYLE_MACRO):
        try:
            ROOT.gROOT.LoadMacro(STYLE_MACRO)
            ROOT.SetReleaseNoteStyle()
            ROOT.DefineColors()
            colors = {"blue": int(ROOT.BLUE), "magenta": int(ROOT.MAGENTA)}
        except Exception as exc:  # style is cosmetic; never fail the plot on it
            warn(f"could not apply src/style.C ({exc}); using ROOT defaults")
    return colors


def _find_object(tree, name: str):
    """Locate an object TTree::Draw just created (its home dir varies)."""
    ROOT = import_root()
    for holder in (ROOT.gDirectory, ROOT.gROOT, tree.GetDirectory()):
        if not holder:
            continue
        obj = holder.Get(name) if hasattr(holder, "Get") else None
        if obj:
            return obj
    return ROOT.gROOT.FindObject(name)


def _drop_object(name: str) -> None:
    """Remove a stale same-named histogram so Draw does not append to it."""
    ROOT = import_root()
    ROOT.gROOT.cd()
    for holder in (ROOT.gDirectory, ROOT.gROOT):
        obj = holder.Get(name)
        if obj:
            obj.Delete()


def auto_range(tree, expr: str, cut: str | None):
    """Let ROOT pick a sensible axis range for `expr` (with `cut` applied)."""
    ROOT = import_root()
    _drop_object("htemp")
    n = tree.Draw(expr, cut or "", "goff")
    if n <= 0:
        raise ToolError(f"expression {expr!r} with cut {cut or '(none)'} yields no entries")
    htemp = _find_object(tree, "htemp")
    if htemp:
        axis = htemp.GetXaxis()
        return axis.GetXmin(), axis.GetXmax()
    values = [tree.GetV1()[i] for i in range(min(n, int(tree.GetEstimate())))]
    lo, hi = min(values), max(values)
    if lo == hi:
        pad = abs(lo) * 0.1 or 1.0
        lo, hi = lo - pad, hi + pad
    return lo, hi


def draw_1d(tree, expr: str, cut: str | None, name: str, bins: int,
            xmin: float | None, xmax: float | None):
    """1D histogram of `expr` under `cut`, detached from any file directory."""
    ROOT = import_root()
    if xmin is None or xmax is None:
        xmin, xmax = auto_range(tree, expr, cut)
    _drop_object(name)
    tree.Draw(f"{expr}>>{name}({bins},{xmin},{xmax})", cut or "", "goff")
    hist = _find_object(tree, name)
    if not hist:
        raise ToolError(f"ROOT did not produce a histogram for {expr!r}")
    hist.SetDirectory(ROOT.nullptr)
    return hist


def draw_2d(tree, yexpr: str, xexpr: str, cut: str | None, name: str, bins: int,
            xrange_: tuple[float | None, float | None],
            yrange_: tuple[float | None, float | None]):
    """2D histogram of `yexpr:xexpr` under `cut`."""
    ROOT = import_root()
    xlo, xhi = xrange_
    ylo, yhi = yrange_
    if xlo is None or xhi is None:
        xlo, xhi = auto_range(tree, xexpr, cut)
    if ylo is None or yhi is None:
        ylo, yhi = auto_range(tree, yexpr, cut)
    _drop_object(name)
    spec = f"({bins},{xlo},{xhi},{bins},{ylo},{yhi})"
    tree.Draw(f"{yexpr}:{xexpr}>>{name}{spec}", cut or "", "goff")
    hist = _find_object(tree, name)
    if not hist:
        raise ToolError(f"ROOT did not produce a histogram for {yexpr}:{xexpr}")
    hist.SetDirectory(ROOT.nullptr)
    return hist


def style_hist(hist, color, marker):
    hist.SetLineColor(color)
    hist.SetFillColorAlpha(color, 0.7)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetMarkerSize(0.7)
    hist.SetLineWidth(2)


def normalize(hist):
    integral = hist.Integral("width")
    if integral > 0:
        hist.Scale(1.0 / integral)


def save_canvas(canvas, out_png: str, hists: dict | None = None,
                save_root: bool = False):
    ROOT = import_root()
    canvas.SaveAs(out_png)
    info(f"wrote    : {rel(out_png)}")
    if save_root and hists:
        out_root = os.path.splitext(out_png)[0] + ".root"
        fout = ROOT.TFile(out_root, "RECREATE")
        if not fout.IsZombie():
            for key, hist in hists.items():
                hist.Write(key)
            fout.Close()
            info(f"wrote    : {rel(out_root)}")


def category_counts(tree, branch: str, cut: str | None) -> dict[int, int]:
    """Per-object counts of an integer per-particle branch within a selection."""
    counts: dict[int, int] = {}
    counter = branch_map(tree)[branch]["counter"]
    for entry in selected_entries(tree, cut):
        tree.GetEntry(entry)
        n_obj = int(getattr(tree, counter)) if counter else 1
        values = getattr(tree, branch)
        for j in range(n_obj):
            key = int(values[j]) if counter else int(values)
            counts[key] = counts.get(key, 0) + 1
    return counts


# --------------------------------------------------------------------------- #
# subcommand: plot
# --------------------------------------------------------------------------- #


def plot_outdir(src: str, cut: str | None, label: str | None) -> str:
    base = file_base(src)
    sub = safe_name(label) if label else (cut_tag(cut) if cut else "all")
    return os.path.join(PLOT_DIR, base, sub)


def cmd_plot(args) -> int:
    ROOT = import_root()
    f = open_read(args.file)
    tree = get_tree(f, args.tree, args.cycle)
    cut = validate_formula(args.cut, tree, "cut") if args.cut else None
    colors = load_style()

    is_2d = ":" in args.expr
    if args.compare:
        if is_2d:
            raise ToolError("--compare cannot be combined with a 2D 'y:x' expression")
        # --compare takes the bare variable stem (pt, eta, ...) or either of the
        # pflow_/truth_ branch names, and expands to both.
        stem = args.expr.strip()
        if stem.startswith(("pflow_", "truth_")):
            stem = stem.split("_", 1)[1]
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", stem):
            raise ToolError(
                f"--compare needs a plain variable name (e.g. pt, eta, pflow_pt), got {args.expr!r}"
            )
        reco_expr = validate_formula(f"pflow_{stem}", tree, "expression")
        truth_expr = validate_formula(f"truth_{stem}", tree, "expression")
        expr = reco_expr
    else:
        expr = validate_formula(args.expr, tree, "expression", allow_colon=is_2d)

    normalized = args.norm if args.norm is not None else bool(args.compare)
    n_sel = count_selection(tree, cut)
    info(f"tree '{tree.GetName()}' : {int(tree.GetEntries()):,d} entries; "
         f"cut selects {n_sel:,d}")
    if n_sel == 0:
        raise ToolError("the cut selects 0 entries; nothing to plot")

    outdir = plot_outdir(args.file, cut, args.label)
    if args.out:
        out_png = args.out
    else:
        stem = safe_name(args.expr.replace(":", "_vs_"))
        stem = ("compare_" if args.compare else "") + stem
        out_png = os.path.join(outdir, stem + args.format)
    out_png = resolve_output(out_png, [args.file], [PLOT_DIR],
                             (".png", ".pdf", ".svg"), args.force)

    canvas = ROOT.TCanvas("ri_canvas", "", 1200, 900)
    canvas.cd()
    if args.logx:
        canvas.SetLogx()
    if args.logy:
        canvas.SetLogy()
    hists: dict = {}

    if args.compare:
        xmin, xmax = args.xmin, args.xmax
        if xmin is None or xmax is None:
            lo_r, hi_r = auto_range(tree, reco_expr, cut)
            lo_t, hi_t = auto_range(tree, truth_expr, cut)
            xmin, xmax = min(lo_r, lo_t), max(hi_r, hi_t)
        h_reco = draw_1d(tree, reco_expr, cut, "h_reco", args.bins, xmin, xmax)
        h_truth = draw_1d(tree, truth_expr, cut, "h_truth", args.bins, xmin, xmax)
        style_hist(h_reco, colors["blue"], 20)
        style_hist(h_truth, colors["magenta"], 24)
        if normalized:
            normalize(h_reco)
            normalize(h_truth)
        title = args.title or f"pflow vs truth {args.expr}"
        h_reco.SetTitle(title)
        h_reco.GetXaxis().SetTitle(args.xlabel or args.expr)
        h_reco.GetYaxis().SetTitle(args.ylabel or ("Normalized entries" if normalized else "Entries"))
        h_reco.SetMaximum(1.25 * max(h_reco.GetMaximum(), h_truth.GetMaximum()))
        h_reco.Draw("hist")
        h_truth.Draw("hist same")
        legend = ROOT.TLegend(0.68, 0.78, 0.90, 0.90)
        legend.AddEntry(h_reco, reco_expr, "f")
        legend.AddEntry(h_truth, truth_expr, "f")
        legend.Draw()
        info(f"entries  : {reco_expr}={h_reco.GetEntries():.0f}  "
             f"{truth_expr}={h_truth.GetEntries():.0f}")
        hists = {"reco_hist": h_reco, "truth_hist": h_truth}
        save_canvas(canvas, out_png, hists, args.save_root)
        return 0

    if is_2d:
        ypart, xpart = [p.strip() for p in expr.split(":", 1)]
        hist = draw_2d(tree, ypart, xpart, cut, "ri_h2", args.bins,
                       (args.xmin, args.xmax), (args.ymin, args.ymax))
        hist.SetTitle(args.title or args.expr)
        hist.GetXaxis().SetTitle(args.xlabel or xpart)
        hist.GetYaxis().SetTitle(args.ylabel or ypart)
        if args.logz:
            canvas.SetLogz()
        canvas.SetRightMargin(0.15)
        hist.Draw("colz")
        info(f"entries  : {hist.GetEntries():.0f}  mean=({hist.GetMean(1):.4g}, {hist.GetMean(2):.4g})")
        save_canvas(canvas, out_png, {"hist": hist}, args.save_root)
        return 0

    hist = draw_1d(tree, expr, cut, "ri_h1", args.bins, args.xmin, args.xmax)
    style_hist(hist, colors["blue"], 20)
    if normalized:
        normalize(hist)
    hist.SetTitle(args.title or args.expr)
    hist.GetXaxis().SetTitle(args.xlabel or args.expr)
    hist.GetYaxis().SetTitle(args.ylabel or ("Normalized entries" if normalized else "Entries"))
    hist.Draw("hist")
    info(f"entries  : {hist.GetEntries():.0f}  mean={hist.GetMean():.6g}  "
         f"rms={hist.GetRMS():.6g}  underflow={hist.GetBinContent(0):.0f}  "
         f"overflow={hist.GetBinContent(hist.GetNbinsX() + 1):.0f}")
    save_canvas(canvas, out_png, {"hist": hist}, args.save_root)
    return 0


# --------------------------------------------------------------------------- #
# subcommand: plot-all
# --------------------------------------------------------------------------- #


def macro_required_branches() -> list[str]:
    needed = ["npflow", "ntruth"]
    for var, _, _, _, _ in COMPARE_VARS:
        needed += [f"pflow_{var}", f"truth_{var}"]
    for var in CATEGORY_VARS:
        needed += [f"pflow_{var}", f"truth_{var}"]
    return needed


def run_repo_plot_macro(target: str) -> int:
    """Run macros/plot_prepped_data.C on `target` from the repo root."""
    if not os.path.isfile(PLOT_MACRO):
        raise ToolError(f"missing plotting macro: {PLOT_MACRO}")
    arg = rel(os.path.abspath(target))
    cmd = ["root", "-l", "-b", "-q", f'macros/plot_prepped_data.C("{arg}")']
    info(f"running  : {' '.join(cmd)}   (cwd={REPO_ROOT})")
    proc = subprocess.run(cmd, cwd=REPO_ROOT, check=False)
    outdir = os.path.join(PLOT_DIR, file_base(target))
    if proc.returncode != 0:
        raise ToolError(f"plot_prepped_data.C exited with code {proc.returncode}")
    if os.path.isdir(outdir):
        pngs = sorted(p for p in os.listdir(outdir) if p.endswith(".png"))
        info(f"plots    : {rel(outdir)}/  ({len(pngs)} png)")
        for png in pngs:
            info(f"           {png}")
    return proc.returncode


def builtin_plot_set(src: str, tree_name: str, cycle: int | None,
                     cut: str | None, outdir: str, force: bool) -> None:
    """Fallback for files that lack some branches the repo macro requires."""
    ROOT = import_root()
    f = open_read(src)
    tree = get_tree(f, tree_name, cycle)
    vcut = validate_formula(cut, tree, "cut") if cut else None
    present = set(branch_map(tree))
    colors = load_style()
    os.makedirs(outdir, exist_ok=True)

    made = 0
    for var, xlabel, xmin, xmax, bins in COMPARE_VARS:
        reco, truth = f"pflow_{var}", f"truth_{var}"
        if reco not in present or truth not in present:
            warn(f"skipping '{var}': {reco}/{truth} not in this file")
            continue
        h_reco = draw_1d(tree, reco, vcut, f"h_reco_{var}", bins, xmin, xmax)
        h_truth = draw_1d(tree, truth, vcut, f"h_truth_{var}", bins, xmin, xmax)
        style_hist(h_reco, colors["blue"], 20)
        style_hist(h_truth, colors["magenta"], 24)
        normalize(h_reco)
        normalize(h_truth)
        canvas = ROOT.TCanvas(f"c_{var}", "", 1200, 900)
        h_reco.SetTitle(f"pflow vs truth {var}")
        h_reco.GetXaxis().SetTitle(xlabel)
        h_reco.GetYaxis().SetTitle("Normalized entries")
        h_reco.SetMaximum(1.25 * max(h_reco.GetMaximum(), h_truth.GetMaximum()))
        h_reco.Draw("hist")
        h_truth.Draw("hist same")
        legend = ROOT.TLegend(0.68, 0.78, 0.90, 0.90)
        legend.AddEntry(h_reco, reco, "f")
        legend.AddEntry(h_truth, truth, "f")
        legend.Draw()
        out_png = resolve_output(os.path.join(outdir, f"compare_{var}.png"),
                                 [src], [PLOT_DIR], (".png",), True)
        canvas.SaveAs(out_png)
        made += 1

    if {"npflow", "ntruth"} <= present:
        h_reco = draw_1d(tree, "npflow", vcut, "h_npflow", 20, None, None)
        h_truth = draw_1d(tree, "ntruth", vcut, "h_ntruth", 20,
                          h_reco.GetXaxis().GetXmin(), h_reco.GetXaxis().GetXmax())
        style_hist(h_reco, colors["blue"], 20)
        style_hist(h_truth, colors["magenta"], 24)
        normalize(h_reco)
        normalize(h_truth)
        canvas = ROOT.TCanvas("c_mult", "", 1200, 900)
        h_reco.SetTitle("npflow vs ntruth")
        h_reco.GetXaxis().SetTitle("multiplicity")
        h_reco.GetYaxis().SetTitle("Normalized entries")
        h_reco.SetMaximum(1.25 * max(h_reco.GetMaximum(), h_truth.GetMaximum()))
        h_reco.Draw("hist")
        h_truth.Draw("hist same")
        legend = ROOT.TLegend(0.68, 0.78, 0.90, 0.90)
        legend.AddEntry(h_reco, "npflow", "f")
        legend.AddEntry(h_truth, "ntruth", "f")
        legend.Draw()
        canvas.SaveAs(os.path.join(outdir, "compare_npflow_ntruth.png"))
        made += 1

    for var in CATEGORY_VARS:
        reco, truth = f"pflow_{var}", f"truth_{var}"
        if reco not in present or truth not in present:
            continue
        for branch in (reco, truth):
            counts = category_counts(tree, branch, vcut)
            total = sum(counts.values()) or 1
            info(f"{branch} composition ({total:,d} objects):")
            for key in sorted(counts, key=lambda k: -counts[k]):
                extra = f" ({CLASS_NAMES.get(key, '?')})" if var == "class" else ""
                info(f"    {key:>12d}{extra:<18s} {counts[key]:>10,d}  "
                     f"{100.0 * counts[key] / total:6.2f}%")
    info(f"plots    : {rel(outdir)}/  ({made} png)")


def missing_macro_branches(src: str, tree_name: str, cycle: int | None) -> list[str]:
    f = open_read(src)
    tree = get_tree(f, tree_name, cycle)
    present = set(branch_map(tree))
    missing = [b for b in macro_required_branches() if b not in present]
    f.Close()
    return missing


def cmd_plot_all(args) -> int:
    missing = missing_macro_branches(args.file, args.tree, args.cycle)

    target = args.file
    skim = None
    if args.cut:
        skim = do_select(args.file, args.tree, args.cycle, args.cut, args.out_skim,
                         args.label, args.max_entries, False, args.force)
        target = skim
    else:
        info("no --cut given: plotting the whole file")

    try:
        if missing:
            warn(
                "macros/plot_prepped_data.C needs branches this file does not "
                f"have: {', '.join(missing)}"
            )
            if args.strict:
                raise ToolError("--strict: refusing to run the macro on an incompatible file")
            warn("falling back to the skill's built-in comparison plots")
            # A skim already contains exactly the selected entries, so the cut
            # must not be re-applied; plot into plots/<skim base>/ so the layout
            # matches what plot_prepped_data.C would have produced.
            outdir = (os.path.join(PLOT_DIR, file_base(target)) if skim
                      else plot_outdir(args.file, args.cut, args.label))
            builtin_plot_set(target, args.tree, None if skim else args.cycle,
                             None if skim else args.cut, outdir, args.force)
        else:
            run_repo_plot_macro(target)
    finally:
        if skim and not args.keep_skim:
            os.remove(skim)
            info(f"removed  : {rel(skim)} (skim not kept; pass --keep-skim to retain)")
    return 0


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def add_common(parser) -> None:
    parser.add_argument("file", help="input ROOT file (opened read-only)")
    parser.add_argument("--tree", default=DEFAULT_TREE,
                        help=f"tree name (default: {DEFAULT_TREE})")
    parser.add_argument("--cycle", type=int, default=None,
                        help="pin a specific TKey cycle of the tree, e.g. 1")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="root_inspect.sh",
        description="Inspect, slice and plot segments of large ROOT TTrees.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  root_inspect.sh info prepped_data/f.root --values\n"
            "  root_inspect.sh count prepped_data/f.root --cut 'eventNumber>3392000 && eventNumber<3395200'\n"
            "  root_inspect.sh plot-all prepped_data/f.root --cut 'eventNumber>3392000 && eventNumber<3395200' --label seg1\n"
            "  root_inspect.sh plot prepped_data/f.root --expr pt --compare --cut 'npflow>3'\n"
            "  root_inspect.sh dump prepped_data/f.root --cut 'eventNumber==3392001'\n"
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("info", help="list trees, cycles and branches")
    add_common(p)
    p.add_argument("--values", action="store_true",
                   help="also print min/max of every event-level branch")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("count", help="count entries passing a cut")
    add_common(p)
    p.add_argument("--cut", required=True, help="TTree selection string")
    p.set_defaults(func=cmd_count)

    p = sub.add_parser("select", help="write a new ROOT file with only the selected entries")
    add_common(p)
    p.add_argument("--cut", required=True, help="TTree selection string")
    p.add_argument("--out", default=None,
                   help="output .root path (default: prepped_data/skims/<base>__<tag>.root)")
    p.add_argument("--label", default=None, help="short name used in the default output path")
    p.add_argument("--also-trees", default=None,
                   help="comma-separated extra trees to skim with the same cut")
    p.add_argument("--max-entries", type=int, default=500_000,
                   help="safety cap on entries copied (default: 500000)")
    p.add_argument("--dry-run", action="store_true", help="report the count, write nothing")
    p.add_argument("--force", action="store_true", help="allow overwriting the output file")
    p.set_defaults(func=cmd_select)

    p = sub.add_parser("dump", help="print the selected events in full")
    add_common(p)
    p.add_argument("--cut", default=None, help="TTree selection string")
    p.add_argument("--max-events", type=int, default=20, help="events to print (default: 20)")
    p.add_argument("--max-particles", type=int, default=40,
                   help="per-event particle rows to print (default: 40)")
    p.add_argument("--vars", default=None, help="comma-separated subset of branches to show")
    p.add_argument("--no-particles", action="store_true", help="event scalars only")
    p.set_defaults(func=cmd_dump)

    p = sub.add_parser("plot", help="one histogram of an expression within a cut")
    add_common(p)
    p.add_argument("--expr", required=True,
                   help="branch or formula; 'y:x' makes a 2D colz plot")
    p.add_argument("--cut", default=None, help="TTree selection string")
    p.add_argument("--compare", action="store_true",
                   help="overlay pflow_<var> vs truth_<var> (pass --expr pt, eta, ...)")
    p.add_argument("--bins", type=int, default=100, help="bins per axis (default: 100)")
    p.add_argument("--xmin", type=float, default=None)
    p.add_argument("--xmax", type=float, default=None)
    p.add_argument("--ymin", type=float, default=None)
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--logx", action="store_true")
    p.add_argument("--logy", action="store_true")
    p.add_argument("--logz", action="store_true")
    norm = p.add_mutually_exclusive_group()
    norm.add_argument("--norm", dest="norm", action="store_true", default=None,
                      help="normalize to unit area (default for --compare)")
    norm.add_argument("--no-norm", dest="norm", action="store_false",
                      help="plot raw counts (default for a single histogram)")
    p.add_argument("--title", default=None)
    p.add_argument("--xlabel", default=None)
    p.add_argument("--ylabel", default=None)
    p.add_argument("--label", default=None, help="short name for the plots/ subdirectory")
    p.add_argument("--out", default=None, help="explicit output image path under plots/")
    p.add_argument("--format", default=".png", choices=[".png", ".pdf", ".svg"])
    p.add_argument("--save-root", action="store_true",
                   help="also write the histogram(s) to a sidecar .root")
    p.add_argument("--force", action="store_true", help="allow overwriting the output image")
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser("plot-all",
                       help="skim by cut, then run macros/plot_prepped_data.C on the skim")
    add_common(p)
    p.add_argument("--cut", default=None,
                   help="TTree selection string (omit to plot the whole file)")
    p.add_argument("--label", default=None, help="short name for skim/plot paths")
    p.add_argument("--out-skim", default=None, help="explicit skim .root path")
    p.add_argument("--max-entries", type=int, default=500_000,
                   help="safety cap on entries copied into the skim (default: 500000)")
    p.add_argument("--keep-skim", dest="keep_skim", action="store_true", default=True,
                   help="keep the skim file (default)")
    p.add_argument("--no-keep-skim", dest="keep_skim", action="store_false",
                   help="delete the skim after plotting")
    p.add_argument("--strict", action="store_true",
                   help="fail instead of falling back when the file lacks branches "
                        "the macro needs")
    p.add_argument("--force", action="store_true", help="allow overwriting outputs")
    p.set_defaults(func=cmd_plot_all)

    return parser


def main(argv: list[str]) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except ToolError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
