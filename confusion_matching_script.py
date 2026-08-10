import ROOT
import numpy as np

FILENAME, TREENAME = "prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0306_GEMC_512_bad_neutral.root", "evt_tree"
DR_MAX = 0.10                      # max sqrt(dtheta^2 + dphi^2) to call it a match

RECO_PID, RECO_ETA, RECO_PHI = "pflow_pdgId", "pflow_eta", "pflow_phi"
TRUE_PID, TRUE_ETA, TRUE_PHI = "truth_pdgId", "truth_eta", "truth_phi"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.2f")
ROOT.gStyle.SetPalette(ROOT.kMint)
ROOT.TColor.InvertPalette()         # light = low, dark = high

codes  = [22, 11, -11, 211, -211, 321, -321, 2212, -2212, 2112, 45]
labels = ["#gamma", "e^{-}", "e^{+}", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}",
          "p", "#bar{p}", "n", "d", "other"]
n   = len(labels)
idx = {c: i for i, c in enumerate(codes)}    # anything else -> "other"


# %%
def eta_to_theta(eta):
    """theta = 2*arctan(exp(-eta)), in (0, pi). Clip guards against overflow."""
    return 2.0 * np.arctan(np.exp(-np.clip(eta, -700.0, 700.0)))


def greedy_match(th_r, phi_r, th_t, phi_t, dr_max=DR_MAX):
    """Return list of (reco_index, truth_index) pairs, closest first, one-to-one."""
    if th_r.size == 0 or th_t.size == 0:
        return []

    dth  = th_r[:, None] - th_t[None, :]
    dphi = np.abs(phi_r[:, None] - phi_t[None, :]) % (2.0 * np.pi)
    dphi = np.minimum(dphi, 2.0 * np.pi - dphi)          # wrap into [0, pi]
    dr   = np.sqrt(dth**2 + dphi**2)

    order = np.argsort(dr, axis=None)
    pairs, used_r, used_t = [], set(), set()
    for flat in order:
        ir, it = np.unravel_index(flat, dr.shape)
        if dr[ir, it] > dr_max:
            break                                        # rest are all worse
        if ir in used_r or it in used_t:
            continue
        used_r.add(ir)
        used_t.add(it)
        pairs.append((ir, it))
    return pairs


# %%
f = ROOT.TFile.Open(FILENAME)
t = f.Get(TREENAME)

names  = [RECO_PID, RECO_ETA, RECO_PHI, TRUE_PID, TRUE_ETA, TRUE_PHI]
leaves = {name: t.GetLeaf(name) for name in names}
if any(l is None for l in leaves.values()):
    raise KeyError("missing branch(es): "
                   + ", ".join(k for k, v in leaves.items() if v is None))

def vals(name):
    leaf = leaves[name]
    return np.fromiter((leaf.GetValue(k) for k in range(leaf.GetLen())),
                       dtype=float, count=leaf.GetLen())

h = ROOT.TH2D("h", "CLAS12 Recon;Reco;Truth", n, 0, n, n, 0, n)
for i, lab in enumerate(labels):
    h.GetXaxis().SetBinLabel(i + 1, lab)
    h.GetYaxis().SetBinLabel(n - i, lab)     # reversed: #gamma at the top

n_reco = n_matched = 0
for i in range(t.GetEntries()):
    t.GetEntry(i)
    if i % 100000 ==0:
        print(f"iterating on entriy {i}")
    pid_r, phi_r = vals(RECO_PID), vals(RECO_PHI)
    pid_t, phi_t = vals(TRUE_PID), vals(TRUE_PHI)
    th_r = eta_to_theta(vals(RECO_ETA))
    th_t = eta_to_theta(vals(TRUE_ETA))

    n_reco += pid_r.size
    for ir, it in greedy_match(th_r, phi_r, th_t, phi_t):
        n_matched += 1
        ix = idx.get(int(pid_r[ir]), n - 1)
        iy = idx.get(int(pid_t[it]), n - 1)
        h.Fill(ix + 0.5, (n - 1 - iy) + 0.5)

print(f"matched {n_matched} / {n_reco} pflow particles "
      f"({100.0 * n_matched / max(n_reco, 1):.1f}%), dR < {DR_MAX}")


# %%
# normalize each truth row to 1
for iy in range(1, n + 1):
    tot = sum(h.GetBinContent(ix, iy) for ix in range(1, n + 1))
    if tot:
        for ix in range(1, n + 1):
            h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) / tot)

c = ROOT.TCanvas("c", "c", 700, 660)
c.SetLeftMargin(0.13)
c.SetBottomMargin(0.12)
h.SetMinimum(0.0)
h.SetMaximum(1.0)
h.SetMarkerSize(1.0)
h.GetXaxis().SetLabelSize(0.04)
h.GetYaxis().SetLabelSize(0.04)
h.Draw("COL TEXT")
c.Draw()                             # renders inline in the notebook
c.SaveAs("plots/clasdis_rga_fa18_inb_45nA_10604MeV-0306_GEMC_512_bad_neutral/confusion_matrix.png")
