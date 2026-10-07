// Apply a particle-level beta cut to e+/e- in the evt_tree of a prepped Parnassus file, writing a new file.
//
// A pflow particle is removed if |pflow_pdgId| == 11 and pflow_beta is outside [betaMin, betaMax].
// After that, an event with no pflow electron (pflow_pdgId == 11) left is dropped entirely, unless
// dropNoElectron = false: then the event is kept and its pflow_scatE is set to 0 (no reconstructed
// scattered electron), so the evaluation chain can filter it out later (e.g. test files that keep
// non-electron events). Every other branch (truth_*, ntruth, truth_scatE, eventNumber) is copied unchanged.
//
// Motivation: REC::Particle beta for e+/e- has rare timing failures (beta from -1458 to +1544 in v6),
// while a real electron at the momenta in these files (p > 250 MeV) has beta = 1 to within 1e-5,
// so the cut removes no physical electrons. Photons already get the same window in CutManager::photonBetaCut.
//
// Usage (compile with + for speed; the input is opened READ, the output with CREATE so it is never overwritten):
//   root -l -b -q 'macros/apply_lepton_beta_cut.C+("in.root","out.root")'
//   root -l -b -q 'macros/apply_lepton_beta_cut.C+("in.root","out.root",0.9,1.1)'
//   root -l -b -q 'macros/apply_lepton_beta_cut.C+("in.root","out.root",0.9,1.1,false)'   // keep no-electron events
//
// Returns 0 on success.

#include <TFile.h>
#include <TTree.h>

#include <cstdlib>
#include <iostream>
#include <memory>

namespace {

const int kNmax = 100;  // array size used by ParnassusTree in src/TreeManager.C

}  // namespace

int apply_lepton_beta_cut(const char *inFileName, const char *outFileName, bool dropNoElectron = true, double betaMin = 0.9,
                          double betaMax = 1.1, const char *treeName = "evt_tree") {
    std::unique_ptr<TFile> fin(TFile::Open(inFileName, "READ"));
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Could not open input file: " << inFileName << std::endl;
        return 1;
    }
    TTree *tin = dynamic_cast<TTree *>(fin->Get(treeName));
    if (!tin) {
        std::cerr << "ERROR: No tree '" << treeName << "' in " << inFileName << std::endl;
        return 1;
    }

    const char *arrays[] = {"pflow_pt", "pflow_phi", "pflow_theta", "pflow_vx", "pflow_vy", "pflow_vz", "pflow_beta"};
    for (const char *name : {"npflow", "pflow_pt", "pflow_phi", "pflow_theta", "pflow_vx", "pflow_vy", "pflow_vz",
                             "pflow_pdgId", "pflow_beta", "pflow_scatE"}) {
        if (!tin->GetBranch(name)) {
            std::cerr << "ERROR: " << treeName << " has no '" << name << "' branch in " << inFileName << std::endl;
            return 1;
        }
    }
    if (tin->GetMaximum("npflow") > kNmax) {
        std::cerr << "ERROR: npflow exceeds kNmax=" << kNmax << " in " << inFileName << std::endl;
        return 1;
    }

    // Addresses are set before CloneTree so the clone shares these buffers; all other branches are
    // copied through ROOT's own buffers.
    int npflow = 0, pflowScatE = 0;
    int pdgId[kNmax];
    double buf[7][kNmax];
    tin->SetBranchStatus("*", 1);
    tin->SetBranchAddress("npflow", &npflow);
    tin->SetBranchAddress("pflow_scatE", &pflowScatE);
    tin->SetBranchAddress("pflow_pdgId", pdgId);
    for (int a = 0; a < 7; ++a) tin->SetBranchAddress(arrays[a], buf[a]);

    std::unique_ptr<TFile> fout(TFile::Open(outFileName, "CREATE"));
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: Could not create output file (does it already exist?): " << outFileName << std::endl;
        return 1;
    }
    fout->SetCompressionSettings(fin->GetCompressionSettings());
    fout->cd();
    TTree *tout = tin->CloneTree(0);

    double *beta = buf[6];
    const Long64_t nIn = tin->GetEntries();
    Long64_t nOut = 0, nEvAffected = 0, nRemovedEm = 0, nRemovedEp = 0, nNeg = 0, nLow = 0, nHigh = 0;
    Long64_t nPflowIn = 0, nPflowOut = 0, nNoElectron = 0, nScatEFlipped = 0;
    for (Long64_t i = 0; i < nIn; ++i) {
        tin->GetEntry(i);
        nPflowIn += npflow;

        int n = 0;
        bool hasElectron = false;
        for (int j = 0; j < npflow; ++j) {
            const bool isLepton = std::abs(pdgId[j]) == 11;
            if (isLepton && (beta[j] < betaMin || beta[j] > betaMax)) {
                (pdgId[j] == 11 ? nRemovedEm : nRemovedEp)++;
                if (beta[j] < 0) ++nNeg;
                else if (beta[j] < betaMin) ++nLow;
                else ++nHigh;
                continue;
            }
            if (pdgId[j] == 11) hasElectron = true;
            if (n != j) {
                pdgId[n] = pdgId[j];
                for (int a = 0; a < 7; ++a) buf[a][n] = buf[a][j];
            }
            ++n;
        }
        if (n != npflow) ++nEvAffected;
        npflow = n;
        if (!hasElectron) {
            ++nNoElectron;
            if (dropNoElectron) continue;
            if (pflowScatE != 0) ++nScatEFlipped;
            pflowScatE = 0;
        }

        tout->Fill();
        ++nOut;
        nPflowOut += npflow;
        if ((i + 1) % 5000000 == 0) std::cout << "  processed " << (i + 1) << " / " << nIn << std::endl;
    }

    fout->cd();
    if (tout->Write("", TObject::kOverwrite) <= 0) {
        std::cerr << "ERROR: Failed to write " << treeName << " to " << outFileName << std::endl;
        return 1;
    }
    fout->Close();

    std::cout << "\nCut: remove e+/e- with beta outside [" << betaMin << ", " << betaMax << "]; "
              << (dropNoElectron ? "drop" : "keep (pflow_scatE -> 0)") << " events with no pflow e- left\n"
              << "Input events             : " << nIn << "\n"
              << "Events with e+- removed  : " << nEvAffected << "\n"
              << "Events with no e- left   : " << nNoElectron << "\n"
              << "Events dropped           : " << (nIn - nOut) << "\n"
              << "pflow_scatE set 1 -> 0   : " << nScatEFlipped << "\n"
              << "Output events            : " << nOut << "\n"
              << "e+- removed              : " << (nRemovedEm + nRemovedEp) << " (e- " << nRemovedEm << ", e+ "
              << nRemovedEp << "; beta<0 " << nNeg << ", 0<=beta<" << betaMin << " " << nLow << ", beta>" << betaMax
              << " " << nHigh << ")\n"
              << "pflow particles in / out : " << nPflowIn << " / " << nPflowOut << "\n"
              << "Wrote                    : " << outFileName << std::endl;
    return 0;
}
