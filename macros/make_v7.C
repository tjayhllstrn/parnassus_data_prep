// Build a Parnassus v7 evt_tree directly from the reco_EventTree and mc_EventTree of processed per-hipo files.
//
// For processed files made before hipo2tree_parnassus.C wrote v7 itself (e.g. v5); new hipo runs need no
// conversion. Each entry of the two EventTrees is turned back into the std::vector<part> lists
// hipo2tree_parnassus.C filled them from, and written with ParnassusTree (src/TreeManager.C):
//   1. pflow_p / truth_p (momentum magnitude, GeV) next to pflow_pt / truth_pt (kept, e.g. for Ht).
//   2. truth_theta unfolded to [0, pi], atan2(truept, truepz). Old files' truetheta came from the folding
//      abs(atan(pt/pz)) version of Kinematics::th.
// The pflow e+/e- beta cut is applied here (it is not in old processed files): e+/e- failing
// CutManager::leptonBetaCut (0.8 < beta < 1.2) are removed. evt_tree is not read, so it does not matter
// whether the input has pflow_beta/truth_beta there.
// Unlike a fresh hipo run, old files keep their pt > 2.5 GeV reco cut (now p > 10.6 GeV) and lack truth
// particles above 179 deg (the 1 deg truth cut was applied to the folded angle).
//
// Scattered electrons (ParnassusTree sets pflow_scatE/truth_scatE from part.is_scattered_electron):
//   reco: hipo2tree_parnassus.C takes the highest-energy pid 11 particle and fills the reco DIS variables only if
//         it survives CutManager::filter_particles (otherwise Q2 = -999). So when reco_EventTree.Q2 != -999 it is
//         the highest-E e- of the stored list. pflow_scatE = 1 if it also survives the beta cut.
//   truth: e- with trueparentid 0 or -1, as in CLAS12Ana::load_mc_particles.
// eventNumber = firstEventNumber + output entry index, so it is unique across the merged files.
//
// Modes:
//   "train": keep only events with pflow_scatE == 1
//   "test" : keep every event, including those without a reco scattered electron (pflow_scatE == 0)
//
// Usage (compile with +; the output is opened with CREATE, so it is never overwritten):
//   root -l -b -q 'macros/make_v7.C+("/volatile/clas12/users/tjhellst/parnassus_data/v5/train_files","v7_train.root","train",-1,50000000)'
//   root -l -b -q 'macros/make_v7.C+("/volatile/clas12/users/tjhellst/parnassus_data/v5/test_files","v7_test.root","test",-1,5000000,200000000)'
//
// Arguments:
//   input                 directory (every *.root in it, LC_ALL=C order), text file with one path per line, or one .root file
//   outFile               output ROOT file: tree evt_tree, compression 101 (ZLIB level 1)
//   mode                  "train" or "test"
//   nFiles = -1           use only the first nFiles input files (-1: all)
//   nEvents = -1          stop after writing nEvents events (-1: all)
//   firstEventNumber = 0  eventNumber of the first output entry
//
// Returns 0 on success, 1 on error or if fewer than nEvents events were available (the file is still written).

#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/CutManager.C"
#include "../src/TreeManager.C"

namespace {

const int kNmax = 100;  // array size used by EventTree/ParnassusTree in src/TreeManager.C

// CutManager keeps its cuts protected (production code goes through filter_particles); expose only its e+/e-
// beta window.
class BetaCut : public CutManager {
   public:
    using CutManager::leptonBetaCut;
};

bool ends_with(const std::string &value, const std::string &suffix) {
    return value.size() >= suffix.size() && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool collect_inputs(const std::string &source, std::vector<std::string> &files) {
    FileStat_t st;
    if (source.empty() || gSystem->GetPathInfo(source.c_str(), st) != 0) {
        std::cerr << "ERROR: Input path is empty or does not exist: " << source << std::endl;
        return false;
    }
    if (R_ISDIR(st.fMode)) {
        TSystemDirectory dir(source.c_str(), source.c_str());
        if (TList *entries = dir.GetListOfFiles()) {
            TIter next(entries);
            while (TObject *obj = next()) {
                const std::string name(obj->GetName());
                if (ends_with(name, ".root")) files.push_back(source.back() == '/' ? source + name : source + "/" + name);
            }
            delete entries;
        }
        std::sort(files.begin(), files.end());
    } else if (ends_with(source, ".root")) {
        files.push_back(source);
    } else {
        std::ifstream in(source);
        std::string line;
        while (std::getline(in, line))
            if (!line.empty() && line[0] != '#') files.push_back(line);
    }
    if (files.empty()) std::cerr << "ERROR: No input .root files found in: " << source << std::endl;
    return !files.empty();
}

void enable(TTree *tree, const std::vector<const char *> &names) {
    tree->SetBranchStatus("*", 0);
    for (const char *name : names) tree->SetBranchStatus(name, 1);
    tree->SetCacheSize(64LL << 20);
}

}  // namespace

int make_v7(const char *input, const char *outFile, const char *mode, int nFiles = -1, Long64_t nEvents = -1,
            int firstEventNumber = 0) {
    const std::string modeStr(mode ? mode : "");
    if (modeStr != "train" && modeStr != "test") {
        std::cerr << "ERROR: mode must be \"train\" or \"test\", got \"" << modeStr << "\"" << std::endl;
        return 1;
    }
    const bool train = modeStr == "train";
    std::vector<std::string> files;
    if (!collect_inputs(input ? input : "", files)) return 1;
    if (nFiles >= 0 && nFiles < static_cast<int>(files.size())) files.resize(nFiles);

    std::unique_ptr<TFile> fout(TFile::Open(outFile, "CREATE"));
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: Could not create output file (does it already exist?): " << outFile << std::endl;
        return 1;
    }
    fout->SetCompressionSettings(101);
    fout->cd();
    ParnassusTree treeParnassus("evt_tree");
    treeParnassus.GetTree()->SetAutoSave(0);  // a single evt_tree cycle in the output
    BetaCut betaCut;
    EVENT recoEvent, mcEvent;  // FillTree takes them; the DIS branches are not written
    EVENT_INFO eventInfo;
    std::vector<part> vec_particles, vec_mcparticles;

    // Input buffers.
    int recoN = 0, mcN = 0;
    double Q2 = 0;
    double px[kNmax], py[kNmax], p[kNmax], phi[kNmax], theta[kNmax], vx[kNmax], vy[kNmax], vz[kNmax], beta[kNmax], E[kNmax];
    double truepx[kNmax], truepy[kNmax], truep[kNmax], truept[kNmax], truepz[kNmax], truephi[kNmax], truevx[kNmax], truevy[kNmax],
        truevz[kNmax], truebeta[kNmax];
    int pid[kNmax], truepid[kNmax], trueparentid[kNmax];
    auto recoPart = [&](int j) {
        part par;
        par.pid = pid[j];
        par.px = px[j];
        par.py = py[j];
        par.p = p[j];
        par.E = E[j];
        par.phi = phi[j];
        par.theta = theta[j];
        par.vx = vx[j];
        par.vy = vy[j];
        par.vz = vz[j];
        par.beta = beta[j];
        return par;
    };
    auto isLepton = [&](int j) { return std::abs(pid[j]) == 11; };

    TStopwatch watch;
    Long64_t nRead = 0, nOut = 0, nScatEIn = 0, nScatECut = 0, nNoElectronQ2 = 0, nLeptonsCut = 0;
    Long64_t nPflow = 0, nTruth = 0, nBackward = 0, nEvBackward = 0, nOutScatE = 0;
    std::vector<std::pair<std::string, Long64_t>> used;  // file, entries read from it
    bool done = false, failed = false;
    for (const auto &path : files) {
        if (done || failed) break;
        std::unique_ptr<TFile> fin(TFile::Open(path.c_str(), "READ"));
        TTree *reco = fin && !fin->IsZombie() ? dynamic_cast<TTree *>(fin->Get("reco_EventTree")) : nullptr;
        TTree *mc = reco ? dynamic_cast<TTree *>(fin->Get("mc_EventTree")) : nullptr;
        if (!reco || !mc || reco->GetEntries() != mc->GetEntries()) {
            std::cerr << "ERROR: Missing reco_EventTree/mc_EventTree or entry counts differ in " << path << std::endl;
            failed = true;
            break;
        }
        enable(reco, {"Nmax", "Q2", "px", "py", "p", "phi", "theta", "vx", "vy", "vz", "beta", "E", "pid"});
        enable(mc, {"Nmax", "truepx", "truepy", "truep", "truept", "truepz", "truephi", "truevx", "truevy", "truevz", "truebeta",
                    "truepid", "trueparentid"});
        reco->SetBranchAddress("Nmax", &recoN);
        reco->SetBranchAddress("Q2", &Q2);
        reco->SetBranchAddress("px", px);
        reco->SetBranchAddress("py", py);
        reco->SetBranchAddress("p", p);
        reco->SetBranchAddress("phi", phi);
        reco->SetBranchAddress("theta", theta);
        reco->SetBranchAddress("vx", vx);
        reco->SetBranchAddress("vy", vy);
        reco->SetBranchAddress("vz", vz);
        reco->SetBranchAddress("beta", beta);
        reco->SetBranchAddress("E", E);
        reco->SetBranchAddress("pid", pid);
        mc->SetBranchAddress("Nmax", &mcN);
        mc->SetBranchAddress("truepx", truepx);
        mc->SetBranchAddress("truepy", truepy);
        mc->SetBranchAddress("truep", truep);
        mc->SetBranchAddress("truept", truept);
        mc->SetBranchAddress("truepz", truepz);
        mc->SetBranchAddress("truephi", truephi);
        mc->SetBranchAddress("truevx", truevx);
        mc->SetBranchAddress("truevy", truevy);
        mc->SetBranchAddress("truevz", truevz);
        mc->SetBranchAddress("truebeta", truebeta);
        mc->SetBranchAddress("truepid", truepid);
        mc->SetBranchAddress("trueparentid", trueparentid);
        TBranch *bRecoN = reco->GetBranch("Nmax");
        TBranch *bMcN = mc->GetBranch("Nmax");

        const Long64_t n = reco->GetEntries();
        Long64_t i = 0;
        for (; i < n; ++i) {
            // Counters first, so an oversized entry cannot overrun the buffers.
            bRecoN->GetEntry(i);
            bMcN->GetEntry(i);
            if (recoN < 0 || recoN > kNmax || mcN < 0 || mcN > kNmax) {
                std::cerr << "ERROR: Nmax outside [0, " << kNmax << "] at entry " << i << " of " << path << std::endl;
                failed = true;
                break;
            }
            reco->GetEntry(i);
            ++nRead;

            int scatIdx = -1;
            if (Q2 != -999) {
                for (int j = 0; j < recoN; ++j)
                    if (pid[j] == 11 && (scatIdx < 0 || E[j] > E[scatIdx])) scatIdx = j;
                if (scatIdx < 0) ++nNoElectronQ2;
                else ++nScatEIn;
            }
            const bool scatE = scatIdx >= 0 && betaCut.leptonBetaCut(recoPart(scatIdx));
            if (scatIdx >= 0 && !scatE) ++nScatECut;
            if (train && !scatE) continue;

            vec_particles.clear();
            for (int j = 0; j < recoN; ++j) {
                part par = recoPart(j);
                if (isLepton(j) && !betaCut.leptonBetaCut(par)) {
                    ++nLeptonsCut;
                    continue;
                }
                if (j == scatIdx) par.is_scattered_electron = 1;
                vec_particles.push_back(par);
            }

            mc->GetEntry(i);
            vec_mcparticles.clear();
            bool backward = false;
            for (int j = 0; j < mcN; ++j) {
                part par;
                par.pid = truepid[j];
                par.truepid = truepid[j];
                par.truepx = truepx[j];
                par.truepy = truepy[j];
                par.truep = truep[j];
                par.truept = truept[j];
                par.truepz = truepz[j];
                par.truephi = truephi[j];
                par.truevx = truevx[j];
                par.truevy = truevy[j];
                par.truevz = truevz[j];
                par.truebeta = truebeta[j];
                par.trueparentid = trueparentid[j];
                if (truepid[j] == 11 && (trueparentid[j] == 0 || trueparentid[j] == -1)) par.is_scattered_electron = 1;
                if (truepz[j] < 0) {
                    backward = true;
                    ++nBackward;
                }
                vec_mcparticles.push_back(par);
            }
            if (backward) ++nEvBackward;

            eventInfo.evnum = firstEventNumber + static_cast<int>(nOut);
            treeParnassus.FillTree(vec_particles, vec_mcparticles, recoEvent, mcEvent, eventInfo);
            ++nOut;
            if (scatE) ++nOutScatE;
            nPflow += vec_particles.size();
            nTruth += vec_mcparticles.size();
            if (nOut % 5000000 == 0) {
                std::cout << "  written " << nOut << " (" << watch.RealTime() << " s)" << std::endl;
                watch.Continue();
            }
            if (nEvents > 0 && nOut >= nEvents) {
                ++i;
                done = true;
                break;
            }
        }
        reco->ResetBranchAddresses();
        mc->ResetBranchAddresses();
        used.push_back({path, i});
        std::cout << "  read " << i << " / " << n << " entries of " << path << "; written so far " << nOut << std::endl;
    }

    fout->cd();
    treeParnassus.Write();
    fout->Close();

    auto pct = [](Long64_t a, Long64_t b) { return b > 0 ? 100.0 * a / b : 0.0; };
    std::cout << "\nmake_v7 (" << modeStr << ")\n"
              << "Files used                      : " << used.size() << " of " << files.size() << "\n";
    for (const auto &u : used) std::cout << "  " << u.second << " entries from " << u.first << "\n";
    std::cout << "Events read                     : " << nRead << "\n"
              << "  with a reco scattered e-      : " << nScatEIn << " (" << pct(nScatEIn, nRead) << "%)\n"
              << "  scattered e- fails beta cut   : " << nScatECut << "\n"
              << "  Q2 set but no e- (unexpected) : " << nNoElectronQ2 << "\n"
              << "Events written                  : " << nOut << (nEvents > 0 ? " (requested " + std::to_string(nEvents) + ")" : "")
              << ", with pflow_scatE == 1: " << nOutScatE << " (" << pct(nOutScatE, nOut) << "%)\n"
              << "e+/e- removed by beta cut       : " << nLeptonsCut << "\n"
              << "pflow / truth particles         : " << nPflow << " / " << nTruth << "\n"
              << "truth particles with pz < 0     : " << nBackward << " (" << pct(nBackward, nTruth) << "%), in "
              << nEvBackward << " events (" << pct(nEvBackward, nOut) << "%)\n"
              << "eventNumber                     : " << firstEventNumber << " to " << firstEventNumber + nOut - 1 << "\n"
              << "Wrote                           : " << outFile << " (" << watch.RealTime() << " s)" << std::endl;

    if (failed) return 1;
    if (nEvents > 0 && nOut < nEvents) {
        std::cerr << "WARNING: only " << nOut << " of the requested " << nEvents << " events were available; add input files." << std::endl;
        return 1;
    }
    return 0;
}
