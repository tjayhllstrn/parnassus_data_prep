// Add pflow_beta / truth_beta branches to the evt_tree of already-prepped ROOT files, in place.
//
// Workaround for files produced before ParnassusTree had beta branches: instead of re-running
// hipo2tree_parnassus.C from the hipo, copy the values over from the sibling trees in the same file.
//   evt_tree.pflow_beta[npflow] <- reco_EventTree.beta[Nmax]
//   evt_tree.truth_beta[ntruth] <- mc_EventTree.truebeta[Nmax]
// This is valid because hipo2tree_parnassus.C fills all three trees once per event from the same
// particle vectors, so entry i and particle j line up across the trees. That alignment is checked
// for every entry (reco Nmax == npflow, mc Nmax == ntruth) before anything is written.
//
// Usage:
//   root -l -b -q 'macros/add_beta_branches.C("prepped_data")'                // every *.root in the directory
//   root -l -b -q 'macros/add_beta_branches.C("prepped_data/file.root")'      // a single file
//   root -l -b -q 'macros/add_beta_branches.C("prepped_data",false,true)'     // dry run: validate only
//
// Optional arguments:
//   recursive = false -> search subdirectories when true
//   dryRun    = false -> validate every file but write nothing when true
//
// Files whose evt_tree already has pflow_beta or truth_beta are skipped, so re-running is safe.
// Returns 0 only if no file failed validation or writing.

#include <TBranch.h>
#include <TFile.h>
#include <TList.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

const int kNmax = 100;  // array size used by EventTree/ParnassusTree in src/TreeManager.C

enum class Status { kAdded, kSkipped, kFailed };

bool ends_with(const std::string &value, const std::string &suffix) {
    if (value.size() < suffix.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void collect_root_files(const std::string &dirPath, bool recursive, std::vector<std::string> &files) {
    TSystemDirectory dir(dirPath.c_str(), dirPath.c_str());
    TList *entries = dir.GetListOfFiles();
    if (!entries) return;

    std::vector<std::string> local;
    TIter next(entries);
    while (TObject *obj = next()) {
        std::string name(obj->GetName());
        if (name == "." || name == "..") continue;
        std::string fullPath = dirPath.back() == '/' ? dirPath + name : dirPath + "/" + name;

        FileStat_t st;
        if (gSystem->GetPathInfo(fullPath.c_str(), st) != 0) continue;
        if (R_ISDIR(st.fMode)) {
            if (recursive) collect_root_files(fullPath, recursive, files);
            continue;
        }
        if (ends_with(name, ".root")) local.push_back(fullPath);
    }
    delete entries;

    std::sort(local.begin(), local.end());
    files.insert(files.end(), local.begin(), local.end());
}

TBranch *require_branch(TTree *tree, const char *name, const std::string &filePath) {
    TBranch *b = tree->GetBranch(name);
    if (!b) {
        std::cerr << "ERROR: " << tree->GetName() << " has no '" << name << "' branch in " << filePath << std::endl;
    }
    return b;
}

// Reads the per-entry particle counts through individual branches (rather than TTree::GetEntry) so
// that only the few branches we need are decompressed and no branch status is changed on evt_tree.
Status process_file(const std::string &filePath, bool dryRun) {
    std::unique_ptr<TFile> file(TFile::Open(filePath.c_str(), dryRun ? "READ" : "UPDATE"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Could not open ROOT file: " << filePath << std::endl;
        return Status::kFailed;
    }

    TTree *evt = dynamic_cast<TTree *>(file->Get("evt_tree"));
    TTree *reco = dynamic_cast<TTree *>(file->Get("reco_EventTree"));
    TTree *mc = dynamic_cast<TTree *>(file->Get("mc_EventTree"));
    if (!evt || !reco || !mc) {
        std::cerr << "ERROR: Missing" << (evt ? "" : " evt_tree") << (reco ? "" : " reco_EventTree")
                  << (mc ? "" : " mc_EventTree") << " in " << filePath << std::endl;
        return Status::kFailed;
    }

    if (evt->GetBranch("pflow_beta") || evt->GetBranch("truth_beta")) {
        std::cout << "SKIP: evt_tree already has pflow_beta/truth_beta: " << filePath << std::endl;
        return Status::kSkipped;
    }

    TBranch *bNpflow = require_branch(evt, "npflow", filePath);
    TBranch *bNtruth = require_branch(evt, "ntruth", filePath);
    TBranch *bRecoN = require_branch(reco, "Nmax", filePath);
    TBranch *bRecoBeta = require_branch(reco, "beta", filePath);
    TBranch *bMcN = require_branch(mc, "Nmax", filePath);
    TBranch *bMcBeta = require_branch(mc, "truebeta", filePath);
    if (!bNpflow || !bNtruth || !bRecoN || !bRecoBeta || !bMcN || !bMcBeta) return Status::kFailed;

    const Long64_t nEntries = evt->GetEntries();
    if (reco->GetEntries() != nEntries || mc->GetEntries() != nEntries) {
        std::cerr << "ERROR: Entry counts differ (evt_tree=" << nEntries << ", reco_EventTree=" << reco->GetEntries()
                  << ", mc_EventTree=" << mc->GetEntries() << ") in " << filePath << std::endl;
        return Status::kFailed;
    }

    int npflow = 0, ntruth = 0, recoN = 0, mcN = 0;
    double recoBeta[kNmax], mcBeta[kNmax];
    double pflow_beta[kNmax], truth_beta[kNmax];
    bNpflow->SetAddress(&npflow);
    bNtruth->SetAddress(&ntruth);
    bRecoN->SetAddress(&recoN);
    bRecoBeta->SetAddress(recoBeta);
    bMcN->SetAddress(&mcN);
    bMcBeta->SetAddress(mcBeta);

    // Validate every entry before touching the file.
    for (Long64_t i = 0; i < nEntries; ++i) {
        bNpflow->GetEntry(i);
        bNtruth->GetEntry(i);
        bRecoN->GetEntry(i);
        bMcN->GetEntry(i);
        if (recoN != npflow || mcN != ntruth || npflow > kNmax || ntruth > kNmax) {
            std::cerr << "ERROR: Particle counts do not line up at entry " << i << " (reco Nmax=" << recoN
                      << " vs npflow=" << npflow << ", mc Nmax=" << mcN << " vs ntruth=" << ntruth << ") in "
                      << filePath << std::endl;
            return Status::kFailed;
        }
    }

    if (dryRun) {
        std::cout << "OK (dry run): " << nEntries << " entries would be updated: " << filePath << std::endl;
        return Status::kAdded;
    }

    // Same leaflists as ParnassusTree in src/TreeManager.C, so the result matches freshly prepped files.
    TBranch *bPflowBeta = evt->Branch("pflow_beta", pflow_beta, "pflow_beta[npflow]/D");
    TBranch *bTruthBeta = evt->Branch("truth_beta", truth_beta, "truth_beta[ntruth]/D");

    for (Long64_t i = 0; i < nEntries; ++i) {
        bNpflow->GetEntry(i);
        bNtruth->GetEntry(i);
        bRecoN->GetEntry(i);
        bRecoBeta->GetEntry(i);
        bMcN->GetEntry(i);
        bMcBeta->GetEntry(i);
        std::copy(recoBeta, recoBeta + npflow, pflow_beta);
        std::copy(mcBeta, mcBeta + ntruth, truth_beta);
        bPflowBeta->Fill();
        bTruthBeta->Fill();
    }

    file->cd();
    if (evt->Write("", TObject::kOverwrite) <= 0) {
        std::cerr << "ERROR: Failed to write updated evt_tree to " << filePath << std::endl;
        return Status::kFailed;
    }
    // Detach our stack buffers before the tree goes away with the file.
    evt->ResetBranchAddresses();
    reco->ResetBranchAddresses();
    mc->ResetBranchAddresses();
    file->Close();

    std::cout << "ADDED: pflow_beta/truth_beta to " << nEntries << " entries: " << filePath << std::endl;
    return Status::kAdded;
}

}  // namespace

int add_beta_branches(const char *inputPath, bool recursive = false, bool dryRun = false) {
    const std::string inPath(inputPath ? inputPath : "");
    FileStat_t inStat;
    if (inPath.empty() || gSystem->GetPathInfo(inPath.c_str(), inStat) != 0) {
        std::cerr << "ERROR: Input path is empty or does not exist: " << inPath << std::endl;
        return 1;
    }

    std::vector<std::string> files;
    if (R_ISDIR(inStat.fMode)) {
        collect_root_files(inPath, recursive, files);
    } else {
        files.push_back(inPath);
    }
    if (files.empty()) {
        std::cerr << "ERROR: No .root files found in: " << inPath << std::endl;
        return 1;
    }

    int nAdded = 0, nSkipped = 0;
    std::vector<std::string> failed;
    for (const auto &path : files) {
        switch (process_file(path, dryRun)) {
            case Status::kAdded: ++nAdded; break;
            case Status::kSkipped: ++nSkipped; break;
            case Status::kFailed: failed.push_back(path); break;
        }
    }

    std::cout << "\nSummary: " << files.size() << " file(s) | " << (dryRun ? "would update " : "updated ") << nAdded
              << " | skipped (already had branches) " << nSkipped << " | failed " << failed.size() << std::endl;
    for (const auto &path : failed) std::cout << "  FAILED: " << path << std::endl;

    return failed.empty() ? 0 : 1;
}
