// Rewrite the evt_tree "eventNumber" branch so that every entry gets a unique,
// sequential value. The raw eventNumber copied from each input hipo file is only
// unique within that file, so a merged evt_tree ends up with duplicate event
// numbers. This macro rebuilds the tree with eventNumber = entry index.
//
// Usage:
//   root -l -b -q 'macros/reassign_event_number.C("in.root","out.root")'
//
// The output is written to a separate file (outFile); the caller is responsible
// for moving it over the original if an in-place update is desired.

#include <TBranch.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <memory>
#include <string>

int reassign_event_number(const char *inFile, const char *outFile) {
    if (!inFile || !outFile) {
        std::cerr << "ERROR: inFile and outFile are required." << std::endl;
        return 1;
    }

    const std::string inPath(inFile);
    const std::string outPath(outFile);
    if (inPath.empty() || outPath.empty()) {
        std::cerr << "ERROR: inFile and outFile cannot be empty." << std::endl;
        return 1;
    }
    if (inPath == outPath) {
        std::cerr << "ERROR: inFile and outFile must differ." << std::endl;
        return 1;
    }

    std::unique_ptr<TFile> in(TFile::Open(inPath.c_str(), "READ"));
    if (!in || in->IsZombie()) {
        std::cerr << "ERROR: Could not open input ROOT file: " << inPath << std::endl;
        return 1;
    }

    TTree *told = dynamic_cast<TTree *>(in->Get("evt_tree"));
    if (!told) {
        std::cerr << "ERROR: Missing evt_tree in file: " << inPath << std::endl;
        return 1;
    }
    if (!told->GetBranch("eventNumber")) {
        std::cerr << "ERROR: evt_tree has no eventNumber branch in file: " << inPath << std::endl;
        return 1;
    }

    // Drop the existing eventNumber branch from the clone; we replace it below.
    told->SetBranchStatus("*", 1);
    told->SetBranchStatus("eventNumber", 0);

    TFile out(outPath.c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: Could not create output ROOT file: " << outPath << std::endl;
        return 1;
    }

    out.cd();
    TTree *tnew = told->CloneTree(0);
    if (!tnew) {
        std::cerr << "ERROR: CloneTree returned null." << std::endl;
        return 1;
    }

    int eventNumber = 0;
    tnew->Branch("eventNumber", &eventNumber, "eventNumber/I");

    const Long64_t nEntries = told->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        told->GetEntry(i);
        eventNumber = static_cast<int>(i);
        tnew->Fill();
    }

    tnew->SetName("evt_tree");
    tnew->Write("evt_tree", TObject::kOverwrite);
    out.Close();

    std::cout << "Reassigned eventNumber for " << nEntries
              << " entries into: " << outPath << std::endl;

    return 0;
}
