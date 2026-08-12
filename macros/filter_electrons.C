// filter_electrons.C
//
// Removes every event that contains at least one electron (|pdgId| == 11)
// in the pflow_pdgId branch, and writes the surviving events to a new
// .root file keeping *all* branches of the original tree.
//
// Usage (interpreted):
//   root -l 'filter_electrons.C("in.root", "out.root", "tree")'
//
// Usage (compiled, faster):
//   root -l 'filter_electrons.C+("in.root", "out.root", "tree")'
//
// An RDataFrame-based alternative is provided at the bottom of this file:
//   root -l
//   .L filter_electrons.C+
//   filter_electrons_rdf("in.root", "out.root", "tree")

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

#include <iostream>

void filter_electrons(const char *inFileName  = "input.root",
                      const char *outFileName = "output_noElectrons.root",
                      const char *treeName    = "tree")
{
   // ---- input ------------------------------------------------------------
   TFile *fin = TFile::Open(inFileName, "READ");
   if (!fin || fin->IsZombie()) {
      std::cerr << "ERROR: cannot open input file " << inFileName << std::endl;
      return;
   }

   TTree *tin = nullptr;
   fin->GetObject(treeName, tin);
   if (!tin) {
      std::cerr << "ERROR: tree '" << treeName << "' not found in "
                << inFileName << std::endl;
      fin->Close();
      return;
   }

   // Make sure every branch is read, so every branch gets copied.
   tin->SetBranchStatus("*", 1);

   const Long64_t nIn = tin->GetEntries();

   // ---- output -----------------------------------------------------------
   TFile *fout = TFile::Open(outFileName, "RECREATE");
   if (!fout || fout->IsZombie()) {
      std::cerr << "ERROR: cannot create output file " << outFileName << std::endl;
      fin->Close();
      return;
   }
   fout->cd();

   // ---- the cut ----------------------------------------------------------
   // Sum$(...) counts how many elements of the pflow_pdgId collection satisfy
   // the condition in a given event; requiring 0 keeps only electron-free
   // events. (This expression also behaves correctly if pflow_pdgId happens
   // to be a plain scalar rather than a vector/array.)
   const char *cut = "Sum$(pflow_pdgId == 11) == 0";

   TTree *tout = tin->CopyTree(cut);
   if (!tout) {
      std::cerr << "ERROR: CopyTree failed. Check that the branch "
                   "'pflow_pdgId' exists." << std::endl;
      fout->Close();
      fin->Close();
      return;
   }

   const Long64_t nOut = tout->GetEntries();

   tout->Write("", TObject::kOverwrite);
   fout->Close();
   fin->Close();

   std::cout << "Input events   : " << nIn << "\n"
             << "Output events  : " << nOut << "\n"
             << "Rejected       : " << (nIn - nOut) << "\n"
             << "Wrote          : " << outFileName << std::endl;
}

// ---------------------------------------------------------------------------
// RDataFrame version. Equivalent result, usually faster on large files, and
// it handles std::vector / RVec branches natively.
//
// NOTE: with implicit multi-threading enabled the output event order is not
// guaranteed to match the input; leave it disabled if event order matters.
// ---------------------------------------------------------------------------
#include <ROOT/RDataFrame.hxx>

void filter_electrons_rdf(const char *inFileName  = "input.root",
                          const char *outFileName = "output_noElectrons.root",
                          const char *treeName    = "tree",
                          bool useMT              = false)
{
   if (useMT) ROOT::EnableImplicitMT();

   ROOT::RDataFrame df(treeName, inFileName);

   auto nIn = df.Count();

   auto filtered = df.Filter("Sum(abs(pflow_pdgId) == 11) == 0", "no electron");

   // Passing no column list snapshots *all* columns of the input tree.
   auto snap = filtered.Snapshot(treeName, outFileName);

   std::cout << "Input events   : " << *nIn << "\n"
             << "Output events  : " << snap->Count().GetValue() << "\n"
             << "Wrote          : " << outFileName << std::endl;
}
