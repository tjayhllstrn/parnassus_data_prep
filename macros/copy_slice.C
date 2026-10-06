// copy_slice.C  —  root -l -b -q copy_slice.C
void copy_slice(const char* in  = "input.root",
                const char* out = "output.root",
                Long64_t first  = 28000000,
                Long64_t n      = 1000000)
{
   TFile* fin = TFile::Open(in);
   TTree* tin = (TTree*)fin->Get("evt_tree");

   TFile* fout = TFile::Open(out, "RECREATE");   // must be open before CopyTree
   TTree* tout = tin->CopyTree("", "", n, first);

   printf("copied %lld entries starting at %lld\n", tout->GetEntries(), first);
   tout->Write();
   fout->Close();
   fin->Close();
}
