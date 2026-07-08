// root -l -b -q 'macros/check_classID.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.root")'

void check_classID(const char* filename){
    //look at the truth_class branch of evt_tree and print the pid/class if the class variable is negative
    TFile* f = TFile::Open(filename,"READ");
    TTree* tree = (TTree*)f->Get("evt_tree");
    const int kNmax = 100;
    int ntruth;
    int nreco;
    int truth_class[kNmax];
    int reco_class[kNmax];
    int truth_pid[kNmax];
    int reco_pid[kNmax];
    tree->SetBranchAddress("ntruth",&ntruth);
    tree->SetBranchAddress("nreco",&nreco);
    tree->SetBranchAddress("truth_class",truth_class);
    tree->SetBranchAddress("reco_class",reco_class);
    tree->SetBranchAddress("truth_pid",truth_pid);
    tree->SetBranchAddress("reco_pid",reco_pid);
    int nentries = tree->GetEntries();
    for(int i = 0; i < nentries; i++){
        tree->GetEntry(i);
        for(int j = 0; j < ntruth; j++){
            if(truth_class[j] < 0){
                std::cout << "Truth Entry: " << i << " Particle: " << j << " PID: " << truth_pid[j] << " Class: " << truth_class[j] << std::endl;
            }
        }
        for(int j = 0; j < nreco; j++){
            if(reco_class[j] < 0){
                std::cout << "Reco Entry: " << i << " Particle: " << j << " PID: " << reco_pid[j] << " Class: " << reco_class[j] << std::endl;
            }
        }
    }
}