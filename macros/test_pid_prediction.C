double PION_MASS = 0.13957000;
double KAON_MASS = 0.49367700;
double PROTON_MASS = 0.93827200;
double ELECTRON_MASS = 0.00051099891;

double Pi_k_border = PION_MASS + (KAON_MASS - PION_MASS) / 2;
double k_p_border = KAON_MASS + (PROTON_MASS - KAON_MASS) / 2;


double calc_rest_mass(double beta, double p){

    if(beta >= 1.0 || beta <= 0.0){
        // std::cerr << "Warning: beta being clipped" << std::endl;
        beta = std::max(0.0, std::min(beta, 1.0));
    }

    double gamma = 1.0/sqrt(1.0-beta*beta);
    double rest_mass = p / (beta*gamma);

    return rest_mass;

}

int predict_pid(double rest_mass, int charge, int classID){
    if(classID == 0) { //charged hadron
        if(charge == 0) return -999; //unknown charge for charged hadron
        if(charge == 1){
            if(rest_mass < Pi_k_border) return 211; //pi+
            else if(rest_mass < k_p_border) return 321; //K+
            else return 2212; //p
        } else if(charge == -1){
            if(rest_mass < Pi_k_border) return -211; //pi-
            else if(rest_mass < k_p_border) return -321; //K-
            else return -2212; //anti-p
        }
    } else if (classID == 3) { //neutral hadron
        if(charge != 0) return -999; //unknown charge for neutral hadron
        if(rest_mass < Pi_k_border) return 111; //pi0
        else if(rest_mass < k_p_border) return 130; //K0L
        else return 2112; //neutron

    } else if(classID == 1) { //electron
        if (charge == -1) return 11;
        else if (charge == 1) return -11;
        else return -999; //unknown charge for electron

    } else if(classID == 2) { //muon
        if (charge == -1) return 13;
        else if (charge == 1) return -13;
        else return -999; //unknown charge for muon
        return 13;
    } else if(classID==4){
        return 22; //photon
    }
    else {
        return -999; //unknown
    }

}

double to_p(double pt, double eta){
    return pt*cosh(eta);
}

double getPdgMass(int pid)  {

    switch (pid) {

    case 22 :
      return 0.;
    case 11 :
      return 0.00051099891;
    case 211:
      return 0.13957000;
    case -211:
      return 0.13957000;
    case 2212:
      return 0.93827200;
    case 2112:
      return 0.93956500;
    case 321:
      return 0.49367700;
    case -321:
      return 0.49367700;
     case -11 :
      return 0.00051099891;
    case -2212:
      return 0.93827200;
    case 45:
      return 1.875612;
    default :
      return 0;
    }
}


// ----------------------------------------------------------------------------
//  Confusion-matrix category helpers
//
//  Both the truth PID and the predicted PID get mapped onto the same small set
//  of labelled categories so we can build a square "confusion matrix" (truth on
//  one axis, prediction on the other). The last category ("other") absorbs any
//  PID not explicitly listed (e.g. deuterons), and predict_pid's -999
//  "no prediction" code gets its own "none" category.
// ----------------------------------------------------------------------------
std::vector<int>         CAT_PIDS   = { 22,      11,     -11,   13,       -13,      211,      -211,     111,       321,     -321,    130,          2212, -2212,      2112, -999};
std::vector<std::string> CAT_LABELS = {"#gamma","e^{-}","e^{+}","#mu^{-}","#mu^{+}","#pi^{+}","#pi^{-}","#pi^{0}","K^{+}","K^{-}","K^{0}_{L}","p","#bar{p}","n", "none", "other"};

// index of a PID within the confusion-matrix axis (last bin = "other")
int pid_to_cat(int pid){
    for(size_t k = 0; k < CAT_PIDS.size(); ++k){
        if(CAT_PIDS[k] == pid) return (int)k;
    }
    return (int)CAT_PIDS.size(); // "other"
}


// ----------------------------------------------------------------------------
//  Run the full prediction test over one particle collection (prefix = "pflow"
//  or "truth"). Predicts the PID of every particle from its beta / pt / eta /
//  charge / classID, compares against the stored PID, and writes out plots +
//  statistics.
// ----------------------------------------------------------------------------
void analyze_set(TTree* tree, const char* prefix, const char* outbase, int n_entries=200){

    const int kNmax = 100;
    double beta[kNmax], pt[kNmax], eta[kNmax];
    int    charge[kNmax], cls[kNmax], pdgId[kNmax], n;

    tree->SetBranchAddress(Form("%s_beta",   prefix), beta);
    tree->SetBranchAddress(Form("%s_pt",     prefix), pt);
    tree->SetBranchAddress(Form("%s_eta",    prefix), eta);
    tree->SetBranchAddress(Form("%s_charge", prefix), charge);
    tree->SetBranchAddress(Form("%s_class",  prefix), cls);
    tree->SetBranchAddress(Form("%s_pdgId",  prefix), pdgId);
    tree->SetBranchAddress(Form("n%s",       prefix), &n);

    const int nCat = (int)CAT_LABELS.size();

    // --- confusion matrix: truth category (x) vs predicted category (y) ------
    TH2D* h_conf = new TH2D(Form("h_conf_%s", prefix),
        Form("PID confusion (%s);%s PID;predicted PID", prefix, prefix),
        nCat, 0, nCat, nCat, 0, nCat);
    for(int b = 0; b < nCat; ++b){
        h_conf->GetXaxis()->SetBinLabel(b+1, CAT_LABELS[b].c_str());
        h_conf->GetYaxis()->SetBinLabel(b+1, CAT_LABELS[b].c_str());
    }

    // --- outcome tally: correct / wrong / no-prediction ----------------------
    TH1D* h_outcome = new TH1D(Form("h_outcome_%s", prefix),
        Form("Prediction outcome (%s);%s PID;particles", prefix, prefix), 3, 0, 3);
    h_outcome->GetXaxis()->SetBinLabel(1, "correct");
    h_outcome->GetXaxis()->SetBinLabel(2, "wrong");
    h_outcome->GetXaxis()->SetBinLabel(3, "no prediction");

    // --- rest-mass comparison ------------------------------------------------
    // percent difference (predicted - true)/true * 100, overall...
    TH1D* h_massdiff = new TH1D(Form("h_massdiff_%s", prefix),
        Form("Rest-mass %% difference (%s);100 #times (m_{pred} - m_{%s})/m_{%s};particles", prefix, prefix, prefix),
        200, -100, 100);
    // ...and split by true species (only species with a defined non-zero mass)
    TH2D* h_massdiff_sp = new TH2D(Form("h_massdiff_sp_%s", prefix),
        Form("Rest-mass %% difference by species (%s);%s PID;100 #times (m_{pred} - m_{%s})/m_{%s}", prefix, prefix, prefix, prefix),
        nCat, 0, nCat, 200, -100, 100);
    for(int b = 0; b < nCat; ++b)
        h_massdiff_sp->GetXaxis()->SetBinLabel(b+1, CAT_LABELS[b].c_str());
    // predicted vs true rest mass scatter
    TH2D* h_mass_scatter = new TH2D(Form("h_mass_scatter_%s", prefix),
        Form("Predicted vs %s rest mass (%s);m_{%s} [GeV];m_{pred} [GeV]", prefix, prefix, prefix),
        150, 0, 1.5, 150, 0, 1.5);
    // fraction of particles whose predicted mass is within 10%% of the truth,
    // per species (TProfile mean of a 0/1 flag = the fraction)
    TProfile* p_within10 = new TProfile(Form("p_within10_%s", prefix),
        Form("Fraction of predicted rest mass within 10%% of %s;%s PID;fraction within 10%%", prefix, prefix),
        nCat, 0, nCat, 0, 1);
    for(int b = 0; b < nCat; ++b)
        p_within10->GetXaxis()->SetBinLabel(b+1, CAT_LABELS[b].c_str());

    // --- counters ------------------------------------------------------------
    long nTotal = 0, nCorrect = 0, nWrong = 0, nNoPred = 0;
    long nWithin10 = 0, nMassCompared = 0;
    // per-class breakdown: class -> {correct, wrong, none}
    std::map<int, std::array<long,3>> byClass;

    long nEvents = tree->GetEntries();
    if (n_entries > 0 && n_entries < nEvents) nEvents = n_entries;
    for(long i = 0; i < nEvents; ++i){
        tree->GetEntry(i);
        for(int j = 0; j < n; ++j){

            double p_mom             = to_p(pt[j], eta[j]);
            double predicted_rest_m  = calc_rest_mass(beta[j], p_mom);
            int    predicted_pid     = predict_pid(predicted_rest_m, charge[j], cls[j]);

            int    true_pid          = pdgId[j];
            double true_rest_m       = getPdgMass(true_pid);

            ++nTotal;

            // ---- outcome ----
            int outcome; // 0 correct, 1 wrong, 2 none
            if(predicted_pid == -999)          outcome = 2;
            else if(predicted_pid == true_pid) outcome = 0;
            else                               outcome = 1;

            if(outcome == 0){ ++nCorrect; h_outcome->Fill(0.5); }
            else if(outcome == 1){ ++nWrong; h_outcome->Fill(1.5); }
            else { ++nNoPred; h_outcome->Fill(2.5); }

            byClass[cls[j]][outcome]++;

            // ---- confusion matrix ----
            h_conf->Fill(pid_to_cat(true_pid) + 0.5, pid_to_cat(predicted_pid) + 0.5);

            // ---- rest-mass comparison (needs a valid truth mass reference) ----
            if(true_rest_m > 0.0){
                double pct = (predicted_rest_m - true_rest_m) / true_rest_m * 100.0;
                h_massdiff->Fill(pct);
                h_massdiff_sp->Fill(pid_to_cat(true_pid) + 0.5, pct);
                h_mass_scatter->Fill(true_rest_m, predicted_rest_m);

                bool within10 = std::fabs(pct) <= 10.0;
                p_within10->Fill(pid_to_cat(true_pid) + 0.5, within10 ? 1.0 : 0.0);
                ++nMassCompared;
                if(within10) ++nWithin10;
            }
        }
    }

    // ------------------------------------------------------------------------
    //  Draw & save
    // ------------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".0f");

    // output directory for this collection: plots/<input-stem>/<prefix>/
    TString dir = TString(outbase) + "/" + prefix;
    gSystem->mkdir(dir, true);

    TCanvas* c1 = new TCanvas(Form("c_conf_%s", prefix), "confusion", 900, 800);
    c1->SetRightMargin(0.13);
    c1->SetBottomMargin(0.13);
    h_conf->Draw("COLZ");
    c1->SaveAs(Form("%s/confusion.png", dir.Data()));

    // row-normalized confusion (purity per true species): each true-PID column
    // scaled to 100%, so the diagonal reads as "% of this species predicted right"
    TH2D* h_conf_norm = (TH2D*)h_conf->Clone(Form("h_conf_norm_%s", prefix));
    h_conf_norm->SetTitle(Form("PID confusion, column-normalized to %% (%s);%s PID;predicted PID", prefix, prefix));
    for(int bx = 1; bx <= nCat; ++bx){
        double col = 0;
        for(int by = 1; by <= nCat; ++by) col += h_conf->GetBinContent(bx, by);
        if(col <= 0) continue;
        for(int by = 1; by <= nCat; ++by)
            h_conf_norm->SetBinContent(bx, by, 100.0 * h_conf->GetBinContent(bx, by) / col);
    }
    TCanvas* c1b = new TCanvas(Form("c_confnorm_%s", prefix), "confusion norm", 900, 800);
    c1b->SetRightMargin(0.13);
    c1b->SetBottomMargin(0.13);
    h_conf_norm->SetMarkerSize(0.8);
    h_conf_norm->Draw("COLZ TEXT");
    c1b->SaveAs(Form("%s/confusion_norm.png", dir.Data()));

    TCanvas* c2 = new TCanvas(Form("c_outcome_%s", prefix), "outcome", 700, 600);
    h_outcome->SetFillColor(kAzure+1);
    h_outcome->SetBarWidth(0.6);
    h_outcome->SetBarOffset(0.2);
    h_outcome->Draw("BAR");
    c2->SaveAs(Form("%s/outcome.png", dir.Data()));

    TCanvas* c3 = new TCanvas(Form("c_massdiff_%s", prefix), "mass diff", 800, 600);
    h_massdiff->SetFillColor(kOrange+7);
    h_massdiff->Draw("HIST");
    c3->SaveAs(Form("%s/massdiff.png", dir.Data()));

    TCanvas* c4 = new TCanvas(Form("c_massdiff_sp_%s", prefix), "mass diff by species", 900, 700);
    c4->SetRightMargin(0.13);
    h_massdiff_sp->Draw("COLZ");
    c4->SaveAs(Form("%s/massdiff_by_species.png", dir.Data()));

    TCanvas* c5 = new TCanvas(Form("c_within10_%s", prefix), "within 10pct", 800, 600);
    p_within10->SetMarkerStyle(20);
    p_within10->SetLineColor(kGreen+2);
    p_within10->SetMarkerColor(kGreen+2);
    p_within10->GetYaxis()->SetRangeUser(0, 1.05);
    p_within10->Draw("P");
    c5->SaveAs(Form("%s/within10_by_species.png", dir.Data()));

    TCanvas* c6 = new TCanvas(Form("c_mass_scatter_%s", prefix), "mass scatter", 800, 700);
    c6->SetRightMargin(0.13);
    h_mass_scatter->Draw("COLZ");
    // reference y = x line
    TLine* diag = new TLine(0, 0, 1.5, 1.5);
    diag->SetLineColor(kRed);
    diag->SetLineStyle(2);
    diag->Draw("SAME");
    c6->SaveAs(Form("%s/mass_scatter.png", dir.Data()));

    // persist the histograms for later inspection
    TFile* fout = new TFile(Form("%s/hists.root", dir.Data()), "RECREATE");
    h_conf->Write(); h_conf_norm->Write(); h_outcome->Write();
    h_massdiff->Write(); h_massdiff_sp->Write(); h_mass_scatter->Write();
    p_within10->Write();
    fout->Close();

    // ------------------------------------------------------------------------
    //  Print statistics
    // ------------------------------------------------------------------------
    auto pct = [](long a, long b){ return b > 0 ? 100.0 * a / b : 0.0; };

    std::cout << "\n=================================================================\n";
    std::cout << "  PID prediction test  --  '" << prefix << "' particles\n";
    std::cout << "=================================================================\n";
    std::cout << "  particles processed : " << nTotal << "\n";
    std::cout << "  correct             : " << nCorrect << "  (" << pct(nCorrect, nTotal) << "%)\n";
    std::cout << "  wrong               : " << nWrong   << "  (" << pct(nWrong,   nTotal) << "%)\n";
    std::cout << "  no prediction       : " << nNoPred  << "  (" << pct(nNoPred,  nTotal) << "%)\n";
    long nAttempted = nCorrect + nWrong;
    std::cout << "  accuracy (of those actually predicted) : "
              << pct(nCorrect, nAttempted) << "%\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "  rest-mass comparison (species with a known mass):\n";
    std::cout << "    particles compared        : " << nMassCompared << "\n";
    std::cout << "    predicted mass within 10% : " << nWithin10
              << "  (" << pct(nWithin10, nMassCompared) << "%)\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "  breakdown by classID (0=chg had,1=e,2=mu,3=neu had,4=photon,-1=unk):\n";
    std::cout << "    class    correct     wrong      none     total\n";
    for(auto& kv : byClass){
        long c = kv.second[0], w = kv.second[1], nn = kv.second[2];
        long t = c + w + nn;
        std::cout << "    " << std::setw(5) << kv.first
                  << std::setw(11) << c
                  << std::setw(10) << w
                  << std::setw(10) << nn
                  << std::setw(10) << t << "\n";
    }
    std::cout << "=================================================================\n";
    std::cout << "  plots written to " << dir << "/*.png\n";
    std::cout << "  histograms saved to " << dir << "/hists.root\n\n";

    // detach branch addresses (locals go out of scope after this function)
    tree->ResetBranchAddresses();
}


void test_pid_prediction(const char* input_file, int n_entries=200){
    //load input ROOT file and extract evt_tree
    TFile* f = TFile::Open(input_file);
    if(!f || f->IsZombie()){
        std::cerr << "Could not open input file: " << input_file << std::endl;
        return;
    }
    TTree* tree = (TTree*)f->Get("evt_tree");
    if(!tree){
        std::cerr << "Could not find 'evt_tree' in " << input_file << std::endl;
        return;
    }

    // build an output base directory (plots/<input-stem>) from the input file;
    // each collection then gets its own <prefix> subdirectory inside analyze_set
    gSystem->mkdir("plots", true);
    TString stem = gSystem->BaseName(input_file);
    if(stem.EndsWith(".root")) stem.Remove(stem.Length() - 5);
    TString outbase = TString("plots/test_pid_prediction/") + stem;

    // Predict PID from beta/pt/eta/charge/classID and compare to the stored PID.
    // Run it on the reconstructed particles (pflow) -- the primary test of the
    // prediction -- and on the MC truth particles as an ideal-input cross-check.
    analyze_set(tree, "pflow", outbase.Data(), n_entries);
    analyze_set(tree, "truth", outbase.Data(), n_entries);

    f->Close();
}
