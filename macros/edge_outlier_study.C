#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include <filesystem>
#include <set>
#include <string>

// Write one flat tree entry per reco particle, BEFORE any cut, with the raw
// variables behind the calorimeter-edge, v_z / v_T outlier and chi2pid cuts,
// the calorimeter hit positions, the matched MC truth and three pass flags, so
// edge_outlier_study.ipynb can try other thresholds without rerunning the hipo
// step. A calorimeter layer with no hit has <layer>_e == 0 and x/y/z/lv/lw 0.
//
//   pass_prod   : kept by CutManager::filter_particles (the production cuts)
//   pass_loose  : kept by the production cuts with chi2pid, the v_z / v_T
//                 outlier cuts and the e- / gamma calorimeter edge cuts switched off
//   pass_nobeta : kept by the production cuts with the e+- / gamma beta outlier
//                 cut (0 < beta < 2) switched off
//
// Truth (truepid, truep, trueE) is the MC particle closest in angle
// (CLAS12Ana::match_mc_to_reco: dR = sqrt((dtheta/2deg)^2 + (dphi/4deg)^2) < 1);
// it is -999 for unmatched particles and for real data.
//
// root -l -b -q 'macros/edge_outlier_study.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","plots/edge_outlier_study/study.root",-1)'

int edge_outlier_study(const char* hipoFile  = "",
                       const char* outFile   = "plots/edge_outlier_study/study.root",
                       const int   maxEvents = 100000){

    bool is_real_data = std::string(hipoFile).find("reco") != std::string::npos;

    std::string outDir = std::filesystem::path(outFile).parent_path().string();
    if (!outDir.empty()) gSystem->mkdir(outDir.c_str(), true);

    TFile* fOut = new TFile(outFile, "RECREATE");
    TTree* tree = new TTree("particles", "reco particles before cuts");

    // The branches read straight from q, which is set to each particle before Fill()
    part q;
    int  run, evnum;
    bool pass_prod, pass_loose, pass_nobeta;
    tree->Branch("run", &run, "run/I");
    tree->Branch("evnum", &evnum, "evnum/I");
    tree->Branch("pid", &q.pid, "pid/I");
    tree->Branch("sector", &q.sector, "sector/I");
    tree->Branch("truepid", &q.truepid, "truepid/I");
    std::vector<std::pair<const char*, double*>> doubles = {
        {"p", &q.p}, {"E", &q.E}, {"theta", &q.theta}, {"phi", &q.phi}, {"beta", &q.beta}, {"chi2", &q.chi2},
        {"vx", &q.vx}, {"vy", &q.vy}, {"vz", &q.vz},
        {"pcal_e", &q.pcal_e}, {"ecin_e", &q.ecin_e}, {"ecout_e", &q.ecout_e},
        {"pcal_lv", &q.pcal_lv}, {"pcal_lw", &q.pcal_lw},
        {"ecin_lv", &q.ecin_lv}, {"ecin_lw", &q.ecin_lw},
        {"ecout_lv", &q.ecout_lv}, {"ecout_lw", &q.ecout_lw},
        {"pcal_x", &q.pcal_x}, {"pcal_y", &q.pcal_y}, {"pcal_z", &q.pcal_z},
        {"ecin_x", &q.ecin_x}, {"ecin_y", &q.ecin_y}, {"ecin_z", &q.ecin_z},
        {"ecout_x", &q.ecout_x}, {"ecout_y", &q.ecout_y}, {"ecout_z", &q.ecout_z},
        {"truep", &q.truep}, {"trueE", &q.trueE}};
    for (auto& branch : doubles)
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    tree->Branch("pass_prod", &pass_prod, "pass_prod/O");
    tree->Branch("pass_loose", &pass_loose, "pass_loose/O");
    tree->Branch("pass_nobeta", &pass_nobeta, "pass_nobeta/O");

    // Configure CLAS12 Reader and HipoChain
    // -------------------------------------
    clas12root::HipoChain _chain;
    _chain.Add(hipoFile);
    auto _config_c12 = _chain.GetC12Reader();
    auto& _c12 = _chain.C12ref();

    // Add Analysis Objects
    // -------------------------------------
    CutManager _cm = CutManager();
    CLAS12Ana clas12ana = CLAS12Ana(_c12);
    clas12ana.set_run_config(_c12); // must precede get_event_info

    std::vector<part> vec_particles;
    std::vector<part> vec_mcparticles;
    EVENT_INFO event_info;

    // pindex of every particle in a filter_particles result
    auto kept_pindex = [](const std::vector<part>& kept) {
        std::set<int> pindices;
        for (const auto& p : kept) pindices.insert(p.pindex);
        return pindices;
    };

    int  whileidx = 0;
    long n_particles = 0, n_prod = 0, n_loose = 0, n_nobeta = 0, n_matched = 0;

    while (_chain.Next() == true && (whileidx < maxEvents || maxEvents < 0)) {
        if (whileidx % 10000 == 0 && whileidx != 0) {
            std::cout << whileidx << " events read" << std::endl;
        }

        clas12ana.get_event_info(_c12, event_info);
        whileidx++;
        _cm.set_run(event_info.run);

        vec_particles = clas12ana.load_reco_particles(_c12);

        // Tag the scattered electron as hipo2tree_parnassus.C does -- filter_particles
        // reads the -1 marker on particle 0 to skip the electron-dependent cuts
        int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
        if (idx_scattered_ele != -1) vec_particles[idx_scattered_ele].is_scattered_electron = 1;
        else if (vec_particles.size() > 0) vec_particles[0].is_scattered_electron = -1;

        if (!is_real_data) {
            vec_mcparticles = clas12ana.load_mc_particles(_c12);
            clas12ana.match_mc_to_reco(vec_particles, vec_mcparticles);
        }

        std::set<int> kept_prod = kept_pindex(_cm.filter_particles(vec_particles));

        // Production cuts with chi2pid, the v_z / v_T outlier cuts and the e- / gamma calorimeter
        // edge cuts switched off: give a copy a chi2, vertex and PCAL edge distance those cuts
        // always accept. No other cut that filter_particles applies reads chi2, vx, vy, vz,
        // pcal_lv or pcal_lw.
        std::vector<part> loose = vec_particles;
        for (auto& p : loose) { p.chi2 = 0; p.vx = 0; p.vy = 0; p.vz = 0; p.pcal_lv = 100; p.pcal_lw = 100; }
        std::set<int> kept_loose = kept_pindex(_cm.filter_particles(loose));

        // Production cuts with the e+- / gamma beta outlier cut switched off: betaCut_outlier is
        // the only cut filter_particles applies that reads beta, and it always accepts beta = 1
        std::vector<part> nobeta = vec_particles;
        for (auto& p : nobeta) p.beta = 1;
        std::set<int> kept_nobeta = kept_pindex(_cm.filter_particles(nobeta));

        run   = event_info.run;
        evnum = event_info.evnum;
        for (const auto& particle : vec_particles) {
            q = particle;
            pass_prod   = kept_prod.count(q.pindex) > 0;
            pass_loose  = kept_loose.count(q.pindex) > 0;
            pass_nobeta = kept_nobeta.count(q.pindex) > 0;
            tree->Fill();

            n_particles++;
            n_prod    += pass_prod;
            n_loose   += pass_loose;
            n_nobeta  += pass_nobeta;
            n_matched += q.truepid != -999;
        }
    }

    fOut->cd();
    tree->Write();
    fOut->Close();

    auto pct = [n_particles](long n){ return n_particles > 0 ? 100.0 * n / n_particles : 0.0; };
    std::cout << "  events read         : " << whileidx << "\n"
              << "  reco particles      : " << n_particles << "\n"
              << "  pass_loose          : " << n_loose << "  (" << pct(n_loose) << "%)\n"
              << "  pass_nobeta         : " << n_nobeta << "  (" << pct(n_nobeta) << "%)\n"
              << "  pass_prod           : " << n_prod << "  (" << pct(n_prod) << "%)\n"
              << "  matched to MC truth : " << n_matched << "  (" << pct(n_matched) << "%)"
              << (is_real_data ? "  -- real data, no truth" : "") << "\n"
              << "  wrote " << outFile << std::endl;

    return 0;
}
