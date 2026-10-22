#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/CutStudy.h"
#include <filesystem>
#include <string>

// Measure what each particle-level cut does on its own: for every cut in the
// registry below, apply ONLY that cut and count how many events are left with
// no surviving particles (i.e. would be written with npflow==0).
//
// Must be run on a raw .hipo -- the prepped ROOT files are already filtered
// (hipo2tree_parnassus.C filters before it fills any tree), so the pre-cut
// particles only exist in the hipo.
//
// root -l -b -q 'macros/cut_study.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","plots/cut_study",20000)'
//
// Add a 4th argument to also save the MC truth particles of the emptied events
// (see captureCuts below):
// root -l -b -q 'macros/cut_study.C("prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","plots/cut_study",20000,"ALL")'

// ===========================================================================
//  CUT REGISTRY -- THIS IS THE PLACE TO EDIT
//
//  One line per cut. To study a new cut: add a line. To stop studying one:
//  delete or comment out its line. Nothing else needs to change.
//
//  Arguments are (group, name, pids, region, cut):
//    group  : label used to section the output table
//    name   : how the cut appears in the table and the CSVs
//    pids   : PIDs the cut is meaningful for; ALL_PIDS = every particle
//    region : IN_FD_ONLY  -> only judge particles in the forward detector
//             ALL_REGIONS -> judge them everywhere
//    cut    : returns true when the particle is KEPT
//
//  Use add_simple for cuts that only look at the particle, and
//  add_with_electron for cuts that also need the scattered electron (those are
//  skipped -- particle kept -- when no scattered electron was found, exactly as
//  filter_particles does).
//
//  The gating must mirror CutManager::filter_particles / apply_*_cuts, or the
//  numbers are meaningless: e.g. chi2pid() returns false for anything that is
//  not a charged pion, so it must be restricted to {211,-211}.
// ===========================================================================
namespace {

using CM = ProbeableCutManager;

const std::vector<int> ALL_PIDS  = {};
const bool             IN_FD_ONLY  = true;
const bool             ALL_REGIONS = false;

std::vector<CutProbe> build_cut_registry() {
  std::vector<CutProbe> probes;

  auto add_simple = [&probes](const char* group, const char* name,
                              std::vector<int> pids, bool fd_only,
                              std::function<bool(CM&, const part&)> cut) {
    probes.push_back(CutProbe{name, group, pids, fd_only,
        [cut](CM& cm, const part& p, const part&, bool){ return cut(cm, p); }});
  };

  auto add_with_electron = [&probes](const char* group, const char* name,
                                    std::vector<int> pids, bool fd_only,
                                    std::function<bool(CM&, const part&, const part&)> cut) {
    probes.push_back(CutProbe{name, group, pids, fd_only,
        [cut](CM& cm, const part& p, const part& e, bool found){
          return found ? cut(cm, p, e) : true;  // no electron -> cut not applied
        }});
  };

  // For whole cut groups, which take electron_found themselves and decide
  // internally which of their members to skip.
  auto add_group = [&probes](const char* group, const char* name, std::vector<int> pids,
                             CutProbe::Fn cut) {
    probes.push_back(CutProbe{name, group, pids, ALL_REGIONS, cut});
  };

  // --- applied to every particle, inside filter_particles -------------------
  add_simple("global", "theta_range_2_125", ALL_PIDS, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.thetaRangeCut(p); });
  add_simple("global", "FT_e_gamma_only",   ALL_PIDS, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.FTLeptonPhotonCut(p); });
  add_simple("global", "p_max_10.6",        ALL_PIDS, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.maxPCut(p); });
  add_simple("global", "unknown_pid_reject", {0},     ALL_REGIONS,
             [](CM&, const part&){ return false; });  // `case 0:` in filter_particles
  add_simple("global", "is_in_FD",             ALL_PIDS, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.is_in_FD(p); });

  // --- electron cuts (apply_electron_cuts) ---------------------------------
  add_simple("electron", "DC_fiducial",  {11}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.DC_fiducial_cut(p); });
  add_simple("electron", "EleSampFrac",  {11}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.EleSampFrac_pass2(p); });
  add_simple("electron", "minEpcal",     {11}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.minEpcal(p); });
  add_simple("electron", "caloEdges",    {11}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.caloEdges(p, 0); });
  add_simple("electron", "leptonBetaCut", {11}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.leptonBetaCut(p); });
  add_simple("electron", "VzCut",        {11}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.VzCut(p); });

  // --- positron cuts (case -11 in filter_particles) ------------------------
  add_simple("positron", "leptonBetaCut", {-11}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.leptonBetaCut(p); });

  // --- charged pion cuts (apply_pion_cuts) ---------------------------------
  add_simple("pion", "DC_fiducial",     {211,-211}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.DC_fiducial_cut(p); });
  add_simple("pion", "chi2pid_3sigma",  {211,-211}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.chi2pid(p, 2); });
  add_with_electron("pion", "vz_e_pi",  {211,-211}, ALL_REGIONS,
             [](CM& cm, const part& p, const part& e){ return cm.vz_e_pi(p, e); });

  // --- photon cuts (apply_photon_cuts) ------------------------------------
  add_simple("photon", "photonBetaCut", {22}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.photonBetaCut(p); });
  add_simple("photon", "caloEdges",     {22}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.caloEdges(p, 0); });
  add_simple("photon", "minEpcal",      {22}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.minEpcal(p); });
  add_with_electron("photon", "photonElectronAngle", {22}, ALL_REGIONS,
             [](CM& cm, const part& p, const part& e){ return cm.photonElectronAngle(p, e); });

  // --- neutron cuts (apply_neutron_cuts) ----------------------------------
  add_simple("neutron", "minNeutronEnergy", {2112}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.minNeutronEnergy(p); });
  add_simple("neutron", "caloEdges",        {2112}, IN_FD_ONLY,
             [](CM& cm, const part& p){ return cm.caloEdges(p, 0); });

  // --- whole per-PID groups, for cross-checking the single cuts above ------
  add_simple("composite", "apply_electron_cuts", {11}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.apply_electron_cuts(p); });
  add_group("composite", "apply_pion_cuts", {211,-211},
            [](CM& cm, const part& p, const part& e, bool found){ return cm.apply_pion_cuts(p, e, found); });
  add_group("composite", "apply_photon_cuts", {22},
            [](CM& cm, const part& p, const part& e, bool found){ return cm.apply_photon_cuts(p, e, found); });
  add_simple("composite", "apply_neutron_cuts", {2112}, ALL_REGIONS,
             [](CM& cm, const part& p){ return cm.apply_neutron_cuts(p); });

  // --- defined in CutManager but NOT applied by production ----------------
  //     uncomment to see what turning one on would cost
  // add_simple("whatif", "minEleMomentum", {11},       ALL_REGIONS,
  //            [](CM& cm, const part& p){ return cm.minEleMomentum(p); });
  // add_simple("whatif", "minPiMomentum",  {211,-211}, ALL_REGIONS,
  //            [](CM& cm, const part& p){ return cm.minPiMomentum(p); });
  // add_simple("whatif", "Ele3calo",       {11},       ALL_REGIONS,
  //            [](CM& cm, const part& p){ return cm.Ele3calo(p); });
  // add_simple("whatif", "photonMinEtot",  {22},       ALL_REGIONS,
  //            [](CM& cm, const part& p){ return cm.photonMinEtot(p); });
  // add_simple("whatif", "caloEdges_tight", {11},      IN_FD_ONLY,
  //            [](CM& cm, const part& p){ return cm.caloEdges(p, 1); });
  // add_simple("whatif", "chi2pid_strict", {211,-211}, ALL_REGIONS,
  //            [](CM& cm, const part& p){ return cm.chi2pid(p, 1); });

  return probes;
}

}  // namespace

// ===========================================================================
//  Event loop -- mirrors macros/hipo2tree_parnassus.C with the tree writing
//  stripped out. Only reco particles matter here.
// ===========================================================================
//  captureCuts: leave empty to only print the table. Otherwise it selects the
//  combination of cuts whose emptied events get written, as MC truth particles
//  in EventTree format, to <outDir>/<stem>/empty_events.root. Selectors are
//  comma-separated group or "group/name" entries from the registry above, or
//  "ALL" for the full production stack. Examples:
//      "ALL"                    every event production would write as npflow==0
//      "global,composite"       the registry's equivalent of the full stack
//      "photon/photonBetaCut"   events emptied by that one cut alone
int cut_study(const char* hipoFile = "",
              const char* outDir   = "plots/cut_study",
              const int   maxEvents = 100000,
              const char* captureCuts = ""){

    std::string stem   = std::filesystem::path(hipoFile).stem().string();
    std::string csvDir = (std::filesystem::path(outDir) / stem).string();

    // Configure CLAS12 Reader and HipoChain
    // -------------------------------------
    clas12root::HipoChain _chain;
    _chain.Add(hipoFile);
    auto _config_c12 = _chain.GetC12Reader();
    auto& _c12 = _chain.C12ref();

    // Add Analysis Objects
    // -------------------------------------
    ProbeableCutManager _cm = ProbeableCutManager();
    CLAS12Ana clas12ana = CLAS12Ana(_c12);
    clas12ana.set_run_config(_c12); // must precede get_event_info

    CutStudy study(build_cut_registry());
    study.set_capture_cuts(captureCuts);

    // Only opened when capture is requested, so the default run writes no ROOT
    // file. The TFile must exist before the EventTree, which attaches its TTree
    // to the current directory.
    TFile*     fEmpty     = nullptr;
    EventTree* treeTruth  = nullptr;   // MC truth particles
    EventTree* treeReco    = nullptr;  // reco particles, BEFORE the cuts
    if (study.capture_enabled()) {
        gSystem->mkdir(csvDir.c_str(), true);
        std::string emptyFile = csvDir + "/empty_events.root";
        fEmpty    = new TFile(emptyFile.c_str(), "RECREATE");
        treeTruth = new EventTree("truth_EventTree");
        treeReco  = new EventTree("reco_EventTree");
        std::cout << "  capturing events emptied by '" << captureCuts
                  << "' to " << emptyFile << std::endl;
    }

    std::vector<part> vec_particles;
    std::vector<part> vec_mcparticles;
    EVENT_INFO event_info;
    EVENT reco_event;
    EVENT mc_event;

    int  whileidx = 0;
    long ntrunc   = 0;   // captured events whose particle list overran EventTree

    // EventTree's fixed arrays hold kNmax=100 particles and FillTree does not
    // clamp, so truncate here rather than overrun them.
    auto clamp_to_tree = [&ntrunc](std::vector<part>& v) {
        if (v.size() > 100) { v.resize(100); ntrunc++; }
    };

    while (_chain.Next() == true && (whileidx < maxEvents || maxEvents < 0)) {
        if (whileidx % 10000 == 0 && whileidx != 0) {
            std::cout << whileidx << " events read" << std::endl;
        }

        clas12ana.get_event_info(_c12, event_info);
        event_info.uID = whileidx;
        whileidx++;
        _cm.set_run(event_info.run);

        vec_particles = clas12ana.load_reco_particles(_c12);

        // Clear the kinematics each event: they are only filled when a
        // scattered electron is found, and the captured events frequently have
        // none -- without this they would inherit the previous event's values.
        reco_event = EVENT();
        mc_event   = EVENT();

        // Tag the scattered electron exactly as the production macro does --
        // the is_scattered_electron==-1 marker on particle 0 is the signal
        // filter_particles reads to skip the electron-dependent cuts.
        int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
        if (idx_scattered_ele != -1){
            vec_particles[idx_scattered_ele].is_scattered_electron = 1;
            clas12ana.fill_reco_event_variables(reco_event, vec_particles);
        }
        if (idx_scattered_ele == -1 && vec_particles.size() > 0) {
            vec_particles[0].is_scattered_electron = -1;
        }

        // Save the emptied events. Both trees get one entry per captured event,
        // so they stay row-aligned: entry i of reco_EventTree is the pre-cut
        // reco content of the same event as entry i of truth_EventTree.
        // MC truth is only loaded for the events we actually save, so the
        // default (no capture) run never pays for it.
        if (study.process_event(_cm, vec_particles) && treeTruth) {
            vec_mcparticles = clas12ana.load_mc_particles(_c12);
            clas12ana.fill_mc_event_variables(mc_event, vec_mcparticles);
            clamp_to_tree(vec_mcparticles);
            treeTruth->FillTree(vec_mcparticles, mc_event, event_info);

            // vec_particles is still the unfiltered list -- process_event takes
            // it by const reference and filter_particles by value.
            clamp_to_tree(vec_particles);
            treeReco->FillTree(vec_particles, reco_event, event_info);
        }
    }

    study.report(std::cout);
    study.write_csv(csvDir);

    if (treeTruth) {
        fEmpty->cd();
        treeTruth->Write();
        treeReco->Write();
        fEmpty->Close();
        std::cout << "  wrote " << csvDir << "/empty_events.root  ("
                  << study.n_events_captured()
                  << " events; truth_EventTree + reco_EventTree)" << std::endl;
        if (ntrunc > 0) {
            std::cout << "  WARNING: " << ntrunc
                      << " particle lists had >100 particles and were truncated to 100"
                      << std::endl;
        }
    }

    return 0;
}
