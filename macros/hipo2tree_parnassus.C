#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include <filesystem>
#include <string>

// root -l -b -q 'macros/hipo2tree_parnassus.C("/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","prepped_data",-1)'

int hipo2tree_parnassus(const char* hipoFile = "",
    const char* outDir = "",
    const int maxEvents = 100){

    std::string out_name = std::filesystem::path(hipoFile).stem().string();
    std::string outputFile = (std::filesystem::path(outDir) / (out_name + ".root")).string();
     
    // Create a TFile to save the data
      TFile* fOut = new TFile(outputFile.c_str(), "RECREATE");
    // Create a TTree to store the data
    EventTree* treeReco = new EventTree("reco_EventTree");
    EventTree* treeMC = new EventTree("mc_EventTree");
    ParnassusTree* treeParnassus = new ParnassusTree("evt_tree");
    
        // Configure CLAS12 Reader and HipoChain
    // -------------------------------------
    clas12root::HipoChain _chain;
    _chain.Add(hipoFile);

    //use the config to select only the events in the hipo files that have the particles we are interested in
    auto _config_c12=_chain.GetC12Reader();

    //generate a reference to the object that contains all the information for the desired events
    auto& _c12 = _chain.C12ref();

    // Add Analysis Objects
    // -------------------------------------
    CutManager _cm = CutManager();
    CLAS12Ana clas12ana = CLAS12Ana(_c12);
    clas12ana.set_run_config(_c12); //grabs the event number index, run number index, and torus info index, saving it to the corresponding object private methods for clas12ana

    // Add Analysis Structs - from the Structs.h src file
    // -------------------------------------
    std::vector<part> vec_particles;
    std::vector<part> vec_mcparticles;
    EVENT_INFO event_info;
    EVENT reco_event;
    EVENT mc_event;

    int whileidx = 0;
    int _ievent = 0;

    while (_chain.Next() == true && (whileidx < maxEvents || maxEvents < 0)) {
        if (whileidx % 10000 == 0 && whileidx != 0) {
            std::cout << whileidx << " events read | "
                      << _ievent * 100.0 / whileidx << "% passed event selection "
                      << std::endl;
        }

        //for each event, populate the neccesary information in the clas12ana and cutmanager contexts using info about each event found in the constants class
        //----------------------------------------------------------------------------------------------------------------------------------------------------------
        clas12ana.get_event_info(_c12, event_info); //this line fills event_info with the relevant properties of that event based found in constants.h
        event_info.uID = whileidx; //keeps track of the iteration in the event_info method
        whileidx++;
        _cm.set_run(event_info.run); //sets relevant runinfo for the cutmanager class to use 

        //clear reco_event and mc_event structs to prepare for new event
        reco_event.reset();
        mc_event.reset();

        //colect reconstructed particle info
        vec_particles = clas12ana.load_reco_particles(_c12);
        //determine which particle is the scattered electron
        int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
        //if the event contains a scattered electron, fill the reco_event struct with the relevant information about that event. Otherwise leave it blank
        if (idx_scattered_ele != -1){
            vec_particles[idx_scattered_ele].is_scattered_electron = 1;
            clas12ana.fill_reco_event_variables(reco_event, vec_particles);
        }
        //if the event doesn't contain a scattered electron, then skip the event 
        if (idx_scattered_ele == -1 && vec_particles.size() > 0) {
            vec_particles[0].is_scattered_electron = -1;
        }
        
        vec_particles = _cm.filter_particles(vec_particles);

        //the scattered electron itself may have been removed by the cuts above. If it was, the DIS
        //variables computed from it no longer describe anything in the output particle list, so put
        //them back to -999 to stay consistent with the scattered-electron flag written to the trees
        if (idx_scattered_ele != -1){
            bool scattered_ele_survived = false;
            for (const auto& p : vec_particles){
                if (p.pid == 11 && p.is_scattered_electron == 1) { scattered_ele_survived = true; break; }
            }
            if (!scattered_ele_survived) reco_event.reset();
        }

        //collect MC particle information
        vec_mcparticles = clas12ana.load_mc_particles(_c12);
        clas12ana.fill_mc_event_variables(mc_event,vec_mcparticles);
        
        //then fill the tree with the vector of particles, event, and event_info
        treeReco->FillTree(vec_particles, reco_event, event_info);
        treeMC->FillTree(vec_mcparticles, mc_event, event_info);
        treeParnassus->FillTree(vec_particles, vec_mcparticles, reco_event, mc_event, event_info);
        _ievent++;
    }

    fOut->cd();
    treeReco->Write();
    treeMC->Write();
    treeParnassus->Write();
    fOut->Close();


    return 0;
    }
