#include "Structs.h"

class BaseTree {
    protected: // common variables
        int A, evnum, fgId, run, hel,uID,hwp,target,tSign;
        double x,y,Q2,nu,W,Pol,tPol;
        double truex,truey,trueQ2,truenu,trueW;
        TTree *tree;
    
    public:
        BaseTree(const char * treename) {
            tree=new TTree(treename,treename);
        }
        BaseTree(){}
    
        virtual ~BaseTree() {}
        virtual void CreateBranches() = 0;
        TTree* GetTree() const {
            return tree;
        }
    
        int GetEntries() { return tree->GetEntries(); }
        void Write() {    tree->Write();    }
        void Write(const char * c) {    tree->Write(c);    }
};



class EventTree : public BaseTree {
    private:
        const static int kNmax = 100;
        int Nmax;
        int pindex[kNmax], status[kNmax], pid[kNmax], truepid[kNmax], trueparentid[kNmax], trueparentpid[kNmax], trueparentparentid[kNmax], trueparentparentpid[kNmax];
        double px[kNmax], py[kNmax], pz[kNmax], p[kNmax], E[kNmax],m[kNmax];
        double vx[kNmax], vy[kNmax], vz[kNmax], chi2[kNmax], beta[kNmax];
        double truepx[kNmax], truepy[kNmax], truepz[kNmax], truep[kNmax], trueE[kNmax], truebeta[kNmax];
        int is_CFR[kNmax];
        double theta[kNmax], eta[kNmax], phi[kNmax], truept[kNmax], truem[kNmax], truetheta[kNmax], trueeta[kNmax], truephi[kNmax], truevx[kNmax], truevy[kNmax], truevz[kNmax];
        int pcal_sector[kNmax], ecin_sector[kNmax], ecout_sector[kNmax];
        double pcal_x[kNmax], pcal_y[kNmax], pcal_z[kNmax];
        double ecin_x[kNmax], ecin_y[kNmax], ecin_z[kNmax];
        double ecout_x[kNmax], ecout_y[kNmax], ecout_z[kNmax];
        double pcal_e[kNmax], pcal_lu[kNmax], pcal_lv[kNmax], pcal_lw[kNmax], pcal_m2u[kNmax], pcal_m2v[kNmax], pcal_m2w[kNmax];
        double ecin_e[kNmax], ecin_lu[kNmax], ecin_lv[kNmax], ecin_lw[kNmax], ecin_m2u[kNmax], ecin_m2v[kNmax], ecin_m2w[kNmax];
        double ecout_e[kNmax], ecout_lu[kNmax], ecout_lv[kNmax], ecout_lw[kNmax], ecout_m2u[kNmax], ecout_m2v[kNmax], ecout_m2w[kNmax];
        double nphe_ltcc[kNmax];
        double nphe_htcc[kNmax];
        int sector[kNmax];
        double traj_x1[kNmax], traj_y1[kNmax], traj_z1[kNmax], traj_x2[kNmax], traj_y2[kNmax], traj_z2[kNmax], traj_x3[kNmax], traj_y3[kNmax], traj_z3[kNmax];
        double edge_1[kNmax], edge_2[kNmax],edge_3[kNmax];
        
    
    public:
        EventTree(const char * treename) : BaseTree(treename) {            
            CreateBranches();
        }
    
        ~EventTree() {}

        void CreateBranches() override {
            tree->Branch("A",&A,"A/I");  
            tree->Branch("evnum",&evnum,"evnum/I");  
            tree->Branch("uID",&uID,"uID/I");
            tree->Branch("run",&run,"run/I");
            tree->Branch("Pol",&Pol,"Pol/D");
            tree->Branch("tPol",&tPol,"tPol/D");
            tree->Branch("hwp",&hwp,"hwp/I");
            tree->Branch("target",&target,"target/I");
            tree->Branch("tSign",&tSign,"tSign/I");
            tree->Branch("Nmax",&Nmax,"Nmax/I");
            tree->Branch("x", &x, "x/D");
            tree->Branch("y", &y, "y/D");
            tree->Branch("W", &W, "W/D");
            tree->Branch("Q2",&Q2,"Q2/D");
            tree->Branch("nu", &nu, "nu/D");
            tree->Branch("truex", &truex, "truex/D");
            tree->Branch("truey", &truey, "truey/D");
            tree->Branch("trueQ2",&trueQ2,"trueQ2/D");
            tree->Branch("trueW", &trueW, "trueW/D");
            tree->Branch("truenu", &truenu, "truenu/D");
            tree->Branch("hel", &hel, "hel/I");

            tree->Branch("truex", &truex, "truex/D");
            tree->Branch("truey", &truey, "truey/D");
            tree->Branch("trueW", &trueW, "trueW/D");
            tree->Branch("truenu", &truenu, "truenu/D");

            tree->Branch("pindex", pindex, "pindex[Nmax]/I");
            tree->Branch("status", status, "status[Nmax]/I");
            tree->Branch("px", px, "px[Nmax]/D");
            tree->Branch("py", py, "py[Nmax]/D");
            tree->Branch("pz", pz, "pz[Nmax]/D");
            tree->Branch("p", p, "p[Nmax]/D");
            tree->Branch("E", E, "E[Nmax]/D");
            tree->Branch("pid", pid, "pid[Nmax]/I");
            tree->Branch("vx", vx, "vx[Nmax]/D");
            tree->Branch("vy", vy, "vy[Nmax]/D");
            tree->Branch("vz", vz, "vz[Nmax]/D");
            tree->Branch("chi2", chi2, "chi2[Nmax]/D");
            tree->Branch("beta", beta, "beta[Nmax]/D");
            tree->Branch("m", m, "m[Nmax]/D");
            tree->Branch("theta", theta, "theta[Nmax]/D");
            tree->Branch("eta", eta, "eta[Nmax]/D");
            tree->Branch("phi", phi, "phi[Nmax]/D");
            tree->Branch("truepx", truepx, "truepx[Nmax]/D");
            tree->Branch("truepy", truepy, "truepy[Nmax]/D");
            tree->Branch("truepz", truepz, "truepz[Nmax]/D");
            tree->Branch("truep", truep, "truep[Nmax]/D");
            tree->Branch("truept", truept, "truept[Nmax]/D");
            tree->Branch("truem", truem, "truem[Nmax]/D");
            tree->Branch("truetheta", truetheta, "truetheta[Nmax]/D");
            tree->Branch("trueeta", trueeta, "trueeta[Nmax]/D");
            tree->Branch("truephi", truephi, "truephi[Nmax]/D");
            tree->Branch("truevx", truevx, "truevx[Nmax]/D");
            tree->Branch("truevy", truevy, "truevy[Nmax]/D");
            tree->Branch("truevz", truevz, "truevz[Nmax]/D");
            tree->Branch("trueE", trueE, "trueE[Nmax]/D");
            tree->Branch("truebeta", truebeta, "truebeta[Nmax]/D");
            tree->Branch("is_CFR", is_CFR, "is_CFR[Nmax]/I");
            tree->Branch("truepid", truepid, "truepid[Nmax]/I");
            tree->Branch("trueparentid", trueparentid, "trueparentid[Nmax]/I");
            tree->Branch("trueparentpid", trueparentpid, "trueparentpid[Nmax]/I");
            tree->Branch("trueparentparentid", trueparentparentid, "trueparentparentid[Nmax]/I");
            tree->Branch("trueparentparentpid", trueparentparentpid, "trueparentparentpid[Nmax]/I");

            tree->Branch("pcal_sector", pcal_sector, "pcal_sector[Nmax]/I");
            tree->Branch("pcal_e", pcal_e, "pcal_e[Nmax]/D");
            tree->Branch("pcal_x", pcal_x, "pcal_x[Nmax]/D");
            tree->Branch("pcal_y", pcal_y, "pcal_y[Nmax]/D");
            tree->Branch("pcal_z", pcal_z, "pcal_z[Nmax]/D");
            tree->Branch("pcal_lu", pcal_lu, "pcal_lu[Nmax]/D");
            tree->Branch("pcal_lv", pcal_lv, "pcal_lv[Nmax]/D");
            tree->Branch("pcal_lw", pcal_lw, "pcal_lw[Nmax]/D");
            tree->Branch("pcal_m2u", pcal_m2u, "pcal_m2u[Nmax]/D");
            tree->Branch("pcal_m2v", pcal_m2v, "pcal_m2v[Nmax]/D");
            tree->Branch("pcal_m2w", pcal_m2w, "pcal_m2w[Nmax]/D");

            tree->Branch("ecin_sector", ecin_sector, "ecin_sector[Nmax]/I");
            tree->Branch("ecin_e", ecin_e, "ecin_e[Nmax]/D");
            tree->Branch("ecin_x", ecin_x, "ecin_x[Nmax]/D");
            tree->Branch("ecin_y", ecin_y, "ecin_y[Nmax]/D");
            tree->Branch("ecin_z", ecin_z, "ecin_z[Nmax]/D");
            tree->Branch("ecin_lu", ecin_lu, "ecin_lu[Nmax]/D");
            tree->Branch("ecin_lv", ecin_lv, "ecin_lv[Nmax]/D");
            tree->Branch("ecin_lw", ecin_lw, "ecin_lw[Nmax]/D");
            tree->Branch("ecin_m2u", ecin_m2u, "ecin_m2u[Nmax]/D");
            tree->Branch("ecin_m2v", ecin_m2v, "ecin_m2v[Nmax]/D");
            tree->Branch("ecin_m2w", ecin_m2w, "ecin_m2w[Nmax]/D");

            tree->Branch("ecout_sector", ecout_sector, "ecout_sector[Nmax]/I");
            tree->Branch("ecout_e", ecout_e, "ecout_e[Nmax]/D");
            tree->Branch("ecout_x", ecout_x, "ecout_x[Nmax]/D");
            tree->Branch("ecout_y", ecout_y, "ecout_y[Nmax]/D");
            tree->Branch("ecout_z", ecout_z, "ecout_z[Nmax]/D");
            tree->Branch("ecout_lu", ecout_lu, "ecout_lu[Nmax]/D");
            tree->Branch("ecout_lv", ecout_lv, "ecout_lv[Nmax]/D");
            tree->Branch("ecout_lw", ecout_lw, "ecout_lw[Nmax]/D");
            tree->Branch("ecout_m2u", ecout_m2u, "ecout_m2u[Nmax]/D");
            tree->Branch("ecout_m2v", ecout_m2v, "ecout_m2v[Nmax]/D");
            tree->Branch("ecout_m2w", ecout_m2w, "ecout_m2w[Nmax]/D");

            tree->Branch("sector", sector, "sector[Nmax]/I");
            tree->Branch("traj_x1", traj_x1, "traj_x1[Nmax]/D");
            tree->Branch("traj_y1", traj_y1, "traj_y1[Nmax]/D");
            tree->Branch("traj_z1", traj_z1, "traj_z1[Nmax]/D");
            tree->Branch("traj_x2", traj_x2, "traj_x2[Nmax]/D");
            tree->Branch("traj_y2", traj_y2, "traj_y2[Nmax]/D");
            tree->Branch("traj_z2", traj_z2, "traj_z2[Nmax]/D");
            tree->Branch("traj_x3", traj_x3, "traj_x3[Nmax]/D");
            tree->Branch("traj_y3", traj_y3, "traj_y3[Nmax]/D");
            tree->Branch("traj_z3", traj_z3, "traj_z3[Nmax]/D");
            tree->Branch("edge_1", edge_1, "edge_1[Nmax]/D");
            tree->Branch("edge_2", edge_2, "edge_2[Nmax]/D");
            tree->Branch("edge_3", edge_3, "edge_3[Nmax]/D");

            tree->Branch("nphe_ltcc", nphe_ltcc, "nphe_ltcc[Nmax]/D");
            tree->Branch("nphe_htcc", nphe_htcc, "nphe_htcc[Nmax]/D");
        }
    
    void FillTree(std::vector<part> vec_particles, EVENT &event, EVENT_INFO &event_info){
        A = event_info.A;
        evnum = event_info.evnum;
        uID = event_info.uID;
        fgId = event_info.fgId;
        run = event_info.run;
        tPol = event_info.tPol;
        hwp = event_info.hwp;
        target = event_info.target;
        tSign = event_info.tSign;
        hel = event_info.hel;
        x = event.x;
        y = event.y;
        Q2 = event.Q2;
        nu = event.nu;
        W = event.W;
        Pol = event_info.Pol;
        truex = event.truex;
        truey = event.truey;
        trueQ2 = event.trueQ2;
        truenu = event.truenu;
        trueW = event.trueW;

        Nmax = vec_particles.size();
        for( int i = 0 ; i < Nmax; i++){
            part par = vec_particles[i];

            pindex[i] = par.pindex;
            status[i] = par.status;
            px[i] = par.px;
            py[i] = par.py;
            pz[i] = par.pz;
            p[i] = par.p;
            E[i] = par.E;
            pid[i] = par.pid;
            vx[i] = par.vx;
            vy[i] = par.vy;
            vz[i] = par.vz;
            chi2[i] = par.chi2;
            beta[i] = par.beta;
            m[i] = par.m;
            theta[i] = par.theta;
            eta[i] = par.eta;
            phi[i] = par.phi;
            truepx[i] = par.truepx;
            truepy[i] = par.truepy;
            truepz[i] = par.truepz;
            truep[i] = par.truep;
            truept[i] = par.truept;
            trueE[i] = par.trueE;
            truem[i] = par.truem;
            truebeta[i] = par.truebeta;
            truetheta[i] = par.truetheta;
            trueeta[i] = par.trueeta;
            truephi[i] = par.truephi;
            truevx[i] = par.truevx;
            truevy[i] = par.truevy;
            truevz[i] = par.truevz;
            is_CFR[i] = par.is_CFR;
            truepid[i] = par.truepid;
            trueparentid[i] = par.trueparentid;
            trueparentpid[i] = par.trueparentpid;
            trueparentparentid[i] = par.trueparentparentid;
            trueparentparentpid[i] = par.trueparentparentpid;
            pcal_sector[i] = par.pcal_sector;
            pcal_e[i] = par.pcal_e;
            pcal_x[i] = par.pcal_x;
            pcal_y[i] = par.pcal_y;
            pcal_z[i] = par.pcal_z;
            pcal_lu[i] = par.pcal_lu;
            pcal_lv[i] = par.pcal_lv;
            pcal_lw[i] = par.pcal_lw;
            pcal_m2u[i] = par.pcal_m2u;
            pcal_m2v[i] = par.pcal_m2v;
            pcal_m2w[i] = par.pcal_m2w;
            ecin_sector[i] = par.ecin_sector;
            ecin_e[i] = par.ecin_e;
            ecin_x[i] = par.ecin_x;
            ecin_y[i] = par.ecin_y;
            ecin_z[i] = par.ecin_z;
            ecin_lu[i] = par.ecin_lu;
            ecin_lv[i] = par.ecin_lv;
            ecin_lw[i] = par.ecin_lw;
            ecin_m2u[i] = par.ecin_m2u;
            ecin_m2v[i] = par.ecin_m2v;
            ecin_m2w[i] = par.ecin_m2w;
            ecout_sector[i] = par.ecout_sector;
            ecout_e[i] = par.ecout_e;
            ecout_x[i] = par.ecout_x;
            ecout_y[i] = par.ecout_y;
            ecout_z[i] = par.ecout_z;
            ecout_lu[i] = par.ecout_lu;
            ecout_lv[i] = par.ecout_lv;
            ecout_lw[i] = par.ecout_lw;
            ecout_m2u[i] = par.ecout_m2u;
            ecout_m2v[i] = par.ecout_m2v;
            ecout_m2w[i] = par.ecout_m2w;
            sector[i] = par.sector;
            traj_x1[i] = par.traj_x1;
            traj_y1[i] = par.traj_y1;
            traj_z1[i] = par.traj_z1;
            traj_x2[i] = par.traj_x2;
            traj_y2[i] = par.traj_y2;
            traj_z2[i] = par.traj_z2;
            traj_x3[i] = par.traj_x3;
            traj_y3[i] = par.traj_y3;
            traj_z3[i] = par.traj_z3;
            edge_1[i] = par.edge_1;
            edge_2[i] = par.edge_2;
            edge_3[i] = par.edge_3;
            nphe_ltcc[i] = par.nphe_ltcc;
            nphe_htcc[i] = par.nphe_htcc;
        }
        tree->Fill();
    }
};

class DihadronTree : public BaseTree {
    private:
        double gamma,eps,depolA,depolB,depolC,depolV,depolW;
        double Mh, z, M1, M2, M12, phi_h, phi_R0, phi_R1, th, z1, z2, xF1, xF2, xF, Mx, phi_h1, phi_h2, delta_phi_h;
        double pT1, pT2, pTtot, P1, P2, Ptot;
        double truegamma,trueeps,truedepolA,truedepolB,truedepolC,truedepolV,truedepolW;
        double truex, trueQ2, trueW, truey, trueM1, trueM2, trueM12, trueMh, truephi_h, truephi_R0, truephi_R1, trueth;
        double truez1, truez2, truexF1, truexF2, truez, truexF, trueMx, truephi_h1, truephi_h2, truedelta_phi_h;
        double truepT1, truepT2, truepTtot, trueP1, trueP2, truePtot;
        int truepid_1, truepid_2, truepid_11, truepid_12, truepid_21, truepid_22;
        int trueparentpid_1, trueparentpid_2, trueparentid_1, trueparentid_2;
        int trueparentparentpid_1, trueparentparentpid_2, trueparentparentid_1, trueparentparentid_2;
        int is_CFR_1, is_CFR_2;
        double p_11, p_12, p_21, p_22;
        int MCmatch, isGoodEventWithoutML;
    public:
        DihadronTree(const char * treename) : BaseTree(treename) {            
            CreateBranches();
        }
    
        ~DihadronTree() {}

        void CreateBranches() override {
            // create branches specific to the DihadronTree
            tree->Branch("A", &A, "A/I");
            tree->Branch("evnum", &evnum, "evnum/I");
            tree->Branch("uID", &uID, "uID/I");
            tree->Branch("fgId", &fgId, "fgId/D");
            tree->Branch("run", &run, "run/I");
            tree->Branch("Pol", &Pol, "Pol/D");
            tree->Branch("hel", &hel, "hel/I");
            tree->Branch("MCmatch", &MCmatch, "MCmatch/I");
            tree->Branch("isGoodEventWithoutML", &isGoodEventWithoutML, "isGoodEventWithoutML/I");
            tree->Branch("is_CFR_1",&is_CFR_1, "is_CFR_1/I");
            tree->Branch("is_CFR_2",&is_CFR_2, "is_CFR_2/I");
            tree->Branch("gamma",&gamma,"gamma/D");
            tree->Branch("eps", &eps, "eps/D");
            tree->Branch("depolA", &depolA, "depolA/D");
            tree->Branch("depolB", &depolB, "depolB/D");
            tree->Branch("depolC", &depolC, "depolC/D");
            tree->Branch("depolV", &depolV, "depolV/D");
            tree->Branch("depolW", &depolW, "depolW/D");
            tree->Branch("x", &x, "x/D");
            tree->Branch("Q2", &Q2, "Q2/D");
            tree->Branch("W", &W, "W/D");
            tree->Branch("y", &y, "y/D");
            tree->Branch("M1", &M1, "M1/D");
            tree->Branch("M2", &M2, "M2/D");
            tree->Branch("M12",&M12,"M12/D");
            tree->Branch("Mh", &Mh, "Mh/D");
            tree->Branch("phi_h", &phi_h, "phi_h/D");
            tree->Branch("phi_R0", &phi_R0, "phi_R0/D");
            tree->Branch("phi_R1", &phi_R1, "phi_R1/D");
            tree->Branch("th", &th, "th/D");
            tree->Branch("z1", &z1, "z1/D");
            tree->Branch("z2", &z2, "z2/D");
            tree->Branch("xF1", &xF1, "xF1/D");
            tree->Branch("xF2", &xF2, "xF2/D");
            tree->Branch("z", &z, "z/D");
            tree->Branch("xF", &xF, "xF/D");
            tree->Branch("Mx", &Mx, "Mx/D");
            tree->Branch("phi_h1", &phi_h1, "phi_h1/D");
            tree->Branch("phi_h2", &phi_h2, "phi_h2/D");
            tree->Branch("delta_phi_h", &delta_phi_h, "delta_phi_h/D");
            tree->Branch("pT1", &pT1, "pT1/D");
            tree->Branch("pT2", &pT2, "pT2/D");
            tree->Branch("pTtot", &pTtot, "pTtot/D");
            tree->Branch("P1", &P1, "P1/D");
            tree->Branch("P2", &P2, "P2/D");
            tree->Branch("Ptot", &Ptot, "Ptot/D");
            tree->Branch("truex", &truex, "truex/D");
            tree->Branch("trueQ2", &trueQ2, "trueQ2/D");
            tree->Branch("trueW", &trueW, "trueW/D");
            tree->Branch("truey", &truey, "truey/D");
            tree->Branch("trueM1", &trueM1, "trueM1/D");
            tree->Branch("trueM2", &trueM2, "trueM2/D");
            tree->Branch("trueM12",&trueM12,"trueM12/D");
            tree->Branch("trueMh", &trueMh, "trueMh/D");
            tree->Branch("truephi_h", &truephi_h, "truephi_h/D");
            tree->Branch("truephi_R0", &truephi_R0, "truephi_R0/D");
            tree->Branch("truephi_R1", &truephi_R1, "truephi_R1/D");
            tree->Branch("trueth", &trueth, "trueth/D");
            tree->Branch("truez1", &truez1, "truez1/D");
            tree->Branch("truez2", &truez2, "truez2/D");
            tree->Branch("truexF1", &truexF1, "truexF1/D");
            tree->Branch("truexF2", &truexF2, "truexF2/D");
            tree->Branch("truez", &truez, "truez/D");
            tree->Branch("truexF", &truexF, "truexF/D");
            tree->Branch("trueMx", &trueMx, "trueMx/D");
            tree->Branch("truephi_h1", &truephi_h1, "truephi_h1/D");
            tree->Branch("truephi_h2", &truephi_h2, "truephi_h2/D");
            tree->Branch("truedelta_phi_h", &truedelta_phi_h, "truedelta_phi_h/D");
            tree->Branch("truepT1", &truepT1, "truepT1/D");
            tree->Branch("truepT2", &truepT2, "truepT2/D");
            tree->Branch("truepTtot", &truepTtot, "truepTtot/D");
            tree->Branch("trueP1", &trueP1, "trueP1/D");
            tree->Branch("trueP2", &trueP2, "trueP2/D");
            tree->Branch("truePtot", &truePtot, "truePtot/D");
            tree->Branch("truepid_1", &truepid_1, "truepid_1/I");
            tree->Branch("truepid_2", &truepid_2, "truepid_2/I");
            tree->Branch("truepid_11", &truepid_11, "truepid_11/I");
            tree->Branch("truepid_12", &truepid_12, "truepid_12/I");
            tree->Branch("truepid_21", &truepid_21, "truepid_21/I");
            tree->Branch("truepid_22", &truepid_22, "truepid_22/I");
            tree->Branch("trueparentpid_1", &trueparentpid_1, "trueparentpid_1/I");
            tree->Branch("trueparentpid_2", &trueparentpid_2, "trueparentpid_2/I");
            tree->Branch("trueparentid_1", &trueparentid_1, "trueparentid_1/I");
            tree->Branch("trueparentid_2", &trueparentid_2, "trueparentid_2/I");
            tree->Branch("trueparentparentpid_1", &trueparentparentpid_1, "trueparentparentpid_1/I");
            tree->Branch("trueparentparentpid_2", &trueparentparentpid_2, "trueparentparentpid_2/I");
            tree->Branch("trueparentparentid_1", &trueparentparentid_1, "trueparentparentid_1/I");
            tree->Branch("trueparentparentid_2", &trueparentparentid_2, "trueparentparentid_2/I");
            tree->Branch("p_11", &p_11,"p_11/D");
            tree->Branch("p_12", &p_12,"p_12/D");
            tree->Branch("p_21", &p_21,"p_21/D");
            tree->Branch("p_22", &p_22,"p_22/D");
        }
    
        void FillTree(EVENT &event, EVENT &tevent, EVENT_INFO &event_info){
            A = event_info.A;
            evnum = event_info.evnum;
            uID = event_info.uID;
            fgId = event_info.fgId;
            run = event_info.run;
            hel = event_info.hel;
            x = event.x;
            y = event.y;
            Q2 = event.Q2;
            nu = event.nu;
            W = event.W;
            Pol = event_info.Pol;
            truex = tevent.truex;
            truey = tevent.truey;
            trueQ2 = tevent.trueQ2;
            truenu = tevent.truenu;
            trueW = tevent.trueW;
            
            gamma = event.gamma;
            eps = event.eps;
            depolA = event.depolA;
            depolB = event.depolB;
            depolC = event.depolC;
            depolV = event.depolV;
            depolW = event.depolW;
            Mh = event.Mh;
            z = event.z;
            M1 = event.M1;
            M2 = event.M2;
            M12 = event.M12;
            phi_h = event.phi_h;
            phi_R0 = event.phi_R0;
            phi_R1 = event.phi_R1;
            th = event.th;
            z1 = event.z1;
            z2 = event.z2;
            xF1 = event.xF1;
            xF2 = event.xF2;
            xF = event.xF;
            Mx = event.Mx;
            phi_h1 = event.phi_h1;
            phi_h2 = event.phi_h2;
            delta_phi_h = event.delta_phi_h;
            pT1 = event.pT1;
            pT2 = event.pT2;
            pTtot = event.pTtot;
            P1 = event.P1;
            P2 = event.P2;
            Ptot = event.Ptot;
            truex = tevent.truex;
            trueQ2 = tevent.trueQ2;
            trueW = tevent.trueW;
            truey = tevent.truey;
            trueM1 = tevent.trueM1;
            trueM2 = tevent.trueM2;
            trueM12 = tevent.trueM12;
            trueMh = tevent.trueMh;
            truephi_h = tevent.truephi_h;
            truephi_R0 = tevent.truephi_R0;
            truephi_R1 = tevent.truephi_R1;
            truedelta_phi_h = tevent.truedelta_phi_h;
            truez1 = tevent.truez1;
            truez2 = tevent.truez2;
            truexF1 = tevent.truexF1;
            truexF2 = tevent.truexF2;
            truez = tevent.truez;
            truexF = tevent.truexF;
            trueMx = tevent.trueMx;
            truephi_h1 = tevent.truephi_h1;
            truephi_h2 = tevent.truephi_h2;
            truepT1 = tevent.truepT1;
            truepT2 = tevent.truepT2;
            truepTtot = tevent.truepTtot;
            trueP1 = tevent.trueP1;
            trueP2 = tevent.trueP2;
            truePtot = tevent.truePtot;
            truepid_1 = tevent.truepid_1;
            truepid_2 = tevent.truepid_2;
            truepid_11 = tevent.truepid_11;
            truepid_12 = tevent.truepid_12;
            truepid_21 = tevent.truepid_21;
            truepid_22 = tevent.truepid_22;
            trueparentpid_1 = tevent.trueparentpid_1;
            trueparentpid_2 = tevent.trueparentpid_2;
            trueparentid_1 = tevent.trueparentid_1;
            trueparentid_2 = tevent.trueparentid_2;
            trueparentparentpid_1 = tevent.trueparentparentpid_1;
            trueparentparentpid_2 = tevent.trueparentparentpid_2;
            trueparentparentid_1 = tevent.trueparentparentid_1;
            trueparentparentid_2 = tevent.trueparentparentid_2;
            p_11 = event.p_11;
            p_12 = event.p_12;
            p_21 = event.p_21;
            p_22 = event.p_22;
            tree->Fill();
        }

};


class ParnassusTree : public BaseTree {
    private:
        const static int kNmax = 100;
        int ntruth, npflow, eventNumber;

        double pflow_pt[kNmax], pflow_phi[kNmax], pflow_eta[kNmax], pflow_vx[kNmax], pflow_vy[kNmax], pflow_vz[kNmax],pflow_beta[kNmax];
        int pflow_class[kNmax], pflow_pdgId[kNmax], pflow_charge[kNmax];
        double truth_pt[kNmax], truth_eta[kNmax], truth_phi[kNmax], truth_vx[kNmax], truth_vy[kNmax], truth_vz[kNmax],truth_beta[kNmax];
        int truth_class[kNmax], truth_pdgId[kNmax], truth_charge[kNmax];

        Kinematics kin;

        int GetClassID(double pid){

            if (abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212 || abs(pid) == 45 || abs(pid) == 1000010020) return 0; //charged hadron (45 = deuteron, CLAS12 EB code)
            if (abs(pid) == 11) return 1; //electron
            if (abs(pid) == 13) return 2; //muon
            if (abs(pid) == 111 || abs(pid) == 130 || abs(pid) == 2112 ||abs(pid) == 310 || abs(pid) ==311 ||abs(pid) == 3122 ||abs(pid) ==3212 || abs(pid) ==3322) return 3; //neutron hadron
            if (abs(pid) == 22) return 4; //photon
            return -1; //unknown particle
        }
        int GetCharge(double pid){
            if (pid == 211 || pid == 321 || pid == 2212 || pid == 11 || pid==13 || pid == 45 || pid == 1000010020) return 1; //positive charge
            if (pid == -211 || pid == -321 || pid == -2212 || pid == -11 || pid==-13) return -1; //negative charge
            return 0; //neutral particle
        }
        
    public:
        ParnassusTree(const char * treename) : BaseTree(treename) {            
            CreateBranches();
        }
    
        ~ParnassusTree() {}

        void CreateBranches() override {  
            tree->Branch("eventNumber",&eventNumber,"eventNumber/I");
            tree->Branch("ntruth",&ntruth,"ntruth/I");
            tree->Branch("npflow",&npflow,"npflow/I");
            
            tree->Branch("pflow_pt", pflow_pt, "pflow_pt[npflow]/D");
            tree->Branch("pflow_phi", pflow_phi, "pflow_phi[npflow]/D");
            tree->Branch("pflow_eta", pflow_eta, "pflow_eta[npflow]/D");
            tree->Branch("pflow_vx", pflow_vx, "pflow_vx[npflow]/D");
            tree->Branch("pflow_vy", pflow_vy, "pflow_vy[npflow]/D"); 
            tree->Branch("pflow_vz", pflow_vz, "pflow_vz[npflow]/D"); 
            tree->Branch("pflow_class", pflow_class, "pflow_class[npflow]/I");
            tree->Branch("pflow_pdgId", pflow_pdgId, "pflow_pdgId[npflow]/I");
            tree->Branch("pflow_beta", pflow_beta, "pflow_beta[npflow]/D");
            tree->Branch("pflow_charge", pflow_charge, "pflow_charge[npflow]/I");

            tree->Branch("truth_pt", truth_pt, "truth_pt[ntruth]/D");
            tree->Branch("truth_phi", truth_phi, "truth_phi[ntruth]/D");
            tree->Branch("truth_eta", truth_eta, "truth_eta[ntruth]/D");
            tree->Branch("truth_vx", truth_vx, "truth_vx[ntruth]/D");
            tree->Branch("truth_vy", truth_vy, "truth_vy[ntruth]/D"); 
            tree->Branch("truth_vz", truth_vz, "truth_vz[ntruth]/D"); 
            tree->Branch("truth_class", truth_class, "truth_class[ntruth]/I");
            tree->Branch("truth_pdgId", truth_pdgId, "truth_pdgId[ntruth]/I");
            tree->Branch("truth_beta", truth_beta, "truth_beta[ntruth]/D");
            tree->Branch("truth_charge", truth_charge, "truth_charge[ntruth]/I");
        }
    
        void FillTree(std::vector<part> vec_particles, std::vector<part> vec_mcparticles, 
            EVENT &pflow_event,EVENT &mc_event, EVENT_INFO &event_info){
            eventNumber = event_info.evnum;
            
            npflow = vec_particles.size();
            for( int i = 0 ; i < npflow; i++){
                part par = vec_particles[i];

                //convert CLAS specific PIDs to PDG IDs
                if(par.pid == 45) par.pid = 1000010020; // deuteron

                pflow_pt[i] = kin.Pt(par.px, par.py);
                pflow_phi[i] = par.phi;
                pflow_eta[i] = kin.eta(par.theta);
                pflow_vx[i] = par.vx;
                pflow_vy[i] = par.vy;
                pflow_vz[i] = par.vz;
                pflow_pdgId[i] = par.pid;
                pflow_class[i] = GetClassID(par.pid);
                pflow_beta[i] = par.beta;
                pflow_charge[i] = GetCharge(par.pid);
            }

            ntruth = vec_mcparticles.size();
            for( int i = 0 ; i < ntruth; i++){
                part par = vec_mcparticles[i];

                //convert CLAS specific PIDs to PDG IDs
                if(par.truepid == 45) par.truepid = 1000010020; // deuteron

                truth_pt[i] = kin.Pt(par.truepx, par.truepy);
                truth_phi[i] = par.truephi;
                truth_eta[i] = kin.eta(par.truetheta);
                truth_vx[i] = par.truevx;
                truth_vy[i] = par.truevy;
                truth_vz[i] = par.truevz;
                truth_pdgId[i] = par.truepid;
                truth_class[i] = GetClassID(par.truepid);       
                truth_beta[i] = par.truebeta;
                truth_charge[i] = GetCharge(par.truepid);
            }

            tree->Fill();
        }
};


