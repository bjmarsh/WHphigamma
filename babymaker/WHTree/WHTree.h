
#ifndef WHTree_h
#define WHTree_h

#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"


using namespace std;

class WHTree {
  public:

    int            evt_run;
    int            evt_lumi;
    long           evt_event;
    bool           isGolden;
    float          scale1fb;
    float          genH_pt;
    float          genH_eta;
    float          genH_phi;
    int            genHMeson_pdgId;
    float          genHMeson_pt;
    float          genHMeson_eta;
    float          genHMeson_phi;
    float          genHPhoton_pt;
    float          genHPhoton_eta;
    float          genHPhoton_phi;
    float          genKm_pt;
    float          genKm_eta;
    float          genKm_phi;
    float          genKp_pt;
    float          genKp_eta;
    float          genKp_phi;
    float          genKK_dR;
    float          genLep_pdgId;
    float          genLep_pt;
    float          genLep_eta;
    float          genLep_phi;
    float          genNeutrino_pdgId;
    float          genNeutrino_pt;
    float          genNeutrino_eta;
    float          genNeutrino_phi;
    float          genvhAngles_cosThetaStar;
    float          genvhAngles_cosTheta1;
    float          genvhAngles_cosTheta2;
    float          genvhAngles_phi;
    float          genvhAngles_phi1;
    float          genvhAngles_m1;
    float          genvhAngles_m2;
    int            nphoton;
    float          photon_pt;
    float          photon_eta;
    float          photon_phi;
    bool           photon_isTightCutBased;
    float          photon_relChHadIso;
    float          photon_dR_higgsPhoton;
    bool           photon_isMatchedGen;
    float          photon_drMinParton;
    int            nmcand;
    vector<float>* mcand_Km_pt = 0;
    vector<float>* mcand_Kp_pt = 0;
    vector<float>* mcand_mass_kaon = 0;
    vector<float>* mcand_mass_pion = 0;
    vector<float>* mcand_pt = 0;
    vector<float>* mcand_eta = 0;
    vector<float>* mcand_phi = 0;
    vector<float>* mcand_dR = 0;
    vector<float>* mcand_dR_genMeson = 0;
    vector<float>* mcand_relIso = 0;
    vector<float>* mcand_photon_mass = 0;
    vector<float>* mcand_photon_dR = 0;
    vector<bool>*  mcand_fromPV = 0;
    int            best_phiCand_idx;
    int            best_rhoCand_idx;
    int            nlep;
    float          lepton_pt;
    float          lepton_eta;
    float          lepton_phi;
    int            lepton_pdgId;
    float          lepton_miniRelIso;
    float          lepton_photon_dR;
    float          met_rawPt;
    float          met_rawPhi;
    float          met_pt;
    float          met_phi;
    vector<float>* vhAngles_cosThetaStar = 0;
    vector<float>* vhAngles_cosTheta1 = 0;
    vector<float>* vhAngles_cosTheta2 = 0;
    vector<float>* vhAngles_phi = 0;
    vector<float>* vhAngles_phi1 = 0;
    vector<float>* vhAngles_m1 = 0;
    vector<float>* vhAngles_m2 = 0;
    int            HLT_SingleEl;
    int            HLT_SingleMu;
    bool           passFilters;
    bool           isHEM;

    WHTree(TTree *tree=0);
    void Init(TTree *tree=0);
    void Fill();
    void Reset();
    void Write(TDirectory *d);
    void GetEntry(ULong64_t entry);
    static void progress(int nEventsTotal, int nEventsChain);

  private:
    TTree *tree;

    TBranch *b_evt_run = 0;
    TBranch *b_evt_lumi = 0;
    TBranch *b_evt_event = 0;
    TBranch *b_isGolden = 0;
    TBranch *b_scale1fb = 0;
    TBranch *b_genH_pt = 0;
    TBranch *b_genH_eta = 0;
    TBranch *b_genH_phi = 0;
    TBranch *b_genHMeson_pdgId = 0;
    TBranch *b_genHMeson_pt = 0;
    TBranch *b_genHMeson_eta = 0;
    TBranch *b_genHMeson_phi = 0;
    TBranch *b_genHPhoton_pt = 0;
    TBranch *b_genHPhoton_eta = 0;
    TBranch *b_genHPhoton_phi = 0;
    TBranch *b_genKm_pt = 0;
    TBranch *b_genKm_eta = 0;
    TBranch *b_genKm_phi = 0;
    TBranch *b_genKp_pt = 0;
    TBranch *b_genKp_eta = 0;
    TBranch *b_genKp_phi = 0;
    TBranch *b_genKK_dR = 0;
    TBranch *b_genLep_pdgId = 0;
    TBranch *b_genLep_pt = 0;
    TBranch *b_genLep_eta = 0;
    TBranch *b_genLep_phi = 0;
    TBranch *b_genNeutrino_pdgId = 0;
    TBranch *b_genNeutrino_pt = 0;
    TBranch *b_genNeutrino_eta = 0;
    TBranch *b_genNeutrino_phi = 0;
    TBranch *b_genvhAngles_cosThetaStar = 0;
    TBranch *b_genvhAngles_cosTheta1 = 0;
    TBranch *b_genvhAngles_cosTheta2 = 0;
    TBranch *b_genvhAngles_phi = 0;
    TBranch *b_genvhAngles_phi1 = 0;
    TBranch *b_genvhAngles_m1 = 0;
    TBranch *b_genvhAngles_m2 = 0;
    TBranch *b_nphoton = 0;
    TBranch *b_photon_pt = 0;
    TBranch *b_photon_eta = 0;
    TBranch *b_photon_phi = 0;
    TBranch *b_photon_isTightCutBased = 0;
    TBranch *b_photon_relChHadIso = 0;
    TBranch *b_photon_dR_higgsPhoton = 0;
    TBranch *b_photon_isMatchedGen = 0;
    TBranch *b_photon_drMinParton = 0;
    TBranch *b_nmcand = 0;
    TBranch *b_mcand_Km_pt = 0;
    TBranch *b_mcand_Kp_pt = 0;
    TBranch *b_mcand_mass_kaon = 0;
    TBranch *b_mcand_mass_pion = 0;
    TBranch *b_mcand_pt = 0;
    TBranch *b_mcand_eta = 0;
    TBranch *b_mcand_phi = 0;
    TBranch *b_mcand_dR = 0;
    TBranch *b_mcand_dR_genMeson = 0;
    TBranch *b_mcand_relIso = 0;
    TBranch *b_mcand_photon_mass = 0;
    TBranch *b_mcand_photon_dR = 0;
    TBranch *b_mcand_fromPV = 0;
    TBranch *b_best_phiCand_idx = 0;
    TBranch *b_best_rhoCand_idx = 0;
    TBranch *b_nlep = 0;
    TBranch *b_lepton_pt = 0;
    TBranch *b_lepton_eta = 0;
    TBranch *b_lepton_phi = 0;
    TBranch *b_lepton_pdgId = 0;
    TBranch *b_lepton_miniRelIso = 0;
    TBranch *b_lepton_photon_dR = 0;
    TBranch *b_met_rawPt = 0;
    TBranch *b_met_rawPhi = 0;
    TBranch *b_met_pt = 0;
    TBranch *b_met_phi = 0;
    TBranch *b_vhAngles_cosThetaStar = 0;
    TBranch *b_vhAngles_cosTheta1 = 0;
    TBranch *b_vhAngles_cosTheta2 = 0;
    TBranch *b_vhAngles_phi = 0;
    TBranch *b_vhAngles_phi1 = 0;
    TBranch *b_vhAngles_m1 = 0;
    TBranch *b_vhAngles_m2 = 0;
    TBranch *b_HLT_SingleEl = 0;
    TBranch *b_HLT_SingleMu = 0;
    TBranch *b_passFilters = 0;
    TBranch *b_isHEM = 0;

};

#endif
