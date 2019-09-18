
#include "WHTree.h"
#include "TTree.h"
#include "TDirectory.h"
#include <iostream>
#include <vector>

using namespace std;

WHTree::WHTree(TTree *tree){
    if(tree != NULL)
        Init(tree);
}

void WHTree::Init(TTree *tree){
    if(tree == NULL){
        // if we don't pass a tree, open in "write" mode. pointers have to be initialized,
        // and then create the branches.
        this->tree = new TTree("Events","");

        this->mcand_Km_pt              = new vector<float>;
        this->mcand_Kp_pt              = new vector<float>;
        this->mcand_mass_kaon          = new vector<float>;
        this->mcand_mass_pion          = new vector<float>;
        this->mcand_pt                 = new vector<float>;
        this->mcand_eta                = new vector<float>;
        this->mcand_phi                = new vector<float>;
        this->mcand_dR                 = new vector<float>;
        this->mcand_dR_genMeson        = new vector<float>;
        this->mcand_relIso             = new vector<float>;
        this->mcand_photon_mass        = new vector<float>;
        this->mcand_photon_dR          = new vector<float>;
        this->mcand_fromPV             = new vector<bool>;
        this->vhAngles_cosThetaStar    = new vector<float>;
        this->vhAngles_cosTheta1       = new vector<float>;
        this->vhAngles_cosTheta2       = new vector<float>;
        this->vhAngles_phi             = new vector<float>;
        this->vhAngles_phi1            = new vector<float>;
        this->vhAngles_m1              = new vector<float>;
        this->vhAngles_m2              = new vector<float>;

        b_evt_run                  = this->tree->Branch("evt_run", &evt_run, "evt_run/I");
        b_evt_lumi                 = this->tree->Branch("evt_lumi", &evt_lumi, "evt_lumi/I");
        b_evt_event                = this->tree->Branch("evt_event", &evt_event, "evt_event/L");
        b_isGolden                 = this->tree->Branch("isGolden", &isGolden, "isGolden/O");
        b_scale1fb                 = this->tree->Branch("scale1fb", &scale1fb, "scale1fb/F");
        b_genH_pt                  = this->tree->Branch("genH_pt", &genH_pt, "genH_pt/F");
        b_genH_eta                 = this->tree->Branch("genH_eta", &genH_eta, "genH_eta/F");
        b_genH_phi                 = this->tree->Branch("genH_phi", &genH_phi, "genH_phi/F");
        b_genHMeson_pdgId          = this->tree->Branch("genHMeson_pdgId", &genHMeson_pdgId, "genHMeson_pdgId/I");
        b_genHMeson_pt             = this->tree->Branch("genHMeson_pt", &genHMeson_pt, "genHMeson_pt/F");
        b_genHMeson_eta            = this->tree->Branch("genHMeson_eta", &genHMeson_eta, "genHMeson_eta/F");
        b_genHMeson_phi            = this->tree->Branch("genHMeson_phi", &genHMeson_phi, "genHMeson_phi/F");
        b_genHPhoton_pt            = this->tree->Branch("genHPhoton_pt", &genHPhoton_pt, "genHPhoton_pt/F");
        b_genHPhoton_eta           = this->tree->Branch("genHPhoton_eta", &genHPhoton_eta, "genHPhoton_eta/F");
        b_genHPhoton_phi           = this->tree->Branch("genHPhoton_phi", &genHPhoton_phi, "genHPhoton_phi/F");
        b_genKm_pt                 = this->tree->Branch("genKm_pt", &genKm_pt, "genKm_pt/F");
        b_genKm_eta                = this->tree->Branch("genKm_eta", &genKm_eta, "genKm_eta/F");
        b_genKm_phi                = this->tree->Branch("genKm_phi", &genKm_phi, "genKm_phi/F");
        b_genKp_pt                 = this->tree->Branch("genKp_pt", &genKp_pt, "genKp_pt/F");
        b_genKp_eta                = this->tree->Branch("genKp_eta", &genKp_eta, "genKp_eta/F");
        b_genKp_phi                = this->tree->Branch("genKp_phi", &genKp_phi, "genKp_phi/F");
        b_genKK_dR                 = this->tree->Branch("genKK_dR", &genKK_dR, "genKK_dR/F");
        b_genLep_pdgId             = this->tree->Branch("genLep_pdgId", &genLep_pdgId, "genLep_pdgId/F");
        b_genLep_pt                = this->tree->Branch("genLep_pt", &genLep_pt, "genLep_pt/F");
        b_genLep_eta               = this->tree->Branch("genLep_eta", &genLep_eta, "genLep_eta/F");
        b_genLep_phi               = this->tree->Branch("genLep_phi", &genLep_phi, "genLep_phi/F");
        b_genNeutrino_pdgId        = this->tree->Branch("genNeutrino_pdgId", &genNeutrino_pdgId, "genNeutrino_pdgId/F");
        b_genNeutrino_pt           = this->tree->Branch("genNeutrino_pt", &genNeutrino_pt, "genNeutrino_pt/F");
        b_genNeutrino_eta          = this->tree->Branch("genNeutrino_eta", &genNeutrino_eta, "genNeutrino_eta/F");
        b_genNeutrino_phi          = this->tree->Branch("genNeutrino_phi", &genNeutrino_phi, "genNeutrino_phi/F");
        b_genvhAngles_cosThetaStar = this->tree->Branch("genvhAngles_cosThetaStar", &genvhAngles_cosThetaStar, "genvhAngles_cosThetaStar/F");
        b_genvhAngles_cosTheta1    = this->tree->Branch("genvhAngles_cosTheta1", &genvhAngles_cosTheta1, "genvhAngles_cosTheta1/F");
        b_genvhAngles_cosTheta2    = this->tree->Branch("genvhAngles_cosTheta2", &genvhAngles_cosTheta2, "genvhAngles_cosTheta2/F");
        b_genvhAngles_phi          = this->tree->Branch("genvhAngles_phi", &genvhAngles_phi, "genvhAngles_phi/F");
        b_genvhAngles_phi1         = this->tree->Branch("genvhAngles_phi1", &genvhAngles_phi1, "genvhAngles_phi1/F");
        b_genvhAngles_m1           = this->tree->Branch("genvhAngles_m1", &genvhAngles_m1, "genvhAngles_m1/F");
        b_genvhAngles_m2           = this->tree->Branch("genvhAngles_m2", &genvhAngles_m2, "genvhAngles_m2/F");
        b_nphoton                  = this->tree->Branch("nphoton", &nphoton, "nphoton/I");
        b_photon_pt                = this->tree->Branch("photon_pt", &photon_pt, "photon_pt/F");
        b_photon_eta               = this->tree->Branch("photon_eta", &photon_eta, "photon_eta/F");
        b_photon_phi               = this->tree->Branch("photon_phi", &photon_phi, "photon_phi/F");
        b_photon_isTightCutBased   = this->tree->Branch("photon_isTightCutBased", &photon_isTightCutBased, "photon_isTightCutBased/O");
        b_photon_relChHadIso       = this->tree->Branch("photon_relChHadIso", &photon_relChHadIso, "photon_relChHadIso/F");
        b_photon_dR_higgsPhoton    = this->tree->Branch("photon_dR_higgsPhoton", &photon_dR_higgsPhoton, "photon_dR_higgsPhoton/F");
        b_photon_isMatchedGen      = this->tree->Branch("photon_isMatchedGen", &photon_isMatchedGen, "photon_isMatchedGen/O");
        b_photon_drMinParton       = this->tree->Branch("photon_drMinParton", &photon_drMinParton, "photon_drMinParton/F");
        b_nmcand                   = this->tree->Branch("nmcand", &nmcand, "nmcand/I");
        b_mcand_Km_pt              = this->tree->Branch("mcand_Km_pt", mcand_Km_pt);
        b_mcand_Kp_pt              = this->tree->Branch("mcand_Kp_pt", mcand_Kp_pt);
        b_mcand_mass_kaon          = this->tree->Branch("mcand_mass_kaon", mcand_mass_kaon);
        b_mcand_mass_pion          = this->tree->Branch("mcand_mass_pion", mcand_mass_pion);
        b_mcand_pt                 = this->tree->Branch("mcand_pt", mcand_pt);
        b_mcand_eta                = this->tree->Branch("mcand_eta", mcand_eta);
        b_mcand_phi                = this->tree->Branch("mcand_phi", mcand_phi);
        b_mcand_dR                 = this->tree->Branch("mcand_dR", mcand_dR);
        b_mcand_dR_genMeson        = this->tree->Branch("mcand_dR_genMeson", mcand_dR_genMeson);
        b_mcand_relIso             = this->tree->Branch("mcand_relIso", mcand_relIso);
        b_mcand_photon_mass        = this->tree->Branch("mcand_photon_mass", mcand_photon_mass);
        b_mcand_photon_dR          = this->tree->Branch("mcand_photon_dR", mcand_photon_dR);
        b_mcand_fromPV             = this->tree->Branch("mcand_fromPV", mcand_fromPV);
        b_best_phiCand_idx         = this->tree->Branch("best_phiCand_idx", &best_phiCand_idx, "best_phiCand_idx/I");
        b_best_rhoCand_idx         = this->tree->Branch("best_rhoCand_idx", &best_rhoCand_idx, "best_rhoCand_idx/I");
        b_nlep                     = this->tree->Branch("nlep", &nlep, "nlep/I");
        b_lepton_pt                = this->tree->Branch("lepton_pt", &lepton_pt, "lepton_pt/F");
        b_lepton_eta               = this->tree->Branch("lepton_eta", &lepton_eta, "lepton_eta/F");
        b_lepton_phi               = this->tree->Branch("lepton_phi", &lepton_phi, "lepton_phi/F");
        b_lepton_pdgId             = this->tree->Branch("lepton_pdgId", &lepton_pdgId, "lepton_pdgId/I");
        b_lepton_miniRelIso        = this->tree->Branch("lepton_miniRelIso", &lepton_miniRelIso, "lepton_miniRelIso/F");
        b_lepton_photon_dR         = this->tree->Branch("lepton_photon_dR", &lepton_photon_dR, "lepton_photon_dR/F");
        b_met_rawPt                = this->tree->Branch("met_rawPt", &met_rawPt, "met_rawPt/F");
        b_met_rawPhi               = this->tree->Branch("met_rawPhi", &met_rawPhi, "met_rawPhi/F");
        b_met_pt                   = this->tree->Branch("met_pt", &met_pt, "met_pt/F");
        b_met_phi                  = this->tree->Branch("met_phi", &met_phi, "met_phi/F");
        b_vhAngles_cosThetaStar    = this->tree->Branch("vhAngles_cosThetaStar", vhAngles_cosThetaStar);
        b_vhAngles_cosTheta1       = this->tree->Branch("vhAngles_cosTheta1", vhAngles_cosTheta1);
        b_vhAngles_cosTheta2       = this->tree->Branch("vhAngles_cosTheta2", vhAngles_cosTheta2);
        b_vhAngles_phi             = this->tree->Branch("vhAngles_phi", vhAngles_phi);
        b_vhAngles_phi1            = this->tree->Branch("vhAngles_phi1", vhAngles_phi1);
        b_vhAngles_m1              = this->tree->Branch("vhAngles_m1", vhAngles_m1);
        b_vhAngles_m2              = this->tree->Branch("vhAngles_m2", vhAngles_m2);
        b_HLT_SingleEl             = this->tree->Branch("HLT_SingleEl", &HLT_SingleEl, "HLT_SingleEl/I");
        b_HLT_SingleMu             = this->tree->Branch("HLT_SingleMu", &HLT_SingleMu, "HLT_SingleMu/I");
        b_passFilters              = this->tree->Branch("passFilters", &passFilters, "passFilters/O");
        b_isHEM                    = this->tree->Branch("isHEM", &isHEM, "isHEM/O");

        Reset();
    }else{
        // if we do pass a tree, open in "read" mode
        this->tree = tree;
        this->tree->SetMakeClass(1);

        this->tree->SetBranchAddress("evt_run", &evt_run, &b_evt_run);
        this->tree->SetBranchAddress("evt_lumi", &evt_lumi, &b_evt_lumi);
        this->tree->SetBranchAddress("evt_event", &evt_event, &b_evt_event);
        this->tree->SetBranchAddress("isGolden", &isGolden, &b_isGolden);
        this->tree->SetBranchAddress("scale1fb", &scale1fb, &b_scale1fb);
        this->tree->SetBranchAddress("genH_pt", &genH_pt, &b_genH_pt);
        this->tree->SetBranchAddress("genH_eta", &genH_eta, &b_genH_eta);
        this->tree->SetBranchAddress("genH_phi", &genH_phi, &b_genH_phi);
        this->tree->SetBranchAddress("genHMeson_pdgId", &genHMeson_pdgId, &b_genHMeson_pdgId);
        this->tree->SetBranchAddress("genHMeson_pt", &genHMeson_pt, &b_genHMeson_pt);
        this->tree->SetBranchAddress("genHMeson_eta", &genHMeson_eta, &b_genHMeson_eta);
        this->tree->SetBranchAddress("genHMeson_phi", &genHMeson_phi, &b_genHMeson_phi);
        this->tree->SetBranchAddress("genHPhoton_pt", &genHPhoton_pt, &b_genHPhoton_pt);
        this->tree->SetBranchAddress("genHPhoton_eta", &genHPhoton_eta, &b_genHPhoton_eta);
        this->tree->SetBranchAddress("genHPhoton_phi", &genHPhoton_phi, &b_genHPhoton_phi);
        this->tree->SetBranchAddress("genKm_pt", &genKm_pt, &b_genKm_pt);
        this->tree->SetBranchAddress("genKm_eta", &genKm_eta, &b_genKm_eta);
        this->tree->SetBranchAddress("genKm_phi", &genKm_phi, &b_genKm_phi);
        this->tree->SetBranchAddress("genKp_pt", &genKp_pt, &b_genKp_pt);
        this->tree->SetBranchAddress("genKp_eta", &genKp_eta, &b_genKp_eta);
        this->tree->SetBranchAddress("genKp_phi", &genKp_phi, &b_genKp_phi);
        this->tree->SetBranchAddress("genKK_dR", &genKK_dR, &b_genKK_dR);
        this->tree->SetBranchAddress("genLep_pdgId", &genLep_pdgId, &b_genLep_pdgId);
        this->tree->SetBranchAddress("genLep_pt", &genLep_pt, &b_genLep_pt);
        this->tree->SetBranchAddress("genLep_eta", &genLep_eta, &b_genLep_eta);
        this->tree->SetBranchAddress("genLep_phi", &genLep_phi, &b_genLep_phi);
        this->tree->SetBranchAddress("genNeutrino_pdgId", &genNeutrino_pdgId, &b_genNeutrino_pdgId);
        this->tree->SetBranchAddress("genNeutrino_pt", &genNeutrino_pt, &b_genNeutrino_pt);
        this->tree->SetBranchAddress("genNeutrino_eta", &genNeutrino_eta, &b_genNeutrino_eta);
        this->tree->SetBranchAddress("genNeutrino_phi", &genNeutrino_phi, &b_genNeutrino_phi);
        this->tree->SetBranchAddress("genvhAngles_cosThetaStar", &genvhAngles_cosThetaStar, &b_genvhAngles_cosThetaStar);
        this->tree->SetBranchAddress("genvhAngles_cosTheta1", &genvhAngles_cosTheta1, &b_genvhAngles_cosTheta1);
        this->tree->SetBranchAddress("genvhAngles_cosTheta2", &genvhAngles_cosTheta2, &b_genvhAngles_cosTheta2);
        this->tree->SetBranchAddress("genvhAngles_phi", &genvhAngles_phi, &b_genvhAngles_phi);
        this->tree->SetBranchAddress("genvhAngles_phi1", &genvhAngles_phi1, &b_genvhAngles_phi1);
        this->tree->SetBranchAddress("genvhAngles_m1", &genvhAngles_m1, &b_genvhAngles_m1);
        this->tree->SetBranchAddress("genvhAngles_m2", &genvhAngles_m2, &b_genvhAngles_m2);
        this->tree->SetBranchAddress("nphoton", &nphoton, &b_nphoton);
        this->tree->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
        this->tree->SetBranchAddress("photon_eta", &photon_eta, &b_photon_eta);
        this->tree->SetBranchAddress("photon_phi", &photon_phi, &b_photon_phi);
        this->tree->SetBranchAddress("photon_isTightCutBased", &photon_isTightCutBased, &b_photon_isTightCutBased);
        this->tree->SetBranchAddress("photon_relChHadIso", &photon_relChHadIso, &b_photon_relChHadIso);
        this->tree->SetBranchAddress("photon_dR_higgsPhoton", &photon_dR_higgsPhoton, &b_photon_dR_higgsPhoton);
        this->tree->SetBranchAddress("photon_isMatchedGen", &photon_isMatchedGen, &b_photon_isMatchedGen);
        this->tree->SetBranchAddress("photon_drMinParton", &photon_drMinParton, &b_photon_drMinParton);
        this->tree->SetBranchAddress("nmcand", &nmcand, &b_nmcand);
        this->tree->SetBranchAddress("mcand_Km_pt", &mcand_Km_pt, &b_mcand_Km_pt);
        this->tree->SetBranchAddress("mcand_Kp_pt", &mcand_Kp_pt, &b_mcand_Kp_pt);
        this->tree->SetBranchAddress("mcand_mass_kaon", &mcand_mass_kaon, &b_mcand_mass_kaon);
        this->tree->SetBranchAddress("mcand_mass_pion", &mcand_mass_pion, &b_mcand_mass_pion);
        this->tree->SetBranchAddress("mcand_pt", &mcand_pt, &b_mcand_pt);
        this->tree->SetBranchAddress("mcand_eta", &mcand_eta, &b_mcand_eta);
        this->tree->SetBranchAddress("mcand_phi", &mcand_phi, &b_mcand_phi);
        this->tree->SetBranchAddress("mcand_dR", &mcand_dR, &b_mcand_dR);
        this->tree->SetBranchAddress("mcand_dR_genMeson", &mcand_dR_genMeson, &b_mcand_dR_genMeson);
        this->tree->SetBranchAddress("mcand_relIso", &mcand_relIso, &b_mcand_relIso);
        this->tree->SetBranchAddress("mcand_photon_mass", &mcand_photon_mass, &b_mcand_photon_mass);
        this->tree->SetBranchAddress("mcand_photon_dR", &mcand_photon_dR, &b_mcand_photon_dR);
        this->tree->SetBranchAddress("mcand_fromPV", &mcand_fromPV, &b_mcand_fromPV);
        this->tree->SetBranchAddress("best_phiCand_idx", &best_phiCand_idx, &b_best_phiCand_idx);
        this->tree->SetBranchAddress("best_rhoCand_idx", &best_rhoCand_idx, &b_best_rhoCand_idx);
        this->tree->SetBranchAddress("nlep", &nlep, &b_nlep);
        this->tree->SetBranchAddress("lepton_pt", &lepton_pt, &b_lepton_pt);
        this->tree->SetBranchAddress("lepton_eta", &lepton_eta, &b_lepton_eta);
        this->tree->SetBranchAddress("lepton_phi", &lepton_phi, &b_lepton_phi);
        this->tree->SetBranchAddress("lepton_pdgId", &lepton_pdgId, &b_lepton_pdgId);
        this->tree->SetBranchAddress("lepton_miniRelIso", &lepton_miniRelIso, &b_lepton_miniRelIso);
        this->tree->SetBranchAddress("lepton_photon_dR", &lepton_photon_dR, &b_lepton_photon_dR);
        this->tree->SetBranchAddress("met_rawPt", &met_rawPt, &b_met_rawPt);
        this->tree->SetBranchAddress("met_rawPhi", &met_rawPhi, &b_met_rawPhi);
        this->tree->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
        this->tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
        this->tree->SetBranchAddress("vhAngles_cosThetaStar", &vhAngles_cosThetaStar, &b_vhAngles_cosThetaStar);
        this->tree->SetBranchAddress("vhAngles_cosTheta1", &vhAngles_cosTheta1, &b_vhAngles_cosTheta1);
        this->tree->SetBranchAddress("vhAngles_cosTheta2", &vhAngles_cosTheta2, &b_vhAngles_cosTheta2);
        this->tree->SetBranchAddress("vhAngles_phi", &vhAngles_phi, &b_vhAngles_phi);
        this->tree->SetBranchAddress("vhAngles_phi1", &vhAngles_phi1, &b_vhAngles_phi1);
        this->tree->SetBranchAddress("vhAngles_m1", &vhAngles_m1, &b_vhAngles_m1);
        this->tree->SetBranchAddress("vhAngles_m2", &vhAngles_m2, &b_vhAngles_m2);
        this->tree->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl, &b_HLT_SingleEl);
        this->tree->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
        this->tree->SetBranchAddress("passFilters", &passFilters, &b_passFilters);
        this->tree->SetBranchAddress("isHEM", &isHEM, &b_isHEM);

    }
}

void WHTree::Fill(){
    tree->Fill();
}

void WHTree::Reset(){
    evt_run = -999;
    evt_lumi = -999;
    evt_event = -999;
    isGolden = -999;
    scale1fb = -999;
    genH_pt = -999;
    genH_eta = -999;
    genH_phi = -999;
    genHMeson_pdgId = -999;
    genHMeson_pt = -999;
    genHMeson_eta = -999;
    genHMeson_phi = -999;
    genHPhoton_pt = -999;
    genHPhoton_eta = -999;
    genHPhoton_phi = -999;
    genKm_pt = -999;
    genKm_eta = -999;
    genKm_phi = -999;
    genKp_pt = -999;
    genKp_eta = -999;
    genKp_phi = -999;
    genKK_dR = -999;
    genLep_pdgId = -999;
    genLep_pt = -999;
    genLep_eta = -999;
    genLep_phi = -999;
    genNeutrino_pdgId = -999;
    genNeutrino_pt = -999;
    genNeutrino_eta = -999;
    genNeutrino_phi = -999;
    genvhAngles_cosThetaStar = -999;
    genvhAngles_cosTheta1 = -999;
    genvhAngles_cosTheta2 = -999;
    genvhAngles_phi = -999;
    genvhAngles_phi1 = -999;
    genvhAngles_m1 = -999;
    genvhAngles_m2 = -999;
    nphoton = -999;
    photon_pt = -999;
    photon_eta = -999;
    photon_phi = -999;
    photon_isTightCutBased = -999;
    photon_relChHadIso = -999;
    photon_dR_higgsPhoton = -999;
    photon_isMatchedGen = -999;
    photon_drMinParton = -999;
    nmcand = -999;
    mcand_Km_pt->clear();
    mcand_Kp_pt->clear();
    mcand_mass_kaon->clear();
    mcand_mass_pion->clear();
    mcand_pt->clear();
    mcand_eta->clear();
    mcand_phi->clear();
    mcand_dR->clear();
    mcand_dR_genMeson->clear();
    mcand_relIso->clear();
    mcand_photon_mass->clear();
    mcand_photon_dR->clear();
    mcand_fromPV->clear();
    best_phiCand_idx = -999;
    best_rhoCand_idx = -999;
    nlep = -999;
    lepton_pt = -999;
    lepton_eta = -999;
    lepton_phi = -999;
    lepton_pdgId = -999;
    lepton_miniRelIso = -999;
    lepton_photon_dR = -999;
    met_rawPt = -999;
    met_rawPhi = -999;
    met_pt = -999;
    met_phi = -999;
    vhAngles_cosThetaStar->clear();
    vhAngles_cosTheta1->clear();
    vhAngles_cosTheta2->clear();
    vhAngles_phi->clear();
    vhAngles_phi1->clear();
    vhAngles_m1->clear();
    vhAngles_m2->clear();
    HLT_SingleEl = -999;
    HLT_SingleMu = -999;
    passFilters = -999;
    isHEM = -999;

}

void WHTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}

void WHTree::GetEntry(ULong64_t i){
    this->tree->GetEntry(i);
}
