// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TTreeCache.h"

// CMS3
#include "CORE/CMS3.h"
#include "CORE/Tools/utils.h"
#include "CORE/Tools/goodrun.h"
#include "CORE/Base.h"
#include "CORE/ElectronSelections.h"
#include "CORE/MuonSelections.h"
#include "CORE/PhotonSelections.h"
#include "CORE/MetSelections.h"
#include "CORE/TriggerSelections.h"
#include "CORE/IsolationTools.h"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"
#include "WHTree/WHTree.h"
#include "vhAngles.h"


using namespace std;
using namespace tas;

TLorentzVector lvconv(LorentzVector p4){
    TLorentzVector k;
    k.SetPtEtaPhiE(p4.pt(), p4.eta(), p4.phi(), p4.E());
    return k;       
}

int ScanChain( TChain* chain, string outdir, string outid, int year=2018, bool isData=false, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    gconf.year = year;
    gconf.ea_version = year==2016 ? 1 : 3;
    float lumi = year==2016 ? 35.92 : year==2017 ? 41.53 : 59.74;
    cout << "isData? " << isData << endl;
    cout << "Year: " << gconf.year << endl;
    cout << "EA version: " << gconf.ea_version << endl;

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    TFile *fout = new TFile(Form("output_files/%s/%s.root",outdir.c_str(),outid.c_str()),"RECREATE");

    // output tree
    WHTree t;
    t.Init();

    // Loop over events to Analyze
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    if( nEvents >= 0 ) nEventsChain = nEvents;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    if(isData){
        if(year==2016)
            set_goodrun_file("jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_snt.txt");
        if(year==2017)
            set_goodrun_file("jsons/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1_snt.txt");
        if(year==2018)
            set_goodrun_file("jsons/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_snt.txt");
    }

    DatasetInfoFromFile datasetInfoFromFile;

    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) {

        // Get File Content
        TString currentFileName = currentFile->GetTitle();
        cout << currentFileName << endl;
        TFile file( currentFile->GetTitle() );
        if(file.IsZombie()){
            cout << "WARNING: couldn't open file! Skipping. " << currentFile->GetTitle() << endl;
            continue;
        }

        bool isSignal = currentFileName.Contains("WH");

        // setup jet energy corrections
        std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
        FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3(0);
        string jecTag = "";
        if(isData){
            if(currentFileName.Contains("2018A")) jecTag = "Autumn18_RunA_V8_DATA";
            if(currentFileName.Contains("2018B")) jecTag = "Autumn18_RunB_V8_DATA";
            if(currentFileName.Contains("2018C")) jecTag = "Autumn18_RunC_V8_DATA";
            if(currentFileName.Contains("2018D")) jecTag = "Autumn18_RunD_V8_DATA";
            if(currentFileName.Contains("2017B")) jecTag = "Fall17_17Nov2017B_V32_DATA";
            if(currentFileName.Contains("2017C")) jecTag = "Fall17_17Nov2017C_V32_DATA";
            if(currentFileName.Contains("2017D")) jecTag = "Fall17_17Nov2017DE_V32_DATA";
            if(currentFileName.Contains("2017E")) jecTag = "Fall17_17Nov2017DE_V32_DATA";
            if(currentFileName.Contains("2017F")) jecTag = "Fall17_17Nov2017F_V32_DATA";
            if(currentFileName.Contains("2016B")) jecTag = "Summer16_07Aug2017BCD_V11_DATA";
            if(currentFileName.Contains("2016C")) jecTag = "Summer16_07Aug2017BCD_V11_DATA";
            if(currentFileName.Contains("2016D")) jecTag = "Summer16_07Aug2017BCD_V11_DATA";
            if(currentFileName.Contains("2016E")) jecTag = "Summer16_07Aug2017EF_V11_DATA";
            if(currentFileName.Contains("2016F")) jecTag = "Summer16_07Aug2017EF_V11_DATA";
            if(currentFileName.Contains("2016G")) jecTag = "Summer16_07Aug2017GH_V11_DATA";
            if(currentFileName.Contains("2016H")) jecTag = "Summer16_07Aug2017GH_V11_DATA";
        }else{
            if (year==2016) jecTag = "Summer16_07Aug2017_V11_MC";
            if (year==2017) jecTag = "Fall17_17Nov2017_V32_MC";
            if (year==2018) jecTag = "Autumn18_V8_MC";
        }
        cout << "Using jecTag: " << jecTag << endl;
        jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/"+jecTag+"_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/"+jecTag+"_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/"+jecTag+"_L3Absolute_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/"+jecTag+"_L2L3Residual_AK4PFchs.txt");        
        jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);


        TTree *tree = (TTree*)file.Get("Events");
        if(fast) TTreeCache::SetLearnEntries(10);
        if(fast) tree->SetCacheSize(128*1024*1024);
        cms3.Init(tree);
    
        // Loop over Events in current file
        if( nEventsTotal >= nEventsChain ) continue;
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {
    
            // Get Event Content
            if( nEventsTotal >= nEventsChain ) continue;
            if(fast) tree->LoadTree(event);
            cms3.GetEntry(event);
            ++nEventsTotal;
    
            // reset output tree variables
            t.Reset();

            // // Progress
            // CMS3::progress( nEventsTotal, nEventsChain );

            if(nEventsTotal%100000==0)
                cout << nEventsTotal << " / " << nEventsChain << endl;

            t.evt_run = evt_run();
            t.evt_lumi = evt_lumiBlock();
            t.evt_event = evt_event();

            if(isData)
                t.isGolden = bool(goodrun(t.evt_run, t.evt_lumi));
            else
                t.isGolden = true;

            // t.HLT_SingleMu = 1;
            // t.HLT_SingleEl = 1;
            if(isData){
                t.HLT_SingleEl = 
                    passHLTTriggerPattern("HLT_Ele27_eta2p1_WPTight_Gsf_v") || // 2016
                    passHLTTriggerPattern("HLT_Ele32_eta2p1_WPTight_Gsf_v") || // 2016
                    passHLTTriggerPattern("HLT_Ele27_WPTight_Gsf_v") ||
                    passHLTTriggerPattern("HLT_Ele32_WPTight_Gsf_v") ||  // 2017,18
                    passHLTTriggerPattern("HLT_Ele35_WPTight_Gsf_v") ||  // 2017,18
                    passHLTTriggerPattern("HLT_Ele38_WPTight_Gsf_v") ||  // 2017,18
                    passHLTTriggerPattern("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v"); // 2016
                t.HLT_SingleMu = 
                    passHLTTriggerPattern("HLT_IsoMu17_eta2p1_v") ||
                    passHLTTriggerPattern("HLT_IsoMu20_v") || passHLTTriggerPattern("HLT_IsoMu20_eta2p1_v") ||
                    passHLTTriggerPattern("HLT_IsoTkMu20_v") || passHLTTriggerPattern("HLT_IsoTkMu20_eta2p1_v") ||
                    passHLTTriggerPattern("HLT_IsoMu24_v") || passHLTTriggerPattern("HLT_IsoTkMu24_v") || 
                    passHLTTriggerPattern("HLT_IsoMu24_eta2p1_v") ||
                    passHLTTriggerPattern("HLT_IsoMu27_v") || passHLTTriggerPattern("HLT_IsoTkMu27_v");
            }

            t.passFilters = 
                filt_ecalTP() &&
                filt_goodVertices() &&
                filt_hbheNoise() &&
                filt_hbheNoiseIso() &&
                (year==2016 || filt_ecalBadCalibFilterUpdate()) &&
                filt_BadPFMuonFilter() &&
                // filt_BadChargedCandidateFilter() &&
                filt_globalSuperTightHalo2016() &&
                (!isData || filt_eeBadSc());

            t.scale1fb = 1.0;
            if(!isData && !isSignal){
                const string cms3_version = evt_CMS3tag().at(0).Data();
                const string dataset_name = evt_dataset().at(0).Data();
                t.scale1fb = datasetInfoFromFile.getScale1fbFromFile(dataset_name, cms3_version);
                if(genps_weight() < 0.0)
                    t.scale1fb *= -1;
                t.scale1fb *= lumi;
            }

            if(!isData){
                LorentzVector p4hmeson, p4hphoton, p4lep, p4neutrino;
                uint ngen = genps_p4().size();
                for(uint i=0; i<ngen; i++){
                    if(!genps_isLastCopy().at(i))
                        continue;
                
                    // higgs
                    if(genps_id().at(i) == 25){
                        t.genH_pt = genps_p4().at(i).pt();
                        t.genH_eta = genps_p4().at(i).eta();
                        t.genH_phi = genps_p4().at(i).phi();
                    }

                    // higgs daughters
                    if(genps_id_mother().at(i) == 25){
                        if(genps_id().at(i) == 22){
                            if(genps_isHardProcess().at(i)){
                                p4hphoton = genps_p4().at(i);
                                t.genHPhoton_pt = genps_p4().at(i).pt();
                                t.genHPhoton_eta = genps_p4().at(i).eta();
                                t.genHPhoton_phi = genps_p4().at(i).phi();
                            }
                        }else if(genps_id().at(i) == 333 || genps_id().at(i) == 113){
                            p4hmeson = genps_p4().at(i);
                            t.genHMeson_pdgId = genps_id().at(i);
                            t.genHMeson_pt = genps_p4().at(i).pt();
                            t.genHMeson_eta = genps_p4().at(i).eta();
                            t.genHMeson_phi = genps_p4().at(i).phi();
                        }else{
                            cout << "WARNING: found higgs daughter that's not 22/113/333. run:lumi:evt " << evt_run() << ":" << evt_lumiBlock() << ":" << evt_event() << ", id=" << genps_id().at(i) << endl;
                        }
                    }

                    // higgs granddaughters from phi (333)
                    if(genps_id_mother().at(i) == 333 && genps_id_mother().at(genps_idx_mother().at(i)) == 25){
                        // only save K+/K-
                        if(genps_id().at(i) == 321){
                            t.genKp_pt = genps_p4().at(i).pt();
                            t.genKp_eta = genps_p4().at(i).eta();
                            t.genKp_phi = genps_p4().at(i).phi();                       
                        }
                        if(genps_id().at(i) == -321){
                            t.genKm_pt = genps_p4().at(i).pt();
                            t.genKm_eta = genps_p4().at(i).eta();
                            t.genKm_phi = genps_p4().at(i).phi();                       
                        }
                    }
                    
                    // W daughters
                    if(abs(genps_id_mother().at(i)) == 24 && 
                       genps_status().at(genps_idx_mother().at(i)) == 62 &&
                       genps_isLastCopy().at(genps_idx_mother().at(i)) ) {

                        if(abs(genps_id().at(i)) % 2 == 1){
                            p4lep = genps_p4().at(i);
                            t.genLep_pdgId = genps_id().at(i);
                            t.genLep_pt = genps_p4().at(i).pt();
                            t.genLep_eta = genps_p4().at(i).eta();
                            t.genLep_phi = genps_p4().at(i).phi();
                        }else{
                            p4neutrino = genps_p4().at(i);
                            t.genNeutrino_pdgId = genps_id().at(i);
                            t.genNeutrino_pt = genps_p4().at(i).pt();
                            t.genNeutrino_eta = genps_p4().at(i).eta();
                            t.genNeutrino_phi = genps_p4().at(i).phi();
                        }
                    }
                
                }

                if(t.genKm_pt>0 && t.genKp_pt>0){
                    t.genKK_dR = DeltaR(t.genKm_eta, t.genKp_eta, t.genKm_phi, t.genKp_phi);
                }

                if(t.genHPhoton_pt > 0 && t.genHMeson_pt > 0 && t.genLep_pt > 0 && t.genNeutrino_pt > 0){
                    VHAngles angles = getAngles(p4neutrino.pt(), p4neutrino.phi(), t.genLep_pdgId, lvconv(p4lep), lvconv(p4hmeson), lvconv(p4hphoton), p4neutrino.pz());
                    t.genvhAngles_cosThetaStar = angles.angles[0];
                    t.genvhAngles_cosTheta1 = angles.angles[1];
                    t.genvhAngles_cosTheta2 = angles.angles[2];
                    t.genvhAngles_phi = angles.angles[3];
                    t.genvhAngles_phi1 = angles.angles[4];
                    t.genvhAngles_m1 = angles.angles[5];
                    t.genvhAngles_m2 = angles.angles[6];                    
                }
            }
            
            // get MET info
            t.met_rawPt = evt_pfmet_raw();
            t.met_rawPhi = evt_pfmetPhi_raw();
            std::pair <float, float> t1met;
            t1met = getT1CHSMET_fromMINIAOD(jet_corrector_pfL1FastJetL2L3, 0, 0, false, false);
            t.met_pt = t1met.first;
            t.met_phi = t1met.second;
            
            t.nphoton = 0;
            LorentzVector p4_photon;
            for(uint ip=0; ip<photons_p4().size(); ip++){
                LorentzVector p4 = photons_p4().at(ip);
                if(p4.pt() < 20.0 || fabs(p4.eta()) > 2.5)
                    continue;
                if(!isMediumPhotonPOG_Fall17V2(ip))
                    continue;
                if(photons_recoChargedHadronIso().at(ip) / p4.pt() > 0.06)
                    continue;

                // check if overlapping a pt>10 electron
                bool overlapEl = false;
                for (unsigned int j = 0; j < els_p4().size(); j++) {
                    LorentzVector el_p4 = els_p4().at(j);
                    if (el_p4.pt() > 10) {
                        overlapEl = DeltaR(p4.eta(), el_p4.eta(), p4.phi(), el_p4.phi()) < 0.2;
                        if (overlapEl) break;
                    }
                }
                if(overlapEl)
                    continue;

                // count all photons that pass the above cuts, but only save the highest pT one
                t.nphoton++;
                if(p4.pt() < t.photon_pt)
                    continue;
                
                t.photon_pt = p4.pt();
                t.photon_eta = p4.eta();
                t.photon_phi = p4.phi();
                t.photon_isTightCutBased = isTightPhotonPOG_Fall17V2(ip);
                t.photon_relChHadIso = photons_recoChargedHadronIso().at(ip) / p4.pt();
                if(t.genHPhoton_pt > 0)
                    t.photon_dR_higgsPhoton = DeltaR(p4.eta(), t.genHPhoton_eta, p4.phi(), t.genHPhoton_phi);
                p4_photon = p4;
            }
            // get genmatching info
            t.photon_isMatchedGen = false;
            t.photon_drMinParton = 999;
            if(!isData && t.photon_pt > 0){
                int bestIdx = -1;
                float bestDR = 0.1;
                float bestEta=0, bestPhi=0;
                for(uint ig=0; ig < genps_p4().size(); ig++){
                    LorentzVector p4 = genps_p4().at(ig);
                    // match to gen photon
                    if(genps_id().at(ig) != 22 || genps_status().at(ig)!=1) continue;
                    if(fabs(genps_id_simplemother().at(ig)) > 22  && genps_id_simplemother().at(ig)!=2212) continue;
                    if(p4.pt() > 2*p4_photon.pt() || p4.pt() < 0.5*p4_photon.pt()) continue;
                    float dR = DeltaR(p4.eta(), p4_photon.eta(), p4.phi(), p4_photon.phi());
                    if(dR < bestDR){
                        bestIdx = ig;
                        bestDR = dR;
                        bestEta = p4.eta();
                        bestPhi = p4.phi();
                    }
                }
                if(bestIdx != -1){
                    t.photon_isMatchedGen = true;
                    float minDR = 999;
                    for(uint ig=0; ig < genps_p4().size(); ig++){
                        LorentzVector p4 = genps_p4().at(ig);
                        if (genps_status().at(ig) != 22 && genps_status().at(ig) != 23) continue;
                        if (fabs(genps_id().at(ig)) > 21) continue;
                        float dr = DeltaR(p4.eta(), bestEta, p4.phi(), bestPhi);
                        if (dr < minDR) minDR = dr;
                    }
                    t.photon_drMinParton = minDR;
                }
            }

            t.nlep = 0;
            LorentzVector p4_lepton;
            for(uint ie=0; ie<els_p4().size(); ie++){
                LorentzVector p4 = els_p4().at(ie);
                if(p4.pt() < 20 || fabs(p4.eta()) > 2.4)
                    continue;
                // id/iso
                if(!electronID(ie,id_level_t::HAD_medium_noiso_v5))
                    continue;
                if(elMiniRelIsoCMS3_EA(ie, gconf.ea_version) > 0.1)
                    continue;

                t.nlep++;

                if(p4.pt() < t.lepton_pt)
                    continue;
                
                p4_lepton = p4;
                t.lepton_pt = p4.pt();
                t.lepton_eta = p4.eta();
                t.lepton_phi = p4.phi();                
                t.lepton_pdgId = -11*els_charge().at(ie);
                t.lepton_miniRelIso = elMiniRelIsoCMS3_EA(ie, gconf.ea_version);
            }
            for(uint im=0; im<mus_p4().size(); im++){
                LorentzVector p4 = mus_p4().at(im);
                if(p4.pt() < 20 || fabs(p4.eta()) > 2.4)
                    continue;
                // id/iso
                // if(!muonID(im,id_level_t::HAD_tight_v4))
                if(!isMediumMuonPOG(im))
                    continue;
                if(muMiniRelIsoCMS3_EA(im, gconf.ea_version) > 0.2)
                    continue;

                t.nlep++;

                if(p4.pt() < t.lepton_pt)
                    continue;
                
                p4_lepton = p4;
                t.lepton_pt = p4.pt();
                t.lepton_eta = p4.eta();
                t.lepton_phi = p4.phi();
                t.lepton_pdgId = -13*mus_charge().at(im);
                t.lepton_miniRelIso = muMiniRelIsoCMS3_EA(im, gconf.ea_version);
            }

            if(t.nlep > 0 && t.nphoton > 0)
                t.lepton_photon_dR = DeltaR(p4_lepton.eta(), p4_photon.eta(), p4_lepton.phi(), p4_photon.phi());

            t.isHEM = false;
            if(year == 2018){
                if((isData && t.evt_run >= 319077) || (!isData && t.evt_event % 1961 < 1286)){
                    if(t.nphoton > 0 && t.photon_eta > -4.7 && t.photon_eta < -1.4 && t.photon_phi > -1.6 && t.photon_phi < -0.8)
                        t.isHEM = true;
                    if(t.nlep > 0 && t.lepton_eta > -4.7 && t.lepton_eta < -1.4 && t.lepton_phi > -1.6 && t.lepton_phi < -0.8)
                        t.isHEM = true;
                }
            }

            // loop over pfcands
            vector<uint> km_idxs, kp_idxs;
            for(uint i=0; i<pfcands_p4().size(); i++){
                if(pfcands_p4().at(i).pt() > 10.0 && fabs(pfcands_p4().at(i).eta()) < 2.4){
                    if(pfcands_particleId().at(i) == 211)
                        kp_idxs.push_back(i);
                    if(pfcands_particleId().at(i) == -211)
                        km_idxs.push_back(i);
                }
            }

            t.nmcand = 0;
            float best_phi_dm = 0.0, best_rho_dm = 0.0;
            for(uint ip=0; ip<kp_idxs.size(); ip++){
                for(uint im=0; im<km_idxs.size(); im++){
                    uint idxp = kp_idxs.at(ip);
                    uint idxm = km_idxs.at(im);
                    
                    // these work under the pion hypothesis, since particle flow assumes pion mass
                    LorentzVector p4m = pfcands_p4().at(idxm);
                    LorentzVector p4p = pfcands_p4().at(idxp);
                    // for the kaon hypothesis, need to modify the energy term to account for kaon mass
                    LorentzVector p4m_K, p4p_K;
                    p4m_K.SetPxPyPzE(p4m.px(), p4m.py(), p4m.pz(), sqrt(p4m.P2() + pow(0.493677,2)));
                    p4p_K.SetPxPyPzE(p4p.px(), p4p.py(), p4p.pz(), sqrt(p4p.P2() + pow(0.493677,2)));

                    LorentzVector p4 = p4m+p4p;
                    LorentzVector p4_K = p4m_K+p4p_K;

                    float dr = DeltaR(p4m.eta(), p4p.eta(), p4m.phi(), p4p.phi());
                    if(dr < 0.1){
                        // compute isolation. Make sure to subtract off the other track candidate
                        float isop = pfcands_trackIso().at(idxp);
                        float isom = pfcands_trackIso().at(idxm);
                        if(fabs(pfcands_dz().at(idxp)) < 0.1 || pfcands_fromPV().at(idxp) > 1)
                            isom -= p4p.pt();
                        if(fabs(pfcands_dz().at(idxm)) < 0.1 || pfcands_fromPV().at(idxm) > 1)
                            isop -= p4m.pt();

                        // don't consider if it isn't isolated
                        if(max(isom, isop) / p4.pt() > 0.06)
                            continue;
                        
                        t.nmcand++;
                        t.mcand_relIso->push_back(max(isom, isop) / p4.pt());

                        float mass = p4.M();
                        float mass_K = p4_K.M();
                        t.mcand_Km_pt->push_back(p4m.pt());
                        t.mcand_Kp_pt->push_back(p4p.pt());
                        t.mcand_mass_pion->push_back(mass);
                        t.mcand_mass_kaon->push_back(mass_K);
                        t.mcand_pt->push_back(p4.pt());
                        t.mcand_eta->push_back(p4.eta());
                        t.mcand_phi->push_back(p4.phi());
                        t.mcand_dR->push_back(dr);

                        bool fromPV = true;
                        // fromPV &= fabs(pfcands_dz().at(idxp)) < 0.1;
                        // fromPV &= fabs(pfcands_dz().at(idxm)) < 0.1;
                        // fromPV &= pfcands_fromPV().at(idxp) > 1;
                        // fromPV &= pfcands_fromPV().at(idxm) > 1;
                        fromPV &= (fabs(pfcands_dz().at(idxp)) < 0.1 || pfcands_fromPV().at(idxp) > 1);
                        fromPV &= (fabs(pfcands_dz().at(idxm)) < 0.1 || pfcands_fromPV().at(idxm) > 1);
                        t.mcand_fromPV->push_back(fromPV);

                        if(t.genHMeson_pt > 0)
                            t.mcand_dR_genMeson->push_back(DeltaR(p4.eta(), t.genHMeson_eta, p4.phi(), t.genHMeson_phi));
                        
                        if(t.best_phiCand_idx < 0 || fabs(mass_K - 1.0195) < best_phi_dm){
                            t.best_phiCand_idx = t.mcand_mass_kaon->size()-1;
                            best_phi_dm = fabs(mass_K - 1.0195);
                        }
                        if(t.best_rhoCand_idx < 0 || fabs(mass - 0.7755) < best_rho_dm){
                            t.best_rhoCand_idx = t.mcand_mass_pion->size()-1;
                            best_rho_dm = fabs(mass - 0.7755);
                        }

                        if(t.nphoton > 0){
                            t.mcand_photon_mass->push_back((p4+p4_photon).M());
                            t.mcand_photon_dR->push_back(DeltaR(p4.eta(), p4_photon.eta(), p4.phi(), p4_photon.phi()));
                        }else{
                            t.mcand_photon_mass->push_back(-999);
                            t.mcand_photon_dR->push_back(-999);
                        }

                        if(t.nphoton > 0 && t.nlep > 0){
                            VHAngles angles = getAngles(t.met_pt, t.met_phi, t.lepton_pdgId, lvconv(p4_lepton), lvconv(p4), lvconv(p4_photon));
                            t.vhAngles_cosThetaStar->push_back(angles.angles[0]);
                            t.vhAngles_cosTheta1->push_back(angles.angles[1]);
                            t.vhAngles_cosTheta2->push_back(angles.angles[2]);
                            t.vhAngles_phi->push_back(angles.angles[3]);
                            t.vhAngles_phi1->push_back(angles.angles[4]);
                            t.vhAngles_m1->push_back(angles.angles[5]);
                            t.vhAngles_m2->push_back(angles.angles[6]);
                        }
                    }
                }
            }

            t.Fill();
        }
  
        // Clean Up
        delete tree;
        file.Close();
    }
    if ( nEventsChain != nEventsTotal ) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  
    t.Write(fout);

    fout->Close();

    // return
    bmark->Stop("benchmark");
    cout << endl;
    cout << nEventsTotal << " Events Processed" << endl;
    cout << "------------------------------" << endl;
    cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
    cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
    cout << endl;
    delete bmark;
    return 0;
}
