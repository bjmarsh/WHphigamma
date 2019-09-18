{
    gSystem->Load("CORE/CMS3_CORE.so");
    // gSystem->Load("CMS3_CORE.so");
    gROOT->ProcessLine(".L CORE/Tools/utils.cc+");
    gROOT->ProcessLine(".L WHTree/WHTree.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");
    
    TChain *ch = new TChain("Events"); 

    // 80X
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/ZJetsToNuNu_HT-600To800_13TeV-madgraph_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/ZJetsToNuNu_HT-800To1200_13TeV-madgraph_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v3/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/V08-00-01/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v2/V08-00-01/merged_ntuple_*.root");

    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/V08-00-01/merged_ntuple_*1.root");

    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/V08-00-01/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv1/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v2/V08-00-01/merged_ntuple_*1.root");

    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/V08-00-05/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_*.root");

    //// 76X
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_76MiniAODv2/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/V07-06-03_MC/merged_ntuple_*.root");

    //// 74X
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");

    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-Asympt25ns_74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-Asympt25ns_74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-Asympt25ns_74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/V07-04-11/merged_ntuple_*.root");

    // ch->Add("/hadoop/cms/store/group/snt/run2_25ns_MiniAODv2/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/V07-04-11/merged_ntuple_*");

    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");
    // ch->Add("/hadoop/cms/store/group/snt/run2_mc2017/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2_MINIAODSIM_CMS4_V09-04-19/merged_ntuple_*1.root");

    // ScanChain(ch, "wjets_ht100to200_94x");
    // ScanChain(ch, "wjets_ht200to400_94x");
    // ScanChain(ch, "wjets_ht400to600_94x");
    // ScanChain(ch, "wjets_ht600to800_94x");
    // ScanChain(ch, "wjets_ht800to1200_94x");
    // ScanChain(ch, "wjets_ht1200to2500_94x");
    // ScanChain(ch, "wjets_ht2500toInf_94x");
    // ScanChain(ch, "wjets_incl_94x");
            
    // ScanChain(ch, "zinv_ht600to800_80x"); 
    // ScanChain(ch, "zinv_ht800to1200_80x"); 
    // ScanChain(ch, "zinv_ht1200to2500_80x"); 
    // ScanChain(ch, "zinv_ht2500toInf_80x"); 
    // ScanChain(ch, "zinv_ht600toInf_74x"); 
}
