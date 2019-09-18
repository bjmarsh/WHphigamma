import os,sys
import dis_client
import glob

TAG = "v1"

datasets = {
    2016: [
        ("/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM", "wg_ext1"),
        ("/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM", "wg_ext3"),
        ("/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM", "wjets"),
        ("/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "ttgjets_ext1"),
        ("/TTGamma_SingleLeptFromT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "ttgamma_top"),
        ("/TTGamma_SingleLeptFromTbar_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "ttgamma_tbar"),
        ("/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "ttsl_top"),
        ("/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "ttsl_top_ext1"),
        ("/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "ttsl_tbar"),
        ("/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "ttsl_tbar_ext1"),
        ("/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "ttjets_incl"),
        ("/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", "tt_powheg"),
        ("/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "ttjets_amcatnlo"),
        ("/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM", "dyjetsll50_amcatnlo_ext2"),
        ("/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM", "dyjetsll50_mglo_ext2"),
        ("/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "dyjetsll50_mglo_ext1"),
        ("/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "dyjetsll10_amcatnlo_ext1"),
        ("/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "ttdl"),
        ("/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "ttdl_ext1"),

        ("/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD" , "data_Run2016B_SingleMuon_17Jul2018_ver2-v1"),
        ("/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD" , "data_Run2016C_SingleMuon_17Jul2018-v1"),
        ("/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD" , "data_Run2016D_SingleMuon_17Jul2018-v1"),
        ("/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD" , "data_Run2016E_SingleMuon_17Jul2018-v1"),
        ("/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD" , "data_Run2016F_SingleMuon_17Jul2018-v1"),
        ("/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD" , "data_Run2016G_SingleMuon_17Jul2018-v1"),
        ("/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD" , "data_Run2016H_SingleMuon_17Jul2018-v1"),
        
        ("/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD" , "data_Run2016B_SingleElectron_17Jul2018_ver2-v1"),
        ("/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD" , "data_Run2016C_SingleElectron_17Jul2018-v1"),
        ("/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD" , "data_Run2016D_SingleElectron_17Jul2018-v1"),
        ("/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD" , "data_Run2016E_SingleElectron_17Jul2018-v1"),
        ("/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD" , "data_Run2016F_SingleElectron_17Jul2018-v1"),
        ("/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD" , "data_Run2016G_SingleElectron_17Jul2018-v1"),
        ("/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD" , "data_Run2016H_SingleElectron_17Jul2018-v1"),

        ],

    2017: [
        ("/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "wg"),
        ("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM", "wjets"),
        ("/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttgjets"),
        ("/TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttgamma_top"),
        ("/TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8/RunIIFall17MiniAOD-PU2017_94X_mc2017_realistic_v11-v1/MINIAODSIM", "ttgamma_tbar"),
        ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_top"),
        ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_tbar"),
        ("/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_incl"),
        ("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_amcatnlo"),
        ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_powheg_newpmx"),
        ("/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttsl_powheg"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "dyjetsll50_mglo"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM", "dyjetsll50_mglo_ext1"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "dyjetsll50_amcatnlo"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM", "dyjetsll50_amcatnlo_ext1"),
        ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "dyjetsll10_mglo"),
        ("/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "ttdl"),

        ("/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"     , "data_Run2017B_SingleMuon_31Mar2018"),
        ("/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD" , "data_Run2017B_SingleElectron_31Mar2018"),
        ("/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"     , "data_Run2017C_SingleMuon_31Mar2018"),
        ("/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD" , "data_Run2017C_SingleElectron_31Mar2018"),
        ("/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"     , "data_Run2017D_SingleMuon_31Mar2018"),
        ("/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD" , "data_Run2017D_SingleElectron_31Mar2018"),
        ("/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"     , "data_Run2017E_SingleMuon_31Mar2018"),
        ("/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD" , "data_Run2017E_SingleElectron_31Mar2018"),
        ("/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"     , "data_Run2017F_SingleMuon_31Mar2018"),
        ("/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD" , "data_Run2017F_SingleElectron_31Mar2018"),
        ],

    2018: [
        ("/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "wg"),
        ("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", "wjets"),
        ("/TTGamma_SingleLeptFromT_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "ttgamma_top"),
        ("/TTGamma_SingleLeptFromTbar_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "ttgamma_tbar"),
        ("/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "ttgjets"),
        ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "ttsl_top"),
        ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "ttsl_tbar"),
        ("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "ttjets_amcatnlo"),
        ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "ttsl_powheg"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "dyjetsll50_amcatnlo"),
        ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "dyjetsll50_mglo"),
        ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", "dyjetsll10_mglo"),
        ("/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "ttdl"),
        ("/WH_HtoRhoGammaPhiGamma/privateMC_102x/MINIAOD", "wh"),

        ("/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD" , "data_Run2018A_SingleMuon_17Sep2018"),
        ("/EGamma/Run2018A-17Sep2018-v2/MINIAOD"     , "data_Run2018A_EGamma_17Sep2018"),
        ("/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD" , "data_Run2018B_SingleMuon_17Sep2018"),
        ("/EGamma/Run2018B-17Sep2018-v1/MINIAOD"     , "data_Run2018B_EGamma_17Sep2018"),
        ("/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD" , "data_Run2018C_SingleMuon_17Sep2018"),
        ("/EGamma/Run2018C-17Sep2018-v1/MINIAOD"     , "data_Run2018C_EGamma_17Sep2018"),
        ("/SingleMuon/Run2018D-PromptReco-v2/MINIAOD", "data_Run2018D_SingleMuon_PromptReco"),
        ("/EGamma/Run2018D-22Jan2019-v2/MINIAOD"     , "data_Run2018D_EGamma_22Jan2019"),
        ],
}


os.system("mkdir -p configs/{0}".format(TAG))
        
for year in datasets:
    print year
    for ds,short in datasets[year]:
        if "WH_HtoRhoGammaPhiGamma" in ds:
            loc = "/hadoop/cms/store/group/snt/run2_mc2018_private/WH_HtoRhoGammaPhiGamma_privateMC_102x_MINIAOD_v1"
        else:
            info = dis_client.query(ds, "snt")["response"]["payload"]
            loc = info[0]["location"]
        print loc

        outdir = "/hadoop/cms/store/user/bemarsh/WHphigamma/{0}/{1}/{2}".format(TAG, year, short)

        fout = open("configs/{0}/config_{1}_{2}.cmd".format(TAG,year,short), 'w')
        fout.write("""
universe=vanilla
when_to_transfer_output = ON_EXIT
+DESIRED_Sites="T2_US_UCSD"
+remote_DESIRED_Sites="T2_US_UCSD"
+Owner = undefined
log=logs/condor_submit.log
output=logs/{0}/{1}/1e.$(Cluster).$(Process).out
error =logs/{0}/{1}/1e.$(Cluster).$(Process).err
notification=Never
x509userproxy=/tmp/x509up_u31592

executable=wrapper.sh
transfer_input_files=input.tar.xz
transfer_executable=True

""".format(year, short))

        os.system("mkdir -p logs/{0}/{1}".format(year, short))

        for f in glob.glob(loc+"/*.root"):
            bn = f.split("/")[-1].split(".")[0]
            fout.write("arguments={0} {1} {2}\nqueue\n\n".format(bn, loc, outdir))

        fout.close()

