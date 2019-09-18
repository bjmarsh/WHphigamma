import ROOT as r
import pyRootPlotMaker as ppm
r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(1)

r.gROOT.ProcessLine(".L ../looper_btree/InvMass.C")

indir = "/nfs-7/userdata/jguiang/rare-higgs/YEAR/v4-0-0/"
years = ["2016","2017","2018"]

sig = "wh"
bkgs = ["dyjets", "wg", "wjets","ttjets","ttg"]
data = ["data*SingleMuon","data*EGamma","data*SingleElectron"]
# data = ["data*SingleMuon","data*EGamma"]

baseline = r.TCut("recoGamma_pt>0 && recoWLepton_pt>35 && recoPhi_pt>0 && isGold && !isHEM && passFilters && (recoHiggs_mass<120||recoHiggs_mass>130||genHiggs_pt>0) && recoWLepton_nLep==1 && recoWLepton_nLepVeto==1")

# isLog, binning
plots = {
    "recoWLepton_pt" : (None, 0, (25,0,200),0.5),
    "recoWLepton_eta" : (None, 0, (30,-3,3),0.5),
    "recoWLepton_phi" : (None, 0, (30,-3.1416,3.1416),0.5),
    "recoGamma_pt" : (None, 0, (25,0,200),0.5),
    "recoGamma_eta" : (None, 0, (30,-3,3),0.5),
    "recoGamma_phi" : (None, 0, (30,-3.1416,3.1416),0.5),
    "recoPhi_pt" : (None, 0, (25,0,200),0.5),
    "recoPhi_eta" : (None, 0, (30,-3,3),0.5),
    "recoPhi_phi" : (None, 0, (30,-3.1416,3.1416),0.5),
    "recoKpKm_dR" : (None, 0, (50,0,0.1),0.4),
    "recoPhi_relIso" : (None, 0, (50,0,0.1),0.5),
    "recoPhi_mass" : ("recoPhi_mass", 0, (25,0.98,1.05),0.15),
    "recoHiggs_massAll" : ("recoHiggs_mass", 0, (36,0,180),0.1),
    "recoHiggs_mass" : ("recoHiggs_mass", 0, (18,90,180),0.1),
    "recoGammaWLepton_dR" : (None, 0, (30,0,4),0.6),
    "recoPhiGamma_dR" : (None, 0, (30,0,4),0.8),
    "recoMagAng_cosThetaStar" : (None, 0, (30,-1,1), 1.0),
    "recoMagAng_cosTheta1" : (None, 0, (30,-1,1), 1.0),
    "recoMagAng_cosTheta2" : (None, 0, (30,-1,1), 1.0),
    "recoMagAng_Phi" : (None, 0, (30,-3.1416,3.1416), 1.0),
    "recoMagAng_Phi1" : (None, 0, (30,-3.1416,3.1416), 1.0),
    "recoMagAng_m1" : ("recoMagAng_m1", 0, (50,0,500), 0.5),
    "recoMet_pt" : ("met_pt", 0, (50,0,200), 0.5),
    "recoMet_phi" : ("met_phi", 0, (30,-3.1416,3.1416), 0.5),
    "recoLepPhi_invMass" : ("InvMass(recoWLepton_pt,recoWLepton_eta,recoWLepton_phi,0,recoGamma_pt,recoGamma_eta,recoGamma_phi,0)", 0, (50,0,300), 0.5),    
}

hists = {}
for proc in [sig]+bkgs+data:
    print proc
    c = r.TChain("tree")
    for year in years:
        c.Add("{0}/output_{1}*.root".format(indir.replace("YEAR",year), proc))
    hists[proc] = {}
    extra = r.TCut()
    if proc in ["wg", "ttg"]:
        extra = r.TCut("genRecoGamma_isMatch && minGammaParton_dR>0.2")
    if proc in ["wjets", "ttsl"]:
        extra = r.TCut("(!genRecoGamma_isMatch || minGammaParton_dR<0.2)")
    if proc in ["wh"]:
        extra = r.TCut("genKp_pt>0")
    if "data" in proc and "SingleMuon" in proc:
        extra = r.TCut("HLT_singleMu && abs(recoWLepton_id)==13")
    if "data" in proc and ("EGamma" in proc or "SingleElectron" in proc):
        extra = r.TCut("HLT_singleEl && abs(recoWLepton_id)==11")
    weight = "scale1fb"
    if proc in ["ttjets","wg"]:
        weight = "scale1fb*(35.9+41.5+59.8)/(35.9*2+41.5+59.8)"
    if proc in ["dyjets"]:
        weight = "scale1fb*(35.9+41.5+59.8)/(35.9*2+41.5*2+59.8)"
    for plot in plots:
        toplot, _, binning, _ = plots[plot]
        hists[proc][plot] = r.TH1D("h_{0}_{1}".format(proc, plot), plot, *binning)
        toplot = toplot if toplot else plot
        c.Draw("{0}>>h_{1}_{2}".format(toplot, proc, plot), "({0})*({1})".format(str(baseline+extra),weight), "goff")



for plot in plots:
    h_data = hists[data[0]][plot]
    for p in data[1:]:
        h_data.Add(hists[p][plot])
    toplot, isLog, binning, sigScale = plots[plot]
    h_bkg_vec = [hists[proc][plot] for proc in bkgs]
    h_sig_vec = [hists[sig][plot]]
    subLegText = "MC scaled by {datamcsf} #pm {datamcsferr}"
    saveAs = "/home/users/bemarsh/public_html/WHphigamma/crplots/v4/{0}".format(plot)
    for ext in [".png",".pdf"]:
        ppm.plotDataMC(h_bkg_vec, bkgs, h_data, h_sig_vec=h_sig_vec, sig_names=[sig], isLog=isLog, 
                       scaleSigToMC=sigScale, scaleMCtoData=True, doOverflow=False, doBkgError=True, 
                       legCoords=(0.60,0.72,0.87,0.89), legNCol=2, lumi=137,
                       xAxisTitle=plot, subLegText=subLegText, doPause=False, saveAs=saveAs+ext)



