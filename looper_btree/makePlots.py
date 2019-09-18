import os
import ROOT as r
import pyRootPlotMaker as ppm
r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L InvMass.C")

DOPHI = True
best_idx = "best_phiCand_idx" if DOPHI else "best_rhoCand_idx"

indir = "/nfs-7/userdata/bemarsh/rare-higgs/btree/v1/YEAR"
years = ["2016","2017","2018"]
# years = ["2018"]

sig = "wh"
bkgs = ["dyjetsll50_m", "wg","wjets","ttsl_t","ttgamma", "ttdl"]
data = ["data*SingleMuon","data*EGamma","data*SingleElectron"]
# data = ["data*SingleMuon","data*EGamma"]

# baseline = r.TCut("photon_pt>20 && lepton_pt>35 && nmcand>0 && isGold && !isHEM && passFilters && (recoHiggs_mass<120||recoHiggs_mass>130||genHiggs_pt>0) && recoWLepton_nLep==1")

bdt_cut = 0.9
baseline = r.TCut("photon_pt>20 && ((abs(lepton_pdgId)==11 && lepton_pt>35)||(abs(lepton_pdgId)==13 && lepton_pt>30)) && nmcand>0 && isGolden && passFilters && !isHEM && (mcand_photon_mass[{0}]<120||mcand_photon_mass[{0}]>130||genH_pt>0) && nlep==1 && mcand_photon_dR[{0}]>0.1".format(best_idx))
# baseline += r.TCut("fabs(lepton_eta)<1.5")
if bdt_cut > 0:
    baseline += r.TCut("bdt_score > {0}".format(bdt_cut))

# isLog, binning
plots = {
    "recoWLepton_pt" : ("lepton_pt", 0, (25,0,200),0.5),
    "recoWLepton_eta" : ("lepton_eta", 0, (30,-3,3),0.5),
    "recoWLepton_phi" : ("lepton_phi", 0, (30,-3.1416,3.1416),0.5),
    "recoGamma_pt" : ("photon_pt", 0, (25,0,200),0.5),
    "recoGamma_eta" : ("photon_eta", 0, (30,-3,3),0.5),
    "recoGamma_phi" : ("photon_phi", 0, (30,-3.1416,3.1416),0.5),
    "recoPhi_pt" : ("mcand_pt[{0}]".format(best_idx), 0, (25,0,200),0.5),
    "recoPhi_eta" : ("mcand_eta[{0}]".format(best_idx), 0, (30,-3,3),0.5),
    "recoPhi_phi" : ("mcand_phi[{0}]".format(best_idx), 0, (30,-3.1416,3.1416),0.5),
    "recoKpKm_dR" : ("mcand_dR[{0}]".format(best_idx), 0, (50,0,0.1),0.4),
    "recoPhi_relIso" : ("mcand_relIso[{0}]".format(best_idx), 0, (50,0,0.1),0.5),
    "recoPhi_mass" : ("mcand_mass_kaon[{0}]".format(best_idx), 0, (25,0.98,1.05),0.15),
    "recoRho_mass" : ("mcand_mass_pion[{0}]".format(best_idx), 0, (25,0.5,1.0),0.15),
    "recoHiggs_massAll" : ("mcand_photon_mass[{0}]".format(best_idx), 0, (36,0,180),0.1),
    "recoHiggs_mass" : ("mcand_photon_mass[{0}]".format(best_idx), 0, (18,90,180),0.1),
    "recoHiggs_mass_fitWindow" : ("mcand_photon_mass[{0}]".format(best_idx), 0, (80,100,180),0.1),
    "recoGammaWLepton_dR" : ("lepton_photon_dR", 0, (30,0,4),0.6),
    "recoPhiGamma_dR" : ("mcand_photon_dR[{0}]".format(best_idx), 0, (30,0,4),0.8),
    "recoMagAng_cosThetaStar" : ("vhAngles_cosThetaStar[{0}]".format(best_idx), 0, (30,-1,1), 1.0),
    "recoMagAng_cosTheta1" : ("vhAngles_cosTheta1[{0}]".format(best_idx), 0, (30,-1,1), 1.0),
    "recoMagAng_cosTheta2" : ("vhAngles_cosTheta2[{0}]".format(best_idx), 0, (30,-1,1), 1.0),
    "recoMagAng_Phi" : ("vhAngles_phi[{0}]".format(best_idx), 0, (30,-3.1416,3.1416), 1.0),
    "recoMagAng_Phi1" : ("vhAngles_phi1[{0}]".format(best_idx), 0, (30,-3.1416,3.1416), 1.0),
    "recoMagAng_m1" : ("vhAngles_m1[{0}]".format(best_idx), 0, (50,0,500), 0.5),
    "recoMet_pt" : ("met_pt", 0, (50,0,200), 0.5),
    "recoLepPho_invMass" : ("InvMass(lepton_pt,lepton_eta,lepton_phi,0,photon_pt,photon_eta,photon_phi,0)", 0, (50,0,300), 0.5),
    "bdt_score" : ("bdt_score", 1, (20,0,1), 1.0),
}

hists = {}
for proc in [sig]+bkgs+data:
    print proc
    c = r.TChain("Events")
    for year in years:
        c.Add("{0}/{1}*.root".format(indir.replace("YEAR",year), proc))
    hists[proc] = {}
    extra = r.TCut()
    if proc in ["wg", "ttg"]:
        extra = r.TCut("photon_isMatchedGen && photon_drMinParton>0.2")
    if proc in ["wjets", "ttsl_t"]:
        extra = r.TCut("(!photon_isMatchedGen || photon_drMinParton<0.2)")
    if proc in ["wh"]:
        if DOPHI:
            extra = r.TCut("genHMeson_pdgId == 333")
        else:
            extra = r.TCut("genHMeson_pdgId == 113")
    if "data" in proc and "SingleMuon" in proc:
        extra = r.TCut("HLT_SingleMu && abs(lepton_pdgId)==13")
    if "data" in proc and ("EGamma" in proc or "SingleElectron" in proc):
        extra = r.TCut("HLT_SingleEl && abs(lepton_pdgId)==11")
    weight = "scale1fb"
    if proc in ["ttsl_t","wg"]:
        weight = "scale1fb*(35.9+41.5+59.8)/(35.9*2+41.5+59.8)"
    if proc in ["dyjetsll50_m"]:
        weight = "scale1fb*(35.9+41.5+59.8)/(35.9*2+41.5*2+59.8)"
    if proc in ["ttdl"]:
        weight = "scale1fb*(35.9+41.5+59.8)/(35.9*2+41.5+59.8*0)"
    # if proc == "wh":
    #     weight = 1.380 * 1000 / 1000000 * 137.2
    for plot in plots:        
        if "fitWindow" in plot:
            extra += r.TCut("fabs(InvMass(lepton_pt,lepton_eta,lepton_phi,0,photon_pt,photon_eta,photon_phi,0)-87.5)>7.5")
            extra += r.TCut("mcand_mass_kaon[best_phiCand_idx]>1.00 && mcand_mass_kaon[best_phiCand_idx]<1.04")
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
    saveAs = "/home/users/bemarsh/public_html/WHphigamma/crplots/btree/v1{0}{1}/{2}".format("_phi" if DOPHI else "_rho", "" if bdt_cut==0.0 else "_bdt_{0}".format(str(bdt_cut).replace(".","p")), plot)
    colors = [r.kAzure-2, r.kOrange-2, r.kSpring-6, r.kRed+1, r.kGray+1, r.kViolet+5]
    os.system("mkdir -p "+os.path.split(saveAs)[0])
    if "Higgs_mass" in plot:
        sigScale = None
        h_sig_vec[0].Scale(1.380 * 1000 * 0.216 / 500000 * 137.2 * 5e-3)
    userMin, userMax = None, None
    if "fitWindow" in plot:
        userMax = 10
    for ext in [".png",".pdf"]:
        ppm.plotDataMC(h_bkg_vec, bkgs, h_data, h_sig_vec=h_sig_vec, sig_names=[sig], isLog=isLog, 
                       scaleSigToMC=sigScale, scaleMCtoData=True, doOverflow=False, doBkgError=True, 
                       legCoords=(0.60,0.72,0.87,0.89), legNCol=2, lumi=137, customColors=colors,
                       userMin=userMin, userMax=userMax,
                       xAxisTitle=plot, subLegText=subLegText, doPause=False, saveAs=saveAs+ext)



