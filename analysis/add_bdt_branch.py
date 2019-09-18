import os,sys
import numpy as np
import ROOT as r
r.gROOT.ProcessLine(".L bdt_phi.C+")

fnames = sys.argv[1:]

for fname in fnames:
    print fname
    fid = r.TFile(fname, "UPDATE")
    
    torig = fid.Get("Events")
    torig.SetBranchStatus("*", 1)
    if torig.GetListOfBranches().FindObject("bdt_score"):
        torig.SetBranchStatus("bdt*", 0)

    t = torig.CloneTree(-1, "fast")

    np_bdt = np.array([0.0], dtype=np.float32)
    b_bdt = t.Branch("bdt_score", np_bdt, "bdt_score/F")

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        j = t.best_phiCand_idx
        s = r.get_prediction(
            t.photon_pt / t.mcand_photon_mass[j],
            t.photon_eta,
            t.photon_phi,
            t.photon_isTightCutBased,
            t.photon_relChHadIso,
            t.mcand_mass_kaon[j],
            t.mcand_pt[j] / t.mcand_photon_mass[j],
            t.mcand_eta[j],
            t.mcand_phi[j],
            t.mcand_dR[j],
            t.mcand_relIso[j],
            t.mcand_photon_dR[j],
            t.lepton_pt,
            t.lepton_eta,
            t.lepton_phi,
            t.lepton_miniRelIso,
            t.lepton_photon_dR,
            t.met_pt,
            t.met_phi,
            t.vhAngles_cosThetaStar[j],
            t.vhAngles_cosTheta1[j],
            t.vhAngles_cosTheta2[j],
            t.vhAngles_phi[j],
            t.vhAngles_phi1[j],
        )

        np_bdt[0] = s
        b_bdt.Fill()

    t.Write("Events", r.TObject.kWriteDelete)
    fid.Close()
