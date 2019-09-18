import ROOT as r
r.gROOT.ProcessLine(".L crystalBall.C")

fid = r.TFile("/nfs-7/userdata/bemarsh/rare-higgs/btree/v1/2018/wh.root")
t = fid.Get("Events")

h_sig = r.TH1D("h_sig","",160,100,180)
t.Draw("mcand_photon_mass[best_phiCand_idx]>>h_sig", "fabs(mcand_mass_kaon[best_phiCand_idx]-1.02)<0.02 && lepton_pt>35", "goff")

h_sig.Scale(1.0 / h_sig.GetMaximum())

f = r.TF1("sig_cb", r.CrystalBall, 100, 180, 5)
f.SetParameter(0, 1)
f.SetParameter(1, 125)
f.SetParameter(2, 2)
f.SetParameter(3, 2)
f.SetParameter(4, 1)
f.SetNpx(300)

h_sig.Fit("sig_cb", "N", "goff", 100, 180)

h_sig.Draw()
f.Draw("SAME")

raw_input()
