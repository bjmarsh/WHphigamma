import scipy.stats
import numpy as np
import ROOT as r
r.gStyle.SetOptStat(0)

SEED = 1
np.random.seed(SEED)
r.gRandom.SetSeed(SEED)

r.gROOT.ProcessLine(".L crystalBall.C")

inject_signal = False

N_bkg = 40
x_scale = 80
N_sig = 5

bkg_fnc = r.TF1("bkg_fnc", r.Exp, 100, 180, 2)
bkg_fnc.SetParameter(0, N_bkg)
bkg_fnc.SetParameter(1, x_scale)

# sig_fnc = r.TF1("sig_fnc", r.CrystalBall, 100, 180, 5)
# sig_fnc.SetParameter(0, 1)
# sig_fnc.SetParameter(1, 125)
# sig_fnc.SetParameter(2, 2)
# sig_fnc.SetParameter(3, 2)
# sig_fnc.SetParameter(4, 1)
# fid = r.TFile("/nfs-7/userdata/bemarsh/rare-higgs/btree/v1/2018/wh.root")
# t = fid.Get("Events")
# h_sig =r.TH1D("h_sig","",160,100,180)
# t.Draw("mcand_photon_mass[best_phiCand_idx]>>h_sig", "fabs(mcand_mass_kaon[best_phiCand_idx]-1.02)<0.02 && lepton_pt>35", "goff")
# h_sig.Scale(1.0/ h_sig.GetMaximum())
# h_sig.Fit("sig_fnc", "N", "goff", 100, 180)
# sig_fnc.SetParameter(0, N_sig * sig_fnc.GetParameter(0) / sig_fnc.Integral(100,180))
# sig_mean = sig_fnc.GetParameter(1)
# sig_width = sig_fnc.GetParameter(2)
# sig_alpha = sig_fnc.GetParameter(3)
# sig_n = sig_fnc.GetParameter(4)

bkg_fnc.SetNpx(300)
# sig_fnc.SetNpx(300)

h_data = r.TH1D("h_data", ";m_{K^{+}K^{#minus}#gamma}", 80,100,180)

if inject_signal:
    Ns = [N_bkg,N_sig]
    fs = [bkg_fnc,sig_fnc]
else:
    Ns = [N_bkg]
    fs = [bkg_fnc]
data_points = []
for N,f in zip(Ns,fs):
    n = np.random.poisson(N)
    for i in range(n):
        m = f.GetRandom()
        data_points.append(m)
        h_data.Fill(m)

h_data.SetMarkerStyle(20)
h_data.SetMarkerColor(r.kBlack)
h_data.SetLineColor(r.kBlack)
h_data.GetXaxis().SetTitleSize(0.04)

h_blinded = h_data.Clone("h_blinded")
for i in range(21,31):
    h_blinded.SetBinContent(i, 0)

bkg_fits = []
bkg_fits.append(r.TF1("bkg_exp", r.Exp, 100, 180, 2))
bkg_fits[-1].SetParameter(0, N_bkg)
bkg_fits[-1].SetParameter(1, x_scale)
bkg_fits.append(r.TF1("bkg_lin", r.Quadratic, 100, 180, 3))
bkg_fits[-1].FixParameter(2, 0)
bkg_fits.append(r.TF1("bkg_quad", r.Quadratic, 100, 180, 3))
bkg_fits.append(r.TF1("bkg_pow1", r.Power, 100, 180, 4))
bkg_fits[-1].SetParameter(0,2e5)
bkg_fits[-1].SetParameter(1,-2)
bkg_fits[-1].FixParameter(2, 0)
bkg_fits[-1].FixParameter(3, 0)
bkg_fits.append(r.TF1("bkg_pow2", r.Power, 100, 180, 4))
bkg_fits[-1].SetParameter(0,2e5)
bkg_fits[-1].SetParameter(1,-2)
bkg_fits[-1].SetParameter(2,2e5)
bkg_fits[-1].SetParameter(3,-2)

r.REJECT=True
for fit in bkg_fits[::-1]:
    h_blinded.Fit(fit, "NLQ", "goff", 100, 180)
r.REJECT=False

print [h_blinded.GetBinContent(i) for i in range(h_blinded.GetNbinsX())]
for f in [bkg_fnc]+bkg_fits:
    NLL = -sum(scipy.stats.poisson.logpmf(h_blinded.GetBinContent(i), f.Eval(h_blinded.GetBinCenter(i))) for i in range(h_blinded.GetNbinsX()))
    print f.GetName(), f.Integral(120,130), NLL

gr_fit = r.TGraphErrors()
gr_err1 = r.TGraphErrors()
gr_err2 = r.TGraphErrors()
for i,x in enumerate(np.linspace(100, 180, 200)):
    gr_err1.SetPoint(i, x, bkg_fits[0].Eval(x))
    gr_err2.SetPoint(i, x, bkg_fits[0].Eval(x))
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err1, 0.68)
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err2, 0.95)

c2 = r.TCanvas("c2","c2", 700, 600)

cols = [r.kRed, r.kBlue, r.kGreen+2, r.kMagenta+2, r.kAzure+8]
for fit,col in zip(bkg_fits,cols):
    fit.SetLineWidth(2)
    fit.SetLineStyle(1)
    fit.SetLineColor(col)

bkg_fnc.SetLineColor(r.kBlack)
bkg_fnc.SetLineStyle(7)
gr_err1.SetFillColor(r.kGreen+1)
gr_err2.SetFillColor(r.kOrange-2)


h_blinded.Draw("PE")
# gr_err2.Draw("3 SAME")
# gr_err1.Draw("3 SAME")
for fit in bkg_fits:
    fit.Draw("L SAME")
bkg_fnc.Draw("L SAME")
h_blinded.Draw("PE SAME")

raw_input()


