from math import exp, sqrt
import numpy as np
import ROOT as r
r.gStyle.SetOptStat(0)

r.gRandom.SetSeed(42)

r.gROOT.ProcessLine(".L ../crystalBall.C")

inject_signal = False

N_bkg = 38
x_scale = 60
N_sig = 5

bkg_fnc = r.TF1("bkg_fnc", r.Exp, 100, 180, 2)
bkg_fnc.SetParameter(0, N_bkg)
bkg_fnc.SetParameter(1, x_scale)

sig_fnc = r.TF1("sig_fnc", r.CrystalBall, 100, 180, 5)
sig_fnc.SetParameter(0, 1)
sig_fnc.SetParameter(1, 125)
sig_fnc.SetParameter(2, 2)
sig_fnc.SetParameter(3, 2)
sig_fnc.SetParameter(4, 1)
fid = r.TFile("/nfs-7/userdata/bemarsh/rare-higgs/btree/v1/2018/wh.root")
t = fid.Get("Events")
h_sig =r.TH1D("h_sig","",160,100,180)
t.Draw("mcand_photon_mass[best_phiCand_idx]>>h_sig", "fabs(mcand_mass_kaon[best_phiCand_idx]-1.02)<0.02 && lepton_pt>35", "goff")
h_sig.Scale(1.0/ h_sig.GetMaximum())
h_sig.Fit("sig_fnc", "N", "goff", 100, 180)
sig_fnc.SetParameter(0, N_sig * sig_fnc.GetParameter(0) / sig_fnc.Integral(100,180))
sig_mean = sig_fnc.GetParameter(1)
sig_width = sig_fnc.GetParameter(2)
sig_alpha = sig_fnc.GetParameter(3)
sig_n = sig_fnc.GetParameter(4)

#tot_fnc = r.TF1("tot_fnc", "bkg_fnc+sig_fnc", 100, 180)

bkg_fnc.SetNpx(300)
sig_fnc.SetNpx(300)
#tot_fnc.SetNpx(300)

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

bkg_fnc.SetLineWidth(3)
bkg_fnc.SetLineColor(r.kRed)
sig_fnc.SetLineWidth(3)
sig_fnc.SetLineColor(r.kBlue)
# tot_fnc.SetLineWidth(3)
# tot_fnc.SetLineColor(r.kBlue)

# c1 = r.TCanvas("c1","c1", 700, 600)

# h_data.Draw("PE")
# # tot_fnc.Draw("SAME")
# bkg_fnc.Draw("SAME")
# h_data.Draw("PE SAME")

# raw_input()

h_blinded = h_data.Clone("h_blinded")
for i in range(21,31):
    h_blinded.SetBinContent(i, 0)
    h_blinded.SetBinError(i, 100000)

bkg_fit = r.TF1("bkg_fit", r.Exp, 100, 180, 2)
bkg_fit.SetParameter(0, N_bkg)
bkg_fit.SetParameter(1, x_scale)

r.REJECT = True
result_ptr = h_blinded.Fit(bkg_fit, "NL", "goff", 100, 180)
fit_result = result_ptr.Get()
r.REJECT = False

gr_fit = r.TGraphErrors()
gr_err1 = r.TGraphErrors()
gr_err2 = r.TGraphErrors()
for i,x in enumerate(np.linspace(100, 180, 200)):
    gr_err1.SetPoint(i, x, bkg_fit.Eval(x))
    gr_err2.SetPoint(i, x, bkg_fit.Eval(x))
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err1, 0.68)
r.TVirtualFitter.GetFitter().GetConfidenceIntervals(gr_err2, 0.95)

c2 = r.TCanvas("c2","c2", 700, 600)

bkg_fit.SetLineWidth(2)
bkg_fit.SetLineStyle(1)
bkg_fit.SetLineColor(r.kRed)
gr_err1.SetFillColor(r.kGreen+1)
gr_err2.SetFillColor(r.kOrange-2)

bkg_fnc.SetLineColor(r.kBlack)
bkg_fnc.SetLineStyle(7)

h_data.Draw("PE")
gr_err2.Draw("3 SAME")
gr_err1.Draw("3 SAME")
sig_fnc.Draw("L SAME")
bkg_fit.Draw("L SAME")
bkg_fnc.Draw("L SAME")
h_data.Draw("PE SAME")

raw_input()


ws = r.RooWorkspace('w', "test_workspace")

M = r.RooRealVar("M","M",100,180)
data = r.RooDataSet("data_obs","data_obs",r.RooArgSet(M))

for i in range(len(data_points)):
    M.setVal(data_points[i])
    data.add(r.RooArgSet(M))
data.Print()

M.setBins(80)
data_binned = r.RooDataHist("data_obs_binned","data_obs_binned", r.RooArgSet(M), data)

x_scale = r.RooRealVar("x_scale","x_scale",-1,-0.01)
print "BLAH", -1.0 / bkg_fit.GetParameter(1)
x_scale.setVal(-1.0 / bkg_fit.GetParameter(1))
bkg_exp = r.RooExponential("bkg_exp", "bkg_exp", M, x_scale)

bkg_norm = bkg_fit.Integral(100,180)
bkg_exp_norm = r.RooRealVar("bkg_exp_norm","bkg_exp_norm", 0, 2*bkg_norm)
bkg_exp_norm.setVal(bkg_norm)

sig_width = r.RooRealVar("sig_width","sig_width", 0.1, 10)
sig_mean = r.RooRealVar("sig_mean","sig_mean", 120, 130)
sig_pdf = r.RooGaussian("sig_pdf","sig_pdf", M, sig_mean, sig_width)
sig_mean.setVal(sig_fnc.GetParameter(1))
sig_width.setVal(sig_fnc.GetParameter(2))
sig_norm = sig_fnc.Integral(100,180)
sig_pdf_norm = r.RooRealVar("sig_pdf_norm", "sig_pdf_norm", 0, 2*sig_norm)
sig_pdf_norm.setVal(sig_norm)
sig_pdf_norm.setConstant(True)

getattr(ws, 'import')(data)
getattr(ws, 'import')(data_binned)
getattr(ws, 'import')(bkg_exp)
getattr(ws, 'import')(bkg_exp_norm)
getattr(ws, 'import')(sig_pdf)
getattr(ws, 'import')(sig_pdf_norm)

ws.writeToFile("test_workspace.root")





