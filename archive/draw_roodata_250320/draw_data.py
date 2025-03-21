import ROOT
from ROOT import TFile, TCanvas, RooWorkspace, TH2D, RooFit, RooRealVar

# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# don't pop up a canvas and also remove stdout -> faster
ROOT.gROOT.SetBatch(True)

# ===== define functions ===== #
def draw_1d(var, title='', save_name='', label1='', label2='J/#psi#rightarrow#mu^{+}#mu^{-}'):
    canvas = TCanvas("canvas", "2D RooDataset Plot", 800, 600)

    hist = ds.createHistogram("hist", var)
    hist.Draw()
    hist.SetTitle(title)

    # # draw labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(0.89, 0.85, label1)
    latex.DrawLatexNDC(0.89, 0.80, label2)
    

    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.Update()
    canvas.Draw()
    if save_name != '':
        canvas.SaveAs('figs/' + save_name)

def draw_2d(hist2D, title='', save_name='', label=''):
    canvas = ROOT.TCanvas("canvas", "", 1200, 900)

    hist2D.Draw('colz')
    hist2D.SetTitle(title)

    # [hist, "Title", "saving name", "Kinematics label"

    texts = []  # python need it to store the objects

    # for i in range(1, hist2D.GetNbinsX() + 1):
    #     for j in range(1, hist2D.GetNbinsY() + 1):
    #         binContent = hist2D.GetBinContent(i, j)
    #         if binContent > 0:  # when content > 0
    #             x = hist2D.GetXaxis().GetBinCenter(i)
    #             y = hist2D.GetYaxis().GetBinCenter(j)
                
    #             text = ROOT.TText(x, y, f"{binContent:.2f}")
    #             text.SetTextSize(0.02)
    #             text.SetTextColor(ROOT.kBlack)
    #             text.SetTextAlign(22)  # draw value in the center of bin
    #             text.Draw()
    #             texts.append(text)  # python need it! Or text object will be overwritten by call.

    # draw labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(0.89, 0.80, 'J/#psi#rightarrow#mu^{+}#mu^{-}')
    latex.DrawLatexNDC(0.89, 0.85, label)

    canvas.SetRightMargin(0.15)
    canvas.Update()
    canvas.Draw()
    if save_name != '':
        canvas.SaveAs('figs/' + save_name)

# get root file and dataset
in_file = TFile('../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root')
ds = in_file.Get('dataset')
# ds.Print('v')


# ===== define local variables ===== #
mass = RooRealVar("mass", "mass (#mu^{+}#mu^{-})", 2.6, 3.5, "GeV/c^{2}")
pt = RooRealVar("pt", "p_{T} (#mu^{+}#mu^{-})", 0, 50, "GeV/c")
y = RooRealVar("y", "rapidity (#mu^{+}#mu^{-})", -3, 3, "")
pt1 = RooRealVar("pt1", "pt of #mu^{+}", 0, 50, "GeV/c")
eta1 = RooRealVar("eta1", "eta of #mu^{+}", -3, 3, "")
pt2 = RooRealVar("pt2", "pt of #mu^{-}", 0, 50, "GeV/c")
eta2 = RooRealVar("eta2", "eta of #mu^{-}", -3, 3, "")
cBin = RooRealVar("cBin", "Centrality bin", 0, 200, "")
# ep2 = RooRealVar("ep2", "2nd order event plane", -100, 100, "")
weight = RooRealVar("weight", "corr weight", 0, 10000, "")
recoQQsign = RooRealVar("recoQQsign", "qq sign", -1, 3, "")
ctau3D = RooRealVar("ctau3D", "c_{#tau}", -5.0, 5.0, "mm")
ctau3DErr = RooRealVar("ctau3DErr", "#sigma_{c#tau}", -0.1, 0.7, "mm")
ctau3DRes = RooRealVar("ctau3DRes", "c_{#tau} Res", -10, 10.0, "")

cos_theta = RooRealVar("cos_theta", "cos#theta", -1.0, 1.0, "")
cos_theta1 = RooRealVar("cos_theta1", "cos#theta (#mu^{+})", -1.0, 1.0, "")
phi = RooRealVar("phi", "", -5.0, 5.0, "#phi")
phi1 = RooRealVar("phi1", "", -5.0, 5.0, "#phi (#mu^{+})")

cos_cs = RooRealVar("cos_cs", "cos#theta_{CS}", -1.0, 1.0, "")
phi_cs = RooRealVar("phi_cs", "#phi_{CS}", -5.0, 5.0, "")
cos_hx = RooRealVar("cos_hx", "cos#theta_{HX}", -1.0, 1.0, "")
phi_hx = RooRealVar("phi_hx", "#phi_{HX}", -5.0, 5.0, "")
cos_ep = RooRealVar("cos_ep", "cos#theta_{EP}", -1.0, 1.0, "")
phi_ep = RooRealVar("phi_ep", "#phi_{EP}", -5.0, 5.0, "")
# phi_cs = RooRealVar("phi_cs", "", -5.0, 5.0, "")

# ===== List of variables ===== #
#     1)        mass = 3.0916  L(1 - 6) // [GeV/c^{2}] "mass variable"
#     2)          pt = 8.305  L(0 - 100) // [GeV/c] "pt variable"
#     3)           y = -1.77866  L(-5 - 5)  "rapidity of the dimuon pair"
#     4)         pt1 = 6.19354  L(0 - 500) // [GeV/c] "pt of muon+"
#     5)         pt2 = 2.41948  L(0 - 500) // [GeV/c] "pt of muon+"
#     6)        eta1 = -1.62806  L(-4 - 4)  "eta of muon+"
#     7)        eta2 = -2.15754  L(-4 - 4)  "eta of muon+"
#     8)      weight = 2.39843  L(0 - 10000)  "corr weight"
#     9)        cBin = 71  L(-100 - 500)  "Centrality bin"
#    10)  recoQQsign = 0  L(-1 - 3)  "qq sign"
#    11)     NumDimu = 1  L(0 - 100)  "number of dimuon"
#    12)      ctau3D = 0.0493783  L(-1000 - 1000) // [mm] "c_{#tau}"
#    13)   ctau3DErr = 0.0304145  L(-1000 - 1000) // [mm] "#sigma_{c#tau}"
#    14)   ctau3DRes = 1.62351  L(-1000 - 1000)  "c_{#tau} res"
#    15)   cos_theta = -0.95081  L(-1 - 1)  ""
#    16)  cos_theta1 = -0.925785  L(-1 - 1)  ""
#    17)         phi = -2.59855  L(-5 - 5)  ""
#    18)        phi1 = -2.43364  L(-5 - 5)  ""
#    19)      cos_cs = 0.607746  L(-1 - 1)  ""
#    20)      phi_cs = -2.16214  L(-5 - 5)  ""
#    21)      cos_hx = 0.216052  L(-1 - 1)  ""
#    22)      phi_hx = -0.741275  L(-5 - 5)  ""
#    23)      cos_ep = -0.813494  L(-1 - 1)  ""
#    24)      phi_ep = -2.29342  L(-5 - 5)  ""


# ===== draw 1d plots ===== #
draw_1d(mass, '', 'mass.png', '')
draw_1d(pt, '', 'pt.png', '')
draw_1d(y, '', 'y.png', '')
draw_1d(pt1, '', 'pt1.png', '', '#mu^{+}')
draw_1d(pt2, '', 'pt2.png', '', '#mu^{-}')
draw_1d(cBin, '', 'cBin.png', '')
draw_1d(weight, '', 'weight.png', '')
draw_1d(ctau3D, '', 'ctau3D.png', '')
draw_1d(ctau3DErr, '', 'ctau3DErr.png', '')
draw_1d(ctau3DRes, '', 'ctau3DRes.png', '')
draw_1d(cos_theta, '', 'cos_theta.png', '')
draw_1d(cos_theta1, '', 'cos_theta1.png', '', '#mu^{+}')
draw_1d(phi1, '', 'phi1.png', '', '#mu^{+}')
draw_1d(cos_cs, '', 'cos_cs.png', '', '')
draw_1d(phi_cs, '', 'phi_cs.png', '', '')
draw_1d(cos_hx, '', 'cos_hx.png', '', '')
draw_1d(phi_hx, '', 'phi_hx.png', '', '')
draw_1d(cos_ep, '', 'cos_ep.png', '', '')
draw_1d(phi_ep, '', 'phi_ep.png', '', '')


# ===== prepare 2d hists ===== #
h2_mass_pt = ds.createHistogram("h2_mass_pt", mass, YVar=dict(var=pt))
h2_mass_ctau3D = ds.createHistogram("h2_mass_ctau3D", mass, YVar=dict(var=ctau3D))
h2_ctau3D_ctauErr = ds.createHistogram("h2_ctau3D_ctauErr", ctau3D, YVar=dict(var=ctau3DErr))
h2_phi_cosTheta = ds.createHistogram("h2_phi_cosTheta", phi, YVar=dict(var=cos_theta))
h2_cos_ep_phi_ep = ds.createHistogram("h2_cos_ep_phi_ep", cos_ep, YVar=dict(var=phi_ep))
h2_cos_cs_phi_cs = ds.createHistogram("h2_cos_cs_phi_cs", cos_cs, YVar=dict(var=phi_cs))
h2_cos_hx_phi_hx = ds.createHistogram("h2_cos_hx_phi_hx", cos_hx, YVar=dict(var=phi_hx))
# h2_phi1_cosTheta1_CS = ds.createHistogram("h2_phi1_cosTheta1_CS", phi_cs, YVar=dict(var=cos_cs))
# h2_phi1_cosTheta1_HX = ds.createHistogram("h2_phi1_cosTheta1_HX", phi_hx, YVar=dict(var=cos_hx))

# ===== draw 2d plots ===== #
draw_2d(h2_mass_pt, '', '2d_mass_pt.png', label = '')
draw_2d(h2_mass_ctau3D, '', '2d_mass_ctau3D.png', label = '')
draw_2d(h2_ctau3D_ctauErr, '', '2d_ctau3D_ctauErr.png', label = '')
draw_2d(h2_phi_cosTheta, '', '2d_phi_cosTheta.png', label = '')
draw_2d(h2_cos_ep_phi_ep, '', '2d_cosEP_phiEP.png', label = '')
draw_2d(h2_cos_cs_phi_cs, '', '2d_cosCS_phiCS.png', label = '')
draw_2d(h2_cos_hx_phi_hx, '', '2d_cosHX_phiHX.png', label = '')

# draw_2d(h_fwd_cent0_20, 'Acc x Eff', 'fwd_cent0_20.png', 
#     label = '1.6 < |y| < 2.4, 0 < cent.(%) < 10')

# draw_2d(h_mid_cent20_60, 'Acc x Eff', 'mid_cent20_60.png', 
#     label = '0 < |y| < 1.6, 10 < cent.(%) < 30')
# draw_2d(h_fwd_cent20_60, 'Acc x Eff', 'fwd_cent20_60.png', 
#     label = '1.6 < |y| < 2.4, 10 < cent.(%) < 30')

# draw_2d(h_mid_cent60_180, 'Acc x Eff', 'mid_cent60_180.png', 
#     label = '0 < |y| < 1.6, 30 < cent.(%) < 90')
# draw_2d(h_fwd_cent60_180, 'Acc x Eff', 'fwd_cent60_180.png', 
#     label = '1.6 < |y| < 2.4, 30 < cent.(%) < 90')




# # ===== 2D graphs ===== #
# h2_mass_pt.Draw("LEGO2")
# c.Draw()
# c.SaveAs(f"figs/{sample}_RooDataset_2D_mass_pt.png")

# h2_mass_ctau3D.Draw("LEGO2")
# c.Draw()
# c.SaveAs(f"figs/{sample}_RooDataset_2D_mass_ctau3D.png")

# h2_phi_cosTheta.Draw("LEGO2")
# c.Draw()
# c.SaveAs(f"figs/{sample}_RooDataset_2D_phi_cosTheta.png")

# h2_phi1_cosTheta1_CS.Draw("LEGO2")
# c.Draw()
# c.SaveAs(f"figs/{sample}_RooDataset_2D_phi1_cosTheta1_CS.png")

# h2_phi1_cosTheta1_HX.Draw("LEGO2")
# c.Draw()
# c.SaveAs(f"figs/{sample}_RooDataset_2D_phi1_cosTheta1_HX.png")