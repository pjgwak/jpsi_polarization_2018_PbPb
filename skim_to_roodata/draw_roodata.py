import ROOT
from ROOT import TFile, TCanvas, RooWorkspace, TH2D, RooFit, RooRealVar

# is_draw_all False일 때 sample, 1D(2D), var1, var2는 prompt에서 전댈하는 게 더 편할지도?

is_draw_all = False
sample = 'mc_pr'
# sample = 'mc_pr'
# sample = 'mc_np'

in_file = TFile()
if sample == 'data':
    in_file = TFile('../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw0_Accw0_PtW0_TnP0_HFHFNom_250201.root')
elif sample == 'mc_pr':
    in_file = TFile('../files_roodata/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root')
elif sample == 'mc_np':
    in_file = TFile('../files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_180_Effw0_Accw0_PtW0_TnP0_HFHFNom_250201.root')

ds = in_file.Get('dataset')

mass = RooRealVar("mass", "mass variable", 2.6, 3.5, "GeV/c^{2}")
pt = RooRealVar("pt", "pt variable", 0, 50, "GeV/c")
y = RooRealVar("y", "rapidity of the dimuon pair", -3, 3, "")
pt1 = RooRealVar("pt1", "pt of muon+", 0, 50, "GeV/c")
eta1 = RooRealVar("eta1", "eta of muon+", -3, 3, "")
pt2 = RooRealVar("pt2", "pt of muon+", 0, 50, "GeV/c")
eta2 = RooRealVar("eta2", "eta of muon+", -3, 3, "")
cBin = RooRealVar("cBin", "Centrality bin", 0, 200, "")
# ep2 = RooRealVar("ep2", "2nd order event plane", -100, 100, "")
weight = RooRealVar("weight", "corr weight", 0, 10000, "")
recoQQsign = RooRealVar("recoQQsign", "qq sign", -1, 3, "")
ctau3D = RooRealVar("ctau3D", "c_{#tau}", -5.0, 5.0, "mm")
ctau3DErr = RooRealVar("ctau3DErr", "#sigma_{c#tau}", -5.0, 5.0, "mm")
ctau3DRes = RooRealVar("ctau3DRes", "c_{#tau}", -0.5, 2.0, "")
evtWeight = RooRealVar("evtWeightVar", "c_{#tau}", -1000, 1000, "")

cos_theta = RooRealVar("cos_theta", "", -1.0, 1.0, "")
cos_theta1 = RooRealVar("cos_theta1", "", -1.0, 1.0, "")
phi = RooRealVar("phi", "", -5.0, 5.0, "")
phi1 = RooRealVar("phi1", "", -5.0, 5.0, "")

cos_cs = RooRealVar("cos_cs", "", -1.0, 1.0, "")
phi_cs = RooRealVar("phi_cs", "", -5.0, 5.0, "")
cos_hx = RooRealVar("cos_hx", "", -1.0, 1.0, "")
phi_hx = RooRealVar("phi_hx", "", -5.0, 5.0, "")
cos_ep = RooRealVar("cos_ep", "", -1.0, 1.0, "")
phi_ep = RooRealVar("phi_ep", "", -5.0, 5.0, "")
# phi_cs = RooRealVar("phi_cs", "", -5.0, 5.0, "")

if is_draw_all:
    ROOT.gROOT.SetBatch(True)

c = TCanvas("c", "2D RooDataset Plot", 800, 600)

if not is_draw_all:
    is1D = True
    
    if is1D: # 1D
        my_hist = ds.createHistogram("my_hist", phi_ep)
        my_hist.Draw()
    else: # 2D
        my_list = [mass, ctau3DRes]
        my_hist = ds.createHistogram("my_hist", my_list[0], YVar=dict(var=my_list[1]))
        my_hist.SetTitle(f';{my_list[0].GetName()};{my_list[1].GetName()}')
        my_hist.Draw("LEGO2")

    c.Draw()
    # c.SaveAs(f"figs/{sample}_RooDataset_1D.png")


# Draw all interesting graphs
if is_draw_all:
    # ===== 1D graphs ===== #
    variable_list = [mass, pt, y, pt1, ctau3D, ctau3DErr, ctau3DRes,
    cos_theta, cos_theta1, phi, phi1, cos_cs,
    phi_cs, cos_hx, phi_hx]

    for var in variable_list:
        name = var.GetName()
        my_hist = ds.createHistogram('h', var)
        # my_hist = ds.createHistogram('h', var, RooFit::WeightVar("evtWeight"))
        my_hist.Draw()
        c.Draw()
        c.SaveAs(f"figs/{sample}_RooDataset_1D_{name}.png")


    # ===== 2D graphs ===== #
    h2_mass_pt = ds.createHistogram("hist", mass, YVar=dict(var=pt))
    h2_mass_ctau3D = ds.createHistogram("hist", ctau3D, YVar=dict(var=mass))
    h2_phi_cosTheta = ds.createHistogram("hist", phi, YVar=dict(var=cos_theta))
    h2_phi1_cosTheta1_CS = ds.createHistogram("hist", phi_cs, YVar=dict(var=cos_cs))
    h2_phi1_cosTheta1_HX = ds.createHistogram("hist", phi_hx, YVar=dict(var=cos_hx))

    h2_mass_pt.Draw("LEGO2")
    c.Draw()
    c.SaveAs(f"figs/{sample}_RooDataset_2D_mass_pt.png")

    h2_mass_ctau3D.Draw("LEGO2")
    c.Draw()
    c.SaveAs(f"figs/{sample}_RooDataset_2D_mass_ctau3D.png")

    h2_phi_cosTheta.Draw("LEGO2")
    c.Draw()
    c.SaveAs(f"figs/{sample}_RooDataset_2D_phi_cosTheta.png")

    h2_phi1_cosTheta1_CS.Draw("LEGO2")
    c.Draw()
    c.SaveAs(f"figs/{sample}_RooDataset_2D_phi1_cosTheta1_CS.png")

    h2_phi1_cosTheta1_HX.Draw("LEGO2")
    c.Draw()
    c.SaveAs(f"figs/{sample}_RooDataset_2D_phi1_cosTheta1_HX.png")