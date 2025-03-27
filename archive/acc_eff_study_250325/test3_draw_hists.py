import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D

#  List of histograms
#   KEY: TH2D	mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4_1;1	Acc x Eff
#   KEY: TH2D	mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4_2;1	Acc x Eff
#   KEY: TH2D	c60_180_fwd_num1;1
#   KEY: TH2D	c60_180_fwd_den1;1
#   KEY: TH2D	c60_180_fwd_num2;1	no pT weight
#   KEY: TH2D	c60_180_fwd_den2;1	no pT weight
#   KEY: TH1D	c60_180_fwd_pt_num;1
#   KEY: TH1D	c60_180_fwd_cos_ep_num1;1
#   KEY: TH1D	c60_180_fwd_cos_ep_num2;1	no pT weight
#   KEY: TH1D	c60_180_fwd_pt_den;1
#   KEY: TH1D	c60_180_fwd_cos_ep_den1;1
#   KEY: TH1D	c60_180_fwd_cos_ep_den2;1	no pT weight
#   KEY: TH2D	hpt_tnp_trig_fwd;1	hpt_tnp_trig_fwd


# ===== set macro configure ===== #
# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# batch mode - don't show the plot
ROOT.gROOT.SetBatch(True)


# ===== read inputs ===== #
# bring input root file
in_file = TFile('roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_test2.root')
# in_file.Print('V')

# ===== get hists ===== #
# 2D - correction
h2_fwd_cent60_180 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4')

# 2D - numbers. 1 means pT weight on, 2 means pT weight off
h2_c60_180_fwd_num1 = in_file.Get('c60_180_fwd_num1')
h2_c60_180_fwd_num2 = in_file.Get('c60_180_fwd_num2')
h2_c60_180_fwd_den1 = in_file.Get('c60_180_fwd_den1')
h2_c60_180_fwd_den2 = in_file.Get('c60_180_fwd_den2')


# 1D - numbers. 1 means pT weight on, 2 means pT weight off
h_c60_180_fwd_pt_num = in_file.Get('c60_180_fwd_pt_num') # pt - always turn on pT weight for pT
h_c60_180_fwd_cos_ep_num1 = in_file.Get('c60_180_fwd_cos_ep_num1') #cos1
h_c60_180_fwd_cos_ep_num2 = in_file.Get('c60_180_fwd_cos_ep_num2') #cos2
h_c60_180_fwd_pt_den = in_file.Get('c60_180_fwd_pt_den') # pt
h_c60_180_fwd_cos_ep_den1 = in_file.Get('c60_180_fwd_cos_ep_den1') #cos1
h_c60_180_fwd_cos_ep_den2 = in_file.Get('c60_180_fwd_cos_ep_den2') #cos2


# =============== draw 1D w/ and w/o number of events btw 1D - num and den separately  =============== #
# ===== numerator ===== #
c_cos_ep_num = TCanvas('c_cos_ep_num', '', 600, 600) # overlab cos_ep1, 2
c_cos_ep_num.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_num1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_num2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_num1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_num2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_num1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_num2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_num1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_num2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_num1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_num2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_num1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_num1.GetYaxis().SetTitle("# of dimuon (numerator)")

c_cos_ep_num.Update()
c_cos_ep_num.Draw()
c_cos_ep_num.SaveAs('figs/1d_pT_weight_comparison/cos_ep_num.png')


# ===== denominator ===== #
c_cos_ep_den = TCanvas('c_cos_ep_den', '', 600, 600) # overlab cos_ep1, 2
c_cos_ep_den.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_den1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_den2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_den1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_den2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_den1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_den2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_den1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_den2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_den1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_den2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_den1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_den1.GetYaxis().SetTitle("# of dimuon (Denominator)")

c_cos_ep_den.Update()
c_cos_ep_den.Draw()
c_cos_ep_den.SaveAs('figs/1d_pT_weight_comparison/cos_ep_den.png')


# ===== correction ===== #
c_cos_ep_cor = TCanvas('c_cos_ep_cor', '', 600, 600) # overlab cos_ep1, 2
c_cos_ep_cor.cd()

# make correction histograms
h_c60_180_fwd_cos_ep_cor1 = h_c60_180_fwd_cos_ep_num1.Clone("h_c60_180_fwd_cos_ep_cor1")  # clone and assign proper name
# h_c60_180_fwd_cos_ep_cor1.Sumw2() # applied when cloning the num.
h_c60_180_fwd_cos_ep_cor1.Divide(h_c60_180_fwd_cos_ep_den1) # num / den. Divdie: divide by (parameter)

h_c60_180_fwd_cos_ep_cor2 = h_c60_180_fwd_cos_ep_num2.Clone("h_c60_180_fwd_cos_ep_cor2")
# h_c60_180_fwd_cos_ep_cor2.Sumw2()
h_c60_180_fwd_cos_ep_cor2.Divide(h_c60_180_fwd_cos_ep_den2)

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_cor1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_cor2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_cor1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_cor2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_cor1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_cor2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_cor1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_cor2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_cor1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_cor2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_cor1.SetMinimum(0)
h_c60_180_fwd_cos_ep_cor1.SetMaximum(1)
h_c60_180_fwd_cos_ep_cor1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_cor1.GetYaxis().SetTitle("Acc x Eff")

c_cos_ep_cor.Update()
c_cos_ep_cor.Draw()
c_cos_ep_cor.SaveAs('figs/1d_pT_weight_comparison/cos_ep_cor.png')


# =============== 1D vs projected 1D - pT on, off =============== #


# =============== projected 1D vs projected 1D - pT on, off =============== #
h_c60_180_fwd_cos_ep_proj_num1 = h2_c60_180_fwd_num1.ProjectionX('h_c60_180_fwd_cos_ep_proj_num1')
h_c60_180_fwd_cos_ep_proj_den1 = h2_c60_180_fwd_den1.ProjectionX('h_c60_180_fwd_cos_ep_proj_den1')
h_c60_180_fwd_cos_ep_proj_num2 = h2_c60_180_fwd_num2.ProjectionX('h_c60_180_fwd_cos_ep_proj_num2')
h_c60_180_fwd_cos_ep_proj_den2 = h2_c60_180_fwd_den2.ProjectionX('h_c60_180_fwd_cos_ep_proj_den2')
h_c60_180_fwd_pt_proj_num1 = h2_c60_180_fwd_num1.ProjectionY('h_c60_180_fwd_pt_proj_num1')
h_c60_180_fwd_pt_proj_den1 = h2_c60_180_fwd_den1.ProjectionY('h_c60_180_fwd_pt_proj_den1')
h_c60_180_fwd_pt_proj_num2 = h2_c60_180_fwd_num2.ProjectionY('h_c60_180_fwd_pt_proj_num2')
h_c60_180_fwd_pt_proj_den2 = h2_c60_180_fwd_den2.ProjectionY('h_c60_180_fwd_pt_proj_den2')

# ===== cosEP - numerator ===== #
c_cos_ep_proj_num = TCanvas('c_cos_ep_proj_num', '', 600, 600) # overlab cos_ep_proj1, 2
c_cos_ep_proj_num.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_num1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_num2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_proj_num1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_num2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_num1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_num2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_num1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_proj_num2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_proj_num1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_proj_num2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_proj_num1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_proj_num1.GetYaxis().SetTitle("# of dimuon (numerator)")

c_cos_ep_proj_num.Update()
c_cos_ep_proj_num.Draw()
c_cos_ep_proj_num.SaveAs('figs/1d_pT_weight_comparison/cos_ep_proj_num.png')


# ===== cosEP - denominator ===== #
c_cos_ep_proj_den = TCanvas('c_cos_ep_proj_den', '', 600, 600) # overlab cos_ep_proj1, 2
c_cos_ep_proj_den.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_den1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_den2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_proj_den1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_den2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_den1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_den2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_den1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_proj_den2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_proj_den1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_proj_den2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_proj_den1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_proj_den1.GetYaxis().SetTitle("# of dimuon (Denominator)")

c_cos_ep_proj_den.Update()
c_cos_ep_proj_den.Draw()
c_cos_ep_proj_den.SaveAs('figs/1d_pT_weight_comparison/cos_ep_proj_den.png')


# ===== cosEP - correction ===== #
c_cos_ep_proj_cor = TCanvas('c_cos_ep_proj_cor', '', 600, 600) # overlab cos_ep_proj1, 2
c_cos_ep_proj_cor.cd()

# make correction histograms
h_c60_180_fwd_cos_ep_proj_cor1 = h_c60_180_fwd_cos_ep_proj_num1.Clone("h_c60_180_fwd_cos_ep_proj_cor1")  # clone and assign proper name
# h_c60_180_fwd_cos_ep_proj_cor1.Sumw2() # applied when cloning the num.
h_c60_180_fwd_cos_ep_proj_cor1.Divide(h_c60_180_fwd_cos_ep_proj_den1) # num / den. Divdie: divide by (parameter)

h_c60_180_fwd_cos_ep_proj_cor2 = h_c60_180_fwd_cos_ep_proj_num2.Clone("h_c60_180_fwd_cos_ep_proj_cor2")
# h_c60_180_fwd_cos_ep_proj_cor2.Sumw2()
h_c60_180_fwd_cos_ep_proj_cor2.Divide(h_c60_180_fwd_cos_ep_proj_den2)

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_cor1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_cos_ep_proj_cor2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_cos_ep_proj_cor1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_cor2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_cor1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_cos_ep_proj_cor2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_cos_ep_proj_cor1.SetMarkerStyle(21)
h_c60_180_fwd_cos_ep_proj_cor2.SetMarkerStyle(22)

h_c60_180_fwd_cos_ep_proj_cor1.Draw('E1') # E1: error bar
h_c60_180_fwd_cos_ep_proj_cor2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_cos_ep_proj_cor1.SetMinimum(0)
h_c60_180_fwd_cos_ep_proj_cor1.SetMaximum(1)
h_c60_180_fwd_cos_ep_proj_cor1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_cos_ep_proj_cor1.GetYaxis().SetTitle("Acc x Eff")

c_cos_ep_proj_cor.Update()
c_cos_ep_proj_cor.Draw()
c_cos_ep_proj_cor.SaveAs('figs/1d_pT_weight_comparison/cos_ep_proj_cor.png')



# ===== pT - numerator ===== #
c_pt_proj_num = TCanvas('c_pt_proj_num', '', 600, 600) # overlab pt_proj1, 2
c_pt_proj_num.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_proj_num1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_pt_proj_num2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_proj_num1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_num2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_num1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_num2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_num1.SetMarkerStyle(21)
h_c60_180_fwd_pt_proj_num2.SetMarkerStyle(22)

h_c60_180_fwd_pt_proj_num1.Draw('E1') # E1: error bar
h_c60_180_fwd_pt_proj_num2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_pt_proj_num1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_proj_num1.GetYaxis().SetTitle("# of dimuon (numerator)")

c_pt_proj_num.Update()
c_pt_proj_num.Draw()
c_pt_proj_num.SaveAs('figs/1d_pT_weight_comparison/pt_proj_num.png')


# ===== pT - denominator ===== #
c_pt_proj_den = TCanvas('c_pt_proj_den', '', 600, 600) # overlab pt_proj1, 2
c_pt_proj_den.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_proj_den1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_pt_proj_den2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_proj_den1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_den2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_den1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_den2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_den1.SetMarkerStyle(21)
h_c60_180_fwd_pt_proj_den2.SetMarkerStyle(22)

h_c60_180_fwd_pt_proj_den1.Draw('E1') # E1: error bar
h_c60_180_fwd_pt_proj_den2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_pt_proj_den1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_proj_den1.GetYaxis().SetTitle("# of dimuon (Denominator)")

c_pt_proj_den.Update()
c_pt_proj_den.Draw()
c_pt_proj_den.SaveAs('figs/1d_pT_weight_comparison/pt_proj_den.png')


# ===== pT - correction ===== #
c_pt_proj_cor = TCanvas('c_pt_proj_cor', '', 600, 600) # overlab pt_proj1, 2
c_pt_proj_cor.cd()

# make correction histograms
h_c60_180_fwd_pt_proj_cor1 = h_c60_180_fwd_pt_proj_num1.Clone("h_c60_180_fwd_pt_proj_cor1")  # clone and assign proper name
# h_c60_180_fwd_pt_proj_cor1.Sumw2() # applied when cloning the num.
h_c60_180_fwd_pt_proj_cor1.Divide(h_c60_180_fwd_pt_proj_den1) # num / den. Divdie: divide by (parameter)

h_c60_180_fwd_pt_proj_cor2 = h_c60_180_fwd_pt_proj_num2.Clone("h_c60_180_fwd_pt_proj_cor2")
# h_c60_180_fwd_pt_proj_cor2.Sumw2()
h_c60_180_fwd_pt_proj_cor2.Divide(h_c60_180_fwd_pt_proj_den2)

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_proj_cor1, 'w/ p_{T} weight', 'ep')
legend.AddEntry(h_c60_180_fwd_pt_proj_cor2, 'w/o p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_proj_cor1.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_cor2.SetMarkerColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_cor1.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_proj_cor2.SetLineColor(ROOT.kRed)
h_c60_180_fwd_pt_proj_cor1.SetMarkerStyle(21)
h_c60_180_fwd_pt_proj_cor2.SetMarkerStyle(22)

h_c60_180_fwd_pt_proj_cor1.Draw('E1') # E1: error bar
h_c60_180_fwd_pt_proj_cor2.Draw('E1 same')
legend.Draw()

h_c60_180_fwd_pt_proj_cor1.SetMinimum(0)
h_c60_180_fwd_pt_proj_cor1.SetMaximum(1)
h_c60_180_fwd_pt_proj_cor1.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_proj_cor1.GetYaxis().SetTitle("Acc x Eff")

c_pt_proj_cor.Update()
c_pt_proj_cor.Draw()
c_pt_proj_cor.SaveAs('figs/1d_pT_weight_comparison/pt_proj_cor.png')

# =============== draw 1D vs projected 1D - pT weight off (num, den, correction)=============== #



# =============== draw 2D plots with values - digit .3f =============== #
# pt weight on, off


# =============== Start pT corrections =============== #
# ===== numerator ===== #
c_pt_num = TCanvas('c_pt_num', '', 600, 600) # overlab pt1, 2
c_pt_num.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_num, 'w/ p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_num.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_num.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_num.SetMarkerStyle(21)

h_c60_180_fwd_pt_num.Draw('E1') # E1: error bar
legend.Draw()

h_c60_180_fwd_pt_num.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_num.GetYaxis().SetTitle("# of dimuon (numerator)")

c_pt_num.Update()
c_pt_num.Draw()
c_pt_num.SaveAs('figs/1d_pT_weight_comparison/pt_num.png')


# ===== denominator ===== #
c_pt_den = TCanvas('c_pt_den', '', 600, 600) # overlab pt1, 2
c_pt_den.cd()
ROOT.gPad.SetLogy()

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_den, 'w/ p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_den.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_den.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_den.SetMarkerStyle(21)

h_c60_180_fwd_pt_den.Draw('E1') # E1: error bar
legend.Draw()

h_c60_180_fwd_pt_den.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_den.GetYaxis().SetTitle("# of dimuon (Denominator)")

c_pt_den.Update()
c_pt_den.Draw()
c_pt_den.SaveAs('figs/1d_pT_weight_comparison/pt_den.png')


# ===== correction ===== #
c_pt_cor = TCanvas('c_pt_cor', '', 600, 600) # overlab pt1, 2
c_pt_cor.cd()

# make correction histograms
h_c60_180_fwd_pt_cor = h_c60_180_fwd_pt_num.Clone("h_c60_180_fwd_pt_cor")  # clone and assign proper name
# h_c60_180_fwd_pt_cor.Sumw2() # applied when cloning the num.
h_c60_180_fwd_pt_cor.Divide(h_c60_180_fwd_pt_den) # num / den. Divdie: divide by (parameter)

# legend
legend = TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_c60_180_fwd_pt_cor, 'w/ p_{T} weight', 'ep')

# cosmetics
h_c60_180_fwd_pt_cor.SetMarkerColor(ROOT.kBlue)
h_c60_180_fwd_pt_cor.SetLineColor(ROOT.kBlue)
h_c60_180_fwd_pt_cor.SetMarkerStyle(21)

h_c60_180_fwd_pt_cor.Draw('E1') # E1: error bar
legend.Draw()

h_c60_180_fwd_pt_cor.SetMinimum(0)
h_c60_180_fwd_pt_cor.SetMaximum(1)
h_c60_180_fwd_pt_cor.GetXaxis().SetTitle("cos#theta_{EP}")
h_c60_180_fwd_pt_cor.GetYaxis().SetTitle("Acc x Eff")

c_pt_cor.Update()
c_pt_cor.Draw()
c_pt_cor.SaveAs('figs/1d_pT_weight_comparison/pt_cor.png')