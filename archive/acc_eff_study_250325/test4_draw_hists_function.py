import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D

# ===== List of histograms ===== #
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


# ===== List of comparison plots ===== #
# 1D - pT on vs off -> cos_ep 1 vs 2
# cos_ep num
# cos_ep den
# cos_ep correction

# projection comparison 1 - 1D vs projected 1D
# cos_ep1 vs projected cos_ep1
# cos_ep num
# cos_ep den
# cos_ep correction

# cos_ep2 vs projected cos_ep2
# cos_ep num
# cos_ep den
# cos_ep correction

# pT num
# pT den
# pT correction

# 2D comparison - pT on vs off
# cos_ep num
# cos_ep den
# cos_ep correction


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


h_c60_180_fwd_cos_ep_proj_num1 = h2_c60_180_fwd_num1.ProjectionX('h_c60_180_fwd_cos_ep_proj_num1')
h_c60_180_fwd_cos_ep_proj_den1 = h2_c60_180_fwd_den1.ProjectionX('h_c60_180_fwd_cos_ep_proj_den1')
h_c60_180_fwd_cos_ep_proj_num2 = h2_c60_180_fwd_num2.ProjectionX('h_c60_180_fwd_cos_ep_proj_num2')
h_c60_180_fwd_cos_ep_proj_den2 = h2_c60_180_fwd_den2.ProjectionX('h_c60_180_fwd_cos_ep_proj_den2')
h_c60_180_fwd_pt_proj_num1 = h2_c60_180_fwd_num1.ProjectionY('h_c60_180_fwd_pt_proj_num1')
h_c60_180_fwd_pt_proj_den1 = h2_c60_180_fwd_den1.ProjectionY('h_c60_180_fwd_pt_proj_den1')
h_c60_180_fwd_pt_proj_num2 = h2_c60_180_fwd_num2.ProjectionY('h_c60_180_fwd_pt_proj_num2')
h_c60_180_fwd_pt_proj_den2 = h2_c60_180_fwd_den2.ProjectionY('h_c60_180_fwd_pt_proj_den2')

def draw_comparison_plot(hist1, hist2, legend1='', legend2='', x_title='', y_title='', save_name=''):
    # numerator
    c_comparison = TCanvas('c_comparison', '', 600, 600)
    c_comparison.cd()
    ROOT.gPad.SetLogy()

    # legend
    legend = TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(hist1, legend1, 'ep')
    legend.AddEntry(hist2, legend2, 'ep')

    # cosmetics
    hist1.SetMarkerColor(ROOT.kBlue)
    hist2.SetMarkerColor(ROOT.kRed)
    hist1.SetLineColor(ROOT.kBlue)
    hist2.SetLineColor(ROOT.kRed)
    hist1.SetMarkerStyle(21)
    hist2.SetMarkerStyle(22)

    hist1.Draw('E1') # E1: error bar
    hist2.Draw('E1 same')
    legend.Draw()

    hist1.GetXaxis().SetTitle(x_title)
    hist1.GetYaxis().SetTitle(y_title)

    c_comparison.Update()
    c_comparison.Draw()
    if save_name != '':
        c_comparison.SaveAs('figs/1d_pT_weight_comparison/' + save_name + '.png')

def draw_correction_plot(h1_num, h1_den, h2_num, h2_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    hist1 = h1_num.Clone("hist1")  # clone and assign proper name
    # hist1.Sumw2() # Don't need applied when cloning the num.
    hist1.Divide(h1_den) # num / den. Divdie: divide by (parameter)

    hist2 = h2_num.Clone("hist2")
    hist2.Divide(h2_den)

    # numerator
    c_comparison = TCanvas('c_comparison', '', 600, 600)
    c_comparison.cd()

    # legend
    legend = TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(hist1, legend1, 'ep')
    legend.AddEntry(hist2, legend2, 'ep')

    # cosmetics
    hist1.SetMarkerColor(ROOT.kBlue)
    hist2.SetMarkerColor(ROOT.kRed)
    hist1.SetLineColor(ROOT.kBlue)
    hist2.SetLineColor(ROOT.kRed)
    hist1.SetMarkerStyle(21)
    hist2.SetMarkerStyle(22)

    hist1.Draw('E1') # E1: error bar
    hist2.Draw('E1 same')
    legend.Draw()
    
    hist1.SetMinimum(0)
    hist1.SetMaximum(1)
    hist1.GetXaxis().SetTitle(x_title)
    hist1.GetYaxis().SetTitle(y_title)

    c_comparison.Update()
    c_comparison.Draw()
    if save_name != '':
        c_comparison.SaveAs('figs/1d_pT_weight_comparison/' + save_name + '.png')


# 1D - pT on vs off -> cos_ep 1 vs 2
# cos_ep num
draw_comparison_plot(h_c60_180_fwd_cos_ep_num1, h_c60_180_fwd_cos_ep_num2,
    legend1='w/ p_{T} weight', legend2='w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (numerator)',
    save_name='cos_ep_pT_on_off_num')

# cos_ep den
draw_comparison_plot(h_c60_180_fwd_cos_ep_den1, h_c60_180_fwd_cos_ep_den2,
    legend1='w/ p_{T} weight', legend2='w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (denominator)',
    save_name='cos_ep_pT_on_off_den')

# cos_ep correction
draw_correction_plot(
    h1_num=h_c60_180_fwd_cos_ep_num1, h1_den=h_c60_180_fwd_cos_ep_den1,
    h2_num=h_c60_180_fwd_cos_ep_num2, h2_den=h_c60_180_fwd_cos_ep_den2,
    legend1='w/ p_{T} weight', legend2='w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='Acc x Eff',
    save_name='cos_ep_pT_on_off_correction')


# projection comparison 1 - 1D vs projected 1D
# cos_ep1 vs projected cos_ep1
# cos_ep num
draw_comparison_plot(h_c60_180_fwd_cos_ep_num1, h_c60_180_fwd_cos_ep_proj_num1,
    legend1='1D hist w/ p_{T} weight', legend2='Projected 1D hist w/ p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (numerator)',
    save_name='cos_ep_1d_vs_proj_1d_pT_on_num')
# cos_ep den
draw_comparison_plot(h_c60_180_fwd_cos_ep_den1, h_c60_180_fwd_cos_ep_proj_den1,
    legend1='1D hist w/ p_{T} weight', legend2='Projected 1D hist w/ p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (denominator)',
    save_name='cos_ep_1d_vs_proj_1d_pT_on_den')
# cos_ep correction
draw_correction_plot(
    h1_num=h_c60_180_fwd_cos_ep_num1, h1_den=h_c60_180_fwd_cos_ep_den1,
    h2_num=h_c60_180_fwd_cos_ep_proj_num1, h2_den=h_c60_180_fwd_cos_ep_proj_den1,
    legend1='1D hist w/ p_{T} weight', legend2='Projected 1D hist w/ p_{T} weight',
    x_title='cos#theta_{EP}', y_title='Acc x Eff',
    save_name='cos_ep_1d_vs_proj_1d_pT_on_correction')


# cos_ep2 vs projected cos_ep2
# cos_ep num
draw_comparison_plot(h_c60_180_fwd_cos_ep_num2, h_c60_180_fwd_cos_ep_proj_num2,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (numerator)',
    save_name='cos_ep_1d_vs_proj_1d_pT_off_num')
# cos_ep den
draw_comparison_plot(h_c60_180_fwd_cos_ep_den2, h_c60_180_fwd_cos_ep_proj_den2,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='# of dimuon (denominator)',
    save_name='cos_ep_1d_vs_proj_1d_pT_off_den')
# cos_ep correction
draw_correction_plot(
    h1_num=h_c60_180_fwd_cos_ep_num2, h1_den=h_c60_180_fwd_cos_ep_den2,
    h2_num=h_c60_180_fwd_cos_ep_proj_num2, h2_den=h_c60_180_fwd_cos_ep_proj_den2,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='cos#theta_{EP}', y_title='Acc x Eff',
    save_name='cos_ep_1d_vs_proj_1d_pT_off_correction')


# pt num
draw_comparison_plot(h_c60_180_fwd_pt_num, h_c60_180_fwd_pt_proj_num1,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='p_{T}', y_title='# of dimuon (numerator)',
    save_name='pt_1d_vs_proj_1d_num')
# pt den
draw_comparison_plot(h_c60_180_fwd_pt_den, h_c60_180_fwd_pt_proj_den1,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='p_{T}', y_title='# of dimuon (denominator)',
    save_name='pt_1d_vs_proj_1d_den')
# pt correction
draw_correction_plot(
    h1_num=h_c60_180_fwd_pt_num, h1_den=h_c60_180_fwd_pt_den,
    h2_num=h_c60_180_fwd_pt_proj_num1, h2_den=h_c60_180_fwd_pt_proj_den1,
    legend1='1D hist w/o p_{T} weight', legend2='Projected 1D hist w/o p_{T} weight',
    x_title='p_{T}', y_title='Acc x Eff',
    save_name='pt_1d_vs_proj_1d_correction')


# 2D comparison - pT on vs off -> Not need.
# differences btw projected 1D hists can be used for comparison