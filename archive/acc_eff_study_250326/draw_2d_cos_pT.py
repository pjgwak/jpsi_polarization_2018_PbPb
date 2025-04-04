import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D


def draw_comparison_plot(hist1, hist2, legend1='', legend2='', x_title='', y_title='', save_name=''):
    # numerator
    c_comparison = TCanvas('c_comparison', '', 600, 600)
    c_comparison.cd()
    ROOT.gPad.SetLogy()

    # set proper Y range
    maxHist1 = hist1.GetMaximum()
    minHist1 = hist1.GetMinimum()
    maxHist2 = hist2.GetMaximum()
    minHist2 = hist2.GetMinimum()

    maxY = ROOT.TMath.Max(maxHist1, maxHist2)
    minY = ROOT.TMath.Min(minHist1, minHist2)

    hist1.GetYaxis().SetRangeUser(1, 5 * maxY)

    # legend
    legend = TLegend(0.6, 0.7, 0.9, 0.9)
    legend.SetBorderSize(0) # no boarder
    legend.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    legend.SetTextSize(0.03) # text size
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
        c_comparison.SaveAs('figs/' + save_name + '.png')

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
        c_comparison.SaveAs('figs/' + save_name + '.png')

def draw_2d_plot(hist2D, title='', x_title='', y_title='', save_name='', label1='', label2=''):
    canvas = ROOT.TCanvas("canvas", "", 1200, 900)

    hist2D.Draw('colz')
    hist2D.SetTitle(title)
    hist2D.GetXaxis().SetTitle(x_title)
    hist2D.GetYaxis().SetTitle(y_title)

    # draw labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(0.89, 0.80, label1)
    latex.DrawLatexNDC(0.89, 0.85, label2)

    canvas.Update()
    canvas.Draw()
    if save_name != '':
        canvas.SaveAs('figs/' + save_name + '.png')


def draw_2d_correction_plot(h2_num, h2_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    canvas = ROOT.TCanvas("canvas", "", 1200, 900)

    # calculate correction
    h_eff = h2_num.Clone("h_eff")
    h_eff.Reset() # set 0

    for i in range(1, h2_num.GetNbinsX() + 1):
        for j in range(1, h2_num.GetNbinsY() + 1):
            content1 = h2_num.GetBinContent(i, j)
            content2 = h2_den.GetBinContent(i, j)
            
            if content2 != 0:
                eff = content1 / content2
                h_eff.SetBinContent(i, j, eff)
            else:
                h_eff.SetBinContent(i, j, 0)

    
    # h_eff.SetTitle(title)
    h_eff.GetXaxis().SetTitle(x_title)
    h_eff.GetYaxis().SetTitle(y_title)
    h_eff.SetMinimum(0.01) # don't draw zero bin
    h_eff.Draw('colz')

    # draw labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(0.89, 0.80, legend1)
    latex.DrawLatexNDC(0.89, 0.85, legend2)

    canvas.Update()
    canvas.Draw()
    if save_name != '':
        canvas.SaveAs('figs/' + save_name + '.png')

# ===== List of histograms ===== #
# // fwd
# fwd_ep_num->Sumw2();
# fwd_ep_den->Sumw2();
# fwd_hx_num->Sumw2();
# fwd_hx_den->Sumw2();
# fwd_cs_num->Sumw2();
# fwd_cs_den->Sumw2();

# // mid
# mid_ep_num->Sumw2();
# mid_ep_den->Sumw2();
# mid_hx_num->Sumw2();
# mid_hx_den->Sumw2();
# mid_cs_num->Sumw2();
# mid_cs_den->Sumw2();


# ===== set macro configure ===== #
# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# batch mode - don't show the plot
ROOT.gROOT.SetBatch(True)


# ===== read inputs -PbPb ===== #
# bring input root file
in_file_pb = TFile('roots/acc_PbPb_Jpsi_PtW1_tnp1_all_event.root')

# ===== get hists ===== #
# 3D map
h3_fwd_ep_num_pb = in_file_pb.Get('fwd_ep_num')
h3_fwd_ep_den_pb = in_file_pb.Get('fwd_ep_den')
h3_fwd_hx_num_pb = in_file_pb.Get('fwd_hx_num')
h3_fwd_hx_den_pb = in_file_pb.Get('fwd_hx_den')
h3_fwd_cs_num_pb = in_file_pb.Get('fwd_cs_num')
h3_fwd_cs_den_pb = in_file_pb.Get('fwd_cs_den')

h3_mid_ep_num_pb = in_file_pb.Get('mid_ep_num')
h3_mid_ep_den_pb = in_file_pb.Get('mid_ep_den')
h3_mid_hx_num_pb = in_file_pb.Get('mid_hx_num')
h3_mid_hx_den_pb = in_file_pb.Get('mid_hx_den')
h3_mid_cs_num_pb = in_file_pb.Get('mid_cs_num')
h3_mid_cs_den_pb = in_file_pb.Get('mid_cs_den')

# 1d map - projection
# cos, phi, pT
h_fwd_cos_ep_num_pb = h3_fwd_ep_num_pb.ProjectionX()
h_fwd_phi_ep_num_pb = h3_fwd_ep_num_pb.ProjectionY()
h_fwd_pT_ep_num_pb = h3_fwd_ep_num_pb.ProjectionX()
h_fwd_cos_ep_den_pb = h3_fwd_ep_den_pb.ProjectionX()
h_fwd_phi_ep_den_pb = h3_fwd_ep_den_pb.ProjectionY()
h_fwd_pT_ep_den_pb = h3_fwd_ep_den_pb.ProjectionX()

h_fwd_cos_hx_num_pb = h3_fwd_hx_num_pb.ProjectionX()
h_fwd_phi_hx_num_pb = h3_fwd_hx_num_pb.ProjectionY()
h_fwd_pT_hx_num_pb = h3_fwd_hx_num_pb.ProjectionX()
h_fwd_cos_hx_den_pb = h3_fwd_hx_den_pb.ProjectionX()
h_fwd_phi_hx_den_pb = h3_fwd_hx_den_pb.ProjectionY()
h_fwd_pT_hx_den_pb = h3_fwd_hx_den_pb.ProjectionX()

h_fwd_cos_cs_num_pb = h3_fwd_cs_num_pb.ProjectionX()
h_fwd_phi_cs_num_pb = h3_fwd_cs_num_pb.ProjectionY()
h_fwd_pT_cs_num_pb = h3_fwd_cs_num_pb.ProjectionX()
h_fwd_cos_cs_den_pb = h3_fwd_cs_den_pb.ProjectionX()
h_fwd_phi_cs_den_pb = h3_fwd_cs_den_pb.ProjectionY()
h_fwd_pT_cs_den_pb = h3_fwd_cs_den_pb.ProjectionX()

# 2d map
h2_fwd_cos_phi_ep_num_pb = h3_fwd_ep_num_pb.Project3D("zx")
h2_fwd_cos_phi_hx_num_pb = h3_fwd_hx_num_pb.Project3D("zx")
h2_fwd_cos_phi_cs_num_pb = h3_fwd_cs_num_pb.Project3D("zx")
h2_fwd_cos_phi_ep_den_pb = h3_fwd_ep_den_pb.Project3D("zx")
h2_fwd_cos_phi_hx_den_pb = h3_fwd_hx_den_pb.Project3D("zx")
h2_fwd_cos_phi_cs_den_pb = h3_fwd_cs_den_pb.Project3D("zx")



# ===== read inputs and hists - pp Acc ===== #
# bring input root file
in_file_pp = TFile('roots/acc_pp_Jpsi_PtW1_tnp1_all_event.root')

# ===== get hists ===== #
# 3D map
h3_fwd_hx_num_pp = in_file_pp.Get('fwd_hx_num')
h3_fwd_hx_den_pp = in_file_pp.Get('fwd_hx_den')
h3_fwd_cs_num_pp = in_file_pp.Get('fwd_cs_num')
h3_fwd_cs_den_pp = in_file_pp.Get('fwd_cs_den')

h3_mid_hx_num_pp = in_file_pp.Get('mid_hx_num')
h3_mid_hx_den_pp = in_file_pp.Get('mid_hx_den')
h3_mid_cs_num_pp = in_file_pp.Get('mid_cs_num')
h3_mid_cs_den_pp = in_file_pp.Get('mid_cs_den')

# 1d map - projection
# 2d map
h2_fwd_cos_pt_hx_num_pp = h3_fwd_hx_num_pp.Project3D("zx")
h2_fwd_cos_pt_cs_num_pp = h3_fwd_cs_num_pp.Project3D("zx")
h2_fwd_cos_pt_hx_den_pp = h3_fwd_hx_den_pp.Project3D("zx")
h2_fwd_cos_pt_cs_den_pp = h3_fwd_cs_den_pp.Project3D("zx")


# ===== draw 2d plots ===== #
# h2_fwd_cos_phi_ep_num_pb = h3_fwd_ep_num_pb.Project3D("xy")
# h2_fwd_cos_phi_hx_num_pb = h3_fwd_hx_num_pb.Project3D("xy")
# h2_fwd_cos_phi_cs_num_pb = h3_fwd_cs_num_pb.Project3D("xy")

# cos vs phi - hx


# vector<double> fwd_pt_bin = {0, 3, 6.5, 9, 12, 50};
# c60_180_fwd_ep_num->GetZaxis()->SetRange(2,2);
# c60_180_fwd_ep_num->Project3D("zx")->Draw("colz");

draw_2d_plot(hist2D=h2_fwd_cos_pt_hx_num_pp,
    title='', x_title='cos#theta_{HX}', y_title='p_{T}', save_name='2d_cos_hx_vs_pT_pp_num', label1='pp', label2='')
draw_2d_plot(hist2D=h2_fwd_cos_pt_hx_den_pp,
    title='', x_title='cos#theta_{HX}', y_title='p_{T}', save_name='2d_cos_hx_vs_pT_pp_den', label1='pp', label2='')
draw_2d_correction_plot(h2_num=h2_fwd_cos_pt_hx_num_pp, h2_den=h2_fwd_cos_pt_hx_den_pp,
    legend1='pp', legend2='',
    x_title='cos#theta_{HX}', y_title='p_{T}',
    save_name='2d_cos_hx_vs_pT_pp')