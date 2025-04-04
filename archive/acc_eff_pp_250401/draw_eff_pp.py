import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D, TPad

def draw_1d_hist(hist1, legend1='', legend2='', x_title='', y_title='', save_name=''):
    # numerator
    c_comparison = TCanvas('c_1d', '', 600, 600)
    c_comparison.cd()
    ROOT.gPad.SetLogy()

    # legend
    legend = TLegend(0.6, 0.7, 0.9, 0.9)
    legend.SetBorderSize(0) # no boarder
    legend.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    legend.SetTextSize(0.03) # text size
    legend.AddEntry(hist1, legend1, 'ep')

    # cosmetics
    hist1.SetMarkerColor(ROOT.kBlue)
    hist1.SetLineColor(ROOT.kBlue)
    hist1.SetMarkerStyle(21)

    hist1.Draw('E1') # E1: error bar
    legend.Draw()

    hist1.GetXaxis().SetTitle(x_title)
    hist1.GetYaxis().SetTitle(y_title)

    c_comparison.Update()
    c_comparison.Draw()
    if save_name != '':
        c_comparison.SaveAs('figs/' + save_name + '.png')


def draw_1d_hist_maps(h_num, h_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    # numerator
    c_hist = TCanvas('c_1d_hist_maps', '', 1800, 600)
    c_hist.Divide(3)

    # numerator
    c_hist.cd(1)
    h_num.SetTitle('')
    h_num.GetXaxis().SetTitle(x_title)
    h_num.GetYaxis().SetTitle('')
    # h_num.GetYaxis().SetLabelSize(0.05)
    # h_num.GetYaxis().SetTitleOffset(1.5)
    h_num.Draw()

    leg1 = TLegend(0.6, 0.8, 0.9, 0.9)
    leg1.SetBorderSize(0) # no boarder
    leg1.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    leg1.SetTextSize(0.03) # text size
    leg1.AddEntry(h_num, 'Numerator', 'ep')
    leg1.Draw()


    # denominator
    c_hist.cd(2)
    h_den.SetTitle('')
    h_den.GetXaxis().SetTitle(x_title)
    h_den.GetYaxis().SetTitle('')
    # h_den.GetYaxis().SetLabelSize(0.05)
    # h_den.GetYaxis().SetTitleOffset(1.5)
    h_den.Draw()

    leg2 = TLegend(0.6, 0.8, 0.9, 0.9)
    leg2.SetBorderSize(0) # no boarder
    leg2.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    leg2.SetTextSize(0.03) # text size
    leg2.AddEntry(h_den, 'Denominator', 'ep')
    leg2.Draw()


    # correction
    c_hist.cd(3)
    h_cor = h_num.Clone("h_cor")
    # h_cor.Sumw2() # Don't need
    h_cor.Divide(h_den) # num / den. Divdie: divide by (parameter)
    
    h_cor.SetTitle('')
    h_cor.GetXaxis().SetTitle(x_title)
    h_cor.GetYaxis().SetTitle('')
    # h_cor.GetYaxis().SetLabelSize(0.05)
    # h_cor.GetYaxis().SetTitleOffset(1.5)
    h_cor.Draw()

    leg3 = TLegend(0.6, 0.8, 0.9, 0.9)
    leg3.SetBorderSize(0) # no boarder
    leg3.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    leg3.SetTextSize(0.03) # text size
    leg3.AddEntry(h_den, 'Efficiency', 'ep')
    leg3.Draw()

    c_hist.Update()
    c_hist.Draw()
    if save_name != '':
        c_hist.SaveAs('figs/' + save_name + '.png')

def draw_2d_hist_maps(h_num, h_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    # numerator
    c_hist = TCanvas('c_2d_hist_maps', '', 2500, 600)
    c_hist.Divide(3)

    # numerator
    c_hist.cd(1)
    h_num.SetTitle('')
    h_num.GetXaxis().SetTitle(x_title)
    h_num.GetYaxis().SetTitle(y_title)
    h_num.Draw('colz')

    lt1 = ROOT.TLatex()
    lt1.SetTextSize(0.03)
    lt1.SetTextAlign(31)
    lt1.DrawLatexNDC(0.89, 0.80, 'Numerator')
    lt1.DrawLatexNDC(0.89, 0.85, legend2)

    # denominator
    c_hist.cd(2)
    h_den.SetTitle('')
    h_den.GetXaxis().SetTitle(x_title)
    h_den.GetYaxis().SetTitle(y_title)
    h_den.Draw('colz')

    lt2 = ROOT.TLatex()
    lt2.SetTextSize(0.03)
    lt2.SetTextAlign(31)
    lt2.DrawLatexNDC(0.89, 0.80, 'Denominator')
    lt2.DrawLatexNDC(0.89, 0.85, legend2)

    # correction
    c_hist.cd(3)
    h_cor = h_num.Clone("h_cor")
    h_cor.Reset() # set 0

    for i in range(1, h_num.GetNbinsX() + 1):
        for j in range(1, h_num.GetNbinsY() + 1):
            content1 = h_num.GetBinContent(i, j)
            content2 = h_den.GetBinContent(i, j)
            
            if content2 != 0:
                eff = content1 / content2
                h_cor.SetBinContent(i, j, eff)
            else:
                h_cor.SetBinContent(i, j, 0)

    
    # h_cor.SetTitle(title)
    h_cor.GetXaxis().SetTitle(x_title)
    h_cor.GetYaxis().SetTitle(y_title)
    h_cor.SetMinimum(0.000001) # don't draw zero bin
    h_cor.Draw('colz')

    # draw labels
    lt3 = ROOT.TLatex()
    lt3.SetTextSize(0.03)
    lt3.SetTextAlign(31)
    lt3.DrawLatexNDC(0.89, 0.80, 'Efficiency')
    lt3.DrawLatexNDC(0.89, 0.85, legend2)


    c_hist.Update()
    c_hist.Draw()
    if save_name != '':
        c_hist.SaveAs('figs/' + save_name + '.png')

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


# ===== set macro configure ===== #
# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# batch mode - don't show the plot
ROOT.gROOT.SetBatch(True)


# ===== read inputs ===== #
# bring input root file
in_file_pp = TFile('roots/eff_pp_Jpsi_PtW1_tnp1_all_event.root')

# ===== get hists ===== #
# 3D map
h3_fwd_lab_num_pp = in_file_pp.Get('fwd_lab_num')
h3_fwd_lab_den_pp = in_file_pp.Get('fwd_lab_den')
h3_fwd_hx_num_pp = in_file_pp.Get('fwd_hx_num')
h3_fwd_hx_den_pp = in_file_pp.Get('fwd_hx_den')
h3_fwd_cs_num_pp = in_file_pp.Get('fwd_cs_num')
h3_fwd_cs_den_pp = in_file_pp.Get('fwd_cs_den')

h3_mid_lab_num_pp = in_file_pp.Get('mid_lab_num')
h3_mid_lab_den_pp = in_file_pp.Get('mid_lab_den')
h3_mid_hx_num_pp = in_file_pp.Get('mid_hx_num')
h3_mid_hx_den_pp = in_file_pp.Get('mid_hx_den')
h3_mid_cs_num_pp = in_file_pp.Get('mid_cs_num')
h3_mid_cs_den_pp = in_file_pp.Get('mid_cs_den')

# 1d map - y
h_y_lab_num_pp = in_file_pp.Get('y_lab_num')
h_y_lab_den_pp = in_file_pp.Get('y_lab_den')

# 1d map - cos
h_fwd_cos_lab_num_pp = h3_fwd_lab_num_pp.ProjectionX()
h_fwd_cos_lab_den_pp = h3_fwd_lab_den_pp.ProjectionX()
h_fwd_cos_hx_num_pp = h3_fwd_hx_num_pp.ProjectionX()
h_fwd_cos_hx_den_pp = h3_fwd_hx_den_pp.ProjectionX()
h_fwd_cos_cs_num_pp = h3_fwd_cs_num_pp.ProjectionX()
h_fwd_cos_cs_den_pp = h3_fwd_cs_den_pp.ProjectionX()

h_mid_cos_lab_num_pp = h3_mid_lab_num_pp.ProjectionX()
h_mid_cos_lab_den_pp = h3_mid_lab_den_pp.ProjectionX()
h_mid_cos_hx_num_pp = h3_mid_hx_num_pp.ProjectionX()
h_mid_cos_hx_den_pp = h3_mid_hx_den_pp.ProjectionX()
h_mid_cos_cs_num_pp = h3_mid_cs_num_pp.ProjectionX()
h_mid_cos_cs_den_pp = h3_mid_cs_den_pp.ProjectionX()

# 1d map - phi
h_fwd_phi_lab_num_pp = h3_fwd_lab_num_pp.ProjectionY()
h_fwd_phi_lab_den_pp = h3_fwd_lab_den_pp.ProjectionY()
h_fwd_phi_hx_num_pp = h3_fwd_hx_num_pp.ProjectionY()
h_fwd_phi_hx_den_pp = h3_fwd_hx_den_pp.ProjectionY()
h_fwd_phi_cs_num_pp = h3_fwd_cs_num_pp.ProjectionY()
h_fwd_phi_cs_den_pp = h3_fwd_cs_den_pp.ProjectionY()

h_mid_phi_lab_num_pp = h3_mid_lab_num_pp.ProjectionY()
h_mid_phi_lab_den_pp = h3_mid_lab_den_pp.ProjectionY()
h_mid_phi_hx_num_pp = h3_mid_hx_num_pp.ProjectionY()
h_mid_phi_hx_den_pp = h3_mid_hx_den_pp.ProjectionY()
h_mid_phi_cs_num_pp = h3_mid_cs_num_pp.ProjectionY()
h_mid_phi_cs_den_pp = h3_mid_cs_den_pp.ProjectionY()

# 1d map - pT
h_fwd_pT_lab_num_pp = h3_fwd_lab_num_pp.ProjectionZ()
h_fwd_pT_lab_den_pp = h3_fwd_lab_den_pp.ProjectionZ()
h_fwd_pT_hx_num_pp = h3_fwd_hx_num_pp.ProjectionZ()
h_fwd_pT_hx_den_pp = h3_fwd_hx_den_pp.ProjectionZ()
h_fwd_pT_cs_num_pp = h3_fwd_cs_num_pp.ProjectionZ()
h_fwd_pT_cs_den_pp = h3_fwd_cs_den_pp.ProjectionZ()

h_mid_pT_lab_num_pp = h3_mid_lab_num_pp.ProjectionZ()
h_mid_pT_lab_den_pp = h3_mid_lab_den_pp.ProjectionZ()
h_mid_pT_hx_num_pp = h3_mid_hx_num_pp.ProjectionZ()
h_mid_pT_hx_den_pp = h3_mid_hx_den_pp.ProjectionZ()
h_mid_pT_cs_num_pp = h3_mid_cs_num_pp.ProjectionZ()
h_mid_pT_cs_den_pp = h3_mid_cs_den_pp.ProjectionZ()

# 2d map - cos vs pT
h2_fwd_cos_pT_lab_num_pp = h3_fwd_lab_num_pp.Project3D("zx")
h2_fwd_cos_pT_lab_den_pp = h3_fwd_lab_den_pp.Project3D("zx")
h2_fwd_cos_pT_hx_num_pp = h3_fwd_hx_num_pp.Project3D("zx")
h2_fwd_cos_pT_hx_den_pp = h3_fwd_hx_den_pp.Project3D("zx")
h2_fwd_cos_pT_cs_num_pp = h3_fwd_cs_num_pp.Project3D("zx")
h2_fwd_cos_pT_cs_den_pp = h3_fwd_cs_den_pp.Project3D("zx")

h2_mid_cos_pT_lab_num_pp = h3_mid_lab_num_pp.Project3D("zx")
h2_mid_cos_pT_lab_den_pp = h3_mid_lab_den_pp.Project3D("zx")
h2_mid_cos_pT_hx_num_pp = h3_mid_hx_num_pp.Project3D("zx")
h2_mid_cos_pT_hx_den_pp = h3_mid_hx_den_pp.Project3D("zx")
h2_mid_cos_pT_cs_num_pp = h3_mid_cs_num_pp.Project3D("zx")
h2_mid_cos_pT_cs_den_pp = h3_mid_cs_den_pp.Project3D("zx")


# ===== draw 1d plots ===== #
# y
draw_1d_hist_maps(h_num=h_y_lab_num_pp, h_den=h_y_lab_den_pp, legend1='', legend2='', x_title='|y|', y_title='', save_name='eff_y_lab_pp')

# cos
draw_1d_hist_maps(h_num=h_fwd_cos_lab_num_pp, h_den=h_fwd_cos_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_fwd/cos_lab_pp')
draw_1d_hist_maps(h_num=h_fwd_cos_hx_num_pp, h_den=h_fwd_cos_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_fwd/cos_hx_pp')
draw_1d_hist_maps(h_num=h_fwd_cos_cs_num_pp, h_den=h_fwd_cos_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_fwd/cos_cs_pp')

draw_1d_hist_maps(h_num=h_mid_cos_lab_num_pp, h_den=h_mid_cos_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_mid/cos_lab_pp')
draw_1d_hist_maps(h_num=h_mid_cos_hx_num_pp, h_den=h_mid_cos_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_mid/cos_hx_pp')
draw_1d_hist_maps(h_num=h_mid_cos_cs_num_pp, h_den=h_mid_cos_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='', save_name='eff_mid/cos_cs_pp')

# phi
draw_1d_hist_maps(h_num=h_fwd_phi_lab_num_pp, h_den=h_fwd_phi_lab_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_fwd/phi_lab_pp')
draw_1d_hist_maps(h_num=h_fwd_phi_hx_num_pp, h_den=h_fwd_phi_hx_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_fwd/phi_hx_pp')
draw_1d_hist_maps(h_num=h_fwd_phi_cs_num_pp, h_den=h_fwd_phi_cs_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_fwd/phi_cs_pp')

draw_1d_hist_maps(h_num=h_mid_phi_lab_num_pp, h_den=h_mid_phi_lab_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_mid/phi_lab_pp')
draw_1d_hist_maps(h_num=h_mid_phi_hx_num_pp, h_den=h_mid_phi_hx_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_mid/phi_hx_pp')
draw_1d_hist_maps(h_num=h_mid_phi_cs_num_pp, h_den=h_mid_phi_cs_den_pp, legend1='', legend2='', x_title='#phi', y_title='', save_name='eff_mid/phi_cs_pp')

# pt
draw_1d_hist_maps(h_num=h_fwd_pT_lab_num_pp, h_den=h_fwd_pT_lab_den_pp, legend1='', legend2='', x_title='p_{T}', y_title='', save_name='eff_fwd/pT_lab_pp')
draw_1d_hist_maps(h_num=h_mid_pT_lab_num_pp, h_den=h_mid_pT_lab_den_pp, legend1='', legend2='', x_title='p_{T}', y_title='', save_name='eff_mid/pT_lab_pp')


# ===== draw 2d plots ===== #
# cos vs pT
draw_2d_hist_maps(h_num=h2_fwd_cos_pT_lab_num_pp, h_den=h2_fwd_cos_pT_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_fwd/2d_cos_pT_lab_pp')
draw_2d_hist_maps(h_num=h2_fwd_cos_pT_hx_num_pp, h_den=h2_fwd_cos_pT_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_fwd/2d_cos_pT_hx_pp')
draw_2d_hist_maps(h_num=h2_fwd_cos_pT_cs_num_pp, h_den=h2_fwd_cos_pT_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_fwd/2d_cos_pT_cs_pp')

draw_2d_hist_maps(h_num=h2_mid_cos_pT_lab_num_pp, h_den=h2_mid_cos_pT_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_mid/2d_cos_pT_lab_pp')
draw_2d_hist_maps(h_num=h2_mid_cos_pT_hx_num_pp, h_den=h2_mid_cos_pT_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_mid/2d_cos_pT_hx_pp')
draw_2d_hist_maps(h_num=h2_mid_cos_pT_cs_num_pp, h_den=h2_mid_cos_pT_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='eff_mid/2d_cos_pT_cs_pp')


# ===== Fwd cos vs phi - accroding to pT bins ===== # - 6
# pT integrated
h2_fwd_cos_phi_lab_num_pp = h3_fwd_lab_num_pp.Project3D("yx")
h2_fwd_cos_phi_lab_den_pp = h3_fwd_lab_den_pp.Project3D("yx")
h2_fwd_cos_phi_hx_num_pp = h3_fwd_hx_num_pp.Project3D("yx")
h2_fwd_cos_phi_hx_den_pp = h3_fwd_hx_den_pp.Project3D("yx")
h2_fwd_cos_phi_cs_num_pp = h3_fwd_cs_num_pp.Project3D("yx")
h2_fwd_cos_phi_cs_den_pp = h3_fwd_cs_den_pp.Project3D("yx")

draw_2d_hist_maps(h_num=h2_fwd_cos_phi_lab_num_pp, h_den=h2_fwd_cos_phi_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_lab_pp')
draw_2d_hist_maps(h_num=h2_fwd_cos_phi_hx_num_pp, h_den=h2_fwd_cos_phi_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_hx_pp')
draw_2d_hist_maps(h_num=h2_fwd_cos_phi_cs_num_pp, h_den=h2_fwd_cos_phi_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_cs_pp')


for idx in (1,2,3,4,5,6):
    h3_fwd_lab_num_pp.GetZaxis().SetRange(idx, idx)
    h3_fwd_lab_den_pp.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_num_pp.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_den_pp.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_num_pp.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_den_pp.GetZaxis().SetRange(idx, idx)
    h2_fwd_cos_phi_lab_num_pp = h3_fwd_lab_num_pp.Project3D("yx")
    h2_fwd_cos_phi_lab_den_pp = h3_fwd_lab_den_pp.Project3D("yx")
    h2_fwd_cos_phi_hx_num_pp = h3_fwd_hx_num_pp.Project3D("yx")
    h2_fwd_cos_phi_hx_den_pp = h3_fwd_hx_den_pp.Project3D("yx")
    h2_fwd_cos_phi_cs_num_pp = h3_fwd_cs_num_pp.Project3D("yx")
    h2_fwd_cos_phi_cs_den_pp = h3_fwd_cs_den_pp.Project3D("yx")


    fwd_pt_bins = ('pT3_6', 'pT6_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')
    name_tag = fwd_pt_bins[idx-1]

    draw_2d_hist_maps(h_num=h2_fwd_cos_phi_lab_num_pp, h_den=h2_fwd_cos_phi_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_lab_pp_' + name_tag)
    draw_2d_hist_maps(h_num=h2_fwd_cos_phi_hx_num_pp, h_den=h2_fwd_cos_phi_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_hx_pp_' + name_tag)
    draw_2d_hist_maps(h_num=h2_fwd_cos_phi_cs_num_pp, h_den=h2_fwd_cos_phi_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_fwd/2d_cos_phi_cs_pp_' + name_tag)


# ===== Mid cos vs phi - accroding to pT bins ===== # - 5
# pT integrated
# pT integrated
h2_mid_cos_phi_lab_num_pp = h3_mid_lab_num_pp.Project3D("yx")
h2_mid_cos_phi_lab_den_pp = h3_mid_lab_den_pp.Project3D("yx")
h2_mid_cos_phi_hx_num_pp = h3_mid_hx_num_pp.Project3D("yx")
h2_mid_cos_phi_hx_den_pp = h3_mid_hx_den_pp.Project3D("yx")
h2_mid_cos_phi_cs_num_pp = h3_mid_cs_num_pp.Project3D("yx")
h2_mid_cos_phi_cs_den_pp = h3_mid_cs_den_pp.Project3D("yx")

draw_2d_hist_maps(h_num=h2_mid_cos_phi_lab_num_pp, h_den=h2_mid_cos_phi_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_lab_pp')
draw_2d_hist_maps(h_num=h2_mid_cos_phi_hx_num_pp, h_den=h2_mid_cos_phi_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_hx_pp')
draw_2d_hist_maps(h_num=h2_mid_cos_phi_cs_num_pp, h_den=h2_mid_cos_phi_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_cs_pp')

# pT integrated
for idx in (1,2,3,4,5):
    h3_mid_lab_num_pp.GetZaxis().SetRange(idx, idx)
    h3_mid_lab_den_pp.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_num_pp.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_den_pp.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_num_pp.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_den_pp.GetZaxis().SetRange(idx, idx)
    h2_mid_cos_phi_lab_num_pp = h3_mid_lab_num_pp.Project3D("yx")
    h2_mid_cos_phi_lab_den_pp = h3_mid_lab_den_pp.Project3D("yx")
    h2_mid_cos_phi_hx_num_pp = h3_mid_hx_num_pp.Project3D("yx")
    h2_mid_cos_phi_hx_den_pp = h3_mid_hx_den_pp.Project3D("yx")
    h2_mid_cos_phi_cs_num_pp = h3_mid_cs_num_pp.Project3D("yx")
    h2_mid_cos_phi_cs_den_pp = h3_mid_cs_den_pp.Project3D("yx")

    mid_pt_bins = ('pT6p5_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')
    name_tag = mid_pt_bins[idx-1]

    draw_2d_hist_maps(h_num=h2_mid_cos_phi_lab_num_pp, h_den=h2_mid_cos_phi_lab_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_lab_pp_' + name_tag)
    draw_2d_hist_maps(h_num=h2_mid_cos_phi_hx_num_pp, h_den=h2_mid_cos_phi_hx_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_hx_pp_' + name_tag)
    draw_2d_hist_maps(h_num=h2_mid_cos_phi_cs_num_pp, h_den=h2_mid_cos_phi_cs_den_pp, legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='eff_mid/2d_cos_phi_cs_pp_' + name_tag)