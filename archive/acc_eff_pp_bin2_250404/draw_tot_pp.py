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


def draw_1d_cor(h_acc_num, h_acc_den, h_eff_num, h_eff_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    c_hist = TCanvas('c_1d_hist', '', 600, 600)
    c_hist.cd()
    # ROOT.gPad.SetLogy()

    # Acceptance
    h_acc = h_acc_num.Clone("h_acc")
    h_acc.Divide(h_acc_den)

    # Efficiency
    h_eff = h_eff_num.Clone("h_eff")
    h_eff.Divide(h_eff_den)
    
    # Acc x Eff
    h_cor = h_acc_num.Clone("h_cor")
    h_cor.Reset()

    for bin in range(1, h_acc.GetNbinsX() + 1):
        bin_content1 = h_acc.GetBinContent(bin)
        bin_content2 = h_eff.GetBinContent(bin)
        prod = bin_content1 * bin_content2
        h_cor.SetBinContent(bin, prod)
    
    # h_cor.SetTitle('')
    # h_cor.GetXaxis().SetTitle(x_title)
    # h_cor.GetYaxis().SetTitle('')
    # h_cor.GetYaxis().SetLabelSize(0.05)
    # h_cor.GetYaxis().SetTitleOffset(1.5)
    
    h_cor.SetMarkerColor(ROOT.kBlue)
    h_cor.SetLineColor(ROOT.kBlue)
    h_cor.SetMarkerStyle(2)

    h_acc.SetMarkerColor(ROOT.kBlack)
    h_eff.SetMarkerColor(ROOT.kRed)
    h_acc.SetLineColor(ROOT.kBlack)
    h_eff.SetLineColor(ROOT.kRed)
    h_acc.SetMarkerStyle(21)
    h_eff.SetMarkerStyle(22)
    
    h_acc.SetMaximum(1.2)
    h_acc.SetMinimum(1e-5)

    h_acc.Draw('ep')
    h_eff.Draw('same ep')
    h_cor.Draw('same ep')

    leg = TLegend(0.6, 0.8, 0.9, 0.9)
    leg.SetBorderSize(0) # no boarder
    leg.SetFillColorAlpha(ROOT.kWhite, 0) # transperency
    leg.SetTextSize(0.03) # text size
    leg.AddEntry(h_acc, 'Acc', 'ep')
    leg.AddEntry(h_eff, 'Eff', 'ep')
    leg.AddEntry(h_cor, 'Acc x Eff', 'ep')
    leg.Draw()

    c_hist.Update()
    c_hist.Draw()
    if save_name != '':
        c_hist.SaveAs('figs/' + save_name + '.png')

def draw_2d_cor(h_acc_num, h_acc_den, h_eff_num, h_eff_den, legend1='', legend2='', x_title='', y_title='', save_name=''):
    c_hist = TCanvas('c_2d_hist', '', 800, 600)
    
    # Acc
    h_acc = h_acc_num.Clone("h_acc")
    h_acc.Reset() # set 0
    for i in range(1, h_acc_num.GetNbinsX() + 1):
        for j in range(1, h_acc_num.GetNbinsY() + 1):
            content1 = h_acc_num.GetBinContent(i, j)
            content2 = h_acc_den.GetBinContent(i, j)
            
            if content2 != 0:
                acc = content1 / content2
                h_acc.SetBinContent(i, j, acc)
            else:
                h_acc.SetBinContent(i, j, 0)

    # Eff
    h_eff = h_eff_num.Clone("h_eff")
    h_eff.Reset()
    for i in range(1, h_eff_num.GetNbinsX() + 1):
        for j in range(1, h_eff_num.GetNbinsY() + 1):
            content1 = h_eff_num.GetBinContent(i, j)
            content2 = h_eff_den.GetBinContent(i, j)
            
            if content2 != 0:
                eff = content1 / content2
                h_eff.SetBinContent(i, j, eff)
            else:
                h_eff.SetBinContent(i, j, 0)    

    # Acc x Eff
    h_cor = h_acc_num.Clone("h_cor")
    h_cor.Reset() # set 0
    for x_bin in range(1, h_acc.GetNbinsX() + 1):
        for y_bin in range(1, h_acc.GetNbinsY() + 1):
            bin_content_1 = h_acc.GetBinContent(x_bin, y_bin)
            bin_content_2 = h_eff.GetBinContent(x_bin, y_bin)
            prod = bin_content_1 * bin_content_2
            # comment-out only for plotting -> 안 쓰는 영역이 주로 잘려서 켜도 괜찮은데 정확한 그림은 아니니까
            # 두 장을 그려야 한다.
            # if prod > 1.2:
            #     prod = 0
            h_cor.SetBinContent(x_bin, y_bin, prod)
    
    # h_cor.SetTitle(title)
    h_cor.GetXaxis().SetTitle(x_title)
    h_cor.GetYaxis().SetTitle(y_title)
    h_cor.SetMinimum(0.001) # don't draw zero bin
    # h_cor.SetMaximum(1.2) # don't draw zero bin
    h_cor.Draw('colz')

    # draw labels
    lt3 = ROOT.TLatex()
    lt3.SetTextSize(0.03)
    lt3.SetTextAlign(31)
    lt3.DrawLatexNDC(0.89, 0.80, 'Acc x Eff')
    lt3.DrawLatexNDC(0.89, 0.85, legend2)



    texts = []  # python need it to store the objects

    for i in range(1, h_cor.GetNbinsX() + 1):
        for j in range(1, h_cor.GetNbinsY() + 1):
            binContent = h_cor.GetBinContent(i, j)
            if binContent > 0:  # when content > 0
                x = h_cor.GetXaxis().GetBinCenter(i)
                y = h_cor.GetYaxis().GetBinCenter(j)
                
                text = ROOT.TText(x, y, f"{binContent:.2f}")
                text.SetTextSize(0.02)
                text.SetTextColor(ROOT.kBlack)
                text.SetTextAlign(22)  # draw value in the center of bin
                text.Draw()
                texts.append(text)  # python need it! Or text object will be overwritten by call.

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
    canvas = ROOT.TCanvas("canvas", "", 600, 600)

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


# ===== read inputs - Acc ===== #
# bring input root file
in_file_acc = TFile('roots/acc_pp_Jpsi_PtW1_tnp1_all_event.root')

# ===== get hists ===== #
# 3D map
h3_fwd_lab_num_acc = in_file_acc.Get('fwd_lab_num')
h3_fwd_lab_den_acc = in_file_acc.Get('fwd_lab_den')
h3_fwd_hx_num_acc = in_file_acc.Get('fwd_hx_num')
h3_fwd_hx_den_acc = in_file_acc.Get('fwd_hx_den')
h3_fwd_cs_num_acc = in_file_acc.Get('fwd_cs_num')
h3_fwd_cs_den_acc = in_file_acc.Get('fwd_cs_den')

h3_mid_lab_num_acc = in_file_acc.Get('mid_lab_num')
h3_mid_lab_den_acc = in_file_acc.Get('mid_lab_den')
h3_mid_hx_num_acc = in_file_acc.Get('mid_hx_num')
h3_mid_hx_den_acc = in_file_acc.Get('mid_hx_den')
h3_mid_cs_num_acc = in_file_acc.Get('mid_cs_num')
h3_mid_cs_den_acc = in_file_acc.Get('mid_cs_den')

# 1d map - y
h_y_lab_num_acc = in_file_acc.Get('y_lab_num')
h_y_lab_den_acc = in_file_acc.Get('y_lab_den')

# 1d map - cos
h_fwd_cos_lab_num_acc = h3_fwd_lab_num_acc.ProjectionX()
h_fwd_cos_lab_den_acc = h3_fwd_lab_den_acc.ProjectionX()
h_fwd_cos_hx_num_acc = h3_fwd_hx_num_acc.ProjectionX()
h_fwd_cos_hx_den_acc = h3_fwd_hx_den_acc.ProjectionX()
h_fwd_cos_cs_num_acc = h3_fwd_cs_num_acc.ProjectionX()
h_fwd_cos_cs_den_acc = h3_fwd_cs_den_acc.ProjectionX()

h_mid_cos_lab_num_acc = h3_mid_lab_num_acc.ProjectionX()
h_mid_cos_lab_den_acc = h3_mid_lab_den_acc.ProjectionX()
h_mid_cos_hx_num_acc = h3_mid_hx_num_acc.ProjectionX()
h_mid_cos_hx_den_acc = h3_mid_hx_den_acc.ProjectionX()
h_mid_cos_cs_num_acc = h3_mid_cs_num_acc.ProjectionX()
h_mid_cos_cs_den_acc = h3_mid_cs_den_acc.ProjectionX()

# 1d map - phi
h_fwd_phi_lab_num_acc = h3_fwd_lab_num_acc.ProjectionY()
h_fwd_phi_lab_den_acc = h3_fwd_lab_den_acc.ProjectionY()
h_fwd_phi_hx_num_acc = h3_fwd_hx_num_acc.ProjectionY()
h_fwd_phi_hx_den_acc = h3_fwd_hx_den_acc.ProjectionY()
h_fwd_phi_cs_num_acc = h3_fwd_cs_num_acc.ProjectionY()
h_fwd_phi_cs_den_acc = h3_fwd_cs_den_acc.ProjectionY()

h_mid_phi_lab_num_acc = h3_mid_lab_num_acc.ProjectionY()
h_mid_phi_lab_den_acc = h3_mid_lab_den_acc.ProjectionY()
h_mid_phi_hx_num_acc = h3_mid_hx_num_acc.ProjectionY()
h_mid_phi_hx_den_acc = h3_mid_hx_den_acc.ProjectionY()
h_mid_phi_cs_num_acc = h3_mid_cs_num_acc.ProjectionY()
h_mid_phi_cs_den_acc = h3_mid_cs_den_acc.ProjectionY()

# 1d map - pT
h_fwd_pT_lab_num_acc = h3_fwd_lab_num_acc.ProjectionZ()
h_fwd_pT_lab_den_acc = h3_fwd_lab_den_acc.ProjectionZ()
h_fwd_pT_hx_num_acc = h3_fwd_hx_num_acc.ProjectionZ()
h_fwd_pT_hx_den_acc = h3_fwd_hx_den_acc.ProjectionZ()
h_fwd_pT_cs_num_acc = h3_fwd_cs_num_acc.ProjectionZ()
h_fwd_pT_cs_den_acc = h3_fwd_cs_den_acc.ProjectionZ()

h_mid_pT_lab_num_acc = h3_mid_lab_num_acc.ProjectionZ()
h_mid_pT_lab_den_acc = h3_mid_lab_den_acc.ProjectionZ()
h_mid_pT_hx_num_acc = h3_mid_hx_num_acc.ProjectionZ()
h_mid_pT_hx_den_acc = h3_mid_hx_den_acc.ProjectionZ()
h_mid_pT_cs_num_acc = h3_mid_cs_num_acc.ProjectionZ()
h_mid_pT_cs_den_acc = h3_mid_cs_den_acc.ProjectionZ()

# 2d map - cos vs pT
h2_fwd_cos_pT_lab_num_acc = h3_fwd_lab_num_acc.Project3D("zx")
h2_fwd_cos_pT_lab_den_acc = h3_fwd_lab_den_acc.Project3D("zx")
h2_fwd_cos_pT_hx_num_acc = h3_fwd_hx_num_acc.Project3D("zx")
h2_fwd_cos_pT_hx_den_acc = h3_fwd_hx_den_acc.Project3D("zx")
h2_fwd_cos_pT_cs_num_acc = h3_fwd_cs_num_acc.Project3D("zx")
h2_fwd_cos_pT_cs_den_acc = h3_fwd_cs_den_acc.Project3D("zx")

h2_mid_cos_pT_lab_num_acc = h3_mid_lab_num_acc.Project3D("zx")
h2_mid_cos_pT_lab_den_acc = h3_mid_lab_den_acc.Project3D("zx")
h2_mid_cos_pT_hx_num_acc = h3_mid_hx_num_acc.Project3D("zx")
h2_mid_cos_pT_hx_den_acc = h3_mid_hx_den_acc.Project3D("zx")
h2_mid_cos_pT_cs_num_acc = h3_mid_cs_num_acc.Project3D("zx")
h2_mid_cos_pT_cs_den_acc = h3_mid_cs_den_acc.Project3D("zx")


# ===== read inputs - Eff ===== #
# bring input root file
in_file_eff = TFile('roots/eff_pp_Jpsi_PtW1_tnp1_all_event.root')

# ===== get hists ===== #
# 3D map
h3_fwd_lab_num_eff = in_file_eff.Get('fwd_lab_num')
h3_fwd_lab_den_eff = in_file_eff.Get('fwd_lab_den')
h3_fwd_hx_num_eff = in_file_eff.Get('fwd_hx_num')
h3_fwd_hx_den_eff = in_file_eff.Get('fwd_hx_den')
h3_fwd_cs_num_eff = in_file_eff.Get('fwd_cs_num')
h3_fwd_cs_den_eff = in_file_eff.Get('fwd_cs_den')

h3_mid_lab_num_eff = in_file_eff.Get('mid_lab_num')
h3_mid_lab_den_eff = in_file_eff.Get('mid_lab_den')
h3_mid_hx_num_eff = in_file_eff.Get('mid_hx_num')
h3_mid_hx_den_eff = in_file_eff.Get('mid_hx_den')
h3_mid_cs_num_eff = in_file_eff.Get('mid_cs_num')
h3_mid_cs_den_eff = in_file_eff.Get('mid_cs_den')

# 1d map - y
h_y_lab_num_eff = in_file_eff.Get('y_lab_num')
h_y_lab_den_eff = in_file_eff.Get('y_lab_den')

# 1d map - cos
h_fwd_cos_lab_num_eff = h3_fwd_lab_num_eff.ProjectionX()
h_fwd_cos_lab_den_eff = h3_fwd_lab_den_eff.ProjectionX()
h_fwd_cos_hx_num_eff = h3_fwd_hx_num_eff.ProjectionX()
h_fwd_cos_hx_den_eff = h3_fwd_hx_den_eff.ProjectionX()
h_fwd_cos_cs_num_eff = h3_fwd_cs_num_eff.ProjectionX()
h_fwd_cos_cs_den_eff = h3_fwd_cs_den_eff.ProjectionX()

h_mid_cos_lab_num_eff = h3_mid_lab_num_eff.ProjectionX()
h_mid_cos_lab_den_eff = h3_mid_lab_den_eff.ProjectionX()
h_mid_cos_hx_num_eff = h3_mid_hx_num_eff.ProjectionX()
h_mid_cos_hx_den_eff = h3_mid_hx_den_eff.ProjectionX()
h_mid_cos_cs_num_eff = h3_mid_cs_num_eff.ProjectionX()
h_mid_cos_cs_den_eff = h3_mid_cs_den_eff.ProjectionX()

# 1d map - phi
h_fwd_phi_lab_num_eff = h3_fwd_lab_num_eff.ProjectionY()
h_fwd_phi_lab_den_eff = h3_fwd_lab_den_eff.ProjectionY()
h_fwd_phi_hx_num_eff = h3_fwd_hx_num_eff.ProjectionY()
h_fwd_phi_hx_den_eff = h3_fwd_hx_den_eff.ProjectionY()
h_fwd_phi_cs_num_eff = h3_fwd_cs_num_eff.ProjectionY()
h_fwd_phi_cs_den_eff = h3_fwd_cs_den_eff.ProjectionY()

h_mid_phi_lab_num_eff = h3_mid_lab_num_eff.ProjectionY()
h_mid_phi_lab_den_eff = h3_mid_lab_den_eff.ProjectionY()
h_mid_phi_hx_num_eff = h3_mid_hx_num_eff.ProjectionY()
h_mid_phi_hx_den_eff = h3_mid_hx_den_eff.ProjectionY()
h_mid_phi_cs_num_eff = h3_mid_cs_num_eff.ProjectionY()
h_mid_phi_cs_den_eff = h3_mid_cs_den_eff.ProjectionY()

# 1d map - pT
h_fwd_pT_lab_num_eff = h3_fwd_lab_num_eff.ProjectionZ()
h_fwd_pT_lab_den_eff = h3_fwd_lab_den_eff.ProjectionZ()
h_fwd_pT_hx_num_eff = h3_fwd_hx_num_eff.ProjectionZ()
h_fwd_pT_hx_den_eff = h3_fwd_hx_den_eff.ProjectionZ()
h_fwd_pT_cs_num_eff = h3_fwd_cs_num_eff.ProjectionZ()
h_fwd_pT_cs_den_eff = h3_fwd_cs_den_eff.ProjectionZ()

h_mid_pT_lab_num_eff = h3_mid_lab_num_eff.ProjectionZ()
h_mid_pT_lab_den_eff = h3_mid_lab_den_eff.ProjectionZ()
h_mid_pT_hx_num_eff = h3_mid_hx_num_eff.ProjectionZ()
h_mid_pT_hx_den_eff = h3_mid_hx_den_eff.ProjectionZ()
h_mid_pT_cs_num_eff = h3_mid_cs_num_eff.ProjectionZ()
h_mid_pT_cs_den_eff = h3_mid_cs_den_eff.ProjectionZ()

# 2d map - cos vs pT
h2_fwd_cos_pT_lab_num_eff = h3_fwd_lab_num_eff.Project3D("zx")
h2_fwd_cos_pT_lab_den_eff = h3_fwd_lab_den_eff.Project3D("zx")
h2_fwd_cos_pT_hx_num_eff = h3_fwd_hx_num_eff.Project3D("zx")
h2_fwd_cos_pT_hx_den_eff = h3_fwd_hx_den_eff.Project3D("zx")
h2_fwd_cos_pT_cs_num_eff = h3_fwd_cs_num_eff.Project3D("zx")
h2_fwd_cos_pT_cs_den_eff = h3_fwd_cs_den_eff.Project3D("zx")

h2_mid_cos_pT_lab_num_eff = h3_mid_lab_num_eff.Project3D("zx")
h2_mid_cos_pT_lab_den_eff = h3_mid_lab_den_eff.Project3D("zx")
h2_mid_cos_pT_hx_num_eff = h3_mid_hx_num_eff.Project3D("zx")
h2_mid_cos_pT_hx_den_eff = h3_mid_hx_den_eff.Project3D("zx")
h2_mid_cos_pT_cs_num_eff = h3_mid_cs_num_eff.Project3D("zx")
h2_mid_cos_pT_cs_den_eff = h3_mid_cs_den_eff.Project3D("zx")


# ===== draw 1d plots ===== #
# y
draw_1d_cor(h_acc_num=h_y_lab_num_acc, h_acc_den=h_y_lab_den_acc,
            h_eff_num=h_y_lab_num_eff, h_eff_den=h_y_lab_den_eff, 
            legend1='', legend2='', x_title='|y|', y_title='', save_name='tot_y_lab_tot')

# cos
draw_1d_cor(h_acc_num=h_fwd_cos_lab_num_acc, h_acc_den=h_fwd_cos_lab_den_acc,
            h_eff_num=h_fwd_cos_lab_num_eff, h_eff_den=h_fwd_cos_lab_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_fwd/cos_lab_tot')
draw_1d_cor(h_acc_num=h_fwd_cos_hx_num_acc, h_acc_den=h_fwd_cos_hx_den_acc,
            h_eff_num=h_fwd_cos_hx_num_eff, h_eff_den=h_fwd_cos_hx_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_fwd/cos_hx_tot')
draw_1d_cor(h_acc_num=h_fwd_cos_cs_num_acc, h_acc_den=h_fwd_cos_cs_den_acc,
            h_eff_num=h_fwd_cos_cs_num_eff, h_eff_den=h_fwd_cos_cs_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_fwd/cos_cs_tot')

draw_1d_cor(h_acc_num=h_mid_cos_lab_num_acc, h_acc_den=h_mid_cos_lab_den_acc,
            h_eff_num=h_mid_cos_lab_num_eff, h_eff_den=h_mid_cos_lab_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_mid/cos_lab_tot')
draw_1d_cor(h_acc_num=h_mid_cos_hx_num_acc, h_acc_den=h_mid_cos_hx_den_acc,
            h_eff_num=h_mid_cos_hx_num_eff, h_eff_den=h_mid_cos_hx_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_mid/cos_hx_tot')
draw_1d_cor(h_acc_num=h_mid_cos_cs_num_acc, h_acc_den=h_mid_cos_cs_den_acc,
            h_eff_num=h_mid_cos_cs_num_eff, h_eff_den=h_mid_cos_cs_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='', save_name='tot_mid/cos_cs_tot')

# # phi
draw_1d_cor(h_acc_num=h_fwd_phi_lab_num_acc, h_acc_den=h_fwd_phi_lab_den_acc,
            h_eff_num=h_fwd_phi_lab_num_eff, h_eff_den=h_fwd_phi_lab_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_fwd/phi_lab_tot')
draw_1d_cor(h_acc_num=h_fwd_phi_hx_num_acc, h_acc_den=h_fwd_phi_hx_den_acc,
            h_eff_num=h_fwd_phi_hx_num_eff, h_eff_den=h_fwd_phi_hx_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_fwd/phi_hx_tot')
draw_1d_cor(h_acc_num=h_fwd_phi_cs_num_acc, h_acc_den=h_fwd_phi_cs_den_acc,
            h_eff_num=h_fwd_phi_cs_num_eff, h_eff_den=h_fwd_phi_cs_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_fwd/phi_cs_tot')

draw_1d_cor(h_acc_num=h_mid_phi_lab_num_acc, h_acc_den=h_mid_phi_lab_den_acc,
            h_eff_num=h_mid_phi_lab_num_eff, h_eff_den=h_mid_phi_lab_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_mid/phi_lab_tot')
draw_1d_cor(h_acc_num=h_mid_phi_hx_num_acc, h_acc_den=h_mid_phi_hx_den_acc,
            h_eff_num=h_mid_phi_hx_num_eff, h_eff_den=h_mid_phi_hx_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_mid/phi_hx_tot')
draw_1d_cor(h_acc_num=h_mid_phi_cs_num_acc, h_acc_den=h_mid_phi_cs_den_acc,
            h_eff_num=h_mid_phi_cs_num_eff, h_eff_den=h_mid_phi_cs_den_eff, 
            legend1='', legend2='', x_title='#phi', y_title='', save_name='tot_mid/phi_cs_tot')

# # pt
draw_1d_cor(h_acc_num=h_fwd_pT_lab_num_acc, h_acc_den=h_fwd_pT_lab_den_acc,
            h_eff_num=h_fwd_pT_lab_num_eff, h_eff_den=h_fwd_pT_lab_den_eff, 
            legend1='', legend2='', x_title='p_{T}', y_title='', save_name='tot_fwd/pT_lab_tot')
draw_1d_cor(h_acc_num=h_mid_pT_cs_num_acc, h_acc_den=h_mid_pT_cs_den_acc,
            h_eff_num=h_mid_pT_cs_num_eff, h_eff_den=h_mid_pT_cs_den_eff, 
            legend1='', legend2='', x_title='p_{T}', y_title='', save_name='tot_mid/pT_cs_tot')


# ===== draw 2d plots ===== #
# cos vs pT
draw_2d_cor(h_acc_num=h2_fwd_cos_pT_lab_num_acc, h_acc_den=h2_fwd_cos_pT_lab_den_acc,
            h_eff_num=h2_fwd_cos_pT_lab_num_eff, h_eff_den=h2_fwd_cos_pT_lab_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_fwd/2d_cos_pT_lab')
draw_2d_cor(h_acc_num=h2_fwd_cos_pT_hx_num_acc, h_acc_den=h2_fwd_cos_pT_hx_den_acc,
            h_eff_num=h2_fwd_cos_pT_hx_num_eff, h_eff_den=h2_fwd_cos_pT_hx_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_fwd/2d_cos_pT_hx')
draw_2d_cor(h_acc_num=h2_fwd_cos_pT_cs_num_acc, h_acc_den=h2_fwd_cos_pT_cs_den_acc,
            h_eff_num=h2_fwd_cos_pT_cs_num_eff, h_eff_den=h2_fwd_cos_pT_cs_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_fwd/2d_cos_pT_cs')

draw_2d_cor(h_acc_num=h2_mid_cos_pT_lab_num_acc, h_acc_den=h2_mid_cos_pT_lab_den_acc,
            h_eff_num=h2_mid_cos_pT_lab_num_eff, h_eff_den=h2_mid_cos_pT_lab_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_mid/2d_cos_pT_lab')
draw_2d_cor(h_acc_num=h2_mid_cos_pT_hx_num_acc, h_acc_den=h2_mid_cos_pT_hx_den_acc,
            h_eff_num=h2_mid_cos_pT_hx_num_eff, h_eff_den=h2_mid_cos_pT_hx_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_mid/2d_cos_pT_hx')
draw_2d_cor(h_acc_num=h2_mid_cos_pT_cs_num_acc, h_acc_den=h2_mid_cos_pT_cs_den_acc,
            h_eff_num=h2_mid_cos_pT_cs_num_eff, h_eff_den=h2_mid_cos_pT_cs_den_eff, 
            legend1='', legend2='', x_title='cos#theta', y_title='p_{T}', save_name='tot_mid/2d_cos_pT_cs')


# ===== Fwd cos vs phi - accroding to pT bins ===== # - 6
# pT integrated
h2_fwd_cos_phi_lab_num_acc = h3_fwd_lab_num_acc.Project3D("yx")
h2_fwd_cos_phi_lab_den_acc = h3_fwd_lab_den_acc.Project3D("yx")
h2_fwd_cos_phi_hx_num_acc = h3_fwd_hx_num_acc.Project3D("yx")
h2_fwd_cos_phi_hx_den_acc = h3_fwd_hx_den_acc.Project3D("yx")
h2_fwd_cos_phi_cs_num_acc = h3_fwd_cs_num_acc.Project3D("yx")
h2_fwd_cos_phi_cs_den_acc = h3_fwd_cs_den_acc.Project3D("yx")

h2_fwd_cos_phi_lab_num_eff = h3_fwd_lab_num_eff.Project3D("yx")
h2_fwd_cos_phi_lab_den_eff = h3_fwd_lab_den_eff.Project3D("yx")
h2_fwd_cos_phi_hx_num_eff = h3_fwd_hx_num_eff.Project3D("yx")
h2_fwd_cos_phi_hx_den_eff = h3_fwd_hx_den_eff.Project3D("yx")
h2_fwd_cos_phi_cs_num_eff = h3_fwd_cs_num_eff.Project3D("yx")
h2_fwd_cos_phi_cs_den_eff = h3_fwd_cs_den_eff.Project3D("yx")

draw_2d_cor(h_acc_num=h2_fwd_cos_phi_lab_num_acc, h_acc_den=h2_fwd_cos_phi_lab_den_acc,
            h_eff_num=h2_fwd_cos_phi_lab_num_eff, h_eff_den=h2_fwd_cos_phi_lab_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_lab')
draw_2d_cor(h_acc_num=h2_fwd_cos_phi_hx_num_acc, h_acc_den=h2_fwd_cos_phi_hx_den_acc,
            h_eff_num=h2_fwd_cos_phi_hx_num_eff, h_eff_den=h2_fwd_cos_phi_hx_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_hx')
draw_2d_cor(h_acc_num=h2_fwd_cos_phi_cs_num_acc, h_acc_den=h2_fwd_cos_phi_cs_den_acc,
            h_eff_num=h2_fwd_cos_phi_cs_num_eff, h_eff_den=h2_fwd_cos_phi_cs_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_cs')



for idx in (1,2,3,4,5,6):
    h3_fwd_lab_num_acc.GetZaxis().SetRange(idx, idx)
    h3_fwd_lab_den_acc.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_num_acc.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_den_acc.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_num_acc.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_den_acc.GetZaxis().SetRange(idx, idx)
    h2_fwd_cos_phi_lab_num_acc = h3_fwd_lab_num_acc.Project3D("yx")
    h2_fwd_cos_phi_lab_den_acc = h3_fwd_lab_den_acc.Project3D("yx")
    h2_fwd_cos_phi_hx_num_acc = h3_fwd_hx_num_acc.Project3D("yx")
    h2_fwd_cos_phi_hx_den_acc = h3_fwd_hx_den_acc.Project3D("yx")
    h2_fwd_cos_phi_cs_num_acc = h3_fwd_cs_num_acc.Project3D("yx")
    h2_fwd_cos_phi_cs_den_acc = h3_fwd_cs_den_acc.Project3D("yx")

    h3_fwd_lab_num_eff.GetZaxis().SetRange(idx, idx)
    h3_fwd_lab_den_eff.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_num_eff.GetZaxis().SetRange(idx, idx)
    h3_fwd_hx_den_eff.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_num_eff.GetZaxis().SetRange(idx, idx)
    h3_fwd_cs_den_eff.GetZaxis().SetRange(idx, idx)
    h2_fwd_cos_phi_lab_num_eff = h3_fwd_lab_num_eff.Project3D("yx")
    h2_fwd_cos_phi_lab_den_eff = h3_fwd_lab_den_eff.Project3D("yx")
    h2_fwd_cos_phi_hx_num_eff = h3_fwd_hx_num_eff.Project3D("yx")
    h2_fwd_cos_phi_hx_den_eff = h3_fwd_hx_den_eff.Project3D("yx")
    h2_fwd_cos_phi_cs_num_eff = h3_fwd_cs_num_eff.Project3D("yx")
    h2_fwd_cos_phi_cs_den_eff = h3_fwd_cs_den_eff.Project3D("yx")


    fwd_pt_bins = ('pT3_6', 'pT6_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')
    name_tag = fwd_pt_bins[idx-1]

    draw_2d_cor(h_acc_num=h2_fwd_cos_phi_lab_num_acc, h_acc_den=h2_fwd_cos_phi_lab_den_acc,
                h_eff_num=h2_fwd_cos_phi_lab_num_eff, h_eff_den=h2_fwd_cos_phi_lab_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_lab_' + name_tag)
    draw_2d_cor(h_acc_num=h2_fwd_cos_phi_hx_num_acc, h_acc_den=h2_fwd_cos_phi_hx_den_acc,
                h_eff_num=h2_fwd_cos_phi_hx_num_eff, h_eff_den=h2_fwd_cos_phi_hx_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_hx_' + name_tag)
    draw_2d_cor(h_acc_num=h2_fwd_cos_phi_cs_num_acc, h_acc_den=h2_fwd_cos_phi_cs_den_acc,
                h_eff_num=h2_fwd_cos_phi_cs_num_eff, h_eff_den=h2_fwd_cos_phi_cs_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_fwd/2d_cos_phi_cs_' + name_tag)


# # ===== Mid cos vs phi - accroding to pT bins ===== # - 5
# pT integrated
h2_mid_cos_phi_lab_num_acc = h3_mid_lab_num_acc.Project3D("yx")
h2_mid_cos_phi_lab_den_acc = h3_mid_lab_den_acc.Project3D("yx")
h2_mid_cos_phi_hx_num_acc = h3_mid_hx_num_acc.Project3D("yx")
h2_mid_cos_phi_hx_den_acc = h3_mid_hx_den_acc.Project3D("yx")
h2_mid_cos_phi_cs_num_acc = h3_mid_cs_num_acc.Project3D("yx")
h2_mid_cos_phi_cs_den_acc = h3_mid_cs_den_acc.Project3D("yx")

h2_mid_cos_phi_lab_num_eff = h3_mid_lab_num_eff.Project3D("yx")
h2_mid_cos_phi_lab_den_eff = h3_mid_lab_den_eff.Project3D("yx")
h2_mid_cos_phi_hx_num_eff = h3_mid_hx_num_eff.Project3D("yx")
h2_mid_cos_phi_hx_den_eff = h3_mid_hx_den_eff.Project3D("yx")
h2_mid_cos_phi_cs_num_eff = h3_mid_cs_num_eff.Project3D("yx")
h2_mid_cos_phi_cs_den_eff = h3_mid_cs_den_eff.Project3D("yx")

draw_2d_cor(h_acc_num=h2_mid_cos_phi_lab_num_acc, h_acc_den=h2_mid_cos_phi_lab_den_acc,
            h_eff_num=h2_mid_cos_phi_lab_num_eff, h_eff_den=h2_mid_cos_phi_lab_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_lab')
draw_2d_cor(h_acc_num=h2_mid_cos_phi_hx_num_acc, h_acc_den=h2_mid_cos_phi_hx_den_acc,
            h_eff_num=h2_mid_cos_phi_hx_num_eff, h_eff_den=h2_mid_cos_phi_hx_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_hx')
draw_2d_cor(h_acc_num=h2_mid_cos_phi_cs_num_acc, h_acc_den=h2_mid_cos_phi_cs_den_acc,
            h_eff_num=h2_mid_cos_phi_cs_num_eff, h_eff_den=h2_mid_cos_phi_cs_den_eff,
            legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_cs')


for idx in (1,2,3,4,5):
    h3_mid_lab_num_acc.GetZaxis().SetRange(idx, idx)
    h3_mid_lab_den_acc.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_num_acc.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_den_acc.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_num_acc.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_den_acc.GetZaxis().SetRange(idx, idx)
    h2_mid_cos_phi_lab_num_acc = h3_mid_lab_num_acc.Project3D("yx")
    h2_mid_cos_phi_lab_den_acc = h3_mid_lab_den_acc.Project3D("yx")
    h2_mid_cos_phi_hx_num_acc = h3_mid_hx_num_acc.Project3D("yx")
    h2_mid_cos_phi_hx_den_acc = h3_mid_hx_den_acc.Project3D("yx")
    h2_mid_cos_phi_cs_num_acc = h3_mid_cs_num_acc.Project3D("yx")
    h2_mid_cos_phi_cs_den_acc = h3_mid_cs_den_acc.Project3D("yx")

    h3_mid_lab_num_eff.GetZaxis().SetRange(idx, idx)
    h3_mid_lab_den_eff.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_num_eff.GetZaxis().SetRange(idx, idx)
    h3_mid_hx_den_eff.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_num_eff.GetZaxis().SetRange(idx, idx)
    h3_mid_cs_den_eff.GetZaxis().SetRange(idx, idx)
    h2_mid_cos_phi_lab_num_eff = h3_mid_lab_num_eff.Project3D("yx")
    h2_mid_cos_phi_lab_den_eff = h3_mid_lab_den_eff.Project3D("yx")
    h2_mid_cos_phi_hx_num_eff = h3_mid_hx_num_eff.Project3D("yx")
    h2_mid_cos_phi_hx_den_eff = h3_mid_hx_den_eff.Project3D("yx")
    h2_mid_cos_phi_cs_num_eff = h3_mid_cs_num_eff.Project3D("yx")
    h2_mid_cos_phi_cs_den_eff = h3_mid_cs_den_eff.Project3D("yx")


    mid_pt_bins = ('pT6p5_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')
    name_tag = mid_pt_bins[idx-1]

    draw_2d_cor(h_acc_num=h2_mid_cos_phi_lab_num_acc, h_acc_den=h2_mid_cos_phi_lab_den_acc,
                h_eff_num=h2_mid_cos_phi_lab_num_eff, h_eff_den=h2_mid_cos_phi_lab_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_lab_' + name_tag)
    draw_2d_cor(h_acc_num=h2_mid_cos_phi_hx_num_acc, h_acc_den=h2_mid_cos_phi_hx_den_acc,
                h_eff_num=h2_mid_cos_phi_hx_num_eff, h_eff_den=h2_mid_cos_phi_hx_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_hx_' + name_tag)
    draw_2d_cor(h_acc_num=h2_mid_cos_phi_cs_num_acc, h_acc_den=h2_mid_cos_phi_cs_den_acc,
                h_eff_num=h2_mid_cos_phi_cs_num_eff, h_eff_den=h2_mid_cos_phi_cs_den_eff,
                legend1='', legend2='', x_title='cos#theta', y_title='#phi', save_name='tot_mid/2d_cos_phi_cs_' + name_tag)