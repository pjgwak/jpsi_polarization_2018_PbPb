import ROOT
from ROOT import TFile
# TFile, TCanvas, RooWorkspace, TH2D, RooFit, RooRealVar

# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# batch mode - dont' draw the plot
ROOT.gROOT.SetBatch(True)

# bring input root file
in_file = TFile('../../eff_acc/roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_250221.root')
# in_file.Print('V')

# bring histograms
h_mid_cent0_20 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy0_1p6')
h_fwd_cent0_20 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy1p6_2p4')

h_mid_cent20_60 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy0_1p6')
h_fwd_cent20_60 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy1p6_2p4')

h_mid_cent60_180 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy0_1p6')
h_fwd_cent60_180 = in_file.Get('mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4')


# define a fcn for drawing and saving
def draw_plot(hist2D, title='', save_name='', label=''):
    canvas = ROOT.TCanvas("canvas", "", 1200, 900)

    hist2D.Draw('colz')
    hist2D.SetTitle(title)

    # [hist, "Title", "saving name", "Kinematics label"

    texts = []  # python need it to store the objects

    for i in range(1, hist2D.GetNbinsX() + 1):
        for j in range(1, hist2D.GetNbinsY() + 1):
            binContent = hist2D.GetBinContent(i, j)
            if binContent > 0:  # when content > 0
                x = hist2D.GetXaxis().GetBinCenter(i)
                y = hist2D.GetYaxis().GetBinCenter(j)
                
                text = ROOT.TText(x, y, f"{binContent:.2f}")
                text.SetTextSize(0.02)
                text.SetTextColor(ROOT.kBlack)
                text.SetTextAlign(22)  # draw value in the center of bin
                text.Draw()
                texts.append(text)  # python need it! Or text object will be overwritten by call.

    # draw labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(0.89, 0.80, 'J/#psi#rightarrow#mu^{+}#mu^{-}')
    latex.DrawLatexNDC(0.89, 0.85, label)

    canvas.Update()
    canvas.Draw()
    if save_name != '':
        canvas.SaveAs('figs/' + save_name)


# ===== drawing ===== #
draw_plot(h_mid_cent0_20, 'Acc x Eff', 'mid_cent0_20.png', 
    label = '0 < |y| < 1.6, 0 < cent.(%) < 20')
draw_plot(h_fwd_cent0_20, 'Acc x Eff', 'fwd_cent0_20.png', 
    label = '1.6 < |y| < 2.4, 0 < cent.(%) < 10')

draw_plot(h_mid_cent20_60, 'Acc x Eff', 'mid_cent20_60.png', 
    label = '0 < |y| < 1.6, 10 < cent.(%) < 30')
draw_plot(h_fwd_cent20_60, 'Acc x Eff', 'fwd_cent20_60.png', 
    label = '1.6 < |y| < 2.4, 10 < cent.(%) < 30')

draw_plot(h_mid_cent60_180, 'Acc x Eff', 'mid_cent60_180.png', 
    label = '0 < |y| < 1.6, 30 < cent.(%) < 90')
draw_plot(h_fwd_cent60_180, 'Acc x Eff', 'fwd_cent60_180.png', 
    label = '1.6 < |y| < 2.4, 30 < cent.(%) < 90')