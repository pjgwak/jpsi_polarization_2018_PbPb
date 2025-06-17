import os, ROOT
from Plot import PlotOverlaid
from ROOT import TFile, TCanvas, TLegend, kGreen, kBlack, TCanvas, TH1


def plotOverlaid(fileName='', baseName='Data2023', formats=['png'], sampleInfo='', isWeight=False, isEP=False):
    """ draw overlaid plots according to 3 rapidity regions
    """

    # create saving folders
    os.makedirs(f'figs/{baseName}_overlaid', exist_ok=True)

    f = TFile.Open(fileName)

    canvas = TCanvas('c_massOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)
    h_epAngHFm2 = f.Get('h_epAngHFm2')
    h_epAngHFm2.SetTitle(';#Psi_{2} (radian);')
    h_epAngHFm2.SetLineColor(ROOT.kBlue)
    h_epAngHFm2.SetLineWidth(2)
    h_epAngHFm2.SetMarkerColor(ROOT.kBlue)
    h_epAngHFm2.SetMarkerStyle(20)

    # draw and save
    h_epAngHFm2.Draw()
    canvas.SaveAs(f'figs/{baseName}_overlaid/epAng_HFm2.png')

    #
    h_epAngHFm2 = f.Get('h_epAngHFp2')
    h_epAngHFm2.SetTitle(';#Psi_{2} (radian);')
    h_epAngHFm2.SetLineColor(ROOT.kBlue)
    h_epAngHFm2.SetLineWidth(2)
    h_epAngHFm2.SetMarkerColor(ROOT.kBlue)
    h_epAngHFm2.SetMarkerStyle(20)

    # draw and save
    h_epAngHFm2.Draw()
    canvas.SaveAs(f'figs/{baseName}_overlaid/epAng_HFp2.png')

    #
    h_epAngHFm2 = f.Get('h_epAngQxQyHFp2')
    h_epAngHFm2.SetTitle(';#Psi_{2} (radian);')
    h_epAngHFm2.SetLineColor(ROOT.kBlue)
    h_epAngHFm2.SetLineWidth(2)
    h_epAngHFm2.SetMarkerColor(ROOT.kBlue)
    h_epAngHFm2.SetMarkerStyle(20)

    # draw and save
    h_epAngHFm2.Draw()
    canvas.SaveAs(f'figs/{baseName}_overlaid/QxQy_HFp2.png')

    #
    h_epAngHFm2 = f.Get('h_epAngQxQyHFm2')
    h_epAngHFm2.SetTitle(';#Psi_{2} (radian);')
    h_epAngHFm2.SetLineColor(ROOT.kBlue)
    h_epAngHFm2.SetLineWidth(2)
    h_epAngHFm2.SetMarkerColor(ROOT.kBlue)
    h_epAngHFm2.SetMarkerStyle(20)

    # draw and save
    h_epAngHFm2.Draw()
    canvas.SaveAs(f'figs/{baseName}_overlaid/QxQy_HFm2.png')


    # cleaning
    f.Close()