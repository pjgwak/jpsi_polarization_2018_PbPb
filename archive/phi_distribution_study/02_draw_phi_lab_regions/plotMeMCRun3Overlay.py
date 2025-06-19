import os, ROOT
from Plot import PlotOverlaid
from ROOT import TFile, TCanvas, TLegend, kGreen, kBlack, TCanvas, TH1


def plotOverlaid(fileName='', baseName='Data2023', formats=['png'], sampleInfo='', isWeight=False, isEP=False):
    """ draw overlaid plots according to 3 rapidity regions
    """

    # create saving folders
    os.makedirs(f'figs/{baseName}_overlaid', exist_ok=True)

    f = TFile.Open(fileName)
    myTree = f.Get('myTree')

    cutSignalMass = '&& (mass > 2.9 && mass < 3.3)'
    cutBaseFwd = 'pt > 3 && pt < 50 && fabs(y) > 1.6 && fabs(y) < 2.4' #+ cutSignalMass
    cutBaseMid = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 1.6'  #+ cutSignalMass
    cutBaseAllY = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 2.4' #+ cutSignalMass

    for idx in range (-32, 32, 2):
        idx_low = idx / 10
        idx_high = idx_low + 0.2

        cutPhiLab = f'&& (phi>{idx_low}) && phi<{idx_high}'
        cutBaseFwd = 'pt > 3 && pt < 50 && (fabs(y) > 1.6 && fabs(y) < 2.4)' + cutPhiLab
        cutBaseMid = 'pt > 6.5 && pt < 50 && (fabs(y) > 0 && fabs(y) < 1.6)'  + cutPhiLab
        cutBaseAllY = 'pt > 6.5 && pt < 50 && (fabs(y) > 0 && fabs(y) < 2.4)' + cutPhiLab

        canvas = TCanvas(f'c_{idx_low}_{idx_high}', '', 800, 600)
        leg = TLegend(0.45, 0.80, 0.9, 0.92)

        h_fwd = ROOT.TH1D(f'h_fwd_{idx_low}_{idx_high}', ';#varphi_{HX} (radian);', 40, -3.141592, 3.141592)
        h_fwd.SetLineColor(ROOT.kBlue)
        h_fwd.SetLineWidth(2)
        h_fwd.SetMarkerColor(ROOT.kBlue)
        h_fwd.SetMarkerStyle(20)

        h_mid = ROOT.TH1D(f'h_mid_{idx_low}_{idx_high}', '', 40, -3.141592, 3.141592)
        h_mid.SetLineColor(ROOT.kRed)
        h_mid.SetLineWidth(2)
        h_mid.SetMarkerColor(ROOT.kRed)
        h_mid.SetMarkerStyle(21)

        h_all = ROOT.TH1D(f'h_all_{idx_low}_{idx_high}', '', 40, -3.141592, 3.141592)
        h_all.SetLineColor(ROOT.kGreen + 2)
        h_all.SetLineWidth(2)
        h_all.SetMarkerColor(ROOT.kGreen + 2)
        h_all.SetMarkerStyle(22)


        # draw and save
        myPlot = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

        myPlot.add(tree=myTree, var='phiHX', cut=cutBaseFwd, hist=h_fwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50, '+ f'{idx_low:.1f}<'+'#varphi_{Lab}'+f'<{idx_high:.1f}', drawOpt='ep')
        myPlot.add(tree=myTree, var='phiHX', cut=cutBaseMid, hist=h_mid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50, '+ f'{idx_low:.1f}<'+'#varphi_{Lab}'+f'<{idx_high:.1f}', drawOpt='ep')
        myPlot.add(tree=myTree, var='phiHX', cut=cutBaseAllY, hist=h_all, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50, '+ f'{idx_low:.1f}<'+'#varphi_{Lab}'+f'<{idx_high:.1f}', drawOpt='ep')

        myPlot.drawRun3MC()
        myPlot.save(f'figs/{baseName}_overlaid/{idx_low:.1f}_{idx_high:.1f}', formats=formats)


    # cleaning
    f.Close()