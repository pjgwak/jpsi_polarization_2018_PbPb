import os, ROOT
from Plot import PlotOverlaid
from ROOT import TFile, TCanvas, TLegend, kGreen, kBlack, TCanvas, TH1


def plotOverlaid(fileName='', baseName='Data2023', formats=['png'], sampleInfo='', isWeight=False):
    """ draw overlaid plots according to 3 rapidity regions
    """

    # create saving folders
    os.makedirs(f'figs/{baseName}_overlaid', exist_ok=True)


    # get data TTree
    f = TFile.Open(fileName)
    myTree = f.Get('myTree')

    cutBaseFwd = 'pt > 3 && pt < 50 && fabs(y) > 1.6 && fabs(y) < 2.4'
    cutBaseMid = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 1.6'
    cutBaseAllY = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 2.4'

    # ===== mass ===== #    
    canvas = TCanvas('c_massOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_massFwd = ROOT.TH1D('h_'+baseName+'Fwd_Mass', ';mass (GeV/c^{2});# of dimuon', 30, 2.6, 3.5)
    h_massFwd.SetLineColor(ROOT.kBlue)
    h_massFwd.SetLineWidth(2)
    h_massFwd.SetMarkerColor(ROOT.kBlue)
    h_massFwd.SetMarkerStyle(20)

    h_massMid = ROOT.TH1D('h_'+baseName+'Mid_Mass', '', 30, 2.6, 3.5)
    h_massMid.SetLineColor(ROOT.kRed)
    h_massMid.SetLineWidth(2)
    h_massMid.SetMarkerColor(ROOT.kRed)
    h_massMid.SetMarkerStyle(21)

    h_massAllY = ROOT.TH1D('h_'+baseName+'AllY_Mass', '', 30, 2.6, 3.5)
    h_massAllY.SetLineColor(ROOT.kGreen + 2)
    h_massAllY.SetLineWidth(2)
    h_massAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_massAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotMass = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotMass.add(tree=myTree, var='mass', cut=cutBaseFwd, hist=h_massFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotMass.add(tree=myTree, var='mass', cut=cutBaseMid, hist=h_massMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotMass.add(tree=myTree, var='mass', cut=cutBaseAllY, hist=h_massAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotMass.drawRun3MC()
    plotMass.save(f'figs/{baseName}_overlaid/mass', formats=formats)


    # ===== pt ===== #
    canvas = TCanvas('c_PtOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_ptFwd = ROOT.TH1D('h_'+baseName+'Fwd_Pt', ';#it{p}_{T} (GeV/c);# of dimuon', 100, 0, 50)
    h_ptFwd.SetLineColor(ROOT.kBlue)
    h_ptFwd.SetLineWidth(2)
    h_ptFwd.SetMarkerColor(ROOT.kBlue)
    h_ptFwd.SetMarkerStyle(20)

    h_ptMid = ROOT.TH1D('h_'+baseName+'Mid_Pt', '', 100, 0, 50)
    h_ptMid.SetLineColor(ROOT.kRed)
    h_ptMid.SetLineWidth(2)
    h_ptMid.SetMarkerColor(ROOT.kRed)
    h_ptMid.SetMarkerStyle(21)

    h_ptAllY = ROOT.TH1D('h_'+baseName+'AllY_Pt', '', 100, 0, 50)
    h_ptAllY.SetLineColor(ROOT.kGreen + 2)
    h_ptAllY.SetLineWidth(2)
    h_ptAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_ptAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotPt = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotPt.add(tree=myTree, var='pt', cut=cutBaseFwd, hist=h_ptFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotPt.add(tree=myTree, var='pt', cut=cutBaseMid, hist=h_ptMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotPt.add(tree=myTree, var='pt', cut=cutBaseAllY, hist=h_ptAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotPt.drawRun3MC()
    plotPt.save(f'figs/{baseName}_overlaid/pt', formats=formats)


    # ===== y ===== #
    canvas = TCanvas('c_YOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_Fwd = ROOT.TH1D('h_'+baseName+'Fwd_Y', ';Rapidity;# of dimuon', 48, -2.4, 2.4)
    h_Fwd.SetLineColor(ROOT.kBlue)
    h_Fwd.SetLineWidth(2)
    h_Fwd.SetMarkerColor(ROOT.kBlue)
    h_Fwd.SetMarkerStyle(20)

    h_Mid = ROOT.TH1D('h_'+baseName+'Mid_Y', '', 48, -2.4, 2.4)
    h_Mid.SetLineColor(ROOT.kRed)
    h_Mid.SetLineWidth(2)
    h_Mid.SetMarkerColor(ROOT.kRed)
    h_Mid.SetMarkerStyle(21)

    h_AllY = ROOT.TH1D('h_'+baseName+'AllY_Y', '', 48, -2.4, 2.4)
    h_AllY.SetLineColor(ROOT.kGreen + 2)
    h_AllY.SetLineWidth(2)
    h_AllY.SetMarkerColor(ROOT.kGreen + 2)
    h_AllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotY = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotY.add(tree=myTree, var='y', cut=cutBaseFwd, hist=h_Fwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotY.add(tree=myTree, var='y', cut=cutBaseMid, hist=h_Mid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotY.add(tree=myTree, var='y', cut=cutBaseAllY, hist=h_AllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotY.drawRun3MC()
    plotY.save(f'figs/{baseName}_overlaid/y', formats=formats)


    # ===== phi ===== #
    canvas = TCanvas('c_phiOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)


    # set histograms
    h_phiFwd = ROOT.TH1D('h_'+baseName+'Fwd_Phi', ';#varphi_{Lab}^{#mu^{+}#mu^{-}} (radian);# of dimuon', 40, -3.141592, 3.141592)
    h_phiFwd.SetLineColor(ROOT.kBlue)
    h_phiFwd.SetLineWidth(2)
    h_phiFwd.SetMarkerColor(ROOT.kBlue)
    h_phiFwd.SetMarkerStyle(20)

    h_phiMid = ROOT.TH1D('h_'+baseName+'Mid_Phi', '', 40, -3.141592, 3.141592)
    h_phiMid.SetLineColor(ROOT.kRed)
    h_phiMid.SetLineWidth(2)
    h_phiMid.SetMarkerColor(ROOT.kRed)
    h_phiMid.SetMarkerStyle(21)

    h_phiAllY = ROOT.TH1D('h_'+baseName+'AllY_Phi', '', 40, -3.141592, 3.141592)
    h_phiAllY.SetLineColor(ROOT.kGreen + 2)
    h_phiAllY.SetLineWidth(2)
    h_phiAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_phiAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotPhi = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotPhi.add(tree=myTree, var='phi', cut=cutBaseFwd, hist=h_phiFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotPhi.add(tree=myTree, var='phi', cut=cutBaseMid, hist=h_phiMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotPhi.add(tree=myTree, var='phi', cut=cutBaseAllY, hist=h_phiAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotPhi.drawRun3MC()
    plotPhi.save(f'figs/{baseName}_overlaid/phi', formats=formats)

    # ===== phi - mupl ===== #
    canvas = TCanvas('c_phi1Overlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_phi1Fwd = ROOT.TH1D('h_'+baseName+'Fwd_Phi1', ';#varphi_{Lab}^{#mu^{+}} (radian);# of dimuon', 40, -3.141592, 3.141592)
    h_phi1Fwd.SetLineColor(ROOT.kBlue)
    h_phi1Fwd.SetLineWidth(2)
    h_phi1Fwd.SetMarkerColor(ROOT.kBlue)
    h_phi1Fwd.SetMarkerStyle(20)
    h_phi1Mid = ROOT.TH1D('h_'+baseName+'Mid_Phi1', '', 40, -3.141592, 3.141592)
    h_phi1Mid.SetLineColor(ROOT.kRed)
    h_phi1Mid.SetLineWidth(2)
    h_phi1Mid.SetMarkerColor(ROOT.kRed)
    h_phi1Mid.SetMarkerStyle(21)

    h_phi1AllY = ROOT.TH1D('h_'+baseName+'AllY_Phi1', '', 40, -3.141592, 3.141592)
    h_phi1AllY.SetLineColor(ROOT.kGreen + 2)
    h_phi1AllY.SetLineWidth(2)
    h_phi1AllY.SetMarkerColor(ROOT.kGreen + 2)
    h_phi1AllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotPhi1 = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotPhi1.add(tree=myTree, var='phi1', cut=cutBaseFwd, hist=h_phi1Fwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotPhi1.add(tree=myTree, var='phi1', cut=cutBaseMid, hist=h_phi1Mid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotPhi1.add(tree=myTree, var='phi1', cut=cutBaseAllY, hist=h_phi1AllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotPhi1.drawRun3MC()
    plotPhi1.save(f'figs/{baseName}_overlaid/phi1', formats=formats)


    # ===== phi HX ===== #
    canvas = TCanvas('c_phiHXOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_phiHXFwd = ROOT.TH1D('h_'+baseName+'Fwd_PhiHX', ';#varphi_{HX} (radian);# of dimuon', 40, -3.141592, 3.141592)
    h_phiHXFwd.SetLineColor(ROOT.kBlue)
    h_phiHXFwd.SetLineWidth(2)
    h_phiHXFwd.SetMarkerColor(ROOT.kBlue)
    h_phiHXFwd.SetMarkerStyle(20)

    h_phiHXMid = ROOT.TH1D('h_'+baseName+'Mid_PhiHX', '', 40, -3.141592, 3.141592)
    h_phiHXMid.SetLineColor(ROOT.kRed)
    h_phiHXMid.SetLineWidth(2)
    h_phiHXMid.SetMarkerColor(ROOT.kRed)
    h_phiHXMid.SetMarkerStyle(21)

    h_phiHXAllY = ROOT.TH1D('h_'+baseName+'AllY_PhiHX', '', 40, -3.141592, 3.141592)
    h_phiHXAllY.SetLineColor(ROOT.kGreen + 2)
    h_phiHXAllY.SetLineWidth(2)
    h_phiHXAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_phiHXAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotPhiHX = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotPhiHX.add(tree=myTree, var='phiHX', cut=cutBaseFwd, hist=h_phiHXFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotPhiHX.add(tree=myTree, var='phiHX', cut=cutBaseMid, hist=h_phiHXMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotPhiHX.add(tree=myTree, var='phiHX', cut=cutBaseAllY, hist=h_phiHXAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotPhiHX.drawRun3MC()
    plotPhiHX.save(f'figs/{baseName}_overlaid/phiHX', formats=formats)


    # ===== phi CS ===== #
    canvas = TCanvas('c_phiCSOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_phiCSFwd = ROOT.TH1D('h_'+baseName+'Fwd_PhiCS', ';#varphi_{CS} (radian);# of dimuon', 40, -3.141592, 3.141592)
    h_phiCSFwd.SetLineColor(ROOT.kBlue)
    h_phiCSFwd.SetLineWidth(2)
    h_phiCSFwd.SetMarkerColor(ROOT.kBlue)
    h_phiCSFwd.SetMarkerStyle(20)

    h_phiCSMid = ROOT.TH1D('h_'+baseName+'Mid_PhiCS', '', 40, -3.141592, 3.141592)
    h_phiCSMid.SetLineColor(ROOT.kRed)
    h_phiCSMid.SetLineWidth(2)
    h_phiCSMid.SetMarkerColor(ROOT.kRed)
    h_phiCSMid.SetMarkerStyle(21)

    h_phiCSAllY = ROOT.TH1D('h_'+baseName+'AllY_PhiCS', '', 40, -3.141592, 3.141592)
    h_phiCSAllY.SetLineColor(ROOT.kGreen + 2)
    h_phiCSAllY.SetLineWidth(2)
    h_phiCSAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_phiCSAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotPhiCS = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotPhiCS.add(tree=myTree, var='phiCS', cut=cutBaseFwd, hist=h_phiCSFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotPhiCS.add(tree=myTree, var='phiCS', cut=cutBaseMid, hist=h_phiCSMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotPhiCS.add(tree=myTree, var='phiCS', cut=cutBaseAllY, hist=h_phiCSAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotPhiCS.drawRun3MC()
    plotPhiCS.save(f'figs/{baseName}_overlaid/phiCS', formats=formats)


    # ===== cos HX ===== #
    canvas = TCanvas('c_cosHXOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_cosHXFwd = ROOT.TH1D('h_'+baseName+'Fwd_CosHX', ';cos#theta_{HX};# of dimuon', 30, -1, 1)
    h_cosHXFwd.SetLineColor(ROOT.kBlue)
    h_cosHXFwd.SetLineWidth(2)
    h_cosHXFwd.SetMarkerColor(ROOT.kBlue)
    h_cosHXFwd.SetMarkerStyle(20)

    h_cosHXMid = ROOT.TH1D('h_'+baseName+'Mid_CosHX', '', 30, -1, 1)
    h_cosHXMid.SetLineColor(ROOT.kRed)
    h_cosHXMid.SetLineWidth(2)
    h_cosHXMid.SetMarkerColor(ROOT.kRed)
    h_cosHXMid.SetMarkerStyle(21)
    h_cosHXAllY = ROOT.TH1D('h_'+baseName+'AllY_CosHX', '', 30, -1, 1)
    h_cosHXAllY.SetLineColor(ROOT.kGreen + 2)
    h_cosHXAllY.SetLineWidth(2)
    h_cosHXAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_cosHXAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotCosHX = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotCosHX.add(tree=myTree, var='cosHX', cut=cutBaseFwd, hist=h_cosHXFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotCosHX.add(tree=myTree, var='cosHX', cut=cutBaseMid, hist=h_cosHXMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotCosHX.add(tree=myTree, var='cosHX', cut=cutBaseAllY, hist=h_cosHXAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotCosHX.drawRun3MC()
    plotCosHX.save(f'figs/{baseName}_overlaid/cosHX', formats=formats)


    # ===== cos CS ===== #
    canvas = TCanvas('c_cosCSOverlay', '', 800, 600)
    leg = TLegend(0.65, 0.83, 0.9, 0.95)

    # set histograms
    h_cosCSFwd = ROOT.TH1D('h_'+baseName+'Fwd_CosCS', ';cos#theta_{CS};# of dimuon', 30, -1, 1)
    h_cosCSFwd.SetLineColor(ROOT.kBlue)
    h_cosCSFwd.SetLineWidth(2)
    h_cosCSFwd.SetMarkerColor(ROOT.kBlue)
    h_cosCSFwd.SetMarkerStyle(20)

    h_cosCSMid = ROOT.TH1D('h_'+baseName+'Mid_CosCS', '', 30, -1, 1)
    h_cosCSMid.SetLineColor(ROOT.kRed)
    h_cosCSMid.SetLineWidth(2)
    h_cosCSMid.SetMarkerColor(ROOT.kRed)
    h_cosCSMid.SetMarkerStyle(21)

    h_cosCSAllY = ROOT.TH1D('h_'+baseName+'AllY_CosCS', '', 30, -1, 1)
    h_cosCSAllY.SetLineColor(ROOT.kGreen + 2)
    h_cosCSAllY.SetLineWidth(2)
    h_cosCSAllY.SetMarkerColor(ROOT.kGreen + 2)
    h_cosCSAllY.SetMarkerStyle(22)

    # create PlotOverlaid instance
    plotCosCS = PlotOverlaid(canvas, leg, sampleInfo=sampleInfo, isWeight=isWeight)

    # add hists
    plotCosCS.add(tree=myTree, var='cosCS', cut=cutBaseFwd, hist=h_cosCSFwd, legEntry='1.6<|y|<2.4, 3<#it{p}_{T}<50', drawOpt='ep')
    plotCosCS.add(tree=myTree, var='cosCS', cut=cutBaseMid, hist=h_cosCSMid, legEntry='|y|<1.6, 6.5<#it{p}_{T}<50', drawOpt='ep')
    plotCosCS.add(tree=myTree, var='cosCS', cut=cutBaseAllY, hist=h_cosCSAllY, legEntry='|y|<2.4, 6.5<#it{p}_{T}<50', drawOpt='ep')

    # draw and save
    plotCosCS.drawRun3MC()
    plotCosCS.save(f'figs/{baseName}_overlaid/cosCS', formats=formats)

    # cleaning
    f.Close()