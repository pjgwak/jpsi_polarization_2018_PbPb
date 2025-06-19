import os, ROOT
from plotMeDataOverlay import *
import plotMeMCRun3Overlay
from ROOT import TFile, TCanvas, TLegend, kGreen, kBlack, TCanvas, TH1

# set batch mode
ROOT.gROOT.SetBatch(True)

# set CMS plot style
rootlogon_path = "../../../rootlogon.C"
if os.path.exists(rootlogon_path):
    print('Use rootlogon.C')
    ROOT.gROOT.ProcessLine(f".X {rootlogon_path}")
else:
    print(f"Not exist: {rootlogon_path}")

# Draw
plotOverlaid(fileName='../../code_organization_reco_skim_250529/13_MC_gen/OniaFlowSkim.Data2018All.AlgRun2Cent.root', baseName = 'Data2018', formats=['png'], sampleInfo='PbPb Data 2018', isWeight=False)
plotOverlaid(fileName='../../code_organization_reco_skim_250529/13_MC_gen/OniaFlowSkim.MC2018Pr.AlgRun2MCReco.root', baseName = 'MC2018', formats=['png'], sampleInfo='PbPb MC 2018', isWeight=True)
plotOverlaid(fileName='../../code_organization_reco_skim_250529/13_MC_gen/OniaFlowSkim.Data2023.AlgDefault.root', baseName = 'Data2023', formats=['png'], sampleInfo='PbPb Data 2023', isWeight=False)
plotMeMCRun3Overlay.plotOverlaid(fileName='../../code_organization_reco_skim_250529/13_MC_gen/OniaFlowSkim.MC2023.AlgRun3MCReco.root', baseName = 'MC2023', formats=['png'], sampleInfo='PbPb MC 2023', isWeight=True)


#OniaFlowSkim.MC2018Pr.AlgRun2MCReco