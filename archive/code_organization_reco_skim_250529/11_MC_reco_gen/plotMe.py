import os, ROOT
from Plot import Plot1D, PlotOverlaid
from plotMeDataOverlay import *
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


# ===== Data ===== #
# DataRun2 - Cent + Peri. Applied full reco cuts
# plotOverlaid(fileName='OniaFlowSkim.Data2018All.AlgRun2Cent.root', baseName = 'DataRun2', formats=['png'])

# DataRun3 - No triggers
# plotOverlaid(fileName='OniaFlowSkim.Data2023.AlgDefault.root', baseName = 'DataRun3', formats=['png'])


# ===== MC ===== #
# DataRun2 - Cent + Peri. Applied full reco cuts
plotOverlaid(fileName='OniaFlowSkim.MC2018Pr.AlgRun2Peri.root', baseName = 'MCRun2PbPbPrompt', formats=['png'])



print('\n===== Finished plotMe.py =====\n')