import os, ROOT
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

# Draw
plotOverlaid(fileName='Run2MC.root', baseName = 'MC', formats=['png'], sampleInfo='PbPb MC 2018 Prompt', isWeight=True)
plotOverlaid(fileName='Run2DataCent.root', baseName = 'DataCent', formats=['png'], sampleInfo='PbPb Data Cent.', isWeight=True)
plotOverlaid(fileName='Run2DataPeri.root', baseName = 'DataPeri', formats=['png'], sampleInfo='PbPb Data Peri', isWeight=True)