from ROOT import gSystem
gSystem.Load('Analysis.so')

from ROOT import EventLoop
eventLoop = EventLoop()

# set the attributes
eventLoop.treeName = 'hionia/myTree'
eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root')

# create algorithms
from ROOT import Algorithm
alg = Algorithm()

# create histograms
from ROOT import TH1D
alg.h_runNb_m = TH1D('h_runNb', 'run number', 150000, 350000, 500000)
alg.h_runNb_m.Sumw2()

# add the algorithm into the event loop
eventLoop.algorithms.push_back(alg)

# create TLorentzVectors
# from ROOT import TLorentzVector
# recoJpsi = TLorentzVector()
# recoMupl = TLorentzVector()
# recoMumi = TLorentzVector()

# intialize and execute the event loop
eventLoop.initialize()
eventLoop.execute()

# Saving and cleaning
from ROOT import TFile
f = TFile.Open('histograms.root', 'recreate')
f.cd()

alg.h_runNb_m.Write()
f.Close()