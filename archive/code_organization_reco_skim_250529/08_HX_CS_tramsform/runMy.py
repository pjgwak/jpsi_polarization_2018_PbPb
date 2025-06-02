from ROOT import gSystem
gSystem.Load('Analysis.so')

# create instance s of the class EventLoop
from Samples import *
eventLoops = []
eventLoops += [Data2023()]
eventLoops += [MC2023()]
# eventLoops += [MC2023GlbTrk()]

# add the algorithm into the event loop
from Algorithms import * 
for eventLoop in eventLoops:
  algs = []
  algs += [AlgDefault()]
  # algs += [AlgOffSoftMuonCut()]
  eventLoop.addAlgorithms(algs)

# execute event loop
for eventLoop in eventLoops:
  eventLoop.execute()
  eventLoop.save()

# ===== old - saving TFile ===== #
# from ROOT import TFile
# f = TFile.Open('histograms.root', 'recreate')
# f.cd()

# alg.h_runNb_m.Write()
# f.Close()