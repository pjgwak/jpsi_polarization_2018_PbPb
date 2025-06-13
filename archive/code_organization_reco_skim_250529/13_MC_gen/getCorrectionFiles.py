from ROOT import gSystem
gSystem.Load('Analysis.so')

# create instance s of the class EventLoop
from Samples import *
eventLoops = []
# eventLoops += [MC2023()]
eventLoops += [MC2018Pr()] # Jpsi Prompt
eventLoops += [MC2018GenOnly()]
# eventLoops += [MC2023GlbTrk()]


# add the algorithm into the event loop
from Algorithms import * 
for eventLoop in eventLoops:
  # common algorithms
  # algs = []
  # algs += [AlgDefault()]
  # # algs += [AlgNoCut()]
  # eventLoop.addAlgorithms(algs)

  # # special algorithms
  # # Run3 Data
  # if isinstance(eventLoop, Data2023):
  #   eventLoop.addAlgorithm(AlgRun2Peri()) # Trigger test

  # # Run3 MC
  # if isinstance(eventLoop, MC2023):
  #   eventLoop.addAlgorithm(AlgRun3MCReco())
  #   eventLoop.addAlgorithm(AlgRun3MCRecoTriggers())

  # # Run2 Data
  # if isinstance(eventLoop, Data2018Cent):
  #   eventLoop.addAlgorithm(AlgRun2Cent())
  # if isinstance(eventLoop, Data2018Peri):
  #   eventLoop.addAlgorithm(AlgRun2Peri())
  
  # Run2 MC
  if isinstance(eventLoop, MC2018Pr):
    eventLoop.addAlgorithm(AlgRun2MCRecoGen())

  if isinstance(eventLoop, MC2018GenOnly):
    eventLoop.addAlgorithm(AlgRun2MCRecoGen())
    eventLoop.addAlgorithm(AlgRun2MCRecoGen())


# execute event loop
for eventLoop in eventLoops:
  eventLoop.executeGen()
  eventLoop.save()

# Merge run2 Data samples
# hadd OniaFlowSkim.Data2018All.AlgRun2.root OniaFlowSkim.Data2018Cent.AlgRun2Cent.root OniaFlowSkim.Data2018Peri.AlgRun2Peri.root