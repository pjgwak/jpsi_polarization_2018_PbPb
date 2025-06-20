from ROOT import gSystem
gSystem.Load('Analysis.so')

# create instance s of the class EventLoop
from Samples import *
eventLoops = []
# eventLoops += [Data2023()]
# eventLoops += [Data2018Cent()]
# eventLoops += [Data2018Peri()]
# eventLoops += [MC2023()]
eventLoops += [MC2018Pr()] # Jpsi Prompt
# eventLoops += [MC2023GlbTrk()]


# add the algorithm into the event loop
from Algorithms import * 
for eventLoop in eventLoops:
  # common algorithms
  algs = []
  algs += [AlgDefault()]
  # algs += [AlgNoCut()]
  eventLoop.addAlgorithms(algs)

  # special algorithms
  # Run2 Data
  if isinstance(eventLoop, Data2018Cent):
    eventLoop.addAlgorithm(AlgRun2Cent())
  if isinstance(eventLoop, Data2018Peri):
    eventLoop.addAlgorithm(AlgRun2Peri())
  
  # Run2 MC
  if isinstance(eventLoop, MC2018Pr):
    eventLoop.addAlgorithm(AlgRun2Peri())


# execute event loop
for eventLoop in eventLoops:
  eventLoop.execute()
  eventLoop.save()