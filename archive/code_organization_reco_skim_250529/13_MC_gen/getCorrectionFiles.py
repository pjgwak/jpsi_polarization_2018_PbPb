from ROOT import gSystem
gSystem.Load('Analysis.so')

# create instance s of the class EventLoop
from Samples import *
eventLoops = []
# run3
eventLoops += [MC2023()]
eventLoops += [MC2023GenOnly()]
eventLoops += [MC2023GenOnlyShower()]

# run2
eventLoops += [MC2018Pr()] # Jpsi Prompt
eventLoops += [MC2018GenOnly()]




# add the algorithm into the event loop
from Algorithms import * 
for eventLoop in eventLoops:
  # run3 algorithms
  if isinstance(eventLoop, MC2023):
    eventLoop.addAlgorithm(AlgRun3MCRecoGen())

  if isinstance(eventLoop, MC2023GenOnly):
    eventLoop.addAlgorithm(AlgRun3MCGenNum())
    eventLoop.addAlgorithm(AlgRun3MCGenDen())

  if isinstance(eventLoop, MC2023GenOnlyShower):
    eventLoop.addAlgorithm(AlgRun3MCGenNum())
    eventLoop.addAlgorithm(AlgRun3MCGenDen())


  # Run2 MC
  if isinstance(eventLoop, MC2018Pr):
    eventLoop.addAlgorithm(AlgRun2MCRecoGen())

  if isinstance(eventLoop, MC2018GenOnly):
    eventLoop.addAlgorithm(AlgRun2MCGenNum())
    eventLoop.addAlgorithm(AlgRun2MCGenDen())


# execute event loop
for eventLoop in eventLoops:
  eventLoop.executeGen()
  eventLoop.save()