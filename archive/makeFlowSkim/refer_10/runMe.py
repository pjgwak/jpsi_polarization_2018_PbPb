# load shared library
from ROOT import gSystem, TH1, kFALSE, TStopwatch
gSystem.Load('Analysis.so')

# avoid auto-registration of the THist in the root system.
TH1.AddDirectory(kFALSE)

# create instance of the class EventLoop
from Samples import *
eventLoops = []
eventLoops += [ SampleGGH() ]
eventLoops += [ SampleVBFH() ]

# create algorithm and add them into the event loops
from Algorithms import *
for eventLoop in eventLoops:
    algs = []
    algs += [AlgDefault()]
    algs += [AlgSF()]
    algs += [AlgDF()]
    eventLoop.addAlgorithms(algs)

# initialize and execute the event loop

# execute parallel jobs
timer = TStopwatch()
from Jobs import Jobs
jobs = Jobs(eventLoops)
jobs.execute()
timer.Stop()
print('The processing took {} s'.format(timer.RealTime()))