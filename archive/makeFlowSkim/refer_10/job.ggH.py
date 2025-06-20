from ROOT import gSystem, TH1
TH1.AddDirectory(False)
gSystem.Load("Analysis")
from Samples import *
from Algorithms import *
eventLoop = SampleGGH()
algs = []
algs += [ AlgDefault() ]
eventLoop.addAlgorithms( algs) 
algs += [ AlgSF() ]
eventLoop.addAlgorithms( algs) 
algs += [ AlgDF() ]
eventLoop.addAlgorithms( algs) 
eventLoop.execute()
eventLoop.save()
print("all ok")
