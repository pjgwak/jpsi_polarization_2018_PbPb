from ROOT import gSystem
gSystem.Load('Analysis.so')

from ROOT import EventLoop
eventLoop = EventLoop()

# set the attributes
eventLoop.treeName = 'hionia/myTree'
eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root')

# intialize and execute the event loop
eventLoop.initialize()
eventLoop.execute()