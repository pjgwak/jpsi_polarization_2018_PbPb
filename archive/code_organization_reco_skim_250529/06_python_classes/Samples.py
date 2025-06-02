from EventLoop import EventLoop

# ---------------------------------
class Data2023(EventLoop):
  """ event loop over the gluon-gluon fusion smaple
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, 'Data2023')

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root')
# --------------------------------
class MC2023(EventLoop):
  """ event loop over the vector boson fusion sample
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, 'MC2023')


    # add the ggH samples into the event loop
    self.eventLoop.isMC = True
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_MC/OniaTree_MC1_Run3_PbPb_no_track_250519/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/OniaTree_MC1_Run3_PbPb_no_track_250519/250519_112255/0000/OniaTree_MC1_Run3_PbPb2023_250520.root')
# ---------------------------------
# ---------------------------------