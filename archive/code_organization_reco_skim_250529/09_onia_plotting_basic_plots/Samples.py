from EventLoop import EventLoop

# ---------------------------------
class Data2023(EventLoop):
  """ event loop over 2023 PbPb Data
  Muon: GlbTrk
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='Data2023', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root')
# --------------------------------
class MC2023(EventLoop):
  """ event loop over 2023 PbPb MC
  Muon: All
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.isMC = True
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_MC/OniaTree_MC1_Run3_PbPb_no_track_250519/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/OniaTree_MC1_Run3_PbPb_no_track_250519/250519_112255/0000/OniaTree_MC1_Run3_PbPb2023_250520.root')
# ---------------------------------
class MC2023GlbTrk(EventLoop): # Not used
  """ event loop over the vector boson fusion sample
  Muon: GlbTrk.
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023GlbTrk', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.isMC = True
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/work/pjgwak/oniatree_5p36/2023_MC/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/OniaTree_MC1_Run3_PbPb_no_track_GlbTrk_250531/250531_092232/Oniatree_MC_miniAOD_GlbTrk_250601.root')
# ---------------------------------