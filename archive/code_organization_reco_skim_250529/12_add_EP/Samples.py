from EventLoop import EventLoop

# ---------------------------------
class Data2023(EventLoop):
  """ event loop over 2023 PbPb Data
  Muon: GlbTrk
  Currently, it's only possible to calculate EP angles for Run3 data samples.
  99.8 % of the miniAOD format was succesfully processed to make an Oniatree
  Run2 - EP angles are not recentered and not flattened
  Run3 MC - Oniatree was not produced due to lack of my EOS folder.
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='Data2023', nEvt=-1)

    # push samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_track_vector_250530/Oniatree_MC0_PbPb_2023_EP_angles_250610.root') # inlcudes EP angles  
    # self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root') # no EP branches

    # set sample flags
    self.eventLoop.isEP = True
# --------------------------------
class Data2018Cent(EventLoop):
  """ event loop over 2018 PbPb Data
  Double muon central
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='Data2018Cent', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.useRun2 = True
    self.eventLoop.treeName = 'myTree'
    
    for i in range(1, 6):
      self.eventLoop.inputFiles.push_back(f'/disk1/Oniatree/polarization/run2_oniatree_PbPb_data_AOD/DM/ReReco_Oniatree_addvn_part{i}.root')
# --------------------------------
class Data2018Peri(EventLoop):
  """ event loop over 2018 PbPb Data
  Double muon peripheral
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='Data2018Peri', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.useRun2 = True
    self.eventLoop.treeName = 'myTree'
    
    for i in range(1, 4):
      self.eventLoop.inputFiles.push_back(f'/disk1/Oniatree/polarization/run2_oniatree_PbPb_data_AOD/DMPeri/ReReco_Oniatree_addvn_part{i}.root')
# --------------------------------
class MC2018Pr(EventLoop):
  """ event loop over 2018 MC - Jpsi Prompt
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2018Pr', nEvt=-1)

    # add the ggH samples into the event loop
    self.eventLoop.useRun2 = True
    self.eventLoop.isMC = True
    self.eventLoop.treeName = 'myTree'
    
    for i in range(1, 3):
      self.eventLoop.inputFiles.push_back(f'/disk1/Oniatree/Jpsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part{i}.root')
# --------------------------------
class MC2023(EventLoop):
  """ event loop over 2023 PbPb MC
  Muon: All
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023', nEvt=1000000)

    # add the ggH samples into the event loop
    self.eventLoop.isMC = True
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/OniaTree_MC1_Run3_PbPb_no_track_250519/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/OniaTree_MC1_Run3_PbPb_no_track_250519/250519_112255/0000/OniaTree_MC1_Run3_PbPb2023_250520.root')
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