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
    EventLoop.__init__(self, name='Data2023', nEvt=-1, isEP=True)

    # push samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_PbPb/Oniatree_MC0_PbPb_2023_EP_angles.root')
# --------------------------------
class Data2018Cent(EventLoop):
  """ event loop over 2018 PbPb Data
  Double muon central
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='Data2018Cent', nEvt=-1, isRun2=True)

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
    EventLoop.__init__(self, name='Data2018Peri', nEvt=-1, isRun2=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'myTree'
    
    for i in range(1, 4):
      self.eventLoop.inputFiles.push_back(f'/disk1/Oniatree/polarization/run2_oniatree_PbPb_data_AOD/DMPeri/ReReco_Oniatree_addvn_part{i}.root')
# --------------------------------
class MC2018Pr(EventLoop):
  """ event loop over 2018 MC - Jpsi Prompt
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2018Pr', nEvt=-1, isMC=True, isRun2=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'myTree'
    
    for i in range(1, 3):
      self.eventLoop.inputFiles.push_back(f'/disk1/Oniatree/Jpsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part{i}.root')
      # --------------------------------
class MC2018GenOnly(EventLoop):
  """ event loop over 2018 MC - Gen only
  It's used for acceptance study
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2018GenOnly', nEvt=-1, isMC=True, isRun2=True, isGenOnly=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/Jpsi/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root')
# --------------------------------
class MC2023(EventLoop):
  """ event loop over 2023 PbPb MC
  Muon: All
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023', nEvt=-1, isMC=True, isEP=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC1_miniAOD_PbPb_2023.root')
# ---------------------------------
# --------------------------------
class MC2023GenOnly(EventLoop):
  """ I think it's not an official MC, isn't it?
  As usually, pythia was used to proudce this sample
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023GenOnly', nEvt=-1, isMC=True, isGenOnly=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_GenOnly/Oniatree_PromptJpsi_5p36TeV_GenOnly.root')
# ---------------------------------
# --------------------------------
class MC2023GenOnlyShower(EventLoop):
  """ I think it's not an official MC, isn't it?
  New functionality, shower mode, of pythia was used.
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, name='MC2023GenOnlyShower', nEvt=-1, isMC=True, isGenOnly=True)

    # add the ggH samples into the event loop
    self.eventLoop.treeName = 'hionia/myTree'
    self.eventLoop.inputFiles.push_back('/disk1/Oniatree/polarization/oniatree_5p36/2023_GenOnly/Oniatree_ShowerJpsi_pTHatMin10_5p36TeV_GenOnly.root')
# ---------------------------------