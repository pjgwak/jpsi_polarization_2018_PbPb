from Algorithm import Algorithm
from ROOT import TH1D, TTree

# ---------------------------------
class AlgDefault(Algorithm):
  """ Default set of oniaFlowSkim
  Apply very coarse selection cuts
  It's good enough for Run3 "Data" samples when consdiering minBias study.
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, "AlgDefault")

    # turn on cuts
    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    # create histograms
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', 'Noimial')
# ---------------------------------
class AlgNoCut(Algorithm):
  """ Turn off all Reco cuts
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgNoCut')

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', 'No reconstruction cut applied')
# ---------------------------------
class AlgRun3MCReco(Algorithm):
  """ Run3 MC Reco
  No triggers, filters
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun3MCReco')

    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    # self.alg.cut_HLTriggerPbPbJpsi2018 = True
    # self.alg.cut_recoQQTrigger = True
    # self.alg.cut_tnp = True

    self.alg.cut_whichGen = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun3MCRecoTriggers(Algorithm):
  """ Run3 MC Reco
  Aplly triggers, filters
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun3MCRecoTriggers')

    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    self.alg.cut_HLTriggerPbPbJpsi2018 = True
    self.alg.cut_recoQQTrigger = True
    # self.alg.cut_tnp = True # it makes TnP errors

    self.alg.cut_whichGen = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
class AlgRun2Cent(Algorithm):
  """ Run2 Data central
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2Cent')

    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    self.alg.cut_HLTriggerPbPbJpsi2018 = True
    self.alg.cut_recoQQTrigger = True
    self.alg.cut_runNb327123 = True
    # self.alg.cut_L2L3FilterPbPbJpsi2018 = True # only MC for TnP

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun2Peri(Algorithm):
  """ Run2 Data Periheral
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2Peri')

    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    self.alg.cut_HLTriggerPbPbJpsi2018 = True
    self.alg.cut_recoQQTrigger = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun2MCReco(Algorithm):
  """ Run2 MC Reco
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2MCReco')

    self.alg.cut_centrality0_180 = True
    self.alg.cut_jpsiMass = True # 2.6 - 3.5
    self.alg.cut_softMuons = True
    self.alg.cut_vtxProbability = True
    self.alg.cut_oppositeSign = True
    self.alg.cut_singleMuonAcc = True

    self.alg.cut_HLTriggerPbPbJpsi2018 = True
    self.alg.cut_recoQQTrigger = True

    self.alg.cut_tnp = True
    self.alg.cut_whichGen = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun2MCRecoGen(Algorithm):
  """ Run2 MC Reco
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2MCRecoGen', isGen=True)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    self.alg.cut_jpsiRapidity = True # 2.6 - 3.5
    self.alg.cut_muChargeGen = True
    self.alg.cut_singleMuonAcc = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun2MCGenNum(Algorithm):
  """ Run2 Gen-Only
  Acceptance numerator
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2MCGenNum', isGen=True, isGenOnly=True)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    
    self.alg.cut_jpsiRapidity = True
    # self.alg.cut_muChargeGen = True
    self.alg.cut_singleMuonAcc = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun2MCGenDen(Algorithm):
  """ Run2 Gen-Only
  Acceptance denominator
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun2MCGenDen', isGen=True, isGenOnly=True)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    self.alg.cut_muChargeGen = False
    self.alg.cut_singleMuonAcc = False

    self.alg.cut_jpsiRapidity = True # |y| < 2.4

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun3MCRecoGen(Algorithm):
  """ Run2 MC Reco
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun3MCRecoGen', isGen=True, isGenOnly=False)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    self.alg.cut_jpsiRapidity = True # 2.6 - 3.5
    self.alg.cut_muChargeGen = True
    self.alg.cut_singleMuonAcc = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun3MCGenNum(Algorithm):
  """ Run3 Gen-Only
  Acceptance numerator
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun3MCGenNum', isGen=True, isGenOnly=True)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    
    self.alg.cut_jpsiRapidity = True
    # self.alg.cut_muChargeGen = True
    self.alg.cut_singleMuonAcc = True

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------
class AlgRun3MCGenDen(Algorithm):
  """ Run3 Gen-Only
  Acceptance denominator
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgRun3MCGenDen', isGen=True, isGenOnly=True)

    self.alg.cut_centrality0_180 = False
    self.alg.cut_jpsiMass = False # 2.6 - 3.5
    self.alg.cut_muChargeGen = False
    self.alg.cut_singleMuonAcc = False

    self.alg.cut_jpsiRapidity = True # |y| < 2.4

    # create histograms: fine binning
    # self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree', '')
# ---------------------------------