from Algorithm import Algorithm
from ROOT import TH1D, TTree

# ---------------------------------
class AlgDefault(Algorithm):
  """ Default set of oniaFlowSkim
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
# ---------------------------------