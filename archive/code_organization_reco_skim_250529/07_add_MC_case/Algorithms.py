from Algorithm import Algorithm
from ROOT import TH1D, TTree

# ---------------------------------
class AlgDefault(Algorithm):
  """ Default set of histograms
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, "AlgDefault")

    # create histograms
    self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree2', '')
  
# ---------------------------------
class AlgFineBinning(Algorithm):
  """ set of histograms with fine binning
  """
  def __init__(self):
    # call inherited constructor
    Algorithm.__init__(self, 'AlgFineBinning')

    # create histograms: fine binning
    self.alg.h_runNb_m = TH1D('h_runNb', 'run number', 40000, 360000, 400000)
    self.alg.m_myTree = TTree('myTree2', '')
# ---------------------------------
# ---------------------------------