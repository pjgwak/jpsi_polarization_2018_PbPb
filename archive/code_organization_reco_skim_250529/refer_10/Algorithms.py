from Algorithm import Algorithm
from ROOT import TH1D

# ----------------------------------
class AlgDefault(Algorithm):
  """ Default set of histograms
  """
  def __init__(self, name='AlgDefault'):
    # call inherited constructor
    Algorithm.__init__(self, name)

    # create histograms
    self.alg.h_ditau_m = TH1D('ditau_m', ';#it{m}_{#tau#tau} [GeV];Events', 30, 0, 300)
    self.alg.h_ditau_visM = TH1D('ditau_visM', ';#it{m}_{ll} [GeV];Events', 30, 0, 150)

# ----------------------------------
class AlgSF(AlgDefault):
  """ Select only same-flavor leptons
  """
  def __init__(self, name='AlgSF'):
    # call inherited constructor
    AlgDefault.__init__(self, name)

    # select same-flavor
    self.alg.cut_selectSF = True
    
    # invariant mass cut < 75 GeV (for visible tau)
    self.alg.cut_maxVisMass.set(75.)

# ----------------------------------
class AlgDF(AlgDefault):
  """ Select only different-flavor leptons
  """
  def __init__(self, name='AlgDF'):
    # call inherited constructor
    AlgDefault.__init__(self, name)

    # select differnt-flavor
    self.alg.cut_selectDF = True