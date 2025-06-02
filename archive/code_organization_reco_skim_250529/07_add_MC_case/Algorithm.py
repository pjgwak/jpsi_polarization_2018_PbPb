from ROOT import Algorithm as AlgorithmCpp
from ROOT import TFile

class Algorithm(object):
  """ wrapper around the cpp Algorithm class
  """
  def __init__(self, name):
    self.name = name
    self.alg = AlgorithmCpp()
  
  def setSumw2(self):
    """ call Sumw2 for all the histograms
    """
    for attrName in dir(self.alg):
      attr = getattr(self.alg, attrName)
      if hasattr(attr, "Sumw2"):
        attr.Sumw2()
    
  def save(self, prefix):
    """ save the histograms into the output file
    """
    # save the everything into the ROOT file
    f = TFile.Open(f'histograms.{prefix}.{self.name}.root', 'recreate')
    f.cd()

    # printout
    print(f'Histograms saved into file: {f.GetName()}')

    # loop through all attributes of self.alg class
    # and save all classes that have 'Write' method
    for attrName in dir(self.alg):
      attr = getattr(self.alg, attrName)
      if hasattr(attr, "Write"):
        attr.Write()
    
    # close the file
    f.Close()