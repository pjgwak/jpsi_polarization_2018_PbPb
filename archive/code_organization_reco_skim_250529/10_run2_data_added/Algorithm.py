from ROOT import Algorithm as AlgorithmCpp
from ROOT import TFile, TH1

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
      if attr and isinstance(attr, TH1):
          attr.Sumw2()
    
  def save(self, prefix):
    """ save the histograms into the output file
    """
    # save the everything into the ROOT file
    # FlowSkim means Skimmed oniatree for flow measuremnet which was done a few years ago.
    # I just keep this name convention. (Be careful: Skim could mean process from AOD to miniAOD )
    f = TFile.Open(f'OniaFlowSkim.{prefix}.{self.name}.root', 'recreate')
    f.cd()

    # printout
    print(f'Histograms saved into file: {f.GetName()}')

    # loop through all attributes of self.alg class
    # and save all classes that have 'Write' method
    for attrName in dir(self.alg):
      attr = getattr(self.alg, attrName)
      if attr and hasattr(attr, "Write"):
        attr.Write()

    # close the file
    f.Close()