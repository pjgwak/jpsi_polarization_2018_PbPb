from ROOT import Algorithm as AlgorithmCpp
from ROOT import TFile

class Algorithm(object):
  """ Wrapper around the cpp Algorithm class
  """
  def __init__(self, name):
    self.name = name

    self.alg = AlgorithmCpp()
  
  def setSumw2(self):
    """ Calls Sumw2 for all the histograms
    """
    for attrName in dir(self.alg):
      attr = getattr(self.alg, attrName)
      if hasattr(attr, 'Sumw2'):
        attr.Sumw2()
  
  def save(self, prefix):
    """ Save the histograms into the output file
    """
    # save everything into the ROOT file
    f = TFile.Open('histograms.{}.{}.root'.format(prefix,self.name), 'recreate')
    f.cd()

    # printout
    print('Histograms saved into file: {}'.format(f.GetName()))

    # loop through all attributes of self.alg class
    # and save all classes that have 'Write' method
    for attrName in dir(self.alg):
      attr = getattr(self.alg, attrName)
      if hasattr(attr, 'Write'):
        attr.Write()
    
    # close the file
    f.Close()