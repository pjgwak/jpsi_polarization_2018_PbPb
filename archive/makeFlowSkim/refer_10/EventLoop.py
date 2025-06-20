# my
from ROOT import EventLoop as EventLoopCpp
from ROOT import TFile

class EventLoop(object):
  """ Wrapper around the cpp class event loop
  """
  def __init__(self, name):
    # set the event loop name
    self.name = name

    # now we can create instance of the cpp class event loop
    self.eventLoop = EventLoopCpp()

    # set the tree name attribute
    self.eventLoop.treeName = 'NOMINAL'

    # keep the list of python Algorithm class instances
    self.algs = [] # for loop 이용하려고 저장.

    # sample meta-data
    self.isMC = False
    self.xSection = 1.
    self.luminosity = 1.
    self.sumOfWeightsGenerated = 0

  def addAlgorithm(self, algorithm):
    """ add an algorithm into this event loop
    """
    algorithm.setSumw2()
    self.algs += [algorithm]
    self.eventLoop.algorithms.push_back(algorithm.alg)
  
  def addAlgorithms(self, algorithms):
    """ add multiple algorithms into this event loop
    """
    for alg in algorithms:
      self.addAlgorithm(alg)
  
  def execute(self):
    """ initialize and execute the event loop
    """
    # load sumOfWeightsGenerated and initialize algorithms
    if self.isMC:
      self.loadSumOfWeightsGenerated()
      self.setSampleWeight()
    
    self.eventLoop.initialize()
    self.eventLoop.execute()
  
  def save(self):
    """ save histograms from all algorithms
    """
    for alg in self.algs:
      alg.save(self.name)

  def loadSumOfWeightsGenerated(self):
    """ Open all input files and adds up the sum of weights for all
    generated events. Only for MC.
    """
    if self.isMC:
      self.sumOfWeightsGenerated = 0
      for inputFile in self.eventLoop.inputFiles:
        # open each file and load the n_metadata histogram
        f = TFile.Open(str(inputFile))
        h_metadata = f.Get('h_metadata')
        # the generated number of events is stored in the 8th bin in the Higgs samples
        self.sumOfWeightsGenerated += h_metadata.GetBinContent(8)

  def setSampleWeight(self):
    """ Calculate the sample weight and pass to the algorithms. Only for MC
    """
    if self.isMC:
      sampleWeight = self.xSection * self.luminosity / self.sumOfWeightsGenerated
      for alg in self.algs:
        alg.alg.p_isMC = True
        alg.alg.p_sampleWeight = sampleWeight