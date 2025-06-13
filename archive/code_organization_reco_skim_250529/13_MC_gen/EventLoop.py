from ROOT import EventLoop as EventLoopCpp

class EventLoop(object):
  """ wrapper around the cpp class event loop
  """
  def __init__(self, name, nEvt=-1, isMC=False, isGenOnly=False,
    isEP=False, isRun2=False):
    # set the event loop name
    self.name = name

    # now we can create instance of the cpp class EventLoop
    self.eventLoop = EventLoopCpp()

    # set the tree name attribute
    self.eventLoop.treeName = 'hionia/myTree'

    # set a number of events
    self.eventLoop.nEvent = nEvt

    # set the flags
    self.eventLoop.isMC = isMC
    self.eventLoop.isGenOnly = isGenOnly
    self.eventLoop.isEP = isEP
    self.eventLoop.useRun2 = isRun2

    # keep the list of python Algorithm class instances
    self.algs = []

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
    self.eventLoop.initialize()
    self.eventLoop.execute()

  def executeGen(self):
    """ initialize and execute the event loop
    """
    self.eventLoop.initialize()
    self.eventLoop.executeGen()

  def save(self):
    """ save histograms from all algorithms
    """
    for alg in self.algs:
      alg.save(self.name)