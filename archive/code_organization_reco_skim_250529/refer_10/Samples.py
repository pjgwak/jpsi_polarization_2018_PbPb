from EventLoop import EventLoop

#---------------------------------
class SampleGGH(EventLoop):
  """ event loop over the gluon-gluon fusion sample
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, 'ggH')

    # add the ggH samples into the event loop
    self.eventLoop.inputFiles.push_back('../../ggH.root')

    # set sample meta-data
    self.isMC = True
    self.xSection = 16.3789 # pb
    self.luminosity = 36074.16 # pb^-1

#---------------------------------
class SampleVBFH(EventLoop):
  """ event loop over the vector boson fusion sample
  """
  def __init__(self):
    # call the inherited constructor
    EventLoop.__init__(self, 'VBFH')

    # add the ggH samples into the event loop
    self.eventLoop.inputFiles.push_back('../../vbfH.root')

    # set sample meta-data
    self.isMC = True
    self.xSection = 1.3741 # pb
    self.luminosity = 36074.16 # pb^-1