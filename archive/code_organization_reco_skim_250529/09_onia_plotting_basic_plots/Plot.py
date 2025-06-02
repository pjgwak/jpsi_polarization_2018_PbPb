from ROOT import TCanvas, TFile, gPad, kWhite, kBlue

class Plot1D(object):
  """ 1D plot helper class
  """

  def __init__(self, tree=None, var='', cut='', hist=None, histTitles='', leg=None, legEntry=['','','',''], canvas=None):
    # self.hist = None # an instance inheriting TH1 class
    self.canvas = canvas # TCanvas
    self.tree = tree # TTree
    self.cut = cut # selection cut
    self.leg = leg # TLegend
    self.legEntry = legEntry # Centent of legend
    self.var = var # branch name
    self.drawOpt = ''

    # used when one doesn't provide a histogram.
    self.hist = hist
    self.histTitles = histTitles

  def draw(self):
    self.canvas.cd()
    self.tree.Draw(f'{self.var}>>{self.hist.GetName()}', self.cut, 'EP')
    
    # get htemp - TTree's Draw() makes a htemp
    # self.hist = gPad.GetPrimitive("htemp")
    self.setHistStyle()

    if self.leg:
      self.drawLegend()
    
    # udpate a pad
    gPad.RedrawAxis()
  
  def setHistStyle(self):
    self.hist.SetTitle(self.histTitles)
    self.hist.SetMarkerColor(kBlue)
    self.hist.SetLineColor(kBlue)
    self.hist.SetMarkerStyle(21)
    self.hist.GetXaxis().CenterTitle(True)
    self.hist.GetYaxis().CenterTitle(True)
    self.hist.GetYaxis().SetRangeUser(0, self.hist.GetMaximum() * 1.4)
  
  def save(self, fPath='', isPng=True, isPdf=False):
    if isPng:
      self.canvas.SaveAs(f'{fPath}{self.canvas.GetName()}.png')
    if isPdf:
      self.canvas.SaveAs(f'{fPath}{self.canvas.GetName()}.pdf')
  
  def drawLegend(self):
    self.leg.SetBorderSize(0) # no boarder
    self.leg.SetFillColorAlpha(kWhite, 0) # transperency
    self.leg.SetTextSize(0.03)
    self.leg.AddEntry(self.hist, self.legEntry[0], 'ep')
  
    # Additional describtions
    if self.legEntry[1]:
      self.leg.AddEntry(self.legEntry[1], self.legEntry[1], '') 
    if self.legEntry[2]:
      self.leg.AddEntry(self.legEntry[2], self.legEntry[2], '')
    if self.legEntry[3]:
      self.leg.AddEntry(self.legEntry[3], self.legEntry[3], '')

    self.leg.Draw()