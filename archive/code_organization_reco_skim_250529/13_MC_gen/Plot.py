import ROOT
from ROOT import TCanvas, TFile, gPad, kWhite, kBlue, TLegend

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

class PlotOverlaid(object):
  """ overlay three plots on one canvas
  Fwd, Mid, AllY
  """

  def __init__(self, canvas, legend, sampleInfo='', isWeight=False):
    self.canvas = canvas
    self.legend = legend
    self.plots = [] # tree, var, cut, hist, legEntries, draw_opt
    self.isWeight = isWeight
    
    # set basic legend style
    self.legend.SetBorderSize(0)
    self.legend.SetFillColorAlpha(ROOT.kWhite, 0)
    self.legend.SetTextSize(0.035)

    self.legSample = TLegend(0.15, 0.83, 0.25, 0.95)
    self.legSample.SetBorderSize(0)
    self.legSample.SetFillColorAlpha(ROOT.kWhite, 0)
    self.legSample.SetTextSize(0.035)
    self.sampleInfo = sampleInfo

  def add(self, tree, var, cut, hist, legEntry, drawOpt='EP'):
    """
    tree: TTre instance, var: branch name
    cut: selection cut, hist: TH1 instance
    legEntry: text to draw as legends
    """
    self.plots.append({
      'tree': tree,
      'var': var,
      'cut': cut,
      'hist': hist,
      'legEntry': legEntry,
      'drawOpt': drawOpt
    })
  
  def addLegSample(self, leg):
    self.legSample = leg

  def draw(self):
    self.canvas.cd()

    # fill the histogram without drawing (GOFF option)
    hists = []
    for plot in self.plots:
      hist = plot['hist']
      drawCommand = f"{plot['var']}>>{hist.GetName()}"
      
      if self.isWeight:
        plot['tree'].Draw(drawCommand, f"{plot['cut']} * weight * TnPweight", 'GOFF')
      else:
        plot['tree'].Draw(drawCommand, plot['cut'], 'GOFF')
      hists.append(hist)

    # find Y max value
    maxY = 0
    for hist in hists:
      if hist.GetEntries() > 0: # check when a hist is not empty
        currentMax = hist.GetMaximum()
        if currentMax > maxY:
          maxY = currentMax
    
    # set Y max
    firstPlot = self.plots[0]
    firstHist = firstPlot['hist']
    firstHist.GetYaxis().SetRangeUser(0, maxY * 1.4 if maxY > 0 else 1)

    # draw first hist
    firstHist.Draw(firstPlot['drawOpt'])
    self.legend.AddEntry(firstHist, firstPlot['legEntry'], 'ep')
    self.legSample.AddEntry(self.sampleInfo, self.sampleInfo, '')
    firstHist.GetXaxis().CenterTitle(True)
    firstHist.GetYaxis().CenterTitle(True)

    # draw other hists
    for idx in range(1, len(self.plots)):
      plot = self.plots[idx]
      hist = plot['hist']
      drawOpt = plot['drawOpt']
      if 'SAME' not in drawOpt.upper():
        drawOpt += ' SAME'
      hist.Draw(drawOpt)
      self.legend.AddEntry(hist, plot['legEntry'], 'ep')

    # draw legend
    self.legend.Draw()
    self.legSample.Draw()

    # update canvas
    self.canvas.Update()
    
    # draw axes
    gPad.RedrawAxis()
  
  def drawRun3MC(self):
    self.canvas.cd()

    # fill the histogram without drawing (GOFF option)
    hists = []
    for plot in self.plots:
      hist = plot['hist']
      drawCommand = f"{plot['var']}>>{hist.GetName()}"
      
      if self.isWeight:
        plot['tree'].Draw(drawCommand, f"{plot['cut']} * weight", 'GOFF')
      else:
        plot['tree'].Draw(drawCommand, plot['cut'], 'GOFF')
      hists.append(hist)

    # find Y max value
    maxY = 0
    for hist in hists:
      if hist.GetEntries() > 0: # check when a hist is not empty
        currentMax = hist.GetMaximum()
        if currentMax > maxY:
          maxY = currentMax
    
    # set Y max
    firstPlot = self.plots[0]
    firstHist = firstPlot['hist']
    firstHist.GetYaxis().SetRangeUser(0, maxY * 1.4 if maxY > 0 else 1)

    # draw first hist
    firstHist.Draw(firstPlot['drawOpt'])
    self.legend.AddEntry(firstHist, firstPlot['legEntry'], 'ep')
    self.legSample.AddEntry(self.sampleInfo, self.sampleInfo, '')
    firstHist.GetXaxis().CenterTitle(True)
    firstHist.GetYaxis().CenterTitle(True)

    # draw other hists
    for idx in range(1, len(self.plots)):
      plot = self.plots[idx]
      hist = plot['hist']
      drawOpt = plot['drawOpt']
      if 'SAME' not in drawOpt.upper():
        drawOpt += ' SAME'
      hist.Draw(drawOpt)
      self.legend.AddEntry(hist, plot['legEntry'], 'ep')

    # draw legend
    self.legend.Draw()
    self.legSample.Draw()

    # update canvas
    self.canvas.Update()
    
    # draw axes
    gPad.RedrawAxis()
  
  def save(self, fPath='', formats=('png')):
    for fmt in formats:
      fullPath = f'{fPath}.{fmt}'
      self.canvas.SaveAs(fullPath)