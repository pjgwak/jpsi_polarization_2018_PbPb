from Plot import Plot
from ROOT import TCanvas, TLegend, kRed, kGreen, kBlack, kMagenta, kYellow

# set ATLAS plot style
from AtlasStyle import setStyle
setStyle()

# draw plots
histNames = ['ditau_m', 'ditau_visM']
flavorNames = ['SF', 'DF']
canvases = []
plots = []
legends = []

for histName in histNames:
  # load individual plots
  plots_ggH = []
  plots_VBFH = []
  for flavorName in flavorNames:
    plots_ggH += [Plot(f'histograms.ggH.Alg{flavorName}.root', histName)]
    plots_VBFH += [Plot(f'histograms.VBFH.Alg{flavorName}.root', histName)]
  
  # add plots from two samples into a single plot
  plot_ggH = Plot(plots_ggH)
  plot_VBFH = Plot(plots_VBFH)

  # set style
  plot_ggH.setStyleSolid(kMagenta)
  plot_VBFH.setStyleSolid(kYellow)

  # create a THStack
  plot = Plot([plot_ggH, plot_VBFH], True)

  # create a canvas
  c = TCanvas(histName, histName, 800, 600)

  # draw
  plot.draw()

  # create the same plot, but this time merging the histograms
  plotErr = Plot([plot_ggH, plot_VBFH])
  plotErr.setStyleErrorbar(kBlack)
  plotErr.draw('same')

  # draw legend
  legend = TLegend(0.6, 0.9, 0.9, 0.65)
  legend.AddEntry(plot_ggH.hist, 'gluon-gluon fusion', 'f')
  legend.AddEntry(plot_VBFH.hist, 'vector-boson fusion', 'f')
  legend.AddEntry(plotErr.hist, 'Stat. uncertainty', 'f')
  legend.Draw('same')

  # save the canvas
  c.Print(f'{histName}.pdf')

  # save the instances so they are not deleted by the garbage collector
  canvases += [c]
  plots += [plot, plotErr]
  legends += [legend]