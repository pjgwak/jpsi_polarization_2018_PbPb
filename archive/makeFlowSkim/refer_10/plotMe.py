from Plot import Plot
from ROOT import TCanvas, TLegend, kRed, kGreen, kBlack

# set ATLAS plot style
from AtlasStyle import setStyle
setStyle()

# draw plots
histNames = ['ditau_m', 'ditau_visM']
sampleNames = ['ggH', 'VBFH']
canvases = []
plots = []
legends = []

for histName in histNames:
  # load individual plots
  plots_sf = []
  plots_df = []
  for sampleName in sampleNames:
    plots_sf += [Plot(f'histograms.{sampleName}.AlgSF.root', histName)]
    plots_df += [Plot(f'histograms.{sampleName}.AlgDF.root', histName)]
  
  # add plots from two samples into a single plot
  plot_sf = Plot(plots_sf)
  plot_df = Plot(plots_df)

  # set style
  plot_sf.setStyleSolid(kRed)
  plot_df.setStyleSolid(kGreen)

  # create a THStack
  plot = Plot([plot_sf, plot_df], True)

  # create a canvas
  c = TCanvas(histName, histName, 800, 600)

  # draw
  plot.draw()

  # create the same plot, but this time merging the histograms
  plotErr = Plot([plot_sf, plot_df])
  plotErr.setStyleErrorbar(kBlack)
  plotErr.draw('same')

  # draw legend
  legend = TLegend(0.6, 0.9, 0.9, 0.65)
  legend.AddEntry(plot_sf.hist, 'SF leptons', 'f')
  legend.AddEntry(plot_df.hist, 'DF leptons', 'f')
  legend.AddEntry(plotErr.hist, 'Stat. uncertainty', 'f')
  legend.Draw('same')

  # save the canvas
  c.Print(f'{histName}.pdf')

  # save the instances so they are not deleted by the garbage collector
  canvases += [c]
  plots += [plot, plotErr]
  legends += [legend]