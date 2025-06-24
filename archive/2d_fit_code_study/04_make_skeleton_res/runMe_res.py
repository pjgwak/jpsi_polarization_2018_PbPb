from JpsiFitter import JpsiFitter
import ROOT
ROOT.gROOT.SetBatch(True)
# import os

mass_result_file = 'roots/mass.root'
splot_result_file = 'roots/splot.root'
output_plot_path = 'figs/res.png'

analysis = JpsiFitter(config_module_name='mc_mass')
analysis.fitter.loadMassResult(mass_result_file) 
analysis.fitter.loadSPlotResult(splot_result_file)
analysis.model.buildCtauResModel(2)
analysis.fitter.performResFit()
analysis.fitter.drawResPlot(output_plot_path)

print('===== Finish runMe_res.py =====')