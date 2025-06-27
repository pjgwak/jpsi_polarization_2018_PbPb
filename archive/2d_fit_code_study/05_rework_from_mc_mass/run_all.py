import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('./Analysis.so')

from ROOT import McMassFitter, MassFitter, SPlotter


# set kinematics
ptLow = 3.0
ptHigh = 6.5
yLow = 1.6
yHigh = 2.4
cLow = 60
cHigh = 180
cos_low = 0.5
cos_high = 0.6


print('===== Start run_fitter.py =====')
print('===== Start MC Fitter =====')
# mc_fitter = McMassFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
# mc_fitter.run()
print('===== Finish MC Fitter =====')


print('===== Start Mass Fitter =====')
# mass_fitter = MassFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
# mass_fitter.run()
print('===== Finish Mass Fitter =====')


print('===== Start ctau3DErr Fitter =====')
splotter = SPlotter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
splotter.run()
print('===== Finish ctau3DErr Fitter =====')


print('===== Finish run_fitter.py =====')