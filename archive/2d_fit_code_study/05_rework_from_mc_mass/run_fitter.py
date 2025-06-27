import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('./McFitter.so')

from ROOT import McMassFitter

print('===== Start run_fitter.py =====')

# set kinematics
ptLow = 3.0
ptHigh = 6.5
yLow = 1.6
yHigh = 2.4
cLow = 60
cHigh = 180
cos_low = 0.5
cos_high = 0.6

# make a McMassFitter instance
print(f'Creating McMassFitter for pt: {ptLow}-{ptHigh}, y: {yLow}-{yHigh}')
fitter = McMassFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)

# call run()
print('Running the fit...')
fitter.run()

print('===== Finish run_fitter.py =====')