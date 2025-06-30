from ROOT import gROOT, gSystem
gROOT.SetBatch(True)
gSystem.Load('./Analysis.so')

from ROOT import McMassFit, MassFit


# from ROOT import McMassFitter, MassFitter, SPlotter, ResolutionFitter, BkgFitter, TrueFitter
# from ROOT import FinalFitter


is_mc_mass = False
is_mass = True
is_err = False
is_res = False
is_bkg = False
is_true_fit = False
is_final = False


# set kinematics
ptLow = 3.0
ptHigh = 6.5
yLow = 1.6
yHigh = 2.4
cLow = 60
cHigh = 180
cos_low = 0.5
cos_high = 0.6

PR = 0; PRw = 1; fEffW = False; fAccW = False; isPtW = False; isTnP = False


print('===== Start run_fitter.py =====')
if is_mc_mass:
  print('===== Start MC Fitter =====')
  mc_fitter = McMassFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, True)
  mc_fitter.run()
  print('===== Finish MC Fitter =====\n\n\n')

if is_mass:
  print('===== Start Mass Fitter =====')
  mass_fitter = MassFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, True)
  mass_fitter.run()
  # mass_fitter.allInOne()
  print('===== Finish Mass Fitter =====\n\n\n')

# if is_err:
#   print('===== Start ctau3DErr Fitter =====')
#   # splotter = SPlotter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   # splotter.run()
#   print('===== Finish ctau3DErr Fitter =====')

# if is_res:
#   print('===== Start ctau3DRes Fitter =====')
#   res_fitter = ResolutionFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   res_fitter.run()
#   print('===== Finish ctau3DRes Fitter =====')

# if is_bkg:
#   print('===== Start ctauBkg Fitter =====')
#   bkg_fitter = BkgFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   bkg_fitter.run()
#   print('===== Finish ctauBkg Fitter =====')

# if is_true_fit:
#   print('===== Start ctauTrue Fitter =====')
#   true_fitter = TrueFitter(ptLow, ptHigh, yLow, yHigh)
#   true_fitter.run()
#   print('===== Finish ctauTrue Fitter =====')

# if is_final:
#   print('===== Start ctauTrue Fitter =====')
#   final_fitter = FinalFitter(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   final_fitter.run()
#   print('===== Finish ctauTrue Fitter =====')

print('===== Finish run_fitter.py =====')