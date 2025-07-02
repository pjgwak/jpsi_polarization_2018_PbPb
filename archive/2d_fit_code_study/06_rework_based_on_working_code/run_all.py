from ROOT import gROOT, gSystem
gROOT.SetBatch(True)
gSystem.Load('./Analysis.so')

from ROOT import McMassFit, MassFit, CtauErrFit, CtauResFit


# from ROOT import McMassFitter, MassFitter, SPlotter, ResolutionFitter, BkgFitter, TrueFitter
# from ROOT import FinalFitter


is_mc_mass = False
is_mass = False
is_err = False
is_res = True
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

PR = 0; PRw = 1; fEffW = False; fAccW = False; isPtW = False; isTnP = True


print('===== Start run_all.py =====')
if is_mc_mass:
  print('===== Start MC Fit =====')
  mc_mass = McMassFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  mc_mass.run()
  print('===== Finish MC Fit =====\n\n\n')

if is_mass:
  print('===== Start Mass Fit =====')
  mass = MassFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  mass.run()
  print('===== Finish Mass Fit =====\n\n\n')

if is_err:
  print('===== Start ctau3DErr Fit =====')
  ctau_err = CtauErrFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  ctau_err.run()
  print('===== Finish ctau3DErr Fit =====')

if is_res:
  print('===== Start ctau3DRes Fit =====')
  res = CtauResFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  res.run()
  print('===== Finish ctau3DRes Fit =====')

# if is_bkg:
#   print('===== Start ctauBkg Fit =====')
#   bkg_Fit = BkgFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   bkg_Fit.run()
#   print('===== Finish ctauBkg Fit =====')

# if is_true_fit:
#   print('===== Start ctauTrue Fit =====')
#   true_Fit = TrueFit(ptLow, ptHigh, yLow, yHigh)
#   true_Fit.run()
#   print('===== Finish ctauTrue Fit =====')

# if is_final:
#   print('===== Start 2D Final Fit =====')
#   final_Fit = FinalFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high)
#   final_Fit.run()
#   print('===== Finish 2D Final Fit =====')

print('===== Finish run_all.py =====')