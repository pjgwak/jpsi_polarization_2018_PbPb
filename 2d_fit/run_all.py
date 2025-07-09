from ROOT import gROOT, gSystem
gROOT.SetBatch(True)
gSystem.Load('./Analysis.so')

from ROOT import Final2DFit, CtauTrueFit, CtauBkgFit, CtauResFit, CtauErrFit, MassFit, McMassFit

is_mc_mass = False
is_mass = False
is_err = False
is_res = False
is_bkg = False
is_gen_fit = False
is_final = True


# set kinematics
ptLow = 3.0
ptHigh = 6.5
yLow = 1.6
yHigh = 2.4
cLow = 60
cHigh = 180
cos_low = 0
cos_high = 0.1

PR = 0; PRw = 1; fEffW = False; fAccW = False; isPtW = False; isTnP = False


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

if is_bkg:
  print('===== Start ctauBkg Fit =====')
  bkg = CtauBkgFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  bkg.run()
  print('===== Finish ctauBkg Fit =====')

if is_gen_fit:
  print('===== Start ctauTrue Fit =====')
  ctau_true = CtauTrueFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  ctau_true.run()
  print('===== Finish ctauTrue Fit =====')

if is_final:
  print('===== Start 2D Final Fit =====')
  final = Final2DFit(ptLow, ptHigh, yLow, yHigh, cLow, cHigh, cos_low, cos_high, PR, PRw, fEffW, fAccW, isPtW, isTnP)
  final.run()
  print('===== Finish 2D Final Fit =====')

print('===== Finish run_all.py =====')