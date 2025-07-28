# --- please change those items
# 1. kinematic values such as ptLow, ptHigh ... ")
# 2. ctauResMin and ctauResMax if you need
# 3. fit parameters according to nExp
# ---

import os, sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent)) # to call local_config
from local_config import get_common_config

# set working directory path - parent of this script
script_path = os.path.abspath(sys.argv[0])
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

import shutil

from ROOT import gSystem, gROOT
gROOT.SetBatch(True)

# load local cofnig
cfg = get_common_config()

# load .so
gSystem.Load(cfg['shared_lib_path'])

# import cpp class
from ROOT import CtauBkgFit

# make an instance
fit = CtauBkgFit(
    ptLow=6.5, ptHigh=50,
    yLow=0, yHigh=1.6,
    cLow=0, cHigh=180,
    cosLow=-1.0, cosHigh=1.0,
    PR=cfg["default_PR"],
    PRw=cfg["default_PRw"],
    fEffW=cfg['default_fEffW'],
    fAccW=cfg['default_fAccW'],
    isPtW=cfg['default_isPtW'],
    isTnP=cfg['default_isTnP']
)

# # user name tag and pdf info
fit.DATE = cfg["date_tag"] # you can change it - usually to distinguish wegiht vs no weight
fit.nGauss = cfg['defalut_true_n_exp'] # nGauss of resolution - it must same with nGauss of ctauRes fit
fit.nGauss = 2
fit.nExp_L = 1; fit.nExp_R = 1

fit.init()

# --- change here ---
# set custom ctauMin, ctauMax - if you need
fit.ctauMin = -2; fit.ctauMax = 2

# pdf parameters
fit.initVar('N_BkgCtau', 300, 1, 10000)
# fit.initVar('b_Bkg', 0.8, 0.3, 0.95)
fit.initVar('fDecayP', 0.2, 0.01, 0.9)
fit.initVar('fDecayM', 0.6, 0.01, 0.9)

# make gauss sigma free - when the peak is too sharp
# fit.initVar('s1_CtauRes', 0.93, 0.8, 3)
# fit.setConstVar('s1_CtauRes', False)


if fit.nExp_L == 1:
    fit.initVar('lambdaDF_Bkg1', 0.2, 0.001, 1.0)
    # fit.setConstVar('lambdaDF_Bkg1', True, 0.4)

if fit.nExp_L == 2:
    fit.initVar('lambdaDF_Bkg1', 1.5, 0.001, 2.0)
    fit.initVar('rDF12', 1.5, 0.001, 5.0)
    fit.initVar('fDF12', 0.4, 0.01, 1.0)

if fit.nExp_R == 1:
    fit.initVar('lambdaDSS_Bkg1', 0.4, 0.001, 2.0)
if fit.nExp_R == 2:
    fit.initVar('lambdaDSS_Bkg1', 0.4, 0.001, 5.0)
    fit.initVar('rDSS12', 4, 1, 10.0)
    fit.initVar('fDSS12', 0.4, 0.01, 1.0)

# if fit.nExp_C == 1:
#     fit.initVar('lambdaDDS_Bkg1', 0.4, 0.01, 2.0)
# if fit.nExp_C == 2:
#     fit.initVar('lambdaDDS_Bkg1', 0.4, 0.01, 2.0)
#     fit.initVar('lambdaDDS_Bkg2', 1.5, 0.001, 2.0)
#     fit.initVar('fDDS12', 0.4, 0.01, 1.0)


fit.run()

print('--- Finish ctauBkg fit ---')


# -----------------------
# --- free user notes ---
# -----------------------
# if fit.nExp == 1:
#     fit.initVar('N_Jpsi_MC', 500000, 400000, 1000000)
#     fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
# if fit.nExp == 2:
#     fit.initVar('N_Jpsi_MC', 500000, 400000, 1000000)
#     fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
#     fit.initVar('r_lambda2', 0.5, -2, 2)  # ln(lambda2 / lambda1)
#     fit.initVar('fDSS', 0.5, 0.01, 0.99)
# if fit.nExp == 3:
#     fit.initVar('N_Jpsi_MC', 500000, 400000, 1000000)
#     fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
#     fit.initVar('r_lambda2', 0.5, 0.001, 2)  # ln(lambda2 / lambda1)
#     fit.initVar('r_lambda3', 0.5, 0.001, 2)  
#     fit.initVar('fDSS', 0.5, 0.01, 0.99)
#     fit.initVar('f2_DSS', 0.5, 0.01, 0.99)
# if fit.nExp == 4:
#     fit.initVar('N_Jpsi_MC', 500000, 400000, 1000000)
#     fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
#     fit.initVar('r_lambda2', 0.5, 0.001, 2)  # ln(lambda2 / lambda1) > 0
#     fit.initVar('r_lambda3', 0.5, 0.001, 2)
#     fit.initVar('r_lambda4', 0.5, 0.001, 2)
#     fit.initVar('fDSS', 0.5, 0.01, 0.99)
#     fit.initVar('f2_DSS', 0.5, 0.01, 0.99)
#     fit.initVar('f3_DSS', 0.5, 0.01, 0.99)