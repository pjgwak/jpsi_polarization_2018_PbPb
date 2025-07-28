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
from ROOT import CtauTrueFit

# make an instance
fit = CtauTrueFit(
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
fit.inputFilePath = cfg["mc_gen_only_input_path"]
fit.nExp = cfg['defalut_true_n_exp'] # nExp in ctauTrue decay pdf
fit.nExp = 2 # for test

fit.init()

# --- change here ---
# set custom ctau3DMin, ctau3DMax - if you need
# fit.ctau3DMax = 6

# pdf parameters
if fit.nExp == 2:
    fit.initVar('N_Jpsi_MC', 1000000, 400000, 5000000)
    fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
    fit.initVar('r_lambda2', 0.5, -2, 2)  # ln(lambda2 / lambda1)
    fit.initVar('fDSS', 0.5, 0.01, 0.99)

fit.run()

print('--- Finish ctauTrue fit ---')


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