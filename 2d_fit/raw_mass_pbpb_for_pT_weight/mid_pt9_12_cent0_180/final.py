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
from ROOT import Final2DFit

# make an instance
fit = Final2DFit(
    ptLow=9, ptHigh=12,
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
fit.inputFilePath = cfg["data_input_path"]
# fit.pdfTypeMassSig = cfg["default_pdf_mass_sig"] # doubleCB, CBG
# fit.pdfTypeMassBkg = cfg["default_pdf_mass_bkg"] # expo, cheby1, ..., cheby6
# fit.nGauss = cfg['defalut_res_n_gauss'] # 1 - 4
# fit.nExpTrue = cfg['defalut_true_n_exp'] # 1 - 3
# fit.nExpBkg_L = cfg['defalut_bkg_n_exp_L'] # 1 - 2
# fit.nExpBkg_R = cfg['defalut_bkg_n_exp_R'] # 1 - 2

fit.init()

# --- change here ---
# set custom ctau3DMin, ctau3DMax - if you need
# fit.ctau3DMax = 6

fit.initVar('b_Jpsi', 0.1, 0.01, 0.9)
# fit.initVar('sigma_1_A_smass_sig', 1.2, 0.8, 2)
fit.initVar('nSig2D', 1000, 800, 3000)
fit.initVar('nBkg2D', 1000, 400, 3000)
fit.setConstVar('sl1', True)
fit.setConstVar('sl2', True)

# fit.setConstVar('nSig2D', True, 10000)
# fit.setConstVar('fDecayM', False)
# fit.setConstVar('b_Jpsi', True, 0.32)
# fit.setConstVar('mean_Jpsi', True, 3.096)
# fit.setConstVar('nSig2D', True, 10000)
# fit.setConstVar('nBkg2D', True, 5000)
# fit.setConstVar('sigma_1_A', True, 0.2)
# fit.setConstVar('fDecayM', True)
# fit.setConstVar('fDecayP', False)

# pdf parameters
# if fit.nExp == 2:
#     fit.initVar('N_Jpsi_MC', 500000, 400000, 1000000)
#     fit.initVar('lambdaDSS', 0.4, 0.001, 1) # > 0
#     fit.initVar('r_lambda2', 0.5, -2, 2)  # ln(lambda2 / lambda1)
#     fit.initVar('fDSS', 0.5, 0.01, 0.99)

fit.run()

print('--- Finish Final2DFit fit ---')


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