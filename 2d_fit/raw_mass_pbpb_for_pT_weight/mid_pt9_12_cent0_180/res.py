# --- please change those items
# 1. kinematic values such as ptLow, ptHigh ... ")
# 2. ctauResMin and ctauResMax if you need
# 3. fit parameters according to nGauss
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
from ROOT import CtauResFit

# make an instance
fit = CtauResFit(
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

# user name tag and pdf info
fit.DATE = cfg["date_tag"] # you can change it - usually to distinguish wegiht vs no weight
fit.nGauss = cfg['defalut_res_n_gauss'] # nGauss in Res pdf
fit.nGauss = 2 # for test

fit.init()

# --- change here ---
# Todo: apply cut
# set user custom ctuRes min range - should min < max - See the printed log
# fit.ctauResMin = -2.5

# pdf parameters
if fit.nGauss == 2:
    fit.initVar('N_Jpsi', 200, 10, 1000)
    fit.initVar('s1_CtauRes', 0.4, 0.001, 2)
    fit.initVar('rS21_CtauRes', 1.5, 1.0, 5.0) # ratio should > 1.0
    fit.initVar('f_CtauRes', 0.5, 0.01, 0.99)


fit.run()

print('--- Finish ctauRes fit ---')


# -----------------------
# --- free user notes ---
# -----------------------
# if fit.nGauss == 1:
#     fit.initVar('N_Jpsi', 30000, 10, 100000)
#     fit.initVar('s1_CtauRes', 0.4, 0.001, 0.9)
# if fit.nGauss == 2:
#     fit.initVar('N_Jpsi', 30000, 10, 100000)
#     fit.initVar('s1_CtauRes', 0.4, 0.001, 0.9)
#     fit.initVar('rS21_CtauRes', 1.5, 1.0, 3.0) # ratio should > 1.0
#     fit.initVar('f_CtauRes', 0.5, 0.01, 0.99)
# if fit.nGauss == 3:
#     fit.initVar('N_Jpsi', 30000, 10, 100000)
#     fit.initVar('s1_CtauRes', 0.4, 0.001, 0.9)
#     fit.initVar('rS21_CtauRes', 1.5, 1.0, 3.0) # ratio should > 1.0
#     fit.initVar('rS32_CtauRes', 1.5, 1.0, 3.0)
#     fit.initVar('f_CtauRes', 0.3, 0.01, 0.99)
#     fit.initVar('f2_CtauRes', 0.3, 0.01, 0.99)
# if fit.nGauss == 4:
#     fit.initVar('N_Jpsi', 30000, 10, 100000)
#     fit.initVar('s1_CtauRes', 0.4, 0.001, 0.9)
#     fit.initVar('rS21_CtauRes', 1.5, 1.0, 3.0) # ratio should > 1.0
#     fit.initVar('rS32_CtauRes', 1.5, 1.0, 3.0)
#     fit.initVar('rS43_CtauRes', 1.5, 1.0, 3.0)
#     fit.initVar('f_CtauRes', 0.3, 0.01, 0.99)
#     fit.initVar('f2_CtauRes', 0.3, 0.01, 0.99)
#     fit.initVar('f3_CtauRes', 0.3, 0.01, 0.99)