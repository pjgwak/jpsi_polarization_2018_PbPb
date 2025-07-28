# --- please change those items
# 1. kinematic values such as ptLow, ptHigh ... (search "changeThisNumber")
# 2. signal and bkg paramters (see below the "fit parameters")
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
from ROOT import CtauErrFit

# make an instance
fit = CtauErrFit(
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

# # input, user name tag and pdf type
fit.inputFilePath = cfg["data_input_path"]
fit.DATE = cfg["date_tag"] # you can change it - usually to distinguish wegiht vs no weight

# fit.init() # input, lablels, pdfs

# set user custom ctuErr max range - should max > min - See the printed log
# fit.ctauErrMax = 0.045


fit.run()

print('--- Finish ctauErr fit ---')


# -----------------------
# --- free user notes ---
# -----------------------