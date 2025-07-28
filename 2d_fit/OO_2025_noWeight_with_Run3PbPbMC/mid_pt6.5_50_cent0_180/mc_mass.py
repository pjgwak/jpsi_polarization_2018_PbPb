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
from ROOT import McMassFit

# make an instance
fit = McMassFit(
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

# input, user name tag and pdf type
fit.inputFilePath = cfg["mc_input_path"]
fit.DATE = cfg["date_tag"] # you can change it - usually to distinguish wegiht vs no weight
fit.pdfType = cfg["default_pdf_mass_sig"] # doubleCB, CBG
# fit.pdfType = 'CBG' # for test
fit.isWeighted = False # should it go to local_config?

# mc mass fit range - you can change it too!
# fit.massMin = 2.6
fit.massMax = 3.21

fit.init() # set lablels, pdfs

# fit parameters
if fit.pdfType == 'doubleCB':
  fit.initVar('N_Jpsi', 2000000, 1000000, 10000000)
  # fit.initVar('mean', 3.096, 3.096, 3.096) # use default value
  fit.initVar('sigma_1_A', 0.02, 0.001, 0.05)
  fit.initVar('x_A', 2, 1.0, 5.0)
  fit.initVar('alpha_1_A', 1.5, 0.01, 3)
  fit.initVar('n_1_A', 1.5, 1, 3)
  fit.initVar('f', 0.5, 0.01, 0.98)

fit.run()

print('--- Finish mc_mass fit ---')


# -----------------------
# --- free user notes ---
# -----------------------
# if pdfType == "DoubleCB":
#   fit.initVar('N_Jpsi', 500000, 200000, 600000)
#   fit.initVar('mean', 3.096, 3.086, 3.106)
#   fit.initVar('sigma_1_A', 0.01, 0.001, 0.1)
#   fit.initVar('x_A', 1.1, 1, 3)
#   fit.initVar('alpha_1_A', 1.5, 0.2, 5)
#   fit.initVar('n_1_A', 1, 0.2, 5)
#   fit.initVar('f', 0.6, 0.05, 0.95)
# elif pdfType == "CBG":
#   fit.initVar('N_Jpsi', 500000, 200000, 600000)
#   fit.initVar('mean', 3.096, 3.086, 3.106)
#   fit.initVar('sigma_cb', 0.01, 0.001, 0.1)
#   fit.initVar('x_A', 1.1, 1, 3)
#   fit.initVar('alpha_cb', 1.5, 0.2, 5)
#   fit.initVar('n_cb', 1, 0.2, 5)
#   fit.initVar('f', 0.6, 0.05, 0.95)
# else:
#     raise ValueError(f"Unsupported pdfType: {pdfType}")