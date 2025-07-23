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
from ROOT import MassFit

# # make an instance
fit = MassFit(
    ptLow=5, ptHigh=6.5,
    yLow=1.6, yHigh=2.4,
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
fit.pdfTypeSig = cfg["default_pdf_mass_sig"] # doubleCB, CBG
fit.pdfTypeBkg = cfg["default_pdf_mass_bkg"] # expo, cheby1, ..., cheby6
# fit.pdfTypeSig = 'CBG' # for test
# fit.pdfTypeBkg = 'expo' # for test
# fit.pdfTypeBkg = 'cheby1' #for test
fit.isWeighted = False

fit.init() # input, lablels, pdfs

# # fit parameters
if fit.pdfTypeSig == 'doubleCB': # signal
  fit.initVar('N_Jpsi', 500, 100, 10000)
  fit.initVar('sigma_1_A', 0.02, 0.0001, 0.2)
  # fixed: x_A, alpha_1_A, n_1_A, f
  # fit.initVar('mean', 3.096, 3.086, 3.106) # use default value
  # fit.initVar('x_A', 1.1, 1, 5) # fixed
  # fit.initVar('alpha_1_A', 1.5, 0.2, 5) # fixed
  # fit.initVar('n_1_A', 1.5, 1, 100) # fixed
  # fit.initVar('f', 0.6, 0.05, 0.95) # fixed
elif fit.pdfTypeSig == "CBG":
  fit.initVar('N_Jpsi', 10000, 2000, 40000)
  fit.initVar('sigma_cb', 0.01, 0.001, 0.1)
  # fit.initVar('mean', 3.096, 3.086, 3.106) # use default value
  # fixed: x_A, alpha_cb, n_cb, f

# bkg
if fit.pdfTypeBkg == 'expo': 
  fit.initVar('N_Bkg', 2000, 10, 20000)
  fit.initVar('lambda', -0.01, -1., 0.)
elif fit.pdfTypeBkg == 'cheby1': 
  fit.initVar('N_Bkg', 2000, 10, 20000)
  fit.initVar('sl1', 0.01, -1, 1)
elif fit.pdfTypeBkg == 'cheby2': 
  # sl1~6 have default values inside MassFit.cpp
  fit.initVar('N_Bkg', 20000, 1000, 50000)
  fit.initVar('sl1', 0.03, -1, 1) # cheby always -1 ~ 1
  fit.initVar('sl2', 0.05, -1, 1)
  # fit.initVar('sl3', 0.01, -1, 1)
  # fit.initVar('sl4', 0.01, -1, 1)
  # fit.initVar('sl5', 0.01, -1, 1)
  # fit.initVar('sl6', 0.01, -1, 1)


fit.run()

print('--- Finish mass fit ---')


# -----------------------
# --- free user notes ---
# -----------------------
# if fit.pdfTypeSig == 'doubleCB': # signal
#   fit.initVar('N_Jpsi', 50000, 10000, 100000)
#   fit.initVar('sigma_1_A', 0.01, 0.001, 0.1)
  # fit.initVar('mean', 3.096, 3.086, 3.106) # use default value
  # fit.initVar('x_A', 1.1, 1, 3) # fixed
  # fit.initVar('alpha_1_A', 1.5, 0.2, 5) # fixed
  # fit.initVar('n_1_A', 1, 0.2, 5) # fixed
  # fit.initVar('f', 0.6, 0.05, 0.95) # fixed
# elif fit.pdfTypeSig == "CBG":
#   fit.initVar('N_Jpsi', 500000, 200000, 600000)
#   fit.initVar('mean', 3.096, 3.086, 3.106)
#   fit.initVar('sigma_cb', 0.01, 0.001, 0.1)
#   fit.initVar('x_A', 1.1, 1, 3)
#   fit.initVar('alpha_cb', 1.5, 0.2, 5)
#   fit.initVar('n_cb', 1, 0.2, 5)
#   fit.initVar('f', 0.6, 0.05, 0.95)

# if fit.pdfTypeBkg == 'expo': 
#   fit.initVar('N_Bkg', 100000, 50000, 200000)
#   fit.initVar('lambda', -0.01, -1., 0.)
# elif fit.pdfTypeBkg == 'cheby2': 
#   # bkg - please change pdf name and comment in (out) according to order
#   fit.initVar('N_Bkg', 500000, 200000, 600000)
#   fit.initVar('sl1', 0.01, -1, 1) # cheby always -1 ~ 1
#   fit.initVar('sl2', 0.01, -1, 1)
#   fit.initVar('sl3', 0.01, -1, 1)
#   fit.initVar('sl4', 0.01, -1, 1)
#   fit.initVar('sl5', 0.01, -1, 1)
#   fit.initVar('sl6', 0.01, -1, 1)