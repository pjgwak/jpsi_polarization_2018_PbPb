import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent)) # to call local_config
from local_config import get_common_config

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
    ptLow=3.0, ptHigh=6.5,
    yLow=1.6, yHigh=2.4,
    cLow=0, cHigh=180,
    cosLow=-1.0, cosHigh=1.0,
    PR=cfg["default_PR"], # 이하 local_config로 이동?
    PRw=cfg["default_PRw"],
    fEffW=False,
    fAccW=False,
    isPtW=False,
    isTnP=False
)

# input, user name tag and pdf type
fit.inputFilePath = cfg["input_path"]
fit.DATE = cfg["date_tag"] # you can change it - usually to distinguish wegiht vs no weight
fit.pdfType = cfg["default_pdf_type"] # doubleCB, CBG
fit.isWeighted = False

# mc mass fit range
fit.massMax = 3.27

fit.init() # set lablels, pdfs

# fit parameters
if fit.pdfType == "doubleCB":
  fit.initVar('N_Jpsi', 500000, 200000, 600000)
  fit.initVar('mean', 3.096, 3.086, 3.106)
  fit.initVar('sigma_1_A', 0.01, 0.001, 0.1)
  fit.initVar('x_A', 1.1, 1, 3)
  fit.initVar('alpha_1_A', 1.5, 0.2, 5)
  fit.initVar('n_1_A', 1, 0.2, 5)
  fit.initVar('f', 0.6, 0.05, 0.95)

fit.run()

# --- set plotting configuration ---
# drawing options etc
# positon
# legend
# etc.
# 그냥 root 파일에 필요한 거 저장해놓고 한 번 더 돌릴까?

print('--- Finish mc_mass fit ---')


# ---- Not used ----
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