import os
import sys

# set working directory path - parent of this script
script_path = os.path.abspath(sys.argv[0])
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

from ROOT import gROOT, gSystem
gROOT.SetBatch(True)

# load library
ret = gSystem.Load('../RooMaker.so')
# print("Load status:", ret)
from ROOT import RooMakerRun3PbPb as Maker

# set configurations
input_path = '../../input_flowskim/roots/OniaFlowSkim_isMC0_Run3PbPb_Minbias.root'

# make an instance
maker = Maker(
    cLow = 0, cHigh = 180,
    massLow = 2.6, massHigh = 3.5,
    isMC = False,
    isAccW = False, isEffW = False,
    isTnP = False, isPtW = False,
    hiHFBinEdge = 0,
    mcType = 0, # 0: PR, 1: NP. Ignore for data
    fInputPath = input_path,
    inputChainName = "myTree",
    userTag = "Run3_PbPb_ptWeightFit"
    # dimusign = True, # we always use OS
)

# run
maker.init()
maker.run()
# maker.inChain.Print()

print('--- Finish test_run.py ---')

# --- comments ---
# Written: 250717
# The output will be used for mass fit w.r.t. mass fit for pT weighting
# no weight mc mass fit
# -> no weight mass fit (output of this code)
# -> dn/dpT