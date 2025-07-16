from ROOT import gROOT, gSystem
gROOT.SetBatch(True)

# load library
ret = gSystem.Load('RooMaker.so')
# print("Load status:", ret)
from ROOT import RooMakerRun3DataPbPb as Maker

# set configurations
input_path = '../archive/makeFlowSkim/13_MC_gen/OniaFlowSkim.Data2023.AlgDefault.root'

# make an instance
maker = Maker(
    cLow = 0, cHigh = 200,
    massLow = 2.6, massHigh = 3.5,
    isMC = False,
    isAccW = False, isEffW = False,
    isTnP = False, isPtW = False,
    hiHFBinEdge = 0,
    mcType = 0, # 0: PR, 1: NP. Ignore for data
    fInputPath = input_path,
    inputChainName = "myTree"
    # dimusign = True, # we always use OS
)

# run
maker.init()
maker.run()
# maker.inChain.Print()

print('--- Finish test_run.py ---')