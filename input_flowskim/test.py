import os
import sys

# set working directory path - parent of this script
# script_path = os.path.abspath(sys.argv[0])
# script_dir = os.path.dirname(script_path)
# os.chdir(script_dir)

from ROOT import gROOT, gSystem
gROOT.SetBatch(True)

# # load library
ret = gSystem.Load('FlowSkim.so')
# print("Load status:", ret)
from ROOT import FlowSkimRun3DataPbPb as FlowSkim

# set io paths
input_path = '/disk1/Oniatree/polarization/oniatree_5p36/2023_PbPb/Oniatree_MC0_PbPb_2023_EP_angles.root'
input_tree = 'hionia/myTree'
output_path = 'dummy_now' # for test
output_tree = 'myTree'
user_tag = 'Run3PbPb_singleMuonMinbias16'

# make an instance
skimmer = FlowSkim(
  isMC = False, kTrigSel = 12, # should be changed
  hiHFBinEdge = 0, 
  PDtype = 1, # 1: DB, or 2: DBPeri -> 그냥 string으로 전달하는 게?
  inputFilePath = input_path, inputTreeName = input_tree,
  outputFilePath = output_path, outputTreeName = output_tree,
  userTag = user_tag
  # Todo: 
  #   filters
  #   std::string outputFilePath, outputTreeName;
)

# run
skimmer.run(nevt = 10) # -1: all event

print('--- Finish test.py ---')

# # --- comments ---
# # Written: 250717
# # The output will be used for mass fit w.r.t. mass fit for pT weighting
# # no weight mc mass fit
# # -> no weight mass fit (output of this code)
# # -> dn/dpT


# # make an instance
# maker = Maker(
#     cLow = 0, cHigh = 180,
#     massLow = 2.6, massHigh = 3.5,
#     isMC = False,
#     isAccW = False, isEffW = False,
#     isTnP = False, isPtW = False,
#     hiHFBinEdge = 0,
#     mcType = 0, # 0: PR, 1: NP. Ignore for data
#     fInputPath = input_path,
#     inputChainName = "myTree",
#     userTag = "Run3_PbPb_ptWeightFit"
#     # dimusign = True, # we always use OS
# )