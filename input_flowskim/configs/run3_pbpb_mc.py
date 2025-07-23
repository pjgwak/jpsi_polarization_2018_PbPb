import os
import sys

# set working directory path - parent of this script
script_path = os.path.abspath(sys.argv[0])
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)

from ROOT import gROOT, gSystem
gROOT.SetBatch(True)

# # load library
ret = gSystem.Load('../FlowSkim.so')
# print("Load status:", ret)
from ROOT import FlowSkimRun3DataPbPb as FlowSkim

# HLT and filters - see headers "git_main/headers/cutsAndBins.h"
# Todo: filter?
trig_sel = -1 # -1: dummy
# trig_sel = 24 # Run3 PbPb Jpsi, MinBias

# set io paths
input_path = '/disk1/Oniatree/polarization/oniatree_5p36/2023_MC/Oniatree_MC1_miniAOD_PbPb_2023.root'
input_tree = 'hionia/myTree'
output_path = 'dummy_now' # for test
output_tree = 'myTree'
user_tag = 'Run3PbPb_Minbias'

# make an instance
skimmer = FlowSkim(
  isMC = True, kTrigSel = trig_sel,
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
skimmer.run(nevt = -1) # -1: all event

print('--- Finish test.py ---')

# --- comments ---
# 