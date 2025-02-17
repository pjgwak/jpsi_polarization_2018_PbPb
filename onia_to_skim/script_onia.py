import ROOT, os
# script가 현재 프롬프트가 위치한 폴더에서 실행되도록 보장
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(script_dir)

# Load macro1 - data
ROOT.gROOT.LoadMacro('onia_to_skim.C')

# what smaples to skim?
is_data_double_muon = True
is_data_peripheral = True
is_mc_prompt = True
is_mc_nonprompt = True
# You have to merge data DB and Peri samples

# pp MC
is_mc_prompt_pp = False
is_mc_nonprompt_pp = False

# parameters
date_label = '250217'
n_evt = -1 # -1 for all events

if is_data_double_muon:
    ROOT.onia_to_skim(
        date_label, # Date label
        n_evt, # number of event to use. -1:all event
        int(False), # isMC
        1, # MCtype (MC only): 1(Signal) or 2(NP)
        1, # PDtype (Data only): DoubleMuon(1) or Peripheral(2)
        0 # hiHFBinEdge: HF systematics. Norminal(0), Up(1), Down(2)
    )

if is_data_peripheral:
    ROOT.onia_to_skim(
        date_label, # Date label
        n_evt, # number of event to use. -1:all event
        int(False), # isMC
        1, # MCtype (MC only): 1(Signal) or 2(NP)
        2, # PDtype (Data only): DoubleMuon(1) or Peripheral(2)
        0 # hiHFBinEdge: HF systematics. Norminal(0), Up(1), Down(2)
    )

if is_mc_prompt:
    ROOT.onia_to_skim(
        date_label, # Date label
        n_evt, # number of event to use. -1:all event
        int(True), # isMC
        1, # MCtype (MC only): 1(Signal) or 2(NP)
        1, # PDtype (Data only): DoubleMuon(1) or Peripheral(2)
        0 # hiHFBinEdge: HF systematics. Norminal(0), Up(1), Down(2)
    )

if is_mc_nonprompt:
    ROOT.onia_to_skim(
        date_label, # Date label
        n_evt, # number of event to use. -1:all event
        int(True), # isMC
        2, # MCtype (MC only): 1(Signal) or 2(NP)
        1, # PDtype (Data only): DoubleMuon(1) or Peripheral(2)
        0 # hiHFBinEdge: HF systematics. Norminal(0), Up(1), Down(2)
    )