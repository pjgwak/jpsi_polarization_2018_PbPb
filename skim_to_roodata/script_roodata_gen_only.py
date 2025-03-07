import ROOT, os
# script가 현재 프롬프트가 위치한 폴더에서 실행되도록 보장
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(script_dir)

# Load macro1 - data
ROOT.gROOT.LoadMacro('skim_to_roodata_gen_only.C')

# what smaples to skim?
# nominal
# is_mc_prompt = False # Basically, we don't need it.
is_mc_nonprompt = True

# pp MC
is_mc_prompt_pp = False
is_mc_nonprompt_pp = False

# syst.


# parameters
date_label = '250307'
n_evt = -1

# Warning: centrality range is not applied in practice - just to keep naming sense
c_low = 0; c_high = 180 # centrality
mass_low = 2.6; mass_high = 3.5 # Jpsi mass


# ===== nominal ===== #
if is_mc_prompt:
        ROOT.skim_to_roodata_gen_only(
        date_label,  # Date label
        n_evt, # how many event to process? -1 for all
        c_low, c_high, # centrality range
        mass_low, mass_high, # centrality range
        True, # dimusign - true: opposite sign muons, false: same sign
        True, # isMC
        1,  # mc_type: 1(Prompt), 2(Nonprompt)
        True, True, True, True, # isAccW, isEffW, isTnP, isPtW,
        0 # HF syst: 0 nominal, 1 HFUp, 2 HFDown
    )
if is_mc_nonprompt:
        ROOT.skim_to_roodata_gen_only(
        date_label,  # Date label
        n_evt, # how many event to process? -1 for all
        c_low, c_high, # centrality range
        mass_low, mass_high, # centrality range
        True, # dimusign - true: opposite sign muons, false: same sign
        True, # isMC
        2,  # mc_type: 1(Prompt), 2(Nonprompt)
        True, True, True, True, # isAccW, isEffW, isTnP, isPtW,
        0 # HF syst: 0 nominal, 1 HFUp, 2 HFDown
    )


# ===== syst ===== #