import ROOT, os
# script가 현재 프롬프트가 위치한 폴더에서 실행되도록 보장
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(script_dir)

# Load macro1 - data
ROOT.gROOT.LoadMacro('skim_to_roodata.C')

# what smaples to skim?
# nominal
is_data = True
is_mc_prompt = True
is_mc_nonprompt = True

# pp MC
is_mc_prompt_pp = False
is_mc_nonprompt_pp = False

# syst.


# parameters
date_label = '250221'
n_evt = -1

c_low = 0; c_high = 180 # centality
mass_low = 2.6; mass_high = 3.5 # Jpsi mass


# ===== nominal ===== #
if is_data:
    ROOT.skim_to_roodata(
        date_label,  # Date label
        n_evt, # how many event to process? -1 for all
        c_low, c_high, # centrality range
        mass_low, mass_high, # centrality range
        True, # dimusign - true: opposite sign muons, false: same sign
        False, # isMC
        1,  # mc_type: 1(Prompt), 2(Nonprompt)
        True, True, True, True, # isAccW, isEffW, isTnP, isPtW,
        0 # HF syst: 0 nominal, 1 HFUp, 2 HFDown
    )
    # isAccW, isEffW, isTnP, isPtW 이 넷은 항상 같이 움직인다. 위에서 처리
    # mc_type, weight_PR 연관 있는지 확인. Data 관련 없으면 .C 코드에서 이름에 조건문 추가.
    # weight_PR (1 prompt, 2 nonprompt) is not used

if is_mc_prompt:
        ROOT.skim_to_roodata(
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
        ROOT.skim_to_roodata(
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

# # ===== syst ===== #