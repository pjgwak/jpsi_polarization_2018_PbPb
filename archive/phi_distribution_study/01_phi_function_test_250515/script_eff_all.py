import os, shutil
import subprocess as sp

# macros should be run from the folder where the prompt is located
current_dir = os.getcwd()
print(f"current directory: {current_dir}")


# ===== turn on/off flag ===== #
is_init = True # Be careful - it removes figs and roots folders
is_eff = False # also comment out the codes below
is_drawing = True


# ===== recreate folders ===== #
if is_init:
    
    # remove old folders
    if os.path.exists("figs"):
        shutil.rmtree("figs")

    os.makedirs("figs/eff_fwd")
    os.makedirs("figs/eff_mid")
    os.makedirs("figs/eff_all_y")
    
    # if os.path.exists("roots"):
    #     shutil.rmtree("roots")
    if not os.path.exists("roots"):
        os.makedirs("roots")


# ===== calculate correction ===== #
if is_eff:
    def run_eff_macro(data_label="all_event", nevt=-1, mc_type=1, is_tnp=1, is_pt_weight=1):
        cmd = f"root -b -q eff_pb.C'(\"{data_label}\", {int(nevt)}, {int(mc_type)}, {is_tnp}, {is_pt_weight})'"
        sp.run([cmd], shell=True, cwd=current_dir)
    
    run_eff_macro(data_label='all_event', nevt=-1, mc_type=1)


# ===== draw plots ===== #
if is_drawing:    
    sp.run(["python3", "draw_eff_three_rapidities_pb.py"], cwd=current_dir)
    