import os, shutil
import subprocess as sp

# macros should be run from the folder where the prompt is located
current_dir = os.getcwd()
print(f"current directory: {current_dir}")


# ===== turn on/off flag ===== #
is_init = False # Be careful - it removes figs and roots folders
is_eff = False
is_drawing = True


# ===== recreate folders ===== #
if is_init:
    # remove old folders
    if os.path.exists("figs"):
        shutil.rmtree("figs")

    if os.path.exists("roots"):
        shutil.rmtree("roots")

    # make folders
    os.makedirs("figs/eff_fwd")
    os.makedirs("figs/eff_mid")
    os.makedirs("figs/eff_all_y")
    os.makedirs("figs/stacked_plots")
    os.makedirs("roots")


# ===== calculate correction ===== #
if is_eff:
    def run_eff_macro(data_label="all_event", nevt=-1, mc_type=1):
        cmd = f"root -b -q eff_pb.C'(\"{data_label}\", {int(nevt)}, {int(mc_type)})'"
        sp.run([cmd], shell=True, cwd=current_dir)
    
    run_eff_macro(data_label='all_event', nevt=-1, mc_type=1)


# ===== draw plots ===== #
if is_drawing:    
    sp.run(["python3", "draw_eff_three_rapidities_pb.py"], cwd=current_dir)
    sp.run(["python3", "draw_eff_stacked_plots_pb.py"], cwd=current_dir)
    