import ROOT
import array
import os
import random
import math

# set batch mode
ROOT.gROOT.SetBatch(True)

# set CMS plot style
rootlogon_path = '../../../rootlogon.C'
if os.path.exists(rootlogon_path):
    # print('Use rootlogon.C')
    ROOT.gROOT.ProcessLine(f'.X {rootlogon_path}')
else:
    print(f'Not exist: {rootlogon_path}')


# get inputs
input_file_name = 'PbPb_ep_angle.root'
tree_name = 'hionia/myTree'
output_prefix = 'plots_' 

file = ROOT.TFile.Open(input_file_name, 'READ')
tree = file.Get(tree_name)

# draw ep, ep_rec and ep_flat on one canvas
c_ep_group = ROOT.TCanvas('c_ep_group', 'ep, ep_rec, ep_flat Distributions', 1500, 500)
c_ep_group.Divide(3, 1) 

# ep_raw
c_ep_group.cd(1)
tree.Draw('ep>>h_ep(100, -2, 2)', '', 'HIST')
h_ep = ROOT.gPad.GetPrimitive('h_ep')

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.03)

if h_ep:
    h_ep.GetXaxis().SetTitle('#Psi_{2,trk} (rad)')
    h_ep.GetXaxis().CenterTitle(True)

    current_max_y = h_ep.GetMaximum()
    if current_max_y > 0: 
        h_ep.SetMaximum(current_max_y * 1.4)
    else:
        h_ep.SetMaximum(1)

    latex.DrawLatex(0.75, 0.85, 'Raw EP Angle')
else:
    print('Error: Histogram h_ep not found in pad 1!')

# ep_rec
c_ep_group.cd(2)
tree.Draw('ep_rec>>h_ep_rec(100, -2, 2)', '', 'HIST')
h_ep_rec = ROOT.gROOT.FindObject('h_ep_rec')
if h_ep_rec:
    h_ep_rec.GetXaxis().SetTitle('#Psi_{2,trk} (rad)')
    h_ep_rec.GetXaxis().CenterTitle(True)

    current_max_y = h_ep_rec.GetMaximum()
    if current_max_y > 0:
        h_ep_rec.SetMaximum(current_max_y * 1.4)
    else:
        h_ep_rec.SetMaximum(1)

    latex.DrawLatex(0.75, 0.85, 'Recentering')
else:
    print('Error: Histogram h_ep_rec not found!')

# ep_flat
c_ep_group.cd(3)
tree.Draw('ep_flat>>h_ep_flat(100, -2, 2)', '', 'HIST')
h_ep_flat = ROOT.gROOT.FindObject('h_ep_flat')
if h_ep_flat:
    h_ep_flat.GetXaxis().SetTitle('#Psi_{2,trk} (rad)')
    h_ep_flat.GetXaxis().CenterTitle(True)

    current_max_y = h_ep_flat.GetMaximum()
    if current_max_y > 0:
        h_ep_flat.SetMaximum(current_max_y * 1.4)
    else:
        h_ep_flat.SetMaximum(1)

    latex.DrawLatex(0.6, 0.85, 'Recentering + Flattening')
else:
    print('Error: Histogram h_ep_flat not found!')

c_ep_group.Update()

# save the plot
ep_group_filename = f'{output_prefix}ep_group'
c_ep_group.SaveAs(f'figs/{ep_group_filename}.png')

# cleaning
file.Close()