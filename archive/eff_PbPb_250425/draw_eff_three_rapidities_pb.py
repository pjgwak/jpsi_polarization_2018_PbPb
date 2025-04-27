from module_draw_fcns import *
import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D, TPad
from enum import IntEnum

# ==== Please check the pt_bins ==== #
fwd_pt_bins = ('pT3_6', 'pT6_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')
mid_pt_bins = ('pT6p5_9', 'pT9_12', 'pT12_15', 'pT15_20', 'pT20-50')

class Axis(IntEnum):
    cent = 0
    cos = 1
    phi = 2
    pt = 3

# ===== set macro configure ===== #
# no stat box for histogram
ROOT.gStyle.SetOptStat(0)

# batch mode - don't show the plot
ROOT.gROOT.SetBatch(True)


# ===== read inputs ===== #
# bring input root file
in_file_pb = TFile('roots/eff_PbPb_Jpsi_PtW1_tnp1_all_event.root')

# ===== bring hists ===== #
# fwd: 1.6 < |y| < 2.4, 3 < pT < 50
h4_fwd_lab_num_pb = in_file_pb.Get('fwd_lab_num')
h4_fwd_lab_den_pb = in_file_pb.Get('fwd_lab_den')
h4_fwd_hx_num_pb = in_file_pb.Get('fwd_hx_num')
h4_fwd_hx_den_pb = in_file_pb.Get('fwd_hx_den')
h4_fwd_cs_num_pb = in_file_pb.Get('fwd_cs_num')
h4_fwd_cs_den_pb = in_file_pb.Get('fwd_cs_den')
h4_fwd_ep_num_pb = in_file_pb.Get('fwd_ep_num')
h4_fwd_ep_den_pb = in_file_pb.Get('fwd_ep_den')

# mid: |y| < 1.6, 6.5 < pT < 50
h4_mid_lab_num_pb = in_file_pb.Get('mid_lab_num')
h4_mid_lab_den_pb = in_file_pb.Get('mid_lab_den')
h4_mid_hx_num_pb = in_file_pb.Get('mid_hx_num')
h4_mid_hx_den_pb = in_file_pb.Get('mid_hx_den')
h4_mid_cs_num_pb = in_file_pb.Get('mid_cs_num')
h4_mid_cs_den_pb = in_file_pb.Get('mid_cs_den')
h4_mid_ep_num_pb = in_file_pb.Get('mid_ep_num')
h4_mid_ep_den_pb = in_file_pb.Get('mid_ep_den')

# all_y: |y| < 2.4, 6.5 < pT < 50
h4_all_y_lab_num_pb = in_file_pb.Get('all_y_lab_num')
h4_all_y_lab_den_pb = in_file_pb.Get('all_y_lab_den')
h4_all_y_hx_num_pb = in_file_pb.Get('all_y_hx_num')
h4_all_y_hx_den_pb = in_file_pb.Get('all_y_hx_den')
h4_all_y_cs_num_pb = in_file_pb.Get('all_y_cs_num')
h4_all_y_cs_den_pb = in_file_pb.Get('all_y_cs_den')
h4_all_y_ep_num_pb = in_file_pb.Get('all_y_ep_num')
h4_all_y_ep_den_pb = in_file_pb.Get('all_y_ep_den')

# 1d map - y
h_y_lab_num_pb = in_file_pb.Get('y_lab_num')
h_y_lab_den_pb = in_file_pb.Get('y_lab_den')


# ===== draw 1d plots ===== #
# rapidity
draw_1d_y_maps(h_num=h_y_lab_num_pb, h_den=h_y_lab_den_pb, legend1='', legend2='', x_title='y', y_title='# of dimuon', cor_y_max=1.6, save_name='eff_y_lab_pb')
# print('Test run finished')
# exit(1)

# cos
draw_1d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='# of dimuon', save_name='eff_fwd/cos_lab_pb')
draw_1d_hist_maps(h4_num=h4_fwd_hx_num_pb, h4_den=h4_fwd_hx_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='# of dimuon', save_name='eff_fwd/cos_hx_pb')
draw_1d_hist_maps(h4_num=h4_fwd_cs_num_pb, h4_den=h4_fwd_cs_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='# of dimuon', 
save_name='eff_fwd/cos_cs_pb')
draw_1d_hist_maps(h4_num=h4_fwd_ep_num_pb, h4_den=h4_fwd_ep_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='# of dimuon', 
save_name='eff_fwd/cos_ep_pb')

draw_1d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='# of dimuon', save_name='eff_mid/cos_lab_pb')
draw_1d_hist_maps(h4_num=h4_mid_hx_num_pb, h4_den=h4_mid_hx_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='# of dimuon', save_name='eff_mid/cos_hx_pb')
draw_1d_hist_maps(h4_num=h4_mid_cs_num_pb, h4_den=h4_mid_cs_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='# of dimuon', 
save_name='eff_mid/cos_cs_pb')
draw_1d_hist_maps(h4_num=h4_mid_ep_num_pb, h4_den=h4_mid_ep_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='# of dimuon', 
save_name='eff_mid/cos_ep_pb')

draw_1d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='# of dimuon', save_name='eff_all_y/cos_lab_pb')
draw_1d_hist_maps(h4_num=h4_all_y_hx_num_pb, h4_den=h4_all_y_hx_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='# of dimuon', save_name='eff_all_y/cos_hx_pb')
draw_1d_hist_maps(h4_num=h4_all_y_cs_num_pb, h4_den=h4_all_y_cs_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='# of dimuon', 
save_name='eff_all_y/cos_cs_pb')
draw_1d_hist_maps(h4_num=h4_all_y_ep_num_pb, h4_den=h4_all_y_ep_den_pb, axis=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='# of dimuon', 
save_name='eff_all_y/cos_ep_pb')

# phi
draw_1d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{Lab}', y_title='# of dimuon', save_name='eff_fwd/phi_lab_pb')
draw_1d_hist_maps(h4_num=h4_fwd_hx_num_pb, h4_den=h4_fwd_hx_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{hx}', y_title='# of dimuon', save_name='eff_fwd/phi_hx_pb')
draw_1d_hist_maps(h4_num=h4_fwd_cs_num_pb, h4_den=h4_fwd_cs_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{cs}', y_title='# of dimuon', 
save_name='eff_fwd/phi_cs_pb')
draw_1d_hist_maps(h4_num=h4_fwd_ep_num_pb, h4_den=h4_fwd_ep_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{ep}', y_title='# of dimuon', 
save_name='eff_fwd/phi_ep_pb')

draw_1d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{Lab}', y_title='# of dimuon', save_name='eff_mid/phi_lab_pb')
draw_1d_hist_maps(h4_num=h4_mid_hx_num_pb, h4_den=h4_mid_hx_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{hx}', y_title='# of dimuon', save_name='eff_mid/phi_hx_pb')
draw_1d_hist_maps(h4_num=h4_mid_cs_num_pb, h4_den=h4_mid_cs_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{cs}', y_title='# of dimuon', 
save_name='eff_mid/phi_cs_pb')
draw_1d_hist_maps(h4_num=h4_mid_ep_num_pb, h4_den=h4_mid_ep_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{ep}', y_title='# of dimuon', 
save_name='eff_mid/phi_ep_pb')

draw_1d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{Lab}', y_title='# of dimuon', save_name='eff_all_y/phi_lab_pb')
draw_1d_hist_maps(h4_num=h4_all_y_hx_num_pb, h4_den=h4_all_y_hx_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{hx}', y_title='# of dimuon', save_name='eff_all_y/phi_hx_pb')
draw_1d_hist_maps(h4_num=h4_all_y_cs_num_pb, h4_den=h4_all_y_cs_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{cs}', y_title='# of dimuon', 
save_name='eff_all_y/phi_cs_pb')
draw_1d_hist_maps(h4_num=h4_all_y_ep_num_pb, h4_den=h4_all_y_ep_den_pb, axis=Axis.phi, legend1='', legend2='', x_title='#phi_{ep}', y_title='# of dimuon', 
save_name='eff_all_y/phi_ep_pb')

# pt
draw_1d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis=Axis.pt, legend1='', legend2='', x_title='p_{T} (GeV/c)', y_title='# of dimuon', save_name='eff_fwd/pT_lab_pb')
draw_1d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis=Axis.pt, legend1='', legend2='', x_title='p_{T} (GeV/c)', y_title='# of dimuon', save_name='eff_mid/pT_lab_pb')
draw_1d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis=Axis.pt, legend1='', legend2='', x_title='p_{T} (GeV/c)', y_title='# of dimuon', save_name='eff_all_y/pT_lab_pb')

# cent
draw_1d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis=Axis.cent, legend1='', legend2='', x_title='centrality', y_title='# of dimuon', save_name='eff_fwd/cent_lab_pb')
draw_1d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis=Axis.cent, legend1='', legend2='', x_title='centrality', y_title='# of dimuon', save_name='eff_mid/cent_lab_pb')
draw_1d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis=Axis.cent, legend1='', legend2='', x_title='centrality', y_title='# of dimuon', save_name='eff_all_y/cent_lab_pb')


# ===== draw 2d plots ===== #
# cos vs pT
draw_2d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='p_{T} (GeV/c)', save_name='eff_fwd/2d_cos_pT_lab_pb')
draw_2d_hist_maps(h4_num=h4_fwd_hx_num_pb, h4_den=h4_fwd_hx_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='p_{T} (GeV/c)', save_name='eff_fwd/2d_cos_pT_hx_pb')
draw_2d_hist_maps(h4_num=h4_fwd_cs_num_pb, h4_den=h4_fwd_cs_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='p_{T} (GeV/c)', save_name='eff_fwd/2d_cos_pT_cs_pb')
draw_2d_hist_maps(h4_num=h4_fwd_ep_num_pb, h4_den=h4_fwd_ep_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='p_{T} (GeV/c)', save_name='eff_fwd/2d_cos_pT_ep_pb')

draw_2d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='p_{T} (GeV/c)', save_name='eff_mid/2d_cos_pT_lab_pb')
draw_2d_hist_maps(h4_num=h4_mid_hx_num_pb, h4_den=h4_mid_hx_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='p_{T} (GeV/c)', save_name='eff_mid/2d_cos_pT_hx_pb')
draw_2d_hist_maps(h4_num=h4_mid_cs_num_pb, h4_den=h4_mid_cs_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='p_{T} (GeV/c)', save_name='eff_mid/2d_cos_pT_cs_pb')
draw_2d_hist_maps(h4_num=h4_mid_ep_num_pb, h4_den=h4_mid_ep_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='p_{T} (GeV/c)', save_name='eff_mid/2d_cos_pT_ep_pb')

draw_2d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='p_{T} (GeV/c)', save_name='eff_all_y/2d_cos_pT_lab_pb')
draw_2d_hist_maps(h4_num=h4_all_y_hx_num_pb, h4_den=h4_all_y_hx_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='p_{T} (GeV/c)', save_name='eff_all_y/2d_cos_pT_hx_pb')
draw_2d_hist_maps(h4_num=h4_all_y_cs_num_pb, h4_den=h4_all_y_cs_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='p_{T} (GeV/c)', save_name='eff_all_y/2d_cos_pT_cs_pb')
draw_2d_hist_maps(h4_num=h4_all_y_ep_num_pb, h4_den=h4_all_y_ep_den_pb, axis1=Axis.pt, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='p_{T} (GeV/c)', save_name='eff_all_y/2d_cos_pT_ep_pb')


# ===== Fwd cos vs phi - accroding to pT bins ===== # - 6
# pT integrated
draw_2d_hist_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_fwd/2d_cos_phi_lab_pb')
draw_2d_hist_maps(h4_num=h4_fwd_hx_num_pb, h4_den=h4_fwd_hx_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='#phi_{HX} (rad)', save_name='eff_fwd/2d_cos_phi_hx_pb')
draw_2d_hist_maps(h4_num=h4_fwd_cs_num_pb, h4_den=h4_fwd_cs_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='#phi_{CS} (rad)', save_name='eff_fwd/2d_cos_phi_cs_pb')
draw_2d_hist_maps(h4_num=h4_fwd_ep_num_pb, h4_den=h4_fwd_ep_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='#phi_{EP} (rad)', save_name='eff_fwd/2d_cos_phi_ep_pb')

draw_2d_hist_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_mid/2d_cos_phi_lab_pb')
draw_2d_hist_maps(h4_num=h4_mid_hx_num_pb, h4_den=h4_mid_hx_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='#phi_{HX} (rad)', save_name='eff_mid/2d_cos_phi_hx_pb')
draw_2d_hist_maps(h4_num=h4_mid_cs_num_pb, h4_den=h4_mid_cs_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='#phi_{CS} (rad)', save_name='eff_mid/2d_cos_phi_cs_pb')
draw_2d_hist_maps(h4_num=h4_mid_ep_num_pb, h4_den=h4_mid_ep_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='#phi_{EP} (rad)', save_name='eff_mid/2d_cos_phi_ep_pb')

draw_2d_hist_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_all_y/2d_cos_phi_lab_pb')
draw_2d_hist_maps(h4_num=h4_all_y_hx_num_pb, h4_den=h4_all_y_hx_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{HX}', y_title='#phi_{HX} (rad)', save_name='eff_all_y/2d_cos_phi_hx_pb')
draw_2d_hist_maps(h4_num=h4_all_y_cs_num_pb, h4_den=h4_all_y_cs_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{CS}', y_title='#phi_{CS} (rad)', save_name='eff_all_y/2d_cos_phi_cs_pb')
draw_2d_hist_maps(h4_num=h4_all_y_ep_num_pb, h4_den=h4_all_y_ep_den_pb, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{EP}', y_title='#phi_{EP} (rad)', save_name='eff_all_y/2d_cos_phi_ep_pb')


# according to different pT bins
draw_cos_phi_maps(h4_num=h4_fwd_lab_num_pb, h4_den=h4_fwd_lab_den_pb, pt_bins=fwd_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_fwd/2d_cos_phi_lab_pb')
draw_cos_phi_maps(h4_num=h4_fwd_hx_num_pb, h4_den=h4_fwd_hx_den_pb, pt_bins=fwd_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='#phi_{hx} (rad)', save_name='eff_fwd/2d_cos_phi_hx_pb')
draw_cos_phi_maps(h4_num=h4_fwd_cs_num_pb, h4_den=h4_fwd_cs_den_pb, pt_bins=fwd_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='#phi_{cs} (rad)', save_name='eff_fwd/2d_cos_phi_cs_pb')
draw_cos_phi_maps(h4_num=h4_fwd_ep_num_pb, h4_den=h4_fwd_ep_den_pb, pt_bins=fwd_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='#phi_{ep} (rad)', save_name='eff_fwd/2d_cos_phi_ep_pb')

draw_cos_phi_maps(h4_num=h4_mid_lab_num_pb, h4_den=h4_mid_lab_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_mid/2d_cos_phi_lab_pb')
draw_cos_phi_maps(h4_num=h4_mid_hx_num_pb, h4_den=h4_mid_hx_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='#phi_{hx} (rad)', save_name='eff_mid/2d_cos_phi_hx_pb')
draw_cos_phi_maps(h4_num=h4_mid_cs_num_pb, h4_den=h4_mid_cs_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='#phi_{cs} (rad)', save_name='eff_mid/2d_cos_phi_cs_pb')
draw_cos_phi_maps(h4_num=h4_mid_ep_num_pb, h4_den=h4_mid_ep_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='#phi_{ep} (rad)', save_name='eff_mid/2d_cos_phi_ep_pb')

draw_cos_phi_maps(h4_num=h4_all_y_lab_num_pb, h4_den=h4_all_y_lab_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{Lab}', y_title='#phi_{Lab} (rad)', save_name='eff_all_y/2d_cos_phi_lab_pb')
draw_cos_phi_maps(h4_num=h4_all_y_hx_num_pb, h4_den=h4_all_y_hx_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{hx}', y_title='#phi_{hx} (rad)', save_name='eff_all_y/2d_cos_phi_hx_pb')
draw_cos_phi_maps(h4_num=h4_all_y_cs_num_pb, h4_den=h4_all_y_cs_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{cs}', y_title='#phi_{cs} (rad)', save_name='eff_all_y/2d_cos_phi_cs_pb')
draw_cos_phi_maps(h4_num=h4_all_y_ep_num_pb, h4_den=h4_all_y_ep_den_pb, pt_bins=mid_pt_bins, axis_pt=Axis.pt, axis1=Axis.phi, axis2=Axis.cos, legend1='', legend2='', x_title='cos#theta_{ep}', y_title='#phi_{ep} (rad)', save_name='eff_all_y/2d_cos_phi_ep_pb')