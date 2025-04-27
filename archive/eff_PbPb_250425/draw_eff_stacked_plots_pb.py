from module_draw_stacked_fcns import *
import ROOT
from ROOT import TFile, TCanvas, TLegend, TH1D, TPad
from enum import IntEnum


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
# cos
draw_comparison_plot(hists1=(h4_fwd_lab_num_pb, h4_fwd_lab_den_pb), hists2=(h4_mid_lab_num_pb, h4_mid_lab_den_pb), hists3=(h4_all_y_lab_num_pb, h4_all_y_lab_den_pb), axis=Axis.cos, title='Efficiency', x_title='cos#theta_{Lab}', y_title='Correction', save_name='stacked_plots/cos_lab_pb')
draw_comparison_plot(hists1=(h4_fwd_hx_num_pb, h4_fwd_hx_den_pb), hists2=(h4_mid_hx_num_pb, h4_mid_hx_den_pb), hists3=(h4_all_y_hx_num_pb, h4_all_y_hx_den_pb), axis=Axis.cos, title='Efficiency', x_title='cos#theta_{HX}', y_title='Correction', save_name='stacked_plots/cos_hx_pb')
draw_comparison_plot(hists1=(h4_fwd_cs_num_pb, h4_fwd_cs_den_pb), hists2=(h4_mid_cs_num_pb, h4_mid_cs_den_pb), hists3=(h4_all_y_cs_num_pb, h4_all_y_cs_den_pb), axis=Axis.cos, title='Efficiency', x_title='cos#theta_{CS}', y_title='Correction', save_name='stacked_plots/cos_cs_pb')
draw_comparison_plot(hists1=(h4_fwd_ep_num_pb, h4_fwd_ep_den_pb), hists2=(h4_mid_ep_num_pb, h4_mid_ep_den_pb), hists3=(h4_all_y_ep_num_pb, h4_all_y_ep_den_pb), axis=Axis.cos, title='Efficiency', x_title='cos#theta_{EP}', y_title='Correction', save_name='stacked_plots/cos_ep_pb')

# phi
draw_comparison_plot(hists1=(h4_fwd_lab_num_pb, h4_fwd_lab_den_pb), hists2=(h4_mid_lab_num_pb, h4_mid_lab_den_pb), hists3=(h4_all_y_lab_num_pb, h4_all_y_lab_den_pb), axis=Axis.phi, title='Efficiency', x_title='#phi_{Lab}', y_title='Correction', save_name='stacked_plots/phi_lab_pb')
draw_comparison_plot(hists1=(h4_fwd_hx_num_pb, h4_fwd_hx_den_pb), hists2=(h4_mid_hx_num_pb, h4_mid_hx_den_pb), hists3=(h4_all_y_hx_num_pb, h4_all_y_hx_den_pb), axis=Axis.phi, title='Efficiency', x_title='#phi_{HX}', y_title='Correction', save_name='stacked_plots/phi_hx_pb')
draw_comparison_plot(hists1=(h4_fwd_cs_num_pb, h4_fwd_cs_den_pb), hists2=(h4_mid_cs_num_pb, h4_mid_cs_den_pb), hists3=(h4_all_y_cs_num_pb, h4_all_y_cs_den_pb), axis=Axis.phi, title='Efficiency', x_title='#phi_{CS}', y_title='Correction', save_name='stacked_plots/phi_cs_pb')
draw_comparison_plot(hists1=(h4_fwd_ep_num_pb, h4_fwd_ep_den_pb), hists2=(h4_mid_ep_num_pb, h4_mid_ep_den_pb), hists3=(h4_all_y_ep_num_pb, h4_all_y_ep_den_pb), axis=Axis.phi, title='Efficiency', x_title='#phi_{EP}', y_title='Correction', save_name='stacked_plots/phi_ep_pb')

# pT
draw_comparison_plot(hists1=(h4_fwd_lab_num_pb, h4_fwd_lab_den_pb), hists2=(h4_mid_lab_num_pb, h4_mid_lab_den_pb), hists3=(h4_all_y_lab_num_pb, h4_all_y_lab_den_pb), axis=Axis.pt, title='Efficiency', x_title='p_{T} (GeV/c)', y_title='Correction', save_name='stacked_plots/pt_lab_pb')

# centrality
draw_comparison_plot(hists1=(h4_fwd_lab_num_pb, h4_fwd_lab_den_pb), hists2=(h4_mid_lab_num_pb, h4_mid_lab_den_pb), hists3=(h4_all_y_lab_num_pb, h4_all_y_lab_den_pb), axis=Axis.cent, title='Efficiency', x_title='Centrality', y_title='Correction', save_name='stacked_plots/cent_lab_pb')