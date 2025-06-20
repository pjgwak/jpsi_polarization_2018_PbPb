//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 24 15:07:26 2025 by ROOT version 6.32.02
// from TTree NOMINAL/NOMINAL
// found on file: ../../ggH.root
//////////////////////////////////////////////////////////

#ifndef tree_h
#define tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "TVector3.h"

class tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          HLT_2e17_lhvloose_nod0_L12EM15VHI;
   Int_t           HTXS_Njets_pTjet25;
   Int_t           HTXS_Njets_pTjet30;
   Int_t           HTXS_Stage0_Category;
   Int_t           HTXS_Stage1_Category_pTjet25GeV;
   Int_t           HTXS_Stage1_Category_pTjet30GeV;
   Int_t           HTXS_errorMode;
   Int_t           HTXS_prodMode;
   Float_t         NOMINAL_pileup_combined_weight;
   UInt_t          NOMINAL_pileup_random_lb_number;
   UInt_t          NOMINAL_pileup_random_run_number;
   Float_t         PRW_DATASF_1down_pileup_combined_weight;
   Float_t         PRW_DATASF_1up_pileup_combined_weight;
   TLorentzVector  *boson_0_truth_p4;
   Int_t           boson_0_truth_pdgId;
   Float_t         boson_0_truth_q;
   Int_t           boson_0_truth_status;
   Int_t           channel_index;
   TLorentzVector  *dijet_p4;
   Int_t           ditau_coll_approx;
   Float_t         ditau_coll_approx_m;
   Float_t         ditau_coll_approx_x0;
   Float_t         ditau_coll_approx_x1;
   Float_t         ditau_cosalpha;
   Float_t         ditau_deta;
   Float_t         ditau_dphi;
   Float_t         ditau_dr;
   Double_t        ditau_higgspt;
   Int_t           ditau_matched;
   Float_t         ditau_matched_cosalpha;
   Float_t         ditau_matched_deta;
   Float_t         ditau_matched_dphi;
   Float_t         ditau_matched_dr;
   TLorentzVector  *ditau_matched_p4;
   Float_t         ditau_matched_qxq;
   Float_t         ditau_matched_scal_sum_pt;
   Float_t         ditau_matched_vis_cosalpha;
   Float_t         ditau_matched_vis_deta;
   Float_t         ditau_matched_vis_dphi;
   Float_t         ditau_matched_vis_dr;
   Float_t         ditau_matched_vis_mass;
   Float_t         ditau_matched_vis_scal_sum_pt;
   Float_t         ditau_matched_vis_vect_sum_pt;
   Float_t         ditau_met_centrality;
   Float_t         ditau_met_lep0_cos_dphi;
   Float_t         ditau_met_lep1_cos_dphi;
   Float_t         ditau_met_min_dphi;
   Float_t         ditau_met_sum_cos_dphi;
   Int_t           ditau_mmc_mlm_fit_status;
   Float_t         ditau_mmc_mlm_m;
   Float_t         ditau_mt_lep0_met;
   Float_t         ditau_mt_lep1_met;
   TLorentzVector  *ditau_p4;
   Float_t         ditau_qxq;
   Float_t         ditau_scal_sum_pt;
   UInt_t          event_clean_EC_TightBad;
   ULong64_t       event_number;
   UInt_t          is_dijet_centrality;
   Int_t           jet_0_b_tagged;
   Float_t         jet_0_fjvt;
   Int_t           jet_0_flavorlabel_part;
   Int_t           jet_0_is_Jvt_HS;
   Float_t         jet_0_jvt;
   TLorentzVector  *jet_0_p4;
   Float_t         jet_0_width;
   TLorentzVector  *jet_0_wztruth_p4;
   Float_t         jet_0_wztruth_pdgid;
   Int_t           jet_1_b_tagged;
   Float_t         jet_1_fjvt;
   Int_t           jet_1_flavorlabel_part;
   Int_t           jet_1_is_Jvt_HS;
   Float_t         jet_1_jvt;
   TLorentzVector  *jet_1_p4;
   Float_t         jet_1_width;
   TLorentzVector  *jet_1_wztruth_p4;
   Float_t         jet_1_wztruth_pdgid;
   Int_t           jet_2_b_tagged;
   Float_t         jet_2_fjvt;
   Int_t           jet_2_flavorlabel_part;
   Int_t           jet_2_is_Jvt_HS;
   Float_t         jet_2_jvt;
   TLorentzVector  *jet_2_p4;
   Float_t         jet_2_width;
   TLorentzVector  *jet_2_wztruth_p4;
   Float_t         jet_2_wztruth_pdgid;
   Float_t         jet_FT_EFF_Eigen_B_0_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_0_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_0_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_0_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_1_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_1_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_1_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_1_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_2_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_2_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_2_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_B_2_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_0_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_0_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_0_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_0_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_1_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_1_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_1_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_1_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_2_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_2_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_2_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_2_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_3_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_3_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_3_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_C_3_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_0_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_0_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_0_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_0_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_1_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_1_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_1_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_1_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_2_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_2_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_2_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_2_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_3_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_3_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_3_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_Eigen_Light_3_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_1up_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_from_charm_1down_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_from_charm_1down_global_ineffSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_from_charm_1up_global_effSF_MV2c10;
   Float_t         jet_FT_EFF_extrapolation_from_charm_1up_global_ineffSF_MV2c10;
   Float_t         jet_JET_JvtEfficiency_1down_central_jets_global_effSF_JVT;
   Float_t         jet_JET_JvtEfficiency_1down_central_jets_global_ineffSF_JVT;
   Float_t         jet_JET_JvtEfficiency_1up_central_jets_global_effSF_JVT;
   Float_t         jet_JET_JvtEfficiency_1up_central_jets_global_ineffSF_JVT;
   Float_t         jet_JET_fJvtEfficiency_1down_forward_jets_global_effSF_JVT;
   Float_t         jet_JET_fJvtEfficiency_1down_forward_jets_global_ineffSF_JVT;
   Float_t         jet_JET_fJvtEfficiency_1up_forward_jets_global_effSF_JVT;
   Float_t         jet_JET_fJvtEfficiency_1up_forward_jets_global_ineffSF_JVT;
   Float_t         jet_NOMINAL_central_jets_global_effSF_JVT;
   Float_t         jet_NOMINAL_central_jets_global_ineffSF_JVT;
   Float_t         jet_NOMINAL_forward_jets_global_effSF_JVT;
   Float_t         jet_NOMINAL_forward_jets_global_ineffSF_JVT;
   Float_t         jet_NOMINAL_global_effSF_MV2c10;
   Float_t         jet_NOMINAL_global_ineffSF_MV2c10;
   Int_t           leplep_fake_tf_bin;
   Float_t         lepton_eta_centrality;
   UInt_t          mc_channel_number;
   TLorentzVector  *met_hpto_p4;
   Float_t         met_more_met_et_ele;
   Float_t         met_more_met_et_jet;
   Float_t         met_more_met_et_muon;
   Float_t         met_more_met_et_pho;
   Float_t         met_more_met_et_soft;
   Float_t         met_more_met_et_tau;
   Float_t         met_more_met_phi_ele;
   Float_t         met_more_met_phi_jet;
   Float_t         met_more_met_phi_muon;
   Float_t         met_more_met_phi_pho;
   Float_t         met_more_met_phi_soft;
   Float_t         met_more_met_phi_tau;
   Float_t         met_more_met_sumet_ele;
   Float_t         met_more_met_sumet_jet;
   Float_t         met_more_met_sumet_muon;
   Float_t         met_more_met_sumet_pho;
   Float_t         met_more_met_sumet_soft;
   Float_t         met_more_met_sumet_tau;
   TLorentzVector  *met_p4;
   Float_t         met_sign_met_over_sqrt_ht;
   Float_t         met_sign_met_over_sqrt_sumet;
   Float_t         met_sign_met_rho;
   Float_t         met_sign_met_rho_ttdir;
   Float_t         met_sign_met_sig_directional;
   Float_t         met_sign_met_sig_directional_ttdir;
   Float_t         met_sign_met_significance;
   Float_t         met_sign_met_significance_ttdir;
   Float_t         met_sign_met_valL;
   Float_t         met_sign_met_valL_ttdir;
   Float_t         met_sign_met_varT;
   Float_t         met_sign_met_varT_ttdir;
   Float_t         met_sumet;
   TLorentzVector  *met_truth_p4;
   Float_t         met_truth_sumet;
   Double_t        mva_random_number;
   Float_t         n_actual_int;
   Float_t         n_actual_int_cor;
   Float_t         n_avg_int;
   Float_t         n_avg_int_cor;
   Int_t           n_bjets;
   Int_t           n_electrons;
   Int_t           n_jets;
   Int_t           n_jets_30;
   Int_t           n_jets_40;
   Int_t           n_jets_central;
   Int_t           n_jets_central_30;
   Int_t           n_jets_central_40;
   Int_t           n_jets_forward;
   Int_t           n_jets_forward_30;
   Int_t           n_jets_forward_40;
   Int_t           n_jets_l1j25;
   Int_t           n_jets_mc_hs;
   Int_t           n_muons;
   Int_t           n_photons;
   Int_t           n_pvx;
   Int_t           n_taus;
   Int_t           n_taus_loose;
   Int_t           n_taus_medium;
   Int_t           n_taus_tight;
   Int_t           n_taus_veryloose;
   UInt_t          n_truth_gluon_jets;
   UInt_t          n_truth_jets;
   UInt_t          n_truth_jets_pt20_eta45;
   UInt_t          n_truth_quark_jets;
   Int_t           n_vx;
   Int_t           primary_vertex;
   TVector3        *primary_vertex_v;
   Float_t         pt_total;
   UInt_t          run_number;
   Float_t         scalar_sum_pt;
   UInt_t          tau_0;
   Float_t         tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk;
   Float_t         tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_0_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_0_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_0_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_0_NOMINAL_EleEffSF_offline_RecoTrk;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_0_NOMINAL_MuEffSF_IsoFCLoose;
   Float_t         tau_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_0_NOMINAL_MuEffSF_Reco_QualMedium;
   Float_t         tau_0_cluster_eta;
   Float_t         tau_0_cluster_eta_be2;
   Float_t         tau_0_cluster_pt;
   UInt_t          tau_0_electron_trig_HLT_2e12_lhloose_L12EM10VH;
   UInt_t          tau_0_electron_trig_HLT_2e17_lhvloose_nod0;
   UInt_t          tau_0_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI;
   UInt_t          tau_0_electron_trig_HLT_e120_lhloose;
   UInt_t          tau_0_electron_trig_HLT_e140_lhloose_nod0;
   UInt_t          tau_0_electron_trig_HLT_e24_lhmedium_L1EM20VH;
   UInt_t          tau_0_electron_trig_HLT_e26_lhtight_nod0_ivarloose;
   UInt_t          tau_0_electron_trig_HLT_e60_lhmedium;
   UInt_t          tau_0_electron_trig_HLT_e60_lhmedium_nod0;
   UInt_t          tau_0_electron_trig_trigger_matched;
   UInt_t          tau_0_emu_trig_HLT_e17_lhloose_mu14;
   UInt_t          tau_0_emu_trig_HLT_e17_lhloose_nod0_mu14;
   UInt_t          tau_0_emu_trig_trigger_matched;
   Int_t           tau_0_id_bad;
   Int_t           tau_0_id_charge;
   Int_t           tau_0_id_loose;
   Int_t           tau_0_id_medium;
   Int_t           tau_0_id_tight;
   Int_t           tau_0_id_veryloose;
   UInt_t          tau_0_iso_FCHighPtCaloOnly;
   UInt_t          tau_0_iso_FCLoose;
   UInt_t          tau_0_iso_FCTight;
   UInt_t          tau_0_iso_FCTightTrackOnly;
   UInt_t          tau_0_iso_FixedCutHighPtTrackOnly;
   UInt_t          tau_0_iso_Gradient;
   Float_t         tau_0_iso_ptvarcone20;
   Float_t         tau_0_iso_ptvarcone30;
   UInt_t          tau_0_matched;
   Int_t           tau_0_matched_classifierParticleOrigin;
   Int_t           tau_0_matched_classifierParticleType;
   Int_t           tau_0_matched_isHad;
   UInt_t          tau_0_matched_leptonic_tau;
   Int_t           tau_0_matched_leptonic_tau_classifierParticleOrigin;
   Int_t           tau_0_matched_leptonic_tau_classifierParticleType;
   TLorentzVector  *tau_0_matched_leptonic_tau_invis_p4;
   Int_t           tau_0_matched_leptonic_tau_mother_pdgId;
   Int_t           tau_0_matched_leptonic_tau_mother_status;
   Int_t           tau_0_matched_leptonic_tau_origin;
   TLorentzVector  *tau_0_matched_leptonic_tau_p4;
   Int_t           tau_0_matched_leptonic_tau_pdgId;
   Float_t         tau_0_matched_leptonic_tau_pz;
   Float_t         tau_0_matched_leptonic_tau_q;
   Int_t           tau_0_matched_leptonic_tau_status;
   Int_t           tau_0_matched_leptonic_tau_type;
   TLorentzVector  *tau_0_matched_leptonic_tau_vis_p4;
   Int_t           tau_0_matched_mother_pdgId;
   Int_t           tau_0_matched_mother_status;
   TLorentzVector  *tau_0_matched_p4;
   Int_t           tau_0_matched_pdgId;
   Float_t         tau_0_matched_q;
   Int_t           tau_0_matched_status;
   Int_t           tau_0_matched_type;
   Int_t           tau_0_muonAuthor;
   Int_t           tau_0_muonType;
   UInt_t          tau_0_muon_trig_HLT_mu18_mu8noL1;
   UInt_t          tau_0_muon_trig_HLT_mu20_iloose_L1MU15;
   UInt_t          tau_0_muon_trig_HLT_mu20_mu8noL1;
   UInt_t          tau_0_muon_trig_HLT_mu22_mu8noL1;
   UInt_t          tau_0_muon_trig_HLT_mu26_ivarmedium;
   UInt_t          tau_0_muon_trig_HLT_mu50;
   UInt_t          tau_0_muon_trig_trigger_matched;
   Int_t           tau_0_origin;
   TLorentzVector  *tau_0_p4;
   Float_t         tau_0_q;
   Float_t         tau_0_trk_d0;
   Float_t         tau_0_trk_d0_sig;
   Float_t         tau_0_trk_pt;
   Float_t         tau_0_trk_pt_error;
   Float_t         tau_0_trk_pvx_z0;
   Float_t         tau_0_trk_pvx_z0_sig;
   Float_t         tau_0_trk_pvx_z0_sintheta;
   Int_t           tau_0_trk_vx;
   TVector3        *tau_0_trk_vx_v;
   Float_t         tau_0_trk_z0;
   Float_t         tau_0_trk_z0_sig;
   Float_t         tau_0_trk_z0_sintheta;
   Int_t           tau_0_type;
   UInt_t          tau_1;
   Float_t         tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk;
   Float_t         tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_1_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;
   Float_t         tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;
   Float_t         tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;
   Float_t         tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;
   Float_t         tau_1_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;
   Float_t         tau_1_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13;
   Float_t         tau_1_NOMINAL_EleEffSF_offline_RecoTrk;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;
   Float_t         tau_1_NOMINAL_MuEffSF_IsoFCLoose;
   Float_t         tau_1_NOMINAL_MuEffSF_IsoFCTightTrackOnly;
   Float_t         tau_1_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly;
   Float_t         tau_1_NOMINAL_MuEffSF_Reco_QualMedium;
   Float_t         tau_1_cluster_eta;
   Float_t         tau_1_cluster_eta_be2;
   Float_t         tau_1_cluster_pt;
   UInt_t          tau_1_electron_trig_HLT_2e12_lhloose_L12EM10VH;
   UInt_t          tau_1_electron_trig_HLT_2e17_lhvloose_nod0;
   UInt_t          tau_1_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI;
   UInt_t          tau_1_electron_trig_HLT_e120_lhloose;
   UInt_t          tau_1_electron_trig_HLT_e140_lhloose_nod0;
   UInt_t          tau_1_electron_trig_HLT_e24_lhmedium_L1EM20VH;
   UInt_t          tau_1_electron_trig_HLT_e26_lhtight_nod0_ivarloose;
   UInt_t          tau_1_electron_trig_HLT_e60_lhmedium;
   UInt_t          tau_1_electron_trig_HLT_e60_lhmedium_nod0;
   UInt_t          tau_1_electron_trig_trigger_matched;
   UInt_t          tau_1_emu_trig_HLT_e17_lhloose_mu14;
   UInt_t          tau_1_emu_trig_HLT_e17_lhloose_nod0_mu14;
   UInt_t          tau_1_emu_trig_trigger_matched;
   Int_t           tau_1_id_bad;
   Int_t           tau_1_id_charge;
   Int_t           tau_1_id_loose;
   Int_t           tau_1_id_medium;
   Int_t           tau_1_id_tight;
   Int_t           tau_1_id_veryloose;
   UInt_t          tau_1_iso_FCHighPtCaloOnly;
   UInt_t          tau_1_iso_FCLoose;
   UInt_t          tau_1_iso_FCTight;
   UInt_t          tau_1_iso_FCTightTrackOnly;
   UInt_t          tau_1_iso_FixedCutHighPtTrackOnly;
   UInt_t          tau_1_iso_Gradient;
   Float_t         tau_1_iso_ptvarcone20;
   Float_t         tau_1_iso_ptvarcone30;
   UInt_t          tau_1_matched;
   Int_t           tau_1_matched_classifierParticleOrigin;
   Int_t           tau_1_matched_classifierParticleType;
   Int_t           tau_1_matched_isHad;
   TLorentzVector  *tau_1_matched_leptonic_tau_invis_p4;
   TLorentzVector  *tau_1_matched_leptonic_tau_p4;
   TLorentzVector  *tau_1_matched_leptonic_tau_vis_p4;
   Int_t           tau_1_matched_mother_pdgId;
   Int_t           tau_1_matched_mother_status;
   TLorentzVector  *tau_1_matched_p4;
   Int_t           tau_1_matched_pdgId;
   Float_t         tau_1_matched_q;
   Int_t           tau_1_matched_status;
   Int_t           tau_1_muonAuthor;
   Int_t           tau_1_muonType;
   UInt_t          tau_1_muon_trig_HLT_mu18_mu8noL1;
   UInt_t          tau_1_muon_trig_HLT_mu20_iloose_L1MU15;
   UInt_t          tau_1_muon_trig_HLT_mu20_mu8noL1;
   UInt_t          tau_1_muon_trig_HLT_mu22_mu8noL1;
   UInt_t          tau_1_muon_trig_HLT_mu26_ivarmedium;
   UInt_t          tau_1_muon_trig_HLT_mu50;
   UInt_t          tau_1_muon_trig_trigger_matched;
   Int_t           tau_1_origin;
   TLorentzVector  *tau_1_p4;
   Float_t         tau_1_q;
   Float_t         tau_1_trk_d0;
   Float_t         tau_1_trk_d0_sig;
   Float_t         tau_1_trk_pt;
   Float_t         tau_1_trk_pt_error;
   Float_t         tau_1_trk_pvx_z0;
   Float_t         tau_1_trk_pvx_z0_sig;
   Float_t         tau_1_trk_pvx_z0_sintheta;
   Int_t           tau_1_trk_vx;
   TVector3        *tau_1_trk_vx_v;
   Float_t         tau_1_trk_z0;
   Float_t         tau_1_trk_z0_sig;
   Float_t         tau_1_trk_z0_sintheta;
   Int_t           tau_1_type;
   Float_t         tau_eta_centrality;
   Float_t         theory_weights_alphaS_down;
   Float_t         theory_weights_alphaS_up;
   Float_t         theory_weights_nominal;
   Float_t         theory_weights_pdf_signal_weight_0;
   Float_t         theory_weights_pdf_signal_weight_1;
   Float_t         theory_weights_pdf_signal_weight_10;
   Float_t         theory_weights_pdf_signal_weight_11;
   Float_t         theory_weights_pdf_signal_weight_12;
   Float_t         theory_weights_pdf_signal_weight_13;
   Float_t         theory_weights_pdf_signal_weight_14;
   Float_t         theory_weights_pdf_signal_weight_15;
   Float_t         theory_weights_pdf_signal_weight_16;
   Float_t         theory_weights_pdf_signal_weight_17;
   Float_t         theory_weights_pdf_signal_weight_18;
   Float_t         theory_weights_pdf_signal_weight_19;
   Float_t         theory_weights_pdf_signal_weight_2;
   Float_t         theory_weights_pdf_signal_weight_20;
   Float_t         theory_weights_pdf_signal_weight_21;
   Float_t         theory_weights_pdf_signal_weight_22;
   Float_t         theory_weights_pdf_signal_weight_23;
   Float_t         theory_weights_pdf_signal_weight_24;
   Float_t         theory_weights_pdf_signal_weight_25;
   Float_t         theory_weights_pdf_signal_weight_26;
   Float_t         theory_weights_pdf_signal_weight_27;
   Float_t         theory_weights_pdf_signal_weight_28;
   Float_t         theory_weights_pdf_signal_weight_29;
   Float_t         theory_weights_pdf_signal_weight_3;
   Float_t         theory_weights_pdf_signal_weight_4;
   Float_t         theory_weights_pdf_signal_weight_5;
   Float_t         theory_weights_pdf_signal_weight_6;
   Float_t         theory_weights_pdf_signal_weight_7;
   Float_t         theory_weights_pdf_signal_weight_8;
   Float_t         theory_weights_pdf_signal_weight_9;
   Float_t         theory_weights_qcd_weight_0;
   Float_t         theory_weights_qcd_weight_1;
   Float_t         theory_weights_qcd_weight_2;
   Float_t         theory_weights_qcd_weight_3;
   Float_t         theory_weights_qcd_weight_4;
   Float_t         theory_weights_qcd_weight_5;
   Float_t         theory_weights_qcd_weight_6;
   Float_t         theory_weights_qcd_weight_7;
   Float_t         theory_weights_qcd_weight_8;
   UInt_t          truth_passedVBFFilter;
   Double_t        weight_mc;

   // List of branches
   TBranch        *b_HLT_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_HTXS_Njets_pTjet25;   //!
   TBranch        *b_HTXS_Njets_pTjet30;   //!
   TBranch        *b_HTXS_Stage0_Category;   //!
   TBranch        *b_HTXS_Stage1_Category_pTjet25GeV;   //!
   TBranch        *b_HTXS_Stage1_Category_pTjet30GeV;   //!
   TBranch        *b_HTXS_errorMode;   //!
   TBranch        *b_HTXS_prodMode;   //!
   TBranch        *b_NOMINAL_pileup_combined_weight;   //!
   TBranch        *b_NOMINAL_pileup_random_lb_number;   //!
   TBranch        *b_NOMINAL_pileup_random_run_number;   //!
   TBranch        *b_PRW_DATASF_1down_pileup_combined_weight;   //!
   TBranch        *b_PRW_DATASF_1up_pileup_combined_weight;   //!
   TBranch        *b_boson_0_truth_p4;   //!
   TBranch        *b_boson_0_truth_pdgId;   //!
   TBranch        *b_boson_0_truth_q;   //!
   TBranch        *b_boson_0_truth_status;   //!
   TBranch        *b_channel_index;   //!
   TBranch        *b_dijet_p4;   //!
   TBranch        *b_ditau_coll_approx;   //!
   TBranch        *b_ditau_coll_approx_m;   //!
   TBranch        *b_ditau_coll_approx_x0;   //!
   TBranch        *b_ditau_coll_approx_x1;   //!
   TBranch        *b_ditau_cosalpha;   //!
   TBranch        *b_ditau_deta;   //!
   TBranch        *b_ditau_dphi;   //!
   TBranch        *b_ditau_dr;   //!
   TBranch        *b_ditau_higgspt;   //!
   TBranch        *b_ditau_matched;   //!
   TBranch        *b_ditau_matched_cosalpha;   //!
   TBranch        *b_ditau_matched_deta;   //!
   TBranch        *b_ditau_matched_dphi;   //!
   TBranch        *b_ditau_matched_dr;   //!
   TBranch        *b_ditau_matched_p4;   //!
   TBranch        *b_ditau_matched_qxq;   //!
   TBranch        *b_ditau_matched_scal_sum_pt;   //!
   TBranch        *b_ditau_matched_vis_cosalpha;   //!
   TBranch        *b_ditau_matched_vis_deta;   //!
   TBranch        *b_ditau_matched_vis_dphi;   //!
   TBranch        *b_ditau_matched_vis_dr;   //!
   TBranch        *b_ditau_matched_vis_mass;   //!
   TBranch        *b_ditau_matched_vis_scal_sum_pt;   //!
   TBranch        *b_ditau_matched_vis_vect_sum_pt;   //!
   TBranch        *b_ditau_met_centrality;   //!
   TBranch        *b_ditau_met_lep0_cos_dphi;   //!
   TBranch        *b_ditau_met_lep1_cos_dphi;   //!
   TBranch        *b_ditau_met_min_dphi;   //!
   TBranch        *b_ditau_met_sum_cos_dphi;   //!
   TBranch        *b_ditau_mmc_mlm_fit_status;   //!
   TBranch        *b_ditau_mmc_mlm_m;   //!
   TBranch        *b_ditau_mt_lep0_met;   //!
   TBranch        *b_ditau_mt_lep1_met;   //!
   TBranch        *b_ditau_p4;   //!
   TBranch        *b_ditau_qxq;   //!
   TBranch        *b_ditau_scal_sum_pt;   //!
   TBranch        *b_event_clean_EC_TightBad;   //!
   TBranch        *b_event_number;   //!
   TBranch        *b_is_dijet_centrality;   //!
   TBranch        *b_jet_0_b_tagged;   //!
   TBranch        *b_jet_0_fjvt;   //!
   TBranch        *b_jet_0_flavorlabel_part;   //!
   TBranch        *b_jet_0_is_Jvt_HS;   //!
   TBranch        *b_jet_0_jvt;   //!
   TBranch        *b_jet_0_p4;   //!
   TBranch        *b_jet_0_width;   //!
   TBranch        *b_jet_0_wztruth_p4;   //!
   TBranch        *b_jet_0_wztruth_pdgid;   //!
   TBranch        *b_jet_1_b_tagged;   //!
   TBranch        *b_jet_1_fjvt;   //!
   TBranch        *b_jet_1_flavorlabel_part;   //!
   TBranch        *b_jet_1_is_Jvt_HS;   //!
   TBranch        *b_jet_1_jvt;   //!
   TBranch        *b_jet_1_p4;   //!
   TBranch        *b_jet_1_width;   //!
   TBranch        *b_jet_1_wztruth_p4;   //!
   TBranch        *b_jet_1_wztruth_pdgid;   //!
   TBranch        *b_jet_2_b_tagged;   //!
   TBranch        *b_jet_2_fjvt;   //!
   TBranch        *b_jet_2_flavorlabel_part;   //!
   TBranch        *b_jet_2_is_Jvt_HS;   //!
   TBranch        *b_jet_2_jvt;   //!
   TBranch        *b_jet_2_p4;   //!
   TBranch        *b_jet_2_width;   //!
   TBranch        *b_jet_2_wztruth_p4;   //!
   TBranch        *b_jet_2_wztruth_pdgid;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_0_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_0_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_0_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_0_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_1_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_1_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_1_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_1_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_2_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_2_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_2_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_B_2_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_0_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_0_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_0_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_0_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_1_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_1_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_1_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_1_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_2_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_2_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_2_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_2_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_3_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_3_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_3_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_C_3_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_0_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_0_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_0_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_0_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_1_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_1_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_1_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_1_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_2_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_2_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_2_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_2_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_3_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_3_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_3_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_Eigen_Light_3_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_from_charm_1down_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_from_charm_1down_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_from_charm_1up_global_effSF_MV2c10;   //!
   TBranch        *b_jet_FT_EFF_extrapolation_from_charm_1up_global_ineffSF_MV2c10;   //!
   TBranch        *b_jet_JET_JvtEfficiency_1down_central_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_JET_JvtEfficiency_1down_central_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_JET_JvtEfficiency_1up_central_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_JET_JvtEfficiency_1up_central_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_JET_fJvtEfficiency_1down_forward_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_JET_fJvtEfficiency_1down_forward_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_JET_fJvtEfficiency_1up_forward_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_JET_fJvtEfficiency_1up_forward_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_NOMINAL_central_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_NOMINAL_central_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_NOMINAL_forward_jets_global_effSF_JVT;   //!
   TBranch        *b_jet_NOMINAL_forward_jets_global_ineffSF_JVT;   //!
   TBranch        *b_jet_NOMINAL_global_effSF_MV2c10;   //!
   TBranch        *b_jet_NOMINAL_global_ineffSF_MV2c10;   //!
   TBranch        *b_leplep_fake_tf_bin;   //!
   TBranch        *b_lepton_eta_centrality;   //!
   TBranch        *b_mc_channel_number;   //!
   TBranch        *b_met_hpto_p4;   //!
   TBranch        *b_met_more_met_et_ele;   //!
   TBranch        *b_met_more_met_et_jet;   //!
   TBranch        *b_met_more_met_et_muon;   //!
   TBranch        *b_met_more_met_et_pho;   //!
   TBranch        *b_met_more_met_et_soft;   //!
   TBranch        *b_met_more_met_et_tau;   //!
   TBranch        *b_met_more_met_phi_ele;   //!
   TBranch        *b_met_more_met_phi_jet;   //!
   TBranch        *b_met_more_met_phi_muon;   //!
   TBranch        *b_met_more_met_phi_pho;   //!
   TBranch        *b_met_more_met_phi_soft;   //!
   TBranch        *b_met_more_met_phi_tau;   //!
   TBranch        *b_met_more_met_sumet_ele;   //!
   TBranch        *b_met_more_met_sumet_jet;   //!
   TBranch        *b_met_more_met_sumet_muon;   //!
   TBranch        *b_met_more_met_sumet_pho;   //!
   TBranch        *b_met_more_met_sumet_soft;   //!
   TBranch        *b_met_more_met_sumet_tau;   //!
   TBranch        *b_met_p4;   //!
   TBranch        *b_met_sign_met_over_sqrt_ht;   //!
   TBranch        *b_met_sign_met_over_sqrt_sumet;   //!
   TBranch        *b_met_sign_met_rho;   //!
   TBranch        *b_met_sign_met_rho_ttdir;   //!
   TBranch        *b_met_sign_met_sig_directional;   //!
   TBranch        *b_met_sign_met_sig_directional_ttdir;   //!
   TBranch        *b_met_sign_met_significance;   //!
   TBranch        *b_met_sign_met_significance_ttdir;   //!
   TBranch        *b_met_sign_met_valL;   //!
   TBranch        *b_met_sign_met_valL_ttdir;   //!
   TBranch        *b_met_sign_met_varT;   //!
   TBranch        *b_met_sign_met_varT_ttdir;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_truth_p4;   //!
   TBranch        *b_met_truth_sumet;   //!
   TBranch        *b_mva_random_number;   //!
   TBranch        *b_n_actual_int;   //!
   TBranch        *b_n_actual_int_cor;   //!
   TBranch        *b_n_avg_int;   //!
   TBranch        *b_n_avg_int_cor;   //!
   TBranch        *b_n_bjets;   //!
   TBranch        *b_n_electrons;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_n_jets_30;   //!
   TBranch        *b_n_jets_40;   //!
   TBranch        *b_n_jets_central;   //!
   TBranch        *b_n_jets_central_30;   //!
   TBranch        *b_n_jets_central_40;   //!
   TBranch        *b_n_jets_forward;   //!
   TBranch        *b_n_jets_forward_30;   //!
   TBranch        *b_n_jets_forward_40;   //!
   TBranch        *b_n_jets_l1j25;   //!
   TBranch        *b_n_jets_mc_hs;   //!
   TBranch        *b_n_muons;   //!
   TBranch        *b_n_photons;   //!
   TBranch        *b_n_pvx;   //!
   TBranch        *b_n_taus;   //!
   TBranch        *b_n_taus_loose;   //!
   TBranch        *b_n_taus_medium;   //!
   TBranch        *b_n_taus_tight;   //!
   TBranch        *b_n_taus_veryloose;   //!
   TBranch        *b_n_truth_gluon_jets;   //!
   TBranch        *b_n_truth_jets;   //!
   TBranch        *b_n_truth_jets_pt20_eta45;   //!
   TBranch        *b_n_truth_quark_jets;   //!
   TBranch        *b_n_vx;   //!
   TBranch        *b_primary_vertex;   //!
   TBranch        *b_primary_vertex_v;   //!
   TBranch        *b_pt_total;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_scalar_sum_pt;   //!
   TBranch        *b_tau_0;   //!
   TBranch        *b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_0_NOMINAL_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_NOMINAL_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_0_cluster_eta;   //!
   TBranch        *b_tau_0_cluster_eta_be2;   //!
   TBranch        *b_tau_0_cluster_pt;   //!
   TBranch        *b_tau_0_electron_trig_HLT_2e12_lhloose_L12EM10VH;   //!
   TBranch        *b_tau_0_electron_trig_HLT_2e17_lhvloose_nod0;   //!
   TBranch        *b_tau_0_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e120_lhloose;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e60_lhmedium;   //!
   TBranch        *b_tau_0_electron_trig_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_tau_0_electron_trig_trigger_matched;   //!
   TBranch        *b_tau_0_emu_trig_HLT_e17_lhloose_mu14;   //!
   TBranch        *b_tau_0_emu_trig_HLT_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_tau_0_emu_trig_trigger_matched;   //!
   TBranch        *b_tau_0_id_bad;   //!
   TBranch        *b_tau_0_id_charge;   //!
   TBranch        *b_tau_0_id_loose;   //!
   TBranch        *b_tau_0_id_medium;   //!
   TBranch        *b_tau_0_id_tight;   //!
   TBranch        *b_tau_0_id_veryloose;   //!
   TBranch        *b_tau_0_iso_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_0_iso_FCLoose;   //!
   TBranch        *b_tau_0_iso_FCTight;   //!
   TBranch        *b_tau_0_iso_FCTightTrackOnly;   //!
   TBranch        *b_tau_0_iso_FixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_0_iso_Gradient;   //!
   TBranch        *b_tau_0_iso_ptvarcone20;   //!
   TBranch        *b_tau_0_iso_ptvarcone30;   //!
   TBranch        *b_tau_0_matched;   //!
   TBranch        *b_tau_0_matched_classifierParticleOrigin;   //!
   TBranch        *b_tau_0_matched_classifierParticleType;   //!
   TBranch        *b_tau_0_matched_isHad;   //!
   TBranch        *b_tau_0_matched_leptonic_tau;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_classifierParticleOrigin;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_classifierParticleType;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_invis_p4;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_mother_pdgId;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_mother_status;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_origin;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_p4;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_pdgId;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_pz;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_q;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_status;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_type;   //!
   TBranch        *b_tau_0_matched_leptonic_tau_vis_p4;   //!
   TBranch        *b_tau_0_matched_mother_pdgId;   //!
   TBranch        *b_tau_0_matched_mother_status;   //!
   TBranch        *b_tau_0_matched_p4;   //!
   TBranch        *b_tau_0_matched_pdgId;   //!
   TBranch        *b_tau_0_matched_q;   //!
   TBranch        *b_tau_0_matched_status;   //!
   TBranch        *b_tau_0_matched_type;   //!
   TBranch        *b_tau_0_muonAuthor;   //!
   TBranch        *b_tau_0_muonType;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu18_mu8noL1;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu20_mu8noL1;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu22_mu8noL1;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu26_ivarmedium;   //!
   TBranch        *b_tau_0_muon_trig_HLT_mu50;   //!
   TBranch        *b_tau_0_muon_trig_trigger_matched;   //!
   TBranch        *b_tau_0_origin;   //!
   TBranch        *b_tau_0_p4;   //!
   TBranch        *b_tau_0_q;   //!
   TBranch        *b_tau_0_trk_d0;   //!
   TBranch        *b_tau_0_trk_d0_sig;   //!
   TBranch        *b_tau_0_trk_pt;   //!
   TBranch        *b_tau_0_trk_pt_error;   //!
   TBranch        *b_tau_0_trk_pvx_z0;   //!
   TBranch        *b_tau_0_trk_pvx_z0_sig;   //!
   TBranch        *b_tau_0_trk_pvx_z0_sintheta;   //!
   TBranch        *b_tau_0_trk_vx;   //!
   TBranch        *b_tau_0_trk_vx_v;   //!
   TBranch        *b_tau_0_trk_z0;   //!
   TBranch        *b_tau_0_trk_z0_sig;   //!
   TBranch        *b_tau_0_trk_z0_sintheta;   //!
   TBranch        *b_tau_0_type;   //!
   TBranch        *b_tau_1;   //!
   TBranch        *b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13;   //!
   TBranch        *b_tau_1_NOMINAL_EleEffSF_offline_RecoTrk;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_IsoFCLoose;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_IsoFCTightTrackOnly;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_NOMINAL_MuEffSF_Reco_QualMedium;   //!
   TBranch        *b_tau_1_cluster_eta;   //!
   TBranch        *b_tau_1_cluster_eta_be2;   //!
   TBranch        *b_tau_1_cluster_pt;   //!
   TBranch        *b_tau_1_electron_trig_HLT_2e12_lhloose_L12EM10VH;   //!
   TBranch        *b_tau_1_electron_trig_HLT_2e17_lhvloose_nod0;   //!
   TBranch        *b_tau_1_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e120_lhloose;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e60_lhmedium;   //!
   TBranch        *b_tau_1_electron_trig_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_tau_1_electron_trig_trigger_matched;   //!
   TBranch        *b_tau_1_emu_trig_HLT_e17_lhloose_mu14;   //!
   TBranch        *b_tau_1_emu_trig_HLT_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_tau_1_emu_trig_trigger_matched;   //!
   TBranch        *b_tau_1_id_bad;   //!
   TBranch        *b_tau_1_id_charge;   //!
   TBranch        *b_tau_1_id_loose;   //!
   TBranch        *b_tau_1_id_medium;   //!
   TBranch        *b_tau_1_id_tight;   //!
   TBranch        *b_tau_1_id_veryloose;   //!
   TBranch        *b_tau_1_iso_FCHighPtCaloOnly;   //!
   TBranch        *b_tau_1_iso_FCLoose;   //!
   TBranch        *b_tau_1_iso_FCTight;   //!
   TBranch        *b_tau_1_iso_FCTightTrackOnly;   //!
   TBranch        *b_tau_1_iso_FixedCutHighPtTrackOnly;   //!
   TBranch        *b_tau_1_iso_Gradient;   //!
   TBranch        *b_tau_1_iso_ptvarcone20;   //!
   TBranch        *b_tau_1_iso_ptvarcone30;   //!
   TBranch        *b_tau_1_matched;   //!
   TBranch        *b_tau_1_matched_classifierParticleOrigin;   //!
   TBranch        *b_tau_1_matched_classifierParticleType;   //!
   TBranch        *b_tau_1_matched_isHad;   //!
   TBranch        *b_tau_1_matched_leptonic_tau_invis_p4;   //!
   TBranch        *b_tau_1_matched_leptonic_tau_p4;   //!
   TBranch        *b_tau_1_matched_leptonic_tau_vis_p4;   //!
   TBranch        *b_tau_1_matched_mother_pdgId;   //!
   TBranch        *b_tau_1_matched_mother_status;   //!
   TBranch        *b_tau_1_matched_p4;   //!
   TBranch        *b_tau_1_matched_pdgId;   //!
   TBranch        *b_tau_1_matched_q;   //!
   TBranch        *b_tau_1_matched_status;   //!
   TBranch        *b_tau_1_muonAuthor;   //!
   TBranch        *b_tau_1_muonType;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu18_mu8noL1;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu20_mu8noL1;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu22_mu8noL1;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu26_ivarmedium;   //!
   TBranch        *b_tau_1_muon_trig_HLT_mu50;   //!
   TBranch        *b_tau_1_muon_trig_trigger_matched;   //!
   TBranch        *b_tau_1_origin;   //!
   TBranch        *b_tau_1_p4;   //!
   TBranch        *b_tau_1_q;   //!
   TBranch        *b_tau_1_trk_d0;   //!
   TBranch        *b_tau_1_trk_d0_sig;   //!
   TBranch        *b_tau_1_trk_pt;   //!
   TBranch        *b_tau_1_trk_pt_error;   //!
   TBranch        *b_tau_1_trk_pvx_z0;   //!
   TBranch        *b_tau_1_trk_pvx_z0_sig;   //!
   TBranch        *b_tau_1_trk_pvx_z0_sintheta;   //!
   TBranch        *b_tau_1_trk_vx;   //!
   TBranch        *b_tau_1_trk_vx_v;   //!
   TBranch        *b_tau_1_trk_z0;   //!
   TBranch        *b_tau_1_trk_z0_sig;   //!
   TBranch        *b_tau_1_trk_z0_sintheta;   //!
   TBranch        *b_tau_1_type;   //!
   TBranch        *b_tau_eta_centrality;   //!
   TBranch        *b_theory_weights_alphaS_down;   //!
   TBranch        *b_theory_weights_alphaS_up;   //!
   TBranch        *b_theory_weights_nominal;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_0;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_1;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_10;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_11;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_12;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_13;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_14;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_15;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_16;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_17;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_18;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_19;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_2;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_20;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_21;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_22;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_23;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_24;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_25;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_26;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_27;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_28;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_29;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_3;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_4;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_5;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_6;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_7;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_8;   //!
   TBranch        *b_theory_weights_pdf_signal_weight_9;   //!
   TBranch        *b_theory_weights_qcd_weight_0;   //!
   TBranch        *b_theory_weights_qcd_weight_1;   //!
   TBranch        *b_theory_weights_qcd_weight_2;   //!
   TBranch        *b_theory_weights_qcd_weight_3;   //!
   TBranch        *b_theory_weights_qcd_weight_4;   //!
   TBranch        *b_theory_weights_qcd_weight_5;   //!
   TBranch        *b_theory_weights_qcd_weight_6;   //!
   TBranch        *b_theory_weights_qcd_weight_7;   //!
   TBranch        *b_theory_weights_qcd_weight_8;   //!
   TBranch        *b_truth_passedVBFFilter;   //!
   TBranch        *b_weight_mc;   //!

   tree(TTree *tree=0);
   virtual ~tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree_cxx
tree::tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../ggH.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../ggH.root");
      }
      f->GetObject("NOMINAL",tree);

   }
   Init(tree);
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   boson_0_truth_p4 = 0;
   dijet_p4 = 0;
   ditau_matched_p4 = 0;
   ditau_p4 = 0;
   jet_0_p4 = 0;
   jet_0_wztruth_p4 = 0;
   jet_1_p4 = 0;
   jet_1_wztruth_p4 = 0;
   jet_2_p4 = 0;
   jet_2_wztruth_p4 = 0;
   met_hpto_p4 = 0;
   met_p4 = 0;
   met_truth_p4 = 0;
   primary_vertex_v = 0;
   tau_0_matched_leptonic_tau_invis_p4 = 0;
   tau_0_matched_leptonic_tau_p4 = 0;
   tau_0_matched_leptonic_tau_vis_p4 = 0;
   tau_0_matched_p4 = 0;
   tau_0_p4 = 0;
   tau_0_trk_vx_v = 0;
   tau_1_matched_leptonic_tau_invis_p4 = 0;
   tau_1_matched_leptonic_tau_p4 = 0;
   tau_1_matched_leptonic_tau_vis_p4 = 0;
   tau_1_matched_p4 = 0;
   tau_1_p4 = 0;
   tau_1_trk_vx_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("HLT_2e17_lhvloose_nod0_L12EM15VHI", &HLT_2e17_lhvloose_nod0_L12EM15VHI, &b_HLT_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("HTXS_Njets_pTjet25", &HTXS_Njets_pTjet25, &b_HTXS_Njets_pTjet25);
   fChain->SetBranchAddress("HTXS_Njets_pTjet30", &HTXS_Njets_pTjet30, &b_HTXS_Njets_pTjet30);
   fChain->SetBranchAddress("HTXS_Stage0_Category", &HTXS_Stage0_Category, &b_HTXS_Stage0_Category);
   fChain->SetBranchAddress("HTXS_Stage1_Category_pTjet25GeV", &HTXS_Stage1_Category_pTjet25GeV, &b_HTXS_Stage1_Category_pTjet25GeV);
   fChain->SetBranchAddress("HTXS_Stage1_Category_pTjet30GeV", &HTXS_Stage1_Category_pTjet30GeV, &b_HTXS_Stage1_Category_pTjet30GeV);
   fChain->SetBranchAddress("HTXS_errorMode", &HTXS_errorMode, &b_HTXS_errorMode);
   fChain->SetBranchAddress("HTXS_prodMode", &HTXS_prodMode, &b_HTXS_prodMode);
   fChain->SetBranchAddress("NOMINAL_pileup_combined_weight", &NOMINAL_pileup_combined_weight, &b_NOMINAL_pileup_combined_weight);
   fChain->SetBranchAddress("NOMINAL_pileup_random_lb_number", &NOMINAL_pileup_random_lb_number, &b_NOMINAL_pileup_random_lb_number);
   fChain->SetBranchAddress("NOMINAL_pileup_random_run_number", &NOMINAL_pileup_random_run_number, &b_NOMINAL_pileup_random_run_number);
   fChain->SetBranchAddress("PRW_DATASF_1down_pileup_combined_weight", &PRW_DATASF_1down_pileup_combined_weight, &b_PRW_DATASF_1down_pileup_combined_weight);
   fChain->SetBranchAddress("PRW_DATASF_1up_pileup_combined_weight", &PRW_DATASF_1up_pileup_combined_weight, &b_PRW_DATASF_1up_pileup_combined_weight);
   fChain->SetBranchAddress("boson_0_truth_p4", &boson_0_truth_p4, &b_boson_0_truth_p4);
   fChain->SetBranchAddress("boson_0_truth_pdgId", &boson_0_truth_pdgId, &b_boson_0_truth_pdgId);
   fChain->SetBranchAddress("boson_0_truth_q", &boson_0_truth_q, &b_boson_0_truth_q);
   fChain->SetBranchAddress("boson_0_truth_status", &boson_0_truth_status, &b_boson_0_truth_status);
   fChain->SetBranchAddress("channel_index", &channel_index, &b_channel_index);
   fChain->SetBranchAddress("dijet_p4", &dijet_p4, &b_dijet_p4);
   fChain->SetBranchAddress("ditau_coll_approx", &ditau_coll_approx, &b_ditau_coll_approx);
   fChain->SetBranchAddress("ditau_coll_approx_m", &ditau_coll_approx_m, &b_ditau_coll_approx_m);
   fChain->SetBranchAddress("ditau_coll_approx_x0", &ditau_coll_approx_x0, &b_ditau_coll_approx_x0);
   fChain->SetBranchAddress("ditau_coll_approx_x1", &ditau_coll_approx_x1, &b_ditau_coll_approx_x1);
   fChain->SetBranchAddress("ditau_cosalpha", &ditau_cosalpha, &b_ditau_cosalpha);
   fChain->SetBranchAddress("ditau_deta", &ditau_deta, &b_ditau_deta);
   fChain->SetBranchAddress("ditau_dphi", &ditau_dphi, &b_ditau_dphi);
   fChain->SetBranchAddress("ditau_dr", &ditau_dr, &b_ditau_dr);
   fChain->SetBranchAddress("ditau_higgspt", &ditau_higgspt, &b_ditau_higgspt);
   fChain->SetBranchAddress("ditau_matched", &ditau_matched, &b_ditau_matched);
   fChain->SetBranchAddress("ditau_matched_cosalpha", &ditau_matched_cosalpha, &b_ditau_matched_cosalpha);
   fChain->SetBranchAddress("ditau_matched_deta", &ditau_matched_deta, &b_ditau_matched_deta);
   fChain->SetBranchAddress("ditau_matched_dphi", &ditau_matched_dphi, &b_ditau_matched_dphi);
   fChain->SetBranchAddress("ditau_matched_dr", &ditau_matched_dr, &b_ditau_matched_dr);
   fChain->SetBranchAddress("ditau_matched_p4", &ditau_matched_p4, &b_ditau_matched_p4);
   fChain->SetBranchAddress("ditau_matched_qxq", &ditau_matched_qxq, &b_ditau_matched_qxq);
   fChain->SetBranchAddress("ditau_matched_scal_sum_pt", &ditau_matched_scal_sum_pt, &b_ditau_matched_scal_sum_pt);
   fChain->SetBranchAddress("ditau_matched_vis_cosalpha", &ditau_matched_vis_cosalpha, &b_ditau_matched_vis_cosalpha);
   fChain->SetBranchAddress("ditau_matched_vis_deta", &ditau_matched_vis_deta, &b_ditau_matched_vis_deta);
   fChain->SetBranchAddress("ditau_matched_vis_dphi", &ditau_matched_vis_dphi, &b_ditau_matched_vis_dphi);
   fChain->SetBranchAddress("ditau_matched_vis_dr", &ditau_matched_vis_dr, &b_ditau_matched_vis_dr);
   fChain->SetBranchAddress("ditau_matched_vis_mass", &ditau_matched_vis_mass, &b_ditau_matched_vis_mass);
   fChain->SetBranchAddress("ditau_matched_vis_scal_sum_pt", &ditau_matched_vis_scal_sum_pt, &b_ditau_matched_vis_scal_sum_pt);
   fChain->SetBranchAddress("ditau_matched_vis_vect_sum_pt", &ditau_matched_vis_vect_sum_pt, &b_ditau_matched_vis_vect_sum_pt);
   fChain->SetBranchAddress("ditau_met_centrality", &ditau_met_centrality, &b_ditau_met_centrality);
   fChain->SetBranchAddress("ditau_met_lep0_cos_dphi", &ditau_met_lep0_cos_dphi, &b_ditau_met_lep0_cos_dphi);
   fChain->SetBranchAddress("ditau_met_lep1_cos_dphi", &ditau_met_lep1_cos_dphi, &b_ditau_met_lep1_cos_dphi);
   fChain->SetBranchAddress("ditau_met_min_dphi", &ditau_met_min_dphi, &b_ditau_met_min_dphi);
   fChain->SetBranchAddress("ditau_met_sum_cos_dphi", &ditau_met_sum_cos_dphi, &b_ditau_met_sum_cos_dphi);
   fChain->SetBranchAddress("ditau_mmc_mlm_fit_status", &ditau_mmc_mlm_fit_status, &b_ditau_mmc_mlm_fit_status);
   fChain->SetBranchAddress("ditau_mmc_mlm_m", &ditau_mmc_mlm_m, &b_ditau_mmc_mlm_m);
   fChain->SetBranchAddress("ditau_mt_lep0_met", &ditau_mt_lep0_met, &b_ditau_mt_lep0_met);
   fChain->SetBranchAddress("ditau_mt_lep1_met", &ditau_mt_lep1_met, &b_ditau_mt_lep1_met);
   fChain->SetBranchAddress("ditau_p4", &ditau_p4, &b_ditau_p4);
   fChain->SetBranchAddress("ditau_qxq", &ditau_qxq, &b_ditau_qxq);
   fChain->SetBranchAddress("ditau_scal_sum_pt", &ditau_scal_sum_pt, &b_ditau_scal_sum_pt);
   fChain->SetBranchAddress("event_clean_EC_TightBad", &event_clean_EC_TightBad, &b_event_clean_EC_TightBad);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("is_dijet_centrality", &is_dijet_centrality, &b_is_dijet_centrality);
   fChain->SetBranchAddress("jet_0_b_tagged", &jet_0_b_tagged, &b_jet_0_b_tagged);
   fChain->SetBranchAddress("jet_0_fjvt", &jet_0_fjvt, &b_jet_0_fjvt);
   fChain->SetBranchAddress("jet_0_flavorlabel_part", &jet_0_flavorlabel_part, &b_jet_0_flavorlabel_part);
   fChain->SetBranchAddress("jet_0_is_Jvt_HS", &jet_0_is_Jvt_HS, &b_jet_0_is_Jvt_HS);
   fChain->SetBranchAddress("jet_0_jvt", &jet_0_jvt, &b_jet_0_jvt);
   fChain->SetBranchAddress("jet_0_p4", &jet_0_p4, &b_jet_0_p4);
   fChain->SetBranchAddress("jet_0_width", &jet_0_width, &b_jet_0_width);
   fChain->SetBranchAddress("jet_0_wztruth_p4", &jet_0_wztruth_p4, &b_jet_0_wztruth_p4);
   fChain->SetBranchAddress("jet_0_wztruth_pdgid", &jet_0_wztruth_pdgid, &b_jet_0_wztruth_pdgid);
   fChain->SetBranchAddress("jet_1_b_tagged", &jet_1_b_tagged, &b_jet_1_b_tagged);
   fChain->SetBranchAddress("jet_1_fjvt", &jet_1_fjvt, &b_jet_1_fjvt);
   fChain->SetBranchAddress("jet_1_flavorlabel_part", &jet_1_flavorlabel_part, &b_jet_1_flavorlabel_part);
   fChain->SetBranchAddress("jet_1_is_Jvt_HS", &jet_1_is_Jvt_HS, &b_jet_1_is_Jvt_HS);
   fChain->SetBranchAddress("jet_1_jvt", &jet_1_jvt, &b_jet_1_jvt);
   fChain->SetBranchAddress("jet_1_p4", &jet_1_p4, &b_jet_1_p4);
   fChain->SetBranchAddress("jet_1_width", &jet_1_width, &b_jet_1_width);
   fChain->SetBranchAddress("jet_1_wztruth_p4", &jet_1_wztruth_p4, &b_jet_1_wztruth_p4);
   fChain->SetBranchAddress("jet_1_wztruth_pdgid", &jet_1_wztruth_pdgid, &b_jet_1_wztruth_pdgid);
   fChain->SetBranchAddress("jet_2_b_tagged", &jet_2_b_tagged, &b_jet_2_b_tagged);
   fChain->SetBranchAddress("jet_2_fjvt", &jet_2_fjvt, &b_jet_2_fjvt);
   fChain->SetBranchAddress("jet_2_flavorlabel_part", &jet_2_flavorlabel_part, &b_jet_2_flavorlabel_part);
   fChain->SetBranchAddress("jet_2_is_Jvt_HS", &jet_2_is_Jvt_HS, &b_jet_2_is_Jvt_HS);
   fChain->SetBranchAddress("jet_2_jvt", &jet_2_jvt, &b_jet_2_jvt);
   fChain->SetBranchAddress("jet_2_p4", &jet_2_p4, &b_jet_2_p4);
   fChain->SetBranchAddress("jet_2_width", &jet_2_width, &b_jet_2_width);
   fChain->SetBranchAddress("jet_2_wztruth_p4", &jet_2_wztruth_p4, &b_jet_2_wztruth_p4);
   fChain->SetBranchAddress("jet_2_wztruth_pdgid", &jet_2_wztruth_pdgid, &b_jet_2_wztruth_pdgid);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_0_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_0_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_0_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_0_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_0_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_0_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_0_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_0_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_0_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_0_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_0_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_0_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_1_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_1_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_1_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_1_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_1_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_1_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_1_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_1_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_1_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_1_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_1_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_1_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_2_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_2_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_2_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_2_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_2_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_2_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_2_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_B_2_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_B_2_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_B_2_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_B_2_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_B_2_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_0_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_0_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_0_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_0_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_0_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_0_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_0_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_0_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_0_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_0_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_0_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_0_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_1_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_1_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_1_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_1_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_1_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_1_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_1_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_1_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_1_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_1_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_1_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_1_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_2_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_2_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_2_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_2_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_2_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_2_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_2_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_2_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_2_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_2_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_2_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_2_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_3_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_3_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_3_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_3_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_3_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_3_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_3_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_C_3_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_C_3_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_C_3_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_C_3_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_C_3_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_0_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_0_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_0_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_0_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_0_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_0_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_0_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_0_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_0_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_0_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_0_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_0_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_1_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_1_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_1_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_1_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_1_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_1_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_1_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_1_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_1_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_1_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_1_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_1_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_2_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_2_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_2_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_2_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_2_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_2_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_2_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_2_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_2_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_2_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_2_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_2_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_3_1down_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_3_1down_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_3_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_3_1down_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_3_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_3_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_3_1up_global_effSF_MV2c10", &jet_FT_EFF_Eigen_Light_3_1up_global_effSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_3_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_Eigen_Light_3_1up_global_ineffSF_MV2c10", &jet_FT_EFF_Eigen_Light_3_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_Eigen_Light_3_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_1down_global_effSF_MV2c10", &jet_FT_EFF_extrapolation_1down_global_effSF_MV2c10, &b_jet_FT_EFF_extrapolation_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_1down_global_ineffSF_MV2c10", &jet_FT_EFF_extrapolation_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_extrapolation_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_1up_global_effSF_MV2c10", &jet_FT_EFF_extrapolation_1up_global_effSF_MV2c10, &b_jet_FT_EFF_extrapolation_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_1up_global_ineffSF_MV2c10", &jet_FT_EFF_extrapolation_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_extrapolation_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_from_charm_1down_global_effSF_MV2c10", &jet_FT_EFF_extrapolation_from_charm_1down_global_effSF_MV2c10, &b_jet_FT_EFF_extrapolation_from_charm_1down_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_from_charm_1down_global_ineffSF_MV2c10", &jet_FT_EFF_extrapolation_from_charm_1down_global_ineffSF_MV2c10, &b_jet_FT_EFF_extrapolation_from_charm_1down_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_from_charm_1up_global_effSF_MV2c10", &jet_FT_EFF_extrapolation_from_charm_1up_global_effSF_MV2c10, &b_jet_FT_EFF_extrapolation_from_charm_1up_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_FT_EFF_extrapolation_from_charm_1up_global_ineffSF_MV2c10", &jet_FT_EFF_extrapolation_from_charm_1up_global_ineffSF_MV2c10, &b_jet_FT_EFF_extrapolation_from_charm_1up_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("jet_JET_JvtEfficiency_1down_central_jets_global_effSF_JVT", &jet_JET_JvtEfficiency_1down_central_jets_global_effSF_JVT, &b_jet_JET_JvtEfficiency_1down_central_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_JET_JvtEfficiency_1down_central_jets_global_ineffSF_JVT", &jet_JET_JvtEfficiency_1down_central_jets_global_ineffSF_JVT, &b_jet_JET_JvtEfficiency_1down_central_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_JET_JvtEfficiency_1up_central_jets_global_effSF_JVT", &jet_JET_JvtEfficiency_1up_central_jets_global_effSF_JVT, &b_jet_JET_JvtEfficiency_1up_central_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_JET_JvtEfficiency_1up_central_jets_global_ineffSF_JVT", &jet_JET_JvtEfficiency_1up_central_jets_global_ineffSF_JVT, &b_jet_JET_JvtEfficiency_1up_central_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_JET_fJvtEfficiency_1down_forward_jets_global_effSF_JVT", &jet_JET_fJvtEfficiency_1down_forward_jets_global_effSF_JVT, &b_jet_JET_fJvtEfficiency_1down_forward_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_JET_fJvtEfficiency_1down_forward_jets_global_ineffSF_JVT", &jet_JET_fJvtEfficiency_1down_forward_jets_global_ineffSF_JVT, &b_jet_JET_fJvtEfficiency_1down_forward_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_JET_fJvtEfficiency_1up_forward_jets_global_effSF_JVT", &jet_JET_fJvtEfficiency_1up_forward_jets_global_effSF_JVT, &b_jet_JET_fJvtEfficiency_1up_forward_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_JET_fJvtEfficiency_1up_forward_jets_global_ineffSF_JVT", &jet_JET_fJvtEfficiency_1up_forward_jets_global_ineffSF_JVT, &b_jet_JET_fJvtEfficiency_1up_forward_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_NOMINAL_central_jets_global_effSF_JVT", &jet_NOMINAL_central_jets_global_effSF_JVT, &b_jet_NOMINAL_central_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_NOMINAL_central_jets_global_ineffSF_JVT", &jet_NOMINAL_central_jets_global_ineffSF_JVT, &b_jet_NOMINAL_central_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_NOMINAL_forward_jets_global_effSF_JVT", &jet_NOMINAL_forward_jets_global_effSF_JVT, &b_jet_NOMINAL_forward_jets_global_effSF_JVT);
   fChain->SetBranchAddress("jet_NOMINAL_forward_jets_global_ineffSF_JVT", &jet_NOMINAL_forward_jets_global_ineffSF_JVT, &b_jet_NOMINAL_forward_jets_global_ineffSF_JVT);
   fChain->SetBranchAddress("jet_NOMINAL_global_effSF_MV2c10", &jet_NOMINAL_global_effSF_MV2c10, &b_jet_NOMINAL_global_effSF_MV2c10);
   fChain->SetBranchAddress("jet_NOMINAL_global_ineffSF_MV2c10", &jet_NOMINAL_global_ineffSF_MV2c10, &b_jet_NOMINAL_global_ineffSF_MV2c10);
   fChain->SetBranchAddress("leplep_fake_tf_bin", &leplep_fake_tf_bin, &b_leplep_fake_tf_bin);
   fChain->SetBranchAddress("lepton_eta_centrality", &lepton_eta_centrality, &b_lepton_eta_centrality);
   fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
   fChain->SetBranchAddress("met_hpto_p4", &met_hpto_p4, &b_met_hpto_p4);
   fChain->SetBranchAddress("met_more_met_et_ele", &met_more_met_et_ele, &b_met_more_met_et_ele);
   fChain->SetBranchAddress("met_more_met_et_jet", &met_more_met_et_jet, &b_met_more_met_et_jet);
   fChain->SetBranchAddress("met_more_met_et_muon", &met_more_met_et_muon, &b_met_more_met_et_muon);
   fChain->SetBranchAddress("met_more_met_et_pho", &met_more_met_et_pho, &b_met_more_met_et_pho);
   fChain->SetBranchAddress("met_more_met_et_soft", &met_more_met_et_soft, &b_met_more_met_et_soft);
   fChain->SetBranchAddress("met_more_met_et_tau", &met_more_met_et_tau, &b_met_more_met_et_tau);
   fChain->SetBranchAddress("met_more_met_phi_ele", &met_more_met_phi_ele, &b_met_more_met_phi_ele);
   fChain->SetBranchAddress("met_more_met_phi_jet", &met_more_met_phi_jet, &b_met_more_met_phi_jet);
   fChain->SetBranchAddress("met_more_met_phi_muon", &met_more_met_phi_muon, &b_met_more_met_phi_muon);
   fChain->SetBranchAddress("met_more_met_phi_pho", &met_more_met_phi_pho, &b_met_more_met_phi_pho);
   fChain->SetBranchAddress("met_more_met_phi_soft", &met_more_met_phi_soft, &b_met_more_met_phi_soft);
   fChain->SetBranchAddress("met_more_met_phi_tau", &met_more_met_phi_tau, &b_met_more_met_phi_tau);
   fChain->SetBranchAddress("met_more_met_sumet_ele", &met_more_met_sumet_ele, &b_met_more_met_sumet_ele);
   fChain->SetBranchAddress("met_more_met_sumet_jet", &met_more_met_sumet_jet, &b_met_more_met_sumet_jet);
   fChain->SetBranchAddress("met_more_met_sumet_muon", &met_more_met_sumet_muon, &b_met_more_met_sumet_muon);
   fChain->SetBranchAddress("met_more_met_sumet_pho", &met_more_met_sumet_pho, &b_met_more_met_sumet_pho);
   fChain->SetBranchAddress("met_more_met_sumet_soft", &met_more_met_sumet_soft, &b_met_more_met_sumet_soft);
   fChain->SetBranchAddress("met_more_met_sumet_tau", &met_more_met_sumet_tau, &b_met_more_met_sumet_tau);
   fChain->SetBranchAddress("met_p4", &met_p4, &b_met_p4);
   fChain->SetBranchAddress("met_sign_met_over_sqrt_ht", &met_sign_met_over_sqrt_ht, &b_met_sign_met_over_sqrt_ht);
   fChain->SetBranchAddress("met_sign_met_over_sqrt_sumet", &met_sign_met_over_sqrt_sumet, &b_met_sign_met_over_sqrt_sumet);
   fChain->SetBranchAddress("met_sign_met_rho", &met_sign_met_rho, &b_met_sign_met_rho);
   fChain->SetBranchAddress("met_sign_met_rho_ttdir", &met_sign_met_rho_ttdir, &b_met_sign_met_rho_ttdir);
   fChain->SetBranchAddress("met_sign_met_sig_directional", &met_sign_met_sig_directional, &b_met_sign_met_sig_directional);
   fChain->SetBranchAddress("met_sign_met_sig_directional_ttdir", &met_sign_met_sig_directional_ttdir, &b_met_sign_met_sig_directional_ttdir);
   fChain->SetBranchAddress("met_sign_met_significance", &met_sign_met_significance, &b_met_sign_met_significance);
   fChain->SetBranchAddress("met_sign_met_significance_ttdir", &met_sign_met_significance_ttdir, &b_met_sign_met_significance_ttdir);
   fChain->SetBranchAddress("met_sign_met_valL", &met_sign_met_valL, &b_met_sign_met_valL);
   fChain->SetBranchAddress("met_sign_met_valL_ttdir", &met_sign_met_valL_ttdir, &b_met_sign_met_valL_ttdir);
   fChain->SetBranchAddress("met_sign_met_varT", &met_sign_met_varT, &b_met_sign_met_varT);
   fChain->SetBranchAddress("met_sign_met_varT_ttdir", &met_sign_met_varT_ttdir, &b_met_sign_met_varT_ttdir);
   fChain->SetBranchAddress("met_sumet", &met_sumet, &b_met_sumet);
   fChain->SetBranchAddress("met_truth_p4", &met_truth_p4, &b_met_truth_p4);
   fChain->SetBranchAddress("met_truth_sumet", &met_truth_sumet, &b_met_truth_sumet);
   fChain->SetBranchAddress("mva_random_number", &mva_random_number, &b_mva_random_number);
   fChain->SetBranchAddress("n_actual_int", &n_actual_int, &b_n_actual_int);
   fChain->SetBranchAddress("n_actual_int_cor", &n_actual_int_cor, &b_n_actual_int_cor);
   fChain->SetBranchAddress("n_avg_int", &n_avg_int, &b_n_avg_int);
   fChain->SetBranchAddress("n_avg_int_cor", &n_avg_int_cor, &b_n_avg_int_cor);
   fChain->SetBranchAddress("n_bjets", &n_bjets, &b_n_bjets);
   fChain->SetBranchAddress("n_electrons", &n_electrons, &b_n_electrons);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("n_jets_30", &n_jets_30, &b_n_jets_30);
   fChain->SetBranchAddress("n_jets_40", &n_jets_40, &b_n_jets_40);
   fChain->SetBranchAddress("n_jets_central", &n_jets_central, &b_n_jets_central);
   fChain->SetBranchAddress("n_jets_central_30", &n_jets_central_30, &b_n_jets_central_30);
   fChain->SetBranchAddress("n_jets_central_40", &n_jets_central_40, &b_n_jets_central_40);
   fChain->SetBranchAddress("n_jets_forward", &n_jets_forward, &b_n_jets_forward);
   fChain->SetBranchAddress("n_jets_forward_30", &n_jets_forward_30, &b_n_jets_forward_30);
   fChain->SetBranchAddress("n_jets_forward_40", &n_jets_forward_40, &b_n_jets_forward_40);
   fChain->SetBranchAddress("n_jets_l1j25", &n_jets_l1j25, &b_n_jets_l1j25);
   fChain->SetBranchAddress("n_jets_mc_hs", &n_jets_mc_hs, &b_n_jets_mc_hs);
   fChain->SetBranchAddress("n_muons", &n_muons, &b_n_muons);
   fChain->SetBranchAddress("n_photons", &n_photons, &b_n_photons);
   fChain->SetBranchAddress("n_pvx", &n_pvx, &b_n_pvx);
   fChain->SetBranchAddress("n_taus", &n_taus, &b_n_taus);
   fChain->SetBranchAddress("n_taus_loose", &n_taus_loose, &b_n_taus_loose);
   fChain->SetBranchAddress("n_taus_medium", &n_taus_medium, &b_n_taus_medium);
   fChain->SetBranchAddress("n_taus_tight", &n_taus_tight, &b_n_taus_tight);
   fChain->SetBranchAddress("n_taus_veryloose", &n_taus_veryloose, &b_n_taus_veryloose);
   fChain->SetBranchAddress("n_truth_gluon_jets", &n_truth_gluon_jets, &b_n_truth_gluon_jets);
   fChain->SetBranchAddress("n_truth_jets", &n_truth_jets, &b_n_truth_jets);
   fChain->SetBranchAddress("n_truth_jets_pt20_eta45", &n_truth_jets_pt20_eta45, &b_n_truth_jets_pt20_eta45);
   fChain->SetBranchAddress("n_truth_quark_jets", &n_truth_quark_jets, &b_n_truth_quark_jets);
   fChain->SetBranchAddress("n_vx", &n_vx, &b_n_vx);
   fChain->SetBranchAddress("primary_vertex", &primary_vertex, &b_primary_vertex);
   fChain->SetBranchAddress("primary_vertex_v", &primary_vertex_v, &b_primary_vertex_v);
   fChain->SetBranchAddress("pt_total", &pt_total, &b_pt_total);
   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("scalar_sum_pt", &scalar_sum_pt, &b_scalar_sum_pt);
   fChain->SetBranchAddress("tau_0", &tau_0, &b_tau_0);
   fChain->SetBranchAddress("tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_0_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_0_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk", &tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk, &b_tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk", &tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk, &b_tau_0_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose", &tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose, &b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly", &tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly, &b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_0_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose", &tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose, &b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly", &tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly, &b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_0_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose", &tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose, &b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly", &tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly, &b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_0_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose", &tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose, &b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly", &tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly, &b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_0_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium", &tau_0_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium, &b_tau_0_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_0_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_0_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_0_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_0_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_0_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_0_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_0_NOMINAL_EleEffSF_offline_RecoTrk", &tau_0_NOMINAL_EleEffSF_offline_RecoTrk, &b_tau_0_NOMINAL_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_0_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_0_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_IsoFCLoose", &tau_0_NOMINAL_MuEffSF_IsoFCLoose, &b_tau_0_NOMINAL_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly", &tau_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly, &b_tau_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_NOMINAL_MuEffSF_Reco_QualMedium", &tau_0_NOMINAL_MuEffSF_Reco_QualMedium, &b_tau_0_NOMINAL_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_0_cluster_eta", &tau_0_cluster_eta, &b_tau_0_cluster_eta);
   fChain->SetBranchAddress("tau_0_cluster_eta_be2", &tau_0_cluster_eta_be2, &b_tau_0_cluster_eta_be2);
   fChain->SetBranchAddress("tau_0_cluster_pt", &tau_0_cluster_pt, &b_tau_0_cluster_pt);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_2e12_lhloose_L12EM10VH", &tau_0_electron_trig_HLT_2e12_lhloose_L12EM10VH, &b_tau_0_electron_trig_HLT_2e12_lhloose_L12EM10VH);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_2e17_lhvloose_nod0", &tau_0_electron_trig_HLT_2e17_lhvloose_nod0, &b_tau_0_electron_trig_HLT_2e17_lhvloose_nod0);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI", &tau_0_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI, &b_tau_0_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e120_lhloose", &tau_0_electron_trig_HLT_e120_lhloose, &b_tau_0_electron_trig_HLT_e120_lhloose);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e140_lhloose_nod0", &tau_0_electron_trig_HLT_e140_lhloose_nod0, &b_tau_0_electron_trig_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e24_lhmedium_L1EM20VH", &tau_0_electron_trig_HLT_e24_lhmedium_L1EM20VH, &b_tau_0_electron_trig_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e26_lhtight_nod0_ivarloose", &tau_0_electron_trig_HLT_e26_lhtight_nod0_ivarloose, &b_tau_0_electron_trig_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e60_lhmedium", &tau_0_electron_trig_HLT_e60_lhmedium, &b_tau_0_electron_trig_HLT_e60_lhmedium);
   fChain->SetBranchAddress("tau_0_electron_trig_HLT_e60_lhmedium_nod0", &tau_0_electron_trig_HLT_e60_lhmedium_nod0, &b_tau_0_electron_trig_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("tau_0_electron_trig_trigger_matched", &tau_0_electron_trig_trigger_matched, &b_tau_0_electron_trig_trigger_matched);
   fChain->SetBranchAddress("tau_0_emu_trig_HLT_e17_lhloose_mu14", &tau_0_emu_trig_HLT_e17_lhloose_mu14, &b_tau_0_emu_trig_HLT_e17_lhloose_mu14);
   fChain->SetBranchAddress("tau_0_emu_trig_HLT_e17_lhloose_nod0_mu14", &tau_0_emu_trig_HLT_e17_lhloose_nod0_mu14, &b_tau_0_emu_trig_HLT_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("tau_0_emu_trig_trigger_matched", &tau_0_emu_trig_trigger_matched, &b_tau_0_emu_trig_trigger_matched);
   fChain->SetBranchAddress("tau_0_id_bad", &tau_0_id_bad, &b_tau_0_id_bad);
   fChain->SetBranchAddress("tau_0_id_charge", &tau_0_id_charge, &b_tau_0_id_charge);
   fChain->SetBranchAddress("tau_0_id_loose", &tau_0_id_loose, &b_tau_0_id_loose);
   fChain->SetBranchAddress("tau_0_id_medium", &tau_0_id_medium, &b_tau_0_id_medium);
   fChain->SetBranchAddress("tau_0_id_tight", &tau_0_id_tight, &b_tau_0_id_tight);
   fChain->SetBranchAddress("tau_0_id_veryloose", &tau_0_id_veryloose, &b_tau_0_id_veryloose);
   fChain->SetBranchAddress("tau_0_iso_FCHighPtCaloOnly", &tau_0_iso_FCHighPtCaloOnly, &b_tau_0_iso_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_0_iso_FCLoose", &tau_0_iso_FCLoose, &b_tau_0_iso_FCLoose);
   fChain->SetBranchAddress("tau_0_iso_FCTight", &tau_0_iso_FCTight, &b_tau_0_iso_FCTight);
   fChain->SetBranchAddress("tau_0_iso_FCTightTrackOnly", &tau_0_iso_FCTightTrackOnly, &b_tau_0_iso_FCTightTrackOnly);
   fChain->SetBranchAddress("tau_0_iso_FixedCutHighPtTrackOnly", &tau_0_iso_FixedCutHighPtTrackOnly, &b_tau_0_iso_FixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_0_iso_Gradient", &tau_0_iso_Gradient, &b_tau_0_iso_Gradient);
   fChain->SetBranchAddress("tau_0_iso_ptvarcone20", &tau_0_iso_ptvarcone20, &b_tau_0_iso_ptvarcone20);
   fChain->SetBranchAddress("tau_0_iso_ptvarcone30", &tau_0_iso_ptvarcone30, &b_tau_0_iso_ptvarcone30);
   fChain->SetBranchAddress("tau_0_matched", &tau_0_matched, &b_tau_0_matched);
   fChain->SetBranchAddress("tau_0_matched_classifierParticleOrigin", &tau_0_matched_classifierParticleOrigin, &b_tau_0_matched_classifierParticleOrigin);
   fChain->SetBranchAddress("tau_0_matched_classifierParticleType", &tau_0_matched_classifierParticleType, &b_tau_0_matched_classifierParticleType);
   fChain->SetBranchAddress("tau_0_matched_isHad", &tau_0_matched_isHad, &b_tau_0_matched_isHad);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau", &tau_0_matched_leptonic_tau, &b_tau_0_matched_leptonic_tau);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_classifierParticleOrigin", &tau_0_matched_leptonic_tau_classifierParticleOrigin, &b_tau_0_matched_leptonic_tau_classifierParticleOrigin);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_classifierParticleType", &tau_0_matched_leptonic_tau_classifierParticleType, &b_tau_0_matched_leptonic_tau_classifierParticleType);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_invis_p4", &tau_0_matched_leptonic_tau_invis_p4, &b_tau_0_matched_leptonic_tau_invis_p4);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_mother_pdgId", &tau_0_matched_leptonic_tau_mother_pdgId, &b_tau_0_matched_leptonic_tau_mother_pdgId);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_mother_status", &tau_0_matched_leptonic_tau_mother_status, &b_tau_0_matched_leptonic_tau_mother_status);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_origin", &tau_0_matched_leptonic_tau_origin, &b_tau_0_matched_leptonic_tau_origin);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_p4", &tau_0_matched_leptonic_tau_p4, &b_tau_0_matched_leptonic_tau_p4);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_pdgId", &tau_0_matched_leptonic_tau_pdgId, &b_tau_0_matched_leptonic_tau_pdgId);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_pz", &tau_0_matched_leptonic_tau_pz, &b_tau_0_matched_leptonic_tau_pz);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_q", &tau_0_matched_leptonic_tau_q, &b_tau_0_matched_leptonic_tau_q);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_status", &tau_0_matched_leptonic_tau_status, &b_tau_0_matched_leptonic_tau_status);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_type", &tau_0_matched_leptonic_tau_type, &b_tau_0_matched_leptonic_tau_type);
   fChain->SetBranchAddress("tau_0_matched_leptonic_tau_vis_p4", &tau_0_matched_leptonic_tau_vis_p4, &b_tau_0_matched_leptonic_tau_vis_p4);
   fChain->SetBranchAddress("tau_0_matched_mother_pdgId", &tau_0_matched_mother_pdgId, &b_tau_0_matched_mother_pdgId);
   fChain->SetBranchAddress("tau_0_matched_mother_status", &tau_0_matched_mother_status, &b_tau_0_matched_mother_status);
   fChain->SetBranchAddress("tau_0_matched_p4", &tau_0_matched_p4, &b_tau_0_matched_p4);
   fChain->SetBranchAddress("tau_0_matched_pdgId", &tau_0_matched_pdgId, &b_tau_0_matched_pdgId);
   fChain->SetBranchAddress("tau_0_matched_q", &tau_0_matched_q, &b_tau_0_matched_q);
   fChain->SetBranchAddress("tau_0_matched_status", &tau_0_matched_status, &b_tau_0_matched_status);
   fChain->SetBranchAddress("tau_0_matched_type", &tau_0_matched_type, &b_tau_0_matched_type);
   fChain->SetBranchAddress("tau_0_muonAuthor", &tau_0_muonAuthor, &b_tau_0_muonAuthor);
   fChain->SetBranchAddress("tau_0_muonType", &tau_0_muonType, &b_tau_0_muonType);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu18_mu8noL1", &tau_0_muon_trig_HLT_mu18_mu8noL1, &b_tau_0_muon_trig_HLT_mu18_mu8noL1);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu20_iloose_L1MU15", &tau_0_muon_trig_HLT_mu20_iloose_L1MU15, &b_tau_0_muon_trig_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu20_mu8noL1", &tau_0_muon_trig_HLT_mu20_mu8noL1, &b_tau_0_muon_trig_HLT_mu20_mu8noL1);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu22_mu8noL1", &tau_0_muon_trig_HLT_mu22_mu8noL1, &b_tau_0_muon_trig_HLT_mu22_mu8noL1);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu26_ivarmedium", &tau_0_muon_trig_HLT_mu26_ivarmedium, &b_tau_0_muon_trig_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("tau_0_muon_trig_HLT_mu50", &tau_0_muon_trig_HLT_mu50, &b_tau_0_muon_trig_HLT_mu50);
   fChain->SetBranchAddress("tau_0_muon_trig_trigger_matched", &tau_0_muon_trig_trigger_matched, &b_tau_0_muon_trig_trigger_matched);
   fChain->SetBranchAddress("tau_0_origin", &tau_0_origin, &b_tau_0_origin);
   fChain->SetBranchAddress("tau_0_p4", &tau_0_p4, &b_tau_0_p4);
   fChain->SetBranchAddress("tau_0_q", &tau_0_q, &b_tau_0_q);
   fChain->SetBranchAddress("tau_0_trk_d0", &tau_0_trk_d0, &b_tau_0_trk_d0);
   fChain->SetBranchAddress("tau_0_trk_d0_sig", &tau_0_trk_d0_sig, &b_tau_0_trk_d0_sig);
   fChain->SetBranchAddress("tau_0_trk_pt", &tau_0_trk_pt, &b_tau_0_trk_pt);
   fChain->SetBranchAddress("tau_0_trk_pt_error", &tau_0_trk_pt_error, &b_tau_0_trk_pt_error);
   fChain->SetBranchAddress("tau_0_trk_pvx_z0", &tau_0_trk_pvx_z0, &b_tau_0_trk_pvx_z0);
   fChain->SetBranchAddress("tau_0_trk_pvx_z0_sig", &tau_0_trk_pvx_z0_sig, &b_tau_0_trk_pvx_z0_sig);
   fChain->SetBranchAddress("tau_0_trk_pvx_z0_sintheta", &tau_0_trk_pvx_z0_sintheta, &b_tau_0_trk_pvx_z0_sintheta);
   fChain->SetBranchAddress("tau_0_trk_vx", &tau_0_trk_vx, &b_tau_0_trk_vx);
   fChain->SetBranchAddress("tau_0_trk_vx_v", &tau_0_trk_vx_v, &b_tau_0_trk_vx_v);
   fChain->SetBranchAddress("tau_0_trk_z0", &tau_0_trk_z0, &b_tau_0_trk_z0);
   fChain->SetBranchAddress("tau_0_trk_z0_sig", &tau_0_trk_z0_sig, &b_tau_0_trk_z0_sig);
   fChain->SetBranchAddress("tau_0_trk_z0_sintheta", &tau_0_trk_z0_sintheta, &b_tau_0_trk_z0_sintheta);
   fChain->SetBranchAddress("tau_0_type", &tau_0_type, &b_tau_0_type);
   fChain->SetBranchAddress("tau_1", &tau_1, &b_tau_1);
   fChain->SetBranchAddress("tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_1_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_1_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk", &tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk, &b_tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk", &tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk, &b_tau_1_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1down_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR_1up_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose", &tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose, &b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly", &tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly, &b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_1_MUON_EFF_ISO_STAT_1down_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose", &tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose, &b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly", &tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly, &b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_1_MUON_EFF_ISO_STAT_1up_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose", &tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose, &b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly", &tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly, &b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_1_MUON_EFF_ISO_SYS_1down_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose", &tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose, &b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly", &tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly, &b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_1_MUON_EFF_ISO_SYS_1up_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_STAT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_STAT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_STAT_LOWPT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_STAT_LOWPT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_SYS_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_SYS_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_SYS_LOWPT_1down_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium", &tau_1_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium, &b_tau_1_MUON_EFF_RECO_SYS_LOWPT_1up_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigStatUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1down_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_1_MUON_EFF_TrigSystUncertainty_1up_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient", &tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_NOMINAL_EleEffSF_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0_2017_2018_e17_lhvloose_nod0_L1EM15VHI_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly", &tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly, &b_tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient", &tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient, &b_tau_1_NOMINAL_EleEffSF_Isolation_MediumLLH_d0z0_v13_Gradient);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_NOMINAL_EleEffSF_MULTI_L_2015_e17_lhloose_2016_2018_e17_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly", &tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly, &b_tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolFCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient", &tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient, &b_tau_1_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_MediumLLH_d0z0_v13_isolGradient);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13", &tau_1_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13, &b_tau_1_NOMINAL_EleEffSF_offline_LooseAndBLayerLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13", &tau_1_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13, &b_tau_1_NOMINAL_EleEffSF_offline_MediumLLH_d0z0_v13);
   fChain->SetBranchAddress("tau_1_NOMINAL_EleEffSF_offline_RecoTrk", &tau_1_NOMINAL_EleEffSF_offline_RecoTrk, &b_tau_1_NOMINAL_EleEffSF_offline_RecoTrk);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu14_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu18_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu22_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone", &tau_1_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone, &b_tau_1_NOMINAL_MuEffSF_HLT_mu8noL1_QualMedium_IsoNone);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_IsoFCLoose", &tau_1_NOMINAL_MuEffSF_IsoFCLoose, &b_tau_1_NOMINAL_MuEffSF_IsoFCLoose);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_IsoFCTightTrackOnly", &tau_1_NOMINAL_MuEffSF_IsoFCTightTrackOnly, &b_tau_1_NOMINAL_MuEffSF_IsoFCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly", &tau_1_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly, &b_tau_1_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_NOMINAL_MuEffSF_Reco_QualMedium", &tau_1_NOMINAL_MuEffSF_Reco_QualMedium, &b_tau_1_NOMINAL_MuEffSF_Reco_QualMedium);
   fChain->SetBranchAddress("tau_1_cluster_eta", &tau_1_cluster_eta, &b_tau_1_cluster_eta);
   fChain->SetBranchAddress("tau_1_cluster_eta_be2", &tau_1_cluster_eta_be2, &b_tau_1_cluster_eta_be2);
   fChain->SetBranchAddress("tau_1_cluster_pt", &tau_1_cluster_pt, &b_tau_1_cluster_pt);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_2e12_lhloose_L12EM10VH", &tau_1_electron_trig_HLT_2e12_lhloose_L12EM10VH, &b_tau_1_electron_trig_HLT_2e12_lhloose_L12EM10VH);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_2e17_lhvloose_nod0", &tau_1_electron_trig_HLT_2e17_lhvloose_nod0, &b_tau_1_electron_trig_HLT_2e17_lhvloose_nod0);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI", &tau_1_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI, &b_tau_1_electron_trig_HLT_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e120_lhloose", &tau_1_electron_trig_HLT_e120_lhloose, &b_tau_1_electron_trig_HLT_e120_lhloose);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e140_lhloose_nod0", &tau_1_electron_trig_HLT_e140_lhloose_nod0, &b_tau_1_electron_trig_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e24_lhmedium_L1EM20VH", &tau_1_electron_trig_HLT_e24_lhmedium_L1EM20VH, &b_tau_1_electron_trig_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e26_lhtight_nod0_ivarloose", &tau_1_electron_trig_HLT_e26_lhtight_nod0_ivarloose, &b_tau_1_electron_trig_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e60_lhmedium", &tau_1_electron_trig_HLT_e60_lhmedium, &b_tau_1_electron_trig_HLT_e60_lhmedium);
   fChain->SetBranchAddress("tau_1_electron_trig_HLT_e60_lhmedium_nod0", &tau_1_electron_trig_HLT_e60_lhmedium_nod0, &b_tau_1_electron_trig_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("tau_1_electron_trig_trigger_matched", &tau_1_electron_trig_trigger_matched, &b_tau_1_electron_trig_trigger_matched);
   fChain->SetBranchAddress("tau_1_emu_trig_HLT_e17_lhloose_mu14", &tau_1_emu_trig_HLT_e17_lhloose_mu14, &b_tau_1_emu_trig_HLT_e17_lhloose_mu14);
   fChain->SetBranchAddress("tau_1_emu_trig_HLT_e17_lhloose_nod0_mu14", &tau_1_emu_trig_HLT_e17_lhloose_nod0_mu14, &b_tau_1_emu_trig_HLT_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("tau_1_emu_trig_trigger_matched", &tau_1_emu_trig_trigger_matched, &b_tau_1_emu_trig_trigger_matched);
   fChain->SetBranchAddress("tau_1_id_bad", &tau_1_id_bad, &b_tau_1_id_bad);
   fChain->SetBranchAddress("tau_1_id_charge", &tau_1_id_charge, &b_tau_1_id_charge);
   fChain->SetBranchAddress("tau_1_id_loose", &tau_1_id_loose, &b_tau_1_id_loose);
   fChain->SetBranchAddress("tau_1_id_medium", &tau_1_id_medium, &b_tau_1_id_medium);
   fChain->SetBranchAddress("tau_1_id_tight", &tau_1_id_tight, &b_tau_1_id_tight);
   fChain->SetBranchAddress("tau_1_id_veryloose", &tau_1_id_veryloose, &b_tau_1_id_veryloose);
   fChain->SetBranchAddress("tau_1_iso_FCHighPtCaloOnly", &tau_1_iso_FCHighPtCaloOnly, &b_tau_1_iso_FCHighPtCaloOnly);
   fChain->SetBranchAddress("tau_1_iso_FCLoose", &tau_1_iso_FCLoose, &b_tau_1_iso_FCLoose);
   fChain->SetBranchAddress("tau_1_iso_FCTight", &tau_1_iso_FCTight, &b_tau_1_iso_FCTight);
   fChain->SetBranchAddress("tau_1_iso_FCTightTrackOnly", &tau_1_iso_FCTightTrackOnly, &b_tau_1_iso_FCTightTrackOnly);
   fChain->SetBranchAddress("tau_1_iso_FixedCutHighPtTrackOnly", &tau_1_iso_FixedCutHighPtTrackOnly, &b_tau_1_iso_FixedCutHighPtTrackOnly);
   fChain->SetBranchAddress("tau_1_iso_Gradient", &tau_1_iso_Gradient, &b_tau_1_iso_Gradient);
   fChain->SetBranchAddress("tau_1_iso_ptvarcone20", &tau_1_iso_ptvarcone20, &b_tau_1_iso_ptvarcone20);
   fChain->SetBranchAddress("tau_1_iso_ptvarcone30", &tau_1_iso_ptvarcone30, &b_tau_1_iso_ptvarcone30);
   fChain->SetBranchAddress("tau_1_matched", &tau_1_matched, &b_tau_1_matched);
   fChain->SetBranchAddress("tau_1_matched_classifierParticleOrigin", &tau_1_matched_classifierParticleOrigin, &b_tau_1_matched_classifierParticleOrigin);
   fChain->SetBranchAddress("tau_1_matched_classifierParticleType", &tau_1_matched_classifierParticleType, &b_tau_1_matched_classifierParticleType);
   fChain->SetBranchAddress("tau_1_matched_isHad", &tau_1_matched_isHad, &b_tau_1_matched_isHad);
   fChain->SetBranchAddress("tau_1_matched_leptonic_tau_invis_p4", &tau_1_matched_leptonic_tau_invis_p4, &b_tau_1_matched_leptonic_tau_invis_p4);
   fChain->SetBranchAddress("tau_1_matched_leptonic_tau_p4", &tau_1_matched_leptonic_tau_p4, &b_tau_1_matched_leptonic_tau_p4);
   fChain->SetBranchAddress("tau_1_matched_leptonic_tau_vis_p4", &tau_1_matched_leptonic_tau_vis_p4, &b_tau_1_matched_leptonic_tau_vis_p4);
   fChain->SetBranchAddress("tau_1_matched_mother_pdgId", &tau_1_matched_mother_pdgId, &b_tau_1_matched_mother_pdgId);
   fChain->SetBranchAddress("tau_1_matched_mother_status", &tau_1_matched_mother_status, &b_tau_1_matched_mother_status);
   fChain->SetBranchAddress("tau_1_matched_p4", &tau_1_matched_p4, &b_tau_1_matched_p4);
   fChain->SetBranchAddress("tau_1_matched_pdgId", &tau_1_matched_pdgId, &b_tau_1_matched_pdgId);
   fChain->SetBranchAddress("tau_1_matched_q", &tau_1_matched_q, &b_tau_1_matched_q);
   fChain->SetBranchAddress("tau_1_matched_status", &tau_1_matched_status, &b_tau_1_matched_status);
   fChain->SetBranchAddress("tau_1_muonAuthor", &tau_1_muonAuthor, &b_tau_1_muonAuthor);
   fChain->SetBranchAddress("tau_1_muonType", &tau_1_muonType, &b_tau_1_muonType);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu18_mu8noL1", &tau_1_muon_trig_HLT_mu18_mu8noL1, &b_tau_1_muon_trig_HLT_mu18_mu8noL1);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu20_iloose_L1MU15", &tau_1_muon_trig_HLT_mu20_iloose_L1MU15, &b_tau_1_muon_trig_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu20_mu8noL1", &tau_1_muon_trig_HLT_mu20_mu8noL1, &b_tau_1_muon_trig_HLT_mu20_mu8noL1);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu22_mu8noL1", &tau_1_muon_trig_HLT_mu22_mu8noL1, &b_tau_1_muon_trig_HLT_mu22_mu8noL1);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu26_ivarmedium", &tau_1_muon_trig_HLT_mu26_ivarmedium, &b_tau_1_muon_trig_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("tau_1_muon_trig_HLT_mu50", &tau_1_muon_trig_HLT_mu50, &b_tau_1_muon_trig_HLT_mu50);
   fChain->SetBranchAddress("tau_1_muon_trig_trigger_matched", &tau_1_muon_trig_trigger_matched, &b_tau_1_muon_trig_trigger_matched);
   fChain->SetBranchAddress("tau_1_origin", &tau_1_origin, &b_tau_1_origin);
   fChain->SetBranchAddress("tau_1_p4", &tau_1_p4, &b_tau_1_p4);
   fChain->SetBranchAddress("tau_1_q", &tau_1_q, &b_tau_1_q);
   fChain->SetBranchAddress("tau_1_trk_d0", &tau_1_trk_d0, &b_tau_1_trk_d0);
   fChain->SetBranchAddress("tau_1_trk_d0_sig", &tau_1_trk_d0_sig, &b_tau_1_trk_d0_sig);
   fChain->SetBranchAddress("tau_1_trk_pt", &tau_1_trk_pt, &b_tau_1_trk_pt);
   fChain->SetBranchAddress("tau_1_trk_pt_error", &tau_1_trk_pt_error, &b_tau_1_trk_pt_error);
   fChain->SetBranchAddress("tau_1_trk_pvx_z0", &tau_1_trk_pvx_z0, &b_tau_1_trk_pvx_z0);
   fChain->SetBranchAddress("tau_1_trk_pvx_z0_sig", &tau_1_trk_pvx_z0_sig, &b_tau_1_trk_pvx_z0_sig);
   fChain->SetBranchAddress("tau_1_trk_pvx_z0_sintheta", &tau_1_trk_pvx_z0_sintheta, &b_tau_1_trk_pvx_z0_sintheta);
   fChain->SetBranchAddress("tau_1_trk_vx", &tau_1_trk_vx, &b_tau_1_trk_vx);
   fChain->SetBranchAddress("tau_1_trk_vx_v", &tau_1_trk_vx_v, &b_tau_1_trk_vx_v);
   fChain->SetBranchAddress("tau_1_trk_z0", &tau_1_trk_z0, &b_tau_1_trk_z0);
   fChain->SetBranchAddress("tau_1_trk_z0_sig", &tau_1_trk_z0_sig, &b_tau_1_trk_z0_sig);
   fChain->SetBranchAddress("tau_1_trk_z0_sintheta", &tau_1_trk_z0_sintheta, &b_tau_1_trk_z0_sintheta);
   fChain->SetBranchAddress("tau_1_type", &tau_1_type, &b_tau_1_type);
   fChain->SetBranchAddress("tau_eta_centrality", &tau_eta_centrality, &b_tau_eta_centrality);
   fChain->SetBranchAddress("theory_weights_alphaS_down", &theory_weights_alphaS_down, &b_theory_weights_alphaS_down);
   fChain->SetBranchAddress("theory_weights_alphaS_up", &theory_weights_alphaS_up, &b_theory_weights_alphaS_up);
   fChain->SetBranchAddress("theory_weights_nominal", &theory_weights_nominal, &b_theory_weights_nominal);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_0", &theory_weights_pdf_signal_weight_0, &b_theory_weights_pdf_signal_weight_0);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_1", &theory_weights_pdf_signal_weight_1, &b_theory_weights_pdf_signal_weight_1);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_10", &theory_weights_pdf_signal_weight_10, &b_theory_weights_pdf_signal_weight_10);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_11", &theory_weights_pdf_signal_weight_11, &b_theory_weights_pdf_signal_weight_11);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_12", &theory_weights_pdf_signal_weight_12, &b_theory_weights_pdf_signal_weight_12);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_13", &theory_weights_pdf_signal_weight_13, &b_theory_weights_pdf_signal_weight_13);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_14", &theory_weights_pdf_signal_weight_14, &b_theory_weights_pdf_signal_weight_14);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_15", &theory_weights_pdf_signal_weight_15, &b_theory_weights_pdf_signal_weight_15);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_16", &theory_weights_pdf_signal_weight_16, &b_theory_weights_pdf_signal_weight_16);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_17", &theory_weights_pdf_signal_weight_17, &b_theory_weights_pdf_signal_weight_17);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_18", &theory_weights_pdf_signal_weight_18, &b_theory_weights_pdf_signal_weight_18);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_19", &theory_weights_pdf_signal_weight_19, &b_theory_weights_pdf_signal_weight_19);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_2", &theory_weights_pdf_signal_weight_2, &b_theory_weights_pdf_signal_weight_2);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_20", &theory_weights_pdf_signal_weight_20, &b_theory_weights_pdf_signal_weight_20);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_21", &theory_weights_pdf_signal_weight_21, &b_theory_weights_pdf_signal_weight_21);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_22", &theory_weights_pdf_signal_weight_22, &b_theory_weights_pdf_signal_weight_22);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_23", &theory_weights_pdf_signal_weight_23, &b_theory_weights_pdf_signal_weight_23);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_24", &theory_weights_pdf_signal_weight_24, &b_theory_weights_pdf_signal_weight_24);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_25", &theory_weights_pdf_signal_weight_25, &b_theory_weights_pdf_signal_weight_25);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_26", &theory_weights_pdf_signal_weight_26, &b_theory_weights_pdf_signal_weight_26);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_27", &theory_weights_pdf_signal_weight_27, &b_theory_weights_pdf_signal_weight_27);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_28", &theory_weights_pdf_signal_weight_28, &b_theory_weights_pdf_signal_weight_28);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_29", &theory_weights_pdf_signal_weight_29, &b_theory_weights_pdf_signal_weight_29);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_3", &theory_weights_pdf_signal_weight_3, &b_theory_weights_pdf_signal_weight_3);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_4", &theory_weights_pdf_signal_weight_4, &b_theory_weights_pdf_signal_weight_4);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_5", &theory_weights_pdf_signal_weight_5, &b_theory_weights_pdf_signal_weight_5);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_6", &theory_weights_pdf_signal_weight_6, &b_theory_weights_pdf_signal_weight_6);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_7", &theory_weights_pdf_signal_weight_7, &b_theory_weights_pdf_signal_weight_7);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_8", &theory_weights_pdf_signal_weight_8, &b_theory_weights_pdf_signal_weight_8);
   fChain->SetBranchAddress("theory_weights_pdf_signal_weight_9", &theory_weights_pdf_signal_weight_9, &b_theory_weights_pdf_signal_weight_9);
   fChain->SetBranchAddress("theory_weights_qcd_weight_0", &theory_weights_qcd_weight_0, &b_theory_weights_qcd_weight_0);
   fChain->SetBranchAddress("theory_weights_qcd_weight_1", &theory_weights_qcd_weight_1, &b_theory_weights_qcd_weight_1);
   fChain->SetBranchAddress("theory_weights_qcd_weight_2", &theory_weights_qcd_weight_2, &b_theory_weights_qcd_weight_2);
   fChain->SetBranchAddress("theory_weights_qcd_weight_3", &theory_weights_qcd_weight_3, &b_theory_weights_qcd_weight_3);
   fChain->SetBranchAddress("theory_weights_qcd_weight_4", &theory_weights_qcd_weight_4, &b_theory_weights_qcd_weight_4);
   fChain->SetBranchAddress("theory_weights_qcd_weight_5", &theory_weights_qcd_weight_5, &b_theory_weights_qcd_weight_5);
   fChain->SetBranchAddress("theory_weights_qcd_weight_6", &theory_weights_qcd_weight_6, &b_theory_weights_qcd_weight_6);
   fChain->SetBranchAddress("theory_weights_qcd_weight_7", &theory_weights_qcd_weight_7, &b_theory_weights_qcd_weight_7);
   fChain->SetBranchAddress("theory_weights_qcd_weight_8", &theory_weights_qcd_weight_8, &b_theory_weights_qcd_weight_8);
   fChain->SetBranchAddress("truth_passedVBFFilter", &truth_passedVBFFilter, &b_truth_passedVBFFilter);
   fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   Notify();
}

bool tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx
