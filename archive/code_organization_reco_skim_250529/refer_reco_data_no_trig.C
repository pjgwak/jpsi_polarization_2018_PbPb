#include <iostream>
#include <vector>
#include <string>
#include "THn.h"
#include "TH2D.h"
#include <TLorentzVector.h>
#include <TAttMarker.h>
#include "../../headers/commonUtility.h"
#include "../../headers/JpsiUtility.h"
#include "../../headers/HiEvtPlaneList.h"
#include "../../headers/cutsAndBin.h"
#include "../../headers/Style.h"
#include "../../headers/tnp_weight_lowptPbPb_num_den_new.h"


bool IsAcceptable(double pt, double eta);
THnD *create_h4d(const std::string &name, const std::string &title, const std::vector<double> &bins_1, const std::vector<double> &bins_2, const std::vector<double> &bins_3, const std::vector<double> &bins_4);

void reco_data_no_trig(string data_label_ = "all_event", int nevt = -1, int isMC = 0, bool isTnP = true, bool isPtWeight = true)
{
  TStopwatch *t = new TStopwatch;
  t->Start();
  cout << "\n=================================\n";
  cout << "\n Start PbPb Eff Calculation\n";
  cout << "\n=================================\n\n";

  using namespace std;
  using namespace hi; // HiEvtPlaneList.h

  // ===== basic setting ===== //
  // gStyle->SetOptStat(0);


  float ptLow = 0.0, ptHigh = 50.0;
  float yLow = 0.0, yHigh = 2.4;

  int kTrigSel = 12;  // jpsi=12,Upsilon=13, ppJpsi=3
  // Run3 Trig #12: HLT_HIL3DoubleMu2_M2to4p5_Open_v
  float muPtCut = 0;  // 3.5, 1.8
  int kL2filter = 16; // jpsi=16,Upsilon=38, ppJpsi=4
  int kL3filter = 17; // jpsi=17,Upsilon=39, ppJpsi=5

  // jpsi mass
  float massLow = 2.6;
  float massHigh = 3.5;
  int cLow = 0;
  int cHigh = 180;

  // labeling
  TString date_label = data_label_;

  // ===== import oniatree ===== //
  TChain muon_chain("hionia/myTree");
  // TChain ep_chain("hionia/tree");

  if (isMC == 0)
  {
    // PbPb MC prompt
    TString filename = Form("/work/pjgwak/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_no_track_250520/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_no_track_250520/250520_121505/Oniatree_MC0_PbPb_2023.root");
    std::cout << "Adding PbPb(2023) Data sample: " << filename << std::endl;

    muon_chain.Add(filename);
    // ep_chain.Add(filename);
  }
  else if (isMC == 2)
  {
    // not used
    // PbPb MC nonprompt
    TString filename = "";
    std::cout << "Adding PbPb: " << filename << std::endl;
  }

  // ===== build histograms ===== //
  vector<double> cent_bin = {0, 20, 40, 60, 100, 180}; // 0 - 180. divide by 2 to get the cent (%)
  vector<double> fwd_pt_bin = {3, 6.5, 9, 12, 15, 20, 50};
  vector<double> mid_pt_bin = {6.5, 9.0, 12, 15.0, 20, 50};
  vector<double> cos_bin = {-1.0, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0};
  vector<double> phi_bin = {-3.15, -2.835, -2.52, -2.205, -1.89, -1.575, -1.26, -0.945, -0.63, -0.315, 0.0, 0.315, 0.63, 0.945, 1.26, 1.575, 1.89, 2.205, 2.52, 2.835, 3.15}; // width - 0.45 rad
  vector<double> y_bin = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

  double n_fwd_pt = fwd_pt_bin.size() - 1;
  double n_mid_pt = mid_pt_bin.size() - 1;
  double n_cos = cos_bin.size() - 1;
  double n_phi = phi_bin.size() - 1;
  double n_rapidity = y_bin.size() - 1;

  // 1d histograms
  auto y_lab_num = std::make_shared<TH1D>("y_lab_num", "", n_rapidity, &y_bin[0]);
  auto y_lab_den = std::make_shared<TH1D>("y_lab_den", "", n_rapidity, &y_bin[0]);

  // 4d histograms - cent, cos, phi, pT
  // ===== fwd ===== //
  THnD *fwd_lab_num = create_h4d("fwd_lab_num", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);
  THnD *fwd_lab_den = create_h4d("fwd_lab_den", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);

  THnD *fwd_hx_num = create_h4d("fwd_hx_num", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);
  THnD *fwd_hx_den = create_h4d("fwd_hx_den", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);

  THnD *fwd_cs_num = create_h4d("fwd_cs_num", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);
  THnD *fwd_cs_den = create_h4d("fwd_cs_den", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);

  THnD *fwd_ep_num = create_h4d("fwd_ep_num", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);
  THnD *fwd_ep_den = create_h4d("fwd_ep_den", "", cent_bin, cos_bin, phi_bin, fwd_pt_bin);

  // ===== mid ===== //
  THnD *mid_lab_num = create_h4d("mid_lab_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *mid_lab_den = create_h4d("mid_lab_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *mid_hx_num = create_h4d("mid_hx_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *mid_hx_den = create_h4d("mid_hx_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *mid_cs_num = create_h4d("mid_cs_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *mid_cs_den = create_h4d("mid_cs_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *mid_ep_num = create_h4d("mid_ep_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *mid_ep_den = create_h4d("mid_ep_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  // ===== |y| < 2.4 ===== //
  THnD *all_y_lab_num = create_h4d("all_y_lab_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *all_y_lab_den = create_h4d("all_y_lab_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *all_y_hx_num = create_h4d("all_y_hx_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *all_y_hx_den = create_h4d("all_y_hx_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *all_y_cs_num = create_h4d("all_y_cs_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *all_y_cs_den = create_h4d("all_y_cs_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  THnD *all_y_ep_num = create_h4d("all_y_ep_num", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);
  THnD *all_y_ep_den = create_h4d("all_y_ep_den", "", cent_bin, cos_bin, phi_bin, mid_pt_bin);

  // =====  make hists to save error ===== //
  y_lab_num->Sumw2();
  y_lab_den->Sumw2();

  // fwd
  fwd_lab_num->Sumw2();
  fwd_lab_den->Sumw2();
  fwd_hx_num->Sumw2();
  fwd_hx_den->Sumw2();
  fwd_cs_num->Sumw2();
  fwd_cs_den->Sumw2();
  fwd_ep_num->Sumw2();
  fwd_ep_den->Sumw2();

  // mid
  mid_lab_num->Sumw2();
  mid_lab_den->Sumw2();
  mid_hx_num->Sumw2();
  mid_hx_den->Sumw2();
  mid_cs_num->Sumw2();
  mid_cs_den->Sumw2();
  mid_ep_num->Sumw2();
  mid_ep_den->Sumw2();

  // |y| < 2.4
  all_y_lab_num->Sumw2();
  all_y_lab_den->Sumw2();
  all_y_hx_num->Sumw2();
  all_y_hx_den->Sumw2();
  all_y_cs_num->Sumw2();
  all_y_cs_den->Sumw2();
  all_y_ep_num->Sumw2();
  all_y_ep_den->Sumw2();


  // ===== connect to input tree ===== //
  // muon_chain
  const int maxBranchSize = 1000;
  Int_t Centrality;
  ULong64_t HLTriggers;
  TBranch *b_Centrality;                 //!
  TBranch *b_HLTriggers;                 //!


  //  muon id
  Short_t Reco_mu_charge[maxBranchSize];
  TBranch *b_Reco_mu_charge; //!
  muon_chain.SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);

  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;
  TClonesArray *Reco_QQ_4mom;
  TClonesArray *Reco_mu_4mom;
  ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
  ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_QQ_size]
  Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_QQ_size;                //!
  TBranch *b_Reco_mu_size;                //!
  TBranch *b_Reco_QQ_4mom;                //!
  TBranch *b_Reco_mu_4mom;                //!
  TBranch *b_Reco_QQ_trig;                //!
  TBranch *b_Reco_mu_trig;                //!
  TBranch *b_Reco_QQ_VtxProb;             //!

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  muon_chain.SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  muon_chain.SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  muon_chain.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  muon_chain.SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  muon_chain.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  muon_chain.SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  muon_chain.SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  muon_chain.SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  muon_chain.SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  TBranch *b_Reco_QQ_mupl_idx;
  TBranch *b_Reco_QQ_mumi_idx;
  muon_chain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  muon_chain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

  Float_t Reco_mu_dxy[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dxy;             //!
  muon_chain.SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  Float_t Reco_mu_dz[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dz;             //!
  muon_chain.SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nTrkWMea;           //!
  muon_chain.SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nPixWMea;           //!
  muon_chain.SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Short_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_QQ_sign;           //!
  muon_chain.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t Reco_mu_SelectionType[maxBranchSize];
  TBranch *b_Reco_mu_SelectionType;
  muon_chain.SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);


  // ===== set event plane branches =====
  const int nEP = 29; // number of event planes in the tree
  double epang[nEP];
  TBranch *b_epang;
  // ep_chain.SetBranchAddress("epang", epang, &b_epang);

  // ===== declare Lorentz vectors ===== //
  TLorentzVector *JP_Reco = nullptr;
  TLorentzVector *mupl_Reco = nullptr;
  TLorentzVector *mumi_Reco = nullptr;

  // ===== set beam information ===== //
  double sqrt_S_NN = 5.02; // TeV
  // Assuming p1 = -p2, two hadron beams are symmetric and and E >> m0
  // Mendelstam value reads sqrt_S_NN = 2 * E.
  double beam1_E = 5.02 / 2; // One beam have half of the energy
  beam1_E *= 1000;           // change TeV->GeV
  double p1_mag = beam1_E;   // Approximation for E >> m0.

  double beam2_E = beam1_E; // symmetry
  double p2_mag = -p1_mag;

  // build vectors
  TVector3 p1_3d_lab(0, 0, p1_mag);
  TVector3 p2_3d_lab(0, 0, p2_mag);

  auto *p1_lab = new TLorentzVector(p1_3d_lab, beam1_E);
  auto *p2_lab = new TLorentzVector(p2_3d_lab, beam2_E);


  // ===== dimuon counters ===== //
  long count_reco_loop = 0, count_soft_muon = 0;

  // ===== main loop ===== //
  if (nevt == -1)
    nevt = muon_chain.GetEntries();

  cout << "\nTotal Events : " << nevt << endl;

  // main loop start
  for (int iev = 0; iev < nevt; ++iev)
  {
    if (iev % 100000 == 0)
      cout << ">>>>> EVENT " << iev << " / " << muon_chain.GetEntries() << " (" << (int)(100. * iev / muon_chain.GetEntries()) << "%)" << endl;

    muon_chain.GetEntry(iev);
    // ep_chain.GetEntry(iev);

    if (Centrality > cHigh || Centrality < cLow)
      continue;

    // get rpAngles
    float epHFm2 = epang[HFm2];
    float epHFp2 = epang[HFp2];


    // ===== Reco loop strat - numerator ===== //
    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      ++count_reco_loop;

      JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      if (!(JP_Reco->M() > massLow && JP_Reco->M() < massHigh))
        continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));

      bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) &&
          passMuonTypePl //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
      );

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) &&
          passMuonTypeMi //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
      );

      if (!(muplSoft && mumiSoft))
        continue;
      ++count_soft_muon;
      if (Reco_QQ_VtxProb[irqq] < 0.01)
        continue;
      if (Reco_QQ_sign[irqq] != 0)
        continue;
      // if (Reco_mu_charge[Reco_QQ_mupl_idx[irqq]] * Reco_mu_charge[Reco_QQ_mumi_idx[irqq]] > 0)
      //   continue;

      Double_t Rapidity_reco = fabs(JP_Reco->Rapidity());
      Double_t pt_reco = fabs(JP_Reco->Pt());
      double cos_lab_ = mupl_Reco->CosTheta();
      double phi_lab_ = mupl_Reco->Phi();

      if (!(fabs(JP_Reco->Rapidity()) < 2.4 && IsAcceptable(mupl_Reco->Pt(), fabs(mupl_Reco->Eta())) && IsAcceptable(mumi_Reco->Pt(), fabs(mumi_Reco->Eta()))))
        continue;

      if (!(fabs(mupl_Reco->Eta()) < 2.4 && fabs(mumi_Reco->Eta()) < 2.4 && fabs(JP_Reco->Rapidity()) < 2.4))
        continue;

      // ===== calculate polarization angles - Reco ===== //
      // polarization unit vectors in lab
      TVector3 uy_lab(0, 1, 0);
      TVector3 uz_lab(0, 0, 1);
      TVector3 uz_ep = uy_lab;

      // HX - Upsilon
      TVector3 QQVector_Lab = JP_Reco->Vect();
      TLorentzVector MuPlusLV_QQRestFrame(*mupl_Reco);
      MuPlusLV_QQRestFrame.Boost(-JP_Reco->BoostVector());
      TVector3 MuPlusVec_Boosted = MuPlusLV_QQRestFrame.Vect();
      MuPlusVec_Boosted.RotateZ(-QQVector_Lab.Phi());
      MuPlusVec_Boosted.RotateY(-QQVector_Lab.Theta());

      float cos_hx_ = MuPlusVec_Boosted.CosTheta();
      float phi_hx_ = MuPlusVec_Boosted.Phi();


      // ===== fill the numerator ===== //
      if (pt_reco < 50 && ((Rapidity_reco > 1.6 && Rapidity_reco < 2.4 && pt_reco > 3) || (Rapidity_reco < 1.6 && pt_reco > 6.5)))
        y_lab_num->Fill(JP_Reco->Rapidity());

      // |y| < 2.4
      if (Rapidity_reco < 2.4 && pt_reco > 6.5 && pt_reco < 50)
      {
        all_y_lab_num->Fill(Centrality, cos_lab_, phi_lab_, pt_reco);
        all_y_hx_num->Fill(Centrality, cos_hx_, phi_hx_, pt_reco);
      }
      // fwd vs mid
      if (Rapidity_reco > 1.6 && Rapidity_reco < 2.4 && pt_reco > 3 && pt_reco < 50)
      {
        fwd_lab_num->Fill(Centrality, cos_lab_, phi_lab_, pt_reco);
        fwd_hx_num->Fill(Centrality, cos_hx_, phi_hx_, pt_reco);
      }
      else if (Rapidity_reco < 1.6 && pt_reco > 6.5 && pt_reco < 50)
      {
        mid_lab_num->Fill(Centrality, cos_lab_, phi_lab_, pt_reco);
        mid_hx_num->Fill(Centrality, cos_hx_, phi_hx_, pt_reco);
      }
    } // end of reco loop
  } // end of main loop

  cout << "count " << count << endl;

  // ===== save output ===== //
  TString outFileName = Form("roots/Run3_PbPb_Jpsi_noTrig_PtW%d_tnp%d_%s.root", isPtWeight, isTnP, date_label.Data());
  TFile *outFile = new TFile(outFileName, "RECREATE");

  // histograms
  y_lab_num->Write();

  // fwd
  fwd_lab_num->Write();
  fwd_hx_num->Write();
  fwd_cs_num->Write();
  fwd_ep_num->Write();

  // mid
  mid_lab_num->Write();
  mid_hx_num->Write();
  mid_cs_num->Write();
  mid_ep_num->Write();

  all_y_lab_num->Write();
  all_y_hx_num->Write();
  all_y_cs_num->Write();
  all_y_ep_num->Write();

  outFile->Close();

  cout << "\n=================================\n";
  cout << "\n Finish PbPb Eff Calculation\n";
  cout << "\n=================================\n\n";
  t->Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}

bool IsAcceptable(double pt, double eta)
{
  if ((fabs(eta) <= 1.2 && pt >= 3.5) || ((fabs(eta) >= 1.2 && fabs(eta) <= 2.1) && pt >= (5.47 - 1.89 * fabs(eta))) || ((fabs(eta) >= 2.1 && fabs(eta) <= 2.4) && pt >= 1.5))
    return true;
  //          if ((fabs(eta) < 1.0 && pt>3.4) || ((fabs(eta) >=1.0 && fabs(eta) <1.6) && pt>(5.8-2.4*fabs(eta))) || ((fabs(eta) >=1.6 && fabs(eta) <2.4) && pt>3.3667-7.0/9.0*fabs(eta))) return true;//2010 version for 2.76 TeV

  else
    return false;
}

THnD *create_h4d(const std::string &name, const std::string &title, const std::vector<double> &bins_1, const std::vector<double> &bins_2, const std::vector<double> &bins_3, const std::vector<double> &bins_4)
{
  const int n_dim = 4;

  const std::vector<const std::vector<Double_t> *> vector_bins = {&bins_1, &bins_2, &bins_3, &bins_4};

  int n_bins[n_dim];
  double x_mins[n_dim];
  double x_maxs[n_dim];

  // fill the arrays
  for (int i = 0; i < n_dim; ++i)
  {
    n_bins[i] = vector_bins[i]->size() - 1;
    x_mins[i] = vector_bins[i]->front();
    x_maxs[i] = vector_bins[i]->back();
  }

  THnD *h = new THnD(name.c_str(), title.c_str(), n_dim, n_bins, x_mins, x_maxs);

  // set the custom bin edges
  h->GetAxis(0)->Set(n_bins[0], bins_1.data());
  h->GetAxis(1)->Set(n_bins[1], bins_2.data());
  h->GetAxis(2)->Set(n_bins[2], bins_3.data());
  h->GetAxis(3)->Set(n_bins[3], bins_4.data());

  return h;
}