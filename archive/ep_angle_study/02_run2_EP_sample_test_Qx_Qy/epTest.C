#include <iostream>
#include "TH2D.h"
#include <TLorentzVector.h>
#include <TAttMarker.h>
#include "../../../headers/commonUtility.h"
#include "../../../headers/JpsiUtility.h"
#include "../../../headers/HiEvtPlaneList.h"
#include "../../../headers/cutsAndBin.h"
#include "../../../headers/Style.h"
#include "../../../headers/tnp_weight_lowptPbPb_num_den_new.h"
// #include "../headers/tnp_weight_lowptPbPb.h"

void divide_2d_hist(const std::shared_ptr<TH2D> &h_num, const std::shared_ptr<TH2D> &h_den,
                    std::shared_ptr<TH2D> &h_eff);
void draw_2d_hist(const std::shared_ptr<TH2D> &h_eff, string legend, int MCtype);
bool IsAcceptable(double pt, double eta);

void epTest(string data_label_ = "test3", int nevt = -1, bool isMC = true, bool isPeri = false)
{
  using namespace std;
  using namespace hi; // HiEvtPlaneList.h

  // Example of using event plane namespace
  cout << " Index of " << EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of " << EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of " << EPNames[trackmid2] << " = " << trackmid2 << endl;

  // labeling
  TString date_label = data_label_;

  // ===== import oniatree ===== //
  TChain muon_chain("myTree");
  TChain ep_chain("tree");

  if (isMC == true)
    for (int i = 1; i <= 2; ++i)
    {
      // PbPb MC prompt
      TString filename = Form("/disk1/Oniatree/Jpsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part%d.root", i);
      std::cout << "Adding PbPb MC prompt sample: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }
  else if (isMC == false && isPeri == false)
    for (int i = 1; i <= 5; ++i)
    {
      // Run2 Central collison
      TString filename = Form("/disk1/Oniatree/polarization/run2_oniatree_PbPb_data_AOD/DM/ReReco_Oniatree_addvn_part%d.root", i);
      std::cout << "Adding PbPb MC prompt sample: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }
  else
    for (int i = 1; i <= 3; ++i)
    {
      // Run2 Central collison
      TString filename = Form("/disk1/Oniatree/polarization/run2_oniatree_PbPb_data_AOD/DMPeri/ReReco_Oniatree_addvn_part%d.root", i);
      std::cout << "Adding PbPb MC prompt sample: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }

  // ===== build histograms ===== //
  TH1D *h_epAngHFm2 = new TH1D("h_epAngHFm2", "", 100, -2, 2);
  TH1D *h_epAngHFp2 = new TH1D("h_epAngHFp2", "", 100, -2, 2);
  TH1D *h_epAngQxQyHFm2 = new TH1D("h_epAngQxQyHFm2", "", 100, -2, 2);
  TH1D *h_epAngQxQyHFp2 = new TH1D("h_epAngQxQyHFp2", "", 100, -2, 2);

  h_epAngHFm2->Sumw2();
  h_epAngHFp2->Sumw2();
  h_epAngQxQyHFm2->Sumw2();
  h_epAngQxQyHFp2->Sumw2();

  // ===== connect to input tree ===== //
  // muon_chain
  const int maxBranchSize = 1000;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Int_t Gen_QQ_size;
  Int_t Gen_mu_size;
  TClonesArray *Gen_QQ_4mom;
  TClonesArray *Gen_mu_4mom;
  ULong64_t Gen_QQ_trig[maxBranchSize];  //[Gen_QQ_size]
  Float_t Gen_QQ_VtxProb[maxBranchSize]; //[Gen_QQ_size]
  TBranch *b_Centrality;                 //!
  TBranch *b_HLTriggers;                 //!
  TBranch *b_Gen_QQ_size;                //!
  TBranch *b_Gen_mu_size;                //!
  TBranch *b_Gen_QQ_4mom;                //!
  TBranch *b_Gen_mu_4mom;                //!
  TBranch *b_Gen_QQ_trig;                //!
  TBranch *b_Gen_QQ_VtxProb;             //!


  Int_t Reco_QQ_size;
  Int_t Reco_mu_size;
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

  Int_t Reco_QQ_mupl_idx[maxBranchSize];
  Int_t Reco_QQ_mumi_idx[maxBranchSize];
  TBranch *b_Reco_QQ_mupl_idx;
  TBranch *b_Reco_QQ_mumi_idx;
  muon_chain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  muon_chain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

  // ===== set event plane branches =====
  const int nEP = 29; // number of event planes in the tree
  double epang[nEP];
  TBranch *b_epang;
  ep_chain.SetBranchAddress("epang", epang, &b_epang);

  double qx[nEP];
  TBranch *b_qx;
  ep_chain.SetBranchAddress("qx", qx, &b_qx);

  double qy[nEP];
  TBranch *b_qy;
  ep_chain.SetBranchAddress("qy", qy, &b_qy);

  TLorentzVector *JP_Reco = nullptr;
  TLorentzVector *mupl_Reco = nullptr;
  TLorentzVector *mumi_Reco = nullptr;

  // ===== prepare weighting variables ===== //
  if (nevt == -1)
    nevt = muon_chain.GetEntries();

  cout << "\nTotal Events : " << nevt << endl;

  // ===== main loop start ===== //
  for (int iev = 0; iev < nevt; ++iev)
  {
    if (iev % 100000 == 0)
      cout << ">>>>> EVENT " << iev << " / " << muon_chain.GetEntries() << " (" << (int)(100. * iev / muon_chain.GetEntries()) << "%)" << endl;

    muon_chain.GetEntry(iev);
    ep_chain.GetEntry(iev);

    // get rpAngles
    float epHFm2 = epang[HFm2];
    float epHFp2 = epang[HFp2];

    float qxHFm2 = qx[HFm2];
    float qxHFp2 = qx[HFp2];

    float qyHFm2 = qy[HFm2];
    float qyHFp2 = qy[HFp2];


    // cout << "epAng: " << epAng << endl;
    // cout << "epAng (Qx, Qy): " << (double)1 / 2 * TMath::ATan2(qy_, qx_) << endl;

    h_epAngHFm2->Fill(epHFm2);
    h_epAngHFp2->Fill(epHFp2);
    h_epAngQxQyHFm2->Fill((double)1 / 2 * TMath::ATan2(qyHFm2, qxHFm2));
    h_epAngQxQyHFp2->Fill((double)1 / 2 * TMath::ATan2(qyHFp2, qxHFp2));
  } // end of main loop

  string outname = data_label_ + ".root";

  TFile outfile(outname.c_str(), "recreate");
  outfile.cd();
  h_epAngHFm2->Write();
  h_epAngHFp2->Write();
  h_epAngQxQyHFm2->Write();
  h_epAngQxQyHFp2->Write();
  outfile.Close();
}