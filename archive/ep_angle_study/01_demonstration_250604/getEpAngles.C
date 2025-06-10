#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranchObject.h"
#include "math.h"
#include <vector>
#include <TVector3.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TParticlePDG.h>
#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include <typeinfo>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h> // gSystem
#include <cstdio> // std::remove

#define MAXCAN 10000



void getEpAngles(bool is_mc = false)
{
  // turn off the warning
  // eg. Warning in <TChain::CopyAddresses>: Could not find branch named 'Reco_trk_4mom_pt' in tree named 'tempTree'
  // The reasons is that tempTree doesn't have 3 branches which myChain has. But it was intended so it's not a matter
  // you can do test run while commenting it out
  gErrorIgnoreLevel = kError;

  // ===== many input files ===== //
  // use this part when you use unmerged root files in lxplus
  TString base_name = "/disk1/Oniatree/polarization/oniatree_5p36/2023_PbPb/OniaTree_RawPrime0_Run3_PbPb_track_vector_250530/HIPhysicsRawPrime0/crab_OniaTree_RawPrime0_Run3_PbPb_track_vector_250530/250530_135847/0000/Oniatree_2023PbPbPromptRecoData_132X_miniAOD";

  TString nameOldTree = "hionia/myTree";
  TChain myChain(nameOldTree);

  // add files into the chain
  for (int i = 0; i <= 600; ++i)
  {
    TString filename = TString::Format("%s_%d.root", base_name.Data(), i);
    if (!gSystem->AccessPathName(filename, kFileExists)) {
      myChain.Add(filename);
      cout << "Input added: " << filename << endl;
    }
      
  }

  // ===== get one input file ===== //
  // TString nameOldTree = "hionia/myTree";
  // TChain myChain(nameOldTree);
  // myChain.Add("Oniatree_2023PbPbPromptRecoData_132X_miniAOD.root");


  // ===== branches ===== //
  // variables in input tree. They will be used for calculatin and will be saved
  Short_t Reco_QQ_size;
  Short_t Reco_QQ_sign[17]; //[Reco_QQ_size]
  TClonesArray *Reco_QQ_4mom = nullptr;
  Short_t Reco_QQ_mupl_idx[1000];
  Short_t Reco_QQ_mumi_idx[1000];
  TClonesArray *Reco_mu_4mom = nullptr;
  Short_t Reco_trk_size;

  // variables in input tree. They will be used for calculatin but will not saved
  vector<float> *Reco_trk_4mom_pt = nullptr;
  vector<float> *Reco_trk_4mom_eta = nullptr;
  vector<float> *Reco_trk_4mom_phi = nullptr;

  // new variables
  Float_t Q2x_tree, Q2y_tree, ep, ep_rec, ep_flat;


  // ===== make output file ===== //
  TFile *outTemp = TFile::Open("out_buffer_very_long_strange_file_name_never_be_used.root", "RECREATE");
  outTemp->cd();

  // ===== set branches from old tree ===== //
  myChain.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  myChain.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  myChain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  myChain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  myChain.SetBranchAddress("Reco_trk_size", &Reco_trk_size);

  // TClonesArray Branches
  myChain.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  myChain.SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // std::vector<float>* Branches
  myChain.SetBranchAddress("Reco_trk_4mom_pt", &Reco_trk_4mom_pt);
  myChain.SetBranchAddress("Reco_trk_4mom_eta", &Reco_trk_4mom_eta);
  myChain.SetBranchAddress("Reco_trk_4mom_phi", &Reco_trk_4mom_phi);

  std::set<std::string> dropBranches = {"Reco_trk_4mom_pt", "Reco_trk_4mom_eta", "Reco_trk_4mom_phi"};
  
  // turn off branches before cloning
  for (const std::string &branch : dropBranches)
    if (myChain.GetBranch(branch.c_str()))
      myChain.SetBranchStatus(branch.c_str(), 0);

  // ===== set new TTree ===== //
  TTree *tempTree = myChain.CloneTree(0);
  tempTree->SetName("tempTree");

  tempTree->Branch("Q2x_tree", &Q2x_tree, "Q2x_tree/F");
  tempTree->Branch("Q2y_tree", &Q2y_tree, "Q2y_tree/F");
  tempTree->Branch("ep", &ep, "ep/F");

  // turn on branches to use in the loop
  for (const std::string &branch : dropBranches)
  {
    if (myChain.GetBranch(branch.c_str()))
    {
      myChain.SetBranchStatus(branch.c_str(), 1);
    }
  }

  // ===== additonal variables ===== //
  TH2F *hQ2xy = new TH2F("hQ2xy", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);
  
  float d1Phi[MAXCAN], d1Eta[MAXCAN];
  float d2Phi[MAXCAN], d2Eta[MAXCAN];

  // ===== main loop 1 - Qx, Qy ===== //
  Long64_t nEvent = myChain.GetEntries();
  // nEvent = 100000;
  cout << "Total events: " << nEvent << "\n\n";
  cout << "Start loop 1\n\n";
  for (Long64_t ievt = 0; ievt < nEvent; ++ievt)
  {
    if (ievt % 100000 == 0 && ievt != 0)
      std::cout << "Event " << ievt << "\n";

    myChain.GetEntry(ievt);
    tempTree->GetEntry(ievt);

    int N = 0;
    int nsize = Reco_QQ_size;
    // cout << ievt << endl ;
    float Q2x_num = 0;
    float Q2y_num = 0;
    float Qw_pt_num = 0;

    // ===== get Qx, Qy ===== //
    for (int iQQ = 0; iQQ < nsize; iQQ++)
    {
      if (Reco_QQ_sign[iQQ] != 0)
        continue; // only opposite-sign muon pairs
      //  if (Reco_QQ_VtxProb[iQQ] < 0.01) continue; // good common vertex proba
      TLorentzVector *qq4mom = (TLorentzVector *)Reco_QQ_4mom->At(iQQ);
      if (fabs(qq4mom->Eta()) > 2.4)
        continue;
      int iMuPlus = (Reco_QQ_mupl_idx)[iQQ];
      int iMuMinus = (Reco_QQ_mumi_idx)[iQQ];
      TLorentzVector *Reco_mupl_4mom = (TLorentzVector *)Reco_mu_4mom->At(iMuPlus);
      TLorentzVector *Reco_mumi_4mom = (TLorentzVector *)Reco_mu_4mom->At(iMuMinus);

      d1Phi[N] = Reco_mupl_4mom->Phi();
      d1Eta[N] = Reco_mupl_4mom->Eta();
      d2Phi[N] = Reco_mumi_4mom->Phi();
      d2Eta[N] = Reco_mumi_4mom->Eta();
      N++;
    }

    // ===== calculate EP from Qx, Qy ===== //
    nsize = Reco_trk_size;
    int Noff = 0;
    for (unsigned it = 0; it < nsize; it++)
    {
      // TLorentzVector *trk4mom = (TLorentzVector *)Reco_trk_4mom->At(it);
      float eta = Reco_trk_4mom_eta->at(it);
      float pt = Reco_trk_4mom_pt->at(it);
      float phi = Reco_trk_4mom_phi->at(it);

      if (fabs(eta) <= 2.4 && pt > 0.4)
        Noff++;

      bool reject = false;
      for (unsigned i = 0; i < N; ++i)
      {
        if (fabs(eta - d1Eta[i]) < 0.03 && fabs(phi - d1Phi[i]) < 0.03)
          reject = true;
        if (fabs(eta - d2Eta[i]) < 0.03 && fabs(phi - d2Phi[i]) < 0.03)
          reject = true;
      }
      if (true == reject)
      {
        // cout << "rejected" << endl ;
        continue;
      }

      if (fabs(eta) <= 2.4 && pt > 0.3 && pt < 3)
      {
        Q2x_num += pt * cos(2 * phi);
        Q2y_num += pt * sin(2 * phi);
        Qw_pt_num += pt;
      }
    }

    // cout << Noff << endl ;
    Q2x_tree = Q2x_num / Qw_pt_num;
    Q2y_tree = Q2y_num / Qw_pt_num;

    float ev = 0.5 * atan2(Q2y_tree, Q2x_tree);
    ep = ev;
    // cout << Q2x_tree << "  " << Q2y_tree << endl ;
    hQ2xy->Fill(Q2x_tree, Q2y_tree);
    
    tempTree->Fill();
  }
  cout << "loop 1 done (Qx, Qy)" << endl;


  // ===== main loop 2 - Recentering ===== //
  // mean Q2x, Q2y
  float mean_Q2x = hQ2xy->GetMean(1);
  float mean_Q2y = hQ2xy->GetMean(2);
  cout << mean_Q2x << "  " << mean_Q2y << endl;

  tempTree->SetBranchAddress("Q2x_tree", &Q2x_tree);
  tempTree->SetBranchAddress("Q2y_tree", &Q2y_tree);

  TTree *recTree = tempTree->CloneTree(0);
  recTree->SetName("recTree");
  recTree->Branch("ep_rec", &ep_rec, "ep_rec/F");
  
  const int N = 10;   // Number of histograms
  TH1F *means_sin[N]; // Array of histogram pointers
  TH1F *means_cos[N]; // Array of histogram pointers
  
  // recentering loop
  for (int i = 0; i < N; i++)
  {
    means_sin[i] = new TH1F(Form("means_sin_%d", i), Form("Histogram %d", i), 20, -1, 1);
    means_cos[i] = new TH1F(Form("means_cos_%d", i), Form("Histogram %d", i), 20, -1, 1);
  }

  cout << "\n\nStart loop 2\n\n";
  for (Int_t i = 0; i < nEvent; i++)
  {
    if (i % 100000 == 0 && i != 0)
      std::cout << "Event " << i << "\n";
    tempTree->GetEntry(i);
    recTree->GetEntry(i);
    float ep1 = 0.5 * atan2(Q2y_tree - mean_Q2y, Q2x_tree - mean_Q2x);
    ep_rec = ep1;
    float j = 1.;

    for (int i = 0; i < N; i++)
    {
      means_sin[i]->Fill(sin(2 * j * ep_rec));
      means_cos[i]->Fill(cos(2 * j * ep_rec));
      j += 1;
    }
    recTree->Fill();
  }
  cout << "loop 2 done" << endl;


  // ===== main loop 3 - Flattening ===== //
  TFile *outFile = TFile::Open("PbPb_ep_angle.root", "RECREATE");
  outFile->cd();
  TDirectory *hioniaDir = outFile->mkdir("hionia");
  outFile->cd("hionia");

  TTree *myTree = recTree->CloneTree(0);
  myTree->SetName("myTree");
  myTree->Branch("ep_flat", &ep_flat, "ep_flat/F");
  recTree->SetBranchAddress("ep_rec", &ep_rec);

  // mean sin, cos
  double mean_sin[N];
  double mean_cos[N];

  cout << "\n\nStart loop 3\n\n";
  for (int i = 0; i < N; i++)
  {
    mean_sin[i] = means_sin[i]->GetMean();
    mean_cos[i] = means_cos[i]->GetMean();
  }

  for (Int_t i = 0; i < nEvent; i++)
  {
    if (i % 100000 == 0 && i != 0)
      std::cout << "Event " << i << "\n";

    recTree->GetEntry(i);
    // cout << i << endl ;
    float psi = 0; // Start with ep
    int n_max = 10;
    for (int n = 0; n < n_max; n++)
    {
      double mean_sin_n = mean_sin[n];
      double mean_cos_n = mean_cos[n];
      double term = (2.0 / (n + 1)) * (-mean_sin_n * cos((2. * (n + 1)) * ep_rec) + mean_cos_n * sin((2. * (n + 1)) * ep_rec));

      //  cout << ep << "  " << term << endl ;
      psi += term;
    }
    //   cout << ep_rec << endl ;
    psi = 1. / 2. * psi + ep_rec;
    //   cout << ep_rec << "  " << psi << endl;
    ep_flat = psi;
    myTree->Fill();
  }
  cout << "loop 3 done" << endl;

  // ===== saving and cleaning ===== //
  myTree->Write();
  outFile->cd();
  hQ2xy->Write();

  outTemp->Close();
  outFile->Close();

  // remove buffer output file
  string bufferFileName = "out_buffer_very_long_strange_file_name_never_be_used.root";
  std::remove(bufferFileName.c_str());
}