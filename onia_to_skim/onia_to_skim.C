#include <ctime>
#include <TLorentzVector.h>
#include "../headers/commonUtility.h"
#include "../headers/HiEvtPlaneList.h"
#include "../headers/cutsAndBin.h"
#include "../headers/tnp_weight_lowptPbPb.h"

static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D *h = 0, double pt = 0);
double getEffWeight(TH1D *h = 0, double pt = 0);

void onia_to_skim(string data_label_ = "test", int nevt = 1000, bool isMC = true, int MCtype = 2, int PDtype = 2, int hiHFBinEdge = 0, int kTrigSel = kTrigJpsi)
{
  TStopwatch *t = new TStopwatch;
  t->Start();
  cout << "\n=================================\n";
  cout << "\n Start Oniatree skimming\n";
  cout << "\n=================================\n";

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace
  cout << " Index of " << EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of " << EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of " << EPNames[trackmid2] << " = " << trackmid2 << endl;

  // labeling
  TString date_label = data_label_;

  TString fPD;
  if (PDtype == 1)
    fPD = "DB";
  else if (PDtype == 2)
    fPD = "DBPeri";

  TString fMCtype;
  if (MCtype == 1)
    fMCtype = "Signal";
  else if (MCtype == 2)
    fMCtype = "NPOnly";

  int trigIndx = 0;
  if (kTrigSel == kTrigJpsi)
    trigIndx = 0;
  else if (kTrigSel == kTrigUps)
    trigIndx = 1;
  else if (kTrigSel == kTrigL1DBOS40100)
    trigIndx = 2;
  else if (kTrigSel == kTrigL1DB50100)
    trigIndx = 3;
  
  int kL2filter = 38;
  int kL3filter = 39;

  TString fCentSelHF = "HFNom";
  if (hiHFBinEdge == 1)
    fCentSelHF = "HFUp";
  else if (hiHFBinEdge == -1)
    fCentSelHF = "HFDown";
  

  // ===== import oniatree ===== //
  TChain muon_chain("myTree");
  TChain ep_chain("tree");

  if (isMC && MCtype==1)
    for (int i = 1; i <= 2; ++i)
    {
      // PbPb MC prompt
      TString filename = Form("/disk1/Oniatree/Jpsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part%d.root", i);
      std::cout << "Adding PbPb MC prompt sample: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }
  else if (isMC && MCtype==2)
  {
    // PbPb MC nonprompt
    TString filename = "/disk1/Oniatree/Jpsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root";
    std::cout << "Adding PbPb MC nonprompt sample: " << filename << std::endl;

    muon_chain.Add(filename);
    ep_chain.Add(filename);
  }
  else if (!isMC && PDtype==1)
    for (int i = 1; i <= 5; ++i)
    {
      // PbPb double muon data
      TString filename = Form("/work/pjgwak/oniatree_PbPb_data_AOD/DM/ReReco_Oniatree_addvn_part%d.root", i);
      std::cout << "Adding PbPb double muon data: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }
  else if (!isMC && PDtype==2)
    for (int i = 1; i <= 3; ++i)
    {
      // PbPb double muon peripheral data
      TString filename = Form("/work/pjgwak/oniatree_PbPb_data_AOD/DMPeri/ReReco_Oniatree_addvn_part%d.root", i);
      std::cout << "Adding PbPb double muon peripheral data: " << filename << std::endl;

      muon_chain.Add(filename);
      ep_chain.Add(filename);
    }


  // ===== set dimuon branches =====
  const int maxBranchSize = 1000;

  UInt_t runNb;
  UInt_t eventNb, LS;
  float zVtx;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Float_t SumET_HF;
  Int_t Reco_QQ_size;
  Int_t Reco_mu_size;
  TClonesArray *Reco_QQ_4mom = nullptr;
  TClonesArray *Reco_mu_4mom = nullptr;
  ULong64_t Reco_QQ_trig[maxBranchSize];  //[Reco_QQ_size]
  ULong64_t Reco_mu_trig[maxBranchSize];  //[Reco_QQ_size]
  Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_runNb;                       //!
  TBranch *b_eventNb;                     //!
  TBranch *b_LS;
  TBranch *b_zVtx;            //!
  TBranch *b_Centrality;      //!
  TBranch *b_HLTriggers;      //!
  TBranch *b_SumET_HF;        //!
  TBranch *b_Reco_QQ_size;    //!
  TBranch *b_Reco_mu_size;    //!
  TBranch *b_Reco_QQ_4mom;    //!
  TBranch *b_Reco_mu_4mom;    //!
  TBranch *b_Reco_QQ_trig;    //!
  TBranch *b_Reco_mu_trig;    //!
  TBranch *b_Reco_QQ_VtxProb; //!

  Bool_t Reco_mu_highPurity[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_mu_highPurity;            //!
  muon_chain.SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  
  if (!isMC)
  {
    muon_chain.SetBranchAddress("runNb", &runNb, &b_runNb);
    muon_chain.SetBranchAddress("LS", &LS, &b_LS);
  }
  
  muon_chain.SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  muon_chain.SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  muon_chain.SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  muon_chain.SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  muon_chain.SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  muon_chain.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  muon_chain.SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  muon_chain.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  muon_chain.SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  muon_chain.SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  muon_chain.SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  muon_chain.SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id
  Int_t Reco_QQ_mupl_idx[maxBranchSize];
  Int_t Reco_QQ_mumi_idx[maxBranchSize];
  TBranch *b_Reco_QQ_mupl_idx;
  TBranch *b_Reco_QQ_mumi_idx;
  muon_chain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  muon_chain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

  Int_t Reco_mu_nTrkHits[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nTrkHits;           //!
  muon_chain.SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t Reco_mu_normChi2_global[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_normChi2_global;             //!
  muon_chain.SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t Reco_mu_nMuValHits[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nMuValHits;           //!
  muon_chain.SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t Reco_mu_StationsMatched[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_StationsMatched;           //!
  muon_chain.SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t Reco_mu_dxy[maxBranchSize];    //[Reco_mu_size]
  Float_t Reco_mu_dxyErr[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dxy;                //!
  TBranch *b_Reco_mu_dxyErr;             //!
  muon_chain.SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  muon_chain.SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t Reco_mu_dz[maxBranchSize];    //[Reco_mu_size]
  Float_t Reco_mu_dzErr[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_dz;                //!
  TBranch *b_Reco_mu_dzErr;             //!
  muon_chain.SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  muon_chain.SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nTrkWMea;           //!
  muon_chain.SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t Reco_mu_TMOneStaTight[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_TMOneStaTight;            //!

  muon_chain.SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
  TBranch *b_Reco_mu_nPixWMea;           //!
  muon_chain.SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_QQ_sign;           //!
  muon_chain.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  
  // use ep_chain
  // Float_t rpAng[29]; //[nEP]
  // TBranch *b_rpAng;  //!
  //  muon_chain.SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t Reco_mu_nPixValHits[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_mu_nPixValHits;           //!
  muon_chain.SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t Reco_mu_ptErr_global[maxBranchSize]; //[Reco_QQ_size]
  TBranch *b_Reco_mu_ptErr_global;             //!
  muon_chain.SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t Reco_mu_SelectionType[maxBranchSize];
  TBranch *b_Reco_mu_SelectionType;
  muon_chain.SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctauErr3D[maxBranchSize];
  TBranch *b_Reco_QQ_ctau3D;
  TBranch *b_Reco_QQ_ctauErr3D;
  muon_chain.SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  muon_chain.SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);

  // MC only
  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if (isMC)
  {
    muon_chain.SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
    muon_chain.SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
  }


  // ===== set event plane branches =====
  const int nEP = 29; // number of event planes in the tree
  double epang[nEP];
  TBranch *b_epang;
  ep_chain.SetBranchAddress("epang", epang, &b_epang);

  
  // ===== make output file ===== //
  TFile *newfile;
  // well, I keep the name "Flow"Skim as legacy. you can change it.

  if (isMC)
    newfile = new TFile(Form("../files_skim/OniaFlowSkim_%s_AOD_isMC%d_%s_%s_%s.root", fTrigName[trigIndx].Data(), isMC, fMCtype.Data(), fCentSelHF.Data(), date_label.Data()), "recreate");
  else
    newfile = new TFile(Form("../files_skim/OniaFlowSkim_%sTrig_%sPD_AOD_isMC%d_%s_%s.root", fTrigName[trigIndx].Data(), fPD.Data(), isMC, fCentSelHF.Data(), date_label.Data()), "recreate");


  TTree *mmevttree = new TTree("mmepevt", "dimuonAndEventPlanes in event based");

  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nDimu;
  float vz;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float y[nMaxDimu];
  float phi[nMaxDimu];
  float phi1[nMaxDimu];
  float phi2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float weight0[nMaxDimu];
  float weight1[nMaxDimu];
  int recoQQsign[nMaxDimu];
  float ctau3D[nMaxDimu];
  float ctau3DErr[nMaxDimu];
  double TnPweight[nMaxDimu] = {1.};
  double weight = 1;

  // for polarization
  float cos_theta[nMaxDimu];
  float cos_theta1[nMaxDimu];
  float cos_cs[nMaxDimu];
  float phi_cs[nMaxDimu];
  float cos_hx[nMaxDimu];
  float phi_hx[nMaxDimu];
  float cos_px[nMaxDimu];
  float phi_px[nMaxDimu];
  float cos_ep[nMaxDimu];
  float phi_ep[nMaxDimu];
  float cos_lab[nMaxDimu];
  float phi_lab[nMaxDimu];

  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event", &evt, "event/I");
  mmevttree->Branch("runN", &runN, "runN/I");
  mmevttree->Branch("lumi", &lumi, "lumi/I");
  mmevttree->Branch("cBin", &cBin, "cBin/I");
  mmevttree->Branch("vz", &vz, "vz/F");
  mmevttree->Branch("nDimu", &nDimu, "nDimu/I");
  mmevttree->Branch("mass", mass, "mass[nDimu]/F");
  mmevttree->Branch("y", y, "y[nDimu]/F");
  mmevttree->Branch("pt", pt, "pt[nDimu]/F");
  mmevttree->Branch("pt1", pt1, "pt1[nDimu]/F");
  mmevttree->Branch("pt2", pt2, "pt2[nDimu]/F");
  mmevttree->Branch("phi", phi, "phi[nDimu]/F");
  mmevttree->Branch("phi1", phi1, "phi1[nDimu]/F");
  mmevttree->Branch("phi2", phi2, "phi2[nDimu]/F");
  mmevttree->Branch("eta", eta, "eta[nDimu]/F");
  mmevttree->Branch("eta1", eta1, "eta1[nDimu]/F");
  mmevttree->Branch("eta2", eta2, "eta2[nDimu]/F");
  mmevttree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  mmevttree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  mmevttree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  mmevttree->Branch("weight", &weight, "weight/D");
  mmevttree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");

  // for polarization
  mmevttree->Branch("cos_theta", cos_theta, "cos_theta[nDimu]/F");
  mmevttree->Branch("cos_theta1", cos_theta1, "cos_theta1[nDimu]/F");
  mmevttree->Branch("cos_cs", cos_cs, "cos_cs[nDimu]/F");
  mmevttree->Branch("phi_cs", phi_cs, "phi_cs[nDimu]/F");
  mmevttree->Branch("cos_hx", cos_hx, "cos_hx[nDimu]/F");
  mmevttree->Branch("phi_hx", phi_hx, "phi_hx[nDimu]/F");
  mmevttree->Branch("cos_px", cos_px, "cos_px[nDimu]/F");
  mmevttree->Branch("phi_px", phi_px, "phi_px[nDimu]/F");
  mmevttree->Branch("cos_ep", cos_ep, "cos_ep[nDimu]/F");
  mmevttree->Branch("phi_ep", phi_ep, "phi_ep[nDimu]/F");

  // ===== declare Lorentz vectors ===== //
  TLorentzVector *JP_Reco = new TLorentzVector;
  TLorentzVector *mupl_Reco = new TLorentzVector;
  TLorentzVector *mumi_Reco = new TLorentzVector;


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


  // ===== tnp weight variables and  dimuon counters ===== //
  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;

  int count = 0;
  int counttnp = 0;


  // ===== event loop start ===== //
  if (nevt == -1)
    nevt = muon_chain.GetEntries();

  cout << "Total events = " << nevt << ", : " << ep_chain.GetEntries() << endl;

  for (int iev = 0; iev < nevt; ++iev)
  {
    if (iev % 100000 == 0)
      cout << ">>>>> EVENT " << iev << " / " << muon_chain.GetEntries() << " (" << (int)(100. * iev / muon_chain.GetEntries()) << "%)" << endl;

    muon_chain.GetEntry(iev);
    ep_chain.GetEntry(iev);

    nDimu = 0;

    if (PDtype == 1 && runNb >= 327123)
      continue;

    if (!((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
      continue;

    // polarization unit vectors in lab
    TVector3 uy_lab(0, 1, 0);
    TVector3 uz_lab(0, 0, 1);
    TVector3 uz_ep = uy_lab;

    // get rpAngles
    float epHFm2 = epang[HFm2];
    float epHFp2 = epang[HFp2];


    // ===== start dimuon loop ===== //
    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      runN = runNb;
      evt = eventNb;
      lumi = LS;
      cBin = -999;
      if (hiHFBinEdge == 0)
        cBin = getHiBinFromhiHF(SumET_HF);
      else if (hiHFBinEdge == 1)
        cBin = getHiBinFromhiHF_Up(SumET_HF);
      else if (hiHFBinEdge == -1)
        cBin = getHiBinFromhiHF_Down(SumET_HF);
      if (cBin == -999)
      {
        cout << "ERROR!!! No HF Centrality Matching!!" << endl;
        return;
      }
      vz = zVtx;

      JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      weight = 1.;
      if (isMC)
        weight = findNcoll(Centrality) * Gen_weight;

      if (!((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
        continue;

      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
          continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
          continue;
      }

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

      if (Reco_QQ_VtxProb[irqq] < 0.01)
        continue;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];

      count++;
      if (isMC)
      {
        tnp_weight = 1;
        tnp_trig_weight_mupl = -1;
        tnp_trig_weight_mumi = -1;
        tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); // mu id
        tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0);                                     // inner tracker

        // Trigger part
        if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
        {
          //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
          //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
          continue;
        }
        bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
        bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
        bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
        bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
        if (mupl_L2Filter == false || mumi_L2Filter == false)
        {
          cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl;
          cout << endl;
        }

        bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
        bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
        bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
        bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
        bool SelDone = false;

        if (mupl_isL2 && mumi_isL3)
        {
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
          SelDone = true;
        }
        else if (mupl_isL3 && mumi_isL2)
        {
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
          SelDone = true;
        }
        else if (mupl_isL3 && mumi_isL3)
        {
          int t[2] = {-1, 1}; // mupl, mumi
          int l = rand() % (2);
          // pick up what will be L2
          if (t[l] == -1)
          {
            tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
            tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
          }
          else if (t[l] == 1)
          {
            tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
            tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
          }
          else
          {
            cout << "ERROR :: No random selection done !!!!" << endl;
            continue;
          }
          SelDone = true;
        }
        if (SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1))
        {
          cout << "ERROR :: No muon filter combination selected !!!!" << endl;
          continue;
        }
        tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
        counttnp++;
      } // end of MC tnp loop


      // ===== calculate polarization angles
      // consider auto correleation
      float eta_ = JP_Reco->Eta();
      float epAng = 0;
      if (eta_ >= 0)
        epAng = epHFm2;
      else
        epAng = epHFp2;

      // get z-axis in EP frame
      uz_ep.Rotate(epAng, uz_lab); // uz_ep = y rotated as epAng in lab frame

      // Boost to Quarkonia rest frame
      std::unique_ptr<TLorentzVector> mupl_qrf = std::make_unique<TLorentzVector>(*mupl_Reco);
      mupl_qrf->Boost(-JP_Reco->BoostVector());
      std::unique_ptr<TLorentzVector> mumi_qrf = std::make_unique<TLorentzVector>(*mumi_Reco);
      mumi_qrf->Boost(-JP_Reco->BoostVector());
      std::unique_ptr<TLorentzVector> p1_qrf = std::make_unique<TLorentzVector>(*p1_lab);
      p1_qrf->Boost(-JP_Reco->BoostVector());
      std::unique_ptr<TLorentzVector> p2_qrf = std::make_unique<TLorentzVector>(*p2_lab);
      p2_qrf->Boost(-JP_Reco->BoostVector());

      // u_p1, u_p2
      TVector3 u_p1 = p1_qrf->Vect();
      u_p1 = u_p1.Unit();
      TVector3 u_p2 = p2_qrf->Vect();
      u_p2 = u_p2.Unit();

      // unit vector of mupl
      TVector3 u_mupl = mupl_qrf->Vect(); // get px, py, pz
      u_mupl = u_mupl.Unit();

      // y axis for frames except EP frame
      TVector3 uy_pol = u_p1.Cross(u_p2); // p1 cross p2
      uy_pol = uy_pol.Unit();

      // CS
      TVector3 uz_cs = u_p1 - u_p2;
      uz_cs = uz_cs.Unit();
      TVector3 ux_cs = uy_pol.Cross(uz_cs);
      ux_cs = ux_cs.Unit();
      float cos_cs_ = u_mupl.Dot(uz_cs);
      float phi_cs_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_cs));

      // HX
      TVector3 uz_hx = JP_Reco->Vect();
      uz_hx = uz_hx.Unit();
      TVector3 ux_hx = uy_pol.Cross(uz_hx); // y cross z
      ux_hx = ux_hx.Unit();
      float cos_hx_ = u_mupl.Dot(uz_hx);
      float phi_hx_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_hx));

      // PX
      TVector3 uz_px = u_p1 + u_p2;
      uz_px = uz_px.Unit();
      TVector3 ux_px = uy_pol.Cross(uz_px);
      ux_px = ux_px.Unit();
      float cos_px_ = u_mupl.Dot(uz_px);
      float phi_px_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_px));

      // EP
      TVector3 uy_ep = p1_3d_lab;
      uy_ep = uy_ep.Unit();
      TVector3 ux_ep = uy_ep.Cross(uz_ep);
      ux_ep = ux_ep.Unit();
      float cos_ep_ = u_mupl.Dot(uz_ep);
      float phi_ep_ = TMath::ATan2(u_mupl.Dot(uy_pol), u_mupl.Dot(ux_ep));

      // ===== Fill the output tree
      if (isMC)
        TnPweight[nDimu] = tnp_weight;
      mass[nDimu] = JP_Reco->M();
      phi[nDimu] = JP_Reco->Phi();
      phi1[nDimu] = mupl_Reco->Phi();
      phi2[nDimu] = mumi_Reco->Phi();
      y[nDimu] = JP_Reco->Rapidity();
      pt[nDimu] = JP_Reco->Pt();
      pt1[nDimu] = mupl_Reco->Pt();
      pt2[nDimu] = mumi_Reco->Pt();
      eta[nDimu] = JP_Reco->Eta();
      eta1[nDimu] = mupl_Reco->Eta();
      eta2[nDimu] = mumi_Reco->Eta();
      cos_theta[nDimu] = JP_Reco->CosTheta();
      cos_theta1[nDimu] = mupl_Reco->CosTheta();
      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      cos_cs[nDimu] = cos_cs_;
      phi_cs[nDimu] = phi_cs_;
      cos_hx[nDimu] = cos_hx_;
      phi_hx[nDimu] = phi_hx_;
      cos_px[nDimu] = cos_px_;
      phi_px[nDimu] = phi_px_;
      cos_ep[nDimu] = cos_ep_;
      phi_ep[nDimu] = phi_ep_;
      nDimu++;

    } // end of dimuon loop
    if (nDimu > 0)
      mmevttree->Fill();

  } // end of event loop
  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;


  // ===== save output ===== //
  newfile->cd();
  mmevttree->Write();
  newfile->Close();

  cout << "\n=================================\n";
  cout << "\n Oniatree skimming is done\n";
  cout << "\n=================================\n";
  t->Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}
