#include "Algorithm.h"
#include "Data.h"

Algorithm::Algorithm() {
  // skip
}

void Algorithm::initialize(Data* data) {
  m_data = data;
  std::cout << "Algorithm::initialize: Attempting to set branches.\n";
  if (m_myTree)
  { // m_myTree가 유효한지 확인
    setBranches();
  }
  else
  {
    std::cerr << "Algorithm::initialize: m_myTree is nullptr. Branches not set.\n";
  }
  // std::cout << "Error here\n";
  // setBranches();
}

void Algorithm::execute(Long64_t irqq) {
  Data& data = *m_data; // local reference

  // set 4-vectors
  data.JP_Reco = static_cast<TLorentzVector*>(data.Reco_QQ_4mom->At(irqq));
  data.mupl_Reco = static_cast<TLorentzVector*>(data.Reco_mu_4mom->At(data.Reco_QQ_mupl_idx[irqq]));
  data.mumi_Reco = static_cast<TLorentzVector*>(data.Reco_mu_4mom->At(data.Reco_QQ_mumi_idx[irqq]));

  // apply selections
  if(!passedSelection(irqq)) return;

  

  // fill the histograms
  standByFilling(irqq);

  // count up nDimu
  // must be increased later than standByFilling()
  ++nDimu;
  // fillHists();
}
bool Algorithm::passedSelection(Long64_t irqq) {
  Data& data = *m_data;

  // apply selection
  if (!(data.Centrality > 0 && data.Centrality < 180)) return false;
  if (!(data.JP_Reco->M() > 2.6 && data.JP_Reco->M() < 3.5)) return false;
  if (!passedSoftMuonCut(irqq)) return false; 
  if (data.Reco_QQ_VtxProb[irqq] < 0.01) return false;
  if (data.Reco_QQ_sign[irqq] != 0) return false;
  if (!passedMuonAcc2018()) return false;

  // passed all selection
  return true;
}

bool Algorithm::passedSoftMuonCut(Long64_t irqq) {
  Data& data = *m_data;
  bool passMuonTypePl = true;
  passMuonTypePl = passMuonTypePl && (data.Reco_mu_SelectionType[data.Reco_QQ_mupl_idx[irqq]] & (1 << 1));
  passMuonTypePl = passMuonTypePl && (data.Reco_mu_SelectionType[data.Reco_QQ_mupl_idx[irqq]] & (1 << 3));

  bool passMuonTypeMi = true;
  passMuonTypeMi = passMuonTypeMi && (data.Reco_mu_SelectionType[data.Reco_QQ_mumi_idx[irqq]] & (1 << 1));
  passMuonTypeMi = passMuonTypeMi && (data.Reco_mu_SelectionType[data.Reco_QQ_mumi_idx[irqq]] & (1 << 3));

  bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
      (data.Reco_mu_nTrkWMea[data.Reco_QQ_mupl_idx[irqq]] > 5) &&
      (data.Reco_mu_nPixWMea[data.Reco_QQ_mupl_idx[irqq]] > 0) &&
      (std::fabs(data.Reco_mu_dxy[data.Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
      (std::fabs(data.Reco_mu_dz[data.Reco_QQ_mupl_idx[irqq]]) < 20.) &&
      passMuonTypePl  //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
  );

  bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
      (data.Reco_mu_nTrkWMea[data.Reco_QQ_mumi_idx[irqq]] > 5) &&
      (data.Reco_mu_nPixWMea[data.Reco_QQ_mumi_idx[irqq]] > 0) &&
      (std::fabs(data.Reco_mu_dxy[data.Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
      (std::fabs(data.Reco_mu_dz[data.Reco_QQ_mumi_idx[irqq]]) < 20.) &&
      passMuonTypeMi //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
  );

  return (muplSoft && mumiSoft);
}

bool Algorithm::passedMuonAcc2018() {
  Data& data = *m_data;

  bool muplAccCut = checkMuonAcc2018(data.mupl_Reco->Pt(), data.mupl_Reco->Eta());
  bool mumiAccCut = checkMuonAcc2018(data.mumi_Reco->Pt(), data.mumi_Reco->Eta());

  return (muplAccCut && mumiAccCut);
}

bool Algorithm::checkMuonAcc2018(double muPt, double muEta) {
  return ( 	( (TMath::Abs(muEta) <= 1.2) && (muPt >=3.5) ) || 
    ( (TMath::Abs(muEta) > 1.2)  && (TMath::Abs(muEta) <= 2.1) && (muPt >= 5.47-1.89*(TMath::Abs(muEta))) ) || 
    ( (TMath::Abs(muEta) > 2.1)  && (TMath::Abs(muEta) <= 2.4) && (muPt >= 1.5) )  
    );
}

void Algorithm::setBranches() {
  Data &data = *m_data;
  m_myTree->Branch("runN", &runN, "runN/I");
  m_myTree->Branch("lumi", &lumi, "lumi/I");


  m_myTree->Branch("event", &evt, "event/I");
  
  m_myTree->Branch("cBin", &cBin, "cBin/I");
  m_myTree->Branch("vz", &vz, "vz/F");

  m_myTree->Branch("nDimu", &nDimu, "nDimu/I");
  m_myTree->Branch("mass", mass, "mass[nDimu]/F");
  m_myTree->Branch("y", y, "y[nDimu]/F");
  m_myTree->Branch("pt", pt, "pt[nDimu]/F");
  m_myTree->Branch("pt1", pt1, "pt1[nDimu]/F");
  m_myTree->Branch("pt2", pt2, "pt2[nDimu]/F");
  m_myTree->Branch("phi", phi, "phi[nDimu]/F");
  m_myTree->Branch("phi1", phi1, "phi1[nDimu]/F");
  m_myTree->Branch("phi2", phi2, "phi2[nDimu]/F");
  m_myTree->Branch("eta", eta, "eta[nDimu]/F");
  m_myTree->Branch("eta1", eta1, "eta1[nDimu]/F");
  m_myTree->Branch("eta2", eta2, "eta2[nDimu]/F");
  m_myTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  m_myTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  m_myTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  m_myTree->Branch("ctau3DTrue", ctau3DTrue, "ctau3DTrue[nDimu]/F");
  m_myTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
  m_myTree->Branch("weight", &weight, "weight/D");
  m_myTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");
}

void Algorithm::standByFilling(Long64_t irqq) {
  Data &data = *m_data;

  evt = data.eventNb;
  runN = data.runNb;
  lumi = data.LS;


  cBin = data.Centrality;
  vz = data.zVtx;

  // Reco variables
  mass[nDimu] = data.JP_Reco->M();
  phi[nDimu] = data.JP_Reco->Phi();
  phi1[nDimu] = data.mupl_Reco->Phi();
  phi2[nDimu] = data.mumi_Reco->Phi();
  y[nDimu] = data.JP_Reco->Rapidity();
  pt[nDimu] = data.JP_Reco->Pt();
  pt1[nDimu] = data.mupl_Reco->Pt();
  pt2[nDimu] = data.mumi_Reco->Pt();
  eta[nDimu] = data.JP_Reco->Eta();
  eta1[nDimu] = data.mupl_Reco->Eta();
  eta2[nDimu] = data.mumi_Reco->Eta();
  // cos_theta[nDimu] = data.JP_Reco->CosTheta();
  // cos_theta1[nDimu] = data.mupl_Reco->CosTheta();
  ctau3D[nDimu] = data.Reco_QQ_ctau3D[irqq];
  ctau3DErr[nDimu] = data.Reco_QQ_ctauErr3D[irqq];
  ctau3DRes[nDimu] = (data.Reco_QQ_ctau3D[irqq]) / (data.Reco_QQ_ctauErr3D[irqq]);

  recoQQsign[nDimu] = data.Reco_QQ_sign[irqq];
  // cos_cs[nDimu] = cos_cs_;
  // phi_cs[nDimu] = phi_cs_;
  // cos_hx[nDimu] = cos_hx_;
  // phi_hx[nDimu] = phi_hx_;
  // cos_px[nDimu] = cos_px_;
  // phi_px[nDimu] = phi_px_;
  // cos_ep[nDimu] = cos_ep_;
  // phi_ep[nDimu] = phi_ep_;
}

void Algorithm::fillTree() {
  // save only when it has dimuons
  if (nDimu>0)
    m_myTree->Fill();
  
  // initialize nDimu for next loop
  nDimu = 0;
}

void Algorithm::fillHists() {
  if (h_runNb_m)
    h_runNb_m->Fill(m_data->runNb);  
}

