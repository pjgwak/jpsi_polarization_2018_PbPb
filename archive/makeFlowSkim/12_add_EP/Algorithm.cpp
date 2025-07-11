#include "Algorithm.h"
#include "Data.h"

Algorithm::Algorithm() {
  // skip
}

void Algorithm::initialize(Data* data, bool isRun2) {
  m_data = data;
  m_isRun2 = isRun2;

  std::cout << "Algorithm::initialize: Attempting to set branches.\n";
  if (m_myTree)
    setBranches();
  else
    std::cerr << "Algorithm::initialize: m_myTree is nullptr. Branches not set.\n";
}

void Algorithm::execute(Long64_t irqq) {
  Data& data = *m_data; // local reference

  Int_t mupl_idx = 0;
  Int_t mumi_idx = 0;

  if (m_isRun2) {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    mupl_idx = dr2->Reco_QQ_mupl_idx[irqq];
    mumi_idx = dr2->Reco_QQ_mumi_idx[irqq];
  } 
  else {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    mupl_idx = static_cast<Int_t>(dr3->Reco_QQ_mupl_idx[irqq]);
    mumi_idx = static_cast<Int_t>(dr3->Reco_QQ_mumi_idx[irqq]);
  }

  // set 4-vectors
  data.JP_Reco = static_cast<TLorentzVector*>(data.Reco_QQ_4mom->At(irqq));
  data.mupl_Reco = static_cast<TLorentzVector *>(data.Reco_mu_4mom->At(mupl_idx));
  data.mumi_Reco = static_cast<TLorentzVector*>(data.Reco_mu_4mom->At(mumi_idx));

  // apply selections
  if(!passSelection(irqq)) return;

  // HX, CS transformation
  MuPlusVector_Helicity(*data.JP_Reco, *data.mupl_Reco); // Call HX before CS
  MuPlusVector_CollinsSoper(*data.JP_Reco);

  // EP tranformation
  if (data.isEP) {
    MuPlusVector_EventPlane(*data.JP_Reco, *data.mupl_Reco);
  }

  // fill the histograms
  standByFilling(irqq);

  // count up nDimu
  // must be increased later than standByFilling()
  ++nDimu;
  // fillHists();
}
bool Algorithm::passSelection(Long64_t irqq) {
  // It is executed in the recoQQ loop

  Data& data = *m_data;

  Int_t qq_sign_val = 0;
  if (m_isRun2) {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    qq_sign_val = dr2->Reco_QQ_sign[irqq];
  }
  else {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    qq_sign_val = static_cast<Int_t>(dr3->Reco_QQ_sign[irqq]);
  }

  // apply selection
  if (cut_centrality0_180 && !(data.Centrality > 0 && data.Centrality < 180)) return false;
  if (cut_jpsiMass && !(data.JP_Reco->M() > 2.6 && data.JP_Reco->M() < 3.5)) return false;
  if (cut_runNb327123 && !passRunNb327123()) return false;
  if (cut_HLTriggerPbPbJpsi2018 && !passHLTriggerJpsiPbPb2018()) return false;
  if (cut_recoQQTrigger && !passRecoQQTrigger(irqq)) return false;
  if (cut_L2L3FilterPbPbJpsi2018 && !passHLFilterJpsiPbPb2018(irqq)) return false;
  if (cut_softMuons && !passSoftMuonCut(irqq)) return false;
  if (cut_vtxProbability && data.Reco_QQ_VtxProb[irqq] < 0.01) return false;
  if (cut_oppositeSign && qq_sign_val != 0) return false;
  if (cut_singleMuonAcc && !passMuonAcc2018()) return false;
  if (cut_whichGen && !passewhichGen(irqq)) return false;
  if (cut_tnp && !passTnpLogic(irqq)) return false;

  // pass all selection
  return true;
}

bool Algorithm::passSoftMuonCut(Long64_t irqq) {
  Data& data = *m_data;

  Int_t mupl_idx = 0;
  Int_t mumi_idx = 0;

  if (m_isRun2) {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    mupl_idx = dr2->Reco_QQ_mupl_idx[irqq];
    mumi_idx = dr2->Reco_QQ_mumi_idx[irqq];
  }
  else {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    mupl_idx = static_cast<Int_t>(dr3->Reco_QQ_mupl_idx[irqq]);
    mumi_idx = static_cast<Int_t>(dr3->Reco_QQ_mumi_idx[irqq]);
  }

  bool passMuonTypePl = true;
  passMuonTypePl = passMuonTypePl && (data.Reco_mu_SelectionType[mupl_idx] & (1 << 1));
  passMuonTypePl = passMuonTypePl && (data.Reco_mu_SelectionType[mupl_idx] & (1 << 3));

  bool passMuonTypeMi = true;
  passMuonTypeMi = passMuonTypeMi && (data.Reco_mu_SelectionType[mumi_idx] & (1 << 1));
  passMuonTypeMi = passMuonTypeMi && (data.Reco_mu_SelectionType[mumi_idx] & (1 << 3));

  bool muplSoft = ( //(Reco_mu_TMOneStaTight[mupl_idx]==true) &&
      (data.Reco_mu_nTrkWMea[mupl_idx] > 5) &&
      (data.Reco_mu_nPixWMea[mupl_idx] > 0) &&
      (std::fabs(data.Reco_mu_dxy[mupl_idx]) < 0.3) &&
      (std::fabs(data.Reco_mu_dz[mupl_idx]) < 20.) &&
      passMuonTypePl  //			 &&  (Reco_mu_highPurity[mupl_idx]==true)
  );

  bool mumiSoft = ( //(Reco_mu_TMOneStaTight[mumi_idx]==true) &&
      (data.Reco_mu_nTrkWMea[mumi_idx] > 5) &&
      (data.Reco_mu_nPixWMea[mumi_idx] > 0) &&
      (std::fabs(data.Reco_mu_dxy[mumi_idx]) < 0.3) &&
      (std::fabs(data.Reco_mu_dz[mumi_idx]) < 20.) &&
      passMuonTypeMi //			 &&  (Reco_mu_highPurity[mumi_idx]==true)
  );

  return (muplSoft && mumiSoft);
}

bool Algorithm::passMuonAcc2018() {
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

bool Algorithm::passRunNb327123() {
  // this is event level cut
  // it will be executed in the Reco loop
  // to manange all cuts in the Algorithm()
  Data& data = *m_data;
  
  if (data.runNb < 327123)
    return true;
  else
    return false;
}

bool Algorithm::passHLFilterJpsiPbPb2018(Long64_t irqq) {
  // this is event level cut
  // it will be executed in the Reco loop
  // to manange all cuts in the Algorithm()
  Data& data = *m_data;

  Int_t mupl_idx = 0;
  Int_t mumi_idx = 0;

  if (m_isRun2) {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    mupl_idx = dr2->Reco_QQ_mupl_idx[irqq];
    mumi_idx = dr2->Reco_QQ_mumi_idx[irqq];
  } 
  else {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    mupl_idx = static_cast<Int_t>(dr3->Reco_QQ_mupl_idx[irqq]);
    mumi_idx = static_cast<Int_t>(dr3->Reco_QQ_mumi_idx[irqq]);
  }

  return ((data.Reco_mu_trig[mupl_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018)) && (data.Reco_mu_trig[mumi_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018)));
}

bool Algorithm::passHLTriggerJpsiPbPb2018() {
  // this is event level cut
  // it will be executed in the Reco loop
  // to manange all cuts in the Algorithm()
  Data& data = *m_data;

  return ((data.HLTriggers & ((ULong64_t)pow(2, kTrigJpsi2018))) == ((ULong64_t)pow(2, kTrigJpsi2018)));
}

bool Algorithm::passRecoQQTrigger(Long64_t irqq) {
  Data& data = *m_data;

  return ((data.Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigJpsi2018))) == ((ULong64_t)pow(2, kTrigJpsi2018)));
}

bool Algorithm::passewhichGen(Long64_t irqq) {
  // Data &data = *m_data;

  Int_t mupl_idx = 0;
  Int_t mumi_idx = 0;

  if (m_isRun2)
  {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    mupl_idx = dr2->Reco_QQ_mupl_idx[irqq];
    mumi_idx = dr2->Reco_QQ_mumi_idx[irqq];

    return (dr2->Reco_mu_whichGen[mupl_idx] != -1 && dr2->Reco_mu_whichGen[mumi_idx] != -1);
  }
  else
  {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    mupl_idx = static_cast<Int_t>(dr3->Reco_QQ_mupl_idx[irqq]);
    mumi_idx = static_cast<Int_t>(dr3->Reco_QQ_mumi_idx[irqq]);

    return (dr3->Reco_mu_whichGen[mupl_idx] != -1 && dr3->Reco_mu_whichGen[mumi_idx] != -1);
  }
}

bool Algorithm::passTnpLogic(Long64_t irqq) {
  // I keep naming covention of original codes
  // Maybe we can conver it later...

  Data &data = *m_data;

  Int_t mupl_idx = 0;
  Int_t mumi_idx = 0;

  if (m_isRun2)
  {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    mupl_idx = dr2->Reco_QQ_mupl_idx[irqq];
    mumi_idx = dr2->Reco_QQ_mumi_idx[irqq];
  }
  else
  {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    mupl_idx = static_cast<Int_t>(dr3->Reco_QQ_mupl_idx[irqq]);
    mumi_idx = static_cast<Int_t>(dr3->Reco_QQ_mumi_idx[irqq]);
  }

  m_tnpWeight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  m_tnpWeight *= tnp_weight_muid_pbpb(data.mupl_Reco->Pt(), data.mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(data.mumi_Reco->Pt(), data.mumi_Reco->Eta(), 0); // mu id
  m_tnpWeight *= tnp_weight_trk_pbpb(data.mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(data.mumi_Reco->Eta(), 0);                                     // inner tracker

  // Trigger part
  if (!((data.Reco_mu_trig[mupl_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018)) && (data.Reco_mu_trig[mumi_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018))))
  {
    //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
    //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
    return false;
  }
  bool mupl_L2Filter = ((data.Reco_mu_trig[mupl_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018))) ? true : false;
  bool mupl_L3Filter = ((data.Reco_mu_trig[mupl_idx] & ((ULong64_t)pow(2, kL3filterJpsi2018))) == ((ULong64_t)pow(2, kL3filterJpsi2018))) ? true : false;
  bool mumi_L2Filter = ((data.Reco_mu_trig[mumi_idx] & ((ULong64_t)pow(2, kL2filterJpsi2018))) == ((ULong64_t)pow(2, kL2filterJpsi2018))) ? true : false;
  bool mumi_L3Filter = ((data.Reco_mu_trig[mumi_idx] & ((ULong64_t)pow(2, kL3filterJpsi2018))) == ((ULong64_t)pow(2, kL3filterJpsi2018))) ? true : false;
  if (mupl_L2Filter == false || mumi_L2Filter == false)
  {
    std::cout << "TnP ERROR !!!! ::: No matched L2 filter2\n\n";
    return false;
  }

  bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
  bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
  bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
  bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
  bool SelDone = false;

  if (mupl_isL2 && mumi_isL3)
  {
    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(data.mupl_Reco->Pt(), data.mupl_Reco->Eta(), 2, 0);
    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(data.mumi_Reco->Pt(), data.mumi_Reco->Eta(), 3, 0);
    SelDone = true;
  }
  else if (mupl_isL3 && mumi_isL2)
  {
    tnp_trig_weight_mupl = tnp_weight_trg_pbpb(data.mupl_Reco->Pt(), data.mupl_Reco->Eta(), 3, 0);
    tnp_trig_weight_mumi = tnp_weight_trg_pbpb(data.mumi_Reco->Pt(), data.mumi_Reco->Eta(), 2, 0);
    SelDone = true;
  }
  else if (mupl_isL3 && mumi_isL3)
  {
    int t[2] = {-1, 1}; // mupl, mumi
    int l = rand() % (2);
    // pick up what will be L2
    if (t[l] == -1)
    {
      tnp_trig_weight_mupl = tnp_weight_trg_pbpb(data.mupl_Reco->Pt(), data.mupl_Reco->Eta(), 2, 0);
      tnp_trig_weight_mumi = tnp_weight_trg_pbpb(data.mumi_Reco->Pt(), data.mumi_Reco->Eta(), 3, 0);
    }
    else if (t[l] == 1)
    {
      tnp_trig_weight_mupl = tnp_weight_trg_pbpb(data.mupl_Reco->Pt(), data.mupl_Reco->Eta(), 3, 0);
      tnp_trig_weight_mumi = tnp_weight_trg_pbpb(data.mumi_Reco->Pt(), data.mumi_Reco->Eta(), 2, 0);
    }
    else
    {
      std::cout << "ERROR :: No random selection done !!!!\n\n";
      return false;
    }
    SelDone = true;
  }
  if (SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1))
  {
    std::cout << "ERROR :: No muon filter combination selected !!!!\n\n";
    return false;
  }

  m_tnpWeight *= tnp_trig_weight_mupl * tnp_trig_weight_mumi;
  return true;
}

void Algorithm::setBranches() {
  Data &data = *m_data;
  if (!data.isMC) {
    m_myTree->Branch("runN", &runN, "runN/I");
    m_myTree->Branch("lumi", &lumi, "lumi/I");
  }

  m_myTree->Branch("event", &evt, "event/I");
  m_myTree->Branch("cBin", &cBin, "cBin/I");
  m_myTree->Branch("vz", &vz, "vz/F");
  m_myTree->Branch("nDimu", &nDimu, "nDimu/I");
  m_myTree->Branch("weight", &weight, "weight/D");

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
  m_myTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");
  m_myTree->Branch("cosHX", cosHX, "cosHX[nDimu]/F");
  m_myTree->Branch("phiHX", phiHX, "phiHX[nDimu]/F");
  m_myTree->Branch("cosCS", cosCS, "cosCS[nDimu]/F");
  m_myTree->Branch("phiCS", phiCS, "phiCS[nDimu]/F");

  if (data.isEP) {
    m_myTree->Branch("cosEP", cosEP, "cosEP[nDimu]/F");
    m_myTree->Branch("phiEP", phiEP, "phiEP[nDimu]/F");
  }
}

void Algorithm::standByFilling(Long64_t irqq) {
  Data &data = *m_data;

  Int_t current_qq_sign_val = 0;
  if (m_isRun2) {
    DataRun2 *dr2 = dynamic_cast<DataRun2 *>(m_data);
    current_qq_sign_val = dr2->Reco_QQ_sign[irqq];
  }
  else {
    DataRun3 *dr3 = dynamic_cast<DataRun3 *>(m_data);
    current_qq_sign_val = static_cast<Int_t>(dr3->Reco_QQ_sign[irqq]);
  }

  if (!data.isMC) {
    runN = data.runNb;
    lumi = data.LS;
  }

  evt = data.eventNb;
  cBin = getHiBinFromhiHF(data.SumET_HF);
  vz = data.zVtx;

  if (data.isMC) {
    weight = findNcoll(data.Centrality) * data.Gen_weight;
    TnPweight[nDimu] = m_tnpWeight;
  }

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
  recoQQsign[nDimu] = current_qq_sign_val;
  cosHX[nDimu] = muplHX.CosTheta();
  phiHX[nDimu] = muplHX.Phi();
  cosCS[nDimu] = muplCS.CosTheta();
  phiCS[nDimu] = muplCS.Phi();

  if (data.isEP) {
    cosEP[nDimu] = m_cosEP;
    phiEP[nDimu] = m_phiEP;
  }
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

Double_t Algorithm::findNcoll(int hiBin) {
  // 2018 case.
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
}

Int_t Algorithm::getHiBinFromhiHF(const Double_t hiHF) {
  Int_t binPos = -1;
  for (int i = 0; i < nBins; ++i)
  {
    if (hiHF >= binTable[i] && hiHF < binTable[i + 1])
    {
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200 * ((Double_t)binPos) / ((Double_t)nBins));
}

void Algorithm::MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
  // ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
  TVector3 QQVector_Lab = QQLV_Lab.Vect();
  TLorentzVector MuPlusLV_QQRestFrame(MuPlusLV_Lab);

  //(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
  //(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
  MuPlusLV_QQRestFrame.Boost(-QQLV_Lab.BoostVector());

  // ******** Rotate the coordinates ******** //
  muplHX = MuPlusLV_QQRestFrame.Vect();

  // Note: TVector3.Rotate() rotates the vectors, not the coordinates, so should rotate -phi and -theta
  muplHX.RotateZ(-QQVector_Lab.Phi());
  muplHX.RotateY(-QQVector_Lab.Theta());
}

void Algorithm::MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab) {
  // one must call a MuPlusVector_Helicity() before call it.

  // ******** Set beam energy for the Collins-Soper reference frame ******** //
  double sqrt_S_NN = 5.02;                    //(Center of mass Energy per nucleon pair in TeV)
  double beamEnergy = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)

  // ******** HX to CS (rotation from HX frame to CS frame) ******** //
  // (1. Boost two beams to upsilon's rest frame)
  // (2. Rotate the coordinates)
  // (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
  // (4. Calculate delta (angle btw ZHX and ZCS))

  // ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
  TLorentzVector Beam1LV_Boosted(0., 0., beamEnergy, beamEnergy);
  TLorentzVector Beam2LV_Boosted(0., 0., -beamEnergy, beamEnergy); // mind the sign!!

  Beam1LV_Boosted.Boost(-QQLV_Lab.BoostVector());
  Beam2LV_Boosted.Boost(-QQLV_Lab.BoostVector());

  // ******** Rotate the coordinates ******** //
  TVector3 Beam1Vector_QQRestFrame(Beam1LV_Boosted.Vect());
  TVector3 Beam2Vector_QQRestFrame(Beam2LV_Boosted.Vect());

  TVector3 QQVector_Lab = QQLV_Lab.Vect();
  Beam1Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
  Beam1Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

  Beam2Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
  Beam2Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

  // ******** Calculate the angle between z_HX and z_CS ******** //
  TVector3 ZHXunitVec(0, 0, 1.);                                                 //(define z_HX unit vector)
  double Angle_B1ZHX = Beam1Vector_QQRestFrame.Angle(ZHXunitVec);                //(angle between beam1 and z_HX)
  double Angle_B2ZHX = Beam2Vector_QQRestFrame.Angle(-ZHXunitVec);               //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
  double Angle_B1miB2 = Beam1Vector_QQRestFrame.Angle(-Beam2Vector_QQRestFrame); //(angle between beam1 and -beam2)

  double delta = 0; //(define and initialize the angle between z_HX and z_CS)

  // Maths for calculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
  if (Angle_B1ZHX > Angle_B2ZHX)
    delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
  else if (Angle_B1ZHX < Angle_B2ZHX)
    delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
  else
    std::cout << "beam1PvecBoosted.Pz() = 0?" << std::endl;

  // ******** Rotate the coordinates along the y-axis by the angle between z_HX and z_CS ******** //
  muplCS = TVector3(muplHX);

  muplCS.RotateY(delta);
}

void Algorithm::MuPlusVector_EventPlane(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
  Data &data = *m_data;

  // ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
  TVector3 QQVector_Lab = QQLV_Lab.Vect();
  TLorentzVector MuPlusLV_QQRestFrame(MuPlusLV_Lab);

  //(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
  //(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
  MuPlusLV_QQRestFrame.Boost(-QQLV_Lab.BoostVector());

  // ******** Rotate the coordinates ******** //
  muplEP = MuPlusLV_QQRestFrame.Vect();
  muplEP.Unit();

  // set z-axis of EP
  TVector3 uz_lab(0, 0, 1); // beam direction in lab frame
  TVector3 uy_lab(0, 1, 0); // unit vector of y direction in lab frame
  TVector3 uz_ep = uy_lab;
  // =========================
  // How can we consider auto-correlations?
  // =========================
  uz_ep.Rotate(data.ep_flat, uz_lab); // rotate z axis as EP angles

  // calculate cos and phi in EP frame
  TVector3 uy_ep = uz_lab;

  TVector3 ux_ep = uy_ep.Cross(uz_ep);
  ux_ep = ux_ep.Unit();

  m_cosEP = muplEP.Dot(uz_ep);
  m_phiEP = TMath::ATan2(muplEP.Dot(uy_ep), muplEP.Dot(ux_ep));
}

// ======================= //
// ===== TnP tables ===== //
// ======================= //
double Algorithm::tnp_weight_muid_pbpb(double pt, double eta, int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  // SF for 0 < |eta| < 1.2
  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    if (x >= 3.5 && x <4) {num = 0.946125; den = 0.95062; statUp = 0.00377662; statDown = 0.00386775;}
    else if (x >= 4 && x <4.5) {num = 0.973366; den = 0.975256; statUp = 0.00284959; statDown = 0.00295956;}
    else if (x >= 4.5 && x <5) {num = 0.984622; den = 0.983749; statUp = 0.00263774; statDown = 0.00277141;}
    else if (x >= 5 && x <5.5) {num = 0.985061; den = 0.98948; statUp = 0.00285151; statDown = 0.00301546;}
    else if (x >= 5.5 && x <6.5) {num = 0.988937; den = 0.99051; statUp = 0.0020037; statDown = 0.00211913;}
    else if (x >= 6.5 && x <8) {num = 0.991257; den = 0.992598; statUp = 0.0017347; statDown = 0.00185756;}
    else if (x >= 8 && x <10.5) {num = 0.986505; den = 0.9885; statUp = 0.00213856; statDown = 0.00228346;}
    else if (x >= 10.5 && x <14) {num = 0.975903; den = 0.967524; statUp = 0.00364064; statDown = 0.0038785;}
    else if (x >= 14 && x <18) {num = 0.979888; den = 0.978318; statUp = 0.00488479; statDown = 0.00538685;}
    else {num = 0.983224; den = 0.984461; statUp = 0.0056352; statDown = 0.00635297;}
  }
  // SF for 1.2 < |eta| < 1.8
  if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    if (x >= 2.37 && x <3) {num = 0.983353; den = 0.981734; statUp = 0.00668073; statDown = 0.00681507;}
    else if (x >= 3 && x <3.5) {num = 0.985622; den = 0.984888; statUp = 0.00411601; statDown = 0.00425201;}
    else if (x >= 3.5 && x <4) {num = 0.980436; den = 0.989513; statUp = 0.00390367; statDown = 0.00408228;}
    else if (x >= 4 && x <4.5) {num = 0.986109; den = 0.991777; statUp = 0.00358665; statDown = 0.00382294;}
    else if (x >= 4.5 && x <5) {num = 0.992053; den = 0.993039; statUp = 0.00353145; statDown = 0.00382873;}
    else if (x >= 5 && x <6) {num = 0.989106; den = 0.994324; statUp = 0.00288884; statDown = 0.00311446;}
    else if (x >= 6 && x <7.5) {num = 0.995133; den = 0.997452; statUp = 0.00241041; statDown = 0.00267219;}
    else if (x >= 7.5 && x <10) {num = 0.998692; den = 0.996495; statUp = 0.00130754; statDown = 0.00252507;}
    else if (x >= 10 && x <15) {num = 0.979734; den = 0.971008; statUp = 0.0202658; statDown = 0.00568435;}
    else {num = 0.996433; den = 0.987155; statUp = 0.0035675; statDown = 0.00767196;}
  }
  // SF for 1.8 < |eta| < 2.1
  if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    if (x >= 1.8 && x <2.5) {num = 0.970771; den = 0.955555; statUp = 0.010087; statDown = 0.0102019;}
    else if (x >= 2.5 && x <3) {num = 0.984571; den = 0.983608; statUp = 0.00598521; statDown = 0.00625035;}
    else if (x >= 3 && x <3.5) {num = 0.991532; den = 0.99085; statUp = 0.00485266; statDown = 0;}
    else if (x >= 3.5 && x <4) {num = 0.982929; den = 0.995395; statUp = 0.00477394; statDown = 0.00523748;}
    else if (x >= 4 && x <4.5) {num = 0.997191; den = 0.998022; statUp = 0.00280897; statDown = 0.00374539;}
    else if (x >= 4.5 && x <5.5) {num = 0.996918; den = 0.997783; statUp = 0.00226149; statDown = 0.00267054;}
    else if (x >= 5.5 && x <7) {num = 1; den = 0.999082; statUp = 9.82446e-09; statDown = 0.00120742;}
    else if (x >= 7 && x <9) {num = 0.996645; den = 0.999904; statUp = 0.00294462; statDown = 0.00358784;}
    else if (x >= 9 && x <12) {num = 0.997479; den = 0.997288; statUp = 0.00252069; statDown = 0.00478073;}
    else {num = 0.989358; den = 0.997918; statUp = 0.00514686; statDown = 0.00685444;}
  }
  // SF for 2.1 < |eta| < 2.4
  if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    if (x >= 1.8 && x <2.2) {num = 0.89788; den = 0.929933; statUp = 0.0167719; statDown = 0.0168054;}
    else if (x >= 2.2 && x <2.7) {num = 0.940431; den = 0.968575; statUp = 0.0107689; statDown = 0.0109474;}
    else if (x >= 2.7 && x <3.2) {num = 0.972522; den = 0.980868; statUp = 0.0084141; statDown = 0.00875246;}
    else if (x >= 3.2 && x <3.7) {num = 0.977761; den = 0.988464; statUp = 0.00677253; statDown = 0.00722554;}
    else if (x >= 3.7 && x <4.7) {num = 0.997206; den = 0.994088; statUp = 0.00279445; statDown = 0.00494402;}
    else if (x >= 4.7 && x <8) {num = 0.996334; den = 0.99794; statUp = 0.0025871; statDown = 0.00292285;}
    else if (x >= 8 && x <11) {num = 1; den = 0.999768; statUp = 6.83212e-09; statDown = 0.00252398;}
    else if (x >= 11 && x <14) {num = 0.978986; den = 0.998297; statUp = 0.0099437; statDown = 0.0135461;}
    else {num = 1; den = 0.996328; statUp = 1.46346e-07; statDown = 0.00540197;}
  }

  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    // syst uncertainties
    if (x >= 3.5 && x < 4) syst = 0.000855038;
    if (x >= 4 && x < 4.5) syst = 0.000240768;
    if (x >= 4.5 && x < 5) syst = 0.000899317;
    if (x >= 5 && x < 5.5) syst = 0.00109601;
    if (x >= 5.5 && x < 6.5) syst = 0.000311285;
    if (x >= 6.5 && x < 8) syst = 0.000476892;
    if (x >= 8 && x < 10.5) syst = 0.00122726;
    if (x >= 10.5 && x < 14) syst = 0.00274564;
    if (x >= 14 && x < 18) syst = 0.00133458;
    if (x >= 18 && x < 30) syst = 0.00365302;
  }
  if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    // syst uncertainties
    if (x >= 2.37 && x < 3) syst = 0.00121011;
    if (x >= 3 && x < 3.5) syst = 0.00177203;
    if (x >= 3.5 && x < 4) syst = 0.00181968;
    if (x >= 4 && x < 4.5) syst = 0.00154391;
    if (x >= 4.5 && x < 5) syst = 0.00183399;
    if (x >= 5 && x < 6) syst = 0.000750531;
    if (x >= 6 && x < 7.5) syst = 0.00111009;
    if (x >= 7.5 && x < 10) syst = 0.00123526;
    if (x >= 10 && x < 15) syst = 0.00156111;
    if (x >= 15 && x < 30) syst = 0.00311647;
  }
  if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    // syst uncertainties
    if (x >= 1.8 && x < 2.5) syst = 0.00425782;
    if (x >= 2.5 && x < 3) syst = 0.00100318;
    if (x >= 3 && x < 3.5) syst = 0.00476823;
    if (x >= 3.5 && x < 4) syst = 0.00107436;
    if (x >= 4 && x < 4.5) syst = 0.00280821;
    if (x >= 4.5 && x < 5.5) syst = 0.000743448;
    if (x >= 5.5 && x < 7) syst = 9.14067e-07;
    if (x >= 7 && x < 9) syst = 0.00197499;
    if (x >= 9 && x < 12) syst = 0.000842135;
    if (x >= 12 && x < 20) syst = 0.000789339;
  }
  if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    // syst uncertainties
    if (x >= 1.8 && x < 2.2) syst = 0.0234659;
    if (x >= 2.2 && x < 2.7) syst = 0.00575045;
    if (x >= 2.7 && x < 3.2) syst = 0.00231199;
    if (x >= 3.2 && x < 3.7) syst = 0.00333113;
    if (x >= 3.7 && x < 4.7) syst = 0.00183657;
    if (x >= 4.7 && x < 8) syst = 0.00226061;
    if (x >= 8 && x < 11) syst = 2.45089e-07;
    if (x >= 11 && x < 14) syst = 0.00132322;
    if (x >= 14 && x < 20) syst = 1.46346e-07;
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//              T R G     P b P b                //
///////////////////////////////////////////////////

double Algorithm::tnp_weight_trg_pbpb(double pt, double eta, int filterId,int idx, bool getRawDen) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (filterId==0) { //L2 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.663241; den = 0.612464; statUp = 0.00756585; statDown = 0.00760942;}
      else if (x >= 4 && x <4.5) {num = 0.85917; den = 0.842014; statUp = 0.00553614; statDown = 0.0056399;}
      else if (x >= 4.5 && x <5) {num = 0.916521; den = 0.900241; statUp = 0.00469998; statDown = 0.00484238;}
      else if (x >= 5 && x <5.5) {num = 0.934641; den = 0.919918; statUp = 0.0045738; statDown = 0.00475841;}
      else if (x >= 5.5 && x <6.5) {num = 0.947297; den = 0.940729; statUp = 0.00343701; statDown = 0.00356287;}
      else if (x >= 6.5 && x <8) {num = 0.949896; den = 0.95601; statUp = 0.00340062; statDown = 0.00353618;}
      else if (x >= 8 && x <10.5) {num = 0.950506; den = 0.96243; statUp = 0.00363768; statDown = 0.00379628;}
      else if (x >= 10.5 && x <14) {num = 0.947321; den = 0.964831; statUp = 0.00478944; statDown = 0.00505914;}
      else if (x >= 14 && x <18) {num = 0.940314; den = 0.966093; statUp = 0.00760886; statDown = 0.00820262;}
      else {num = 0.943817; den = 0.957341; statUp = 0.00874496; statDown = 0.00954934;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.698036; den = 0.656794; statUp = 0.0116241; statDown = 0.0117424;}
      else if (x >= 3 && x <3.5) {num = 0.829209; den = 0.793366; statUp = 0.00781434; statDown = 0.00789364;}
      else if (x >= 3.5 && x <4) {num = 0.900878; den = 0.867213; statUp = 0.00646057; statDown = 0;}
      else if (x >= 4 && x <4.5) {num = 0.917489; den = 0.917832; statUp = 0.00630865; statDown = 0.00655031;}
      else if (x >= 4.5 && x <5) {num = 0.935786; den = 0.936854; statUp = 0.00635004; statDown = 0.00668515;}
      else if (x >= 5 && x <6) {num = 0.940635; den = 0.950382; statUp = 0.00519892; statDown = 0.00546059;}
      else if (x >= 6 && x <7.5) {num = 0.942288; den = 0.955739; statUp = 0.00542; statDown = 0.00571559;}
      else if (x >= 7.5 && x <10) {num = 0.953266; den = 0.959009; statUp = 0.00550351; statDown = 0.00587165;}
      else if (x >= 10 && x <15) {num = 0.948259; den = 0.959592; statUp = 0.00686182; statDown = 0.0074032;}
      else {num = 0.913688; den = 0.954037; statUp = 0.0126313; statDown = 0.0139534;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.641279; den = 0.556324; statUp = 0.031358; statDown = 0.0311058;}
      else if (x >= 2 && x <2.5) {num = 0.746463; den = 0.740273; statUp = 0.014056; statDown = 0.0141737;}
      else if (x >= 2.5 && x <3) {num = 0.882257; den = 0.893319; statUp = 0.0101587; statDown = 0.0104338;}
      else if (x >= 3 && x <3.5) {num = 0.908732; den = 0.927004; statUp = 0.00969244; statDown = 0.0101057;}
      else if (x >= 3.5 && x <4) {num = 0.918071; den = 0.943802; statUp = 0.00962387; statDown = 0.010125;}
      else if (x >= 4 && x <4.5) {num = 0.907847; den = 0.946961; statUp = 0.0108044; statDown = 0.0115161;}
      else if (x >= 4.5 && x <5.5) {num = 0.930721; den = 0.954231; statUp = 0.00791393; statDown = 0.00835311;}
      else if (x >= 5.5 && x <6.5) {num = 0.892041; den = 0.952605; statUp = 0.0117125; statDown = 0.0123847;}
      else if (x >= 6.5 && x <8) {num = 0.917612; den = 0.946287; statUp = 0.0198739; statDown = 0;}
      else if (x >= 8 && x <9.5) {num = 0.914016; den = 0.940696; statUp = 0.0135905; statDown = 0.0147506;}
      else if (x >= 9.5 && x <13) {num = 0.915476; den = 0.936928; statUp = 0.0139215; statDown = 0.0151525;}
      else {num = 0.905898; den = 0.924578; statUp = 0.0221738; statDown = 0.0243314;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.789011; den = 0.794511; statUp = 0.0182577; statDown = 0.018358;}
      else if (x >= 2.2 && x <2.7) {num = 0.854198; den = 0.867878; statUp = 0.0140189; statDown = 0.0142627;}
      else if (x >= 2.7 && x <3.2) {num = 0.878207; den = 0.893527; statUp = 0.0129086; statDown = 0.0132743;}
      else if (x >= 3.2 && x <3.7) {num = 0.877411; den = 0.911546; statUp = 0.0121215; statDown = 0.0125934;}
      else if (x >= 3.7 && x <4.7) {num = 0.885503; den = 0.912297; statUp = 0.0108807; statDown = 0.011314;}
      else if (x >= 4.7 && x <6.5) {num = 0.912669; den = 0.928676; statUp = 0.00953789; statDown = 0.010012;}
      else if (x >= 6.5 && x <8.5) {num = 0.900634; den = 0.922926; statUp = 0.014265; statDown = 0.0151886;}
      else if (x >= 8.5 && x <11) {num = 0.869202; den = 0.916807; statUp = 0.0198993; statDown = 0.021416;}
      else {num = 0.891181; den = 0.915654; statUp = 0.0229342; statDown = 0.0252136;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00136539;
      if (x >= 4 && x < 4.5) syst = 0.00105145;
      if (x >= 4.5 && x < 5) syst = 0.00132265;
      if (x >= 5 && x < 5.5) syst = 0.000788367;
      if (x >= 5.5 && x < 6.5) syst = 0.000562329;
      if (x >= 6.5 && x < 8) syst = 0.000155103;
      if (x >= 8 && x < 10.5) syst = 0.000158629;
      if (x >= 10.5 && x < 14) syst = 0.000280705;
      if (x >= 14 && x < 18) syst = 0.00108834;
      if (x >= 18 && x < 30) syst = 0.00305146;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00329732;
      if (x >= 3 && x < 3.5) syst = 0.00520027;
      if (x >= 3.5 && x < 4) syst = 0.00240912;
      if (x >= 4 && x < 4.5) syst = 0.00222157;
      if (x >= 4.5 && x < 5) syst = 0.0017746;
      if (x >= 5 && x < 6) syst = 0.00167846;
      if (x >= 6 && x < 7.5) syst = 0.000384188;
      if (x >= 7.5 && x < 10) syst = 0.00184184;
      if (x >= 10 && x < 15) syst = 0.000586472;
      if (x >= 15 && x < 30) syst = 0.00495382;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0133628;
      if (x >= 2 && x < 2.5) syst = 0.0347588;
      if (x >= 2.5 && x < 3) syst = 0.00227196;
      if (x >= 3 && x < 3.5) syst = 0.00511758;
      if (x >= 3.5 && x < 4) syst = 0.00466125;
      if (x >= 4 && x < 4.5) syst = 0.0105296;
      if (x >= 4.5 && x < 5.5) syst = 0.00242737;
      if (x >= 5.5 && x < 6.5) syst = 0.00472684;
      if (x >= 6.5 && x < 8) syst = 0.00762573;
      if (x >= 8 && x < 9.5) syst = 0.00627021;
      if (x >= 9.5 && x < 13) syst = 0.00231563;
      if (x >= 13 && x < 20) syst = 0.019257;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0690381;
      if (x >= 2.2 && x < 2.7) syst = 0.0073301;
      if (x >= 2.7 && x < 3.2) syst = 0.00969524;
      if (x >= 3.2 && x < 3.7) syst = 0.00401152;
      if (x >= 3.7 && x < 4.7) syst = 0.0131915;
      if (x >= 4.7 && x < 6.5) syst = 0.00669816;
      if (x >= 6.5 && x < 8.5) syst = 0.00311553;
      if (x >= 8.5 && x < 11) syst = 0.00999295;
      if (x >= 11 && x < 20) syst = 0.0151298;
    }
  }
  if (filterId==1) { //L3 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0981413; den = 0.071216; statUp = 0.00475984; statDown = 0.00464341;}
      else if (x >= 4 && x <4.5) {num = 0.309366; den = 0.234954; statUp = 0.00730881; statDown = 0.00724846;}
      else if (x >= 4.5 && x <5) {num = 0.49696; den = 0.427459; statUp = 0.00850388; statDown = 0.00850324;}
      else if (x >= 5 && x <5.5) {num = 0.646567; den = 0.569803; statUp = 0.00897182; statDown = 0.00902901;}
      else if (x >= 5.5 && x <6.5) {num = 0.717727; den = 0.665688; statUp = 0.00690496; statDown = 0.00696443;}
      else if (x >= 6.5 && x <8) {num = 0.770814; den = 0.736851; statUp = 0.00662338; statDown = 0.00670205;}
      else if (x >= 8 && x <10.5) {num = 0.791509; den = 0.777511; statUp = 0.0068552; statDown = 0.00695591;}
      else if (x >= 10.5 && x <14) {num = 0.826768; den = 0.814073; statUp = 0.00832433; statDown = 0.00852041;}
      else if (x >= 14 && x <18) {num = 0.799276; den = 0.820805; statUp = 0.0131407; statDown = 0.0135476;}
      else {num = 0.844468; den = 0.83657; statUp = 0.0140064; statDown = 0.0146299;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.331992; den = 0.302908; statUp = 0.011463; statDown = 0.0112805;}
      else if (x >= 3 && x <3.5) {num = 0.429327; den = 0.424718; statUp = 0.00993062; statDown = 0.00987992;}
      else if (x >= 3.5 && x <4) {num = 0.543359; den = 0.527324; statUp = 0.00926375; statDown = 0.00951398;}
      else if (x >= 4 && x <4.5) {num = 0.590264; den = 0.603494; statUp = 0.0113814; statDown = 0.0114274;}
      else if (x >= 4.5 && x <5) {num = 0.638602; den = 0.645253; statUp = 0.0126665; statDown = 0.0127591;}
      else if (x >= 5 && x <6) {num = 0.671782; den = 0.678571; statUp = 0.0104761; statDown = 0.0105625;}
      else if (x >= 6 && x <7.5) {num = 0.683073; den = 0.719277; statUp = 0.0110239; statDown = 0.0111476;}
      else if (x >= 7.5 && x <10) {num = 0.746848; den = 0.753451; statUp = 0.0114305; statDown = 0.0116301;}
      else if (x >= 10 && x <15) {num = 0.772519; den = 0.800301; statUp = 0.0133498; statDown = 0.013682;}
      else {num = 0.768513; den = 0.831634; statUp = 0.0213549; statDown = 0.0221413;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.0302057; den = 0.0291672; statUp = 0.0103319; statDown = 0.00934636;}
      else if (x >= 2 && x <2.5) {num = 0.173305; den = 0.159501; statUp = 0.0110778; statDown = 0.0107941;}
      else if (x >= 2.5 && x <3) {num = 0.470274; den = 0.464604; statUp = 0.0151192; statDown = 0.0150115;}
      else if (x >= 3 && x <3.5) {num = 0.661007; den = 0.654318; statUp = 0.0158638; statDown = 0.0159609;}
      else if (x >= 3.5 && x <4) {num = 0.707769; den = 0.711729; statUp = 0.0158855; statDown = 0.0161;}
      else if (x >= 4 && x <4.5) {num = 0.726282; den = 0.741631; statUp = 0.0174238; statDown = 0.017772;}
      else if (x >= 4.5 && x <5.5) {num = 0.763392; den = 0.796057; statUp = 0.0131827; statDown = 0.013438;}
      else if (x >= 5.5 && x <6.5) {num = 0.767094; den = 0.82346; statUp = 0.016378; statDown = 0.0168236;}
      else if (x >= 6.5 && x <8) {num = 0.810425; den = 0.833902; statUp = 0.0145812; statDown = 0.0150811;}
      else if (x >= 8 && x <9.5) {num = 0.819571; den = 0.840655; statUp = 0.0190366; statDown = 0.01987;}
      else if (x >= 9.5 && x <13) {num = 0.829677; den = 0.856559; statUp = 0.0192692; statDown = 0.0202277;}
      else {num = 0.862606; den = 0.873057; statUp = 0.0264458; statDown = 0.0282459;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.0980329; den = 0.092124; statUp = 0.0109739; statDown = 0.0105493;}
      else if (x >= 2.2 && x <2.7) {num = 0.264092; den = 0.248237; statUp = 0.0154033; statDown = 0.0150975;}
      else if (x >= 2.7 && x <3.2) {num = 0.417324; den = 0.407684; statUp = 0.0182234; statDown = 0.0179965;}
      else if (x >= 3.2 && x <3.7) {num = 0.486753; den = 0.49461; statUp = 0.0184525; statDown = 0.0183476;}
      else if (x >= 3.7 && x <4.7) {num = 0.569753; den = 0.574351; statUp = 0.0169356; statDown = 0.0169599;}
      else if (x >= 4.7 && x <6.5) {num = 0.674727; den = 0.651274; statUp = 0.0161099; statDown = 0.0162599;}
      else if (x >= 6.5 && x <8.5) {num = 0.720956; den = 0.710909; statUp = 0.0216217; statDown = 0.0220983;}
      else if (x >= 8.5 && x <11) {num = 0.714358; den = 0.74928; statUp = 0.0275205; statDown = 0.0283679;}
      else {num = 0.748962; den = 0.778156; statUp = 0.0324675; statDown = 0.0339396;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000650706;
      if (x >= 4 && x < 4.5) syst = 0.00108325;
      if (x >= 4.5 && x < 5) syst = 0.00298052;
      if (x >= 5 && x < 5.5) syst = 0.00341277;
      if (x >= 5.5 && x < 6.5) syst = 0.000613358;
      if (x >= 6.5 && x < 8) syst = 0.000661917;
      if (x >= 8 && x < 10.5) syst = 0.000756817;
      if (x >= 10.5 && x < 14) syst = 0.000783331;
      if (x >= 14 && x < 18) syst = 0.00218267;
      if (x >= 18 && x < 30) syst = 0.00267133;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00235318;
      if (x >= 3 && x < 3.5) syst = 0.0041557;
      if (x >= 3.5 && x < 4) syst = 0.00434833;
      if (x >= 4 && x < 4.5) syst = 0.00151395;
      if (x >= 4.5 && x < 5) syst = 0.00325726;
      if (x >= 5 && x < 6) syst = 0.00151752;
      if (x >= 6 && x < 7.5) syst = 0.0019946;
      if (x >= 7.5 && x < 10) syst = 0.00321667;
      if (x >= 10 && x < 15) syst = 0.00350192;
      if (x >= 15 && x < 30) syst = 0.00441023;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0108482;
      if (x >= 2 && x < 2.5) syst = 0.00318395;
      if (x >= 2.5 && x < 3) syst = 0.00449944;
      if (x >= 3 && x < 3.5) syst = 0.00272957;
      if (x >= 3.5 && x < 4) syst = 0.00779276;
      if (x >= 4 && x < 4.5) syst = 0.00622458;
      if (x >= 4.5 && x < 5.5) syst = 0.00114197;
      if (x >= 5.5 && x < 6.5) syst = 0.0077565;
      if (x >= 6.5 && x < 8) syst = 0.00157525;
      if (x >= 8 && x < 9.5) syst = 0.00373986;
      if (x >= 9.5 && x < 13) syst = 0.00445251;
      if (x >= 13 && x < 20) syst = 0.0284466;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0294959;
      if (x >= 2.2 && x < 2.7) syst = 0.00108294;
      if (x >= 2.7 && x < 3.2) syst = 0.0152933;
      if (x >= 3.2 && x < 3.7) syst = 0.00300734;
      if (x >= 3.7 && x < 4.7) syst = 0.00349077;
      if (x >= 4.7 && x < 6.5) syst = 0.0159082;
      if (x >= 6.5 && x < 8.5) syst = 0.00321438;
      if (x >= 8.5 && x < 11) syst = 0.00306389;
      if (x >= 11 && x < 20) syst = 0.0125999;
    }
  }
  if (filterId==2) { //L2 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.687436; den = 0.628026; statUp = 0.00741958; statDown = 0.00747367;}
      else if (x >= 4 && x <4.5) {num = 0.879668; den = 0.855425; statUp = 0.00516048; statDown = 0.0052745;}
      else if (x >= 4.5 && x <5) {num = 0.936813; den = 0.91093; statUp = 0.00414716; statDown = 0.00429608;}
      else if (x >= 5 && x <5.5) {num = 0.948744; den = 0.930006; statUp = 0.00405143; statDown = 0.00424179;}
      else if (x >= 5.5 && x <6.5) {num = 0.964058; den = 0.95091; statUp = 0.00286448; statDown = 0.00299563;}
      else if (x >= 6.5 && x <8) {num = 0.967607; den = 0.967695; statUp = 0.0027257; statDown = 0.00286633;}
      else if (x >= 8 && x <10.5) {num = 0.977261; den = 0.975624; statUp = 0.0024585; statDown = 0.00262625;}
      else if (x >= 10.5 && x <14) {num = 0.975199; den = 0.979933; statUp = 0.00326164; statDown = 0.00354865;}
      else if (x >= 14 && x <18) {num = 0.978416; den = 0.983164; statUp = 0.00450526; statDown = 0.00515026;}
      else {num = 0.986566; den = 0.984365; statUp = 0.00413127; statDown = 0.00501923;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.707798; den = 0.663284; statUp = 0.0115376; statDown = 0.0116219;}
      else if (x >= 3 && x <3.5) {num = 0.841787; den = 0.79966; statUp = 0.00756481; statDown = 0.00765363;}
      else if (x >= 3.5 && x <4) {num = 0.915148; den = 0.874984; statUp = 0.0058153; statDown = 0;}
      else if (x >= 4 && x <4.5) {num = 0.932092; den = 0.925267; statUp = 0.00578591; statDown = 0.0060332;}
      else if (x >= 4.5 && x <5) {num = 0.95135; den = 0.943497; statUp = 0.005532; statDown = 0.00587852;}
      else if (x >= 5 && x <6) {num = 0.968164; den = 0.958147; statUp = 0.00392457; statDown = 0.00419381;}
      else if (x >= 6 && x <7.5) {num = 0.962352; den = 0.96425; statUp = 0.00435955; statDown = 0.00466282;}
      else if (x >= 7.5 && x <10) {num = 0.979864; den = 0.969778; statUp = 0.00356325; statDown = 0.00395684;}
      else if (x >= 10 && x <15) {num = 0.983282; den = 0.973169; statUp = 0.00384401; statDown = 0.00442388;}
      else {num = 0.985422; den = 0.979848; statUp = 0.00507078; statDown = 0.00689321;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.651185; den = 0.562127; statUp = 0.0311006; statDown = 0.0309056;}
      else if (x >= 2 && x <2.5) {num = 0.75568; den = 0.744834; statUp = 0.0139774; statDown = 0.0139217;}
      else if (x >= 2.5 && x <3) {num = 0.895072; den = 0.898086; statUp = 0.00967712; statDown = 0.00997311;}
      else if (x >= 3 && x <3.5) {num = 0.911965; den = 0.932775; statUp = 0.0093472; statDown = 0.00978825;}
      else if (x >= 3.5 && x <4) {num = 0.933577; den = 0.949087; statUp = 0.008866; statDown = 0.00937714;}
      else if (x >= 4 && x <4.5) {num = 0.92098; den = 0.952348; statUp = 0.0100317; statDown = 0.0107328;}
      else if (x >= 4.5 && x <5.5) {num = 0.945015; den = 0.959924; statUp = 0.00698531; statDown = 0.00744697;}
      else if (x >= 5.5 && x <6.5) {num = 0.915225; den = 0.961517; statUp = 0.0105111; statDown = 0.0112266;}
      else if (x >= 6.5 && x <8) {num = 0.942614; den = 0.958878; statUp = 0.00845552; statDown = 0.00920477;}
      else if (x >= 8 && x <9.5) {num = 0.943992; den = 0.95865; statUp = 0.0111118; statDown = 0.0123342;}
      else if (x >= 9.5 && x <13) {num = 0.953176; den = 0.956512; statUp = 0.0100647; statDown = 0.0114277;}
      else {num = 0.9711; den = 0.959808; statUp = 0.0112001; statDown = 0.0137596;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.80688; den = 0.804533; statUp = 0.017851; statDown = 0.0179625;}
      else if (x >= 2.2 && x <2.7) {num = 0.864384; den = 0.87583; statUp = 0.0134953; statDown = 0.0137959;}
      else if (x >= 2.7 && x <3.2) {num = 0.892179; den = 0.904052; statUp = 0.0121752; statDown = 0.0125694;}
      else if (x >= 3.2 && x <3.7) {num = 0.896585; den = 0.926574; statUp = 0.0112721; statDown = 0.0117765;}
      else if (x >= 3.7 && x <4.7) {num = 0.910998; den = 0.928332; statUp = 0.00974335; statDown = 0.010226;}
      else if (x >= 4.7 && x <6.5) {num = 0.938687; den = 0.949857; statUp = 0.00795089; statDown = 0.00845509;}
      else if (x >= 6.5 && x <8.5) {num = 0.941825; den = 0.954442; statUp = 0.0109965; statDown = 0.0143609;}
      else if (x >= 8.5 && x <11) {num = 0.931535; den = 0.959109; statUp = 0.0136929; statDown = 0.0156881;}
      else {num = 0.969056; den = 0.967084; statUp = 0.00977435; statDown = 0.0126053;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00119995;
      if (x >= 4 && x < 4.5) syst = 0.000801484;
      if (x >= 4.5 && x < 5) syst = 0.00142786;
      if (x >= 5 && x < 5.5) syst = 0.000859141;
      if (x >= 5.5 && x < 6.5) syst = 0.000855793;
      if (x >= 6.5 && x < 8) syst = 0.000338442;
      if (x >= 8 && x < 10.5) syst = 0.000905661;
      if (x >= 10.5 && x < 14) syst = 0.000193737;
      if (x >= 14 && x < 18) syst = 0.000621028;
      if (x >= 18 && x < 30) syst = 0.0029276;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00301699;
      if (x >= 3 && x < 3.5) syst = 0.0051637;
      if (x >= 3.5 && x < 4) syst = 0.00271564;
      if (x >= 4 && x < 4.5) syst = 0.00128082;
      if (x >= 4.5 && x < 5) syst = 0.00105614;
      if (x >= 5 && x < 6) syst = 0.00120191;
      if (x >= 6 && x < 7.5) syst = 0.000729975;
      if (x >= 7.5 && x < 10) syst = 0.00139352;
      if (x >= 10 && x < 15) syst = 0.00151879;
      if (x >= 15 && x < 30) syst = 0.00138277;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0234362;
      if (x >= 2 && x < 2.5) syst = 0.00781699;
      if (x >= 2.5 && x < 3) syst = 0.0020642;
      if (x >= 3 && x < 3.5) syst = 0.00494294;
      if (x >= 3.5 && x < 4) syst = 0.00372959;
      if (x >= 4 && x < 4.5) syst = 0.0101533;
      if (x >= 4.5 && x < 5.5) syst = 0.00248577;
      if (x >= 5.5 && x < 6.5) syst = 0.00480156;
      if (x >= 6.5 && x < 8) syst = 0.00535204;
      if (x >= 8 && x < 9.5) syst = 0.00407749;
      if (x >= 9.5 && x < 13) syst = 0.000734987;
      if (x >= 13 && x < 20) syst = 0.0152108;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0547508;
      if (x >= 2.2 && x < 2.7) syst = 0.0439035;
      if (x >= 2.7 && x < 3.2) syst = 0.0100721;
      if (x >= 3.2 && x < 3.7) syst = 0.00486924;
      if (x >= 3.7 && x < 4.7) syst = 0.0164241;
      if (x >= 4.7 && x < 6.5) syst = 0.0045128;
      if (x >= 6.5 && x < 8.5) syst = 0.00615735;
      if (x >= 8.5 && x < 11) syst = 0.00521994;
      if (x >= 11 && x < 20) syst = 0.00496602;
    }
  }
  if (filterId==3) { //L3 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0981413; den = 0.0714076; statUp = 0.00475984; statDown = 0.00464341;}
      else if (x >= 4 && x <4.5) {num = 0.309591; den = 0.234967; statUp = 0.00731017; statDown = 0.00724988;}
      else if (x >= 4.5 && x <5) {num = 0.49696; den = 0.427491; statUp = 0.00850388; statDown = 0.00850324;}
      else if (x >= 5 && x <5.5) {num = 0.646567; den = 0.569805; statUp = 0.00897182; statDown = 0.00902901;}
      else if (x >= 5.5 && x <6.5) {num = 0.717727; den = 0.665698; statUp = 0.00690496; statDown = 0.00696443;}
      else if (x >= 6.5 && x <8) {num = 0.771046; den = 0.736859; statUp = 0.00662127; statDown = 0.00670002;}
      else if (x >= 8 && x <10.5) {num = 0.792067; den = 0.777534; statUp = 0.00684886; statDown = 0.00695048;}
      else if (x >= 10.5 && x <14) {num = 0.826589; den = 0.814236; statUp = 0.00832558; statDown = 0.00852162;}
      else if (x >= 14 && x <18) {num = 0.800339; den = 0.820918; statUp = 0.0131166; statDown = 0.0135246;}
      else {num = 0.846856; den = 0.837225; statUp = 0.0139208; statDown = 0.0145458;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.307823; den = 0.284114; statUp = 0.0111767; statDown = 0.0110334;}
      else if (x >= 3 && x <3.5) {num = 0.429139; den = 0.424849; statUp = 0.00992916; statDown = 0.00988034;}
      else if (x >= 3.5 && x <4) {num = 0.54449; den = 0.527662; statUp = 0.45551; statDown = 0.00969908;}
      else if (x >= 4 && x <4.5) {num = 0.591156; den = 0.604174; statUp = 0.0113833; statDown = 0.0114203;}
      else if (x >= 4.5 && x <5) {num = 0.639967; den = 0.645913; statUp = 0.0126515; statDown = 0.0127542;}
      else if (x >= 5 && x <6) {num = 0.673449; den = 0.679156; statUp = 0.0104608; statDown = 0.010556;}
      else if (x >= 6 && x <7.5) {num = 0.685635; den = 0.720417; statUp = 0.0109999; statDown = 0.0111263;}
      else if (x >= 7.5 && x <10) {num = 0.749011; den = 0.754465; statUp = 0.0113948; statDown = 0.0115976;}
      else if (x >= 10 && x <15) {num = 0.773253; den = 0.801728; statUp = 0.0133331; statDown = 0.0136666;}
      else {num = 0.773038; den = 0.833024; statUp = 0.0211603; statDown = 0.0220621;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.00155461; den = 0.000463755; statUp = 0.00211758; statDown = 0.00108585;}
      else if (x >= 2 && x <2.5) {num = 0.00555387; den = 0.00858885; statUp = 0.00315843; statDown = 0.00555387;}
      else if (x >= 2.5 && x <3) {num = 0.451703; den = 0.447277; statUp = 0.0150343; statDown = 0.0149314;}
      else if (x >= 3 && x <3.5) {num = 0.655861; den = 0.650893; statUp = 0.0159254; statDown = 0.0160134;}
      else if (x >= 3.5 && x <4) {num = 0.706533; den = 0.710335; statUp = 0.0158931; statDown = 0.0161213;}
      else if (x >= 4 && x <4.5) {num = 0.726282; den = 0.741482; statUp = 0.0174238; statDown = 0.017772;}
      else if (x >= 4.5 && x <5.5) {num = 0.764391; den = 0.796199; statUp = 0.013155; statDown = 0.0134088;}
      else if (x >= 5.5 && x <6.5) {num = 0.769821; den = 0.824468; statUp = 0.016319; statDown = 0.0167672;}
      else if (x >= 6.5 && x <8) {num = 0.811763; den = 0.834174; statUp = 0.0145482; statDown = 0.0150501;}
      else if (x >= 8 && x <9.5) {num = 0.819571; den = 0.841319; statUp = 0.0190366; statDown = 0.01987;}
      else if (x >= 9.5 && x <13) {num = 0.829677; den = 0.857122; statUp = 0.0192692; statDown = 0.0202277;}
      else {num = 0.874981; den = 0.874474; statUp = 0.0258047; statDown = 0.0277093;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.00366548; den = 0.000413051; statUp = 0.00211949; statDown = 0.00167583;}
      else if (x >= 2.2 && x <2.7) {num = 0.116176; den = 0.109916; statUp = 0.0108872; statDown = 0.0105043;}
      else if (x >= 2.7 && x <3.2) {num = 0.413123; den = 0.401827; statUp = 0.0181733; statDown = 0.017947;}
      else if (x >= 3.2 && x <3.7) {num = 0.482035; den = 0.491334; statUp = 0.0184417; statDown = 0.0183373;}
      else if (x >= 3.7 && x <4.7) {num = 0.568894; den = 0.573223; statUp = 0.0169342; statDown = 0.0169658;}
      else if (x >= 4.7 && x <6.5) {num = 0.675048; den = 0.651616; statUp = 0.0160907; statDown = 0.0162558;}
      else if (x >= 6.5 && x <8.5) {num = 0.722882; den = 0.711847; statUp = 0.0215971; statDown = 0.022054;}
      else if (x >= 8.5 && x <11) {num = 0.714358; den = 0.750096; statUp = 0.0275205; statDown = 0.0283679;}
      else {num = 0.753355; den = 0.779864; statUp = 0.0323455; statDown = 0.0336041;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000650706;
      if (x >= 4 && x < 4.5) syst = 0.0010869;
      if (x >= 4.5 && x < 5) syst = 0.00298052;
      if (x >= 5 && x < 5.5) syst = 0.00341277;
      if (x >= 5.5 && x < 6.5) syst = 0.000613358;
      if (x >= 6.5 && x < 8) syst = 0.000658119;
      if (x >= 8 && x < 10.5) syst = 0.000756756;
      if (x >= 10.5 && x < 14) syst = 0.000662617;
      if (x >= 14 && x < 18) syst = 0.00220571;
      if (x >= 18 && x < 30) syst = 0.00215326;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00406324;
      if (x >= 3 && x < 3.5) syst = 0.00422745;
      if (x >= 3.5 && x < 4) syst = 0.00493964;
      if (x >= 4 && x < 4.5) syst = 0.0015019;
      if (x >= 4.5 && x < 5) syst = 0.00349953;
      if (x >= 5 && x < 6) syst = 0.00165421;
      if (x >= 6 && x < 7.5) syst = 0.00195686;
      if (x >= 7.5 && x < 10) syst = 0.00305233;
      if (x >= 10 && x < 15) syst = 0.00341103;
      if (x >= 15 && x < 30) syst = 0.00425449;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.00154707;
      if (x >= 2 && x < 2.5) syst = 0.00239348;
      if (x >= 2.5 && x < 3) syst = 0.00578521;
      if (x >= 3 && x < 3.5) syst = 0.0019864;
      if (x >= 3.5 && x < 4) syst = 0.00826595;
      if (x >= 4 && x < 4.5) syst = 0.00622458;
      if (x >= 4.5 && x < 5.5) syst = 0.00155048;
      if (x >= 5.5 && x < 6.5) syst = 0.00738518;
      if (x >= 6.5 && x < 8) syst = 0.00155169;
      if (x >= 8 && x < 9.5) syst = 0.00373986;
      if (x >= 9.5 && x < 13) syst = 0.00445251;
      if (x >= 13 && x < 20) syst = 0.028681;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.000519911;
      if (x >= 2.2 && x < 2.7) syst = 0.0288676;
      if (x >= 2.7 && x < 3.2) syst = 0.013137;
      if (x >= 3.2 && x < 3.7) syst = 0.00153582;
      if (x >= 3.7 && x < 4.7) syst = 0.00361851;
      if (x >= 4.7 && x < 6.5) syst = 0.0155374;
      if (x >= 6.5 && x < 8.5) syst = 0.00321391;
      if (x >= 8.5 && x < 11) syst = 0.00306389;
      if (x >= 11 && x < 20) syst = 0.0199929;
    }
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  double toreturn = (getRawDen) ? den : ((num+syst_factor+stat_factor)/den);
  return toreturn;
}

///////////////////////////////////////////////////
//              T R K     P b P b                //
///////////////////////////////////////////////////

double Algorithm::tnp_weight_trk_pbpb(double eta, int idx) {
  double x = eta;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  //SF in eta bins
  if (x >= -2.4 && x < -1.6) {num = 0.994498; den = 0.998413; statUp = 0.00287386; statDown = 0.00290888;}
  if (x >= -1.6 && x < -1.2) {num = 0.973539; den = 0.967322; statUp = 0.00399501; statDown = 0.00407843;}
  if (x >= -1.2 && x < -0.9) {num = 0.964465; den = 0.970816; statUp = 0.00861188; statDown = 0.00869453;}
  if (x >= -0.9 && x < -0.6) {num = 0.96081; den = 0.974407; statUp = 0.0223599; statDown = 0.00682405;}
  if (x >= -0.6 && x < -0.3) {num = 0.964464; den = 0.97802; statUp = 0.00612474; statDown = 0.00621291;}
  if (x >= -0.3 && x < 0.3) {num = 0.963862; den = 0.966583; statUp = 0.00496914; statDown = 0.00503378;}
  if (x >= 0.3 && x < 0.6) {num = 0.956897; den = 0.967248; statUp = 0.00672757; statDown = 0.0068202;}
  if (x >= 0.6 && x < 0.9) {num = 0.964172; den = 0.966882; statUp = 0.00735892; statDown = 0.00754429;}
  if (x >= 0.9 && x < 1.2) {num = 0.961874; den = 0.955987; statUp = 0.00987473; statDown = 0.0099638;}
  if (x >= 1.2 && x < 1.6) {num = 0.964754; den = 0.964653; statUp = 0.0042287; statDown = 0.00430601;}
  if (x >= 1.6 && x < 2.4) {num = 0.999937; den = 0.998771; statUp = 6.33084e-05; statDown = 0.00310832;}

  // syst uncertainties
  if (x >= -2.4 && x < -1.6) syst = 0.00129036;
  if (x >= -1.6 && x < -1.2) syst = 0.00376932;
  if (x >= -1.2 && x < -0.9) syst = 0.00125496;
  if (x >= -0.9 && x < -0.6) syst = 0.00240006;
  if (x >= -0.6 && x < -0.3) syst = 0.00228604;
  if (x >= -0.3 && x < 0.3) syst = 0.00494002;
  if (x >= 0.3 && x < 0.6) syst = 0.00529974;
  if (x >= 0.6 && x < 0.9) syst = 0.0023114;
  if (x >= 0.9 && x < 1.2) syst = 0.012059;
  if (x >= 1.2 && x < 1.6) syst = 0.00278996;
  if (x >= 1.6 && x < 2.4) syst = 0.00104774;

  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}