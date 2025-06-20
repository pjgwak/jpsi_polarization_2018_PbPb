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
  if(!passedSelection(irqq)) return;

  // HX, CS transformation
  MuPlusVector_Helicity(*data.JP_Reco, *data.mupl_Reco); // Call HX before CS
  MuPlusVector_CollinsSoper(*data.JP_Reco);

  // fill the histograms
  standByFilling(irqq);

  // count up nDimu
  // must be increased later than standByFilling()
  ++nDimu;
  // fillHists();
}
bool Algorithm::passedSelection(Long64_t irqq) {
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
  if (cut_runNb327123 && !passedRunNb327123()) return false;
  if (cut_HLTriggerPbPbJpsi2018 && !passedHLTriggerJpsiPbPb2018()) return false;
  if (cut_recoQQTrigger && !passedRecoQQTrigger(irqq)) return false;
  if (cut_L2L3FilterPbPbJpsi2018 && !passedHLFilterJpsiPbPb2018(irqq)) return false;
  if (cut_softMuons && !passedSoftMuonCut(irqq)) return false;
  if (cut_vtxProbability && data.Reco_QQ_VtxProb[irqq] < 0.01) return false;
  if (cut_oppositeSign && qq_sign_val != 0) return false;
  if (cut_singleMuonAcc && !passedMuonAcc2018()) return false;

  // passed all selection
  return true;
}

bool Algorithm::passedSoftMuonCut(Long64_t irqq) {
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

bool Algorithm::passedRunNb327123() {
  // this is event level cut
  // it will be executed in the Reco loop
  // to manange all cuts in the Algorithm()
  Data& data = *m_data;
  
  if (data.runNb < 327123)
    return true;
  else
    return false;
}

bool Algorithm::passedHLFilterJpsiPbPb2018(Long64_t irqq) {
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

bool Algorithm::passedHLTriggerJpsiPbPb2018() {
  // this is event level cut
  // it will be executed in the Reco loop
  // to manange all cuts in the Algorithm()
  Data& data = *m_data;

  return ((data.HLTriggers & ((ULong64_t)pow(2, kTrigJpsi2018))) == ((ULong64_t)pow(2, kTrigJpsi2018)));
}

bool Algorithm::passedRecoQQTrigger(Long64_t irqq) {
  Data& data = *m_data;

  return ((data.Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigJpsi2018))) == ((ULong64_t)pow(2, kTrigJpsi2018)));
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

  recoQQsign[nDimu] = current_qq_sign_val;
  cosHX[nDimu] = muplHX.CosTheta();
  phiHX[nDimu] = muplHX.Phi();
  cosCS[nDimu] = muplCS.CosTheta();
  phiCS[nDimu] = muplCS.Phi();
  // cos_cs[nDimu] = cos_cs_;
  // phi_cs[nDimu] = phi_cs_;
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
