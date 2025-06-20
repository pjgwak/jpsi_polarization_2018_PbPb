#include "Algorithm.h"
#include "Data.h"

Algorithm::Algorithm() {
  // skip
}

void Algorithm::initialize(Data* data) {
  m_data = data;
  
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
  fillHists();
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

void Algorithm::fillHists() {
  if (h_runNb_m)
    h_runNb_m->Fill(m_data->runNb);
}

