#include "Data.h"

Data::Data(TTree *tree, bool mc, bool isEP) : isMC(mc), isEP(isEP), m_tree(tree) {
  // m_tree->SetBranchAddress("XXX", &XXX);
  if (!isMC) {
    m_tree->SetBranchAddress("runNb", &runNb);
    m_tree->SetBranchAddress("LS", &LS);
  }

  if (isMC) {
    // m_tree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
    m_tree->SetBranchAddress("Gen_weight", &Gen_weight);
  }

  if (isEP) {
    m_tree->SetBranchAddress("ep_flat", &ep_flat);
  }

  m_tree->SetBranchAddress("eventNb", &eventNb);
  m_tree->SetBranchAddress("Centrality", &Centrality);
  m_tree->SetBranchAddress("zVtx", &zVtx);
  m_tree->SetBranchAddress("HLTriggers", &HLTriggers);
  m_tree->SetBranchAddress("SumET_HF", &SumET_HF);

  // TClonesArray pointers - need "&""
  m_tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  m_tree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // fixed-size array - no "&"
  m_tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  m_tree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  m_tree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  m_tree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  m_tree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  m_tree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  m_tree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  m_tree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  m_tree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  m_tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
}

// ====== Childern classes ===== //
DataRun2::DataRun2(TTree *tree, bool isMC, bool isEP) : Data(tree, isMC, isEP) {
  m_tree->SetBranchAddress("Reco_QQ_size", &this->Reco_QQ_size);
  m_tree->SetBranchAddress("Reco_QQ_size", &this->Reco_QQ_size);
  m_tree->SetBranchAddress("Reco_QQ_mupl_idx", this->Reco_QQ_mupl_idx);
  m_tree->SetBranchAddress("Reco_QQ_mumi_idx", this->Reco_QQ_mumi_idx);
  m_tree->SetBranchAddress("Reco_QQ_sign", this->Reco_QQ_sign);

  if (isMC) 
    m_tree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
}

DataRun3::DataRun3(TTree *tree, bool isMC, bool isEP) : Data(tree, isMC, isEP) {
  m_tree->SetBranchAddress("Reco_QQ_size", &this->Reco_QQ_size);
  m_tree->SetBranchAddress("Reco_QQ_mupl_idx", this->Reco_QQ_mupl_idx);
  m_tree->SetBranchAddress("Reco_QQ_mumi_idx", this->Reco_QQ_mumi_idx);
  m_tree->SetBranchAddress("Reco_QQ_sign", this->Reco_QQ_sign);
  
  if (isMC)
    m_tree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
}