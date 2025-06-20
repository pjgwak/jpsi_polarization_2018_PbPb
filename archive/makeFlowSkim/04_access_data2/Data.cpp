#include "Data.h"

Data::Data(TTree* tree) : m_tree(tree) {
  // m_tree->SetBranchAddress("XXX", &XXX);
  m_tree->SetBranchAddress("runNb", &runNb);
  m_tree->SetBranchAddress("Centrality", &Centrality);
  m_tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);

  // TClonesArray pointers - need "&""
  m_tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  m_tree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // fixed-size array - no "&"
  m_tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  m_tree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  m_tree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  m_tree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  m_tree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  m_tree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  m_tree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  m_tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  m_tree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
}