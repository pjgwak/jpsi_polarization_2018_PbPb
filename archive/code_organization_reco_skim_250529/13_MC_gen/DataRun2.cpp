#include "DataRun2.h"

DataRun2::DataRun2(TTree *tree, bool isMC, bool isGenOnly, bool isEP, bool isRun2)
  : Data(tree, isMC, isGenOnly, isEP, isRun2) {
  if (isMC && m_tree->GetBranch("Gen_QQ_size")) m_tree->SetBranchAddress("Gen_QQ_size", &this->Gen_QQ_size_run2);
  if (isMC && m_tree->GetBranch("Gen_QQ_mupl_idx")) m_tree->SetBranchAddress("Gen_QQ_mupl_idx", this->Gen_QQ_mupl_idx_run2);
  if (isMC && m_tree->GetBranch("Gen_QQ_mumi_idx")) m_tree->SetBranchAddress("Gen_QQ_mumi_idx", this->Gen_QQ_mumi_idx_run2);
  if (isMC && m_tree->GetBranch("Gen_mu_charge")) m_tree->SetBranchAddress("Gen_mu_charge", this->Gen_mu_charge_run2);

  if (m_tree->GetBranch("Reco_QQ_size")) m_tree->SetBranchAddress("Reco_QQ_size", &this->Reco_QQ_size_run2);
  if (m_tree->GetBranch("Reco_mu_whichGen")) m_tree->SetBranchAddress("Reco_mu_whichGen", this->Reco_mu_whichGen_run2);
  if (m_tree->GetBranch("Reco_QQ_mupl_idx")) m_tree->SetBranchAddress("Reco_QQ_mupl_idx", this->Reco_QQ_mupl_idx_run2);
  if (m_tree->GetBranch("Reco_QQ_mumi_idx")) m_tree->SetBranchAddress("Reco_QQ_mumi_idx", this->Reco_QQ_mumi_idx_run2);
  if (m_tree->GetBranch("Reco_QQ_sign")) m_tree->SetBranchAddress("Reco_QQ_sign", this->Reco_QQ_sign_run2);
}